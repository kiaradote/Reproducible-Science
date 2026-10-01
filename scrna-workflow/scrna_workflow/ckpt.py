"""Checkpoint engine: restartable, file-based stages.

Each stage reads ONLY the previous stage's checkpoint file and writes its own checkpoint, a report and a manifest
entry. State lives in files, not in notebook memory, so any stage can be re-run from disk.
"""
import hashlib, json, os, sys, time, platform
from pathlib import Path


class StageFailed(RuntimeError):
    """Raised when a stage's report is FAIL (the checkpoint is still saved so you can inspect it)."""


def shorten(text, n):
    """Cut text at a word boundary and mark the cut."""
    return text if len(text) <= n else text[:n].rsplit(" ", 1)[0] + " ..."


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(chunk), b""):
            h.update(c)
    return h.hexdigest()


def params_hash(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


def peak_rss_mb():
    """Peak resident memory of this process so far (monotonic), in MB."""
    try:
        import psutil
        p = psutil.Process()
        mi = p.memory_info()
        return round(getattr(mi, "peak_wset", mi.rss) / 2**20)
    except Exception:
        pass
    try:
        import resource
        r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(r / 1024 if sys.platform != "darwin" else r / 2**20)
    except Exception:
        return None


def versions():
    out = {"python": platform.python_version()}
    for m in ("numpy", "pandas", "scipy", "anndata", "scanpy"):
        try:
            out[m] = __import__(m).__version__
        except Exception:
            pass
    return out


class Workflow:
    """stages: ordered list of names. params: {"global": {...}, "<stage>": {...}}.
    kinds: {"<stage>": "h5ad" | "json"} (output file type; default h5ad)."""

    WHAT = {"fetch": "the files are the intended ones and complete", "assemble": "matrix, genes and cell metadata line up and match what was expected", "qc": "damaged or empty cells are removed and the filters are sensible",
            "features": "counts are normalised and informative genes are chosen reproducibly", "dimred": "the main structure is captured without being driven by technical effects",
            "cluster": "groups of cells are stable under resampling and each is distinct", "annotate": "cell-type names are supported by marker genes and say when they are not",
            "compare": "the design can answer the condition question, and the differences are statistically valid", "interpret": "the biological conclusions hold up when analysis choices change"}
    TITLES = {"fetch": "0. Fetch", "assemble": "1. Assemble and data integrity", "qc": "2. QC and filtering",
              "features": "3. Normalization and feature selection", "dimred": "4. Dimensionality reduction", "cluster": "5. Clustering",
              "annotate": "6. Cell type annotation", "interpret": "8. Biological interpretation and validation", "compare": "7. Condition comparison"}

    def __init__(self, root, params, stages, kinds=None):
        self.root, self.params, self.stages = Path(root), params, list(stages)
        self.kinds = kinds or {}
        for d in ("checkpoints", "reports", "figures"):
            (self.root / d).mkdir(parents=True, exist_ok=True)
        self.mpath = self.root / "manifest.json"
        self.m = json.load(open(self.mpath)) if self.mpath.exists() else {}

    # ---- paths and bookkeeping
    def ckpt_path(self, stage):
        n = self.stages.index(stage)
        return self.root / "checkpoints" / f"{n:02d}_{stage}.{self.kinds.get(stage, 'h5ad')}"

    def _save_manifest(self):
        tmp = str(self.mpath) + ".tmp"
        json.dump(self.m, open(tmp, "w"), indent=1, default=str)
        os.replace(tmp, self.mpath)

    def _stage_params(self, stage):
        return {"global": self.params.get("global", {}), "stage": self.params.get(stage, {})}

    def _prev(self, stage):
        i = self.stages.index(stage)
        return self.stages[i - 1] if i > 0 else None

    # ---- state: not_run | current | stale | failed
    def state(self, stage):
        """Returns (state, reason). A stage is stale if its parameters or its input checkpoint changed, or an upstream stage is stale/failed."""
        e = self.m.get(stage)
        if e is None:
            return "not_run", "never run"
        prev = self._prev(stage)
        if prev is not None:
            ps, why = self.state(prev)
            if ps in ("not_run", "stale"):
                return "stale", f"upstream stage '{prev}' is {ps.replace('_', ' ')}"
            if e.get("input_hash") != self.m[prev]["output_hash"]:
                return "stale", f"input checkpoint of '{prev}' changed"
        if e.get("params_hash") != params_hash(self._stage_params(stage)):
            return "stale", "parameters changed"
        if not Path(e["output"]).exists():
            return "stale", "checkpoint file missing"
        return ("failed" if e["status"] == "FAIL" and not e.get("accepted") else "current"), e["status"]

    def status(self):
        import pandas as pd
        rows = []
        for s in self.stages:
            st, why = self.state(s)
            e = self.m.get(s, {})
            rows.append({"stage": s, "state": st, "report": e.get("status", ""), "why": why if st != "current" else "",
                         "accepted": bool(e.get("accepted")), "seconds": e.get("seconds"), "peak_mb": e.get("peak_mb"),
                         "checkpoint_mb": e.get("output_mb")})
        return pd.DataFrame(rows).set_index("stage")

    # ---- loading
    def load(self, stage, backed=None):
        p = self.ckpt_path(stage)
        if not p.exists():
            raise FileNotFoundError(f"no checkpoint for stage '{stage}' ({p})")
        if p.suffix == ".json":
            return json.load(open(p))
        import anndata as ad
        return ad.read_h5ad(p, backed=backed)

    # ---- control
    def invalidate(self, stage):
        """Forget this stage and everything after it (checkpoint files are deleted). Next run recomputes them."""
        for s in self.stages[self.stages.index(stage):]:
            self.m.pop(s, None)
            p = self.ckpt_path(s)
            if p.exists():
                p.unlink()
        self._save_manifest()

    def accept(self, stage, reason):
        """Record a human decision to continue past a FAIL (kept in the manifest, shown in status())."""
        assert reason and reason.strip(), "give a reason"
        self.m[stage]["accepted"] = reason
        self._save_manifest()

    # ---- run
    def run(self, stage, fn, writer, force=False):
        """fn(prev_checkpoint_or_None, stage_params, global_params) -> (output, records).
        writer(records, prefix) -> overall status string. Skips if the stage is current (unless force)."""
        st, why = self.state(stage)
        if st == "current" and not force:
            print(f"[{stage}] up to date ({self.m[stage]['status']}); loaded from checkpoint")
            return self.load(stage)
        prev = self._prev(stage)
        prev_obj = None
        if prev is not None:
            ps, pwhy = self.state(prev)
            if ps != "current":
                raise StageFailed(f"cannot run '{stage}': upstream stage '{prev}' is {ps} ({pwhy}). "
                                  + ("Inspect its report, then fix it or call accept()." if ps == "failed" else "Run it first."))
            prev_obj = self.load(prev)
        t0 = time.time()
        out, records = fn(prev_obj, self.params.get(stage, {}), self.params.get("global", {}))
        path = self.ckpt_path(stage)
        if path.suffix == ".json":
            json.dump(out, open(path, "w"), indent=1, default=str)
        else:
            out.write_h5ad(path)
        overall = writer(records, str(self.root / "reports" / stage))
        self.m[stage] = {
            "status": overall, "params_hash": params_hash(self._stage_params(stage)),
            "input_hash": self.m[prev]["output_hash"] if prev else None,
            "output": str(path), "output_hash": sha256_file(path), "output_mb": round(path.stat().st_size / 2**20, 1),
            "report": str(self.root / "reports" / (stage + ".md")), "seconds": round(time.time() - t0, 1),
            "peak_mb": peak_rss_mb(), "versions": versions(),
            "time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        self._save_manifest()
        try:
            self.summary()
        except Exception as exc:                      # the summary must never break a stage
            print(f"[{stage}] run summary not updated: {exc}")
        print(f"[{stage}] {overall} in {self.m[stage]['seconds']} s, peak {self.m[stage]['peak_mb']} MB; report: {self.m[stage]['report']}")
        if overall == "FAIL":
            raise StageFailed(f"stage '{stage}' reported FAIL. Read {self.m[stage]['report']}; the checkpoint is saved. "
                              "Fix the cause and re-run, or call accept(stage, reason) to continue deliberately.")
        return out

    # ---- one-page record of the whole run
    def summary(self, write=True):
        """Build RUN_SUMMARY.md (and run_summary.json) in the run folder: verdict and reason per stage, every check with what it tests,
        the checks that did not run, parameters, versions, timing, figures with their legends, and how to resume. Rebuilt after every stage."""
        import re
        rank = {"PASS": 0, "WARN": 1, "FAIL": 2}
        lines, rows, not_run_all, totals = [], [], [], {"PASS": 0, "WARN": 0, "FAIL": 0, "INFO": 0}
        ds = self.params.get("global", {})
        done = [s for s in self.stages if s in self.m]
        worst = max((self.m[s]["status"] for s in done), key=lambda v: rank.get(v, 0), default="not run")
        stale = [s for s in self.stages if self.state(s)[0] in ("stale", "failed")]
        ver = self.m[done[-1]].get("versions", {}) if done else {}
        lines += ["# Run summary", "",
                  f"- **Run folder:** `{self.root}`", f"- **Species:** {ds.get('species', 'not set')}",
                  f"- **Generated:** {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}",
                  f"- **Software:** " + ", ".join(f"{k} {v}" for k, v in ver.items()),
                  f"- **Stages finished:** {len(done)} of {len(self.stages)}; **worst verdict so far:** {worst}"
                  + (f"; **needs attention:** {', '.join(stale)}" if stale else ""), ""]
        fetch = self.m.get("fetch")
        if fetch and Path(fetch["output"]).exists() and Path(fetch["output"]).suffix == ".json":
            files = json.load(open(fetch["output"])).get("files", [])
            if files:
                lines += ["**Input files (SHA-256 recorded at download):**", ""] + [f"- `{f['name']}`: {f['bytes']:,} bytes, `{f['sha256']}`" for f in files] + [""]
        lines += ["## Stages at a glance", "", "| stage | state | verdict | why | time (s) | peak MB | checkpoint |", "|---|---|---|---|---|---|---|"]
        body, qual = [], []
        for s in self.stages:
            st, why = self.state(s)
            e = self.m.get(s)
            title = self.TITLES.get(s, s)
            if e is None:
                lines.append(f"| {title} | not run |  |  |  |  |  |")
                qual.append((title, self.WHAT.get(s, ""), "Not done", "this stage has not been run", f"run the notebook cell for `{s}`"))
                continue
            rp = Path(e["report"])
            js = rp.with_suffix(".json")
            info = json.load(open(js)) if js.exists() else {}
            for k in totals:
                totals[k] += info.get("counts", {}).get(k, 0)
            nr = info.get("not_run", [])
            not_run_all += [(title, x) for x in nr]
            reason = re.sub(r"\s+", " ", info.get("why", ""))[:240]
            acc = f" (accepted: {e['accepted']})" if e.get("accepted") else ""
            nwarn = info.get("counts", {}).get("WARN", 0)
            if st in ("stale", "failed"):
                rating, note = "Needs attention", f"stage is {st}: {why}"
            elif e["status"] == "FAIL":
                rating, note = ("Accepted despite a failed check", f"accepted with the reason: {e.get('accepted')}") if e.get("accepted") else ("Not reliable until fixed", shorten(reason, 140))
            elif e["status"] == "WARN":
                rating, note = f"Usable, with {nwarn} caveat(s)", shorten(reason, 140)
            elif info.get("counts") and info["counts"].get("PASS", 0) == 0 and info["counts"].get("INFO", 0) > 0:
                rating, note = "Not applicable", "no scored check ran for this dataset: " + shorten(reason, 100)
            else:
                rating, note = "Good", "every scored check passed"
            qual.append((title, self.WHAT.get(s, ""), rating, note, f"`wf.invalidate(\"{s}\")` then rerun; checkpoint `{Path(e['output']).name}`"))
            lines.append(f"| {title} | {st} | {e['status']}{acc} | {reason} | {e.get('seconds')} | {e.get('peak_mb')} | `{Path(e['output']).name}` ({e.get('output_mb')} MB) |")
            rows.append({"stage": s, "state": st, "verdict": e["status"], "why": info.get("why", ""), "counts": info.get("counts", {}), "not_run": nr,
                         "seconds": e.get("seconds"), "peak_mb": e.get("peak_mb"), "checkpoint": e["output"], "accepted": e.get("accepted")})
            txt = rp.read_text(encoding="utf-8") if rp.exists() else "(report file missing)"
            txt = "\n".join(("### " + l[2:]) if l.startswith("# ") else l for l in txt.splitlines())
            prm = json.dumps(self.params.get(s, {}), indent=1, default=str)
            body += [f"## {title}", "", f"Parameters used (stage-specific): ", "", "```json", prm, "```", "", txt, ""]
        gi = lines.index("## Stages at a glance")
        nq = {k: sum(q[2].startswith(k) for q in qual) for k in ("Good", "Usable", "Not reliable", "Needs attention", "Accepted", "Not done", "Not applicable")}
        if nq["Not reliable"] or nq["Needs attention"]:
            overall_q = "Not reliable yet: fix the stages marked below before using the results"
        elif nq["Not done"]:
            overall_q = "Incomplete: some stages have not been run"
        elif nq["Usable"] or nq["Accepted"]:
            overall_q = "Usable, with caveats: read the caveats before drawing conclusions"
        else:
            overall_q = "Good: every scored check passed"
        qlines = ["## Quality check", "",
                  f"**Overall: {overall_q}.** {nq['Good']} stage(s) good, {nq['Usable']} usable with caveats, {nq['Accepted']} accepted despite a failed check, {nq['Not reliable'] + nq['Needs attention']} not reliable or needing attention, {nq['Not done']} not done, {nq['Not applicable']} not applicable.", "",
                  "How to read it: *Good* means every scored check met its default threshold; *Usable* means at least one check warned and the caveat applies to everything built on that stage; *Not reliable* means a check failed and nothing after it should be trusted until it is fixed. "
                  "Thresholds are defaults chosen for typical data, not guarantees of a correct biological answer, and a check that did not run counts for nothing (see the last section).", "",
                  "| stage | what it makes sure of | quality | why | to go back to this point |", "|---|---|---|---|---|"]
        qlines += [f"| {a} | {b} | {c} | {d.replace('|', '/')} | {e_} |" for a, b, c, d, e_ in qual]
        lines[gi:gi] = qlines + [""]
        lines += ["", f"**Checks across all stages:** {totals['PASS']} passed, {totals['WARN']} warned, {totals['FAIL']} failed, {totals['INFO']} informational. "
                  "A WARN is a point to review, not an automatic fault; a FAIL stops the run until fixed or accepted with a reason.", ""]
        lines += body
        if not_run_all:
            lines += ["## Checks that did not run", "", "A check that did not run has not passed. Each is listed with the stage that owns it.", ""] + [f"- {t}: {x}" for t, x in not_run_all] + [""]
        figs = sorted((self.root / "figures").rglob("*.png"))
        if figs:
            lines += ["## Figures", ""]
            seen = set()
            for f in figs:
                d = f.parent
                if d not in seen:
                    seen.add(d)
                    lines += [f"**{d.relative_to(self.root).as_posix()}/**", ""]
                    lg = d / "figure_legends.md"
                    if lg.exists():
                        lines += ["Legends: see `" + lg.relative_to(self.root).as_posix() + "`", ""]
                rel = f.relative_to(self.root).as_posix()
                capf = d / "figure_captions.json"
                cap = json.load(open(capf, encoding="utf-8")).get(f.name, "") if capf.exists() else ""
                lines += [f"![{f.stem}]({rel})", "", f"*{f.stem.replace('_', ' ')}: {cap}*" if cap else f"*{f.stem.replace('_', ' ')}*", ""]
        lines += ["## Folder map", "", "- `checkpoints/`: one restartable file per stage (`NN_stage.h5ad`, or `.json` for fetch)",
                  "- `reports/`: `stage.md` (read this) and `stage.json` (same content for programs) for every stage",
                  "- `figures/`: plots, with `figure_legends.md` beside each stage's figures", "- `manifest.json`: parameters hash, input and output hashes, timing, memory, versions per stage",
                  "- `RUN_SUMMARY.md`: this file, rebuilt after every stage", "",
                  "## Going back and resuming", "",
                  "- Rerun the notebook: finished stages whose parameters and input checkpoint are unchanged load from disk; the first changed stage and everything after it recompute.",
                  "- Force one stage: `run_stage(\"<stage>\", ..., force=True)`; forget a stage and all later ones: `wf.invalidate(\"<stage>\")`.",
                  "- A stage that reported FAIL stops the run with its checkpoint kept; fix the cause, or continue on purpose with `wf.accept(\"<stage>\", \"reason\")` (the reason is recorded above).",
                  "- Stage state tracks parameters and input files, not code: after editing code, `invalidate` the stage."]
        md = "\n".join(lines) + "\n"
        if write:
            (self.root / "RUN_SUMMARY.md").write_text(md, encoding="utf-8")
            json.dump({"stages": rows, "totals": totals, "worst": worst, "needs_attention": stale}, open(self.root / "run_summary.json", "w"), indent=1, default=str)
        return md
