"""Shared helpers for the general data-analysis workflow: test records with plain-English reasons, quality ratings, chat cards, reports,
a run summary and restart logic. Standard library only."""
import hashlib
import json
import os

STATUS_RANK = {"PASS": 0, "INFO": 0, "WARN": 1, "FAIL": 2}
CHECKPOINTS = ("intake", "qc", "normalize", "structure", "compare", "biology", "interpret")
CHECKPOINT_TITLES = {"intake": "1. Intake and integrity", "qc": "2. QC and filtering", "normalize": "3. Normalization and features",
                     "structure": "4. Structure discovery", "compare": "5. Group comparison", "biology": "6. Biological meaning", "interpret": "7. Interpretation and robustness"}

def da_record(check, status, what, found, why, group=None):
    """One test. what = the plain-English definition of what was checked; found = the numbers; why = the reason the result is this status
    (the rule or threshold and how the data compared). A reason is required so no result is ever given without it."""
    if status not in STATUS_RANK:
        raise ValueError(f"status must be PASS, WARN, FAIL or INFO, got {status!r}")
    for name, text in (("what", what), ("found", found), ("why", why)):
        if not str(text).strip():
            raise ValueError(f"da_record: '{name}' must say something (check: {check})")
    r = {"check": check, "status": status, "what": str(what), "found": str(found), "why": str(why)}
    if group:
        r["group"] = group
    return r

def da_counts(records):
    return {s: sum(r["status"] == s for r in records) for s in ("PASS", "WARN", "FAIL", "INFO")}

def da_overall(records):
    """PASS, WARN or FAIL: the worst status among the records (INFO does not count)."""
    worst = max((STATUS_RANK[r["status"]] for r in records), default=0)
    return ["PASS", "WARN", "FAIL"][worst]

def da_rating(records=None, accepted_reason=None, ran=True):
    """Quality rating in words and a one-line note. Not done / Not applicable / Accepted despite a failed check / Not reliable until fixed /
    Usable, with N caveat(s) / Good."""
    if not ran or records is None:
        return "Not done", "this checkpoint has not been run"
    c = da_counts(records)
    if c["PASS"] + c["WARN"] + c["FAIL"] == 0:
        return "Not applicable", "no scored check could run for this data"
    if c["FAIL"]:
        if accepted_reason:
            return "Accepted despite a failed check", f"continued on purpose: {accepted_reason}"
        return "Not reliable until fixed", f"{c['FAIL']} check(s) failed"
    if c["WARN"]:
        return f"Usable, with {c['WARN']} caveat(s)", f"{c['WARN']} check(s) warned"
    return "Good", "every scored check passed"

def da_card(title, asks, records, conclusion, not_run=None, caveats=None, go_back=None, max_pass=5, accepted_reason=None):
    """The plain-English checkpoint card for the chat, as markdown. Shows every WARN and FAIL and up to max_pass PASS rows (INFO rows after them),
    says how many tests are not shown, then the conclusion, what did not run, carried caveats and how to go back."""
    overall = da_overall(records)
    rating, _ = da_rating(records, accepted_reason)
    order = {"FAIL": 0, "WARN": 1, "PASS": 2, "INFO": 3}
    shown = [r for r in records if r["status"] in ("FAIL", "WARN")]
    shown += [r for r in records if r["status"] == "PASS"][:max_pass]
    shown += [r for r in records if r["status"] == "INFO"][:2]
    shown.sort(key=lambda r: order[r["status"]])
    hidden = len(records) - len(shown)
    L = [f"**{title}: {overall}, {rating[0].lower() + rating[1:]}**", f"*Asks:* {asks}", "",
         "| What it checked | What we found | Result and why |", "|---|---|---|"]
    for r in shown:
        L.append(f"| {r['what']} | {r['found']} | {r['status']}: {r['why']} |".replace("\n", " "))
    if hidden > 0:
        c = da_counts([r for r in records if r not in shown])
        L += ["", f"{hidden} more test(s) not shown ({c['PASS']} passed, {c['INFO']} informational); see the report file."]
    L += ["", f"*Conclusion:* {conclusion}"]
    if not_run:
        L.append("*Did not run:* " + "; ".join(not_run if isinstance(not_run, (list, tuple)) else [not_run]))
    if caveats:
        L.append("*Carried caveats:* " + "; ".join(caveats))
    if go_back:
        L.append(f"*Go back:* {go_back}")
    return "\n".join(L)

def da_write_report(records, catalog, prefix, title, asks="", conclusion="", caveats=None, accepted_reason=None):
    """Write <prefix>.json and <prefix>.md. catalog = {check name: plain definition} of every check this checkpoint can run; checks absent from the
    records are listed as Not run. Returns the overall status."""
    overall = da_overall(records)
    rating, note = da_rating(records, accepted_reason)
    seen = {r.get("group") or r["check"] for r in records}
    not_run = {k: v for k, v in (catalog or {}).items() if k not in seen}
    counts = da_counts(records)
    with open(prefix + ".json", "w", encoding="utf8") as f:
        json.dump({"overall": overall, "rating": rating, "note": note, "counts": counts, "not_run": list(not_run), "checks": records}, f, indent=1)
    L = [f"# {title}: {overall}", "", f"**Quality rating:** {rating}. {note}.", ""]
    if asks:
        L += [f"**Asks:** {asks}", ""]
    L += [f"**Tests run ({len(records)}):** {counts['PASS']} pass, {counts['WARN']} warn, {counts['FAIL']} fail, {counts['INFO']} info", "",
          "| status | test | what it checked | what we found | result and why |", "|---|---|---|---|---|"]
    for r in records:
        L.append(f"| {r['status']} | {r['check']} | {r['what']} | {r['found']} | {r['why']} |".replace("\n", " "))
    if conclusion:
        L += ["", f"**Conclusion:** {conclusion}"]
    if not_run:
        L += ["", "**Not run:**", ""] + [f"- {k}: {v}" for k, v in not_run.items()]
    if caveats:
        L += ["", "**Carried caveats:**", ""] + [f"- {c}" for c in caveats]
    with open(prefix + ".md", "w", encoding="utf8") as f:
        f.write("\n".join(L) + "\n")
    return overall

def da_run_summary(entries, prefix):
    """entries: list of dicts with stage, rating, note, go_back (stages in order; a missing stage should be passed with rating 'Not done').
    Writes <prefix>.md with the quality table at the top and returns the overall run rating."""
    worst_order = ["Good", "Not applicable", "Usable", "Accepted", "Not done", "Not reliable"]

    def rank(rating):
        for i, k in enumerate(worst_order):
            if rating.startswith(k):
                return i
        return len(worst_order)

    ranks = [rank(e["rating"]) for e in entries]
    if any(r == 5 for r in ranks):
        overall = "Not reliable yet: fix the checkpoints marked below before using the results"
    elif any(r == 4 for r in ranks):
        overall = "Incomplete: some checkpoints have not been run"
    elif any(r in (2, 3) for r in ranks):
        overall = "Usable, with caveats: read the caveats before drawing conclusions"
    else:
        overall = "Good: every scored check passed"
    L = ["# Run summary", "", "## Quality check", "", f"**Overall: {overall}.**", "",
         "| checkpoint | quality | why | to go back to this point |", "|---|---|---|---|"]
    for e in entries:
        L.append(f"| {CHECKPOINT_TITLES.get(e['stage'], e['stage'])} | {e['rating']} | {e.get('note', '')} | {e.get('go_back', '')} |".replace("\n", " "))
    with open(prefix + ".md", "w", encoding="utf8") as f:
        f.write("\n".join(L) + "\n")
    return overall

def da_hash_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def da_hash_params(params):
    return hashlib.sha256(json.dumps(params, sort_keys=True, default=str).encode("utf8")).hexdigest()

def da_is_current(manifest, stage, params, input_hash=None):
    """Restart logic. manifest = {stage: {params_hash, input_hash, output, output_hash}}. A checkpoint is current when it was run with the same parameters on the
    same input and its output file still has the recorded hash. Returns (True/False, reason)."""
    m = manifest.get(stage)
    if not m:
        return False, "never run"
    if m.get("params_hash") != da_hash_params(params):
        return False, "parameters changed"
    if m.get("input_hash") != input_hash:
        return False, "input changed"
    out = m.get("output")
    if not out or not os.path.exists(out):
        return False, "output file missing"
    if da_hash_file(out) != m.get("output_hash"):
        return False, "output file changed since it was written"
    return True, "current"
