
import re

SAFE_COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#4D4D4D")
BAD_CMAPS = ("jet", "rainbow", "hsv", "gist_rainbow", "nipy_spectral", "gist_ncar", "turbo")
NEON_NAMES = ("cyan", "aqua", "yellow", "lime", "magenta", "fuchsia")
NEON_HEX = ("#00ffff", "#ffff00", "#00ff00", "#ff00ff", "#0ff", "#ff0", "#0f0", "#f0f")


def apply_figure_rules(font_size=10):
    """Set plotting defaults: 2D only, readable fonts, safe color cycle, no top/right spines.
    Call once per deliverable figure, then build with fig, ax = plt.subplots()."""
    import matplotlib as mpl
    from cycler import cycler
    mpl.rcParams.update({
        "font.size": font_size, "axes.labelsize": font_size + 1, "axes.titlesize": font_size + 3,
        "axes.titleweight": "bold", "axes.labelweight": "normal", "axes.titlepad": 8,
        "figure.titlesize": font_size + 4, "figure.titleweight": "bold",
        "xtick.labelsize": font_size, "ytick.labelsize": font_size, "legend.fontsize": font_size,
        "legend.title_fontsize": font_size,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "figure.dpi": 100, "savefig.dpi": 300, "savefig.bbox": "tight",
        "axes.prop_cycle": cycler(color=list(SAFE_COLORS)),
    })
    return SAFE_COLORS


def title_bottom(fig, text, font_size=13):
    """Caption-style bold title placed at the BOTTOM of the figure (T2). Use instead of a top title, never both."""
    fig.text(0.5, -0.02, text, ha="center", va="top", fontsize=font_size, fontweight="bold")


def categorical_colors(n):
    """Distinct, dark-enough hues for unordered groups. Refuses more than 6 (rules K1, S1)."""
    if n > len(SAFE_COLORS):
        raise ValueError("More than 6 categories: split into small multiples, bin, or highlight a subset (K1, S1).")
    return list(SAFE_COLORS[:n])


def colormap_for(kind):
    """kind: 'sequential' (ordered magnitude), 'density' (counts per bin), 'diverging' (change around a meaningful center)."""
    table = {"sequential": "Blues", "density": "Greys", "diverging": "RdBu_r"}
    if kind not in table:
        raise ValueError("kind must be sequential, density or diverging; unordered groups use categorical_colors (K3)")
    return table[kind]


def symmetric_limits(values, center=0.0):
    """vmin, vmax symmetric around the center so a diverging scale is anchored at the reference (K3, H6)."""
    import numpy as np
    m = float(np.nanmax(np.abs(np.asarray(values, dtype=float) - center)))
    return center - m, center + m


def bubble_area(values, max_area=600.0):
    """Marker sizes for scatter(s=...) so AREA is proportional to value (H4). Values must be >= 0."""
    import numpy as np
    v = np.asarray(values, dtype=float)
    if np.nanmin(v) < 0:
        raise ValueError("Area encoding needs non-negative values")
    return v / np.nanmax(v) * max_area


def break_at_gaps(x, y, max_gap):
    """Insert NaN where consecutive x differ by more than max_gap, so a line is not drawn across unmeasured intervals (H7, H8)."""
    import numpy as np
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    gx, gy = [x[0]], [y[0]]
    for i in range(1, len(x)):
        if x[i] - x[i - 1] > max_gap:
            gx.append((x[i] + x[i - 1]) / 2.0)
            gy.append(float("nan"))
        gx.append(x[i])
        gy.append(y[i])
    return np.array(gx), np.array(gy)


def figure_spec_template():
    """Markdown spec to fill in BEFORE plotting (P1). The reviewer checks the figure against it."""
    return """# Figure spec
- Purpose: exploration | explanation   (complexity budget differs, S1)
- One question / claim (a sentence the data and design support, E7):
- Data source (saved results file, never in-memory variables, P2):
- Variable types (categorical / ordered / continuous / time) and units:
- Chart type and why it fits the question (C1):
- x, y, color, size, facet columns:
- n per group, and what one n is (cells vs cultures vs subjects) (E1, E3):
- Denominators / exclusions / transformations (E1, H6):
- Pairing or grouping structure (E5):
- Error bars or intervals, defined (SD / SE / 95% CI) (E3):
- Reference or control to show (E6):
- Axis ranges and baseline, with reason (H3, H9):
- Missing data handling (H8):
- Color plan: variable type -> palette (K1-K7):
"""


def lint_rules():
    return (
    ("H2", "FAIL", r"projection\s*=\s*['\"]3d['\"]|Axes3D|plot_surface|bar3d|plot_wireframe|plot3D|scatter3D", "3D plot: avoid 3D entirely"),
    ("H1", "FAIL", r"\.pie\([^)]*(explode|shadow)", "Pie with explode/shadow: pies must be flat 2D"),
    ("H3", "FAIL", r"\.bar[h]?\([\s\S]*?\)[\s\S]{0,400}?set_[yx]lim\(\s*[1-9]|set_[yx]lim\(\s*[1-9][\s\S]{0,400}?\.bar[h]?\(", "Bars with nonzero axis start: bar length must start at zero"),
    ("K3", "FAIL", r"cmap\s*=\s*['\"](%s)['\"]|\.cm\.(%s)\b" % ("|".join(BAD_CMAPS), "|".join(BAD_CMAPS)), "Rainbow-type colormap: use sequential or diverging"),
    ("K2", "FAIL", r"['\"](%s)['\"]|['\"](%s)['\"]" % ("|".join(NEON_NAMES), "|".join(NEON_HEX)), "Neon color (cyan/yellow/lime/magenta): hard to read on white"),
    ("H8", "WARN", r"\.fillna\(|\.interpolate\(", "Filling or interpolating values: missing is not measured; justify and mark estimated values"),
    ("P3", "WARN", None, "Random numbers used without a seed"),
    ("S2", "WARN", None, "Missing axis labels (need name and units)"),
    ("E3", "WARN", None, "Error bars without SD/SE/CI definition in the script (put it in the axis, legend or caption text)"),
    ("H9", "WARN", None, "Multiple panels without shared axes (sharex/sharey): comparisons need a common axis"),
    ("P2", "WARN", None, "Output not saved as both a vector file (svg/pdf) and a png"),
    ("H6", "WARN", r"set_[xy]scale\(\s*['\"]log|loglog|semilog|\.set\(.*scale\s*=\s*['\"]log", "Log scale: label it, and consider rendering linear and log for the user to choose"),
    ("H1b", "WARN", r"\.pie\(", "Pie chart: only for a few parts of a whole with a stated total; confirm bars would not be clearer"),
    ("E4", "WARN", None, "Bars or means plotted without raw observations (points, strip, violin) in the script"),
    ("hyg", "WARN", r"(^|\s)plt\.(plot|bar|scatter|hist|savefig)\(", "Use fig, ax = plt.subplots() and fig.savefig(), not bare plt calls (lineage)"),
)


def lint_script(src):
    """Static checks on a plotting script (text or a path to a .py file). Returns a list of dicts {rule, severity, message, line}.
    Code-level only: it cannot judge the rendered image, the data, or whether the chart type fits the question."""
    import os
    if isinstance(src, str) and os.path.exists(src):
        src = open(src, encoding="utf8").read()
    code = "\n".join(l.split("#")[0] if "#" in l and l.strip().startswith("#") else l for l in src.splitlines())
    out = []

    def add(rule, sev, msg, m=None):
        line = code[:m.start()].count("\n") + 1 if m is not None else None
        out.append({"rule": rule, "severity": sev, "message": msg, "line": line})

    for rule, sev, pat, msg in lint_rules():
        if pat:
            for m in re.finditer(pat, code, flags=re.I if rule == "K2" else 0):
                add(rule, sev, msg, m)
                break
    plots = re.search(r"\.(plot|bar|barh|scatter|hist|errorbar|imshow|pie|boxplot|violinplot)\(", code)
    if re.search(r"np\.random\.|random\.|\.sample\(", code) and not re.search(r"seed|default_rng\(\s*\d|RandomState\(\s*\d", code):
        add("P3", "WARN", "Random numbers used without a seed")
    if plots and not (re.search(r"set_xlabel|xlabel\(|set\(.*xlabel", code) and re.search(r"set_ylabel|ylabel\(|set\(.*ylabel", code)) and not re.search(r"\.pie\(", code):
        add("S2", "WARN", "Missing axis labels (need name and units)")
    if re.search(r"errorbar\(|yerr\s*=|xerr\s*=", code) and not re.search(r"\b(SD|SEM|SE|CI|s\.d\.|std|confidence)\b", code, flags=re.I):
        add("E3", "WARN", "Error bars without SD/SE/CI definition in the script")
    if re.search(r"subplots\(\s*[^)]*(ncols|nrows)\s*=\s*[2-9]|subplots\(\s*[1-9]\s*,\s*[2-9]|subplots\(\s*[2-9]\s*,\s*[1-9]", code) and not re.search(r"sharex|sharey|set_xlim|set_ylim", code):
        add("H9", "WARN", "Multiple panels without shared axes (sharex/sharey)")
    saves = re.findall(r"savefig\([^)]*\)", code)
    if saves and not (re.search(r"\.(svg|pdf)", " ".join(saves)) and re.search(r"\.png", " ".join(saves))):
        add("P2", "WARN", "Output not saved as both a vector file (svg/pdf) and a png")
    if re.search(r"\.bar[h]?\(|barplot\(", code) and not re.search(r"scatter\(|stripplot|swarmplot|violin|boxplot|\.plot\(.*marker|jitter", code):
        add("E4", "WARN", "Bars or means plotted without raw observations in the script")
    for m in re.finditer(r"(?<![A-Za-z_])(fontsize|labelsize|titlesize|font_size)\s*=\s*([0-7](?:\.\d+)?)\b", code):
        if float(m.group(2)) < 8:
            add("T1", "FAIL", "Font size %s pt is below the 8 pt floor" % m.group(2), m)
            break
    if plots and not re.search(r"set_title\(|suptitle\(|title_bottom\(|\.title\(|fig\.text\(", code):
        add("T2", "WARN", "No title found (top or bottom, bold, larger than the axis labels); the caption must then carry the claim")
    if re.search(r"set_title\([^)]*fontweight\s*=\s*['\"](normal|light)", code):
        add("T2", "WARN", "Title is not bold")
    seen, res = set(), []
    for f in out:
        k = (f["rule"], f["message"])
        if k not in seen:
            seen.add(k)
            res.append(f)
    return res


def audit_figure(fig, min_font=8.0, max_legend=6):
    """Structural checks on a LIVE matplotlib figure (before saving). Returns list of {rule, severity, message}."""
    import numpy as np
    out = []
    axes = fig.get_axes()
    for ax in axes:
        if getattr(ax, "name", "") == "3d":
            out.append({"rule": "H2", "severity": "FAIL", "message": "3D axes present"})
        bars = [p for p in ax.patches if p.__class__.__name__ == "Rectangle" and p.get_width() != 0 and p.get_height() != 0]
        if len(bars) > 1:
            lo, hi = ax.get_ylim()
            if lo > 0 or hi < 0:
                out.append({"rule": "H3", "severity": "FAIL", "message": "Bar axis does not include zero (y-limits %.3g..%.3g)" % (lo, hi)})
        if ax.get_visible() and ax.has_data():
            if not ax.get_ylabel() or not ax.get_xlabel():
                if not any(p.__class__.__name__ == "Wedge" for p in ax.patches):
                    out.append({"rule": "S2", "severity": "WARN", "message": "Axis without label (x='%s', y='%s')" % (ax.get_xlabel(), ax.get_ylabel())})
        sizes = [t.get_fontsize() for t in ax.get_xticklabels() + ax.get_yticklabels() if t.get_text()]
        if sizes and min(sizes) < min_font:
            out.append({"rule": "S2", "severity": "WARN", "message": "Tick label font %.1f pt is below %.0f pt" % (min(sizes), min_font)})
        leg = ax.get_legend()
        if leg is not None and len(leg.get_texts()) > max_legend:
            out.append({"rule": "S1", "severity": "WARN", "message": "Legend has %d entries (> %d): split, bin or label directly" % (len(leg.get_texts()), max_legend)})
        for im in ax.images + ax.collections:
            cm = getattr(im, "get_cmap", None)
            if cm is not None:
                try:
                    if im.get_cmap().name in BAD_CMAPS:
                        out.append({"rule": "K3", "severity": "FAIL", "message": "Rainbow-type colormap '%s'" % im.get_cmap().name})
                except Exception:
                    pass

    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for ax in axes:
        if not ax.get_visible():
            continue
        t = ax.title
        lab = [x for x in (ax.xaxis.label, ax.yaxis.label) if x.get_text()]
        if t.get_text() and lab:
            ls = max(x.get_fontsize() for x in lab)
            if t.get_fontsize() <= ls:
                out.append({"rule": "T2", "severity": "FAIL", "message": "Title (%.1f pt) is not larger than the axis labels (%.1f pt)" % (t.get_fontsize(), ls)})
            if t.get_fontweight() not in ("bold", "heavy", "semibold", 600, 700, 800, 900):
                out.append({"rule": "T2", "severity": "WARN", "message": "Title is not bold"})
        if t.get_text() and t.get_fontsize() < min_font:
            out.append({"rule": "T1", "severity": "FAIL", "message": "Title font below %.0f pt" % min_font})
        for x in lab:
            if x.get_fontsize() < min_font:
                out.append({"rule": "T1", "severity": "FAIL", "message": "Axis label font %.1f pt is below %.0f pt" % (x.get_fontsize(), min_font)})
                break
        tl = [x.get_fontsize() for x in ax.get_xticklabels() + ax.get_yticklabels() if x.get_text()]
        if tl and lab and max(x.get_fontsize() for x in lab) < max(tl):
            out.append({"rule": "T3", "severity": "WARN", "message": "Axis labels are smaller than tick labels"})
        leg = ax.get_legend()
        if leg is not None:
            fs = [x.get_fontsize() for x in leg.get_texts()]
            if fs and min(fs) < min_font:
                out.append({"rule": "T4", "severity": "FAIL", "message": "Legend text %.1f pt is below %.0f pt" % (min(fs), min_font)})
            bb = leg.get_window_extent(rend)
            hit = False
            for ln in ax.lines:
                pts = ax.transData.transform(np.column_stack([ln.get_xdata(orig=False), ln.get_ydata(orig=False)]).astype(float)) if len(ln.get_xdata()) and not isinstance(ln.get_xdata()[0], str) else []
                if len(pts) and ((pts[:, 0] > bb.x0) & (pts[:, 0] < bb.x1) & (pts[:, 1] > bb.y0) & (pts[:, 1] < bb.y1)).any():
                    hit = True
            for c in ax.collections:
                try:
                    off = np.asarray(c.get_offsets(), dtype=float)
                    pts = c.get_offset_transform().transform(off)
                    if len(pts) and ((pts[:, 0] > bb.x0) & (pts[:, 0] < bb.x1) & (pts[:, 1] > bb.y0) & (pts[:, 1] < bb.y1)).any():
                        hit = True
                except Exception:
                    pass
            if hit:
                out.append({"rule": "T4", "severity": "WARN", "message": "Legend covers plotted data: move it outside the axes or to an empty area"})
    for tx in fig.texts:
        if tx.get_text() and tx.get_fontsize() < min_font:
            out.append({"rule": "T1", "severity": "FAIL", "message": "Figure text below %.0f pt" % min_font})
    if len(axes) > 1:
        xl = {tuple(np.round(a.get_xlim(), 6)) for a in axes if a.has_data()}
        yl = {tuple(np.round(a.get_ylim(), 6)) for a in axes if a.has_data()}
        if len(xl) > 1 or len(yl) > 1:
            out.append({"rule": "H9", "severity": "WARN", "message": "Panels have different axis limits: confirm they are not meant to be compared, otherwise share them"})
    return out


def review_card(findings):
    """Render findings as a plain-English markdown table (rule, verdict, reason)."""
    if not findings:
        return "No code-level findings."
    rows = ["| Rule | Verdict | Finding | Line |", "|---|---|---|---|"]
    for f in findings:
        rows.append("| %s | %s | %s | %s |" % (f["rule"], f["severity"], f["message"], f.get("line") or ""))
    return "\n".join(rows)
