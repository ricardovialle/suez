"""Build the suez documentation pages in doc/ from doc/src/.

    python3 doc/src/build.py                 # doc/index.html, usage.html, performance.html, age_varying_eqtl.html
    python3 doc/src/build.py --artifact DIR  # the same site as one multi-page artifact (index without a skeleton)

Pages are fragments in doc/src/pages/ with placeholders:
    {{EX:name}}    the `## @ex name` block of doc/src/examples.R (run `Rscript doc/src/examples.R` to check them)
    {{DATA:name}}  JSON made by a loader below from the benchmark results in bench/; "usage" is the plot data of
                   the usage page, written by `Rscript doc/src/usage_data.R` (which runs examples.R)
    {{TOC}}        an index of the page's <h2 id> / <h3 id> headings
    {{JS:name}}    doc/src/name.js inlined (charts.js is also {{CHARTS_JS}})
"""
import csv, html, json, os, re, shutil, sys

SRC = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(SRC, "..", ".."))
DOC = os.path.join(ROOT, "doc")

PAGES = [("index", "Overview"), ("usage", "Usage"), ("performance", "Performance"), ("age_varying_eqtl", "Age-varying eQTLs")]
EXTERNAL = [("suez.html", "Vignette"), ("https://github.com/ricardovialle/suez", "GitHub")]
FONTS = ('<link rel="preconnect" href="https://fonts.googleapis.com">\n<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
         'family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=Space+Grotesk:wght@500;600&display=swap">')


# ------------------------------------------------------------------ examples
def examples():
    out, name, buf = {}, None, []
    for line in open(os.path.join(SRC, "examples.R")):
        m = re.match(r"## @ex (\S+)", line)
        if m: name, buf = m.group(1), []; continue
        if line.startswith("## @end"):
            if name: out[name] = "".join(buf).rstrip("\n")
            name = None; continue
        if name: buf.append(line)
    return out

def highlight_r(code):
    lines = []
    for ln in code.split("\n"):
        # comment: a # outside a string literal
        q, cut = None, None
        for i, ch in enumerate(ln):
            if q:
                if ch == "\\": continue
                if ch == q: q = None
            elif ch in "\"'": q = ch
            elif ch == "#": cut = i; break
        code_part, comment = (ln, "") if cut is None else (ln[:cut], ln[cut:])
        lines.append(html.escape(code_part) + (f'<span class="c">{html.escape(comment)}</span>' if comment else ""))
    return "<pre><code>" + "\n".join(lines) + "</code></pre>"


# ------------------------------------------------------------------ data loaders
def tsv(path):
    with open(os.path.join(ROOT, path)) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    def conv(v):
        try: return float(v)
        except (TypeError, ValueError): return v
    return [{k: conv(v) for k, v in r.items()} for r in rows]

def data_age():
    r = json.load(open(os.path.join(ROOT, "bench", "age_varying", "results", "sim_age.json")))
    def thin(a, n=240):   # QQ points: n quantiles plus the 40 largest
        a = sorted(a); m = len(a)
        idx = sorted(set([int(round(i * (m - 1) / (n - 1))) for i in range(n)] + list(range(max(0, m - 40), m))))
        return [[i, round(a[i], 4)] for i in idx]
    out = {k: r[k] for k in ["grid", "mu_true", "n_ind", "n_samples", "reps", "n_snp", "sec_linear", "sec_spline", "calib"]}
    out["shapes"] = {s: {k: v for k, v in x.items() if not k.startswith("qq")} for s, x in r["shapes"].items()}
    c = r["shapes"]["constant"]   # no G x age anywhere: every SNP is null
    out["qq"] = {"m": len(c["qq_linear"]), "linear": thin(c["qq_linear"]), "spline": thin(c["qq_spline"])}
    return out

def data_perf():
    tq = "bench/tensorqtl/results/"
    env = dict(l.rstrip("\n").split(": ", 1) for l in open(os.path.join(ROOT, tq, "environment.txt")) if ": " in l)
    base = {(r["dataset"], r["step"]): r for r in tsv("bench/results/final_summary.tsv") if r["cores"] == 1}
    r12 = [r for r in tsv("bench/results/r12_final_summary.tsv")]
    hist = []
    for r in r12:
        b = base.get((r["dataset"], r["step"]))
        if r["cores"] != 1 or b is None: continue
        hist.append({"dataset": r["dataset"], "step": r["step"], "v04": b["wall_s_base"], "v05": b["wall_s_final"],
                     "now": r["wall_c4"], "rss04": b["maxrss_mb_base"], "rss_now": r["rss_c4"]})
    cores = [{"dataset": r["dataset"], "cores": int(r["cores"]), "wall": r["wall_c4"], "tree_rss": r["tree_rss_c4"]}
             for r in r12 if r["step"] == "step2" and r["dataset"] in ("M", "M_dense", "XL_grm", "R_real")]
    return {"speed": tsv(tq + "speed_wide.tsv"), "calib": tsv(tq + "calibration.tsv"), "power": tsv(tq + "power.tsv"),
            "perm": tsv(tq + "permutations.tsv"), "env": env, "history": hist, "cores": cores,
            "r12_scaling": [r for r in r12 if r["cores"] == 1 and r["step"] == "step2"]}

LOADERS = {"age": data_age, "perf": data_perf,
           "usage": lambda: json.load(open(os.path.join(SRC, "usage_data.json")))}   # from Rscript doc/src/usage_data.R


# ------------------------------------------------------------------ assembly
def toc(body):
    items = []
    for lvl, hid, text in re.findall(r'<h([23]) id="([^"]+)"[^>]*>(.*?)</h\1>', body, flags=re.S):
        text = re.sub(r"<[^>]+>", "", text).strip()
        pad = ' style="padding-left:22px"' if lvl == "3" else ""
        items.append(f'<li><a href="#{hid}"{pad}>{text}</a></li>')
    return '<nav class="toc" aria-label="On this page"><ol>' + "".join(items) + "</ol></nav>"

def nav(current):
    cur = ' aria-current="page"'
    links = [f'<a href="{p}.html"{cur if p == current else ""}>{t}</a>' for p, t in PAGES]
    links += [f'<a href="{h}">{t}</a>' for h, t in EXTERNAL]
    return f'<header class="topnav"><div class="in"><a class="brand" href="index.html">suez</a><nav class="links" aria-label="Documentation">{"".join(links)}</nav></div></header>'

def render(name, ex, css):
    src = open(os.path.join(SRC, "pages", name + ".html")).read()
    title = re.search(r"<title>(.*?)</title>", src).group(1)
    body = re.sub(r"<title>.*?</title>\s*", "", src, count=1)
    def sub_ex(m):
        if m.group(1) not in ex: sys.exit(f"{name}: no example '{m.group(1)}' in examples.R")
        return highlight_r(ex[m.group(1)])
    body = re.sub(r"\{\{EX:(\w+)\}\}", sub_ex, body)
    body = re.sub(r"\{\{DATA:(\w+)\}\}", lambda m: json.dumps(LOADERS[m.group(1)](), separators=(",", ":")), body)
    body = body.replace("{{TOC}}", toc(body)).replace("{{CHARTS_JS}}", open(os.path.join(SRC, "charts.js")).read())
    body = re.sub(r"\{\{JS:(\w+)\}\}", lambda m: open(os.path.join(SRC, m.group(1) + ".js")).read(), body)
    head = f"<title>{title}</title>\n{FONTS}\n<style>\n{css}</style>\n"
    return head, nav(name) + '\n<div class="tip" id="tip" hidden></div>\n<div class="page">' + body + '</div>\n'

def main():
    ex, css = examples(), open(os.path.join(SRC, "style.css")).read()
    art = sys.argv[sys.argv.index("--artifact") + 1] if "--artifact" in sys.argv else None
    if art: os.makedirs(art, exist_ok=True)
    for name, _ in PAGES:
        head, body = render(name, ex, css)
        full = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                '<meta name="viewport" content="width=device-width, initial-scale=1">\n' + head + "</head>\n<body>\n" + body + "</body>\n</html>\n")
        open(os.path.join(DOC, name + ".html"), "w").write(full)
        if art: open(os.path.join(art, name + ".html"), "w").write(head + body if name == "index" else full)
        print("wrote doc/" + name + ".html")
    if art:
        for f in ("suez.html", "suez_vs_tensorqtl.html"): shutil.copy(os.path.join(DOC, f), art)
        print("artifact files in", art)

if __name__ == "__main__":
    main()
