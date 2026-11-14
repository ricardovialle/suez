# Generates man/figures/suez_schematic.svg (the README figure). Run from the package root:
#   python3 data-raw/make_schematic.py
import math, random, re

W, H = 1000, 812
INK, MUTED, LINE = "#1f2933", "#52606d", "#9aa5b1"
S1, S1F, S1L = "#2563a8", "#eef4fb", "#bcd3ee"     # step 1: blue
S2, S2F, S2L = "#1f7a5a", "#edf7f2", "#b9dfcd"     # step 2: green
S3, S3F, S3L = "#7a4a9e", "#f5f0fa", "#d9c8e8"     # local-ancestry mode: purple
AFR, EUR = "#c9822b", "#4f9c9a"                    # ancestry colours in the haplotype sketch
INF = "#f5f6f8"
FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif"
out = []

def esc(s): return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def text(x, y, s, size=13, color=INK, weight="normal", anchor="start", italic=False):
    """s may contain _{sub}; subscripts are drawn smaller and lower."""
    parts = re.split(r"(_\{[^}]*\})", s)
    sub, body, low = size * 0.72, [], False
    for p in parts:   # dy is applied to the next character, so the text after a subscript carries the reset
        if p.startswith("_{"):
            body.append(f'<tspan dy="{size*0.3 if not low else 0:.1f}" font-size="{sub:.1f}">{esc(p[2:-1])}</tspan>')
            low = True
        elif p:
            body.append(f'<tspan dy="{-size*0.3:.1f}">{esc(p)}</tspan>' if low else esc(p)); low = False
    st = ' font-style="italic"' if italic else ""
    out.append(f'<text x="{x}" y="{y}" font-size="{size}" fill="{color}" font-weight="{weight}" '
               f'text-anchor="{anchor}"{st}>{"".join(body)}</text>')

def rect(x, y, w, h, fill, stroke="none", r=8, sw=1.2, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')

def arrow(pts, color=LINE, sw=1.8, label=None, lpos=None):
    d = "M" + " L".join(f"{x},{y}" for x, y in pts)
    mk = {LINE: "a-grey", S1: "a-s1", S2: "a-s2", S3: "a-s3"}[color]
    out.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{sw}" marker-end="url(#{mk})"/>')
    if label: text(*lpos, label, 12, color, "600")

out.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
           f'font-family="{FONT}" role="img" aria-labelledby="t d">')
out.append('<title id="t">suez workflow</title><desc id="d">Step 1 learns a sample covariance shared by '
           'all genes; step 2 tests every cis-SNP with nested linear mixed models that use it.</desc>')
out.append("<defs>" + "".join(
    f'<marker id="{i}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
    f'<path d="M0,0 L10,5 L0,10 z" fill="{c}"/></marker>' for i, c in [("a-grey", LINE), ("a-s1", S1), ("a-s2", S2), ("a-s3", S3)]) + "</defs>")
rect(0, 0, W, H, "#ffffff", r=14)

# ---------------- inputs ----------------
text(30, 34, "Inputs", 15, MUTED, "700")
cards = [  # y, h, title, lines
    (50, 92, "Expression", ["genes × samples (log2)", "e.g. TMM log2 CPM"]),
    (158, 92, "Sample annotation", ["individual (repeated samples)", "condition: factor or numeric"]),
    (266, 72, "Kinship / GRM", ["optional, individuals ×", "individuals"]),
    (420, 112, "Genotypes", ["SNPs × individuals (0/1/2)", "+ SNP positions", "+ gene positions"]),
    (548, 92, "Covariates", ["optional, samples × q:", "PCs, sex, batch, ..."]),
    (672, 112, "Local ancestry", ["optional, admixed samples:", "Tractor files per ancestry", "(dosage + hapcount)"]),
]
for y, h, t, ls in cards:
    rect(20, y, 200, h, INF, "#d9dde3")
    text(34, y + 24, t, 14, INK, "700")
    for i, l in enumerate(ls): text(34, y + 46 + 18 * i, l, 12, MUTED)

# ---------------- step 1 ----------------
X1, Y1, W1, H1 = 262, 22, 718, 290
rect(X1, Y1, W1, H1, S1F, S1L, r=12)
text(X1 + 18, Y1 + 30, "Step 1 · learn_covariance()", 16, S1, "700")
text(X1 + 18, Y1 + 50, "one sample covariance shared by all genes", 12.5, MUTED, italic=True)
fx, fy = X1 + 18, Y1 + 86
text(fx, fy, "for every gene g:   y_{g} = X b_{g} + e_{g},    e_{g} ~ N(0, Σ)", 14)
text(fx, fy + 22, "X = condition design (removed exactly: REML on contrasts ⟂ X)", 12, MUTED)
text(fx, fy + 54, "Σ = σ²_{noise} I + σ²_{ind} S + σ²_{kin} K + Σ_{p} σ²_{p} x_{p}x_{p}ᵀ", 15)
legend = [("S", "same individual (repeated samples)"), ("K", "kinship between individuals"),
          ("x_{p}", "P latent factors (hidden confounders)")]
for i, (a, b) in enumerate(legend):
    text(fx + 8, fy + 82 + 19 * i, a, 12.5, INK, "700"); text(fx + 44, fy + 82 + 19 * i, b, 12, MUTED)
text(fx, fy + 154, "Fitted jointly over all genes by maximum likelihood (Stan, L-BFGS);", 12, MUTED)
text(fx, fy + 170, "run time does not grow with the number of genes.", 12, MUTED)

# heatmap of cov_suez: 6 individuals x 2 conditions, ordered by individual
random.seed(3)
n = 12; ind = [i // 2 for i in range(n)]
kin = [[1.0 if a == b else (0.5 if {a, b} in ({0, 1}, {3, 4}) else 0.0) for b in range(6)] for a in range(6)]
f1 = [random.gauss(0, 1) for _ in range(n)]
f2 = [(1 if i % 2 else -1) * 0.8 + random.gauss(0, 0.4) for i in range(n)]
M = [[0.55 * (ind[i] == ind[j]) + 0.35 * kin[ind[i]][ind[j]] + 0.12 * f1[i] * f1[j] + 0.1 * f2[i] * f2[j]
      for j in range(n)] for i in range(n)]
mx = max(abs(v) for r_ in M for v in r_)
hx, hy, cs = X1 + 520, Y1 + 74, 13
def blend(v):
    t = max(-1, min(1, v / mx))
    c0 = (255, 255, 255); c1 = (37, 99, 168) if t >= 0 else (196, 92, 60)
    return "#%02x%02x%02x" % tuple(round(c0[k] + abs(t) * (c1[k] - c0[k])) for k in range(3))
for i in range(n):
    for j in range(n):
        out.append(f'<rect x="{hx + j*cs}" y="{hy + i*cs}" width="{cs}" height="{cs}" fill="{blend(M[i][j])}"/>')
out.append(f'<rect x="{hx}" y="{hy}" width="{n*cs}" height="{n*cs}" fill="none" stroke="{S1}" stroke-width="1"/>')
text(hx + n*cs/2, hy - 10, "samples × samples", 11.5, MUTED, anchor="middle")
text(hx + n*cs/2, hy + n*cs + 22, "cov_suez", 14, S1, "700", "middle")
text(hx + n*cs/2, hy + n*cs + 39, "= Σ without the noise term", 11.5, MUTED, anchor="middle")

# ---------------- step 2 ----------------
X2, Y2, W2, H2 = 262, 352, 718, 300
rect(X2, Y2, W2, H2, S2F, S2L, r=12)
text(X2 + 18, Y2 + 30, "Step 2 · map_interaction_qtl()", 16, S2, "700")
text(X2 + 18, Y2 + 50, "every gene × every SNP within ± cisdist of it", 12.5, MUTED, italic=True)

# column a: nested models + reaction-norm sketch
ax = X2 + 18
text(ax, Y2 + 82, "Three nested models", 13, INK, "700")
models = [("M0", "y = C b + e"), ("M1", "y = C b + g β + e"), ("M2", "y = C b + g β + (g·c) γ + e")]
for i, (m, f) in enumerate(models):
    text(ax, Y2 + 106 + 21 * i, m, 13, S2, "700"); text(ax + 32, Y2 + 106 + 21 * i, f, 13.5)
text(ax, Y2 + 176, "e ~ N(0, σ²_{k} K + σ² I),   K = cov_suez", 13)
text(ax, Y2 + 196, "C: condition + covariates, g: genotype (std.)", 11.5, MUTED)
# mini plot: expression vs genotype per condition
px, py, pw, ph = ax + 10, Y2 + 212, 150, 54
out.append(f'<path d="M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}" fill="none" stroke="{LINE}" stroke-width="1.2"/>')
for k, gt in enumerate(["0", "1", "2"]):
    gx = px + 25 + 50 * k
    text(gx, py + ph + 13, gt, 10.5, MUTED, anchor="middle")
    for (base, slope, col) in [(38, 2, "#7b8794"), (34, 24, S2)]:
        yv = py + ph - (base * 0.6 + slope * k * 0.9) * 0.8
        for dx in (-6, 0, 6):
            out.append(f'<circle cx="{gx+dx}" cy="{yv + random.uniform(-4, 4):.1f}" r="2.1" fill="{col}" opacity="0.75"/>')
for (base, slope, col) in [(38, 2, "#7b8794"), (34, 24, S2)]:
    y0 = py + ph - base * 0.48; y2 = py + ph - (base * 0.6 + slope * 2 * 0.9) * 0.8
    out.append(f'<line x1="{px+25}" y1="{y0:.1f}" x2="{px+125}" y2="{y2:.1f}" stroke="{col}" stroke-width="1.8"/>')
    text(px + pw - 14, y2 + 4, "treated" if col == S2 else "control", 11, col, "600")
text(px + pw / 2, py + ph + 26, "genotype", 10.5, MUTED, anchor="middle")

# column b: likelihood over h
bx = X2 + 300
text(bx, Y2 + 82, "Exact fit of each model", 13, INK, "700")
for i, l in enumerate(["rotate by the eigenvectors of K: the", "covariance becomes diagonal, b and", "the scale are closed-form, leaving",
                       "h = σ²_{k} / (σ²_{k} + σ²)", "(default K: eigenvectors built exactly)"]):
    text(bx, Y2 + 104 + 17 * i, l, 12, MUTED)
cx, cy, cw, ch = bx + 6, Y2 + 186, 200, 66
def ll(x): return 0.9 * math.exp(-((x + 1.6) / 0.8) ** 2) + 1.0 * math.exp(-((x - 1.2) / 0.6) ** 2)
pts = [(cx + (x + 4) / 8 * cw, cy + ch - ll(x) * ch * 0.9) for x in [i / 50 * 8 - 4 for i in range(51)]]
out.append(f'<path d="M{cx},{cy} L{cx},{cy+ch} L{cx+cw},{cy+ch}" fill="none" stroke="{LINE}" stroke-width="1.2"/>')
out.append('<path d="M' + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts) + f'" fill="none" stroke="{S2}" stroke-width="2"/>')
for x in [i / 16 * 8 - 4 for i in range(17)]:
    out.append(f'<circle cx="{cx + (x + 4) / 8 * cw:.1f}" cy="{cy + ch - ll(x) * ch * 0.9:.1f}" r="2.3" fill="#ffffff" stroke="{MUTED}" stroke-width="1"/>')
for xp in (-1.6, 1.2):
    X = cx + (xp + 4) / 8 * cw; Y = cy + ch - ll(xp) * ch * 0.9
    out.append(f'<path d="M{X:.1f},{Y-11:.1f} l-4,-7 h8 z" fill="{S2}"/>')
text(cx + cw / 2, cy + ch + 15, "log10 δ,  δ = σ²_{k} / σ²", 10.5, MUTED, anchor="middle")
text(cx - 6, cy + 8, "log L", 10.5, MUTED, anchor="end")
text(bx, Y2 + 290, "grid, then Newton at every local peak", 11.5, S2, "600")

# column c: tests
tx = X2 + 540
text(tx, Y2 + 82, "Likelihood-ratio tests", 13, INK, "700")
text(tx, Y2 + 108, "p_geno", 13.5, S2, "700"); text(tx, Y2 + 126, "M1 vs M0,  χ²_{1}", 12, MUTED)
text(tx, Y2 + 154, "p_interact", 13.5, S2, "700"); text(tx, Y2 + 172, "M2 vs M1,  χ²_{df_interact}", 12, MUTED)
text(tx, Y2 + 186, "(response eQTL)", 12, MUTED)
rect(tx - 6, Y2 + 202, 158, 74, "#ffffff", S2L, r=6, dash="4 3")
text(tx + 4, Y2 + 222, "optional bootstrap", 12, INK, "600")
for i, l in enumerate(["simulate under M1 to", "calibrate p_interact", "(l_boot_* columns)"]):
    text(tx + 4, Y2 + 239 + 15 * i, l, 11.5, MUTED)

# ---------------- local-ancestry mode ----------------
X3, Y3, W3, H3 = 262, 668, 718, 124
rect(X3, Y3, W3, H3, S3F, S3L, r=12)
text(X3 + 18, Y3 + 28, "Optional · local ancestry", 15, S3, "700")
text(X3 + 210, Y3 + 28, "map_interaction_qtl(local_ancestry = read_tractor(...))", 12, MUTED, italic=True)
# two haplotypes of one individual, coloured by ancestry, with the tested SNP
hx3, hy3, hw = X3 + 18, Y3 + 50, 196
for r_, segs in enumerate([[(0, 0.55, AFR), (0.55, 1, EUR)], [(0, 0.3, EUR), (0.3, 1, AFR)]]):
    yy = hy3 + 26 * r_
    for a, b, col in segs:
        out.append(f'<rect x="{hx3 + a*hw:.1f}" y="{yy}" width="{(b-a)*hw:.1f}" height="12" fill="{col}" opacity="0.85"/>')
sx = hx3 + 0.66 * hw
out.append(f'<line x1="{sx:.1f}" y1="{hy3-6}" x2="{sx:.1f}" y2="{hy3+44}" stroke="{INK}" stroke-width="1.4" stroke-dasharray="3 2"/>')
out.append(f'<circle cx="{sx:.1f}" cy="{hy3+6}" r="4" fill="#ffffff" stroke="{INK}" stroke-width="1.4"/>')
out.append(f'<circle cx="{sx:.1f}" cy="{hy3+32}" r="4" fill="{INK}"/>')
text(hx3, hy3 + 62, "AFR", 11, AFR, "700"); text(hx3 + 34, hy3 + 62, "EUR", 11, EUR, "700")
text(hx3 + 70, hy3 + 62, "● ALT allele at the SNP", 11, MUTED)
# terms and models
mx3 = X3 + 240
text(mx3, Y3 + 54, "genotype g split by haplotype ancestry:", 12, MUTED)
text(mx3, Y3 + 74, "LA = # AFR haplotypes,  G_{AFR}, G_{EUR} = ALT on each", 12.5)
text(mx3, Y3 + 96, "M1: LA + G_{AFR} + G_{EUR}", 13); text(mx3, Y3 + 114, "M2: + LA·c + G_{AFR}·c + G_{EUR}·c", 13)
# tests
tx3 = X3 + 545
text(tx3, Y3 + 54, "p_AFR, p_EUR", 12.5, S3, "700"); text(tx3 + 92, Y3 + 54, "(= Tractor)", 11.5, MUTED)
text(tx3, Y3 + 74, "p_het", 12.5, S3, "700"); text(tx3 + 44, Y3 + 74, "effects differ", 11.5, MUTED)
text(tx3, Y3 + 94, "p_interact_AFR / _EUR", 12.5, S3, "700")
text(tx3, Y3 + 112, "ancestry-specific response", 11.5, MUTED)

# ---------------- arrows ----------------
arrow([(220, 96), (X1 - 4, 96)])                      # expression -> step 1
arrow([(220, 204), (X1 - 4, 204)])                    # annotation -> step 1
arrow([(220, 302), (X1 - 4, 302)])                    # kinship -> step 1
arrow([(236, 96), (236, 470), (X2 - 4, 470)])         # expression -> step 2
arrow([(244, 204), (244, 486), (X2 - 4, 486)])        # annotation -> step 2
arrow([(220, 520), (X2 - 4, 520)])                    # genotypes -> step 2
arrow([(220, 594), (X2 - 4, 594)])                    # covariates -> step 2
arrow([(220, 728), (X3 - 4, 728)], S3)                # local ancestry -> local-ancestry mode
arrow([(X3 + W3 - 60, Y3 + 4), (X3 + W3 - 60, Y2 + H2 + 4)], S3)   # replaces g in step 2
arrow([(hx + n*cs/2, Y1 + H1 - 4), (hx + n*cs/2, Y2 + 4)], S1, 2.2, "sample_kernel", (hx + n*cs/2 + 10, Y1 + H1 + 26))

out.append("</svg>")
open("man/figures/suez_schematic.svg", "w").write("\n".join(out) + "\n")
