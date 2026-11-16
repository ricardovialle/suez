// Figures of the usage guide. U: doc/src/usage_data.json (made by doc/src/usage_data.R from the examples' results).
const secs = v => v < 10 ? v.toFixed(2) + " s" : v.toFixed(1) + " s";
const GV = ["--g0", "--g1", "--g2"];
const jit = (i, w) => (((i * 9301 + 49297) % 233280) / 233280 - 0.5) * w;   // deterministic jitter
const mean = a => a.reduce((x, y) => x + y, 0) / a.length;
const ext = (a, pad = 0.05) => { const lo = Math.min(...a), hi = Math.max(...a), d = (hi - lo) || 1; return [lo - d * pad, hi + d * pad]; };
function ticks(lo, hi, n = 5) { const step = Math.pow(10, Math.floor(Math.log10((hi - lo) / n))); const s = [1, 2, 5, 10].map(k => k * step).find(k => (hi - lo) / k <= n) || step * 10;
  const out = []; for (let t = Math.ceil(lo / s) * s; t <= hi + 1e-9; t += s) out.push(+t.toFixed(10)); return out; }
const fmtT = t => Math.abs(t) >= 100 ? t.toFixed(0) : +t.toPrecision(3);
function svgIn(id, W, H, label) { return el("svg", {viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": label, style: `max-width:${Math.round(W * 1.3)}px`}, document.getElementById(id)); }
const ptip = p => p < 1e-3 ? p.toExponential(1) : p.toPrecision(2);

// expression by genotype, one panel per group (condition / time point); dots coloured by ALT allele count
function genoPanels(id, expr, geno, group, levels, title) {
  const n = levels.length, W = Math.max(300, 210 * n + 50), H = 240, m = {l: 44, r: 10, t: 26, b: 34}, pw = (W - m.l - m.r) / n;
  const svg = svgIn(id, W, H, title), [y0, y1] = ext(expr), Y = scale(y0, y1, H - m.b, m.t);
  for (const t of ticks(y0, y1)) { el("line", {x1: m.l, x2: W - m.r, y1: Y(t), y2: Y(t), class: "gl"}, svg);
    el("text", {x: m.l - 6, y: Y(t) + 4, "text-anchor": "end"}, svg).textContent = fmtT(t); }
  el("text", {x: m.l - 6, y: m.t - 10}, svg).textContent = "expression";
  levels.forEach((lv, k) => {
    const x0 = m.l + pw * k, X = g => x0 + pw * (0.2 + 0.3 * g);
    el("text", {x: x0 + pw / 2, y: m.t - 10, "text-anchor": "middle", class: "lab"}, svg).textContent = lv;
    if (k) el("line", {x1: x0, x2: x0, y1: m.t, y2: H - m.b, class: "ax"}, svg);
    const means = [];
    for (let g = 0; g <= 2; g++) {
      const idx = expr.map((_, i) => i).filter(i => group[i] === lv && geno[i] === g);
      el("text", {x: X(g), y: H - m.b + 15, "text-anchor": "middle"}, svg).textContent = g;
      idx.forEach(i => { const c = el("circle", {cx: X(g) + jit(i, pw * 0.16), cy: Y(expr[i]), r: 3.2, fill: `var(${GV[g]})`, "fill-opacity": 0.85}, svg);
        hover(c, `${lv}, ${g} ALT allele${g === 1 ? "" : "s"}<br>expression ${expr[i].toFixed(2)}`); });
      if (idx.length) means.push([X(g), mean(idx.map(i => expr[i])), idx.length]);
    }
    el("path", {d: means.map((p, i) => (i ? "L" : "M") + p[0] + "," + Y(p[1])).join(""), fill: "none", stroke: "var(--fg)", "stroke-width": 2}, svg);
    means.forEach(p => { const c = el("rect", {x: p[0] - 9, y: Y(p[1]) - 1.5, width: 18, height: 3, fill: "var(--fg)"}, svg);
      hover(c, `${lv}: mean ${p[1].toFixed(2)} (n = ${p[2]})`); });
  });
  el("text", {x: W - m.r, y: H - 2, "text-anchor": "end"}, svg).textContent = "ALT alleles";
}

// -log10 p along the gene's cis window, eQTL and response tests
function regional(id, d) {
  const W = 520, H = 240, m = {l: 40, r: 10, t: 24, b: 34}, xs = d.pos.map(p => p / 1e6);
  const svg = svgIn(id, W, H, "Regional association"), [x0, x1] = ext(xs.concat([d.body[0] / 1e6, d.body[1] / 1e6])), ymax = Math.max(...d.lg, ...d.li) * 1.08;
  const X = scale(x0, x1, m.l, W - m.r), Y = scale(0, ymax, H - m.b, m.t);
  axes(svg, W, H, m, X, Y, ticks(x0, x1, 5), ticks(0, ymax, 5), {fmtx: t => t.toFixed(2), xlab: "position (Mb)", ylab: "−log10 p"});
  el("rect", {x: X(d.body[0] / 1e6), y: H - m.b - 6, width: Math.max(2, X(d.body[1] / 1e6) - X(d.body[0] / 1e6)), height: 6, fill: "var(--muted)"}, svg);
  [["lg", "--geno", "eQTL (p_geno)"], ["li", "--inter", "response (p_interact)"]].forEach(([k, v, n]) => d[k].forEach((y, i) => {
    const lead = d.snps[i] === d.lead;
    const c = el("circle", {cx: X(xs[i]), cy: Y(y), r: lead ? 6 : 3.5, fill: `var(${v})`, "fill-opacity": lead ? 1 : 0.7, stroke: lead ? "var(--fg)" : "none", "stroke-width": 1.5}, svg);
    hover(c, `<div class="t">${d.snps[i]}</div>${sw(v)}${n}: −log10 p ${y.toFixed(1)}`); }));
}

// generic x/y scatter of -log10 p: two runs of the same tests, highlighted planted SNPs
function pp(id, x, y, hl, xl, yl, hlName) {
  const W = 360, H = 300, m = {l: 40, r: 12, t: 22, b: 34}, mx = Math.max(...x, ...y) * 1.05;
  const svg = svgIn(id, W, H, `${yl} against ${xl}`), X = scale(0, mx, m.l, W - m.r), Y = scale(0, mx, H - m.b, m.t);
  axes(svg, W, H, m, X, Y, ticks(0, mx, 5), ticks(0, mx, 5), {xlab: xl, ylab: yl});
  el("line", {x1: X(0), y1: Y(0), x2: X(mx), y2: Y(mx), stroke: "var(--muted)", "stroke-dasharray": "3 3"}, svg);
  x.map((_, i) => i).sort((a, b) => hl[a] - hl[b]).forEach(i => { const c = el("circle", {cx: X(x[i]), cy: Y(y[i]), r: hl[i] ? 4.5 : 3, fill: hl[i] ? "var(--inter)" : "var(--other)", "fill-opacity": hl[i] ? 0.95 : 0.6}, svg);
    hover(c, `${hl[i] ? hlName : "other SNP"}<br>${xl}: ${x[i].toFixed(1)}<br>${yl}: ${y[i].toFixed(1)}`); });
}

// QQ plot, several series; observed values above 20 are drawn at 20
function qq(id, series) {
  const W = 360, H = 300, m = {l: 40, r: 12, t: 22, b: 34};
  const ex = Math.max(...series.flatMap(s => s.d.e)) * 1.05, oy = Math.min(Math.max(ex, ...series.flatMap(s => s.d.o)) * 1.05, 20);
  const svg = svgIn(id, W, H, "QQ plot"), X = scale(0, ex, m.l, W - m.r), Y = scale(0, oy, H - m.b, m.t);
  axes(svg, W, H, m, X, Y, ticks(0, ex, 4), ticks(0, oy, 5), {xlab: "expected −log10 p", ylab: "observed −log10 p"});
  const dm = Math.min(ex, oy);
  el("line", {x1: X(0), y1: Y(0), x2: X(dm), y2: Y(dm), stroke: "var(--muted)", "stroke-dasharray": "3 3"}, svg);
  series.forEach(s => s.d.e.forEach((e, i) => { const c = el("circle", {cx: X(e), cy: Y(Math.min(s.d.o[i], oy)), r: 2.6, fill: `var(${s.v})`, "fill-opacity": 0.8}, svg);
    if (i < 30) hover(c, `${sw(s.v)}${s.n}<br>expected ${e.toFixed(2)}, observed ${s.d.o[i].toFixed(2)}`); }));
}

// heatmap of a samples x samples covariance block
function heat(id, M, labels) {
  const n = M.length, c = 15, m = {l: 96, t: 8, r: 8, b: 8}, W = m.l + n * c + m.r, H = m.t + n * c + m.b;
  const svg = svgIn(id, W, H, "Sample covariance block"), mx = Math.max(...M.flat().map(Math.abs));
  M.forEach((row, i) => { el("text", {x: m.l - 6, y: m.t + i * c + c * 0.75, "text-anchor": "end", style: "font-size:9px"}, svg).textContent = labels[i];
    row.forEach((v, j) => { const r = el("rect", {x: m.l + j * c, y: m.t + i * c, width: c - 1, height: c - 1, rx: 2, fill: v >= 0 ? "var(--geno)" : "var(--inter)", "fill-opacity": Math.max(0.04, Math.abs(v) / mx)}, svg);
      hover(r, `${labels[i]} × ${labels[j]}<br>covariance ${v.toFixed(3)}`); }); });
}

// horizontal bars with labels
function hbars(id, items, xlab, fv) {
  const W = 420, rowH = 22, m = {l: 120, r: 50, t: 8, b: 30}, H = m.t + m.b + rowH * items.length, mx = Math.max(...items.map(d => d[1])) * 1.05;
  const svg = svgIn(id, W, H, xlab), X = scale(0, mx, m.l, W - m.r);
  axes(svg, W, H, m, X, () => 0, ticks(0, mx, 4), [], {xgrid: true, xlab});
  items.forEach(([lab, v, col], i) => { const y = m.t + i * rowH + 4, w = Math.max(1, X(v) - X(0));
    el("text", {x: m.l - 8, y: y + 11, "text-anchor": "end", class: "lab"}, svg).textContent = lab;
    el("path", {d: `M${X(0)},${y}H${X(0) + w - Math.min(4, w)}Q${X(0) + w},${y} ${X(0) + w},${y + Math.min(4, w)}V${y + 14 - Math.min(4, w)}Q${X(0) + w},${y + 14} ${X(0) + w - Math.min(4, w)},${y + 14}H${X(0)}Z`, fill: `var(${col || "--geno"})`}, svg);
    el("text", {x: X(0) + w + 5, y: y + 11, class: "lab"}, svg).textContent = fv ? fv(v) : v < 0.01 ? v.toExponential(1) : v.toFixed(3); });
}

// ---------------------------------------------------------------- the figures
(function () {
  const s = U.step1, f = s.factors;
  hbars("f-step1-var", [["noise", s.s.isotropic_noise], ["same individual", s.s.same_ind], ["kinship (GRM)", s.s.kinship_matrix]]
    .concat(f.slice(0, 6).map((v, i) => [`factor ${i + 1}`, v, "--g1"])), "variance");
  heat("f-step1-heat", s.cov, s.labels);
  const kept = f.filter(v => v > 0.01).length;
  document.getElementById("t-step1").innerHTML = `Two latent factors carry most of the shared variation (hidden confounders), ${kept} of ${f.length} kept a variance above 0.01${kept < f.length ? " and the rest shrank towards 0" : ""}. The same-individual and kinship terms are small here because the factors already explain much of the between-sample covariance.`;
})();
(function () {
  const d = U.two;
  genoPanels("f-two-geno", d.expr, d.geno, d.cond, ["control", "treated"], "Expression by genotype and condition");
  regional("f-two-reg", d);
  document.getElementById("t-two").innerHTML = `<b>${d.gene}</b>, lead SNP ${d.snp} (planted: ${d.planted || "none"}). The genotype effect is visible only in treated samples: p_interact = ${ptip(d.p_interact)}, p_geno = ${ptip(d.p_geno)}.`;
})();
(function () {
  const d = U.eqtl;
  genoPanels("f-eqtl", d.expr, d.geno, d.expr.map(() => "control samples"), ["control samples"], "Expression by genotype");
  document.getElementById("t-eqtl").innerHTML = `<b>${d.gene}</b> / ${d.snp} (planted: ${d.planted || "none"}): β = ${d.beta} ± ${d.se} per ALT allele, p_geno = ${ptip(d.p)}. The slope through the genotype means is what <code>beta_geno</code> estimates.`;
})();
(function () {
  const d = U.factor;
  genoPanels("f-factor", d.expr, d.geno, d.tp, ["0h", "6h", "24h"], "Expression by genotype at each time point");
  document.getElementById("t-factor").innerHTML = `<b>${d.gene}</b> / ${d.snp}: flat at 0h, an eQTL at 6h and 24h. One joint test with ${d.df} df: p_interact = ${ptip(d.p)}.`;
})();
(function () {
  const d = U.numeric, W = 520, H = 260, m = {l: 44, r: 70, t: 22, b: 34};
  const svg = svgIn("f-numeric", W, H, "Expression against log dose by genotype"), [x0, x1] = ext(d.x), [y0, y1] = ext(d.expr);
  const X = scale(x0, x1, m.l, W - m.r), Y = scale(y0, y1, H - m.b, m.t);
  axes(svg, W, H, m, X, Y, ticks(x0, x1, 5), ticks(y0, y1, 5), {fmtx: fmtT, fmty: fmtT, xlab: "log10(dose + 1)", ylab: "expression"});
  d.x.forEach((x, i) => { const c = el("circle", {cx: X(x), cy: Y(d.expr[i]), r: 3, fill: `var(${GV[d.geno[i]]})`, "fill-opacity": 0.75}, svg);
    hover(c, `${d.geno[i]} ALT alleles<br>log dose ${x.toFixed(2)}, expression ${d.expr[i].toFixed(2)}`); });
  for (let g = 0; g <= 2; g++) {   // least-squares line per genotype group
    const ii = d.x.map((_, i) => i).filter(i => d.geno[i] === g); if (ii.length < 3) continue;
    const mx = mean(ii.map(i => d.x[i])), my = mean(ii.map(i => d.expr[i]));
    const b = ii.reduce((a, i) => a + (d.x[i] - mx) * (d.expr[i] - my), 0) / ii.reduce((a, i) => a + (d.x[i] - mx) ** 2, 0);
    el("line", {x1: X(x0), x2: X(x1), y1: Y(my + b * (x0 - mx)), y2: Y(my + b * (x1 - mx)), stroke: `var(${GV[g]})`, "stroke-width": 2.5}, svg);
    el("text", {x: X(x1) + 6, y: Y(my + b * (x1 - mx)) + 4, class: "lab"}, svg).textContent = `${g} ALT`;
  }
  document.getElementById("t-numeric").innerHTML = `<b>${d.gene}</b> / ${d.snp}: the lines fan out as the dose rises; the response test asks whether their slopes differ (p_interact = ${ptip(d.p)}).`;
})();
(function () {
  const d = U.spline, W = 360, H = 250, m = {l: 40, r: 10, t: 22, b: 34};
  const svg = svgIn("f-spline-beta", W, H, "beta by age"), lo = d.beta.map((b, i) => b - 1.96 * d.se[i]), hi = d.beta.map((b, i) => b + 1.96 * d.se[i]);
  const [y0, y1] = ext(lo.concat(hi, d.beta_lin, d.truth)), X = scale(30, 85, m.l, W - m.r), Y = scale(y0, y1, H - m.b, m.t);
  axes(svg, W, H, m, X, Y, [30, 40, 50, 60, 70, 80], ticks(y0, y1, 5), {fmty: fmtT, xlab: "age", ylab: "β per allele"});
  el("line", {x1: X(30), x2: X(85), y1: Y(0), y2: Y(0), class: "ax"}, svg);
  el("path", {d: path(d.ages, hi, X, Y) + path(d.ages.slice().reverse(), lo.slice().reverse(), X, Y).replace("M", "L") + "Z", fill: "var(--band)"}, svg);
  el("path", {d: path(d.ages, d.truth, X, Y), fill: "none", stroke: "var(--truth)", "stroke-width": 2, "stroke-dasharray": "5 4"}, svg);
  el("path", {d: path(d.ages, d.beta_lin, X, Y), fill: "none", stroke: "var(--linear)", "stroke-width": 2}, svg);
  el("path", {d: path(d.ages, d.beta, X, Y), fill: "none", stroke: "var(--spline)", "stroke-width": 2}, svg);
  const hit = el("rect", {x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent"}, svg);
  hover(hit, ev => { const r = svg.getBoundingClientRect(), a = 30 + ((ev.clientX - r.left) / r.width * W - m.l) / (W - m.l - m.r) * 55, i = Math.max(0, Math.min(d.ages.length - 1, Math.round((a - 30) / 5)));
    return `<div class="t">Age ${d.ages[i]}</div>${sw("--truth")}true ${d.truth[i].toFixed(2)}<br>${sw("--spline")}spline ${d.beta[i].toFixed(2)} ± ${(1.96 * d.se[i]).toFixed(2)}<br>${sw("--linear")}linear ${d.beta_lin[i].toFixed(2)}`; });
  const svg2 = svgIn("f-spline-mu", W, H, "mu by age"), mlo = d.mu.map((b, i) => b - 1.96 * d.se_mu[i]), mhi = d.mu.map((b, i) => b + 1.96 * d.se_mu[i]);
  const [u0, u1] = ext(mlo.concat(mhi)), Y2 = scale(u0, u1, H - m.b, m.t);
  axes(svg2, W, H, m, X, Y2, [30, 40, 50, 60, 70, 80], ticks(u0, u1, 5), {fmty: fmtT, xlab: "age", ylab: "μ(age), relative"});
  el("path", {d: path(d.ages, mhi, X, Y2) + path(d.ages.slice().reverse(), mlo.slice().reverse(), X, Y2).replace("M", "L") + "Z", fill: "var(--band)"}, svg2);
  el("path", {d: path(d.ages, d.mu, X, Y2), fill: "none", stroke: "var(--spline)", "stroke-width": 2}, svg2);
  document.getElementById("t-spline").innerHTML = `<b>${d.gene}</b> / ${d.snp}, planted as a mid-life peak. The spline curve follows it (p_interact = ${ptip(d.p_spline)}, 3 df); a linear age term draws a declining line and misses the peak (p = ${ptip(d.p_linear)}, 1 df).`;
})();
(function () {
  const d = U.unbalanced;
  pp("f-unbal", d.full, d.sub, d.causal, "−log10 p, all samples", "−log10 p, 80% of samples", "planted response SNP");
  const pi = Object.entries(d.per_ind).map(([k, v]) => `${v} with ${k}`).join(", ");
  document.getElementById("t-unbal").innerHTML = `Response test on all ${d.n_full} samples against ${d.n_sub} kept samples (individuals: ${pi}). Null SNPs stay along the diagonal. Planted response SNPs mostly drop below it, the information lost with a fifth of the samples, but the strongest stay well clear of the nulls. No individual had to be dropped.`;
})();
(function () {
  const d = U.covariates;
  pp("f-cov", d.without, d.with, d.causal, "−log10 p_geno, no covariates", "−log10 p_geno, with batch and sex", "planted SNP");
  const bm = Object.entries(d.batch_means).map(([k, v]) => `${k} ${v.toFixed(2)}`).join(", ");
  document.getElementById("t-cov").innerHTML = `Mean expression differs by batch (${bm}). Adding batch as a covariate removes that noise: planted SNPs (orange) move above the diagonal, null SNPs stay on it.`;
})();
(function () {
  heat("f-rel-cov", U.related.cov, U.related.labels);
  heat("f-rel-grm", U.related.grm, U.related.labels);
})();
(function () {
  const d = U.la, W = 380, H = 330, m = {l: 44, r: 12, t: 22, b: 38};
  const COL = {"AFR eQTL": "--g1", "AFR response eQTL": "--g1", "EUR eQTL": "--inter", "shared eQTL": "--tq2", "none planted": "--other"};
  const all = d.b_afr.concat(d.b_eur), [lo, hi] = ext(all), svg = svgIn("f-la", W, H, "Ancestry-specific effects at lead SNPs");
  const X = scale(lo, hi, m.l, W - m.r), Y = scale(lo, hi, H - m.b, m.t);
  axes(svg, W, H, m, X, Y, ticks(lo, hi, 5), ticks(lo, hi, 5), {fmtx: fmtT, fmty: fmtT, xlab: "beta_AFR", ylab: "beta_EUR"});
  el("line", {x1: X(lo), y1: Y(lo), x2: X(hi), y2: Y(hi), stroke: "var(--muted)", "stroke-dasharray": "3 3"}, svg);
  d.b_afr.forEach((x, i) => { const t = d.type[i], ring = /response/.test(t);
    const c = el("circle", {cx: X(x), cy: Y(d.b_eur[i]), r: 5.5, fill: ring ? "var(--surface)" : `var(${COL[t]})`, stroke: `var(${COL[t]})`, "stroke-width": ring ? 2.5 : 0}, svg);
    hover(c, `<div class="t">${t}</div>beta_AFR ${x.toFixed(2)}<br>beta_EUR ${d.b_eur[i].toFixed(2)}`); });
  const types = Object.keys(d.found);
  hbars("f-la-found", types.flatMap(t => [[`${t.replace(" eQTL", "")}: local ancestry`, d.found[t], "--geno"], [`pooled genotype`, d.found_pooled[t], "--other"]]), "share of planted SNPs with p_geno < 10⁻³", v => Math.round(v * 100) + "%");
})();
(function () {
  const d = U.fdr, W = 520, H = 240, m = {l: 40, r: 12, t: 22, b: 34};
  const svg = svgIn("f-fdr", W, H, "Gene-level p-values by rank"), n = Math.min(80, d.inter.length), mx = Math.max(...d.geno.slice(0, n), ...d.inter.slice(0, n)) * 1.05;
  const X = scale(1, n, m.l, W - m.r), Y = scale(0, mx, H - m.b, m.t);
  axes(svg, W, H, m, X, Y, ticks(1, n, 6).filter(t => t >= 1), ticks(0, mx, 5), {xlab: "gene rank", ylab: "−log10 gene-level p"});
  [["geno", "--geno", "eQTL"], ["inter", "--inter", "response"]].forEach(([k, v, nm]) => d[k].slice(0, n).forEach((y, i) => {
    const sig = d[k + "_sig"][i], c = el("circle", {cx: X(i + 1), cy: Y(y), r: sig ? 4 : 3, fill: sig ? `var(${v})` : "var(--surface)", stroke: `var(${v})`, "stroke-width": 1.5}, svg);
    hover(c, `${nm} rank ${i + 1}<br>−log10 p ${y.toFixed(1)}<br>${sig ? "FDR < 5%" : "not significant"}${k === "inter" && d.planted_inter[i] ? "<br>planted response eQTL" : ""}`); }));
  const ns = k => d[k + "_sig"].filter(Boolean).length, pl = d.inter_sig.filter((s, i) => s && d.planted_inter[i]).length;
  document.getElementById("t-fdr").innerHTML = `Filled points pass 5% FDR: <b>${ns("geno")}</b> eGenes and <b>${ns("inter")}</b> response eGenes, ${pl} of which are planted response eQTLs (10 were planted).`;
})();
(function () {
  qq("f-calib", [{d: U.calib.geno, v: "--geno", n: "p_geno"}, {d: U.calib.inter, v: "--inter", n: "p_interact"}]);
  qq("f-perm", [{d: U.perm.geno, v: "--geno", n: "p_geno, permuted"}, {d: U.perm.inter, v: "--inter", n: "p_interact, permuted"}]);
  document.getElementById("t-calib").innerHTML = `Real data: λ = ${U.calib.geno.lambda} (eQTL) and ${U.calib.inter.lambda} (response). The bulk follows the diagonal and true effects rise above it at the top right.`;
  document.getElementById("t-perm").innerHTML = `Permuted genotypes: λ = ${U.perm.geno.lambda} and ${U.perm.inter.lambda}; with no true effects left, every point stays near the diagonal.`;
})();
(function () {
  const ds = [...new Set(P.cores.map(r => r.dataset))], cols = ["--suez", "--tq", "--tq2", "--g2"];
  const W = 420, H = 260, m = {l: 40, r: 16, t: 22, b: 34}, svg = svgIn("f-scale", W, H, "Speed-up against cores");
  const X = scale(1, 8, m.l, W - m.r, true), Y = scale(0, 8, H - m.b, m.t);
  axes(svg, W, H, m, X, Y, [1, 2, 4, 8], [0, 2, 4, 6, 8], {xlab: "cores (log)", ylab: "speed-up"});
  el("line", {x1: X(1), y1: Y(1), x2: X(8), y2: Y(8), stroke: "var(--muted)", "stroke-dasharray": "3 3"}, svg);
  ds.forEach((dname, k) => { const rows = P.cores.filter(r => r.dataset === dname).sort((a, b) => a.cores - b.cores), w1 = rows[0].wall;
    el("path", {d: rows.map((r, i) => (i ? "L" : "M") + X(r.cores) + "," + Y(w1 / r.wall)).join(""), fill: "none", stroke: `var(${cols[k]})`, "stroke-width": 2}, svg);
    rows.forEach(r => { const c = el("circle", {cx: X(r.cores), cy: Y(w1 / r.wall), r: 4, fill: `var(${cols[k]})`, stroke: "var(--surface)", "stroke-width": 1.5}, svg);
      hover(c, `<div class="t">${dname}, ${r.cores} cores</div>${secs(r.wall)} (${(w1 / r.wall).toFixed(1)}×)`); });
  });
  document.getElementById("l-scale").innerHTML = ds.map((d, k) => `<span><i class="sw" style="border-color:var(${cols[k]})"></i>${d}</span>`).join("");
})();
