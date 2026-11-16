// Shared chart helpers for the suez docs: SVG elements, axes on one scale, hover tooltips. Colors come from CSS tokens.
const NS = "http://www.w3.org/2000/svg";
function el(tag, attrs, parent) { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
const tip = document.getElementById("tip");
function showTip(ev, html) { tip.innerHTML = html; tip.hidden = false;
  const r = tip.getBoundingClientRect(); let x = ev.clientX + 14, y = ev.clientY + 14;
  if (x + r.width > innerWidth - 8) x = ev.clientX - r.width - 14; if (y + r.height > innerHeight - 8) y = ev.clientY - r.height - 14;
  tip.style.left = x + "px"; tip.style.top = y + "px"; }
function hideTip() { tip.hidden = true; }
function hover(node, html) { node.addEventListener("pointermove", ev => showTip(ev, typeof html === "function" ? html(ev) : html)); node.addEventListener("pointerleave", hideTip); }
const sw = v => `<i style="background:var(${v})"></i>`;
const pct = x => (100 * x).toFixed(x < 0.1 ? 1 : 0) + "%";
const fmt = (x, d = 2) => x.toFixed(d);
// linear or log10 scale: s(v) -> pixel
function scale(d0, d1, r0, r1, log) { const f = log ? Math.log10 : (v => v); const a = f(d0), b = f(d1); return v => r0 + (f(v) - a) / (b - a) * (r1 - r0); }
// axes: xs, ys are scales; xt, yt tick values; fmtx / fmty label formatters
function axes(svg, W, H, m, xs, ys, xt, yt, o = {}) {
  for (const t of yt) { el("line", {x1: m.l, x2: W - m.r, y1: ys(t), y2: ys(t), class: "gl"}, svg);
    el("text", {x: m.l - 6, y: ys(t) + 4, "text-anchor": "end"}, svg).textContent = o.fmty ? o.fmty(t) : t; }
  for (const t of xt) { if (o.xgrid) el("line", {x1: xs(t), x2: xs(t), y1: m.t, y2: H - m.b, class: "gl"}, svg);
    el("text", {x: xs(t), y: H - m.b + 15, "text-anchor": "middle"}, svg).textContent = o.fmtx ? o.fmtx(t) : t; }
  el("line", {x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: "ax"}, svg);
  if (o.xlab) el("text", {x: W - m.r, y: H - 2, "text-anchor": "end"}, svg).textContent = o.xlab;
  if (o.ylab) el("text", {x: m.l - 6, y: m.t - 8, "text-anchor": "start"}, svg).textContent = o.ylab;
}
const path = (xs, ys, X, Y) => xs.map((x, i) => (i ? "L" : "M") + X(x).toFixed(1) + "," + Y(ys[i]).toFixed(1)).join("");
