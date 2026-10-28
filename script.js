// ========================================================================
// === Constants & Helpers ================================================
// ========================================================================

function getGamma() {
  return parseFloat(document.getElementById("gamma").value) || 1.4;
}

function deg2rad(deg) { return (deg * Math.PI) / 180; }
function rad2deg(rad) { return (rad * 180) / Math.PI; }

function newtonRaphson(func, deriv, guess, target, g, maxIter = 50, tol = 1e-9) {
  let x = guess;
  for (let i = 0; i < maxIter; i++) {
    const f = func(x, g) - target;
    if (Math.abs(f) < tol) return x;
    const df = deriv(x, g);
    if (Math.abs(df) < 1e-15) break;
    x = x - f / df;
    if (x < 0) x = tol;
  }
  return NaN;
}

function updateText(id, value, precision = 4) {
  const el = document.getElementById(id);
  if (el) {
    const text = isNaN(value) ? "–" : value.toFixed(precision);
    if (el.tagName === "INPUT") {
      el.value = text;
    } else {
      el.textContent = text;
    }
  }
}

function setVisible(id, visible) {
  const el = document.getElementById(id);
  if (el) { el.style.display = visible ? "block" : "none"; }
}

function arcPath(cx, cy, r, a1, a2) {
  const x1 = cx + r * Math.cos(a1);
  const y1 = cy - r * Math.sin(a1);
  const x2 = cx + r * Math.cos(a2);
  const y2 = cy - r * Math.sin(a2);
  return `M ${x1},${y1} A ${r},${r} 0 0 0 ${x2},${y2}`;
}

// Helper for flashing new results
function flashResults(panelId) {
    const results = document.querySelectorAll(`#${panelId}-calculator .results .row`);
    results.forEach(el => {
        el.classList.remove('is-new-result-row');
        void el.offsetWidth; // Trigger reflow
        el.classList.add('is-new-result-row');
    });
}

// ========================================================================
// === Gas Dynamics Models (The "Math") ===================================
// ========================================================================

const isentropic = {
  T0_T: (M, g) => (1 + (g - 1) / 2 * M ** 2),
  p0_p: (M, g) => (1 + (g - 1) / 2 * M ** 2) ** (g / (g - 1)),
  rho0_rho: (M, g) => (1 + (g - 1) / 2 * M ** 2) ** (1 / (g - 1)),
  A_Astar: (M, g) => (1 / M) * (((1 + (g - 1) / 2 * M ** 2)) / ((g + 1) / 2)) ** ((g + 1) / (2 * (g - 1))),
  M_from_T: (T0_T, g) => Math.sqrt((T0_T - 1) * 2 / (g - 1)),
  M_from_p: (p0_p, g) => Math.sqrt(((p0_p) ** ((g - 1) / g) - 1) * 2 / (g - 1)),
  M_from_rho: (rho0_rho, g) => Math.sqrt(((rho0_rho) ** (g - 1) - 1) * 2 / (g - 1)),
  A_Astar_func: (M, g) => isentropic.A_Astar(M, g),
  A_Astar_deriv: (M, g) => {
    return isentropic.A_Astar(M, g) * (M ** 2 - 1) / (M * (1 + (g - 1) / 2 * M ** 2));
  },
  M_from_A: (A_Astar, g, supersonic = false) => {
    if (A_Astar < 1) return NaN;
    if (A_Astar === 1) return 1;
    // FIX: Subsonic guess was converging to supersonic root.
    // M_sub ≈ 1/A_Astar is a much better guess.
    const guess = supersonic ? (1 + 0.5 * A_Astar) : (1 / (A_Astar**0.5));
    return newtonRaphson(isentropic.A_Astar_func, isentropic.A_Astar_deriv, guess, A_Astar, g);
  }
};

function normalShockRatios(M1, g) {
  if (M1 <= 1) return {};
  const M1_sq = M1 ** 2;
  const gm1 = g - 1;
  const gp1 = g + 1;

  const M2 = Math.sqrt((1 + gm1 / 2 * M1_sq) / (g * M1_sq - gm1 / 2));
  const p2p1 = 1 + 2 * g / gp1 * (M1_sq - 1);
  const rho2rho1 = gp1 * M1_sq / (gm1 * M1_sq + 2);
  const T2T1 = p2p1 / rho2rho1;
  const p02p01 = isentropic.p0_p(M2, g) * (p2p1) * (1 / isentropic.p0_p(M1, g));
  const p02p1 = isentropic.p0_p(M2, g) * p2p1;

  return { M1, M2, p2p1, rho2rho1, T2T1, p02p01, p02p1 };
}

const normalShockInverse = {
  M1_from_M2: (M2, g) => {
    if (M2 <= 0) return NaN;
    const M2_sq = M2 ** 2;
    const gm1 = g - 1;
    const M1_sq = (1 + gm1 / 2 * M2_sq) / (g * M2_sq - gm1 / 2);
    return M1_sq < 1 ? NaN : Math.sqrt(M1_sq);
  },
  M1_from_p2p1: (p2p1, g) => {
    if (p2p1 <= 1) return NaN;
    const M1_sq = (p2p1 - 1) * (g + 1) / (2 * g) + 1;
    return Math.sqrt(M1_sq);
  },
  M1_from_rho2rho1: (rho2rho1, g) => {
    const gp1 = g + 1;
    const gm1 = g - 1;
    if (rho2rho1 <= 1 || rho2rho1 > gp1 / gm1) return NaN;
    const M1_sq = (2 * rho2rho1) / (gp1 - gm1 * rho2rho1);
    return Math.sqrt(M1_sq);
  },
  M1_from_T2T1: (T2T1, g) => {
    if (T2T1 <= 1) return NaN;
    const gp1 = g + 1;
    const gm1 = g - 1;
    const a = 2 * g * gm1;
    const b = (4 * g - gm1**2 - gp1**2 * T2T1);
    const c = -2 * gm1;
    const discriminant = b**2 - 4*a*c;
    if (discriminant < 0) return NaN;
    const X1 = (-b + Math.sqrt(discriminant)) / (2*a);
    return X1 < 1 ? NaN : Math.sqrt(X1);
  },
  p02p01_func: (M1, g) => normalShockRatios(M1, g).p02p01,
  p02p01_deriv: (M1, g) => {
    const h = 1e-6;
    const f_plus = normalShockRatios(M1 + h, g).p02p01 || 0;
    const f_minus = normalShockRatios(M1 - h, g).p02p01 || 0;
    return (f_plus - f_minus) / (2 * h);
  },
  M1_from_p02p01: (p02p01, g) => {
    if (p02p01 >= 1 || p02p01 <= 0) return NaN;
    const guess = 1 - 2 * (p02p01 - 1);
    return newtonRaphson(normalShockInverse.p02p01_func, normalShockInverse.p02p01_deriv, guess < 1 ? 1.1 : guess, p02p01, g);
  },
  p02p1_func: (M1, g) => normalShockRatios(M1, g).p02p1,
  p02p1_deriv: (M1, g) => {
    const h = 1e-6;
    const f_plus = normalShockRatios(M1 + h, g).p02p1 || 0;
    const f_minus = normalShockRatios(M1 - h, g).p02p1 || 0;
    return (f_plus - f_minus) / (2 * h);
  },
  M1_from_p02p1: (p02p1, g) => {
    if (p02p1 <= 1) return NaN;
    const guess = 1 + (p02p1 / 5);
    return newtonRaphson(normalShockInverse.p02p1_func, normalShockInverse.p02p1_deriv, guess < 1 ? 1.1 : guess, p02p1, g);
  }
};

function obliqueShockRatios(M1, beta, theta, g) {
  const Mn1 = M1 * Math.sin(beta);
  if (Mn1 <= 1) return { M1, theta: rad2deg(theta), beta: rad2deg(beta) };
  
  const { M2: Mn2, p2p1, rho2rho1, T2T1, p02p01 } = normalShockRatios(Mn1, g);
  if (isNaN(Mn2)) return { M1, theta: rad2deg(theta), beta: rad2deg(beta) };
  
  const M2 = Mn2 / Math.sin(beta - theta);
  return { M1, M2, Mn1, Mn2, p2p1, rho2rho1, T2T1, p02p01, theta: rad2deg(theta), beta: rad2deg(beta) };
}

function ob_beta_from_M1_theta(M1, theta, g) {
  const thetaRad = deg2rad(theta);
  const mu = Math.asin(1 / M1);
  const c = Math.tan(mu) ** 2;
  const a = ((g - 1) / 2 + ((g + 1) / 2) * c) * Math.tan(thetaRad);
  const b = ((g + 1) / 2 + ((g + 3) / 2) * c) * Math.tan(thetaRad);
  const num = 4 * (1 - 3 * a * b) ** 3;
  const den = (27 * a * a * c + 9 * a * b - 2) ** 2;
  const d = Math.sqrt(num / den - 1);
  const term1 = (b + 9 * a * c) / (2 * (1 - 3 * a * b));
  const term2 = (d * (27 * a * a * c + 9 * a * b - 2)) / (6 * a * (1 - 3 * a * b));
  const angle = (1 / 3) * Math.atan(1 / d); // n=0 for weak shock
  return Math.atan(term1 - term2 * Math.tan(angle));
}

function ob_theta_from_M1_beta(M1, beta, g) {
  const betaRad = deg2rad(beta);
  const M1_sq = M1 ** 2;
  const num = (M1_sq * Math.sin(betaRad) ** 2 - 1);
  const den = M1_sq * (g + Math.cos(2 * betaRad)) + 2;
  const thetaRad = Math.atan(2 * (1/Math.tan(betaRad)) * num / den);
  return rad2deg(thetaRad);
}

function ob_M1_from_theta_beta(theta, beta, g) {
  const thetaRad = deg2rad(theta);
  const betaRad = deg2rad(beta);
  if (beta <= theta) return NaN;
  const T = Math.tan(thetaRad);
  const C = 1 / Math.tan(betaRad);
  const S2 = Math.sin(betaRad) ** 2;
  const G = g + Math.cos(2 * betaRad);
  const M1_sq = 2 * (C + T) / (2 * C * S2 - T * G);
  return M1_sq < 1 ? NaN : Math.sqrt(M1_sq);
}

const prandtlMeyer = {
  nu: (M, g) => {
    if (M < 1) return 0;
    if (M === Infinity) {
        const gp1 = g + 1;
        const gm1 = g - 1;
        return rad2deg(Math.sqrt(gp1 / gm1) * Math.PI / 2 - Math.PI / 2);
    }
    const gp1 = g + 1;
    const gm1 = g - 1;
    const g_ratio = Math.sqrt(gp1 / gm1);
    const M_sq = M ** 2;
    const term1 = Math.sqrt(gm1 / gp1 * (M_sq - 1));
    const term2 = Math.sqrt(M_sq - 1);
    return rad2deg(g_ratio * Math.atan(term1) - Math.atan(term2));
  },
  nu_rad: (M, g) => deg2rad(prandtlMeyer.nu(M, g)),
  nu_deriv: (M, g) => {
    if (M <= 1) return 0;
    const M_sq = M ** 2;
    const gm1 = g - 1;
    return (Math.sqrt(M_sq - 1) / (1 + gm1 / 2 * M_sq)) * (1 / M);
  },
  M_from_nu: (nu_deg, g) => {
    if (nu_deg < 0) return 1;
    const gp1 = g + 1;
    const gm1 = g - 1;
    const nu_max = rad2deg(Math.sqrt(gp1 / gm1) * Math.PI / 2 - Math.PI / 2);
    if (nu_deg > nu_max) return Infinity;
    const nu_rad_target = deg2rad(nu_deg);
    const guess = 1.5 + (nu_deg / nu_max) * 10;
    return newtonRaphson(prandtlMeyer.nu_rad, prandtlMeyer.nu_deriv, guess, nu_rad_target, g);
  }
};


// ========================================================================
// === Diagram Updaters ===================================================
// ========================================================================

// Figures share the visual language of supersonic-airfoils.html:
// shocks red, expansion fans blue, Mach waves grey dashed, bodies in the
// body color, region data in small square label boxes, all drawn to scale.

const FIG = {
  n(v, d) { return isFinite(v) ? v.toFixed(d).replace("-", "−") : "–"; },
  defs(id) {
    const mk = (name, col, size = 7) =>
      `<marker id="${id}-${name}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="${size}" markerHeight="${size}" orient="auto-start-reverse"><path d="M0,0L10,5L0,10z" fill="${col}"/></marker>`;
    return `<defs>${mk("m", "var(--text-muted)")}${mk("k", "var(--text-main)")}${mk("c", "var(--danger)", 6)}${mk("e", "var(--exp)", 6)}${mk("a", "var(--accent)", 6)}</defs>`;
  },
  line(x1, y1, x2, y2, stroke, w = 1.5, extra = "") {
    return `<line x1="${x1.toFixed(1)}" y1="${y1.toFixed(1)}" x2="${x2.toFixed(1)}" y2="${y2.toFixed(1)}" stroke="${stroke}" stroke-width="${w}" ${extra}/>`;
  },
  text(x, y, str, o = {}) {
    const a = o.anchor || "middle", sz = o.size || 13, fill = o.fill || "var(--text-main)", wt = o.weight || 400;
    const fam = o.num ? "var(--f-num)" : "var(--f-body)";
    return `<text x="${x.toFixed(1)}" y="${y.toFixed(1)}" style="font-family:${fam};font-size:${sz}px;font-weight:${wt};text-anchor:${a};dominant-baseline:middle" fill="${fill}">${str}</text>`;
  },
  // Square label box: first line bold, following lines muted.
  box(x, y, lines, anchor = "middle") {
    const w = Math.max(...lines.map(l => l.replace(/<[^>]+>/g, "").length)) * 6.9 + 16, h = lines.length * 16 + 8;
    const x0 = anchor === "start" ? x : anchor === "end" ? x - w : x - w / 2, y0 = y - h / 2;
    let s = `<rect x="${x0.toFixed(1)}" y="${y0.toFixed(1)}" width="${w.toFixed(1)}" height="${h}" fill="var(--card)" fill-opacity=".92" stroke="var(--border)"/>`;
    lines.forEach((l, i) => {
      s += FIG.text(x0 + 8, y0 + 12 + i * 16, l, { anchor: "start", size: 12.5, num: true, weight: i === 0 ? 600 : 400, fill: i === 0 ? "var(--text-main)" : "var(--text-muted)" });
    });
    return s;
  },
  // Arc centred at (cx,cy) from angle a1 to a2 (radians, counter-clockwise positive, y up)
  arc(cx, cy, r, a1, a2, stroke) {
    const p = a => [cx + r * Math.cos(a), cy - r * Math.sin(a)];
    const [x1, y1] = p(a1), [x2, y2] = p(a2);
    const sweep = a2 > a1 ? 0 : 1;
    return `<path d="M${x1.toFixed(1)},${y1.toFixed(1)} A${r},${r} 0 0 ${sweep} ${x2.toFixed(1)},${y2.toFixed(1)}" fill="none" stroke="${stroke}" stroke-width="1.4"/>`;
  },
  message(svg, msg) {
    const vb = svg.viewBox.baseVal;
    svg.innerHTML = FIG.text(vb.width / 2, vb.height / 2, msg, { fill: "var(--text-muted)", size: 14 });
  }
};

// --- Normal shock: velocity arrows scaled by u2/u1 = rho1/rho2 ---
function updateNormalShockDiagram(M1, r) {
  const svg = document.getElementById("ns-diagram");
  if (!svg) return;
  if (!isFinite(M1) || !r) { FIG.message(svg, "Enter a valid input and press Compute."); return; }
  const id = "ns", xs = 330, top = 40, bot = 200;
  const u21 = 1 / r.rho2rho1, L1 = 230, L2 = Math.max(L1 * u21, 14);
  let s = FIG.defs(id);
  // downstream region tint
  s += `<rect x="${xs}" y="${top}" width="${640 - xs - 20}" height="${bot - top}" fill="var(--danger)" fill-opacity=".05"/>`;
  s += FIG.text(30, 22, "1 · upstream", { anchor: "start", weight: 600 });
  s += FIG.text(xs + 16, 22, "2 · downstream", { anchor: "start", weight: 600 });
  // velocity arrows
  [70, 120, 170].forEach(y => {
    s += FIG.line(xs - 20 - L1, y, xs - 20, y, "var(--text-main)", 2, `marker-end="url(#${id}-k)"`);
    s += FIG.line(xs + 16, y, xs + 16 + L2, y, "var(--text-main)", 2, `marker-end="url(#${id}-k)"`);
  });
  // shock
  s += FIG.line(xs, top - 6, xs, bot + 6, "var(--danger)", 4);
  s += FIG.text(xs - 8, bot + 22, "normal shock", { anchor: "end", fill: "var(--danger)", size: 12.5, weight: 600 });
  // data
  s += FIG.box(30, 272, [`M₁ = ${FIG.n(M1, 3)}`, `u₁ (reference)`, `p₁, T₁, ρ₁, p₀₁`], "start");
  s += FIG.box(xs + 16, 272, [`M₂ = ${FIG.n(r.M2, 3)}`, `u₂/u₁ = ${FIG.n(u21, 3)}`, `p₂/p₁ = ${FIG.n(r.p2p1, 3)}  ·  T₂/T₁ = ${FIG.n(r.T2T1, 3)}`, `ρ₂/ρ₁ = ${FIG.n(r.rho2rho1, 3)}  ·  p₀₂/p₀₁ = ${FIG.n(r.p02p01, 4)}`], "start");
  svg.setAttribute("viewBox", "0 0 640 318");
  svg.innerHTML = s;
}

// --- Oblique shock on a compression corner, to scale ---
function updateObliqueDiagram(thetaDeg, betaDeg, r) {
  const svg = document.getElementById("ob-shockDiagram");
  if (!svg) return;
  const okT = isFinite(thetaDeg) && thetaDeg >= 0, okB = isFinite(betaDeg) && betaDeg > 0;
  if (!okT && !okB) { FIG.message(svg, "Enter two of M₁, θ, β and press Compute."); return; }
  const id = "ob", W = 640, H = 380, cx = 230, cy = 300, x0 = 20;
  const th = okT ? deg2rad(thetaDeg) : 0, be = okB ? deg2rad(betaDeg) : NaN;
  const P = (r_, a) => [cx + r_ * Math.cos(a), cy - r_ * Math.sin(a)];
  let s = FIG.defs(id);
  // body: flat wall then ramp at theta
  const rampLen = 460, re = P(rampLen, th);
  s += `<polygon points="${x0},${cy} ${cx},${cy} ${re[0].toFixed(1)},${re[1].toFixed(1)} ${re[0].toFixed(1)},${H} ${x0},${H}" fill="var(--body)" fill-opacity=".12" stroke="none"/>`;
  s += `<polyline points="${x0},${cy} ${cx},${cy} ${re[0].toFixed(1)},${re[1].toFixed(1)}" fill="none" stroke="var(--body)" stroke-width="3" stroke-linejoin="round"/>`;
  s += FIG.line(cx, cy, cx + 170, cy, "var(--text-muted)", 1, 'stroke-dasharray="4 4"');
  // streamlines: horizontal upstream, turned by theta behind the shock
  if (okB) {
    [45, 95, 145, 195].forEach(h => {
      const yi = cy - h, xi = cx + h / Math.tan(be);
      const xEnd = Math.min(xi, W - 10);
      const xm = x0 + 0.55 * (Math.min(xEnd, cx + 30) - x0);
      s += FIG.line(x0, yi, xm, yi, "var(--text-muted)", 1.3, `marker-end="url(#${id}-m)"`);
      s += FIG.line(xm, yi, xEnd, yi, "var(--text-muted)", 1.3);
      if (xi < W - 10) {
        const L = 400, ex = xi + L * Math.cos(th), ey = yi - L * Math.sin(th);
        s += FIG.line(xi, yi, ex, ey, "var(--text-muted)", 1.3);
      }
    });
    // shock
    const Ls = Math.min((cy - 14) / Math.sin(be), (W - 10 - cx) / Math.cos(be));
    const se = P(Ls, be);
    s += FIG.line(cx, cy, se[0], se[1], "var(--danger)", 3);
    // velocity decomposition at a point on the shock: V1 = tangential + normal
    if (r && isFinite(r.Mn1)) {
      const q = P(0.5 * Ls, be), V = 90;
      const tail = [q[0] - V, q[1]];
      const tEnd = [tail[0] + V * Math.cos(be) * Math.cos(be), tail[1] - V * Math.cos(be) * Math.sin(be)];
      s += FIG.line(tail[0], tail[1], q[0], q[1], "var(--text-main)", 2, `marker-end="url(#${id}-k)"`);
      s += FIG.line(tail[0], tail[1], tEnd[0], tEnd[1], "var(--accent)", 1.6, `stroke-dasharray="4 3" marker-end="url(#${id}-a)"`);
      s += FIG.line(tEnd[0], tEnd[1], q[0], q[1], "var(--accent)", 2, `marker-end="url(#${id}-a)"`);
      s += FIG.text(tail[0] - 6, tail[1], "V₁", { anchor: "end", size: 13, weight: 600 });
      s += FIG.text(tEnd[0] - 6, tEnd[1] - 12, `M₁ₙ = M₁ sin β = ${FIG.n(r.Mn1, 3)}`, { anchor: "end", size: 12.5, fill: "var(--accent)", num: true, weight: 600 });
    }
    // beta arc and label
    s += FIG.arc(cx, cy, 120, 0, be, "var(--danger)");
    const bl = P(134, be * 0.55);
    s += FIG.text(bl[0], bl[1], `β = ${FIG.n(betaDeg, 2)}°`, { anchor: "start", fill: "var(--danger)", num: true, weight: 600 });
  }
  if (okT && th > 0) {
    s += FIG.arc(cx, cy, 70, 0, th, "var(--text-main)");
    const tl = P(80, -0.12);
    s += FIG.text(tl[0], tl[1] + 4, `θ = ${FIG.n(thetaDeg, 2)}°`, { anchor: "start", num: true, weight: 600 });
  }
  // region labels
  s += FIG.box(x0 + 4, 30, r ? ["1 · upstream", `M₁ = ${FIG.n(r.M1, 3)}`] : ["1 · upstream"], "start");
  if (r && okB) {
    // region 2 is a thin sliver between shock and ramp: label it from inside the body with a leader
    const Ls2 = Math.min((cy - 14) / Math.sin(be), (W - 10 - cx) / Math.cos(be));
    const pt = P(Math.min(380, 0.85 * Ls2), (th + be) / 2), bx = W - 14, by = H - 34;
    s += FIG.line(pt[0], pt[1], bx - 70, by - 26, "var(--text-muted)", 1);
    s += `<circle cx="${pt[0].toFixed(1)}" cy="${pt[1].toFixed(1)}" r="3" fill="var(--text-main)"/>`;
    s += FIG.box(bx, by, ["2 · behind the shock", `M₂ = ${FIG.n(r.M2, 3)}`, `p₂/p₁ = ${FIG.n(r.p2p1, 3)}`], "end");
  }
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  svg.innerHTML = s;
}

// --- Wedge airfoil: freestream horizontal, wedge pitched nose-up by alpha ---
function updateWedgeDiagram(g, M1, delta, alpha, lower, upper) {
  const svg = document.getElementById("wa-diagram");
  if (!svg) return;
  const id = "wa", W = 640, H = 400, ox = 210, oy = 200, ch = 260;
  const d = deg2rad(delta), a = deg2rad(alpha), h = Math.tan(d);
  const toS = (x, y) => [ox + ch * (x * Math.cos(a) + y * Math.sin(a)), oy - ch * (-x * Math.sin(a) + y * Math.cos(a))];
  const ray = (p, ang, L) => [p[0] + L * Math.cos(ang), p[1] - L * Math.sin(ang)];
  const LE = toS(0, 0), UT = toS(1, h), LT = toS(1, -h), L = 300;
  let s = FIG.defs(id);
  // freestream
  [60, 130, 270, 340].forEach(y => { s += FIG.line(16, y, 96, y, "var(--text-muted)", 1.5, `marker-end="url(#${id}-m)"`); });
  s += FIG.text(16, 98, `M₁ = ${FIG.n(M1, 2)}`, { anchor: "start", num: true, fill: "var(--text-muted)" });
  // waves from the leading edge
  const fan = (a1, a2) => { let o = ""; for (let k = 0; k < 6; k++) { const t = k / 5, e = ray(LE, a1 + (a2 - a1) * t, L * 0.9); o += FIG.line(LE[0], LE[1], e[0], e[1], "var(--exp)", k === 0 || k === 5 ? 1.6 : 0.9); } return o; };
  if (lower && isFinite(lower.beta)) { const e = ray(LE, -deg2rad(lower.beta), L); s += FIG.line(LE[0], LE[1], e[0], e[1], "var(--danger)", 2.6); }
  if (upper && upper.type === "Shock") { const e = ray(LE, deg2rad(upper.beta), L); s += FIG.line(LE[0], LE[1], e[0], e[1], "var(--danger)", 2.6); }
  else if (upper && upper.type === "Expansion") { s += fan(Math.asin(1 / M1), (d - a) + Math.asin(1 / upper.M2)); }
  else if (upper) { const e = ray(LE, Math.asin(1 / M1), L); s += FIG.line(LE[0], LE[1], e[0], e[1], "var(--text-muted)", 1.2, 'stroke-dasharray="6 5"'); }
  // chord line and alpha
  s += FIG.line(LE[0], LE[1], LE[0] + 130, LE[1], "var(--text-muted)", 1, 'stroke-dasharray="3 4"');
  s += FIG.line(LE[0], LE[1], (UT[0] + LT[0]) / 2, (UT[1] + LT[1]) / 2, "var(--text-muted)", 0.8, 'stroke-dasharray="2 3"');
  if (Math.abs(alpha) > 0.3) {
    const fwd = ray(LE, Math.PI - a, 90), hor = [LE[0] - 90, LE[1]];
    s += FIG.line(LE[0], LE[1], fwd[0], fwd[1], "var(--accent)", 1, 'stroke-dasharray="3 4"');
    s += FIG.line(LE[0], LE[1], hor[0], hor[1], "var(--text-muted)", 1, 'stroke-dasharray="3 4"');
    s += FIG.arc(LE[0], LE[1], 70, Math.PI, Math.PI - a, "var(--accent)");
    const al = ray(LE, Math.PI - a / 2, 78);
    s += FIG.text(al[0] - 12, al[1] + (a > 0 ? -12 : 12), `α = ${FIG.n(alpha, 1)}°`, { anchor: "end", fill: "var(--accent)", num: true, size: 12.5, weight: 600 });
  }
  // body
  s += `<polygon points="${[LE, UT, LT].map(p => p[0].toFixed(1) + "," + p[1].toFixed(1)).join(" ")}" fill="var(--body)" fill-opacity=".14" stroke="var(--body)" stroke-width="2" stroke-linejoin="round"/>`;
  // pressure arrows into each face, length proportional to p/p1
  if (lower && upper && isFinite(lower.p_p1) && isFinite(upper.p_p1)) {
    const pmax = Math.max(lower.p_p1, upper.p_p1, 1), k = 55 / pmax;
    const face = (A, B, p, up, label) => {
      const dx = B[0] - A[0], dy = B[1] - A[1], l = Math.hypot(dx, dy);
      let n = [dy / l, -dx / l]; if ((up && n[1] > 0) || (!up && n[1] < 0)) n = [-n[0], -n[1]]; // outward, screen coords
      const len = Math.max(p * k, 14), col = p > 1.0005 ? "var(--danger)" : p < 0.9995 ? "var(--exp)" : "var(--text-muted)", mk = p > 1.0005 ? "c" : p < 0.9995 ? "e" : "m";
      let o = "";
      [0.3, 0.55, 0.8].forEach(t => {
        const b = [A[0] + dx * t, A[1] + dy * t];
        o += FIG.line(b[0] + n[0] * (len + 3), b[1] + n[1] * (len + 3), b[0] + n[0] * 3, b[1] + n[1] * 3, col, 1.6, `marker-end="url(#${id}-${mk})"`);
      });
      const m = [A[0] + dx * 0.55, A[1] + dy * 0.55], dd = len + 34;
      o += FIG.box(m[0] + n[0] * dd, m[1] + n[1] * dd, label);
      return o;
    };
    s += face(LE, UT, upper.p_p1, true, [`upper · p/p₁ = ${FIG.n(upper.p_p1, 3)}`, `M = ${FIG.n(upper.M2, 3)}`]);
    s += face(LE, LT, lower.p_p1, false, [`lower · p/p₁ = ${FIG.n(lower.p_p1, 3)}`, `M = ${FIG.n(lower.M2, 3)}`]);
  }
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  svg.innerHTML = s;
}

// ========================================================================
// === Calculator Controllers =============================================
// ========================================================================

function clearResults(prefix, fields) {
  fields.forEach(id => updateText(`${prefix}-${id}`, NaN));
}

// --- 1. Isentropic Controller ---
function calculateIsentropic() {
  const g = getGamma();
  const inputType = document.getElementById("is-inputType").value;
  const inputValue = parseFloat(document.getElementById("is-inputValue").value);
  const warning = document.getElementById("is-warning");
  warning.textContent = "";
  let M = NaN, T0_T, p0_p, rho0_rho, A_Astar;
  const resultIDs = ["is-M", "is-p0_p", "is-rho0_rho", "is-T0_T", "is-A_Astar"];

  try {
    switch (inputType) {
      case "M":
        M = inputValue;
        if (M <= 0) throw new Error("M must be > 0.");
        break;
      case "T0_T":
        T0_T = inputValue;
        if (T0_T < 1) throw new Error("T₀/T must be ≥ 1.");
        M = isentropic.M_from_T(T0_T, g);
        break;
      case "p0_p":
        p0_p = inputValue;
        if (p0_p < 1) throw new Error("p₀/p must be ≥ 1.");
        M = isentropic.M_from_p(p0_p, g);
        break;
      case "rho0_rho":
        rho0_rho = inputValue;
        if (rho0_rho < 1) throw new Error("ρ₀/ρ must be ≥ 1.");
        M = isentropic.M_from_rho(rho0_rho, g);
        break;
      case "A_Astar_sub":
      case "A_Astar_super":
        A_Astar = inputValue;
        if (A_Astar < 1) throw new Error("A/A* must be ≥ 1.");
        M = isentropic.M_from_A(A_Astar, g, inputType === "A_Astar_super");
        break;
    }
    if (isNaN(M)) throw new Error("Calculation failed. Check inputs or solver failed.");
    T0_T = isentropic.T0_T(M, g);
    p0_p = isentropic.p0_p(M, g);
    rho0_rho = isentropic.rho0_rho(M, g);
    A_Astar = isentropic.A_Astar(M, g);
    flashResults("isentropic");
  } catch (err) {
    warning.textContent = err.message;
    clearResults("is", resultIDs.map(id => id.split('-')[1]));
  }
  updateText("is-M", M);
  updateText("is-p0_p", p0_p);
  updateText("is-rho0_rho", rho0_rho);
  updateText("is-T0_T", T0_T);
  updateText("is-A_Astar", A_Astar);
}

// --- 2. Normal Shock Controller ---
function calculateNormalShock() {
  const g = getGamma();
  const inputType = document.getElementById("ns-inputType").value;
  const inputValue = parseFloat(document.getElementById("ns-inputValue").value);
  const warning = document.getElementById("ns-warning");
  warning.textContent = "";
  let M1 = NaN;
  const resultIDs = ["M1", "M2", "p2p1", "rho2rho1", "T2T1", "p02p01", "p02p1"];
  
  try {
    switch (inputType) {
      case "M1":
        M1 = inputValue;
        if (M1 <= 1) throw new Error("M₁ must be > 1 for a normal shock.");
        break;
      case "M2":
        M1 = normalShockInverse.M1_from_M2(inputValue, g);
        if (isNaN(M1)) throw new Error("Invalid M₂. Must be < 1 and > sqrt((g-1)/2g).");
        break;
      case "p2_p1":
        M1 = normalShockInverse.M1_from_p2p1(inputValue, g);
        if (isNaN(M1)) throw new Error("Invalid p₂/p₁. Must be > 1.");
        break;
      case "rho2_rho1":
        M1 = normalShockInverse.M1_from_rho2rho1(inputValue, g);
        if (isNaN(M1)) throw new Error(`Invalid ρ₂/ρ₁. Must be 1 < ratio < ${(g+1)/(g-1)}.`);
        break;
      case "T2_T1":
        M1 = normalShockInverse.M1_from_T2T1(inputValue, g);
        if (isNaN(M1)) throw new Error("Invalid T₂/T₁. Must be > 1.");
        break;
      case "p02_p01":
        M1 = normalShockInverse.M1_from_p02p01(inputValue, g);
        if (isNaN(M1)) throw new Error("Invalid p₀₂/p₀₁. Must be < 1 and > 0.");
        break;
      case "p02_p1":
        M1 = normalShockInverse.M1_from_p02p1(inputValue, g);
        if (isNaN(M1)) throw new Error("Invalid p₀₂/p₁. Must be > 1.");
        break;
    }
    if (isNaN(M1) || M1 <= 1) throw new Error("Calculation failed. Ensure input is in a valid range.");
    
    const { M2, p2p1, rho2rho1, T2T1, p02p01, p02p1 } = normalShockRatios(M1, g);
    updateText("ns-M1", M1);
    updateText("ns-M2", M2);
    updateText("ns-p2p1", p2p1);
    updateText("ns-rho2rho1", rho2rho1);
    updateText("ns-T2T1", T2T1);
    updateText("ns-p02p01", p02p01);
    updateText("ns-p02p1", p02p1);
    updateNormalShockDiagram(M1, { M2, p2p1, rho2rho1, T2T1, p02p01 });
    flashResults("normal-shock");
  } catch (err) {
    warning.textContent = err.message;
    clearResults("ns", resultIDs);
    updateNormalShockDiagram(NaN, null);
  }
}

// --- 3. Oblique Shock Controller ---
function previewObliqueDiagram() {
  const g = getGamma();
  const M1 = parseFloat(document.getElementById("ob-M1").value);
  const theta = parseFloat(document.getElementById("ob-theta").value);
  const betaIn = parseFloat(document.getElementById("ob-beta").value);
  if (!(M1 > 1) || !(theta >= 0) || isFinite(betaIn)) { updateObliqueDiagram(theta, betaIn, null); return; }
  const betaRad = ob_beta_from_M1_theta(M1, theta, g);
  if (isNaN(betaRad)) { updateObliqueDiagram(theta, NaN, null); return; }
  const r = obliqueShockRatios(M1, betaRad, deg2rad(theta), g);
  updateObliqueDiagram(theta, rad2deg(betaRad), { M1, M2: r.M2, p2p1: r.p2p1, Mn1: r.Mn1 });
}

function clearObliqueInputs() {
    const inputs = ["ob-M1", "ob-theta", "ob-beta"];
    inputs.forEach(id => {
        document.getElementById(id).value = "";
        document.getElementById(id).classList.remove("is-calculated-input");
    });
    document.getElementById("ob-warning").textContent = "";
}

function calculateObliqueShock() {
  const g = getGamma();
  const M1_in = parseFloat(document.getElementById("ob-M1").value);
  const theta_in = parseFloat(document.getElementById("ob-theta").value);
  const beta_in = parseFloat(document.getElementById("ob-beta").value);
  const warning = document.getElementById("ob-warning");
  warning.textContent = "";
  
  const allInputs = [document.getElementById("ob-M1"), document.getElementById("ob-theta"), document.getElementById("ob-beta")];
  allInputs.forEach(el => el.classList.remove("is-calculated-input"));

  const resultIDs = ["M1n", "M2n", "M2", "p2p1", "rho2rho1", "T2T1", "p02p01"];
  let M1 = M1_in, theta = theta_in, beta = beta_in;
  let betaRad, thetaRad, calculatedField = "";

  try {
    const inputs = [!isNaN(M1), !isNaN(theta), !isNaN(beta)];
    const inputCount = inputs.filter(Boolean).length;
    if (inputCount !== 2) {
      throw new Error("Please provide exactly two inputs (M₁, θ, or β).");
    }

    if (!isNaN(M1) && !isNaN(theta)) { // Case 1: (M1, θ) -> β
      if (M1 <= 1) throw new Error("M₁ must be > 1.");
      if (theta < 0) throw new Error("θ must be ≥ 0.");
      betaRad = ob_beta_from_M1_theta(M1, theta, g);
      beta = rad2deg(betaRad);
      if (isNaN(beta)) throw new Error("No attached shock solution (θ > θ_max).");
      thetaRad = deg2rad(theta);
      calculatedField = "ob-beta";
    } else if (!isNaN(M1) && !isNaN(beta)) { // Case 2: (M1, β) -> θ
      if (M1 <= 1) throw new Error("M₁ must be > 1.");
      if (beta <= 0 || beta > 90) throw new Error("β must be in (0, 90].");
      const mu = rad2deg(Math.asin(1/M1));
      if (beta <= mu) throw new Error(`β (${beta.toFixed(1)}°) must be > μ (${mu.toFixed(1)}°).`);
      theta = ob_theta_from_M1_beta(M1, beta, g);
      if (theta < 0) throw new Error("Detached shock region (strong shock).");
      betaRad = deg2rad(beta);
      thetaRad = deg2rad(theta);
      calculatedField = "ob-theta";
    } else if (!isNaN(theta) && !isNaN(beta)) { // Case 3: (θ, β) -> M1
      if (beta <= theta) throw new Error("β must be > θ.");
      if (theta < 0) throw new Error("θ must be ≥ 0.");
      M1 = ob_M1_from_theta_beta(theta, beta, g);
      if (isNaN(M1)) throw new Error("Invalid (θ, β) pair.");
      const mu = rad2deg(Math.asin(1/M1));
      if (beta <= mu) throw new Error(`Implied Mach wave (β ≤ μ).`);
      betaRad = deg2rad(beta);
      thetaRad = deg2rad(theta);
      calculatedField = "ob-M1";
    }
    
    const { Mn1, Mn2, M2, p2p1, rho2rho1, T2T1, p02p01 } = obliqueShockRatios(M1, betaRad, thetaRad, g);
    if (isNaN(M2)) throw new Error("Calculation failed. Check inputs.");

    updateText("ob-M1", M1, 3);
    updateText("ob-theta", theta, 2);
    updateText("ob-beta", beta, 2);
    if (calculatedField) {
        document.getElementById(calculatedField).classList.add("is-calculated-input");
    }
    
    updateText("ob-M1n", Mn1);
    updateText("ob-M2n", Mn2);
    updateText("ob-M2", M2);
    updateText("ob-p2p1", p2p1);
    updateText("ob-rho2rho1", rho2rho1);
    updateText("ob-T2T1", T2T1);
    updateText("ob-p02p01", p02p01);
    
    updateObliqueDiagram(theta, beta, { M1, M2, p2p1, Mn1 });
    flashResults("oblique-shock");

  } catch (err) {
    warning.textContent = err.message;
    clearResults("ob", resultIDs); 
    if (isNaN(M1_in)) updateText("ob-M1", NaN);
    if (isNaN(theta_in)) updateText("ob-theta", NaN);
    if (isNaN(beta_in)) updateText("ob-beta", NaN);
    updateObliqueDiagram(isNaN(theta) ? theta_in : theta, isNaN(beta) ? beta_in : beta, null);
  }
}

// --- 4. Prandtl-Meyer Controller ---
function calculatePrandtlMeyer() {
  const g = getGamma();
  const M1 = parseFloat(document.getElementById("pm-M1").value);
  const thetaDeg = parseFloat(document.getElementById("pm-theta").value);
  const warning = document.getElementById("pm-warning");
  warning.textContent = "";
  const resultIDs = ["nu1", "mu1", "nu2", "mu2", "M2", "p2p1", "rho2rho1", "T2T1"];
  
  try {
    if (isNaN(M1) || M1 <= 1) throw new Error("M₁ must be > 1 for an expansion.");
    if (isNaN(thetaDeg) || thetaDeg <= 0) throw new Error("θ must be > 0.");

    const nu1 = prandtlMeyer.nu(M1, g);
    const mu1 = rad2deg(Math.asin(1/M1));
    const nu2 = nu1 + thetaDeg;
    
    const nu_max = prandtlMeyer.nu(Infinity, g);
    if (nu2 > nu_max) {
      throw new Error(`Expansion angle is too large. Max ν is ${nu_max.toFixed(1)}°.`);
    }
    
    const M2 = prandtlMeyer.M_from_nu(nu2, g);
    const mu2 = M2 === Infinity ? 0 : rad2deg(Math.asin(1/M2));

    const p2p1 = isentropic.p0_p(M1, g) / isentropic.p0_p(M2, g);
    const rho2rho1 = isentropic.rho0_rho(M1, g) / isentropic.rho0_rho(M2, g);
    const T2T1 = isentropic.T0_T(M1, g) / isentropic.T0_T(M2, g);

    updateText("pm-nu1", nu1, 2);
    updateText("pm-mu1", mu1, 2);
    updateText("pm-nu2", nu2, 2);
    updateText("pm-mu2", mu2, 2);
    updateText("pm-M2", M2);
    updateText("pm-p2p1", p2p1);
    updateText("pm-rho2rho1", rho2rho1);
    updateText("pm-T2T1", T2T1);
    flashResults("prandtl-meyer");

  } catch (err) {
    warning.textContent = err.message;
    clearResults("pm", resultIDs);
  }
}

// --- 5. Wedge Airfoil Controller ---
function calculateWedgeAirfoil() {
    const g = getGamma();
    const M1 = parseFloat(document.getElementById("wa-M1").value);
    const delta = parseFloat(document.getElementById("wa-delta").value);
    const alpha = parseFloat(document.getElementById("wa-alpha").value);
    const warning = document.getElementById("wa-warning");
    warning.textContent = "";

    const lowerResultIDs = ["wa-theta-lower", "wa-p-lower-p1", "wa-beta-lower", "wa-m1n-lower", "wa-m2-lower"];
    const upperShockResultIDs = ["wa-theta-upper-shock", "wa-beta-upper", "wa-m1n-upper", "wa-m2-upper-shock"];
    const upperExpResultIDs = ["wa-theta-upper-exp", "wa-nu1-upper", "wa-nu2-upper", "wa-m2-upper-exp"];
    const upperAlignedResultIDs = ["wa-m2-upper-aligned"];

    let lowerResults = {}, upperResults = {};

    try {
        if (isNaN(M1) || M1 <= 1) throw new Error("M₁ must be > 1.");
        if (isNaN(delta)) throw new Error("Wedge Half-Angle (δ) must be a number.");
        if (isNaN(alpha)) throw new Error("Angle of Attack (α) must be a number.");

        // --- Lower Surface (Compression) ---
        const theta_lower = delta + alpha;
        if (theta_lower < 0) {
            // This case (e.g., negative alpha > delta) becomes an expansion.
            // For simplicity, we'll throw an error, as problem implies compression.
            throw new Error("Lower surface has negative deflection (δ + α < 0). This case is not supported.");
        }
        const beta_lower_rad = ob_beta_from_M1_theta(M1, theta_lower, g);
        if (isNaN(beta_lower_rad)) {
            throw new Error("Lower surface shock is detached (θ_lower > θ_max).");
        }
        const { p2p1: p_lower_p1, Mn1: m1n_lower, M2: m2_lower } = obliqueShockRatios(M1, beta_lower_rad, deg2rad(theta_lower), g);
        
        lowerResults = { type: "Shock", theta: theta_lower, beta: rad2deg(beta_lower_rad), M1n: m1n_lower, M2: m2_lower, p_p1: p_lower_p1 };
        
        updateText("wa-p-lower-p1", p_lower_p1, 4);
        updateText("wa-theta-lower", theta_lower, 2);
        updateText("wa-beta-lower", lowerResults.beta, 2);
        updateText("wa-m1n-lower", m1n_lower, 3);
        updateText("wa-m2-lower", m2_lower, 3);


        // --- Upper Surface (Compression or Expansion) ---
        const theta_eff = delta - alpha;
        let p_upper_p1;

        if (theta_eff > 0) {
            // Case 1/2: Compression
            setVisible("wa-results-upper-shock", true);
            setVisible("wa-results-upper-exp", false);
            setVisible("wa-results-upper-aligned", false);

            const beta_upper_rad = ob_beta_from_M1_theta(M1, theta_eff, g);
            if (isNaN(beta_upper_rad)) {
                throw new Error("Upper surface shock is detached (θ_upper > θ_max).");
            }
            const { p2p1, Mn1: m1n_upper, M2: m2_upper } = obliqueShockRatios(M1, beta_upper_rad, deg2rad(theta_eff), g);
            p_upper_p1 = p2p1;
            
            upperResults = { type: "Shock", theta: theta_eff, beta: rad2deg(beta_upper_rad), M1n: m1n_upper, M2: m2_upper, p_p1: p_upper_p1 };
            
            updateText("wa-theta-upper-shock", theta_eff, 2);
            updateText("wa-beta-upper", upperResults.beta, 2);
            updateText("wa-m1n-upper", m1n_upper, 3);
            updateText("wa-m2-upper-shock", m2_upper, 3);

        } else if (theta_eff < 0) {
            // Case 4: Expansion
            setVisible("wa-results-upper-shock", false);
            setVisible("wa-results-upper-exp", true);
            setVisible("wa-results-upper-aligned", false);

            const theta_exp = -theta_eff;
            const nu1 = prandtlMeyer.nu(M1, g);
            const nu2 = nu1 + theta_exp;
            const nu_max = prandtlMeyer.nu(Infinity, g);
            if (nu2 > nu_max) {
                throw new Error(`Upper expansion angle is too large. Max ν is ${nu_max.toFixed(1)}°.`);
            }
            const M2_upper = prandtlMeyer.M_from_nu(nu2, g);
            p_upper_p1 = isentropic.p0_p(M1, g) / isentropic.p0_p(M2_upper, g);

            upperResults = { type: "Expansion", theta: theta_exp, nu1: nu1, nu2: nu2, M2: M2_upper, p_p1: p_upper_p1 };

            updateText("wa-theta-upper-exp", theta_exp, 2);
            updateText("wa-nu1-upper", nu1, 2);
            updateText("wa-nu2-upper", nu2, 2);
            updateText("wa-m2-upper-exp", M2_upper, 3);
        } else {
            // Case 3: Aligned Flow
            setVisible("wa-results-upper-shock", false);
            setVisible("wa-results-upper-exp", false);
            setVisible("wa-results-upper-aligned", true);

            p_upper_p1 = 1.0;
            upperResults = { type: "Aligned", theta: 0, M2: M1, p_p1: 1.0 };
            updateText("wa-m2-upper-aligned", M1, 3);
        }
        
        updateText("wa-p-upper-p1", p_upper_p1, 4);
        
        flashResults("wedge-airfoil");
        updateWedgeDiagram(g, M1, delta, alpha, lowerResults, upperResults);

    } catch (err) {
        warning.textContent = err.message;
        clearResults("wa", ["p-lower-p1", ...lowerResultIDs]);
        clearResults("wa", ["p-upper-p1"]);
        setVisible("wa-results-upper-shock", false);
        setVisible("wa-results-upper-exp", false);
        setVisible("wa-results-upper-aligned", false);
        document.getElementById("wa-diagram").innerHTML = ''; // Clear diagram
    }
}


// ========================================================================
// === Initialization =====================================================
// ========================================================================

function runAllCalculators() {
  calculateIsentropic();
  calculateNormalShock();
  // Oblique shock is not computed on load (it fills the third input), but its
  // figure shows the default case so the tab never opens on an empty drawing.
  previewObliqueDiagram();
  calculatePrandtlMeyer();
  calculateWedgeAirfoil();
}

document.addEventListener("DOMContentLoaded", () => {
  // --- Main Tab Navigation ---
  const navButtons = document.querySelectorAll(".nav-btn");
  const panels = document.querySelectorAll(".calculator-panel");
  // Each tab has its own address (index.html#oblique-shock, ...), so the
  // Supersonic Airfoils page can link straight to any calculator here.
  const showCalculator = (name, updateHash) => {
    const panel = document.getElementById(`${name}-calculator`);
    if (!panel) return false;
    navButtons.forEach(b => {
      const on = b.dataset.calculator === name;
      b.classList.toggle("active", on);
      if (on) b.setAttribute("aria-current", "page"); else b.removeAttribute("aria-current");
    });
    panels.forEach(p => p.classList.remove("active"));
    panel.classList.add("active");
    if (updateHash) {
      try { history.replaceState(null, "", "#" + name); } catch (e) { location.hash = name; }
    }
    return true;
  };
  navButtons.forEach(btn => {
    if (!btn.dataset.calculator) return; // links to other pages (Supersonic Airfoils)
    btn.addEventListener("click", () => showCalculator(btn.dataset.calculator, true));
  });
  showCalculator(location.hash.slice(1), false) || showCalculator("isentropic", false);
  window.addEventListener("hashchange", () => showCalculator(location.hash.slice(1), false));

  // --- Theme Toggle ---
  const themeToggle = document.getElementById("theme-toggle-checkbox");
  const storedTheme = localStorage.getItem("theme");
  const prefersLight = window.matchMedia("(prefers-color-scheme: light)").matches;
  const applyTheme = (theme) => {
    document.body.dataset.theme = theme;
    themeToggle.checked = (theme === "light");
  };
  if (storedTheme) {
    applyTheme(storedTheme);
  } else {
    applyTheme(prefersLight ? "light" : "dark");
  }
  themeToggle.addEventListener("change", () => {
    const newTheme = themeToggle.checked ? "light" : "dark";
    localStorage.setItem("theme", newTheme);
    applyTheme(newTheme);
  });

  // --- Font Size Control ---
  const fontIncrease = document.getElementById("font-increase");
  const fontDecrease = document.getElementById("font-decrease");
  const fontSizes = ["xs", "sm", "md", "lg", "xl"];
  let currentSizeIndex = fontSizes.indexOf(localStorage.getItem("font-size") || "md");
  if (currentSizeIndex === -1) currentSizeIndex = 2; // Default to 'md'
  
  const setFontSize = (index) => {
    if (index < 0) index = 0;
    if (index >= fontSizes.length) index = fontSizes.length - 1;
    currentSizeIndex = index;
    const size = fontSizes[index];
    document.body.dataset.fontSize = size;
    localStorage.setItem("font-size", size);
  };
  
  fontIncrease.addEventListener("click", () => setFontSize(currentSizeIndex + 1));
  fontDecrease.addEventListener("click", () => setFontSize(currentSizeIndex - 1));
  
  // Set initial font size
  setFontSize(currentSizeIndex);


  // --- Setup Calculator Event Listeners ---
  document.getElementById("is-computeBtn").addEventListener("click", calculateIsentropic);
  document.getElementById("is-inputType").addEventListener("change", calculateIsentropic);
  
  document.getElementById("ns-computeBtn").addEventListener("click", calculateNormalShock);
  document.getElementById("ns-inputType").addEventListener("change", calculateNormalShock);
  
  document.getElementById("ob-computeBtn").addEventListener("click", calculateObliqueShock);
  document.getElementById("ob-clearBtn").addEventListener("click", clearObliqueInputs);
  
  document.getElementById("pm-computeBtn").addEventListener("click", calculatePrandtlMeyer);

  document.getElementById("wa-computeBtn").addEventListener("click", calculateWedgeAirfoil);
  
  document.getElementById("gamma").addEventListener("change", runAllCalculators);

  // --- Run calculators on initial load ---
  runAllCalculators();
});

