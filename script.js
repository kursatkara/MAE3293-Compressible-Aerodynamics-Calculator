\
"use strict";

// ---------------------------
// Utility helpers
// ---------------------------
function getGamma() {
  var gammaInput = document.getElementById("gamma");
  var g = parseFloat(gammaInput && gammaInput.value);
  if (!isFinite(g) || g <= 1) {
    return 1.4;
  }
  return g;
}

function readNumber(id, opts) {
  opts = opts || {};
  var el = document.getElementById(id);
  if (!el) return NaN;
  var v = parseFloat(el.value);
  if (!isFinite(v)) return NaN;
  if (opts.min !== undefined && v < opts.min) return NaN;
  if (opts.max !== undefined && v > opts.max) return NaN;
  return v;
}

function fmt(x) {
  if (!isFinite(x)) return "–";
  var ax = Math.abs(x);
  if (ax !== 0 && (ax < 1e-3 || ax > 1e4)) {
    return x.toExponential(4);
  }
  return x.toFixed(4);
}

function setText(id, val) {
  var el = document.getElementById(id);
  if (el) el.textContent = fmt(val);
}

function showWarning(id, msg) {
  var el = document.getElementById(id);
  if (el) el.textContent = msg || "";
}

// ---------------------------
// Basic gas-dynamics relations
// ---------------------------
function tRatioFromM(M, g) {
  return 1 / (1 + 0.5 * (g - 1) * M * M);
}

function pRatioFromM(M, g) {
  var t = tRatioFromM(M, g);
  return Math.pow(t, g / (g - 1));
}

function rhoRatioFromM(M, g) {
  var t = tRatioFromM(M, g);
  return Math.pow(t, 1 / (g - 1));
}

function areaRatioFromM(M, g) {
  var term1 = 1 / M;
  var term2 = (2 / (g + 1)) * (1 + 0.5 * (g - 1) * M * M);
  var exponent = (g + 1) / (2 * (g - 1));
  return term1 * Math.pow(term2, exponent);
}

// Inversion helpers (isentropic)
function MFromTRatio(T_T0, g) {
  if (T_T0 <= 0 || T_T0 >= 1) return NaN;
  var a = 0.5 * (g - 1);
  var M2 = (1 / T_T0 - 1) / a;
  return M2 > 0 ? Math.sqrt(M2) : NaN;
}

function MFromPRatio(P_P0, g) {
  if (P_P0 <= 0 || P_P0 >= 1) return NaN;
  var exponent = (g - 1) / g;
  var t = Math.pow(P_P0, exponent);
  var a = 0.5 * (g - 1);
  var M2 = (1 / t - 1) / a;
  return M2 > 0 ? Math.sqrt(M2) : NaN;
}

function MFromRhoRatio(rho_rho0, g) {
  if (rho_rho0 <= 0 || rho_rho0 >= 1) return NaN;
  var exponent = g - 1;
  var t = Math.pow(rho_rho0, exponent);
  var a = 0.5 * (g - 1);
  var M2 = (1 / t - 1) / a;
  return M2 > 0 ? Math.sqrt(M2) : NaN;
}

function bisection(fn, a, b, maxIter, tol) {
  maxIter = maxIter || 60;
  tol = tol || 1e-8;
  var fa = fn(a);
  var fb = fn(b);
  if (!isFinite(fa) || !isFinite(fb) || fa * fb > 0) return NaN;
  var mid, fmid;
  for (var i = 0; i < maxIter; i++) {
    mid = 0.5 * (a + b);
    fmid = fn(mid);
    if (!isFinite(fmid)) return NaN;
    if (Math.abs(fmid) < tol || (b - a) * 0.5 < tol) return mid;
    if (fa * fmid <= 0) {
      b = mid;
      fb = fmid;
    } else {
      a = mid;
      fa = fmid;
    }
  }
  return mid;
}

function MFromAreaRatio(A_Astar, g, regime) {
  if (A_Astar < 1) return NaN;
  var fn = function (M) {
    return areaRatioFromM(M, g) - A_Astar;
  };
  if (regime === "subsonic") {
    return bisection(fn, 1e-6, 0.999);
  }
  return bisection(fn, 1.0001, 50);
}

// ---------------------------
// Normal shock
// ---------------------------
function normalShock(M1, g) {
  if (M1 <= 1) return null;

  var M1n2 = M1 * M1;
  var numerator = 1 + 0.5 * (g - 1) * M1n2;
  var denominator = g * M1n2 - 0.5 * (g - 1);
  if (denominator <= 0) return null;
  var M2 = Math.sqrt(numerator / denominator);

  var P2_P1 = 1 + (2 * g / (g + 1)) * (M1n2 - 1);
  var rho2_rho1 = ((g + 1) * M1n2) / ((g - 1) * M1n2 + 2);
  var T2_T1 = P2_P1 / rho2_rho1;

  var P01_P1 = Math.pow(1 + 0.5 * (g - 1) * M1n2, g / (g - 1));
  var P02_P2 = Math.pow(1 + 0.5 * (g - 1) * M2 * M2, g / (g - 1));
  var P02_P01 = (P2_P1 * P02_P2) / P01_P1;

  return {
    M2: M2,
    P2_P1: P2_P1,
    T2_T1: T2_T1,
    rho2_rho1: rho2_rho1,
    P02_P01: P02_P01
  };
}

// ---------------------------
// Oblique shock
// ---------------------------
function thetaFromBeta(M1, beta, g) {
  var sinB = Math.sin(beta);
  var Mn1_sq = M1 * M1 * sinB * sinB;
  var num = 2 * (Mn1_sq - 1) / Math.tan(beta);
  var den = M1 * M1 * (g + Math.cos(2 * beta)) + 2;
  var tanTheta = num / den;
  if (!isFinite(tanTheta) || tanTheta < 0) return NaN;
  return Math.atan(tanTheta);
}

function obliqueShockMaxTheta(M1, g) {
  if (M1 <= 1) return { thetaMax: 0, betaAtMax: NaN };
  var mu = Math.asin(1 / M1) + 1e-4;
  var betaMax = 0.5 * Math.PI - 1e-4;
  var N = 800;
  var thetaMax = 0;
  var betaAtMax = mu;

  for (var i = 0; i <= N; i++) {
    var beta = mu + (betaMax - mu) * (i / N);
    var theta = thetaFromBeta(M1, beta, g);
    if (isFinite(theta) && theta > thetaMax) {
      thetaMax = theta;
      betaAtMax = beta;
    }
  }
  return { thetaMax: thetaMax, betaAtMax: betaAtMax };
}

function findObliqueShockBeta(M1, thetaTargetRad, g) {
  var mu = Math.asin(1 / M1) + 1e-4;
  var betaMax = 0.5 * Math.PI - 1e-4;
  var N = 800;
  var brackets = [];

  var prevBeta = mu;
  var prevTheta = thetaFromBeta(M1, prevBeta, g);
  var prevF = prevTheta - thetaTargetRad;

  for (var i = 1; i <= N; i++) {
    var beta = mu + (betaMax - mu) * (i / N);
    var theta = thetaFromBeta(M1, beta, g);
    if (!isFinite(theta)) {
      prevBeta = beta;
      prevTheta = theta;
      prevF = NaN;
      continue;
    }
    var f = theta - thetaTargetRad;
    if (isFinite(prevF) && prevF * f <= 0) {
      brackets.push([prevBeta, beta]);
    }
    prevBeta = beta;
    prevTheta = theta;
    prevF = f;
  }
  return brackets;
}

function obliqueShock(M1, thetaDeg, g, branch) {
  if (M1 <= 1) {
    return { error: "M₁ must be supersonic for an attached oblique shock." };
  }
  var thetaRad = (thetaDeg * Math.PI) / 180;
  var maxInfo = obliqueShockMaxTheta(M1, g);
  var thetaMaxDeg = (maxInfo.thetaMax * 180) / Math.PI;

  if (thetaDeg <= 0 || !isFinite(thetaDeg)) {
    return { error: "Deflection angle θ must be positive.", thetaMaxDeg: thetaMaxDeg };
  }
  if (thetaDeg > thetaMaxDeg + 1e-6) {
    return {
      error: "θ exceeds θ_max; attached oblique shock does not exist (detached shock).",
      thetaMaxDeg: thetaMaxDeg
    };
  }

  var brackets = findObliqueShockBeta(M1, thetaRad, g);
  if (!brackets.length) {
    return {
      error: "Could not find a valid oblique shock solution for this θ.",
      thetaMaxDeg: thetaMaxDeg
    };
  }

  var betaWeak, betaStrong;
  if (brackets.length === 1) {
    betaWeak = bisection(function (b) {
      return thetaFromBeta(M1, b, g) - thetaRad;
    }, brackets[0][0], brackets[0][1]);
    betaStrong = betaWeak;
  } else {
    var beta1 = bisection(function (b) {
      return thetaFromBeta(M1, b, g) - thetaRad;
    }, brackets[0][0], brackets[0][1]);
    var beta2 = bisection(function (b) {
      return thetaFromBeta(M1, b, g) - thetaRad;
    }, brackets[1][0], brackets[1][1]);
    if (!isFinite(beta1) || !isFinite(beta2)) {
      return {
        error: "Failed to converge while solving for β.",
        thetaMaxDeg: thetaMaxDeg
      };
    }
    betaWeak = Math.min(beta1, beta2);
    betaStrong = Math.max(beta1, beta2);
  }

  var beta = branch === "strong" ? betaStrong : betaWeak;
  var sinB = Math.sin(beta);
  var M1n = M1 * sinB;
  var ns = normalShock(M1n, g);
  if (!ns) {
    return {
      error: "Failed to compute normal-shock quantities from M₁n.",
      thetaMaxDeg: thetaMaxDeg
    };
  }
  var thetaRadFinal = thetaFromBeta(M1, beta, g);
  var M2 = ns.M2 / Math.sin(beta - thetaRadFinal);

  return {
    thetaMaxDeg: thetaMaxDeg,
    betaDeg: beta * 180 / Math.PI,
    M2: M2,
    P2_P1: ns.P2_P1,
    T2_T1: ns.T2_T1,
    rho2_rho1: ns.rho2_rho1,
    P02_P01: ns.P02_P01
  };
}

// ---------------------------
// Prandtl-Meyer expansion
// ---------------------------
function nuFromM(M, g) {
  if (M < 1) return NaN;
  var gm1 = g - 1;
  var gp1 = g + 1;
  var a = Math.sqrt(gp1 / gm1);
  var inner = Math.sqrt((gm1 / gp1) * (M * M - 1));
  return a * Math.atan(inner) - Math.atan(Math.sqrt(M * M - 1));
}

function MFromNu(nuTarget, g) {
  if (nuTarget <= 0) return NaN;
  var fn = function (M) {
    return nuFromM(M, g) - nuTarget;
  };
  return bisection(fn, 1.0001, 50);
}

// ---------------------------
// Public DOM helpers
// ---------------------------
function updateIsentropicLabel() {
  var select = document.getElementById("iso-input-type");
  var label = document.getElementById("iso-input-label");
  var regime = document.getElementById("iso-area-regime");
  if (!select || !label || !regime) return;
  var v = select.value;
  if (v === "M") {
    label.textContent = "Mach Number (M)";
    regime.style.display = "none";
  } else if (v === "T_T0") {
    label.textContent = "T/T₀";
    regime.style.display = "none";
  } else if (v === "P_P0") {
    label.textContent = "P/P₀";
    regime.style.display = "none";
  } else if (v === "rho_rho0") {
    label.textContent = "ρ/ρ₀";
    regime.style.display = "none";
  } else {
    label.textContent = "A/A*";
    regime.style.display = "block";
  }
}

function calculateIsentropic() {
  var g = getGamma();
  var known = document.getElementById("iso-input-type").value;
  var val = readNumber("iso-input-value", { min: 0 });
  if (!isFinite(val) || val <= 0) {
    alert("Please enter a positive value for the isentropic input.");
    return;
  }

  var M;
  if (known === "M") {
    M = val;
  } else if (known === "T_T0") {
    M = MFromTRatio(val, g);
  } else if (known === "P_P0") {
    M = MFromPRatio(val, g);
  } else if (known === "rho_rho0") {
    M = MFromRhoRatio(val, g);
  } else if (known === "A_Astar") {
    var radios = document.querySelectorAll("input[name='iso-area-branch']");
    var regimeValue = "supersonic";
    for (var i = 0; i < radios.length; i++) {
      if (radios[i].checked) {
        regimeValue = radios[i].value;
        break;
      }
    }
    M = MFromAreaRatio(val, g, regimeValue);
  }

  if (!isFinite(M) || M <= 0) {
    alert("No valid Mach number could be found for the given input and γ.");
    return;
  }

  var T_T0 = tRatioFromM(M, g);
  var P_P0 = pRatioFromM(M, g);
  var rho_rho0 = rhoRatioFromM(M, g);
  var A_Astar = areaRatioFromM(M, g);

  setText("iso-M", M);
  setText("iso-T_T0", T_T0);
  setText("iso-P_P0", P_P0);
  setText("iso-rho_rho0", rho_rho0);
  setText("iso-A_Astar", A_Astar);
}

function calculateNormalShock() {
  var g = getGamma();
  var M1 = readNumber("ns-M1", { min: 1.0 });
  if (!isFinite(M1) || M1 <= 1) {
    alert("Please enter a supersonic upstream Mach number M₁ > 1.");
    return;
  }
  var res = normalShock(M1, g);
  if (!res) {
    alert("Failed to compute normal-shock relations for this Mach number and γ.");
    return;
  }

  setText("ns-M2", res.M2);
  setText("ns-P2_P1", res.P2_P1);
  setText("ns-T2_T1", res.T2_T1);
  setText("ns-rho2_rho1", res.rho2_rho1);
  setText("ns-P02_P01", res.P02_P01);
}

function calculateObliqueShock() {
  var g = getGamma();
  var M1 = readNumber("os-M1", { min: 1.0 });
  var thetaDeg = readNumber("os-theta", { min: 0 });
  if (!isFinite(M1) || M1 <= 1) {
    alert("Please enter a supersonic upstream Mach number M₁ > 1.");
    return;
  }
  if (!isFinite(thetaDeg) || thetaDeg <= 0) {
    alert("Please enter a positive deflection angle θ in degrees.");
    return;
  }

  var radios = document.querySelectorAll("input[name='os-branch']");
  var branch = "weak";
  for (var i = 0; i < radios.length; i++) {
    if (radios[i].checked) {
      branch = radios[i].value;
      break;
    }
  }

  var result = obliqueShock(M1, thetaDeg, g, branch);
  showWarning("os-warning", "");

  if (result.error) {
    showWarning("os-warning", result.error);
    setText("os-theta-max", result.thetaMaxDeg || NaN);
    setText("os-beta", NaN);
    setText("os-M2", NaN);
    setText("os-P2_P1", NaN);
    setText("os-T2_T1", NaN);
    setText("os-rho2_rho1", NaN);
    setText("os-P02_P01", NaN);
    return;
  }

  setText("os-theta-max", result.thetaMaxDeg);
  setText("os-beta", result.betaDeg);
  setText("os-M2", result.M2);
  setText("os-P2_P1", result.P2_P1);
  setText("os-T2_T1", result.T2_T1);
  setText("os-rho2_rho1", result.rho2_rho1);
  setText("os-P02_P01", result.P02_P01);
}

function calculatePrandtlMeyer() {
  var g = getGamma();
  var type = document.getElementById("pm-input-type").value;
  var val = readNumber("pm-input-value", { min: 0 });
  if (!isFinite(val) || val <= 0) {
    alert("Please enter a positive value.");
    return;
  }

  var M, nu;
  if (type === "M") {
    M = val;
    if (M <= 1) {
      alert("Prandtl–Meyer expansion is defined for supersonic M > 1.");
      return;
    }
    nu = nuFromM(M, g);
  } else {
    var nuTarget = (val * Math.PI) / 180;
    M = MFromNu(nuTarget, g);
    if (!isFinite(M) || M <= 1) {
      alert("Could not find a valid Mach number for this ν and γ.");
      return;
    }
    nu = nuFromM(M, g);
  }

  setText("pm-M", M);
  setText("pm-nu", nu * 180 / Math.PI);
}

function updatePmLabel() {
  var select = document.getElementById("pm-input-type");
  var label = document.getElementById("pm-input-label");
  if (!select || !label) return;
  if (select.value === "M") {
    label.textContent = "Mach Number (M)";
  } else {
    label.textContent = "ν [degrees]";
  }
}

// Initialize labels on load
if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", function () {
    updateIsentropicLabel();
    updatePmLabel();
  });
} else {
  updateIsentropicLabel();
  updatePmLabel();
}
