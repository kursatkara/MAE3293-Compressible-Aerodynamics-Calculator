// =======================================================
// Global Constants & Helper Functions
// =======================================================
let GAMMA = parseFloat(document.getElementById('gamma-input').value); // Global gamma, can be changed by user
const DEG_TO_RAD = Math.PI / 180;
const RAD_TO_DEG = 180 / Math.PI;
const DISPLAY_PRECISION = 5; // Number of decimal places for results

// Utility function to format numbers for display
function formatNumber(num) {
    if (typeof num !== 'number' || isNaN(num) || !isFinite(num)) {
        return 'N/A';
    }
    return num.toFixed(DISPLAY_PRECISION);
}

// Function to display error messages
function displayError(elementId, message) {
    const errorElement = document.getElementById(elementId);
    errorElement.textContent = message;
    errorElement.classList.add('visible');
}

// Function to hide error messages
function hideError(elementId) {
    const errorElement = document.getElementById(elementId);
    errorElement.textContent = '';
    errorElement.classList.remove('visible');
}

// Update GAMMA when the input changes
document.getElementById('gamma-input').addEventListener('input', (event) => {
    const newGamma = parseFloat(event.target.value);
    if (!isNaN(newGamma) && newGamma > 1.0) { // Gamma must be > 1
        GAMMA = newGamma;
        // Re-run calculations that might be auto-updating or suggest recalculation
        // For this project, we'll rely on explicit button clicks after gamma change.
        console.log(`Gamma updated to: ${GAMMA}`);
    } else {
        displayError('isentropic-error', 'Gamma must be a number greater than 1.');
    }
});

// =======================================================
// Aerodynamic Calculation Functions (Verified from previous step)
// =======================================================

// --- Isentropic Flow Relations ---
function isentropicT_T0(M) {
    return 1 / (1 + ((GAMMA - 1) / 2) * M * M);
}
function isentropicP_P0(M) {
    return Math.pow(isentropicT_T0(M), GAMMA / (GAMMA - 1));
}
function isentropicRho_Rho0(M) {
    return Math.pow(isentropicT_T0(M), 1 / (GAMMA - 1));
}
function isentropicA_Astar(M) {
    if (M <= 0) return NaN;
    const term1 = (GAMMA + 1) / 2;
    const term2 = 1 + ((GAMMA - 1) / 2) * M * M;
    const exponent = (GAMMA + 1) / (2 * (GAMMA - 1));
    return (1 / M) * Math.pow(term2 / term1, exponent);
}

function getM_from_T_T0(T_T0_ratio) {
    if (T_T0_ratio <= 0 || T_T0_ratio > 1) return NaN;
    return Math.sqrt((2 / (GAMMA - 1)) * (1 / T_T0_ratio - 1));
}
function getM_from_P_P0(P_P0_ratio) {
    if (P_P0_ratio <= 0 || P_P0_ratio > 1) return NaN;
    return Math.sqrt((2 / (GAMMA - 1)) * (Math.pow(1 / P_P0_ratio, (GAMMA - 1) / GAMMA) - 1));
}
function getM_from_Rho_Rho0(Rho_Rho0_ratio) {
    if (Rho_Rho0_ratio <= 0 || Rho_Rho0_ratio > 1) return NaN;
    return Math.sqrt((2 / (GAMMA - 1)) * (Math.pow(1 / Rho_Rho0_ratio, GAMMA - 1) - 1));
}
function getM_from_A_Astar(A_Astar_ratio, flowType) {
    if (A_Astar_ratio < 1) return NaN;
    let M_guess;
    if (flowType === 'subsonic') { M_guess = 0.5; }
    else if (flowType === 'supersonic') { M_guess = 2.0; }
    else { return NaN; }

    const maxIterations = 100;
    const tolerance = 1e-6;

    for (let i = 0; i < maxIterations; i++) {
        const F_M = isentropicA_Astar(M_guess) - A_Astar_ratio;
        if (Math.abs(F_M) < tolerance) { return M_guess; }
        const h = 1e-6;
        const dF_dM = (isentropicA_Astar(M_guess + h) - isentropicA_Astar(M_guess - h)) / (2 * h);
        if (dF_dM === 0 || isNaN(dF_dM)) return NaN;
        M_guess -= F_M / dF_dM;
        if (M_guess <= 0) M_guess = 0.01;
    }
    return NaN;
}

// --- Normal Shock Relations ---
function normalShockM2(M1) {
    if (M1 < 1) return NaN;
    const numerator = 1 + ((GAMMA - 1) / 2) * M1 * M1;
    const denominator = GAMMA * M1 * M1 - (GAMMA - 1) / 2;
    return Math.sqrt(numerator / denominator);
}
function normalShockP2_P1(M1) {
    if (M1 < 1) return NaN;
    return 1 + (2 * GAMMA / (GAMMA + 1)) * (M1 * M1 - 1);
}
function normalShockT2_T1(M1) {
    if (M1 < 1) return NaN;
    const M2 = normalShockM2(M1);
    if (isNaN(M2)) return NaN;
    return (1 + ((GAMMA - 1) / 2) * M1 * M1) / (1 + ((GAMMA - 1) / 2) * M2 * M2);
}
function normalShockRho2_Rho1(M1) {
    if (M1 < 1) return NaN;
    return normalShockP2_P1(M1) / normalShockT2_T1(M1);
}
function normalShockP02_P01(M1) {
    if (M1 < 1) return NaN;
    const P2_P1 = normalShockP2_P1(M1);
    const M2 = normalShockM2(M1);
    const P01_P1_ratio = Math.pow(1 + ((GAMMA - 1) / 2) * M1 * M1, GAMMA / (GAMMA - 1));
    const P02_P2_ratio = Math.pow(1 + ((GAMMA - 1) / 2) * M2 * M2, GAMMA / (GAMMA - 1));
    return (P02_P2_ratio / P01_P1_ratio) * P2_P1;
}

// --- Oblique Shock Relations ---
function calculateMaxTheta(M1) {
    if (M1 < 1) return NaN;
    let beta_rad_search_start = Math.asin(1 / M1); // Mach angle
    let beta_rad_search_end = Math.PI / 2; // Normal shock
    let max_theta = 0;

    const numPoints = 1000;
    for (let i = 0; i < numPoints; i++) {
        const beta_rad = beta_rad_search_start + (beta_rad_search_end - beta_rad_search_start) * i / numPoints;
        const M1n_test = M1 * Math.sin(beta_rad);
        if (M1n_test < 1) continue;

        const theta_test = Math.atan(
            (2 * (1 / Math.tan(beta_rad)) * (M1 * M1 * Math.sin(beta_rad) * Math.sin(beta_rad) - 1)) /
            (M1 * M1 * (GAMMA + Math.cos(2 * beta_rad)) + 2)
        );
        max_theta = Math.max(max_theta, theta_test);
    }
    return max_theta * RAD_TO_DEG;
}

function calculateBeta(M1, theta_deg, shockType) {
    if (M1 < 1) return NaN;
    if (theta_deg < 0) return NaN;

    const theta_rad = theta_deg * DEG_TO_RAD;
    const maxTheta = calculateMaxTheta(M1);

    if (theta_deg > maxTheta + 0.001) { // Add small tolerance for floating point
        return NaN; // Shock is detached
    }

    let beta_guess_rad;
    if (shockType === 'weak') {
        const mu1 = Math.asin(1 / M1);
        beta_guess_rad = mu1 + 0.05 * DEG_TO_RAD; // Start slightly above Mach angle
    } else if (shockType === 'strong') {
        beta_guess_rad = (85 - theta_deg / 2) * DEG_TO_RAD; // Heuristic guess for strong shock
        if (beta_guess_rad < 5 * DEG_TO_RAD) beta_guess_rad = 5 * DEG_TO_RAD; // Prevent too small
    } else {
        return NaN;
    }

    const maxIterations = 200;
    const tolerance = 1e-7;

    for (let i = 0; i < maxIterations; i++) {
        if (beta_guess_rad <= 0 || beta_guess_rad >= Math.PI / 2) {
            beta_guess_rad = Math.max(0.01 * DEG_TO_RAD, Math.min(89.9 * DEG_TO_RAD, beta_guess_rad));
        }

        const sin_sq_beta = Math.sin(beta_guess_rad) * Math.sin(beta_guess_rad);
        const F_beta = (2 * (1 / Math.tan(beta_guess_rad)) * (M1 * M1 * sin_sq_beta - 1)) / (M1 * M1 * (GAMMA + Math.cos(2 * beta_guess_rad)) + 2) - Math.tan(theta_rad);

        if (Math.abs(F_beta) < tolerance) {
            return beta_guess_rad * RAD_TO_DEG;
        }

        const h = 1e-7;
        const beta_plus_h = beta_guess_rad + h;
        const beta_minus_h = beta_guess_rad - h;

        const sin_sq_beta_plus_h = Math.sin(beta_plus_h) * Math.sin(beta_plus_h);
        const F_beta_plus_h = (2 * (1 / Math.tan(beta_plus_h)) * (M1 * M1 * sin_sq_beta_plus_h - 1)) / (M1 * M1 * (GAMMA + Math.cos(2 * beta_plus_h)) + 2) - Math.tan(theta_rad);

        const sin_sq_beta_minus_h = Math.sin(beta_minus_h) * Math.sin(beta_minus_h);
        const F_beta_minus_h = (2 * (1 / Math.tan(beta_minus_h)) * (M1 * M1 * sin_sq_beta_minus_h - 1)) / (M1 * M1 * (GAMMA + Math.cos(2 * beta_minus_h)) + 2) - Math.tan(theta_rad);

        const dF_dbeta = (F_beta_plus_h - F_beta_minus_h) / (2 * h);

        if (Math.abs(dF_dbeta) < 1e-10) { // Avoid division by very small numbers
            break;
        }
        beta_guess_rad -= F_beta / dF_dbeta;

        beta_guess_rad = Math.max(0.001 * DEG_TO_RAD, Math.min(89.999 * DEG_TO_RAD, beta_guess_rad));
    }
    return NaN;
}

function calculateObliqueShockProperties(M1, beta_deg) {
    if (M1 < 1) return null;
    const beta_rad = beta_deg * DEG_TO_RAD;
    const mu1_rad = Math.asin(1 / M1);

    if (beta_rad < mu1_rad - 0.001 * DEG_TO_RAD || beta_deg > 90) return null; // Small tolerance

    const M1n = M1 * Math.sin(beta_rad);
    if (M1n < 1) return null; // No shock if normal component is subsonic

    const M2n = normalShockM2(M1n);
    const M1t = M1 * Math.cos(beta_rad);
    const M2t = M1t; // Tangential velocity component is unchanged
    const M2 = Math.sqrt(M2n * M2n + M2t * M2t);

    const P2_P1 = normalShockP2_P1(M1n);
    const T2_T1 = normalShockT2_T1(M1n);
    const Rho2_Rho1 = normalShockRho2_Rho1(M1n);
    const P02_P01 = normalShockP02_P01(M1n);

    // Theta can be derived from beta directly, for verification or display
    const theta_rad = Math.atan(
        (2 * (1 / Math.tan(beta_rad)) * (M1 * M1 * Math.sin(beta_rad) * Math.sin(beta_rad) - 1)) /
        (M1 * M1 * (GAMMA + Math.cos(2 * beta_rad)) + 2)
    );

    return {
        M1n: M1n,
        M2: M2,
        M2n: M2n,
        P2_P1: P2_P1,
        T2_T1: T2_T1,
        Rho2_Rho1: Rho2_Rho1,
        P02_P01: P02_P01,
        theta_deg: theta_rad * RAD_TO_DEG
    };
}

// --- Prandtl-Meyer Expansion Fan ---
function prandtlMeyer_nu(M) {
    if (M < 1) return NaN;
    const term1 = Math.sqrt((GAMMA + 1) / (GAMMA - 1));
    const term2 = Math.atan(Math.sqrt(((GAMMA - 1) / (GAMMA + 1)) * (M * M - 1)));
    const term3 = Math.atan(Math.sqrt(M * M - 1));
    return (term1 * term2 - term3) * RAD_TO_DEG; // Result in degrees
}

function getM_from_nu(nu_deg) {
    // Inverse Prandtl-Meyer function requires iteration (Newton-Raphson)
    if (nu_deg < 0) return NaN;

    let M_guess = 1.5; // Initial guess for M > 1
    const maxIterations = 100;
    const tolerance = 1e-6; // for nu_deg

    for (let i = 0; i < maxIterations; i++) {
        const F_M = prandtlMeyer_nu(M_guess) - nu_deg;
        if (Math.abs(F_M) < tolerance) {
            return M_guess;
        }

        const h = 1e-6;
        const dF_dM = (prandtlMeyer_nu(M_guess + h) - prandtlMeyer_nu(M_guess - h)) / (2 * h);
        if (dF_dM === 0 || isNaN(dF_dM)) return NaN;

        M_guess -= F_M / dF_dM;
        if (M_guess < 1) M_guess = 1.0001; // Ensure M remains supersonic
    }
    return NaN;
}


// =======================================================
// UI Interaction and Calculation Logic
// =======================================================

// --- Isentropic Flow Calculator ---
const isentropicInputType = document.getElementById('isentropic-input-type');
const isentropicInputValue = document.getElementById('isentropic-input-value');
const isentropicInputLabel = document.getElementById('isentropic-input-label');
const isentropicFlowTypeGroup = document.getElementById('isentropic-flow-type-group');
const isentropicFlowType = document.getElementById('isentropic-flow-type');

isentropicInputType.addEventListener('change', () => {
    const selectedType = isentropicInputType.value;
    isentropicFlowTypeGroup.classList.add('hidden'); // Hide by default

    switch (selectedType) {
        case 'M':
            isentropicInputLabel.textContent = 'Mach Number (M)';
            isentropicInputValue.min = '0';
            isentropicInputValue.value = '1.0';
            break;
        case 'T_T0':
            isentropicInputLabel.textContent = 'T/T₀ (0 < T/T₀ ≤ 1)';
            isentropicInputValue.min = '0.0001';
            isentropicInputValue.max = '1';
            isentropicInputValue.value = '0.9';
            break;
        case 'P_P0':
            isentropicInputLabel.textContent = 'P/P₀ (0 < P/P₀ ≤ 1)';
            isentropicInputValue.min = '0.0001';
            isentropicInputValue.max = '1';
            isentropicInputValue.value = '0.8';
            break;
        case 'Rho_Rho0':
            isentropicInputLabel.textContent = 'ρ/ρ₀ (0 < ρ/ρ₀ ≤ 1)';
            isentropicInputValue.min = '0.0001';
            isentropicInputValue.max = '1';
            isentropicInputValue.value = '0.85';
            break;
        case 'A_Astar':
            isentropicInputLabel.textContent = 'A/A* (A/A* ≥ 1)';
            isentropicInputValue.min = '1';
            isentropicInputValue.value = '1.5';
            isentropicFlowTypeGroup.classList.remove('hidden'); // Show flow type for A/A*
            break;
    }
    hideError('isentropic-error'); // Clear previous errors
});

function calculateIsentropic() {
    hideError('isentropic-error');
    const inputType = isentropicInputType.value;
    const inputValue = parseFloat(isentropicInputValue.value);
    let M, T_T0_val, P_P0_val, Rho_Rho0_val, A_Astar_val;

    if (isNaN(inputValue)) {
        displayError('isentropic-error', 'Please enter a valid number.');
        return;
    }

    try {
        switch (inputType) {
            case 'M':
                if (inputValue < 0) throw new Error('Mach number cannot be negative.');
                M = inputValue;
                T_T0_val = isentropicT_T0(M);
                P_P0_val = isentropicP_P0(M);
                Rho_Rho0_val = isentropicRho_Rho0(M);
                A_Astar_val = isentropicA_Astar(M);
                break;
            case 'T_T0':
                if (inputValue <= 0 || inputValue > 1) throw new Error('T/T₀ must be between 0 and 1 (exclusive of 0).');
                M = getM_from_T_T0(inputValue);
                T_T0_val = inputValue;
                P_P0_val = isentropicP_P0(M);
                Rho_Rho0_val = isentropicRho_Rho0(M);
                A_Astar_val = isentropicA_Astar(M);
                break;
            case 'P_P0':
                if (inputValue <= 0 || inputValue > 1) throw new Error('P/P₀ must be between 0 and 1 (exclusive of 0).');
                M = getM_from_P_P0(inputValue);
                T_T0_val = isentropicT_T0(M);
                P_P0_val = inputValue;
                Rho_Rho0_val = isentropicRho_Rho0(M);
                A_Astar_val = isentropicA_Astar(M);
                break;
            case 'Rho_Rho0':
                if (inputValue <= 0 || inputValue > 1) throw new Error('ρ/ρ₀ must be between 0 and 1 (exclusive of 0).');
                M = getM_from_Rho_Rho0(inputValue);
                T_T0_val = isentropicT_T0(M);
                P_P0_val = isentropicP_P0(M);
                Rho_Rho0_val = isentropicRho_Rho0(M);
                A_Astar_val = isentropicA_Astar(M);
                break;
            case 'A_Astar':
                if (inputValue < 1) throw new Error('A/A* must be 1 or greater.');
                const flowType = isentropicFlowType.value;
                M = getM_from_A_Astar(inputValue, flowType);
                A_Astar_val = inputValue;
                T_T0_val = isentropicT_T0(M);
                P_P0_val = isentropicP_P0(M);
                Rho_Rho0_val = isentropicRho_Rho0(M);
                break;
        }

        if (isNaN(M)) throw new Error('Could not calculate Mach number from input. Check values.');

        document.getElementById('isentropic-M').textContent = formatNumber(M);
        document.getElementById('isentropic-T_T0').textContent = formatNumber(T_T0_val);
        document.getElementById('isentropic-P_P0').textContent = formatNumber(P_P0_val);
        document.getElementById('isentropic-Rho_Rho0').textContent = formatNumber(Rho_Rho0_val);
        document.getElementById('isentropic-A_Astar').textContent = formatNumber(A_Astar_val);

    } catch (e) {
        displayError('isentropic-error', e.message);
        document.getElementById('isentropic-M').textContent = 'N/A';
        document.getElementById('isentropic-T_T0').textContent = 'N/A';
        document.getElementById('isentropic-P_P0').textContent = 'N/A';
        document.getElementById('isentropic-Rho_Rho0').textContent = 'N/A';
        document.getElementById('isentropic-A_Astar').textContent = 'N/A';
    }
}


// --- Normal Shock Calculator ---
function calculateNormalShock() {
    hideError('normal-shock-error');
    const M1 = parseFloat(document.getElementById('M1-normal-shock').value);

    if (isNaN(M1) || M1 < 1) {
        displayError('normal-shock-error', 'Upstream Mach number (M₁) must be > 1.');
        document.getElementById('M2-normal-shock').textContent = 'N/A';
        document.getElementById('P2-P1-normal-shock').textContent = 'N/A';
        document.getElementById('T2-T1-normal-shock').textContent = 'N/A';
        document.getElementById('Rho2-Rho1-normal-shock').textContent = 'N/A';
        document.getElementById('P02-P01-normal-shock').textContent = 'N/A';
        return;
    }

    const M2 = normalShockM2(M1);
    const P2_P1 = normalShockP2_P1(M1);
    const T2_T1 = normalShockT2_T1(M1);
    const Rho2_Rho1 = normalShockRho2_Rho1(M1);
    const P02_P01 = normalShockP02_P01(M1);

    if (isNaN(M2)) {
        displayError('normal-shock-error', 'Calculation error. Check M₁ value.');
    }

    document.getElementById('M2-normal-shock').textContent = formatNumber(M2);
    document.getElementById('P2-P1-normal-shock').textContent = formatNumber(P2_P1);
    document.getElementById('T2-T1-normal-shock').textContent = formatNumber(T2_T1);
    document.getElementById('Rho2-Rho1-normal-shock').textContent = formatNumber(Rho2_Rho1);
    document.getElementById('P02-P01-normal-shock').textContent = formatNumber(P02_P01);
}

// --- Oblique Shock Calculator ---
function calculateObliqueShock() {
    hideError('oblique-shock-error');
    const M1 = parseFloat(document.getElementById('M1-oblique-shock').value);
    const theta_deg = parseFloat(document.getElementById('theta-oblique-shock').value);
    const shockType = document.getElementById('oblique-shock-type').value;

    document.getElementById('max-theta-oblique-shock').textContent = formatNumber(calculateMaxTheta(M1));


    if (isNaN(M1) || M1 < 1) {
        displayError('oblique-shock-error', 'Upstream Mach number (M₁) must be > 1.');
        // Clear all results
        Object.keys(calculateObliqueShockProperties(2, 45) || {}).forEach(key => {
            if (key !== 'theta_deg' && key !== 'M1n' && key !== 'M2n') { // Don't clear intermediate
                document.getElementById(key.replace(/([A-Z])/g, '-$1').toLowerCase() + '-oblique-shock')?.textContent = 'N/A';
            }
        });
        document.getElementById('beta-oblique-shock').textContent = 'N/A';
        return;
    }

    if (isNaN(theta_deg) || theta_deg < 0) {
        displayError('oblique-shock-error', 'Deflection angle (θ) must be a non-negative number.');
        Object.keys(calculateObliqueShockProperties(2, 45) || {}).forEach(key => {
            if (key !== 'theta_deg' && key !== 'M1n' && key !== 'M2n') { // Don't clear intermediate
                document.getElementById(key.replace(/([A-Z])/g, '-$1').toLowerCase() + '-oblique-shock')?.textContent = 'N/A';
            }
        });
        document.getElementById('beta-oblique-shock').textContent = 'N/A';
        return;
    }

    const beta_deg = calculateBeta(M1, theta_deg, shockType);

    if (isNaN(beta_deg)) {
        if (theta_deg > calculateMaxTheta(M1)) {
            displayError('oblique-shock-error', `Shock is detached for M₁=${M1} and θ=${theta_deg}°. Max θ is ${formatNumber(calculateMaxTheta(M1))}°.`);
        } else {
            displayError('oblique-shock-error', 'Could not calculate Beta. Check M₁ and θ.');
        }
        document.getElementById('beta-oblique-shock').textContent = 'N/A';
        Object.keys(calculateObliqueShockProperties(2, 45) || {}).forEach(key => {
            if (key !== 'theta_deg' && key !== 'M1n' && key !== 'M2n') {
                document.getElementById(key.replace(/([A-Z])/g, '-$1').toLowerCase() + '-oblique-shock')?.textContent = 'N/A';
            }
        });
        return;
    }

    const results = calculateObliqueShockProperties(M1, beta_deg);

    if (!results) {
        displayError('oblique-shock-error', 'Calculation error. Beta value might be invalid or no shock formed.');
        document.getElementById('beta-oblique-shock').textContent = 'N/A';
        Object.keys(calculateObliqueShockProperties(2, 45) || {}).forEach(key => {
            if (key !== 'theta_deg' && key !== 'M1n' && key !== 'M2n') {
                document.getElementById(key.replace(/([A-Z])/g, '-$1').toLowerCase() + '-oblique-shock')?.textContent = 'N/A';
            }
        });
        return;
    }

    document.getElementById('beta-oblique-shock').textContent = formatNumber(beta_deg);
    document.getElementById('M2-oblique-shock').textContent = formatNumber(results.M2);
    document.getElementById('P2-P1-oblique-shock').textContent = formatNumber(results.P2_P1);
    document.getElementById('T2-T1-oblique-shock').textContent = formatNumber(results.T2_T1);
    document.getElementById('Rho2-Rho1-oblique-shock').textContent = formatNumber(results.Rho2_Rho1);
    document.getElementById('P02-P01-oblique-shock').textContent = formatNumber(results.P02_P01);
}

// --- Prandtl-Meyer Expansion Calculator ---
const prandtlMeyerInputType = document.getElementById('prandtl-meyer-input-type');
const prandtlMeyerInputValue = document.getElementById('prandtl-meyer-input-value');
const prandtlMeyerInputLabel = document.getElementById('prandtl-meyer-input-label');

prandtlMeyerInputType.addEventListener('change', () => {
    const selectedType = prandtlMeyerInputType.value;
    if (selectedType === 'M') {
        prandtlMeyerInputLabel.textContent = 'Mach Number (M)';
        prandtlMeyerInputValue.min = '1.0001';
        prandtlMeyerInputValue.value = '2.0';
    } else {
        prandtlMeyerInputLabel.textContent = 'ν [degrees] (ν ≥ 0)';
        prandtlMeyerInputValue.min = '0';
        prandtlMeyerInputValue.value = '20';
    }
    hideError('prandtl-meyer-error');
});

function calculatePrandtlMeyer() {
    hideError('prandtl-meyer-error');
    const inputType = prandtlMeyerInputType.value;
    const inputValue = parseFloat(prandtlMeyerInputValue.value);
    let M_val, nu_val;

    if (isNaN(inputValue)) {
        displayError('prandtl-meyer-error', 'Please enter a valid number.');
        return;
    }

    try {
        if (inputType === 'M') {
            if (inputValue < 1) throw new Error('Mach number (M) must be > 1 for Prandtl-Meyer expansion.');
            M_val = inputValue;
            nu_val = prandtlMeyer_nu(M_val);
        } else { // inputType === 'nu'
            if (inputValue < 0) throw new Error('Prandtl-Meyer function (ν) must be non-negative.');
            nu_val = inputValue;
            M_val = getM_from_nu(nu_val);
        }

        if (isNaN(M_val) || isNaN(nu_val)) throw new Error('Calculation failed. Check input values.');

        document.getElementById('prandtl-meyer-M').textContent = formatNumber(M_val);
        document.getElementById('prandtl-meyer-nu').textContent = formatNumber(nu_val);

    } catch (e) {
        displayError('prandtl-meyer-error', e.message);
        document.getElementById('prandtl-meyer-M').textContent = 'N/A';
        document.getElementById('prandtl-meyer-nu').textContent = 'N/A';
    }
}

// Initial calculations on page load
document.addEventListener('DOMContentLoaded', () => {
    calculateIsentropic();
    calculateNormalShock();
    calculateObliqueShock();
    calculatePrandtlMeyer();
});