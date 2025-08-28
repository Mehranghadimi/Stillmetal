/* ===============================
   Aluminothermic Reduction – app.js
   Split from the original inline <script>
   =============================== */

// ---------- Config ----------
const PY_FILE = 'main.py'; // your Python analyzer (Pyodide will run this)

// ---------- Tiny helpers ----------
const $ = (sel) => document.querySelector(sel);
const round = (x, d = 2) => Number.parseFloat(x).toFixed(d);
const clone = (obj) => JSON.parse(JSON.stringify(obj));

function enhanceNumberInputs() {
  document.querySelectorAll('input[type="number"]').forEach((inp) => {
    inp.addEventListener('focus', () => {
      try {
        inp.select();
      } catch (e) {}
    });
    inp.addEventListener('blur', () => {
      if ((inp.value || '').trim() === '') inp.value = '0';
    });
  });
}

// ---------- Constants ----------
const MOLAR_MASSES = {
  CaO: 56.077,
  SiO2: 60.084,
  FeO: 71.844,
  Fe2O3: 159.69,
  Fe3O4: 231.533,
  MgO: 40.304,
  Al2O3: 101.961,
  P2O5: 141.944,
  Cr2O3: 151.99,
  MnO: 70.937,
  TiO2: 79.866,
  Al: 26.9815,
  Si: 28.0855,
  Fe: 55.845,
  P: 30.9738,
  Cr: 51.9961,
  Mn: 54.938,
  Ti: 47.867,
};

const OXIDE_GROUPS = {
  Main: ['CaO', 'MgO', 'Al2O3', 'SiO2'],
  Iron: ['FeO', 'Fe2O3', 'Fe3O4'],
  Other: ['P2O5', 'Cr2O3', 'MnO', 'TiO2'],
};

const REDUCTION_PRIORITY = ['FeO', 'Fe2O3', 'Fe3O4', 'SiO2', 'P2O5', 'Cr2O3', 'MnO', 'TiO2'];

const REDUCTION_STOICHIOMETRY = {
  SiO2: { Al_coeff: 4 / 3, metal_coeff: 1, metal_key: 'Si', Al2O3_coeff: 2 / 3 },
  FeO: { Al_coeff: 2 / 3, metal_coeff: 1, metal_key: 'Fe', Al2O3_coeff: 1 / 3 },
  Fe2O3: { Al_coeff: 2, metal_coeff: 2, metal_key: 'Fe', Al2O3_coeff: 1 },
  Fe3O4: { Al_coeff: 8 / 3, metal_coeff: 3, metal_key: 'Fe', Al2O3_coeff: 4 / 3 },
  P2O5: { Al_coeff: 10 / 3, metal_coeff: 2, metal_key: 'P', Al2O3_coeff: 5 / 3 },
  Cr2O3: { Al_coeff: 2, metal_coeff: 2, metal_key: 'Cr', Al2O3_coeff: 1 },
  TiO2: { Al_coeff: 4 / 3, metal_coeff: 1, metal_key: 'Ti', Al2O3_coeff: 2 / 3 },
  MnO: { Al_coeff: 2 / 3, metal_coeff: 1, metal_key: 'Mn', Al2O3_coeff: 1 / 3 },
};

// Headers for oxide–temperature table (Temperature at the END; Viscosity last & editable)
const CROSS_HEADERS = [
  'C/A',
  'SiO2 [g]',
  'Al2O3 [g]',
  'CaO [g]',
  'MgO [g]',
  'MnO [g]',
  'Temperature [°C]',
  'Viscosity',
];
// Map table header label -> composition key
const HEADER_TO_KEY = {
  'SiO2 [g]': 'SiO2',
  'Al2O3 [g]': 'Al2O3',
  'CaO [g]': 'CaO',
  'MgO [g]': 'MgO',
};

// ---------- State ----------
let initialRes = null;
let initialTotals = null;
let lastInputs = null;
let lastCARatios = [];
let lastCompositions = []; // post-reduction compositions aligned with lastCARatios

// ---------- Core calculations ----------
function performReductionCalculations(oxideQuantities, al2o3Limit = null) {
  let currentAl = 0,
    totalSi = 0,
    totalFe = 0,
    totalAl2O3 = 0,
    totalP = 0,
    totalCr = 0,
    totalMn = 0,
    totalTi = 0;
  const remaining = clone(oxideQuantities);

  for (const oxide of REDUCTION_PRIORITY) {
    const amt = remaining[oxide] || 0;
    if (amt <= 0 || !(oxide in REDUCTION_STOICHIOMETRY)) continue;

    const molesOx = amt / MOLAR_MASSES[oxide];
    const s = REDUCTION_STOICHIOMETRY[oxide];
    const alNeededFull = molesOx * s.Al_coeff * MOLAR_MASSES.Al;
    const al2o3FromFull = molesOx * s.Al2O3_coeff * MOLAR_MASSES.Al2O3;

    if (al2o3Limit !== null && totalAl2O3 + al2o3FromFull > al2o3Limit) {
      const al2o3Remain = al2o3Limit - totalAl2O3;
      if (al2o3Remain <= 0) break;

      const molesReduce = Math.min(al2o3Remain / MOLAR_MASSES.Al2O3 / s.Al2O3_coeff, molesOx);
      const alCons = molesReduce * s.Al_coeff * MOLAR_MASSES.Al;
      const oxReduced = molesReduce * MOLAR_MASSES[oxide];
      const al2o3Formed = molesReduce * s.Al2O3_coeff * MOLAR_MASSES.Al2O3;

      currentAl += alCons;
      totalAl2O3 += al2o3Formed;

      const metalMoles = molesReduce * s.metal_coeff;
      const metalMass = metalMoles * MOLAR_MASSES[s.metal_key];

      if (s.metal_key === 'Si') totalSi += metalMass;
      else if (s.metal_key === 'Fe') totalFe += metalMass;
      else if (s.metal_key === 'P') totalP += metalMass;
      else if (s.metal_key === 'Cr') totalCr += metalMass;
      else if (s.metal_key === 'Mn') totalMn += metalMass;
      else if (s.metal_key === 'Ti') totalTi += metalMass;

      remaining[oxide] -= oxReduced;
      break;
    } else {
      currentAl += alNeededFull;
      totalAl2O3 += al2o3FromFull;

      const metalMoles = molesOx * s.metal_coeff;
      const metalMass = metalMoles * MOLAR_MASSES[s.metal_key];
      if (s.metal_key === 'Si') totalSi += metalMass;
      else if (s.metal_key === 'Fe') totalFe += metalMass;
      else if (s.metal_key === 'P') totalP += metalMass;
      else if (s.metal_key === 'Cr') totalCr += metalMass;
      else if (s.metal_key === 'Mn') totalMn += metalMass;
      else if (s.metal_key === 'Ti') totalTi += metalMass;

      remaining[oxide] = 0;
    }
  }

  const unreduced = {};
  Object.entries(remaining).forEach(([k, v]) => {
    if (v > 0.001) unreduced[k] = v;
  });

  return {
    total_al_needed_g: currentAl,
    total_si_formed_g: totalSi,
    total_fe_formed_g: totalFe,
    total_al2o3_formed_g: totalAl2O3,
    total_p_formed_g: totalP,
    total_cr_formed_g: totalCr,
    total_mn_formed_g: totalMn,
    total_ti_formed_g: totalTi,
    unreduced_oxides: unreduced,
  };
}

function calculateAndFormatRow(oxideQ, desiredCA, initialTotalAl2O3, initialRes) {
  const CaO = oxideQ.CaO || 0;
  const initAl2O3 = oxideQ.Al2O3 || 0;
  const requiredTotalAl2O3 = CaO / desiredCA;

  let alNeeded = 0;
  let leftover = 'None';
  let leftoverBasicity = 'N/A';
  let fesi25 = 'N/A';
  let extraFe = 'N/A';
  let totalSi = 0;
  let totalFe = 0;
  let extraSiO2Added = 0;

  const initRatio = initialTotalAl2O3 > 0 ? CaO / initialTotalAl2O3 : 0;
  let postComp = { CaO: oxideQ.CaO || 0, MgO: oxideQ.MgO || 0 };

  if (desiredCA > initRatio) {
    // need less alumina -> partial reduction up to target
    let currentAl = 0,
      currentAl2O3 = 0,
      curSi = 0,
      curFe = 0;
    const remaining = clone(oxideQ);

    const priority = OXIDE_GROUPS.Iron.filter((ox) => ox in remaining).concat(['SiO2', 'P2O5', 'Cr2O3', 'MnO', 'TiO2']);
    const S = {
      SiO2: { Al_coeff: 4 / 3, metal_coeff: 1, metal_key: 'Si', Al2O3_coeff: 2 / 3 },
      FeO: { Al_coeff: 2 / 3, metal_coeff: 1, metal_key: 'Fe', Al2O3_coeff: 1 / 3 },
      Fe2O3: { Al_coeff: 2, metal_coeff: 2, metal_key: 'Fe', Al2O3_coeff: 1 },
      Fe3O4: { Al_coeff: 8 / 3, metal_coeff: 3, metal_key: 'Fe', Al2O3_coeff: 4 / 3 },
      P2O5: { Al_coeff: 10 / 3, metal_coeff: 2, metal_key: 'P', Al2O3_coeff: 5 / 3 },
      Cr2O3: { Al_coeff: 2, metal_coeff: 2, metal_key: 'Cr', Al2O3_coeff: 1 },
      MnO: { Al_coeff: 2 / 3, metal_coeff: 1, metal_key: 'Mn', Al2O3_coeff: 1 / 3 },
      TiO2: { Al_coeff: 4 / 3, metal_coeff: 1, metal_key: 'Ti', Al2O3_coeff: 2 / 3 },
    };

    const targetAl2O3FromReduction = Math.max(0, requiredTotalAl2O3 - initAl2O3);
    for (const ox of priority) {
      if ((remaining[ox] || 0) <= 0) continue;

      const molesOx = remaining[ox] / MOLAR_MASSES[ox];
      const s = S[ox];
      const al2o3Full = molesOx * s.Al2O3_coeff * MOLAR_MASSES.Al2O3;
      const alFull = molesOx * s.Al_coeff * MOLAR_MASSES.Al;

      if (currentAl2O3 + al2o3Full > targetAl2O3FromReduction) {
        const al2o3Remain = targetAl2O3FromReduction - currentAl2O3;
        if (al2o3Remain <= 0) break;

        const molesReduce = Math.min(al2o3Remain / MOLAR_MASSES.Al2O3 / s.Al2O3_coeff, molesOx);
        const alCons = molesReduce * s.Al_coeff * MOLAR_MASSES.Al;
        const metalMoles = molesReduce * s.metal_coeff;
        const metalMass = metalMoles * MOLAR_MASSES[s.metal_key];

        currentAl += alCons;
        currentAl2O3 += molesReduce * s.Al2O3_coeff * MOLAR_MASSES.Al2O3;

        if (s.metal_key === 'Si') curSi += metalMass;
        else if (s.metal_key === 'Fe') curFe += metalMass;

        remaining[ox] -= molesReduce * MOLAR_MASSES[ox];
        break;
      } else {
        currentAl += alFull;
        currentAl2O3 += al2o3Full;

        const metalMoles = molesOx * s.metal_coeff;
        const metalMass = metalMoles * MOLAR_MASSES[s.metal_key];

        if (s.metal_key === 'Si') curSi += metalMass;
        else if (s.metal_key === 'Fe') curFe += metalMass;

        remaining[ox] = 0;
      }
    }

    alNeeded = currentAl;
    totalSi = curSi;
    totalFe = curFe;
    extraSiO2Added = 0;

    const ironKey = OXIDE_GROUPS.Iron.find((k) => (oxideQ[k] || 0) > 0);
    const parts = [];
    if (ironKey && (remaining[ironKey] || 0) > 0.001) parts.push(`${ironKey}: ${round(remaining[ironKey])} g`);
    if ((remaining.SiO2 || 0) > 0.001) parts.push(`SiO2: ${round(remaining.SiO2)} g`);
    leftover = parts.length ? parts.join(', ') : 'None';

    if ((remaining.SiO2 || 0) > 0) leftoverBasicity = round((oxideQ.CaO || 0) / remaining.SiO2, 3);

    postComp = { ...postComp, ...remaining, Al2O3: initAl2O3 + currentAl2O3 };
  } else {
    // need more alumina -> full reduction + stoichiometric add of SiO2+Al to make up deficiency
    const full = performReductionCalculations(oxideQ, null);
    const totalAl2O3After = (oxideQ.Al2O3 || 0) + full.total_al2o3_formed_g;
    const deficiency = requiredTotalAl2O3 - totalAl2O3After;

    let addAl = 0,
      addSi = 0,
      addSiO2 = 0;

    if (deficiency > 0) {
      // From stoichiometry: Al2O3 formation per SiO2+Al reaction
      const molesAl2O3 = deficiency / MOLAR_MASSES.Al2O3;
      const molesAl = molesAl2O3 * (4 / 2);
      const molesSi = molesAl2O3 * (3 / 2);
      const molesSiO2 = molesAl2O3 * (3 / 2);

      addAl = molesAl * MOLAR_MASSES.Al;
      addSi = molesSi * MOLAR_MASSES.Si;
      addSiO2 = molesSiO2 * MOLAR_MASSES.SiO2; // added SiO2 fully reacts
    }

    alNeeded = full.total_al_needed_g + addAl;
    totalSi = full.total_si_formed_g + addSi;
    totalFe = full.total_fe_formed_g;
    extraSiO2Added = addSiO2;
    leftover = 'None';
    leftoverBasicity = 'N/A';

    postComp = { ...postComp, ...(full.unreduced_oxides || {}), Al2O3: Math.max(requiredTotalAl2O3, 0) };
    ['SiO2', 'FeO', 'Fe2O3', 'Fe3O4', 'P2O5', 'Cr2O3', 'MnO', 'TiO2'].forEach((k) => {
      if (!(k in postComp)) postComp[k] = 0;
    });
  }

  if (totalSi > 0) {
    const totalFeSi25 = totalSi / 0.25;
    const feRequired = totalFeSi25 * 0.75;
    const addFe = Math.max(0, feRequired - totalFe);
    fesi25 = round(totalFeSi25);
    extraFe = round(addFe);
  }

  const alloyPercents = { FeSi30: 30, FeSi35: 35, FeSi40: 40, FeSi45: 45 };
  const alloy_results = {};
  if (totalSi > 0) {
    Object.entries(alloyPercents).forEach(([name, pct]) => {
      const fraction = pct / 100;
      const totalAlloy = totalSi / fraction;
      const requiredIron = totalAlloy - totalSi;
      const extraIron = Math.max(0, requiredIron - totalFe);
      alloy_results[name] = extraIron.toFixed(2);
    });
  } else {
    Object.keys(alloyPercents).forEach((name) => (alloy_results[name] = 'N/A'));
  }

  return {
    ca_ratio: round(desiredCA, 2),
    al_needed: round(alNeeded),
    extra_silica: round(extraSiO2Added),
    leftover_slag: leftover,
    basicity: leftoverBasicity,
    fesi25,
    extra_fe: extraFe,
    alloy_results,
    post_comp: postComp,
  };
}

// ---------- UI glue ----------
function renderInitialStats(q, res) {
  const SiO2 = q.SiO2 || 0;
  const CaO = q.CaO || 0;
  const cs = SiO2 > 0 ? CaO / SiO2 : null;

  const initAl2O3 = q.Al2O3 || 0;
  const formedAl2O3 = res.total_al2o3_formed_g || 0;
  const totalAl2O3 = initAl2O3 + formedAl2O3;

  const si = res.total_si_formed_g || 0;
  const fe = res.total_fe_formed_g || 0;

  let totalFeSi25 = null,
    extraFe25 = null;
  if (si > 0) {
    totalFeSi25 = si / 0.25;
    const feRequired = totalFeSi25 * 0.75;
    extraFe25 = Math.max(0, feRequired - fe);
  }

  const caAlTotal = totalAl2O3 > 0 && CaO > 0 ? CaO / totalAl2O3 : null;

  $('#statCS').textContent = cs === null ? 'n/a' : round(cs, 3);
  $('#statAl').textContent = round(res.total_al_needed_g);
  $('#statSi').textContent = round(si);
  $('#statFe').textContent = round(fe);
  $('#statFeSi25').textContent = totalFeSi25 === null ? 'n/a' : round(totalFeSi25);
  $('#statExtraFe25').textContent = extraFe25 === null ? 'n/a' : round(extraFe25);
  $('#statAl2O3').textContent = round(formedAl2O3);
  $('#statAl2O3Init').textContent = round(initAl2O3);
  $('#statAl2O3Total').textContent = round(totalAl2O3);
  $('#statCATotal').textContent = caAlTotal === null ? 'n/a' : round(caAlTotal, 3);

  const all = [...OXIDE_GROUPS.Main, ...OXIDE_GROUPS.Iron, ...OXIDE_GROUPS.Other];
  const lines = all.map((k) => `${k}: ${round(q[k] || 0)} g`).join('\n');
  $('#rawList').textContent = lines;

  $('#initialStats').style.display = 'grid';
}

function sweepAndRender() {
  if (!lastInputs || !initialTotals) {
    alert('First click "Run initial calculation".');
    return;
  }

  const start = parseFloat($('#startRatio').value);
  const end = parseFloat($('#endRatio').value);
  const step = parseFloat($('#stepRatio').value);

  if (!(start > 0 && end > 0 && step > 0 && start <= end)) {
    alert('Invalid C/A sweep inputs.');
    return;
  }

  const tbody = $('#resultsTable tbody');
  tbody.innerHTML = '';
  const tbodyAlloy = $('#alloyTable tbody');
  tbodyAlloy.innerHTML = '';

  lastCARatios = [];
  lastCompositions = [];

  for (let r = start; r <= end + step / 2; r = r + step) {
    const row = calculateAndFormatRow(lastInputs.q, r, initialTotals.totalAl2O3, initialRes);
    lastCARatios.push(Number(row.ca_ratio));
    lastCompositions.push(row.post_comp || {});

    const tr = document.createElement('tr');
    [row.ca_ratio, row.al_needed, row.extra_silica, row.leftover_slag, row.basicity, row.fesi25, row.extra_fe].forEach(
      (t) => {
        const td = document.createElement('td');
        td.textContent = t;
        tr.appendChild(td);
      }
    );
    tbody.appendChild(tr);

    const tr2 = document.createElement('tr');
    const vals = [
      row.ca_ratio,
      row.alloy_results.FeSi30,
      row.alloy_results.FeSi35,
      row.alloy_results.FeSi40,
      row.alloy_results.FeSi45,
    ];
    vals.forEach((t) => {
      const td = document.createElement('td');
      td.textContent = t;
      tr2.appendChild(td);
    });
    tbodyAlloy.appendChild(tr2);
  }

  $('#resultsContainer').style.display = 'block';
  $('#resultsCard').style.display = 'block';
  $('#alloyContainer').style.display = 'block';
  $('#alloyCard').style.display = 'block';

  // Enable cross table + CA charge sections
  $('#genCross').disabled = false;
  $('#downloadCross').disabled = true; // will enable after cross is rendered
  $('#crossCard').style.display = 'block';

  const caBtn = $('#genCACharges');
  if (caBtn) caBtn.disabled = false;
  const caSec = $('#caChargeSection');
  if (caSec) caSec.style.display = 'block';
}

function toCSV() {
  const rows = [];
  rows.push(
    [
      'C/A Ratio',
      'Al Needed (g)',
      'Extra SiO2 Added (g)',
      'Leftover Slag (g)',
      'Leftover Basicity (C/S)',
      'FeSi25 Result (g)',
      'Extra Fe Needed (g)',
    ].join(',')
  );
  document.querySelectorAll('#resultsTable tbody tr').forEach((tr) => {
    const cells = Array.from(tr.children).map((td) => '"' + (td.textContent || '').replaceAll('"', '""') + '"');
    rows.push(cells.join(','));
  });
  const csv = rows.join('\n');
  const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'aluminothermic_results.csv';
  a.click();
  URL.revokeObjectURL(url);
}

// ---------- Oxide–Temperature table ----------
function parseTemperatureInput(str) {
  if (!str) return [];
  const parts = str
    .split(',')
    .map((s) => s.trim())
    .filter(Boolean);

  if (parts.length === 1) {
    const v = Number(parts[0]);
    return Number.isFinite(v) ? [v] : [];
  }
  if (parts.length === 3) {
    const start = Number(parts[0]),
      end = Number(parts[1]),
      step = Number(parts[2]);
    if ([start, end, step].every(Number.isFinite) && step > 0 && end >= start) {
      const out = [];
      for (let t = start; t <= end + 1e-9; t += step) out.push(Math.round(t * 1000) / 1000);
      return out;
    }
  }
  const out = [];
  let ok = true;
  parts.forEach((p) => {
    const v = Number(p);
    if (Number.isFinite(v)) out.push(v);
    else ok = false;
  });
  return ok ? out : [];
}

function renderCrossHeader() {
  const tr = $('#crossHeader');
  tr.innerHTML = '';
  CROSS_HEADERS.forEach((h) => {
    const th = document.createElement('th');
    th.textContent = h;
    tr.appendChild(th);
  });
}

function wireViscosityPaste() {
  const inputs = Array.from(document.querySelectorAll('#crossTable tbody input.visc-input'));
  inputs.forEach((inp, startIndex) => {
    inp.addEventListener('paste', (e) => {
      const clip = (e.clipboardData && e.clipboardData.getData('text')) || '';
      if (!clip) return;
      e.preventDefault();
      const values = parsePastedColumn(clip);
      for (let i = 0; i < values.length && startIndex + i < inputs.length; i++) {
        const v = values[i];
        inputs[startIndex + i].value = v === '' ? '' : String(v);
      }
    });
  });
}

function parsePastedColumn(text) {
  if (!text) return [];
  const rows = text.replace(/\r/g, '').split('\n').filter(Boolean);
  const vals = [];
  for (const row of rows) {
    const first = row.split(/\t|,|;|\s+/)[0];
    const v = parseFloat(first);
    vals.push(Number.isFinite(v) ? v : '');
  }
  return vals;
}

function renderCrossTable() {
  const tempStr = $('#tempInput').value.trim();
  const temps = parseTemperatureInput(tempStr);
  if (!temps.length) {
    alert(
      'Enter temperature as a single value (e.g., 1200), a range "start, end, step" (e.g., 1100, 1400, 100), or a comma list (e.g., 1150, 1250, 1350).'
    );
    return;
  }
  if (!lastCARatios.length) {
    alert('Generate the C/A table first.');
    return;
  }

  renderCrossHeader();
  const tbody = document.querySelector('#crossTable tbody');
  tbody.innerHTML = '';

  const q0 = lastInputs ? lastInputs.q : {};
  const tempColIndex = CROSS_HEADERS.indexOf('Temperature [°C]');
  const dataLabels = CROSS_HEADERS.slice(1, tempColIndex);

  lastCARatios.forEach((ca, idx) => {
    const comp = lastCompositions[idx] || {};
    temps.forEach((T) => {
      const tr = document.createElement('tr');

      // C/A
      let cells = [String(Number(ca).toFixed(2))];

      // oxide columns using post-reduction comp; fallback to original if missing
      dataLabels.forEach((label) => {
        const key = HEADER_TO_KEY[label];
        let val = 0;
        if (key) {
          if (key in comp) val = comp[key];
          else if (key in q0) val = q0[key] || 0;
        }
        cells.push(String(round(val)));
      });

      // Temperature value
      cells.push(String(T));

      // Append fixed cells
      cells.forEach((v) => {
        const td = document.createElement('td');
        td.textContent = v;
        tr.appendChild(td);
      });

      // Viscosity input (editable)
      const tdInput = document.createElement('td');
      const input = document.createElement('input');
      input.type = 'number';
      input.step = 'any';
      input.placeholder = '';
      input.className = 'visc-input';
      input.style.width = '100%';
      input.style.boxSizing = 'border-box';
      input.style.padding = '8px 10px';
      input.style.border = '1px solid #2a386d';
      input.style.borderRadius = '8px';
      input.style.background = '#0c1230';
      input.style.color = '#e9f0ff';
      tdInput.appendChild(input);
      tr.appendChild(tdInput);

      tbody.appendChild(tr);
    });
  });

  $('#crossContainer').style.display = 'block';
  $('#crossCard').style.display = 'block';
  $('#downloadCross').disabled = false;

  wireViscosityPaste();
}

function downloadCrossCSV() {
  const rows = [];
  rows.push(CROSS_HEADERS.join(','));

  document.querySelectorAll('#crossTable tbody tr').forEach((tr) => {
    const cells = Array.from(tr.children).map((td, idx) => {
      if (idx === tr.children.length - 1) {
        // viscosity cell
        const inp = td.querySelector('input');
        const v = inp ? inp.value : '';
        return `"${(v || '').replaceAll('"', '""')}"`;
      }
      return `"${(td.textContent || '').replaceAll('"', '""')}"`;
    });
    rows.push(cells.join(','));
  });

  const csv = rows.join('\n');
  const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'oxide_temperature_table.csv';
  a.click();
  URL.revokeObjectURL(url);
}

function copyCrossToClipboard() {
  // Copy numeric values only, skip C/A column, no header
  const lines = [];
  const rows = document.querySelectorAll('#crossTable tbody tr');

  rows.forEach((tr) => {
    const tds = Array.from(tr.children).slice(1); // drop first col (C/A)
    const vals = tds.map((td, i, arr) => {
      if (i === arr.length - 1) {
        // viscosity input
        const inp = td.querySelector('input');
        const num = parseFloat(inp ? inp.value : '');
        return Number.isFinite(num) ? String(num) : '';
      }
      const num = parseFloat(td.textContent || '');
      return Number.isFinite(num) ? String(num) : '';
    });
    lines.push(vals.join('\t'));
  });

  const text = lines.join('\n');
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard
      .writeText(text)
      .then(() => alert('Copied numeric table values (without C/A).'))
      .catch(() => fallbackCopy(text));
  } else {
    fallbackCopy(text);
  }
}

function fallbackCopy(text) {
  const ta = document.createElement('textarea');
  ta.value = text;
  document.body.appendChild(ta);
  ta.select();
  try {
    document.execCommand('copy');
    alert('Copied table (without C/A) to clipboard.');
  } catch (e) {
    alert('Copy failed.');
  } finally {
    document.body.removeChild(ta);
  }
}

// ---------- Per-C/A charge tables ----------
function parseCARange(str) {
  if (!str) return [];
  const parts = str
    .split(',')
    .map((s) => s.trim())
    .filter(Boolean);

  if (parts.length === 1) {
    const v = Number(parts[0]);
    return Number.isFinite(v) ? [v] : [];
  }
  if (parts.length === 3) {
    const s = Number(parts[0]),
      e = Number(parts[1]),
      st = Number(parts[2]);
    if ([s, e, st].every(Number.isFinite) && st > 0 && e >= s) {
      const out = [];
      for (let x = s; x <= e + 1e-9; x += st) out.push(Math.round(x * 1000) / 1000);
      return out;
    }
  }
  const out = [];
  let ok = true;
  parts.forEach((p) => {
    const v = Number(p);
    if (Number.isFinite(v)) out.push(v);
    else ok = false;
  });
  return ok ? out : [];
}

function copyTextLines(lines) {
  const ta = document.createElement('textarea');
  ta.style.position = 'fixed';
  ta.style.opacity = '0';
  ta.value = lines.join('\n');
  document.body.appendChild(ta);
  ta.select();
  try {
    document.execCommand('copy');
    alert('Copied.');
  } catch (e) {
    alert('Copy failed.');
  } finally {
    document.body.removeChild(ta);
  }
}

function generateCAChargeTables() {
  if (!lastInputs || !initialTotals) {
    alert('First click "Run initial calculation".');
    return;
  }
  const caVals = parseCARange($('#caRangeInput').value);
  if (!caVals.length) {
    alert('Invalid C/A input.');
    return;
  }
  const grade = $('#feSiGrade').value; // "25","30","35","40","45"
  const container = $('#caChargesContainer');
  container.innerHTML = '';
  const q0 = lastInputs.q || {};
  const oxList = [...OXIDE_GROUPS.Main, ...OXIDE_GROUPS.Iron, ...OXIDE_GROUPS.Other];

  caVals.forEach((val) => {
    const row = calculateAndFormatRow(q0, Number(val), initialTotals.totalAl2O3, initialRes);
    const items = [];
    const push = (g, sym) => {
      const n = Number(g) || 0;
      if (Math.abs(n) > 1e-9) items.push([n, sym]);
    };

    // initial oxides only if non-zero
    oxList.forEach((k) => {
      const v = Number(q0[k] || 0);
      if (v && v !== 0) push(v, k);
    });

    // extra SiO2 merge
    const extraSiO2 = Math.max(0, Number(row.extra_silica || 0));
    if (extraSiO2 > 0) {
      const idx = items.findIndex(([, s]) => s === 'SiO2');
      if (idx >= 0) items[idx][0] += extraSiO2;
      else push(extraSiO2, 'SiO2');
    }

    // Aluminum to add
    const alToAdd = Math.max(0, Number(row.al_needed || 0));
    if (alToAdd > 0) push(alToAdd, 'Al');

    // Extra Fe for selected FeSi grade
    let feAdd = 0;
    if (grade === '25') feAdd = Number(row.extra_fe || 0);
    else {
      const key = 'FeSi' + grade;
      feAdd = Number((row.alloy_results && row.alloy_results[key]) || 0);
    }
    if (feAdd > 0) push(feAdd, 'Fe');

    // Render: "<grams> <symbol>"
    const lines = items.map(([g, s]) => Number(g).toFixed(2) + ' ' + s);

    const wrap = document.createElement('div');
    wrap.style.marginTop = '12px';
    const label = document.createElement('div');
    label.textContent = 'C/A ' + Number(val).toFixed(2);
    label.style.opacity = '0.85';
    const pre = document.createElement('pre');
    pre.textContent = lines.join('\n');
    const btn = document.createElement('button');
    btn.className = 'btn';
    btn.textContent = '📋 Copy';
    btn.addEventListener('click', () => copyTextLines(lines));

    wrap.appendChild(label);
    wrap.appendChild(pre);
    wrap.appendChild(btn);
    container.appendChild(wrap);
  });

  $('#caChargeSection').style.display = 'block';
  container.style.display = 'block';
}

// ---------- Self tests ----------
function almostEqual(a, b, eps = 1e-3) {
  return Math.abs(a - b) <= eps;
}
function runSelfTests() {
  let pass = true,
    msgs = [];

  let t1 = performReductionCalculations({ Fe2O3: 159.69 }, null);
  pass &= almostEqual(t1.total_al_needed_g, 2 * 26.9815);
  pass &= almostEqual(t1.total_al2o3_formed_g, 101.961);
  pass &= almostEqual(t1.total_fe_formed_g, 2 * 55.845);
  msgs.push('T1 Fe2O3 stoich OK');

  let t2 = performReductionCalculations({ SiO2: 60.084 }, null);
  pass &= almostEqual(t2.total_al_needed_g, (4 / 3) * 26.9815);
  pass &= almostEqual(t2.total_al2o3_formed_g, (2 / 3) * 101.961);
  pass &= almostEqual(t2.total_si_formed_g, 28.0855);
  msgs.push('T2 SiO2 stoich OK');

  let t3 = performReductionCalculations({ FeO: 71.844 }, null);
  pass &= almostEqual(t3.total_al_needed_g, (2 / 3) * 26.9815);
  pass &= almostEqual(t3.total_al2o3_formed_g, (1 / 3) * 101.961);
  pass &= almostEqual(t3.total_fe_formed_g, 55.845);
  msgs.push('T3 FeO stoich OK');

  let t4 = performReductionCalculations({ Fe3O4: 231.533 }, null);
  pass &= almostEqual(t4.total_al_needed_g, (8 / 3) * 26.9815);
  pass &= almostEqual(t4.total_al2o3_formed_g, (4 / 3) * 101.961);
  pass &= almostEqual(t4.total_fe_formed_g, 3 * 55.845);
  msgs.push('T4 Fe3O4 stoich OK');

  // Temperature parsing tests
  const a = parseTemperatureInput('1200');
  pass &= a.length === 1 && a[0] === 1200;
  msgs.push('T5 Temp single OK');

  const b = parseTemperatureInput('1100, 1400, 100');
  pass &= b.length === 4 && b[0] === 1100 && b[b.length - 1] === 1400;
  msgs.push('T6 Temp range OK');

  const c = parseTemperatureInput('1150, 1250, 1350');
  pass &= c.length === 3 && c[1] === 1250;
  msgs.push('T7 Temp list OK');

  // UI constant sanity
  pass &= CROSS_HEADERS[CROSS_HEADERS.length - 1] === 'Viscosity';
  msgs.push('T8 Viscosity header OK');

  alert((pass ? '✅ Tests passed' : '❌ Tests failed') + '\n' + msgs.join('\n'));
}

// ---------- Pyodide (XML → CSV Analyzer) ----------
let pyodideInstance = null;

async function initPyodideAndLoad() {
  const statusEl = $('#xmlStatus');
  try {
    // pyodide.js is loaded via <script src="...pyodide.js"> in index.html
    if (typeof loadPyodide !== 'function') {
      statusEl.textContent = 'Python failed: Pyodide script not loaded.';
      return;
    }
    statusEl.textContent = 'Loading Python...';
    pyodideInstance = await loadPyodide({
      indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/',
    });

    // Load your Python file
    statusEl.textContent = `Loading ${PY_FILE}…`;
    const resp = await fetch(PY_FILE, { cache: 'no-store' });
    if (!resp.ok) {
      statusEl.textContent = `Python failed: could not fetch ${PY_FILE} (${resp.status})`;
      return;
    }
    const pyCode = await resp.text();

    // Execute Python code in Pyodide
    await pyodideInstance.runPythonAsync(pyCode);

    // Probe: ensure analyze_xml_b64 exists
    const hasFunc = await pyodideInstance.runPythonAsync(
      'import builtins\n"analyze_xml_b64" in globals() or hasattr(builtins, "analyze_xml_b64")'
    );
    if (!hasFunc) {
      statusEl.textContent =
        'Python loaded, but analyze_xml_b64() not found. Please ensure main.py defines analyze_xml_b64(b64data).';
      return;
    }

    statusEl.textContent = 'Python ready.';
  } catch (e) {
    console.error(e);
    $('#xmlStatus').textContent = 'Python failed: ' + (e && e.message ? e.message : e);
  }
}

async function handleXMLFile(file) {
  const st = document.getElementById('xmlStatus');
  if (!pyodideInstance) {
    st.textContent = 'Python not loaded yet…';
    return;
  }
  try {
    st.textContent = `Processing ${file.name} …`;

    // 1) Read file as ArrayBuffer (browser side)
    const buf = await file.arrayBuffer();

    // 2) Write to Pyodide FS (no base64, no spreads)
    const xmlPath = 'input.xml';
    pyodideInstance.FS.writeFile(xmlPath, new Uint8Array(buf));

    // 3) Run Python to produce CSV file
    const outFile = 'analysis.csv';
    await pyodideInstance.runPythonAsync(
      `analyze_xml_path_to_file("${xmlPath}", "${outFile}")`
    );

    // 4) Read CSV back from Pyodide FS
    const csvText = pyodideInstance.FS.readFile(outFile, { encoding: 'utf8' });

    // 5) Expose for download + copy
    const blob = new Blob([csvText], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const link = document.getElementById('xmlDownloadLink');
    link.href = url;
    link.download = file.name.replace(/\.xml$/i, '') + '_analysis.csv';

    document.getElementById('xmlDownloadArea').style.display = 'block';
    document.getElementById('xmlCopyBtn').onclick = () => {
      navigator.clipboard.writeText(csvText);
      alert('CSV copied');
    };

    st.textContent = 'Done.';
  } catch (e) {
    console.error(e);
    st.textContent = 'Python failed: ' + (e && e.message ? e.message : e);
  }
}

function wireDropZone() {
  const dz = $('#dropZone');
  const fi = $('#xmlFileInput');
  if (!dz || !fi) return;

  dz.addEventListener('click', () => fi.click());

  dz.addEventListener('dragover', (e) => {
    e.preventDefault();
    dz.style.background = '#222';
  });

  dz.addEventListener('dragleave', () => (dz.style.background = 'transparent'));

  dz.addEventListener('drop', (e) => {
    e.preventDefault();
    dz.style.background = 'transparent';
    if (e.dataTransfer.files.length) handleXMLFile(e.dataTransfer.files[0]);
  });

  fi.addEventListener('change', () => {
    if (fi.files.length) handleXMLFile(fi.files[0]);
  });
}

// ---------- Boot ----------
function onReady() {
  enhanceNumberInputs();

  // Initial calc
  $('#run')?.addEventListener('click', () => {
    const getInputs = () => {
      const q = {};
      OXIDE_GROUPS.Main.forEach((k) => (q[k] = parseFloat($('#' + k).value || '0') || 0));

      const ironKey = $('#ironSelect').value;
      OXIDE_GROUPS.Iron.forEach((k) => (q[k] = 0));
      q[ironKey] = parseFloat($('#IronQty').value || '0') || 0;

      OXIDE_GROUPS.Other.forEach((k) => (q[k] = parseFloat($('#' + k).value || '0') || 0));
      return { q, ironKey };
    };

    lastInputs = getInputs();
    const q = lastInputs.q;
    initialRes = performReductionCalculations(q, null);
    initialTotals = { totalAl2O3: (q.Al2O3 || 0) + initialRes.total_al2o3_formed_g };

    renderInitialStats(q, initialRes);
    $('#sweep').disabled = false;
    $('#downloadCSV').disabled = false;
  });

  // Sweep
  $('#sweep')?.addEventListener('click', sweepAndRender);
  $('#downloadCSV')?.addEventListener('click', toCSV);

  // Cross table
  $('#genCross')?.addEventListener('click', renderCrossTable);
  $('#downloadCross')?.addEventListener('click', downloadCrossCSV);
  $('#copyCross')?.addEventListener('click', copyCrossToClipboard);

  // Sample loader
  $('#loadSample')?.addEventListener('click', () => {
    $('#CaO').value = 45;
    $('#MgO').value = 5;
    $('#Al2O3').value = 5;
    $('#SiO2').value = 25;
    $('#ironSelect').value = 'Fe2O3';
    $('#IronQty').value = 5;
    $('#P2O5').value = 2;
    $('#Cr2O3').value = 3;
    $('#MnO').value = 3;
    $('#TiO2').value = 2;
  });

  // Tests
  $('#runTests')?.addEventListener('click', runSelfTests);

  // Per-C/A charge tables
  $('#genCACharges')?.addEventListener('click', generateCAChargeTables);

  // XML analyzer
  wireDropZone();
  initPyodideAndLoad();
}

if (document.readyState !== 'loading') onReady();
else document.addEventListener('DOMContentLoaded', onReady);
