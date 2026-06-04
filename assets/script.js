const charts = {};

function toggleMenu() {
  document.getElementById('navLinks').classList.toggle('show');
}

function fmt(x, d = 4) {
  if (!isFinite(x)) return 'No calculable';
  const abs = Math.abs(x);
  if (abs !== 0 && (abs < 0.0001 || abs > 1000000)) return Number(x).toExponential(3);
  return Number(x).toFixed(d).replace(/\.0000$/, '').replace(/\.000$/, '').replace(/\.00$/, '');
}

function parseVector(str) {
  return str.split(/[;,\n]+/).map(v => Number(v.trim())).filter(v => !Number.isNaN(v));
}

function parseMatrix(str) {
  return str.trim().split(';').map(row => row.split(',').map(v => Number(v.trim())));
}

function assertSquare(A) {
  if (!A.length || A.some(r => r.length !== A.length)) throw new Error('La matriz debe ser cuadrada.');
}

function zeros(n, m = n) {
  return Array.from({ length: n }, () => Array(m).fill(0));
}

function identity(n) {
  const I = zeros(n);
  for (let i = 0; i < n; i++) I[i][i] = 1;
  return I;
}

function dot(a, b) { return a.reduce((s, v, i) => s + v * b[i], 0); }
function add(a, b) { return a.map((v, i) => v + b[i]); }
function sub(a, b) { return a.map((v, i) => v - b[i]); }
function scale(a, k) { return a.map(v => v * k); }
function norm2(a) { return Math.sqrt(dot(a, a)); }
function normInf(a) { return Math.max(...a.map(v => Math.abs(v))); }
function matVec(A, x) { return A.map(row => dot(row, x)); }
function transpose(A) { return A[0].map((_, j) => A.map(row => row[j])); }
function matMul(A, B) {
  const n = A.length, m = B[0].length, p = B.length;
  const C = zeros(n, m);
  for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) for (let k = 0; k < p; k++) C[i][j] += A[i][k] * B[k][j];
  return C;
}

function residualNorm(A, x, b) { return normInf(sub(matVec(A, x), b)); }

function solveGauss(Ain, bin) {
  const A = Ain.map(r => r.slice());
  const b = bin.slice();
  const n = A.length;
  for (let k = 0; k < n; k++) {
    let maxRow = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[maxRow][k])) maxRow = i;
    if (Math.abs(A[maxRow][k]) < 1e-14) throw new Error('Sistema singular o casi singular.');
    [A[k], A[maxRow]] = [A[maxRow], A[k]];
    [b[k], b[maxRow]] = [b[maxRow], b[k]];
    for (let i = k + 1; i < n; i++) {
      const factor = A[i][k] / A[k][k];
      for (let j = k; j < n; j++) A[i][j] -= factor * A[k][j];
      b[i] -= factor * b[k];
    }
  }
  const x = Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let sum = b[i];
    for (let j = i + 1; j < n; j++) sum -= A[i][j] * x[j];
    x[i] = sum / A[i][i];
  }
  return x;
}

function luSolve(Ain, b) {
  const A = Ain.map(r => r.slice());
  const n = A.length;
  const L = identity(n), U = zeros(n), P = identity(n);
  for (let k = 0; k < n; k++) {
    let pivot = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[pivot][k])) pivot = i;
    if (Math.abs(A[pivot][k]) < 1e-14) throw new Error('No es posible factorizar: pivote nulo.');
    if (pivot !== k) {
      [A[k], A[pivot]] = [A[pivot], A[k]];
      [P[k], P[pivot]] = [P[pivot], P[k]];
      for (let j = 0; j < k; j++) [L[k][j], L[pivot][j]] = [L[pivot][j], L[k][j]];
    }
    for (let j = k; j < n; j++) {
      U[k][j] = A[k][j];
    }
    for (let i = k + 1; i < n; i++) {
      L[i][k] = A[i][k] / U[k][k];
      for (let j = k; j < n; j++) A[i][j] -= L[i][k] * U[k][j];
    }
  }
  const Pb = matVec(P, b);
  const y = Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let sum = Pb[i];
    for (let j = 0; j < i; j++) sum -= L[i][j] * y[j];
    y[i] = sum;
  }
  const x = Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let sum = y[i];
    for (let j = i + 1; j < n; j++) sum -= U[i][j] * x[j];
    x[i] = sum / U[i][i];
  }
  return { x, L, U, iterations: [] };
}

function jacobi(A, b, tol, maxIter) {
  const n = A.length;
  let x = Array(n).fill(0);
  const iterations = [];
  for (let it = 1; it <= maxIter; it++) {
    const xn = Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      let s = b[i];
      for (let j = 0; j < n; j++) if (j !== i) s -= A[i][j] * x[j];
      xn[i] = s / A[i][i];
    }
    const err = normInf(sub(xn, x));
    iterations.push({ it, x: xn.slice(), err });
    x = xn;
    if (err < tol) break;
  }
  return { x, iterations };
}

function gaussSeidel(A, b, tol, maxIter, omega = 1) {
  const n = A.length;
  let x = Array(n).fill(0);
  const iterations = [];
  for (let it = 1; it <= maxIter; it++) {
    const old = x.slice();
    for (let i = 0; i < n; i++) {
      let s = b[i];
      for (let j = 0; j < n; j++) if (j !== i) s -= A[i][j] * x[j];
      const gs = s / A[i][i];
      x[i] = (1 - omega) * x[i] + omega * gs;
    }
    const err = normInf(sub(x, old));
    iterations.push({ it, x: x.slice(), err });
    if (err < tol) break;
  }
  return { x, iterations };
}

function isSymmetric(A) {
  for (let i = 0; i < A.length; i++) for (let j = i + 1; j < A.length; j++) if (Math.abs(A[i][j] - A[j][i]) > 1e-10) return false;
  return true;
}

function conjugateGradient(Ain, bin, tol, maxIter) {
  let A = Ain.map(r => r.slice());
  let b = bin.slice();
  let note = '';
  if (!isSymmetric(A)) {
    const AT = transpose(A);
    A = matMul(AT, A);
    b = matVec(AT, b);
    note = 'La matriz no era simétrica; se resolvieron ecuaciones normales AᵀAx=Aᵀb.';
  }
  const n = A.length;
  let x = Array(n).fill(0);
  let r = sub(b, matVec(A, x));
  let p = r.slice();
  let rsold = dot(r, r);
  const iterations = [];
  for (let it = 1; it <= maxIter; it++) {
    const Ap = matVec(A, p);
    const denom = dot(p, Ap);
    if (Math.abs(denom) < 1e-14) break;
    const alpha = rsold / denom;
    x = add(x, scale(p, alpha));
    r = sub(r, scale(Ap, alpha));
    const rsnew = dot(r, r);
    const err = Math.sqrt(rsnew);
    iterations.push({ it, x: x.slice(), err });
    if (err < tol) break;
    p = add(r, scale(p, rsnew / rsold));
    rsold = rsnew;
  }
  return { x, iterations, note };
}

function inverseMatrix(A) {
  const n = A.length;
  const inv = zeros(n);
  for (let j = 0; j < n; j++) {
    const e = Array(n).fill(0); e[j] = 1;
    const col = solveGauss(A, e);
    for (let i = 0; i < n; i++) inv[i][j] = col[i];
  }
  return inv;
}

function norm1Matrix(A) {
  let max = 0;
  for (let j = 0; j < A[0].length; j++) {
    let s = 0;
    for (let i = 0; i < A.length; i++) s += Math.abs(A[i][j]);
    max = Math.max(max, s);
  }
  return max;
}

function conditionNumber(A) {
  try {
    return norm1Matrix(A) * norm1Matrix(inverseMatrix(A));
  } catch (err) {
    return Infinity;
  }
}

function htmlTable(headers, rows) {
  const head = `<tr>${headers.map(h => `<th>${h}</th>`).join('')}</tr>`;
  const body = rows.map(r => `<tr>${r.map(c => `<td>${c}</td>`).join('')}</tr>`).join('');
  return `<table>${head}${body}</table>`;
}

function createChart(id, type, labels, datasets) {
  if (charts[id]) charts[id].destroy();
  const ctx = document.getElementById(id);
  charts[id] = new Chart(ctx, {
    type,
    data: { labels, datasets },
    options: {
      responsive: true,
      interaction: { mode: 'index', intersect: false },
      plugins: { legend: { position: 'bottom' } },
      scales: { y: { beginAtZero: false } }
    }
  });
}

function runLinear() {
  const out = document.getElementById('linearResult');
  try {
    const A = parseMatrix(document.getElementById('linA').value);
    const b = parseVector(document.getElementById('linB').value);
    assertSquare(A);
    if (b.length !== A.length) throw new Error('El vector b debe tener el mismo tamaño que la matriz.');
    if (A.some((r, i) => Math.abs(r[i]) < 1e-12)) throw new Error('Hay un valor diagonal nulo; los métodos iterativos no pueden aplicarse directamente.');
    const method = document.getElementById('linMethod').value;
    const tol = Number(document.getElementById('linTol').value);
    const maxIter = Number(document.getElementById('linMaxIter').value);
    const omega = Number(document.getElementById('omega').value);
    let result, methodName;
    if (method === 'lu') { result = luSolve(A, b); methodName = 'LU'; }
    if (method === 'jacobi') { result = jacobi(A, b, tol, maxIter); methodName = 'Jacobi'; }
    if (method === 'gs') { result = gaussSeidel(A, b, tol, maxIter, 1); methodName = 'Gauss-Seidel'; }
    if (method === 'sor') { result = gaussSeidel(A, b, tol, maxIter, omega); methodName = `SOR con ω=${omega}`; }
    if (method === 'cg') { result = conjugateGradient(A, b, tol, maxIter); methodName = 'Gradiente conjugado'; }
    const x = result.x;
    const res = residualNorm(A, x, b);
    const cond = conditionNumber(A);
    const stable = cond < 100 ? '<span class="ok">estable</span>' : cond < 10000 ? '<span class="warn">moderadamente sensible</span>' : '<span class="bad">mal condicionado</span>';
    const rows = x.map((v, i) => [`Zona ${i + 1}`, `${fmt(v, 4)} unidades`]);
    const iterRows = (result.iterations || []).slice(0, 12).map(it => [it.it, it.x.map(v => fmt(v, 4)).join(', '), fmt(it.err, 6)]);
    out.innerHTML = `
      <h3>Resultado por ${methodName}</h3>
      <p>${result.note || ''}</p>
      ${htmlTable(['Variable', 'Cantidad recomendada'], rows)}
      <p><strong>Norma del residuo:</strong> ${fmt(res, 6)} | <strong>Número de condición aproximado:</strong> ${fmt(cond, 3)} → ${stable}</p>
      <p><strong>Interpretación:</strong> la solución indica la cantidad que debe asignarse a cada zona para cumplir la demanda bajo las restricciones ingresadas. Si el número de condición es alto, una pequeña variación de demanda o bloqueo puede cambiar bastante el resultado.</p>
      ${iterRows.length ? `<h4>Primeras iteraciones</h4>${htmlTable(['Iteración', 'x aproximado', 'Error'], iterRows)}` : ''}
    `;
    createChart('linearChart', 'bar', ['Zona 1/Norte', 'Zona 2/Centro', 'Zona 3/Sur'], [{ label: 'Cantidad asignada', data: x }]);
  } catch (err) {
    out.innerHTML = `<p class="bad">Error: ${err.message}</p>`;
  }
}

function reserveDerivative(t, R, entrada, consumo, crecimiento) {
  const demand = consumo * (1 + (crecimiento / 100) * t);
  return entrada - demand;
}

function simulateScalarODE(method, R0, entrada, consumo, crecimiento, days, h) {
  const steps = Math.ceil(days / h);
  let t = 0, R = R0;
  const data = [{ t, R }];
  for (let i = 0; i < steps; i++) {
    if (method === 'euler') {
      R = R + h * reserveDerivative(t, R, entrada, consumo, crecimiento);
    } else if (method === 'heun') {
      const k1 = reserveDerivative(t, R, entrada, consumo, crecimiento);
      const pred = R + h * k1;
      const k2 = reserveDerivative(t + h, pred, entrada, consumo, crecimiento);
      R = R + h * (k1 + k2) / 2;
    } else {
      const k1 = reserveDerivative(t, R, entrada, consumo, crecimiento);
      const k2 = reserveDerivative(t + h / 2, R + h * k1 / 2, entrada, consumo, crecimiento);
      const k3 = reserveDerivative(t + h / 2, R + h * k2 / 2, entrada, consumo, crecimiento);
      const k4 = reserveDerivative(t + h, R + h * k3, entrada, consumo, crecimiento);
      R = R + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    t = Number((t + h).toFixed(10));
    data.push({ t, R: Math.max(R, 0) });
    if (R < 0) R = 0;
  }
  return data;
}

function criticalDay(data, crit) {
  const point = data.find(p => p.R <= crit);
  return point ? point.t : null;
}

function runReserves() {
  const out = document.getElementById('reservesResult');
  const R0 = Number(document.getElementById('r0').value);
  const entrada = Number(document.getElementById('entrada').value);
  const consumo = Number(document.getElementById('consumo').value);
  const crec = Number(document.getElementById('crecConsumo').value);
  const crit = Number(document.getElementById('critico').value);
  const days = Number(document.getElementById('diasReservas').value);
  const h = Number(document.getElementById('hReservas').value);
  const euler = simulateScalarODE('euler', R0, entrada, consumo, crec, days, h);
  const heun = simulateScalarODE('heun', R0, entrada, consumo, crec, days, h);
  const rk4 = simulateScalarODE('rk4', R0, entrada, consumo, crec, days, h);
  const rows = [
    ['Euler', criticalDay(euler, crit) ?? 'No llega al crítico', fmt(euler[euler.length - 1].R, 2)],
    ['Heun', criticalDay(heun, crit) ?? 'No llega al crítico', fmt(heun[heun.length - 1].R, 2)],
    ['RK4', criticalDay(rk4, crit) ?? 'No llega al crítico', fmt(rk4[rk4.length - 1].R, 2)]
  ];
  out.innerHTML = `
    <h3>Resultado de la reserva</h3>
    ${htmlTable(['Método', 'Día crítico', 'Reserva final'], rows)}
    <p><strong>Interpretación:</strong> si el consumo supera de forma sostenida la entrada, la reserva disminuye hasta un nivel crítico. RK4 se toma como referencia por su mayor estabilidad numérica.</p>
  `;
  createChart('reservesChart', 'line', rk4.map(p => p.t), [
    { label: 'Euler', data: euler.map(p => p.R), fill: false },
    { label: 'Heun', data: heun.map(p => p.R), fill: false },
    { label: 'RK4', data: rk4.map(p => p.R), fill: false },
    { label: 'Nivel crítico', data: rk4.map(() => crit), borderDash: [6, 6], fill: false }
  ]);
}

function lagrangeEval(x, y, xp) {
  let sum = 0;
  for (let i = 0; i < x.length; i++) {
    let L = 1;
    for (let j = 0; j < x.length; j++) if (j !== i) L *= (xp - x[j]) / (x[i] - x[j]);
    sum += y[i] * L;
  }
  return sum;
}

function newtonCoefficients(x, y) {
  const n = x.length;
  const coef = y.slice();
  for (let j = 1; j < n; j++) for (let i = n - 1; i >= j; i--) coef[i] = (coef[i] - coef[i - 1]) / (x[i] - x[i - j]);
  return coef;
}

function newtonEval(x, coef, xp) {
  let p = coef[coef.length - 1];
  for (let i = coef.length - 2; i >= 0; i--) p = p * (xp - x[i]) + coef[i];
  return p;
}

function splineNatural(x, y) {
  const n = x.length;
  const a = y.slice();
  const b = Array(n - 1).fill(0), d = Array(n - 1).fill(0), h = Array(n - 1).fill(0);
  for (let i = 0; i < n - 1; i++) h[i] = x[i + 1] - x[i];
  const alpha = Array(n).fill(0);
  for (let i = 1; i < n - 1; i++) alpha[i] = 3 / h[i] * (a[i + 1] - a[i]) - 3 / h[i - 1] * (a[i] - a[i - 1]);
  const c = Array(n).fill(0), l = Array(n).fill(0), mu = Array(n).fill(0), z = Array(n).fill(0);
  l[0] = 1;
  for (let i = 1; i < n - 1; i++) {
    l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
    mu[i] = h[i] / l[i];
    z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
  }
  l[n - 1] = 1;
  for (let j = n - 2; j >= 0; j--) {
    c[j] = z[j] - mu[j] * c[j + 1];
    b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
    d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
  }
  return { a, b, c, d, x };
}

function splineEval(s, xp) {
  const x = s.x;
  let i = x.length - 2;
  for (let j = 0; j < x.length - 1; j++) if (xp >= x[j] && xp <= x[j + 1]) { i = j; break; }
  const dx = xp - x[i];
  return s.a[i] + s.b[i] * dx + s.c[i] * dx * dx + s.d[i] * dx * dx * dx;
}

function runInterpolation() {
  const out = document.getElementById('interpResult');
  try {
    const x = parseVector(document.getElementById('interpX').value);
    const y = parseVector(document.getElementById('interpY').value);
    const xp = Number(document.getElementById('interpEval').value);
    if (x.length !== y.length || x.length < 3) throw new Error('Debe ingresar al menos 3 puntos y la misma cantidad de días y precios.');
    const pairs = x.map((v, i) => [v, y[i]]).sort((a, b) => a[0] - b[0]);
    const xs = pairs.map(p => p[0]), ys = pairs.map(p => p[1]);
    const coef = newtonCoefficients(xs, ys);
    const spl = splineNatural(xs, ys);
    const lag = lagrangeEval(xs, ys, xp);
    const neu = newtonEval(xs, coef, xp);
    const spi = splineEval(spl, xp);
    const min = xs[0], max = xs[xs.length - 1];
    const grid = [];
    for (let i = 0; i <= 80; i++) grid.push(min + (max - min) * i / 80);
    const increase = ((ys[ys.length - 1] - ys[0]) / ys[0]) * 100;
    out.innerHTML = `
      <h3>Estimación de precio para el día ${xp}</h3>
      ${htmlTable(['Método', 'Precio estimado Bs'], [['Lagrange', fmt(lag, 2)], ['Newton', fmt(neu, 2)], ['Spline cúbico natural', fmt(spi, 2)]])}
      <p><strong>Incremento total observado:</strong> ${fmt(increase, 2)}% desde el primer hasta el último dato.</p>
      <p><strong>Interpretación:</strong> la interpolación permite estimar precios en días sin registro. Para datos muy dispersos, el spline suele ser más recomendable porque evita oscilaciones fuertes del polinomio global.</p>
    `;
    createChart('interpChart', 'line', grid.map(v => fmt(v, 1)), [
      { label: 'Lagrange', data: grid.map(v => lagrangeEval(xs, ys, v)), fill: false },
      { label: 'Newton', data: grid.map(v => newtonEval(xs, coef, v)), fill: false },
      { label: 'Spline', data: grid.map(v => splineEval(spl, v)), fill: false },
      { label: 'Datos reales', data: grid.map(v => {
        const idx = xs.findIndex(xx => Math.abs(xx - v) < 0.2);
        return idx >= 0 ? ys[idx] : null;
      }), showLine: false, pointRadius: 5 }
    ]);
  } catch (err) {
    out.innerHTML = `<p class="bad">Error: ${err.message}</p>`;
  }
}

function loadDefaultIntegration() {
  const xs = Array.from({ length: 31 }, (_, i) => i);
  const ys = xs.map(d => 82 + 1.15 * d + 3 * Math.sin(d / 4));
  document.getElementById('intX').value = xs.join(',');
  document.getElementById('intY').value = ys.map(v => fmt(v, 2)).join(',');
}

function isEquallySpaced(x) {
  const h = x[1] - x[0];
  return x.every((v, i) => i === 0 || Math.abs((v - x[i - 1]) - h) < 1e-8);
}

function trapezoid(x, y) {
  let area = 0;
  for (let i = 0; i < x.length - 1; i++) area += (x[i + 1] - x[i]) * (y[i] + y[i + 1]) / 2;
  return area;
}

function simpson13(x, y) {
  const n = x.length - 1;
  if (!isEquallySpaced(x) || n % 2 !== 0) return NaN;
  const h = (x[x.length - 1] - x[0]) / n;
  let s = y[0] + y[n];
  for (let i = 1; i < n; i++) s += (i % 2 === 0 ? 2 : 4) * y[i];
  return h * s / 3;
}

function simpson38(x, y) {
  const n = x.length - 1;
  if (!isEquallySpaced(x) || n % 3 !== 0) return NaN;
  const h = (x[x.length - 1] - x[0]) / n;
  let s = y[0] + y[n];
  for (let i = 1; i < n; i++) s += (i % 3 === 0 ? 2 : 3) * y[i];
  return 3 * h * s / 8;
}

function runIntegration() {
  const out = document.getElementById('integrationResult');
  try {
    const x = parseVector(document.getElementById('intX').value);
    const y = parseVector(document.getElementById('intY').value);
    if (x.length !== y.length || x.length < 2) throw new Error('Debe ingresar la misma cantidad de valores x y costos.');
    const trap = trapezoid(x, y);
    const s13 = simpson13(x, y);
    const s38 = simpson38(x, y);
    const base = y[0] * (x[x.length - 1] - x[0]);
    const loss = trap - base;
    out.innerHTML = `
      <h3>Costo acumulado mensual</h3>
      ${htmlTable(['Método', 'Gasto acumulado Bs'], [
        ['Trapecio', fmt(trap, 2)],
        ['Simpson 1/3', Number.isNaN(s13) ? 'No aplica: requiere subintervalos pares e h constante' : fmt(s13, 2)],
        ['Simpson 3/8', Number.isNaN(s38) ? 'No aplica: requiere subintervalos múltiplos de 3 e h constante' : fmt(s38, 2)]
      ])}
      <p><strong>Gasto si el precio no subía:</strong> ${fmt(base, 2)} Bs.</p>
      <p><strong>Pérdida aproximada de poder adquisitivo:</strong> ${fmt(loss, 2)} Bs usando trapecio como referencia.</p>
      <p><strong>Interpretación:</strong> el área bajo la curva representa el gasto acumulado. Cuando la curva se eleva con el tiempo, la diferencia frente al costo constante mide la pérdida económica aproximada.</p>
    `;
    createChart('integrationChart', 'line', x, [{ label: 'Costo diario Bs', data: y, fill: true }]);
  } catch (err) {
    out.innerHTML = `<p class="bad">Error: ${err.message}</p>`;
  }
}

function rootFunction(a, b, c, x) { return a * x * x + b * x + c; }
function rootDerivative(a, b, x) { return 2 * a * x + b; }

function bisection(a, b, c, l, u, tol = 1e-6, maxIter = 50) {
  const rows = [];
  let fl = rootFunction(a, b, c, l), fu = rootFunction(a, b, c, u);
  if (fl * fu > 0) return { root: NaN, rows, error: 'No hay cambio de signo en el intervalo.' };
  let mid = l;
  for (let it = 1; it <= maxIter; it++) {
    mid = (l + u) / 2;
    const fm = rootFunction(a, b, c, mid);
    const err = Math.abs(u - l) / 2;
    rows.push({ it, x: mid, fx: fm, err });
    if (Math.abs(fm) < tol || err < tol) break;
    if (fl * fm < 0) { u = mid; fu = fm; } else { l = mid; fl = fm; }
  }
  return { root: mid, rows };
}

function newtonRoot(a, b, c, x0, tol = 1e-6, maxIter = 30) {
  const rows = [];
  let x = x0;
  for (let it = 1; it <= maxIter; it++) {
    const fx = rootFunction(a, b, c, x);
    const dfx = rootDerivative(a, b, x);
    if (Math.abs(dfx) < 1e-14) break;
    const xn = x - fx / dfx;
    const err = Math.abs(xn - x);
    rows.push({ it, x: xn, fx: rootFunction(a, b, c, xn), err });
    x = xn;
    if (err < tol) break;
  }
  return { root: x, rows };
}

function secantRoot(a, b, c, x0, x1, tol = 1e-6, maxIter = 30) {
  const rows = [];
  let p0 = x0, p1 = x1;
  for (let it = 1; it <= maxIter; it++) {
    const f0 = rootFunction(a, b, c, p0), f1 = rootFunction(a, b, c, p1);
    if (Math.abs(f1 - f0) < 1e-14) break;
    const p = p1 - f1 * (p1 - p0) / (f1 - f0);
    const err = Math.abs(p - p1);
    rows.push({ it, x: p, fx: rootFunction(a, b, c, p), err });
    p0 = p1; p1 = p;
    if (err < tol) break;
  }
  return { root: p1, rows };
}

function estimateOrder(rows, root) {
  if (rows.length < 4 || !isFinite(root)) return 'No estimable';
  const e = rows.slice(-4).map(r => Math.abs(r.x - root)).filter(v => v > 1e-14);
  if (e.length < 3) return 'No estimable';
  const p = Math.log(e[e.length - 1] / e[e.length - 2]) / Math.log(e[e.length - 2] / e[e.length - 3]);
  return fmt(p, 3);
}

function runRoots() {
  const out = document.getElementById('rootsResult');
  const a = Number(document.getElementById('rootA').value);
  const b = Number(document.getElementById('rootB').value);
  const c = Number(document.getElementById('rootC').value);
  const l = Number(document.getElementById('rootL').value);
  const u = Number(document.getElementById('rootU').value);
  const x0 = Number(document.getElementById('rootX0').value);
  const sx0 = Number(document.getElementById('secX0').value);
  const sx1 = Number(document.getElementById('secX1').value);
  const bis = bisection(a, b, c, l, u);
  const newt = newtonRoot(a, b, c, x0);
  const sec = secantRoot(a, b, c, sx0, sx1);
  const rows = [
    ['Bisección', bis.error ? bis.error : fmt(bis.root, 6), estimateOrder(bis.rows, bis.root)],
    ['Newton-Raphson', fmt(newt.root, 6), estimateOrder(newt.rows, newt.root)],
    ['Secante', fmt(sec.root, 6), estimateOrder(sec.rows, sec.root)]
  ];
  const iterRows = (bis.rows.length ? bis.rows : newt.rows).slice(0, 12).map(r => [r.it, fmt(r.x, 6), fmt(r.fx, 6), fmt(r.err, 6)]);
  out.innerHTML = `
    <h3>Raíz o umbral crítico</h3>
    ${htmlTable(['Método', 'Raíz aproximada', 'Orden de convergencia estimado'], rows)}
    ${htmlTable(['Iteración', 'x', 'f(x)', 'Error'], iterRows)}
    <p><strong>Interpretación:</strong> el valor de x encontrado representa el punto crítico. Si x es día, entonces cerca de ese día el gasto acumulado iguala el ingreso disponible; después de ese punto aparece déficit.</p>
  `;
  const labels = [], values = [];
  for (let i = 0; i <= 80; i++) {
    const x = l + (u - l) * i / 80;
    labels.push(fmt(x, 2)); values.push(rootFunction(a, b, c, x));
  }
  createChart('rootsChart', 'line', labels, [
    { label: 'f(x)', data: values, fill: false },
    { label: 'Eje cero', data: labels.map(() => 0), borderDash: [6, 6], fill: false }
  ]);
}

function runRumors() {
  const out = document.getElementById('rumorResult');
  try {
    const A = parseMatrix(document.getElementById('rumorA').value);
    const b = parseVector(document.getElementById('rumorB').value);
    const rumorPct = Number(document.getElementById('rumorPct').value) / 100;
    const stockPct = Number(document.getElementById('stockPct').value) / 100;
    assertSquare(A);
    const baseX = solveGauss(A, b);
    const pertB = b.map(v => v * (1 + rumorPct));
    const pertA = A.map((row, i) => row.map((v, j) => i === j ? v * (1 - stockPct) : v));
    const pertX = solveGauss(pertA, pertB);
    const cond = conditionNumber(A);
    const relB = norm2(sub(pertB, b)) / Math.max(norm2(b), 1e-14);
    const relX = norm2(sub(pertX, baseX)) / Math.max(norm2(baseX), 1e-14);
    const variations = pertX.map((v, i) => Math.abs((v - baseX[i]) / (Math.abs(baseX[i]) || 1)) * 100);
    const vulnIndex = variations.indexOf(Math.max(...variations));
    const rows = baseX.map((v, i) => [`Zona ${i + 1}`, fmt(v, 4), fmt(pertX[i], 4), `${fmt(variations[i], 2)}%`]);
    const level = cond > 1000 || relX > relB * 5 ? '<span class="bad">mal condicionado o muy sensible</span>' : '<span class="ok">estable ante la perturbación ingresada</span>';
    out.innerHTML = `
      <h3>Sensibilidad ante rumores</h3>
      ${htmlTable(['Zona', 'Solución base', 'Con rumor/stock reducido', 'Variación'], rows)}
      <p><strong>Número de condición aproximado:</strong> ${fmt(cond, 3)}. Cambio relativo en demanda: ${fmt(relB * 100, 2)}%. Cambio relativo en solución: ${fmt(relX * 100, 2)}%. Diagnóstico: ${level}.</p>
      <p><strong>Zona más vulnerable:</strong> Zona ${vulnIndex + 1}, porque presenta la mayor variación relativa.</p>
      <p><strong>Interpretación:</strong> cuando el sistema está mal condicionado, rumores o compras impulsivas pueden generar cambios desproporcionados en la asignación requerida.</p>
    `;
    createChart('rumorChart', 'bar', ['Zona 1', 'Zona 2', 'Zona 3'], [
      { label: 'Base', data: baseX },
      { label: 'Con rumor', data: pertX }
    ]);
  } catch (err) {
    out.innerHTML = `<p class="bad">Error: ${err.message}</p>`;
  }
}

function socialDeriv(state, p) {
  const [N, M, D] = state;
  return [
    -p.a * N * M + p.b * D,
    p.a * N * M - p.c * M * D,
    p.k * M - p.r * D
  ];
}

function vectorStepHeun(state, h, p) {
  const k1 = socialDeriv(state, p);
  const pred = add(state, scale(k1, h));
  const k2 = socialDeriv(pred, p);
  return add(state, scale(add(k1, k2), h / 2)).map(v => Math.max(0, v));
}

function vectorStepRK4(state, h, p) {
  const k1 = socialDeriv(state, p);
  const k2 = socialDeriv(add(state, scale(k1, h / 2)), p);
  const k3 = socialDeriv(add(state, scale(k2, h / 2)), p);
  const k4 = socialDeriv(add(state, scale(k3, h)), p);
  const sum = add(add(k1, scale(k2, 2)), add(scale(k3, 2), k4));
  return add(state, scale(sum, h / 6)).map(v => Math.max(0, v));
}

function runSocial() {
  const out = document.getElementById('socialResult');
  const state0 = [Number(document.getElementById('N0').value), Number(document.getElementById('M0').value), Number(document.getElementById('D0').value)];
  const p = {
    a: Number(document.getElementById('coefA').value),
    b: Number(document.getElementById('coefB').value),
    c: Number(document.getElementById('coefC').value),
    k: Number(document.getElementById('coefK').value),
    r: Number(document.getElementById('coefR').value)
  };
  const days = Number(document.getElementById('socialDays').value);
  const h = Number(document.getElementById('socialH').value);
  const method = document.getElementById('socialMethod').value;
  const data = [{ t: 0, N: state0[0], M: state0[1], D: state0[2] }];
  let state = state0.slice();
  for (let t = h; t <= days + 1e-9; t += h) {
    state = method === 'heun' ? vectorStepHeun(state, h, p) : vectorStepRK4(state, h, p);
    data.push({ t: Number(t.toFixed(10)), N: state[0], M: state[1], D: state[2] });
  }
  const firstM = data[0].M;
  const last = data[data.length - 1];
  const trend = last.M > firstM * 1.2 ? '<span class="bad">tiende a masificarse</span>' : last.M < firstM * 0.8 ? '<span class="ok">tiende a disminuir</span>' : '<span class="warn">tiende a estabilizarse</span>';
  const rows = data.filter((_, i) => i % Math.max(1, Math.floor(data.length / 10)) === 0 || i === data.length - 1).map(p => [fmt(p.t, 1), fmt(p.N, 2), fmt(p.M, 2), fmt(p.D, 2)]);
  out.innerHTML = `
    <h3>Dinámica social por ${method === 'heun' ? 'Heun' : 'RK4'}</h3>
    ${htmlTable(['Día', 'N(t)', 'M(t)', 'D(t)'], rows)}
    <p><strong>Diagnóstico:</strong> el número de manifestantes ${trend}. Manifestantes iniciales: ${fmt(firstM, 2)}; manifestantes finales: ${fmt(last.M, 2)}.</p>
    <p><strong>Interpretación:</strong> si aumenta la tasa de diálogo o mediación efectiva, el término cMD reduce M(t). Si no existen mediadores o la influencia a es alta, el conflicto puede crecer con mayor rapidez.</p>
  `;
  createChart('socialChart', 'line', data.map(p => p.t), [
    { label: 'Neutrales N(t)', data: data.map(p => p.N), fill: false },
    { label: 'Manifestantes M(t)', data: data.map(p => p.M), fill: false },
    { label: 'Mediadores D(t)', data: data.map(p => p.D), fill: false }
  ]);
}

window.addEventListener('DOMContentLoaded', () => {
  loadDefaultIntegration();
  runLinear();
  runReserves();
  runInterpolation();
  runIntegration();
  runRoots();
  runRumors();
  runSocial();
});
