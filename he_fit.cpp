// he_fit.cpp - подгонка аналитической модели (вариант B) под аттестованный эталон.
//
//   he_fit данные.csv|- [-o coef.txt] [--init coef.txt] [--dmax 4] [--ridge 1e-10]
//          [--abinitio data/he4_abinitio_bich2007.csv] [--keep-dilute]
//
// Вход: CSV с заголовком; столбцы (регистр не важен, порядок любой):
//   T [К], P [Па], rho [кг/м3] – обязательны;  mu [Па с], lambda [Вт/(м К)] – необязательны.
// Данные получают из аттестованного источника (REFPROP, CoolProp, NIST WebBook) на сетке по всей
// области 50-3000 К, 10 Па - 100 МПа (рекомендуется >= 1500 точек, равномерно в ln T и ln P,
// с плотным покрытием низких T и высоких P).
//
// Что подгоняется (все задачи – линейные по искомым коэффициентам):
//   1) поправочные члены остаточной энергии Гельмгольца  n_i (T/300)^(-t_i) (rho/100)^(d_i)
//      по невязке коэффициента сжимаемости Z (веса 1/Z => минимизируется относительная ошибка rho);
//   2) полиномы ln mu0(T), ln lambda0(T) степени 3 по s = ln(T/300) и добавка плотностной поправки
//      sum d_jk s^k r^j (j=1..3, k=0..2, r = rho/100) к поправке Энскога (последняя берется
//      из подогнанного уравнения состояния).
// Режим привязки к ab initio (--abinitio, Bich, Hellmann, Vogel 2007): B(T) (полином степени 5 по s) и
// ln mu0(T), ln lambda0(T) (степени 3) подгоняются к таблице стандартных значений в диапазоне 50-3000 К.
// Данные CSV можно заменить на '-', тогда выполняется только эта привязка (поправочные члены = 0).
// Если CSV с плотным газом задан, флаг --keep-dilute не дает подгонке переноса менять mu0 и lambda0:
// по плотным данным подгоняется только добавка D(s, rho). Рекомендуемая последовательность:
//   he_fit - --abinitio data/he4_abinitio_bich2007.csv -o dilute.txt
//   he_fit dense.csv --init dilute.txt --keep-dilute -o coef.txt
// Выход: файл коэффициентов для AnalyticModel::load / he_calc --model analytic:файл / he_table_build --coef.
// Утилита печатает ошибки модели относительно данных ДО и ПОСЛЕ подгонки; критерий принятия – по ним.
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include "he_eos.h"
#include "he_console.h"

using namespace he;

struct Row { double T, P, rho, mu, lam; };

// Решение min |A x - y|^2 + ridge*|x|^2 (после нормировки столбцов) методом Хаусхолдера (QR)
static bool lsq(std::vector<double> A, std::vector<double> y, int m, int n, double ridge, std::vector<double>& x) {
  std::vector<double> sc(n, 1.0);
  for (int j = 0; j < n; ++j) {
    double s = 0;
    for (int i = 0; i < m; ++i) s += A[(size_t)i * n + j] * A[(size_t)i * n + j];
    sc[j] = s > 0 ? std::sqrt(s) : 1.0;
    for (int i = 0; i < m; ++i) A[(size_t)i * n + j] /= sc[j];
  }
  const double rr = std::sqrt(ridge);
  for (int j = 0; j < n; ++j) {  // добавочные строки регуляризации
    for (int k = 0; k < n; ++k) A.push_back(k == j ? rr : 0.0);
    y.push_back(0.0);
  }
  const int M = m + n;
  for (int k = 0; k < n; ++k) {
    double nrm = 0;
    for (int i = k; i < M; ++i) nrm += A[(size_t)i * n + k] * A[(size_t)i * n + k];
    nrm = std::sqrt(nrm);
    if (nrm == 0) return false;
    const double alpha = A[(size_t)k * n + k] > 0 ? -nrm : nrm;
    std::vector<double> v(M - k);
    for (int i = k; i < M; ++i) v[i - k] = A[(size_t)i * n + k];
    v[0] -= alpha;
    double vv = 0;
    for (double t : v) vv += t * t;
    if (vv == 0) continue;
    for (int j = k; j < n; ++j) {
      double d = 0;
      for (int i = k; i < M; ++i) d += v[i - k] * A[(size_t)i * n + j];
      d = 2 * d / vv;
      for (int i = k; i < M; ++i) A[(size_t)i * n + j] -= d * v[i - k];
    }
    double d = 0;
    for (int i = k; i < M; ++i) d += v[i - k] * y[i];
    d = 2 * d / vv;
    for (int i = k; i < M; ++i) y[i] -= d * v[i - k];
  }
  x.assign(n, 0.0);
  for (int k = n - 1; k >= 0; --k) {
    double s = y[k];
    for (int j = k + 1; j < n; ++j) s -= A[(size_t)k * n + j] * x[j];
    x[k] = s / A[(size_t)k * n + k];
  }
  for (int j = 0; j < n; ++j) x[j] /= sc[j];
  return true;
}

static std::vector<std::string> split(const std::string& s) {
  std::vector<std::string> o;
  std::string cur;
  for (char c : s) { if (c == ',' || c == ';' || c == '\t') { o.push_back(cur); cur.clear(); } else if (c != '\r' && c != '\n') cur += c; }
  o.push_back(cur);
  return o;
}

struct Stat { double max = 0, rms = 0; int n = 0; void add(double e) { max = std::fmax(max, std::fabs(e)); rms += e * e; ++n; } double r() const { return n ? std::sqrt(rms / n) : 0; } };

static void evaluate(const AnalyticModel& m, const std::vector<Row>& d, const char* tag) {
  Stat sr, sm, sl;
  for (const Row& r : d) {
    Result o;
    if (m.PT(r.P, r.T, HP_RHO | (r.mu > 0 ? (unsigned)HP_MU : 0u) | (r.lam > 0 ? (unsigned)HP_LAM : 0u), o) != HE_OK) { sr.add(1.0); continue; }
    sr.add(o.rho / r.rho - 1);
    if (r.mu > 0) sm.add(o.mu / r.mu - 1);
    if (r.lam > 0) sl.add(o.lambda / r.lam - 1);
  }
  std::printf("  %-7s rho: max %.3e, rms %.3e", tag, sr.max, sr.r());
  if (sm.n) std::printf(" | mu: max %.3e, rms %.3e", sm.max, sm.r());
  if (sl.n) std::printf(" | lambda: max %.3e, rms %.3e", sl.max, sl.r());
  std::printf("\n");
}

// Привязка B(T), mu0(T), lambda0(T) к ab initio-таблице (столбцы: T, B, C, eta[мкПа с], lambda[мВт/(м К)])
static bool stage_abinitio(AnalyticModel& m, const std::string& file) {
  FILE* f = std::fopen(file.c_str(), "r");
  if (!f) { std::fprintf(stderr, "Ошибка: не удалось открыть %s\n", file.c_str()); return false; }
  std::vector<double> T, B, eta, lam;
  char line[512];
  while (std::fgets(line, sizeof line, f)) {
    double t, b, c, e, l;
    if (line[0] == '#' || !(line[0] == '-' || (line[0] >= '0' && line[0] <= '9'))) continue;
    if (std::sscanf(line, "%lf,%lf,%lf,%lf,%lf", &t, &b, &c, &e, &l) == 5 && t >= T_MIN && t <= T_MAX) {
      T.push_back(t); B.push_back(b); eta.push_back(e * 1e-6); lam.push_back(l * 1e-3);
    }
  }
  std::fclose(f);
  const int n = (int)T.size();
  if (n < 12) { std::fprintf(stderr, "Ошибка: в %s меньше 12 точек в диапазоне 50-3000 К\n", file.c_str()); return false; }
  std::vector<double> AB((size_t)n * 6), yb(n, 1.0), AM((size_t)n * 4), ym(n), yl(n), cb, cm, cl;
  for (int i = 0; i < n; ++i) {
    const double s = std::log(T[i] / 300.0);
    for (int k = 0; k < 6; ++k) AB[(size_t)i * 6 + k] = std::pow(s, k) / B[i];   // относительная ошибка B
    for (int k = 0; k < 4; ++k) AM[(size_t)i * 4 + k] = std::pow(s, k);
    ym[i] = std::log(eta[i]);
    yl[i] = std::log(lam[i]);
  }
  if (!lsq(AB, yb, n, 6, 1e-12, cb) || !lsq(AM, ym, n, 4, 1e-12, cm) || !lsq(AM, yl, n, 4, 1e-12, cl)) {
    std::fprintf(stderr, "Ошибка: вырожденная система привязки к ab initio\n");
    return false;
  }
  for (int k = 0; k < 6; ++k) m.coef().Bc[k] = cb[k];
  for (int k = 0; k < 4; ++k) { m.coef().mu0[k] = cm[k]; m.coef().lam0[k] = cl[k]; }
  for (int k = 0; k < 9; ++k) m.coef().mud[k] = m.coef().lamd[k] = 0.0;
  m.finalize();
  double eb = 0, em = 0, el = 0;
  for (int i = 0; i < n; ++i) {
    const double s = std::log(T[i] / 300.0);
    double vb = 0, vm = 0, vl = 0;
    for (int k = 5; k >= 0; --k) vb = vb * s + cb[k];
    for (int k = 3; k >= 0; --k) { vm = vm * s + cm[k]; vl = vl * s + cl[k]; }
    eb = std::fmax(eb, std::fabs(vb / B[i] - 1));
    em = std::fmax(em, std::fabs(std::exp(vm) / eta[i] - 1));
    el = std::fmax(el, std::fabs(std::exp(vl) / lam[i] - 1));
  }
  std::printf("Привязка к ab initio (%d точек, 50-3000 К): макс. ошибка полиномов  B %.3f %%, mu0 %.3f %%, lambda0 %.3f %%\n", n, eb * 100, em * 100, el * 100);
  return true;
}

int main(int argc, char** argv) {
  if (argc < 2) { std::fprintf(stderr, "Использование: %s данные.csv [-o coef.txt] [--init coef.txt] [--dmax 4] [--ridge 1e-10]\n", argv[0]); return 1; }
  std::string in = argv[1], out = "he_eos_coef.txt", init, abfile;
  bool keepDilute = false;
  int dmax = 4;
  double ridge = 1e-10;
  for (int i = 2; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "-o" && i + 1 < argc) out = argv[++i];
    else if (a == "--init" && i + 1 < argc) init = argv[++i];
    else if (a == "--dmax" && i + 1 < argc) dmax = std::atoi(argv[++i]);
    else if (a == "--ridge" && i + 1 < argc) ridge = std::atof(argv[++i]);
    else if (a == "--abinitio" && i + 1 < argc) abfile = argv[++i];
    else if (a == "--keep-dilute") keepDilute = true;
    else { std::fprintf(stderr, "неизвестный параметр %s\n", a.c_str()); return 1; }
  }
  if (in == "-") {
    if (abfile.empty()) { std::fprintf(stderr, "Ошибка: вместо CSV '-' требуется --abinitio\n"); return 1; }
    AnalyticModel m0;
    std::string e0;
    if (!init.empty() && !m0.load(init.c_str(), &e0)) { std::fprintf(stderr, "Ошибка: %s\n", e0.c_str()); return 2; }
    m0.coef().ncorr = 0;
    if (!stage_abinitio(m0, abfile)) return 2;
    if (!m0.save(out.c_str())) { std::fprintf(stderr, "Ошибка: не удалось записать %s\n", out.c_str()); return 3; }
    std::printf("Коэффициенты записаны: %s\nПроверьте: he_check_abinitio data/he4_abinitio_bich2007.csv --model analytic:%s\n", out.c_str(), out.c_str());
    return 0;
  }
  FILE* f = std::fopen(in.c_str(), "r");
  if (!f) { std::fprintf(stderr, "Ошибка: не удалось открыть %s\n", in.c_str()); return 2; }
  char line[2048];
  if (!std::fgets(line, sizeof line, f)) { std::fprintf(stderr, "Ошибка: пустой файл\n"); return 2; }
  int iT = -1, iP = -1, iR = -1, iM = -1, iL = -1;
  auto hd = split(line);
  for (size_t k = 0; k < hd.size(); ++k) {
    std::string h;
    for (char c : hd[k]) if (!std::isspace((unsigned char)c)) h += (char)std::tolower((unsigned char)c);
    if (h == "t") iT = (int)k; else if (h == "p") iP = (int)k; else if (h == "rho") iR = (int)k;
    else if (h == "mu") iM = (int)k; else if (h == "lambda" || h == "lam") iL = (int)k;
  }
  if (iT < 0 || iP < 0 || iR < 0) { std::fprintf(stderr, "Ошибка: нужны столбцы T, P, rho\n"); return 2; }
  std::vector<Row> d;
  while (std::fgets(line, sizeof line, f)) {
    auto c = split(line);
    if ((int)c.size() <= std::max(iT, std::max(iP, iR))) continue;
    Row r{std::atof(c[iT].c_str()), std::atof(c[iP].c_str()), std::atof(c[iR].c_str()), 0, 0};
    if (iM >= 0 && (int)c.size() > iM) r.mu = std::atof(c[iM].c_str());
    if (iL >= 0 && (int)c.size() > iL) r.lam = std::atof(c[iL].c_str());
    if (r.T > 0 && r.P > 0 && r.rho > 0) d.push_back(r);
  }
  std::fclose(f);

  AnalyticModel m;
  std::string err;
  if (!init.empty() && !m.load(init.c_str(), &err)) { std::fprintf(stderr, "Ошибка: %s\n", err.c_str()); return 2; }
  m.coef().ncorr = 0;
  if (!abfile.empty() && !stage_abinitio(m, abfile)) return 2;
  std::printf("Прочитано %zu точек. Ошибки модели относительно данных:\n", d.size());
  evaluate(m, d, "до:");

  // ---- 1) поправочные члены по невязке Z
  std::vector<double> tt, dd;
  const double tl[] = {0, 0.5, 1, 2, 3};
  for (int k = 1; k <= dmax; ++k) for (double t : tl) { tt.push_back(t); dd.push_back(k); }
  const int nt = (int)tt.size(), mrows = (int)d.size();
  if (mrows < 2 * nt) { std::fprintf(stderr, "Ошибка: слишком мало точек (%d) для %d членов\n", mrows, nt); return 2; }
  if (nt > EOS_MAXC) { std::fprintf(stderr, "Ошибка: слишком много членов\n"); return 2; }
  std::vector<double> A((size_t)mrows * nt), y(mrows), x;
  for (int k = 0; k < mrows; ++k) {
    const Row& r = d[k];
    const double Zref = r.P / (r.rho * R_SPEC * r.T), Z0 = m.Z_of(r.T, r.rho), w = 1.0 / Zref;
    y[k] = (Zref - Z0) * w;
    for (int i = 0; i < nt; ++i) A[(size_t)k * nt + i] = w * dd[i] * std::pow(r.T / 300.0, -tt[i]) * std::pow(r.rho / 100.0, dd[i]);
  }
  if (!lsq(A, y, mrows, nt, ridge, x)) { std::fprintf(stderr, "Ошибка: вырожденная система\n"); return 2; }
  for (int i = 0; i < nt; ++i) { m.coef().cn[i] = x[i]; m.coef().ct[i] = tt[i]; m.coef().cd[i] = dd[i]; }
  m.coef().ncorr = nt;
  m.finalize();

  // ---- 2) перенос: ln f = poly3(s) + ln F_Энског(eta) + sum_{j=1..3, k=0..2} d_jk s^k r^j,  r = rho/100
  for (int which = 0; which < 2; ++which) {
    std::vector<int> ids;
    for (int k = 0; k < mrows; ++k) if ((which == 0 ? d[k].mu : d[k].lam) > 0) ids.push_back(k);
    const int nc = keepDilute ? 9 : 13, off = keepDilute ? 0 : 4;
    if ((int)ids.size() < 3 * nc) continue;
    std::vector<double> B(ids.size() * nc), yy(ids.size()), c;
    for (size_t q = 0; q < ids.size(); ++q) {
      const Row& r = d[ids[q]];
      double fm, fl;
      m.enskog(r.T, r.rho, fm, fl);
      const double s = std::log(r.T / 300.0), rr = r.rho / 100.0;
      yy[q] = std::log(which == 0 ? r.mu : r.lam) - std::log(which == 0 ? fm : fl);
      if (keepDilute) {
        const double* p0k = which == 0 ? m.coef().mu0 : m.coef().lam0;
        yy[q] -= p0k[0] + s * (p0k[1] + s * (p0k[2] + s * p0k[3]));
      } else {
        for (int k = 0; k < 4; ++k) B[q * nc + k] = std::pow(s, k);
      }
      for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) B[q * nc + off + 3 * j + k] = std::pow(rr, j + 1) * std::pow(s, k);
    }
    if (!lsq(B, yy, (int)ids.size(), nc, ridge, c)) continue;
    double* p0 = which == 0 ? m.coef().mu0 : m.coef().lam0;
    double* pd = which == 0 ? m.coef().mud : m.coef().lamd;
    if (!keepDilute) for (int k = 0; k < 4; ++k) p0[k] = c[k];
    for (int k = 0; k < 9; ++k) pd[k] = c[off + k];
  }
  std::printf("Ошибки модели после подгонки (%d поправочных членов):\n", nt);
  evaluate(m, d, "после:");
  if (!m.save(out.c_str())) { std::fprintf(stderr, "Ошибка: не удалось записать %s\n", out.c_str()); return 3; }
  std::printf("Коэффициенты записаны: %s\n", out.c_str());
  std::printf("ВНИМАНИЕ: проверьте ошибки на НЕЗАВИСИМОМ наборе точек (не использованном в подгонке),\n"
              "особенно при T < 100 К и P > 30 МПа, и при необходимости измените --dmax / набор данных.\n");
  return 0;
}
