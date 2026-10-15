// he_test.cpp - автоматические проверки обоих вариантов (A: таблица, B: аналитическая модель)
//   he_test table.bin [coef.txt]
// Проверяет: маски, термодинамические тождества, обращение (P,H)/(P,U), граничные случаи
// и экстраполяцию, согласие A с B, а также измеряет скорость.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include "he_eos.h"
#include "he_table.h"
#include "he_console.h"

using namespace he;

static int g_fail = 0;
#define CHECK(cond, ...)                                       \
  do {                                                         \
    if (!(cond)) { ++g_fail; std::printf("  [НЕУДАЧА] "); std::printf(__VA_ARGS__); std::printf("\n"); } \
  } while (0)

static const unsigned kBits[] = {HP_T, HP_RHO, HP_CP, HP_CV, HP_U, HP_H, HP_LAM, HP_MU, HP_PR, HP_DUDT, HP_DRDT, HP_DRDU, HP_DRDH, HP_DRDP};
static double get(const Result& r, unsigned b) {
  switch (b) {
    case HP_T: return r.T; case HP_RHO: return r.rho; case HP_CP: return r.cp; case HP_CV: return r.cv; case HP_U: return r.U;
    case HP_H: return r.H; case HP_LAM: return r.lambda; case HP_MU: return r.mu; case HP_PR: return r.Pr; case HP_DUDT: return r.dUdT;
    case HP_DRDT: return r.drhodT; case HP_DRDU: return r.drhodU; case HP_DRDH: return r.drhodH; default: return r.drhodP;
  }
}
static const char* kNm[] = {"T", "rho", "cp", "cv", "U", "H", "lambda", "mu", "Pr", "dUdT", "drhodT", "drhodU", "drhodH", "drhodP"};

static double pct(std::vector<double>& v, double q) {
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[(size_t)(q * (v.size() - 1))];
}

static void run_model(const Model& m, const char* tag) {
  std::printf("\n=== %s: %s ===\n", tag, m.name());
  std::mt19937_64 rng(2024);
  std::uniform_real_distribution<double> ux(std::log(T_MIN * 1.0001), std::log(T_MAX * 0.9999)), uy(std::log(P_MIN * 1.0001), std::log(P_MAX * 0.9999));
  auto rp = [&](double& P, double& T) { T = std::exp(ux(rng)); P = std::exp(uy(rng)); };

  // 1. маски
  int bad = 0;
  for (int n = 0; n < 500; ++n) {
    double P, T; rp(P, T);
    Result all;
    m.PT(P, T, HP_ALL, all);
    for (unsigned b : kBits) {
      Result r;
      if (m.PT(P, T, b, r) != HE_OK) { ++bad; continue; }
      if (b != HP_T && get(r, b) != get(all, b)) ++bad;               // побитовое совпадение
      if ((r.computed & ~(b | HP_T)) != 0) ++bad;                      // ничего лишнего
      for (unsigned c : kBits) if (c != b && c != HP_T && !std::isnan(get(r, c))) ++bad;  // прочее = NaN
    }
  }
  CHECK(bad == 0, "маски: %d расхождений", bad);
  std::printf("  маски: %s\n", bad ? "ОШИБКА" : "ок (единичная маска = полная маска побитово, лишнее не считается)");

  // 2. термодинамические тождества (конечные разности по самой модели)
  std::vector<double> e[4];
  for (int n = 0; n < 3000; ++n) {
    double P, T; rp(P, T);
    const double hT = 1e-5 * T, hP = 1e-5 * P;
    Result c, a, b, pp, pm;
    m.PT(P, T, HP_ALL, c); m.PT(P, T + hT, HP_ALL, a); m.PT(P, T - hT, HP_ALL, b);
    m.PT(P + hP, T, HP_ALL, pp); m.PT(P - hP, T, HP_ALL, pm);
    e[0].push_back(std::fabs((a.H - b.H) / (2 * hT) / c.cp - 1));
    e[1].push_back(std::fabs((a.rho - b.rho) / (2 * hT) / c.drhodT - 1));
    e[2].push_back(std::fabs((a.U - b.U) / (2 * hT) / c.dUdT - 1));
    e[3].push_back(std::fabs((pp.rho - pm.rho) / (2 * hP) / c.drhodP - 1));
  }
  const char* en[4] = {"cp = dH/dT", "drho/dT = d(rho)/dT", "dU/dT = d(U)/dT", "drho/dP = d(rho)/dP"};
  for (int q = 0; q < 4; ++q) {
    const double p99 = pct(e[q], 0.99), mx = e[q].back();
    std::printf("  согласованность %-22s p99 = %.1e, max = %.1e\n", en[q], p99, mx);
    CHECK(p99 < 1e-6, "тождество %s нарушено (p99=%.1e)", en[q], p99);
  }

  // 3. обращение по (P,H) и (P,U)
  int itmax = 0; double itsum = 0, terr = 0; int cnt = 0;
  Hint hint;
  for (int mode = 0; mode < 2; ++mode)
    for (int n = 0; n < 5000; ++n) {
      double P, T; rp(P, T);
      Result f, r;
      m.PT(P, T, HP_ALL, f);
      const int er = mode == 0 ? m.PH(P, f.H, HP_ALL, r, (n & 1) ? &hint : nullptr) : m.PU(P, f.U, HP_ALL, r, (n & 1) ? &hint : nullptr);
      if (er != HE_OK) { CHECK(false, "обращение не удалось: код %d (P=%g, T=%g)", er, P, T); continue; }
      terr = std::fmax(terr, std::fabs(r.T / T - 1));
      itmax = std::max(itmax, r.iters); itsum += r.iters; ++cnt;
      for (unsigned b : {HP_RHO, HP_CP, HP_LAM, HP_DRDP}) if (std::fabs(get(r, b) / get(f, b) - 1) > 1e-9) CHECK(false, "PH/PU: %s не совпала", kNm[__builtin_ctz(b)]);
    }
  std::printf("  обращение PH/PU: макс. |dT/T| = %.1e, итераций ср. %.2f, макс. %d\n", terr, itsum / cnt, itmax);
  CHECK(terr < 1e-10, "точность обращения по T");

  // 4. граничные случаи и экстраполяция
  Result r;
  CHECK(m.PT(NAN, 300, HP_ALL, r) == HE_BAD_INPUT, "NaN по P");
  CHECK(m.PT(1e5, -5, HP_ALL, r) == HE_BAD_INPUT, "T < 0");
  CHECK(m.PT(1e5, 1e5, HP_ALL, r) == HE_BAD_INPUT, "T вне жестких границ");
  CHECK(m.PH(1e5, 1e12, HP_ALL, r) == HE_NO_SOLUTION, "H слишком велика -> NO_SOLUTION");
  CHECK(m.PH(1e5, 1e3, HP_ALL, r) == HE_NO_SOLUTION, "H слишком мала -> NO_SOLUTION");
  CHECK(m.PU(1e5, -1, HP_ALL, r) == HE_BAD_INPUT, "U < 0");
  struct { double P, T; unsigned fl; } ex[] = {{1e6, 40, FL_LOW_T}, {1e6, 3500, FL_HIGH_T}, {5, 300, FL_LOW_P}, {2e8, 300, FL_HIGH_P}, {5, 3500, FL_LOW_P | FL_HIGH_T}};
  for (auto& x : ex) {
    const int er = m.PT(x.P, x.T, HP_ALL, r);
    CHECK(er == HE_OK && r.flags == x.fl && std::isfinite(r.rho) && r.rho > 0 && std::isfinite(r.mu), "экстраполяция P=%g T=%g: код %d, флаги %u", x.P, x.T, er, r.flags);
  }
  // непрерывность через границу диапазона
  double jump = 0;
  for (double P : {10.0, 1e3, 1e5, 1e7, 1e8}) {
    Result a, b;
    m.PT(P, T_MIN, HP_ALL, a); m.PT(P, T_MIN * 0.999, HP_ALL, b);
    jump = std::fmax(jump, std::fabs(a.rho / b.rho - 1));
    m.PT(P, T_MAX, HP_ALL, a); m.PT(P, T_MAX * 1.001, HP_ALL, b);
    jump = std::fmax(jump, std::fabs(a.rho / b.rho - 1));
  }
  for (double T : {50.0, 300.0, 3000.0}) {
    Result a, b;
    m.PT(P_MIN, T, HP_ALL, a); m.PT(P_MIN * 0.999, T, HP_ALL, b);
    jump = std::fmax(jump, std::fabs(a.rho / b.rho - 1));
    m.PT(P_MAX, T, HP_ALL, a); m.PT(P_MAX * 1.001, T, HP_ALL, b);
    jump = std::fmax(jump, std::fabs(a.rho / b.rho - 1));
  }
  std::printf("  граничные случаи и экстраполяция: проверены; скачок rho на границе при сдвиге 0.1 %% = %.1e (допустимо ~1e-3)\n", jump);
  CHECK(jump < 2e-3, "разрыв при выходе за границу");
}

static void speed(const Model& m, const char* tag) {
  std::mt19937_64 rng(7);
  const int N = 400000;
  std::vector<double> P(N), T(N), H(N), U(N), Ps(N), Ts(N);
  std::uniform_real_distribution<double> ux(std::log(60.0), std::log(2900.0)), uy(std::log(20.0), std::log(9e7));
  Result r;
  for (int i = 0; i < N; ++i) {
    T[i] = std::exp(ux(rng)); P[i] = std::exp(uy(rng));
    m.PT(P[i], T[i], HP_H | HP_U, r); H[i] = r.H; U[i] = r.U;
    Ps[i] = 5e6; Ts[i] = 300.0 + 700.0 * i / N;  // "канал": плавно меняющееся состояние
  }
  struct Case { const char* name; int kind; unsigned mask; bool seq; };
  Case cs[] = {{"PT   все 13 (случайные точки)", 0, HP_ALL, false}, {"PT   только rho", 0, HP_RHO, false},
               {"PT   rho,h,cp (случайные)", 0, HP_RHO | HP_H | HP_CP, false}, {"PT   все 13 (канал, подсказка)", 0, HP_ALL, true},
               {"PH   все 13 (случайные)", 1, HP_ALL, false}, {"PH   только T,rho (канал)", 1, HP_RHO, true},
               {"PU   все 13 (случайные)", 2, HP_ALL, false}, {"PU   только T,rho (канал)", 2, HP_RHO, true}};
  std::printf("\n--- скорость, %s (нс на вызов, 1 поток) ---\n", tag);
  double sink = 0;
  for (auto& c : cs) {
    std::vector<double> Hs(N), Us(N);
    if (c.seq) for (int i = 0; i < N; ++i) { m.PT(Ps[i], Ts[i], HP_H | HP_U, r); Hs[i] = r.H; Us[i] = r.U; }
    const double* pp = c.seq ? Ps.data() : P.data();
    const double* xx = c.kind == 0 ? (c.seq ? Ts.data() : T.data()) : c.kind == 1 ? (c.seq ? Hs.data() : H.data()) : (c.seq ? Us.data() : U.data());
    Hint h;
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < N; ++i) {
      if (c.kind == 0) m.PT(pp[i], xx[i], c.mask, r, &h); else if (c.kind == 1) m.PH(pp[i], xx[i], c.mask, r, &h); else m.PU(pp[i], xx[i], c.mask, r, &h);
      sink += r.rho;
    }
    const double ns = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / N * 1e9;
    std::printf("  %-34s %8.0f нс\n", c.name, ns);
  }
  if (sink == 12345.678) std::printf(" ");
}

int main(int argc, char** argv) {
  if (argc < 2) { std::fprintf(stderr, "Использование: %s table.bin [coef.txt]\n", argv[0]); return 1; }
  std::string err;
  auto tab = TableModel::open(argv[1], err);
  if (!tab) { std::fprintf(stderr, "Ошибка: %s\n", err.c_str()); return 2; }
  AnalyticModel an;
  if (argc > 2 && !an.load(argv[2], &err)) { std::fprintf(stderr, "Ошибка: %s\n", err.c_str()); return 2; }
  std::printf("Таблица: %zu листьев, tol=%.1e, эталон: %s\n", tab->n_leaves(), tab->header().tol, tab->header().ref_name);

  run_model(*tab, "A");
  run_model(an, "B");

  // согласие A и B (A построена по B, поэтому разность = ошибка интерполяции)
  std::printf("\n=== A против B (отклонение в области, 100000 точек) ===\n");
  std::mt19937_64 rng(99);
  std::uniform_real_distribution<double> ux(std::log(T_MIN), std::log(T_MAX)), uy(std::log(P_MIN), std::log(P_MAX));
  double mx[14] = {0};
  for (int n = 0; n < 100000; ++n) {
    const double T = std::exp(ux(rng)), P = std::exp(uy(rng));
    Result a, b;
    tab->PT(P, T, HP_ALL, a); an.PT(P, T, HP_ALL, b);
    for (int k = 1; k < 14; ++k) mx[k] = std::fmax(mx[k], std::fabs(get(a, kBits[k]) / get(b, kBits[k]) - 1));
  }
  double w = 0;
  for (int k = 1; k < 14; ++k) { std::printf("  %-7s %.2e\n", kNm[k], mx[k]); w = std::fmax(w, mx[k]); }
  CHECK(w < 5 * tab->header().tol + 1e-12 || tab->header().ref_name[0] != 'a', "A и B расходятся сильнее допуска (%.2e)", w);

  speed(*tab, "A: таблица");
  speed(an, "B: аналитическая модель");
  std::printf("\n%s\n", g_fail ? "ИТОГ: ЕСТЬ НЕУДАЧИ" : "ИТОГ: все проверки пройдены");
  return g_fail ? 1 : 0;
}
