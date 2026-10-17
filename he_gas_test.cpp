// he_gas_test.cpp - проверка he_gas.* (ASPIG/ASPIGU/ASPIGH/AGHmax/AGUmax).
//   ./he_gas_test [таблица.bin [коэффициенты.txt]]     - основной набор (на реальной модели)
//   ./he_gas_test --fail                               - модель не загружается: err = 1, сообщение непустое
// Эталон - прямые вызовы he::Model (AGGetHeModel()); плюс тождества (beta, k, Pr, обращения).
#include <cmath>
#include <cstdio>
#include <cstring>
#include <thread>
#include <vector>
#include "he_api.h"
#include "he_gas.h"
#ifdef HE_GAS_TEST_MOCK
#include "he_eos.h"
#endif
#include "he_console.h"

static int g_fail = 0, g_n = 0;
#define CHECK(c, ...) do { ++g_n; if (!(c)) { ++g_fail; std::printf("FAIL %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)
static double rel(double a, double b) { return std::fabs(a - b) / std::fmax(std::fabs(b), 1e-300); }

static const double PS[] = {10, 1e3, 1e5, 1e6, 1e7, 1e8};
static const double TS[] = {50, 99.99, 100.01, 300, 1000, 3000};

int main(int argc, char** argv) {
  if (argc > 1 && !std::strcmp(argv[1], "--fail")) {
    ASPIGInit("no_such_table.bin", "no_such_coef.txt");
    double ff[AG_NFF]; int nff[AG_NPROP], err = 0;
    for (int i = 0; i < AG_NPROP; ++i) nff[i] = 1;
    ff[0] = -777;
    ASPIG(4, 1e5, 300, ff, nff, &err);
    CHECK(err == 1, "err=%d", err);
    CHECK(ff[0] == -777, "ff[0] изменен при ошибке");
    CHECK(AGLastError()[0] != 0, "пустое сообщение об ошибке");
    CHECK(ff[AG_MOL] > 0 && ff[AG_R] > 0, "ff[47..49] должны заполняться всегда");
    double a, b; AGHmax(4, &a, &b);
    CHECK(a > 0 && b > a, "запасные границы H: %g %g", a, b);
    std::printf("--fail: %d проверок, ошибок %d (%s)\n", g_n, g_fail, AGLastError());
    return g_fail != 0;
  }
#ifdef HE_GAS_TEST_MOCK
  ASPIGInit("none.bin", "mock-coef");
#else
  if (argc > 1) ASPIGInit(argv[1], argc > 2 ? argv[2] : nullptr);
#endif
  const he::Model* m = AGGetHeModel();
  if (!m) { std::printf("модель не загружена: %s\n", AGLastError()); return 2; }
  std::printf("модель: %s\n", AGModelInfo());

  int all[AG_NPROP]; for (int i = 0; i < AG_NPROP; ++i) all[i] = 1;
  double ff[AG_NFF], g[AG_NFF];
  int err;

  // 1. ASPIG против прямых вызовов Model + тождества
  for (double P : PS) for (double T : TS) {
    for (int i = 0; i < AG_NFF; ++i) ff[i] = -777;
    ASPIG(4, P, T, ff, all, &err);
    he::Result r;
    CHECK(err == 0, "err P=%g T=%g", P, T);
    CHECK(m->PT(P, T, he::HP_ALL, r) == he::HE_OK, "model");
    CHECK(ff[AG_RHO] == r.rho && ff[AG_CP] == r.cp && ff[AG_MU] == r.mu && ff[AG_LAM] == r.lambda, "rho/cp/mu/lam");
    CHECK(ff[AG_U] == r.U && ff[AG_H] == r.H && ff[AG_CV] == r.cv && ff[AG_PR] == r.Pr, "U/H/cv/Pr");
    CHECK(ff[AG_DUDT] == r.dUdT && ff[AG_DRDT] == r.drhodT && ff[AG_DRDP] == r.drhodP, "dUdT/drdT/drdP");
    CHECK(ff[AG_DRDU] == r.drhodU && ff[AG_DRDH] == r.drhodH && ff[AG_T] == T, "drdU/drdH/T");
    CHECK(rel(ff[AG_BETA], -r.drhodT / r.rho) < 1e-14, "beta");
    CHECK(rel(ff[AG_K], r.cp / r.cv) < 1e-14, "k");
    CHECK(rel(ff[AG_R], he::R_SPEC) < 1e-14 && rel(ff[AG_MOL], 4.002602e-3) < 1e-14, "Mol/R");
    CHECK(rel(ff[AG_PR], ff[AG_CP] * ff[AG_MU] / ff[AG_LAM]) < 1e-12, "Pr = cp mu/lambda");
    CHECK(rel(ff[AG_DRDH], ff[AG_DRDT] / ff[AG_CP]) < 1e-12, "drdH = drdT/cp");
    CHECK(rel(ff[AG_DRDU], ff[AG_DRDT] / ff[AG_DUDT]) < 1e-12, "drdU = drdT/dUdT");
    // 2. обращение: (P,H) и (P,U) возвращают T
    ASPIGH(4, P, ff[AG_H], g, all, &err);
    CHECK(err == 0 && rel(g[AG_T], T) < 1e-9, "PH: T=%.12g -> %.12g (err %d)", T, g[AG_T], err);
    CHECK(rel(g[AG_RHO], ff[AG_RHO]) < 1e-9 && rel(g[AG_MU], ff[AG_MU]) < 1e-9, "PH rho/mu");
    ASPIGU(4, P, ff[AG_U], g, all, &err);
    CHECK(err == 0 && rel(g[AG_T], T) < 1e-9, "PU: T=%.12g -> %.12g (err %d)", T, g[AG_T], err);
  }

  // 3. маска nff: пишутся только запрошенные; ff[47..49] - всегда
  {
    int one[AG_NPROP] = {0}; one[AG_MU] = 1;
    for (int i = 0; i < AG_NFF; ++i) g[i] = -777;
    ASPIG(4, 1e5, 300, g, one, &err);
    CHECK(err == 0, "mask err");
    for (int i = 0; i < AG_NPROP; ++i) CHECK((g[i] == -777) == (i != AG_MU), "mask i=%d", i);
    CHECK(g[AG_MOL] > 0 && g[AG_K] > 1 && g[AG_R] > 0, "ff[47..49] при пустой маске");
    int none[AG_NPROP] = {0};
    for (int i = 0; i < AG_NFF; ++i) g[i] = -777;
    ASPIG(4, 1e5, 300, g, none, &err);
    int touched = 0; for (int i = 0; i < AG_NPROP; ++i) touched += (g[i] != -777);
    CHECK(err == 0 && touched == 0 && g[AG_K] > 1, "пустая маска: err %d, тронуто %d", err, touched);
    for (int i = 0; i < AG_NFF; ++i) g[i] = -777;
    ASPIG(4, 1e5, 300, g, nullptr, &err);  // NULL = все
    CHECK(err == 0 && g[0] != -777 && g[AG_T] == 300, "nff = NULL");
  }

  // 4. ошибочные входы: ff[0..14] не меняются, err = 1, ff[47..49] заполнены
  {
    const double bad[][2] = {{-1, 300}, {1e5, -5}, {NAN, 300}, {1e5, NAN}, {0, 300}, {1e5, 1e7}, {1e12, 300}};
    for (auto& b : bad) {
      for (int i = 0; i < AG_NFF; ++i) g[i] = -777;
      ASPIG(4, b[0], b[1], g, all, &err);
      CHECK(err == 1 && g[0] == -777 && g[AG_MOL] > 0, "bad P=%g T=%g err=%d", b[0], b[1], err);
    }
    ASPIGH(4, 1e5, 1e30, g, all, &err); CHECK(err == 1, "H вне диапазона -> err=1 (%d)", err);
    ASPIGU(4, 1e5, -1e30, g, all, &err); CHECK(err == 1, "U вне диапазона -> err=1 (%d)", err);
    ASPIGH(4, NAN, 1e6, g, all, &err); CHECK(err == 1, "P=NaN");
  }

  // 5. экстраполяция и строгий режим
  {
    ASPIG(4, 1e5, 40, g, all, &err);
    CHECK(err == 0 && (AGLastFlags() & he::FL_LOW_T), "T=40 К: err=%d flags=%d", err, AGLastFlags());
    ASPIG(4, 1e5, 300, g, all, &err);
    CHECK(err == 0 && AGLastFlags() == 0, "флаги сбрасываются");
    AGSetStrict(1);
    ASPIG(4, 1e5, 40, g, all, &err); CHECK(err == 1, "strict: T=40 -> err=1");
    ASPIG(4, 1e5, 300, g, all, &err); CHECK(err == 0, "strict: T=300 -> err=0");
    AGSetStrict(0);
  }

  // 6. границы H, U и согласованность с моделью
  {
    double h0, h1, u0, u1; AGHmax(4, &h0, &h1); AGUmax(4, &u0, &u1);
    CHECK(h0 < h1 && u0 < u1, "порядок границ");
    double hmin = 1e300, hmax = -1e300, umin = 1e300, umax = -1e300;
    for (double T : {50.0, 3000.0}) for (int k = 0; k <= 70; ++k) {
      double P = 10 * std::pow(1e7, k / 70.0); he::Result r;
      m->PT(P, T, he::HP_H | he::HP_U, r);
      hmin = std::fmin(hmin, r.H); hmax = std::fmax(hmax, r.H); umin = std::fmin(umin, r.U); umax = std::fmax(umax, r.U);
    }
    CHECK(rel(h0, hmin) < 1e-3 && rel(h1, hmax) < 1e-3, "H: [%g,%g] vs плотная сетка [%g,%g]", h0, h1, hmin, hmax);
    CHECK(rel(u0, umin) < 1e-3 && rel(u1, umax) < 1e-3, "U: [%g,%g] vs [%g,%g]", u0, u1, umin, umax);
    // границы пригодны как вход: ASPIGH на границах при допустимом P находит T внутри 50-3000 К (экстраполяция допустима)
    ASPIGH(4, 1e5, h0 * 0.5 + h1 * 0.5, g, all, &err); CHECK(err == 0, "середина диапазона H");
    std::printf("Hmin=%.6g Hmax=%.6g Umin=%.6g Umax=%.6g\n", h0, h1, u0, u1);
    double a = 1, b = 1; AGHmax(1, &a, &b); CHECK(a == 0 && b == 0, "заглушка AGHmax");
    AGUmax(7, &a, &b); CHECK(a == 0 && b == 0, "заглушка AGUmax");
    AGHmax(4, nullptr, nullptr);  // не падает
  }

  // 7. заглушки других газов
  {
    const int ngs[] = {1, 2, 3, 5, 6, 7, 8, 0, 9, -3, 99};
    const double mol[] = {28.9647e-3, 28.0134e-3, 31.9988e-3, 2.01588e-3, 44.0095e-3, 39.948e-3, 131.293e-3, 28.9647e-3, 28.9647e-3, 28.9647e-3, 28.9647e-3};
    for (size_t i = 0; i < sizeof ngs / sizeof *ngs; ++i) {
      for (int k = 0; k < AG_NFF; ++k) g[k] = -777;
      ASPIG(ngs[i], 1e5, 300, g, all, &err);
      CHECK(err == 1 && g[0] == -777 && rel(g[AG_MOL], mol[i]) < 1e-12 && rel(g[AG_R], 8.314462618 / mol[i]) < 1e-12 && g[AG_K] > 1.2,
            "заглушка Ng=%d: err=%d Mol=%g", ngs[i], err, g[AG_MOL]);
      ASPIGH(ngs[i], 1e5, 1e6, g, all, &err); CHECK(err == 1, "заглушка PH Ng=%d", ngs[i]);
      ASPIGU(ngs[i], 1e5, 1e6, g, all, &err); CHECK(err == 1, "заглушка PU Ng=%d", ngs[i]);
    }
    ASPIG(4, 1e5, 300, nullptr, all, &err); CHECK(err == 1, "ff = NULL");
    ASPIG(4, 1e5, 300, g, all, nullptr);  // err = NULL: не падает
  }

  // 8. ASPIGInit после загрузки
  CHECK(ASPIGInit("x", "y") == -1, "повторный ASPIGInit");

  // 9. потоки: одинаковые результаты, нет гонок
  {
    ASPIG(4, 1e7, 500, g, all, &err);
    std::vector<std::thread> th; std::vector<int> bad(8, 0);
    for (int t = 0; t < 8; ++t) th.emplace_back([&, t] {
      double q[AG_NFF]; int e;
      for (int i = 0; i < 20000; ++i) {
        ASPIG(4, 1e7, 500, q, all, &e);
        if (e || q[AG_RHO] != g[AG_RHO] || q[AG_MU] != g[AG_MU]) { bad[t] = 1; break; }
        ASPIGH(4, 1e7, g[AG_H], q, all, &e);
        if (e || rel(q[AG_T], 500) > 1e-9) { bad[t] = 1; break; }
      }
    });
    for (auto& x : th) x.join();
    int nb = 0; for (int b : bad) nb += b;
    CHECK(nb == 0, "потоки: %d из 8 дали расхождение", nb);
  }

#ifdef HE_GAS_TEST_MOCK
  // 10. замкнутые формулы для стендовой модели (идеальный газ + постоянный B)
  {
    const double R = he::R_SPEC, B = he::MOCK_B, P = 1e7, T = 400;
    ASPIG(4, P, T, g, all, &err);
    const double rho = 1 / (R * T / P + B);
    CHECK(rel(g[AG_RHO], rho) < 1e-14, "rho");
    CHECK(rel(g[AG_CP], 2.5 * R) < 1e-14 && rel(g[AG_CV], 1.5 * R) < 1e-12, "cp, cv");
    CHECK(rel(g[AG_H], 2.5 * R * T + B * P) < 1e-14, "H");
    CHECK(rel(g[AG_U], 2.5 * R * T + B * P - P / rho) < 1e-12, "U");
    CHECK(rel(g[AG_BETA], (R / P) * rho) < 1e-12, "beta = rho R/P");
    CHECK(rel(g[AG_DRDP], rho * rho * R * T / (P * P)) < 1e-12, "drho/dP");
    CHECK(rel(g[AG_K], 5.0 / 3) < 1e-12, "k = 5/3");
    double h0, h1; AGHmax(4, &h0, &h1);
    CHECK(rel(h0, 2.5 * R * 50 + B * 10) < 1e-12 && rel(h1, 2.5 * R * 3000 + B * 1e8) < 1e-12, "Hmin/Hmax замкнуто: %.9g %.9g", h0, h1);
  }
#endif
  std::printf("he_gas_test: %d проверок, ошибок %d\n", g_n, g_fail);
  return g_fail != 0;
}
