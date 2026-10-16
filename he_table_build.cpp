// he_table_build.cpp - Программа 1 (версия 2): построение таблицы остаточного потенциала Гиббса
// G = ln(phi) на квадродереве с адаптивным измельчением, контроль точности и согласованности,
// запись бинарной БД (заголовок с порядком байтов, версией и контрольной суммой).
//
// Алгоритм:
//  1. Область x = ln T, y = ln P делится на rx*ry корневых ячеек; каждая ячейка (лист) хранит
//     биквинтический эрмитов полином G (36 коэффициентов) и бикубические полиномы ln(lambda), ln(mu).
//  2. Данные в узлах: G и ее производные до порядка (2,2). Производные до второго порядка
//     берутся из эталона аналитически (Z, h, cp, drho/dT, drho/dP); третьего и четвертого порядка –
//     разностями этих величин (шаг --step в координатах ln T, ln P).
//  3. Лист проверяется на 9 точках по ВСЕМ 13 выходным величинам (включая производные и Pr);
//     при ошибке > tol делится на 4 (до --maxlevel).
//  4. Итоговый контроль: случайные точки, скачки на стыках листьев, самосогласованность эталона.
//
// Сборка: g++ -O2 -std=c++17 he_table_build.cpp he_ref.cpp he_table.cpp he_eos.cpp he_api.cpp
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <random>
#include <unordered_map>

#include "he_eos.h"
#include "he_ref.h"
#include "he_table.h"
#include "he_console.h"

using namespace he;

namespace {

constexpr int MAXL = 12;  // разрешение целочисленных координат узлов
constexpr int NOUT = 13;
const char* kNames[NOUT] = {"rho", "cp", "cv", "U", "H", "lambda", "mu", "Pr", "dU/dT", "drho/dT", "drho/dU", "drho/dH", "drho/dP"};

struct NodeData { double d[9]; double t[2][4]; };

double X0, X1, Y0, Y1, g_step = 0.02;
uint32_t RX = 2, RY = 8;
std::unordered_map<uint64_t, NodeData> g_cache;

[[noreturn]] void die(const char* msg, double T, double P) {
  std::fprintf(stderr, "ОШИБКА: %s (T=%g К, P=%g Па)\n", msg, T, P);
  std::exit(2);
}

NodeData node_at(uint32_t ix, uint32_t iy) {
  const uint64_t key = ((uint64_t)ix << 32) | iy;
  auto it = g_cache.find(key);
  if (it != g_cache.end()) return it->second;
  const double x = X0 + (X1 - X0) * ix / (double)(RX << MAXL), y = Y0 + (Y1 - Y0) * iy / (double)(RY << MAXL);
  const double R = R_SPEC, h = g_step;
  double gx[3][3], gy[3][3], gxx[3][3], gxy[3][3], gyy[3][3], lm[3][3], ll[3][3];
  double G0 = 0;
  for (int i = -1; i <= 1; ++i)
    for (int j = -1; j <= 1; ++j) {
      const double T = std::exp(x + i * h), P = std::exp(y + j * h);
      RefState r;
      if (!ref_eval(T, P, r) || !(r.rho > 0) || !(r.mu > 0) || !(r.lam > 0)) die("эталонная модель не дала результата", T, P);
      const double RT = R * T, Z = P / (r.rho * RT);
      const int a = i + 1, b = j + 1;
      gy[a][b] = Z - 1.0;
      gx[a][b] = -(r.H - 2.5 * RT) / RT;
      gxx[a][b] = -(r.cp - 2.5 * R) / R - gx[a][b];
      const double dvdT = -r.drdT / (r.rho * r.rho), dvdP = -r.drdP / (r.rho * r.rho);
      gxy[a][b] = (P / R) * dvdT - (1.0 + gy[a][b]);
      gyy[a][b] = (P * P / RT) * dvdP + 1.0 + gy[a][b];
      lm[a][b] = std::log(r.mu);
      ll[a][b] = std::log(r.lam);
      if (i == 0 && j == 0) G0 = r.lnphi;
    }
  NodeData n;
  const double h2 = 2 * h, hh = h * h;
  n.d[0] = G0;
  n.d[1] = gx[1][1];
  n.d[2] = gy[1][1];
  n.d[3] = gxx[1][1];
  n.d[4] = gxy[1][1];
  n.d[5] = gyy[1][1];
  n.d[6] = 0.5 * ((gxx[1][2] - gxx[1][0]) / h2 + (gxy[2][1] - gxy[0][1]) / h2);                    // G_xxy
  n.d[7] = 0.5 * ((gxy[1][2] - gxy[1][0]) / h2 + (gyy[2][1] - gyy[0][1]) / h2);                    // G_xyy
  n.d[8] = (1.0 / 3.0) * ((gxy[2][2] - gxy[2][0] - gxy[0][2] + gxy[0][0]) / (4 * hh) +             // G_xxyy
                          (gyy[2][1] - 2 * gyy[1][1] + gyy[0][1]) / hh +
                          (gxx[1][2] - 2 * gxx[1][1] + gxx[1][0]) / hh);
  const double* tf[2] = {&ll[0][0], &lm[0][0]};
  for (int k = 0; k < 2; ++k) {
    auto F = [&](int a, int b) { return tf[k][a * 3 + b]; };
    n.t[k][0] = F(1, 1);
    n.t[k][1] = (F(2, 1) - F(0, 1)) / h2;
    n.t[k][2] = (F(1, 2) - F(1, 0)) / h2;
    n.t[k][3] = (F(2, 2) - F(2, 0) - F(0, 2) + F(0, 0)) / (4 * hh);
  }
  g_cache.emplace(key, n);
  return n;
}

// Эрмитовы базисы в мономиальной форме (строки: f0, f0', f0'', f1, f1', f1'')
const double M6[6][6] = {{1, 0, 0, -10, 15, -6}, {0, 1, 0, -6, 8, -3}, {0, 0, 0.5, -1.5, 1.5, -0.5},
                         {0, 0, 0, 10, -15, 6},  {0, 0, 0, -4, 7, -3},  {0, 0, 0, 0.5, -1, 0.5}};
const double M4[4][4] = {{1, 0, -3, 2}, {0, 1, -2, 1}, {0, 0, 3, -2}, {0, 0, -1, 1}};
const int IDX3[3][3] = {{0, 2, 5}, {1, 4, 7}, {3, 6, 8}};  // [порядок по x][порядок по y] -> номер в NodeData::d
const int IDX2[2][2] = {{0, 2}, {1, 3}};

struct LeafId { uint32_t l, i, j; };

LeafRec build_leaf(const LeafId& id) {
  const double hx = (X1 - X0) / (double)(RX << id.l), hy = (Y1 - Y0) / (double)(RY << id.l);
  const uint32_t sh = MAXL - id.l;
  NodeData nd[2][2];
  for (int a = 0; a < 2; ++a)
    for (int b = 0; b < 2; ++b) nd[a][b] = node_at((id.i + a) << sh, (id.j + b) << sh);
  LeafRec r;
  std::memset(&r, 0, sizeof r);
  r.level = id.l; r.i = id.i; r.j = id.j;
  double D[6][6], T1[6][6];
  for (int cx = 0; cx < 2; ++cx)
    for (int ox = 0; ox < 3; ++ox)
      for (int cy = 0; cy < 2; ++cy)
        for (int oy = 0; oy < 3; ++oy)
          D[cx * 3 + ox][cy * 3 + oy] = nd[cx][cy].d[IDX3[ox][oy]] * std::pow(hx, ox) * std::pow(hy, oy);
  for (int a = 0; a < 6; ++a)
    for (int j = 0; j < 6; ++j) { double s = 0; for (int b = 0; b < 6; ++b) s += D[a][b] * M6[b][j]; T1[a][j] = s; }
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) { double s = 0; for (int a = 0; a < 6; ++a) s += M6[a][i] * T1[a][j]; r.c[6 * i + j] = s; }
  for (int k = 0; k < 2; ++k) {
    double E[4][4], F1[4][4];
    for (int cx = 0; cx < 2; ++cx)
      for (int ox = 0; ox < 2; ++ox)
        for (int cy = 0; cy < 2; ++cy)
          for (int oy = 0; oy < 2; ++oy)
            E[cx * 2 + ox][cy * 2 + oy] = nd[cx][cy].t[k][IDX2[ox][oy]] * std::pow(hx, ox) * std::pow(hy, oy);
    for (int a = 0; a < 4; ++a)
      for (int j = 0; j < 4; ++j) { double s = 0; for (int b = 0; b < 4; ++b) s += E[a][b] * M4[b][j]; F1[a][j] = s; }
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) { double s = 0; for (int a = 0; a < 4; ++a) s += M4[a][i] * F1[a][j]; r.tr[k][4 * i + j] = s; }
  }
  return r;
}

LeafView view_of(const LeafRec& r) {
  const double hx = (X1 - X0) / (double)(RX << r.level), hy = (Y1 - Y0) / (double)(RY << r.level);
  return {X0 + r.i * hx, Y0 + r.j * hy, 1.0 / hx, 1.0 / hy, r.c, &r.tr[0][0]};
}

void eval_in_leaf(const LeafView& lv, double x, double y, Result& out) {
  const double T = std::exp(x), P = std::exp(y);
  Thermo s;
  leaf_thermo(lv, x, y, 0, 0, T, P, N_RHO | N_H | N_CP | N_DV, s);
  double lam = 0, mu = 0;
  leaf_transport(lv, x, y, 0, 0, HP_LAM | HP_MU | HP_PR, lam, mu);
  assemble(P, T, HP_ALL, s, lam, mu, out);
}

void outs(const Result& r, double a[NOUT]) {
  const double v[NOUT] = {r.rho, r.cp, r.cv, r.U, r.H, r.lambda, r.mu, r.Pr, r.dUdT, r.drhodT, r.drhodU, r.drhodH, r.drhodP};
  std::memcpy(a, v, sizeof v);
}

void ref_outs(double T, double P, double b[NOUT]) {
  RefState q;
  if (!ref_eval(T, P, q)) die("эталонная модель не дала результата", T, P);
  const double v[NOUT] = {q.rho, q.cp, q.cv, q.U, q.H, q.lam, q.mu, q.cp * q.mu / q.lam, q.dUdT, q.drdT, q.drdU, q.drdH, q.drdP};
  std::memcpy(b, v, sizeof v);
}

double max_rel(const double a[NOUT], const double b[NOUT], double e[NOUT] = nullptr) {
  double m = 0;
  for (int k = 0; k < NOUT; ++k) {
    const double d = std::fabs(a[k] - b[k]) / std::fabs(b[k]);
    if (e && d > e[k]) e[k] = d;
    if (d > m) m = d;
  }
  return m;
}

double leaf_error(const LeafRec& rec) {
  static const double uv[9][2] = {{.5, .5}, {.5, 0}, {.5, 1}, {0, .5}, {1, .5}, {.25, .25}, {.75, .25}, {.25, .75}, {.75, .75}};
  const LeafView lv = view_of(rec);
  const double hx = 1.0 / lv.ihx, hy = 1.0 / lv.ihy;
  double worst = 0;
  for (auto& p : uv) {
    const double x = lv.x0 + p[0] * hx, y = lv.y0 + p[1] * hy;
    Result r;
    eval_in_leaf(lv, x, y, r);
    double a[NOUT], b[NOUT];
    outs(r, a);
    ref_outs(std::exp(x), std::exp(y), b);
    worst = std::fmax(worst, max_rel(a, b));
  }
  return worst;
}

// Самосогласованность эталона: тождества термодинамики, проверяемые конечными разностями
void selfcheck(int n) {
  std::mt19937_64 rng(777);
  std::uniform_real_distribution<double> ux(X0, X1), uy(Y0, Y1);
  double e[8] = {0};
  const char* nm[8] = {"cp = dH/dT", "drho/dT", "dU/dT", "drho/dP", "Maxwell dH/dP (от v)", "cv = cp+T dvdT^2/dvdP", "drho/dU", "drho/dH"};
  for (int k = 0; k < n; ++k) {
    const double T = std::exp(ux(rng)), P = std::exp(uy(rng)), hT = 1e-4 * T, hP = 1e-4 * P;
    RefState c, a, b, pp, pm;
    if (!ref_eval(T, P, c) || !ref_eval(T + hT, P, a) || !ref_eval(T - hT, P, b) || !ref_eval(T, P + hP, pp) || !ref_eval(T, P - hP, pm)) die("selfcheck: сбой эталона", T, P);
    auto rel = [](double x, double y) { return std::fabs(x - y) / std::fabs(y); };
    const double dvdT = -c.drdT / (c.rho * c.rho), dvdP = -c.drdP / (c.rho * c.rho);
    const double v[8] = {rel((a.H - b.H) / (2 * hT), c.cp), rel((a.rho - b.rho) / (2 * hT), c.drdT), rel((a.U - b.U) / (2 * hT), c.dUdT),
                         rel((pp.rho - pm.rho) / (2 * hP), c.drdP), std::fabs((pp.H - pm.H) / (2 * hP) - (1.0 / c.rho - T * dvdT)) * c.rho,
                         rel(c.cv, c.cp + T * dvdT * dvdT / dvdP), rel(c.drdU, c.drdT / c.dUdT), rel(c.drdH, c.drdT / c.cp)};
    for (int q = 0; q < 8; ++q) e[q] = std::fmax(e[q], v[q]);
  }
  std::printf("Самосогласованность эталона (%d точек, макс. относительное расхождение):\n", n);
  double w = 0;
  for (int q = 0; q < 8; ++q) { std::printf("  %-24s %.2e\n", nm[q], e[q]); w = std::fmax(w, e[q]); }
  if (w > 1e-4) std::printf("  ПРЕДУПРЕЖДЕНИЕ: эталон термодинамически несогласован (>1e-4); таблица унаследует это.\n");
}

}  // namespace

int main(int argc, char** argv) {
  std::string out = "he_table2.bin", coef;
  double tol = 1e-4;
  int maxlevel = 8, minlevel = 1, nverify = 200000, nself = 0;
  bool verbose = false;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto nx = [&]() -> const char* { if (i + 1 >= argc) { std::fprintf(stderr, "нет значения для %s\n", a.c_str()); std::exit(1); } return argv[++i]; };
    if (a == "-o") out = nx();
    else if (a == "--tol") tol = std::atof(nx());
    else if (a == "--rx") RX = std::atoi(nx());
    else if (a == "--ry") RY = std::atoi(nx());
    else if (a == "--maxlevel") maxlevel = std::atoi(nx());
    else if (a == "--minlevel") minlevel = std::atoi(nx());
    else if (a == "--step") g_step = std::atof(nx());
    else if (a == "--verify") nverify = std::atoi(nx());
    else if (a == "--coef") coef = nx();
    else if (a == "--selfcheck") nself = std::atoi(nx());
    else if (a == "-v") verbose = true;
    else {
      std::printf("Использование: %s [-o файл] [--tol 1e-4] [--rx 2] [--ry 8] [--minlevel 1] [--maxlevel 8]\n"
                  "       [--step 0.02] [--verify 200000] [--selfcheck N] [--coef файл_коэффициентов] [-v]\n", argv[0]);
      return a == "-h" || a == "--help" ? 0 : 1;
    }
  }
  if (maxlevel > 10) maxlevel = 10;
  if (minlevel > maxlevel) minlevel = maxlevel;
  X0 = std::log(T_MIN); X1 = std::log(T_MAX); Y0 = std::log(P_MIN); Y1 = std::log(P_MAX);
  std::string err;
  if (!ref_init(coef.c_str(), &err)) { std::fprintf(stderr, "ОШИБКА: %s\n", err.c_str()); return 1; }
  std::printf("Эталон: %s\nДиапазон: T %.0f-%.0f К, P %.0f Па - %.0f МПа; tol=%.2g, корневая сетка %ux%u, уровни %d..%d\n",
              ref_name(), T_MIN, T_MAX, P_MIN, P_MAX / 1e6, tol, RX, RY, minlevel, maxlevel);
  const auto t0 = std::chrono::steady_clock::now();
  if (nself > 0) selfcheck(nself);

  std::vector<LeafId> work;
  for (uint32_t i = 0; i < RX; ++i) for (uint32_t j = 0; j < RY; ++j) work.push_back({0, i, j});
  std::vector<LeafRec> leaves;
  int unresolved = 0;
  double worst_leaf = 0;
  size_t done = 0;
  while (!work.empty()) {
    const LeafId id = work.back();
    work.pop_back();
    const LeafRec rec = build_leaf(id);
    double e = 0;
    const bool must = (int)id.l < minlevel;
    if (!must) e = leaf_error(rec);
    if ((must || e > tol) && (int)id.l < maxlevel) {
      for (uint32_t a = 0; a < 2; ++a) for (uint32_t b = 0; b < 2; ++b) work.push_back({id.l + 1, 2 * id.i + a, 2 * id.j + b});
    } else {
      if (e > tol) ++unresolved;
      worst_leaf = std::fmax(worst_leaf, e);
      leaves.push_back(rec);
    }
    if (++done % 200 == 0) {
      std::printf("  обработано ячеек: %zu, в очереди: %zu, листьев: %zu\r", 
                  done, work.size(), leaves.size());
      std::fflush(stdout);
    }    
  }
  std::printf("\n");
  TableHeader h;
  std::memset(&h, 0, sizeof h);
  std::memcpy(h.magic, "HE4TBL02", 8);
  h.endian = 0x01020304u; h.version = 2; h.header_size = sizeof(TableHeader); h.leaf_size = sizeof(LeafRec);
  h.rx = RX; h.ry = RY; h.Tmin = T_MIN; h.Tmax = T_MAX; h.Pmin = P_MIN; h.Pmax = P_MAX; h.tol = tol;
  std::snprintf(h.ref_name, sizeof h.ref_name, "%s", ref_name());
  auto tab = TableModel::from_leaves(h, std::move(leaves));
  std::printf("Листьев: %zu (макс. уровень %u), необработанных (> tol при maxlevel): %d, макс. ошибка листьев %.3e\n",
              tab->n_leaves(), tab->header().max_level, unresolved, worst_leaf);

  // ---- контроль: случайные точки через общий API (так же, как в эксплуатации)
  std::mt19937_64 rng(12345);
  std::uniform_real_distribution<double> ux(X0, X1), uy(Y0, Y1);
  double emax[NOUT] = {0}, esum[NOUT] = {0}, worst = 0, wT = 0, wP = 0;
  Hint hint;
  for (int n = 0; n < nverify + 4; ++n) {
    double x, y;
    if (n < 4) { x = (n & 1) ? X1 : X0; y = (n & 2) ? Y1 : Y0; } else { x = ux(rng); y = uy(rng); }
    const double T = std::exp(x), P = std::exp(y);
    Result r;
    if (tab->PT(P, T, HP_ALL, r, &hint) != HE_OK) die("таблица не дала результата", T, P);
    double a[NOUT], b[NOUT], e[NOUT] = {0};
    outs(r, a);
    ref_outs(T, P, b);
    const double m = max_rel(a, b, e);
    for (int k = 0; k < NOUT; ++k) { esum[k] += e[k]; emax[k] = std::fmax(emax[k], e[k]); }
    if (m > worst) { worst = m; wT = T; wP = P; }
  }
  std::printf("\nКонтроль на %d случайных точках (отклонение таблицы от эталона):\n  %-9s %12s %12s\n", nverify + 4, "величина", "максимум", "среднее");
  for (int k = 0; k < NOUT; ++k) std::printf("  %-9s %12.3e %12.3e\n", kNames[k], emax[k], esum[k] / (nverify + 4));
  std::printf("  наибольшее отклонение %.3e при T=%.2f К, P=%.4g Па\n", worst, wT, wP);

  // ---- скачки на стыках листьев (непрерывность)
  double seam = 0;
  for (size_t k = 0; k < tab->n_leaves(); ++k) {
    const LeafView lv = tab->view((int)k);
    const double hx = 1.0 / lv.ihx, hy = 1.0 / lv.ihy;
    for (int side = 0; side < 4; ++side)
      for (double t : {0.25, 0.5, 0.75}) {
        double u = side == 0 ? 0.0 : side == 1 ? 1.0 : t, v = side == 2 ? 0.0 : side == 3 ? 1.0 : t;
        const double x = lv.x0 + u * hx, y = lv.y0 + v * hy;
        const double xo = x + (side == 0 ? -1e-7 * hx : side == 1 ? 1e-7 * hx : 0), yo = y + (side == 2 ? -1e-7 * hy : side == 3 ? 1e-7 * hy : 0);
        if (xo < tab->X0() || xo > tab->X1() || yo < tab->Y0() || yo > tab->Y1()) continue;
        const int k2 = tab->locate(xo, yo, nullptr);
        if (k2 == (int)k) continue;
        Result r1, r2;
        eval_in_leaf(lv, x, y, r1);
        eval_in_leaf(tab->view(k2), x, y, r2);
        double a[NOUT], b[NOUT];
        outs(r1, a); outs(r2, b);
        seam = std::fmax(seam, max_rel(a, b));
      }
  }
  std::printf("  макс. скачок выходных величин на стыках листьев: %.3e\n", seam);
  if (worst > 0.01) std::printf("  ПРЕДУПРЕЖДЕНИЕ: ошибка таблицы превышает 1 %%\n");

  h = tab->header();
  h.verified_err = worst;
  tab = TableModel::from_leaves(h, [&] { std::vector<LeafRec> v; for (size_t k = 0; k < tab->n_leaves(); ++k) v.push_back(tab->rec(k)); return v; }());
  if (!tab->save(out.c_str())) { std::fprintf(stderr, "не удалось записать %s\n", out.c_str()); return 3; }
  const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("\nЗаписано: %s, %zu листьев, %.2f МБ, контрольная сумма %016llx, время %.1f с\n", out.c_str(), tab->n_leaves(),
              tab->bytes() / 1048576.0, (unsigned long long)tab->header().checksum, sec);
  if (verbose) std::printf("узлов эталона: %zu\n", g_cache.size());
  return 0;
}
