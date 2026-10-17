// he_gas.cpp - реализация he_gas.h: гелий-4 через he::Model (таблица или аналитика), остальные газы - заглушки.
#include "he_gas.h"

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <mutex>
#include <string>

#include "he_api.h"
#include "he_eos.h"
#include "he_table.h"

#ifndef HE_GAS_MOL_SCALE
#define HE_GAS_MOL_SCALE 1.0   // ff[47] в кг/моль; 1000.0 - в кг/кмоль
#endif

namespace {

const double RM = 8.314462618;  // Дж/(моль К)

struct GasInfo { double mol; double k; };  // кг/моль; показатель адиабаты идеального газа (только для заглушек)
const GasInfo kGas[9] = {
    {28.9647e-3, 1.4},        // 0: не используется
    {28.9647e-3, 1.4},        // 1 воздух
    {28.0134e-3, 1.4},        // 2 N2
    {31.9988e-3, 1.4},        // 3 O2
    {4.002602e-3, 5.0 / 3},   // 4 He
    {2.01588e-3, 1.41},       // 5 H (принят H2; для атомарного H: 1.00794e-3)
    {44.0095e-3, 1.29},       // 6 CO2
    {39.948e-3, 5.0 / 3},     // 7 Ar
    {131.293e-3, 5.0 / 3}};   // 8 Xe

inline int gidx(int Ng) { return (Ng >= 1 && Ng <= 8) ? Ng : 1; }  // иначе воздух

inline void set_consts(int gi, double k, double* ff) {
  ff[AG_MOL] = kGas[gi].mol * HE_GAS_MOL_SCALE;
  ff[AG_K] = k;
  ff[AG_R] = RM / kGas[gi].mol;
}

// ---------------- гелий: загрузка модели ----------------
struct HeState {
  std::unique_ptr<he::Model> model;
  std::string err, info;
  double Hmin = 0, Hmax = 0, Umin = 0, Umax = 0;
};
HeState g_he;
std::once_flag g_once;
std::string g_table_path, g_coef_path;
std::atomic<int> g_strict{0};
thread_local int t_flags = 0;
thread_local he::Hint t_hint;

std::string env_or(const char* name, const char* dflt) {
  const char* e = std::getenv(name);
  return (e && *e) ? std::string(e) : std::string(dflt);
}

void scan_bounds(HeState& s) {
  // Границы по паспортному диапазону: T = T_MIN, T_MAX; P - логарифмическая сетка (H, U монотонны по T).
  const double Ts[2] = {he::T_MIN, he::T_MAX};
  double hmin = std::numeric_limits<double>::infinity(), hmax = -hmin, umin = hmin, umax = -hmin;
  bool ok = (s.model != nullptr);
  if (ok) {
    he::Hint h;
    const int NP = 14;
    for (int it = 0; it < 2 && ok; ++it)
      for (int ip = 0; ip <= NP; ++ip) {
        const double P = (ip == NP) ? he::P_MAX : he::P_MIN * std::pow(he::P_MAX / he::P_MIN, double(ip) / NP);
        he::Result r;
        if (s.model->PT(P, Ts[it], he::HP_H | he::HP_U, r, &h) != he::HE_OK || !std::isfinite(r.H) || !std::isfinite(r.U)) {
          ok = false;
          break;
        }
        if (r.H < hmin) hmin = r.H;
        if (r.H > hmax) hmax = r.H;
        if (r.U < umin) umin = r.U;
        if (r.U > umax) umax = r.U;
      }
  }
  if (!ok) {  // запасной вариант: идеальный газ (ошибка реального газа до ~2 % по H)
    hmin = 2.5 * he::R_SPEC * he::T_MIN; hmax = 2.5 * he::R_SPEC * he::T_MAX;
    umin = 1.5 * he::R_SPEC * he::T_MIN; umax = 1.5 * he::R_SPEC * he::T_MAX;
  }
  s.Hmin = hmin; s.Hmax = hmax; s.Umin = umin; s.Umax = umax;
}

void init_he() {
  const std::string tp = g_table_path.empty() ? env_or("HE_TABLE_FILE", "he_table2.bin") : g_table_path;
  const std::string cp = g_coef_path.empty() ? env_or("HE_COEF_FILE", "data/he_coef_abinitio.txt") : g_coef_path;
  std::string e1;
  std::unique_ptr<he::TableModel> t = he::TableModel::open(tp.c_str(), e1, he::OPEN_MMAP);
  if (t) {
    g_he.info = "table: " + tp;
    g_he.model = std::move(t);
  } else {
    // Намеренно НЕ откатываемся на встроенные коэффициенты: они не привязаны к эталону (до 7 % по mu, lambda).
    std::unique_ptr<he::AnalyticModel> a(new he::AnalyticModel());
    std::string e2;
    if (a->load(cp.c_str(), &e2)) {
      g_he.info = "analytic: " + cp;
      g_he.model = std::move(a);
    } else {
      g_he.err = "he table '" + tp + "': " + e1 + "; coefficients '" + cp + "': " + e2;
    }
  }
  scan_bounds(g_he);
}

inline bool ensure_he() {
  std::call_once(g_once, init_he);
  return g_he.model != nullptr;
}

// ---------------- гелий: расчет ----------------
const unsigned kBits[AG_NPROP] = {
    he::HP_RHO, he::HP_CP, he::HP_MU, he::HP_LAM, he::HP_U, he::HP_DUDT, he::HP_DRDT, he::HP_DRDP,
    he::HP_DRDU, he::HP_PR, he::HP_RHO | he::HP_DRDT, he::HP_H, he::HP_CV, he::HP_DRDH, he::HP_T};

inline bool want(const int* nff, int i) { return nff == nullptr || nff[i] == 1; }

double pick(const he::Result& r, int i) {
  switch (i) {
    case AG_RHO: return r.rho;
    case AG_CP: return r.cp;
    case AG_MU: return r.mu;
    case AG_LAM: return r.lambda;
    case AG_U: return r.U;
    case AG_DUDT: return r.dUdT;
    case AG_DRDT: return r.drhodT;
    case AG_DRDP: return r.drhodP;
    case AG_DRDU: return r.drhodU;
    case AG_PR: return r.Pr;
    case AG_BETA: return -r.drhodT / r.rho;
    case AG_H: return r.H;
    case AG_CV: return r.cv;
    case AG_DRDH: return r.drhodH;
    default: return r.T;
  }
}

// mode: 0 - (P,T), 1 - (P,H), 2 - (P,U)
void he_calc(int mode, double P, double X, double* ff, const int* nff, int* err) {
  *err = 1;
  set_consts(4, 5.0 / 3, ff);
  t_flags = 0;
  if (!ensure_he()) return;
  unsigned mask = he::HP_CP | he::HP_CV;  // для k = cp/cv
  for (int i = 0; i < AG_NPROP; ++i)
    if (want(nff, i)) mask |= kBits[i];
  if (mode != 0) mask |= he::HP_T;
  he::Result r;
  const he::Model& m = *g_he.model;
  const int e = mode == 0 ? m.PT(P, X, mask, r, &t_hint) : mode == 1 ? m.PH(P, X, mask, r, &t_hint) : m.PU(P, X, mask, r, &t_hint);
  if (e != he::HE_OK) return;
  t_flags = int(r.flags);
  if (g_strict.load() && r.flags) return;
  double tmp[AG_NPROP];
  for (int i = 0; i < AG_NPROP; ++i) {
    if (!want(nff, i)) continue;
    if ((r.computed & kBits[i]) != kBits[i]) return;  // модель не вернула запрошенное
    tmp[i] = pick(r, i);
    if (!std::isfinite(tmp[i])) return;
  }
  for (int i = 0; i < AG_NPROP; ++i)
    if (want(nff, i)) ff[i] = tmp[i];
  if (std::isfinite(r.cp) && std::isfinite(r.cv) && r.cv > 0) ff[AG_K] = r.cp / r.cv;
  *err = 0;
}

// Общая диспетчеризация: 4 - гелий, остальное - заглушка (воздух по умолчанию).
void dispatch(int Ng, int mode, double P, double X, double* ff, int* nff, int* err) {
  int dummy;
  int* e = err ? err : &dummy;
  *e = 1;
  if (!ff) return;
  if (Ng == 4) { he_calc(mode, P, X, ff, nff, e); return; }
  // --- ЗАГЛУШКА: свойства других газов не реализованы ---
  const int gi = gidx(Ng);
  set_consts(gi, kGas[gi].k, ff);
  *e = 1;
}

}  // namespace

// ---------------- публичный интерфейс ----------------
HE_GAS_API void ASPIG(int Ng, double P, double T, double* ff, int* nff, int* err) { dispatch(Ng, 0, P, T, ff, nff, err); }
HE_GAS_API void ASPIGH(int Ng, double P, double H, double* ff, int* nff, int* err) { dispatch(Ng, 1, P, H, ff, nff, err); }
HE_GAS_API void ASPIGU(int Ng, double P, double U, double* ff, int* nff, int* err) { dispatch(Ng, 2, P, U, ff, nff, err); }

HE_GAS_API void AGHmax(int Ng, double* Hmin, double* Hmax) {
  double a = 0, b = 0;
  if (Ng == 4) { ensure_he(); a = g_he.Hmin; b = g_he.Hmax; }  // остальные газы - заглушка (0, 0)
  if (Hmin) *Hmin = a;
  if (Hmax) *Hmax = b;
}

HE_GAS_API void AGUmax(int Ng, double* Umin, double* Umax) {
  double a = 0, b = 0;
  if (Ng == 4) { ensure_he(); a = g_he.Umin; b = g_he.Umax; }
  if (Umin) *Umin = a;
  if (Umax) *Umax = b;
}

HE_GAS_API int ASPIGInit(const char* table_path, const char* coef_path) {
  bool done = true;
  std::call_once(g_once, [&] {
    done = false;
    g_table_path = table_path ? table_path : "";
    g_coef_path = coef_path ? coef_path : "";
    init_he();
  });
  return done ? -1 : 0;
}

HE_GAS_API int AGLastFlags(void) { return t_flags; }
HE_GAS_API void AGSetStrict(int strict) { g_strict.store(strict ? 1 : 0); }
HE_GAS_API const char* AGLastError(void) { ensure_he(); return g_he.err.c_str(); }
HE_GAS_API const char* AGModelInfo(void) { ensure_he(); return g_he.info.c_str(); }
const he::Model* AGGetHeModel() { ensure_he(); return g_he.model.get(); }
