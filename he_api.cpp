// he_api.cpp - общая часть: маски, вызовы PT/PH/PU, обращение по T, пакетные вызовы
#include "he_api.h"

namespace he {

static const double NaN = std::numeric_limits<double>::quiet_NaN();

unsigned thermo_need(unsigned mask) {
  unsigned n = 0;
  if (mask & HP_RHO) n |= N_RHO;
  if (mask & HP_H) n |= N_H;
  if (mask & HP_U) n |= N_H | N_RHO;
  if (mask & (HP_CP | HP_PR)) n |= N_CP;
  if (mask & (HP_CV | HP_DUDT)) n |= N_CP | N_DV;
  if (mask & (HP_DRDT | HP_DRDP)) n |= N_RHO | N_DV;
  if (mask & (HP_DRDU | HP_DRDH)) n |= N_RHO | N_DV | N_CP;
  return n;
}

static void init_result(Result& r) {
  r.P = r.T = r.rho = r.cp = r.cv = r.U = r.H = r.lambda = r.mu = r.Pr = NaN;
  r.dUdT = r.drhodT = r.drhodU = r.drhodH = r.drhodP = NaN;
  r.computed = 0;
  r.flags = 0;
  r.iters = 0;
}

// Все производные выражаются через (rho, h, cp, dv/dT, dv/dP) – поэтому согласованы по построению.
void assemble(double P, double T, unsigned mask, const Thermo& s, double lam, double mu, Result& r) {
  r.P = P;
  r.T = T;
  unsigned c = HP_T;
  if (mask & HP_RHO) { r.rho = s.rho; c |= HP_RHO; }
  if (mask & HP_H) { r.H = s.h; c |= HP_H; }
  if (mask & HP_CP) { r.cp = s.cp; c |= HP_CP; }
  if (mask & HP_U) { r.U = s.h - P / s.rho; c |= HP_U; }
  if (mask & HP_CV) { r.cv = s.cp + T * s.dvdT * s.dvdT / s.dvdP; c |= HP_CV; }
  const double rho2 = s.rho * s.rho;
  if (mask & HP_DRDT) { r.drhodT = -rho2 * s.dvdT; c |= HP_DRDT; }
  if (mask & HP_DRDP) { r.drhodP = -rho2 * s.dvdP; c |= HP_DRDP; }
  if (mask & HP_DUDT) { r.dUdT = s.cp - P * s.dvdT; c |= HP_DUDT; }
  if (mask & HP_DRDU) { r.drhodU = (-rho2 * s.dvdT) / (s.cp - P * s.dvdT); c |= HP_DRDU; }
  if (mask & HP_DRDH) { r.drhodH = (-rho2 * s.dvdT) / s.cp; c |= HP_DRDH; }
  if (mask & HP_LAM) { r.lambda = lam; c |= HP_LAM; }
  if (mask & HP_MU) { r.mu = mu; c |= HP_MU; }
  if (mask & HP_PR) { r.Pr = s.cp * mu / lam; c |= HP_PR; }
  r.computed = c;
}

int Model::PT(double P, double T, unsigned mask, Result& r, Hint* h) const {
  init_result(r);
  if (!(P >= P_HARD_LO && P <= P_HARD_HI && T >= T_HARD_LO && T <= T_HARD_HI)) return HE_BAD_INPUT;
  if (T < T_MIN) r.flags |= FL_LOW_T;
  if (T > T_MAX) r.flags |= FL_HIGH_T;
  if (P < P_MIN) r.flags |= FL_LOW_P;
  if (P > P_MAX) r.flags |= FL_HIGH_P;

  const bool tr = (mask & (HP_LAM | HP_MU | HP_PR)) != 0;
  unsigned need = thermo_need(mask);
  if (tr && transport_needs_rho()) need |= N_RHO;
  Thermo s;
  if (need) {
    const int e = thermo(T, P, need, s, h);
    if (e) return e;
  }
  double lam = NaN, mu = NaN;
  if (tr) {
    const int e = transport(T, P, s.rho, mask & (HP_LAM | HP_MU | HP_PR), lam, mu, h);
    if (e) return e;
  }
  assemble(P, T, mask, s, lam, mu, r);
  return HE_OK;
}

// Решение h(T,P) = H  (или u(T,P) = U) относительно T методом Ньютона с контролем скобки.
int Model::invert(double P, double X, bool isU, unsigned mask, Result& r, Hint* hint) const {
  init_result(r);
  if (!(P >= P_HARD_LO && P <= P_HARD_HI && std::isfinite(X) && X > 0.0)) return HE_BAD_INPUT;

  double T = (hint && hint->T >= T_HARD_LO && hint->T <= T_HARD_HI) ? hint->T
                                                                    : X / ((isU ? 1.5 : 2.5) * R_SPEC);
  double lo = T_HARD_LO, hi = T_HARD_HI;
  T = std::fmin(std::fmax(T, lo * 1.001), hi * 0.999);
  const unsigned need = N_H | N_CP | (isU ? (N_RHO | N_DV) : 0u);

  bool ok = false;
  int it = 0;
  for (; it < 60; ++it) {
    Thermo s;
    const int e = thermo(T, P, need, s, hint);
    if (e) return e;
    const double f = (isU ? s.h - P / s.rho : s.h) - X;
    const double d = isU ? s.cp - P * s.dvdT : s.cp;
    if (f == 0.0) { ok = true; break; }
    if (f > 0.0) hi = T; else lo = T;
    double Tn = T - f / d;
    if (std::fabs(Tn - T) <= 1e-13 * T) { T = Tn; ok = true; break; }
    if (!(Tn > lo && Tn < hi)) Tn = 0.5 * (lo + hi);  // защита: деление пополам
    T = Tn;
  }
  if (!ok || T < T_HARD_LO * (1 + 1e-6) || T > T_HARD_HI * (1 - 1e-6)) return HE_NO_SOLUTION;

  const int e = PT(P, T, mask, r, hint);
  if (e) return e;
  r.iters = it + 1;
  if (isU) { r.U = X; r.computed |= HP_U; } else { r.H = X; r.computed |= HP_H; }
  if (hint) hint->T = T;
  return HE_OK;
}

int Model::PH(double P, double H, unsigned mask, Result& r, Hint* h) const { return invert(P, H, false, mask, r, h); }
int Model::PU(double P, double U, unsigned mask, Result& r, Hint* h) const { return invert(P, U, true, mask, r, h); }

#define HE_BATCH(NAME, CALL)                                                                          \
  int Model::NAME(size_t n, const double* P, const double* X, unsigned mask, Result* out, int* status) const { \
    Hint h;                                                                                           \
    int bad = 0;                                                                                      \
    for (size_t i = 0; i < n; ++i) {                                                                  \
      const int e = CALL;                                                                             \
      if (status) status[i] = e;                                                                      \
      if (e) ++bad;                                                                                   \
    }                                                                                                 \
    return bad;                                                                                       \
  }
HE_BATCH(PT_batch, PT(P[i], X[i], mask, out[i], &h))
HE_BATCH(PH_batch, PH(P[i], X[i], mask, out[i], &h))
HE_BATCH(PU_batch, PU(P[i], X[i], mask, out[i], &h))
#undef HE_BATCH

}  // namespace he
