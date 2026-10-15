// he_ref.cpp - эталон: (1) USE_COOLPROP – CoolProp (Helium, Ortiz-Vega) [не проверялось в этой среде];
//                      (2) по умолчанию – AnalyticModel (he_eos) со своими или загруженными коэффициентами.
#include "he_ref.h"

#include <cmath>

#include "he_eos.h"

#ifdef USE_COOLPROP
#include "CoolProp/CoolProp.h"
#endif

namespace he {

#ifndef USE_COOLPROP
static AnalyticModel g_model;
static std::string g_name = "analytic (встроенные коэффициенты)";

bool ref_init(const char* coef_path, std::string* err) {
  if (coef_path && *coef_path) {
    if (!g_model.load(coef_path, err)) return false;
    g_name = std::string("analytic (коэффициенты из ") + coef_path + ")";
  }
  return true;
}

bool ref_eval(double T, double P, RefState& o) {
  Result r;
  if (g_model.PT(P, T, HP_ALL, r) != HE_OK) return false;
  double g;
  if (g_model.lnphi(T, P, g) != HE_OK) return false;
  o = {r.rho, r.cp, r.cv, r.U, r.H, r.lambda, r.mu, r.dUdT, r.drhodT, r.drhodU, r.drhodH, r.drhodP, g};
  return true;
}
const char* ref_name() { return g_name.c_str(); }

#else
using CoolProp::PropsSI;
static double g_hshift = 0, g_s0 = 0;
static bool g_ok = false;
bool ref_init(const char*, std::string* err) {
  try {
    const double T0 = 300.0, P0 = 1.0;  // при 1 Па остаточные части пренебрежимо малы
    g_hshift = PropsSI("Hmass", "T", T0, "P", P0, "Helium") - 2.5 * R_SPEC * T0;
    g_s0 = PropsSI("Smass", "T", T0, "P", P0, "Helium");  // s_id(T0,P0)
    g_ok = true;
  } catch (...) {
    if (err) *err = "CoolProp: не удалось инициализировать Helium";
    return false;
  }
  return true;
}
bool ref_eval(double T, double P, RefState& o) {
  try {
    if (!g_ok && !ref_init(nullptr, nullptr)) return false;
    const char* F = "Helium";
    o.rho = PropsSI("Dmass", "T", T, "P", P, F);
    o.cp = PropsSI("Cpmass", "T", T, "P", P, F);
    o.cv = PropsSI("Cvmass", "T", T, "P", P, F);
    o.H = PropsSI("Hmass", "T", T, "P", P, F) - g_hshift;
    o.U = PropsSI("Umass", "T", T, "P", P, F) - g_hshift;
    o.lam = PropsSI("conductivity", "T", T, "P", P, F);
    o.mu = PropsSI("viscosity", "T", T, "P", P, F);
    o.dUdT = PropsSI("d(Umass)/d(T)|P", "T", T, "P", P, F);
    o.drdT = PropsSI("d(Dmass)/d(T)|P", "T", T, "P", P, F);
    o.drdU = PropsSI("d(Dmass)/d(Umass)|P", "T", T, "P", P, F);
    o.drdH = PropsSI("d(Dmass)/d(Hmass)|P", "T", T, "P", P, F);
    o.drdP = PropsSI("d(Dmass)/d(P)|T", "T", T, "P", P, F);
    const double s = PropsSI("Smass", "T", T, "P", P, F);
    const double sid = g_s0 + 2.5 * R_SPEC * std::log(T / 300.0) - R_SPEC * std::log(P / 1.0);
    const double hr = o.H - 2.5 * R_SPEC * T;
    o.lnphi = (hr - T * (s - sid)) / (R_SPEC * T);
    return std::isfinite(o.rho) && std::isfinite(o.mu) && std::isfinite(o.lnphi);
  } catch (...) { return false; }
}
const char* ref_name() { return "CoolProp (Helium, Ortiz-Vega EOS)"; }
#endif

}  // namespace he
