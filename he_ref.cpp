// he_ref.cpp - эталон: (1) USE_COOLPROP – CoolProp (Helium, Ortiz-Vega) [не проверялось в этой среде];
//                      (2) по умолчанию – AnalyticModel (he_eos) со своими или загруженными коэффициентами.
#include "he_ref.h"

#include <cmath>

#include "he_eos.h"

#ifdef USE_COOLPROP
#include "CoolProp/CoolProp.h"
#include "CoolProp/AbstractState.h"
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
static std::shared_ptr<CoolProp::AbstractState> g_state;
bool ref_init(const char*, std::string* err) {
  try {
    g_state.reset(CoolProp::AbstractState::factory("HEOS", "Helium"));
    const double T0 = 300.0, P0 = 1.0; // при 1 Па остаточные части пренебрежимо малы
    g_state->update(CoolProp::PT_INPUTS, P0, T0);
    g_hshift = g_state->hmass() - 2.5 * R_SPEC * T0;
    g_s0 = g_state->smass();
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
    // Один раз обновляем состояние по (P, T) — все свойства считываются мгновенно:
    g_state->update(CoolProp::PT_INPUTS, P, T);
    o.rho = g_state->rhomass();
    o.cp = g_state->cpmass();
    o.cv = g_state->cvmass();
    o.H = g_state->hmass() - g_hshift;
    o.U = g_state->umass() - g_hshift;
    o.lam = g_state->conductivity();
    o.mu = g_state->viscosity();
    
    // Производные через first_partial_deriv:
    o.dUdT = g_state->first_partial_deriv(CoolProp::iUmass, CoolProp::iT, CoolProp::iP);
    o.drdT = g_state->first_partial_deriv(CoolProp::iDmass, CoolProp::iT, CoolProp::iP);
    o.drdP = g_state->first_partial_deriv(CoolProp::iDmass, CoolProp::iP, CoolProp::iT);
    o.drdU = o.drdT / o.dUdT;
    o.drdH = o.drdT / o.cp;

    const double s = g_state->smass();
    const double sid = g_s0 + 2.5 * R_SPEC * std::log(T / 300.0) - R_SPEC * std::log(P / 1.0);
    const double hr = o.H - 2.5 * R_SPEC * T;
    o.lnphi = (hr - T * (s - sid)) / (R_SPEC * T);
    return std::isfinite(o.rho) && std::isfinite(o.mu) && std::isfinite(o.lnphi);
  } catch (...) { return false; }
}
const char* ref_name() { return "CoolProp (Helium, Ortiz-Vega EOS)"; }
#endif

}  // namespace he
