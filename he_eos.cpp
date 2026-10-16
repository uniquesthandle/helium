// he_eos.cpp - Вариант B: аналитическое уравнение состояния и корреляции переноса
#include "he_eos.h"

#include <cstdio>
#include <cstring>

namespace he {

static const double MW = 4.002602e-3;  // кг/моль

EosCoef::EosCoef() {
  static const double B0[6] = {11.93825937, -0.9122812673, -0.8253417466, 0.2744727437, -0.002339236121, -0.01473705211};
  for (int k = 0; k < 6; ++k) Bc[k] = B0[k];
  std::memset(cn, 0, sizeof cn);
  std::memset(ct, 0, sizeof ct);
  std::memset(cd, 0, sizeof cd);
  // mu0 = 3.674e-7 T^0.7 (Petersen), lambda0 = (15/4) R mu0
  mu0[0] = std::log(3.674e-7) + 0.7 * std::log(300.0);
  mu0[1] = 0.7;
  mu0[2] = mu0[3] = 0.0;
  for (int k = 0; k < 9; ++k) mud[k] = lamd[k] = 0.0;
  for (int k = 0; k < 4; ++k) lam0[k] = mu0[k];
  lam0[0] += std::log(3.75 * R_SPEC);
}

void AnalyticModel::finalize() const {
  nt_ = nd_ = 0;
  for (int i = 0; i < c_.ncorr; ++i) {
    int k = 0, j = 0;
    while (k < nt_ && std::fabs(tv_[k] - c_.ct[i]) > 1e-12) ++k;
    if (k == nt_) tv_[nt_++] = c_.ct[i];
    while (j < nd_ && std::fabs(dv_[j] - c_.cd[i]) > 1e-12) ++j;
    if (j == nd_) { dv_[nd_++] = c_.cd[i]; for (int q = 0; q < EOS_MAXC; ++q) w_[j][q] = 0.0; }
    w_[j][k] += c_.cn[i];
  }
  dirty_ = false;
}

void AnalyticModel::prep(double T, Prep& p) const {
  if (dirty_) finalize();
  const double s = std::log(T / 300.0);
  double v = 0, dv = 0, ddv = 0;
  for (int k = 5; k >= 0; --k) {
    ddv = ddv * s + 2 * dv;
    dv = dv * s + v;
    v = v * s + c_.Bc[k];
  }
  const double conv = 1e-6 / MW;  // см3/моль -> м3/кг
  p.T = T;
  p.b = v * conv;
  p.b1 = dv / T * conv;
  p.b2 = (ddv - dv) / (T * T) * conv;
  if (nd_) {
    double th[EOS_MAXC];
    for (int k = 0; k < nt_; ++k) th[k] = std::exp(-tv_[k] * s);
    const double iT = 1.0 / T;
    for (int j = 0; j < nd_; ++j) {
      double a0 = 0, a1 = 0, a2 = 0;
      for (int k = 0; k < nt_; ++k) {
        const double q = w_[j][k] * th[k];
        a0 += q;
        a1 += q * (-tv_[k]) * iT;
        a2 += q * tv_[k] * (tv_[k] + 1) * iT * iT;
      }
      p.A0[j] = a0; p.A1[j] = a1; p.A2[j] = a2;
    }
  }
}

double AnalyticModel::b_of(double T) const {
  Prep p;
  prep(T, p);
  return p.b;
}

// lvl 1: Phi, Phi_rho, Phi_rhorho;  lvl 2: + Phi_T, Phi_TT, Phi_rhoT
void AnalyticModel::helm(const Prep& p, double rho, int lvl, Hd& d) const {
  const double eta = p.b * rho * 0.25, w = 1 - eta, w2 = w * w, w3 = w2 * w, w4 = w3 * w;
  const double phi = eta * (4 - 3 * eta) / w2, f1 = (4 - 2 * eta) / w3, f2 = (10 - 4 * eta) / w4;
  d.Phi = phi;
  d.Pr = f1 * p.b * 0.25;
  d.Prr = f2 * p.b * p.b * (1.0 / 16);
  d.PT = d.PTT = d.PrT = 0.0;
  const double rb1 = rho * p.b1 * 0.25;
  if (lvl >= 2) {
    d.PT = f1 * rb1;
    d.PTT = f2 * rb1 * rb1 + f1 * rho * p.b2 * 0.25;
    d.PrT = f2 * (p.b * 0.25) * rb1 + f1 * p.b1 * 0.25;
  }
  if (nd_) {
    const double r = rho / 100.0, ir = 1.0 / r;
    for (int j = 0; j < nd_; ++j) {
      const double dd = dv_[j];
      const int id = (int)dd;
      double rd;
      if (dd == (double)id && id >= 1 && id <= 4) { rd = r; for (int q = 1; q < id; ++q) rd *= r; }
      else rd = std::pow(r, dd);
      const double rd1 = rd * ir, rd2 = rd1 * ir;
      d.Phi += p.A0[j] * rd;
      d.Pr += p.A0[j] * dd * rd1 * 0.01;
      d.Prr += p.A0[j] * dd * (dd - 1) * rd2 * 1e-4;
      if (lvl >= 2) {
        d.PT += p.A1[j] * rd;
        d.PTT += p.A2[j] * rd;
        d.PrT += p.A1[j] * dd * rd1 * 0.01;
      }
    }
  }
}

int AnalyticModel::solve_rho(const Prep& p, double T, double P, double& rho) const {
  const double RT = R_SPEC * T;
  double lo = 0.0, hi = p.b > 0 ? 0.999 * 4.0 / p.b : 1.0e4;  // b <= 0 (T < ~24 K): верхняя граница по плотности – условная
  rho = P / RT;
  if (rho >= hi) rho = 0.5 * hi;
  Hd d;
  for (int it = 0; it < 200; ++it) {
    helm(p, rho, 1, d);
    const double f = rho * RT * (1 + rho * d.Pr) - P;
    const double df = RT * (1 + 2 * rho * d.Pr + rho * rho * d.Prr);
    if (f > 0) hi = rho; else lo = rho;
    double rn = rho - f / df;
    if (std::fabs(rn - rho) <= 1e-14 * rn && rn > 0) { rho = rn; return 0; }
    if (!(rn > lo && rn < hi)) rn = 0.5 * (lo + hi);
    rho = rn;
  }
  return HE_MODEL_FAIL;
}

int AnalyticModel::thermo(double T, double P, unsigned need, Thermo& s, Hint*) const {
  Prep p;
  prep(T, p);
  double rho;
  if (solve_rho(p, T, P, rho)) return HE_MODEL_FAIL;
  const bool full = (need & (N_H | N_CP | N_DV)) != 0;
  Hd d;
  helm(p, rho, full ? 2 : 1, d);
  const double R = R_SPEC, Z = 1 + rho * d.Pr;
  s.rho = rho;
  if (need & N_H) s.h = 1.5 * R * T - R * T * T * d.PT + R * T * Z;
  const double Prho = R * T * (1 + 2 * rho * d.Pr + rho * rho * d.Prr);
  const double PTv = rho * R * (1 + rho * d.Pr + T * rho * d.PrT);
  if (need & N_DV) {
    s.dvdT = PTv / (rho * rho * Prho);
    s.dvdP = -1.0 / (rho * rho * Prho);
  }
  if (need & N_CP) {
    const double cv = 1.5 * R - R * (2 * T * d.PT + T * T * d.PTT);
    s.cp = cv + T * PTv * PTv / (rho * rho * Prho);
  }
  return HE_OK;
}

int AnalyticModel::transport(double T, double, double rho, unsigned which, double& lam, double& mu, Hint*) const {
  Prep p;
  prep(T, p);
  const double eta = p.b * rho * 0.25, w = 1 - eta;
  const double chi = (1 - 0.5 * eta) / (w * w * w), br = 4 * eta;
  const double s = std::log(T / 300.0), r = rho / 100.0;
  auto dcor = [&](const double* d) {
    return r * (d[0] + s * (d[1] + s * d[2])) + r * r * (d[3] + s * (d[4] + s * d[5])) + r * r * r * (d[6] + s * (d[7] + s * d[8]));
  };
  if (which & (HP_MU | HP_PR)) {
    const double m0 = std::exp(c_.mu0[0] + s * (c_.mu0[1] + s * (c_.mu0[2] + s * c_.mu0[3])) + dcor(c_.mud));
    mu = m0 * (1 / chi + 0.8 * br + 0.7614 * br * br * chi);
  }
  if (which & (HP_LAM | HP_PR)) {
    const double l0 = std::exp(c_.lam0[0] + s * (c_.lam0[1] + s * (c_.lam0[2] + s * c_.lam0[3])) + dcor(c_.lamd));
    lam = l0 * (1 / chi + 1.2 * br + 0.7574 * br * br * chi);
  }
  return HE_OK;
}

void AnalyticModel::enskog(double T, double rho, double& fmu, double& flam) const {
  const double eta = b_of(T) * rho * 0.25, w = 1 - eta;
  const double chi = (1 - 0.5 * eta) / (w * w * w), br = 4 * eta;
  fmu = 1 / chi + 0.8 * br + 0.7614 * br * br * chi;
  flam = 1 / chi + 1.2 * br + 0.7574 * br * br * chi;
}

int AnalyticModel::lnphi(double T, double P, double& g) const {
  Prep p;
  prep(T, p);
  double rho;
  if (solve_rho(p, T, P, rho)) return HE_MODEL_FAIL;
  Hd d;
  helm(p, rho, 1, d);
  const double Z = 1 + rho * d.Pr;
  g = d.Phi + Z - 1 - std::log(Z);
  return HE_OK;
}

double AnalyticModel::Z_of(double T, double rho) const {
  Prep p;
  prep(T, p);
  Hd d;
  helm(p, rho, 1, d);
  return 1 + rho * d.Pr;
}

bool AnalyticModel::save(const char* path) const {
  FILE* f = std::fopen(path, "w");
  if (!f) return false;
  std::fprintf(f, "# коэффициенты аналитической модели гелия-4 (he_eos)\nB");
  for (int k = 0; k < 6; ++k) std::fprintf(f, " %.12g", c_.Bc[k]);
  std::fprintf(f, "\nMU0");
  for (int k = 0; k < 4; ++k) std::fprintf(f, " %.12g", c_.mu0[k]);
  std::fprintf(f, "\nLAM0");
  for (int k = 0; k < 4; ++k) std::fprintf(f, " %.12g", c_.lam0[k]);
  std::fprintf(f, "\nMUD");
  for (int k = 0; k < 9; ++k) std::fprintf(f, " %.12g", c_.mud[k]);
  std::fprintf(f, "\nLAMD");
  for (int k = 0; k < 9; ++k) std::fprintf(f, " %.12g", c_.lamd[k]);
  std::fprintf(f, "\n");
  for (int i = 0; i < c_.ncorr; ++i) std::fprintf(f, "CORR %.12g %.12g %.12g\n", c_.cn[i], c_.ct[i], c_.cd[i]);
  return std::fclose(f) == 0;
}

bool AnalyticModel::load(const char* path, std::string* err) {
  FILE* f = std::fopen(path, "r");
  if (!f) { if (err) *err = "не удалось открыть файл коэффициентов"; return false; }
  EosCoef c;
  c.ncorr = 0;
  char line[512];
  bool ok = true;
  while (ok && std::fgets(line, sizeof line, f)) {
    char key[32];
    if (line[0] == '#' || std::sscanf(line, "%31s", key) != 1) continue;
    const char* rest = line + std::strlen(key);
    if (!std::strcmp(key, "B")) {
      ok = std::sscanf(rest, "%lf %lf %lf %lf %lf %lf", &c.Bc[0], &c.Bc[1], &c.Bc[2], &c.Bc[3], &c.Bc[4], &c.Bc[5]) == 6;
    } else if (!std::strcmp(key, "MU0")) {
      ok = std::sscanf(rest, "%lf %lf %lf %lf", &c.mu0[0], &c.mu0[1], &c.mu0[2], &c.mu0[3]) == 4;
    } else if (!std::strcmp(key, "LAM0")) {
      ok = std::sscanf(rest, "%lf %lf %lf %lf", &c.lam0[0], &c.lam0[1], &c.lam0[2], &c.lam0[3]) == 4;
    } else if (!std::strcmp(key, "MUD") || !std::strcmp(key, "LAMD")) {
      double* d = !std::strcmp(key, "MUD") ? c.mud : c.lamd;
      ok = std::sscanf(rest, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8]) == 9;
    } else if (!std::strcmp(key, "CORR")) {
      ok = c.ncorr < EOS_MAXC && std::sscanf(rest, "%lf %lf %lf", &c.cn[c.ncorr], &c.ct[c.ncorr], &c.cd[c.ncorr]) == 3;
      if (ok) ++c.ncorr;
    } else {
      ok = false;
    }
  }
  std::fclose(f);
  if (!ok) { if (err) *err = "ошибка формата файла коэффициентов"; return false; }
  c_ = c;
  return true;
}

}  // namespace he
