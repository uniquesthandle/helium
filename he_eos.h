// he_eos.h - Вариант B: компактное аналитическое уравнение состояния гелия-4.
//
// Остаточная энергия Гельмгольца, Phi = a^r/(R T), в переменных (T, rho):
//   Phi = phi_CS(eta) + sum_i n_i * (T/300)^(-t_i) * (rho/100)^(d_i),   eta = b(T)*rho/4,
//   phi_CS(eta) = eta (4 - 3 eta)/(1 - eta)^2   (твердые сферы, Карнахан-Старлинг),
//   b(T) = B2(T) – второй вириальный коэффициент (полином по ln(T/300)).
// Все свойства – аналитические производные Phi (термодинамически согласованы).
//
// Перенос: mu = mu0(T) * Fmu(eta) * exp(D_mu(s, rho)), lambda = lambda0(T) * Flam(eta) * exp(D_lam(s, rho)),
// где mu0, lambda0 – полиномы ln mu0, ln lambda0 по s = ln(T/300), F – плотностные поправки Энскога,
// D – подгоняемая добавка (по умолчанию 0).
//
// Коэффициенты по умолчанию – встроенные (НЕ аттестованы). Для эксплуатации их нужно
// подогнать под аттестованный эталон (REFPROP/CoolProp/NIST) утилитой he_fit и загрузить
// из файла: AnalyticModel::load("coef.txt").
#pragma once
#include <string>

#include "he_api.h"

namespace he {

constexpr int EOS_MAXC = 40;

struct EosCoef {
  double Bc[6];  // B2(T), см3/моль: sum Bc[k] s^k, s = ln(T/300)
  int ncorr = 0;
  double cn[EOS_MAXC], ct[EOS_MAXC], cd[EOS_MAXC];  // поправочные члены n, t, d
  double mu0[4];   // ln mu0 [Па с] = sum mu0[k] s^k
  double lam0[4];  // ln lambda0 [Вт/(м К)] = sum lam0[k] s^k
  double mud[9], lamd[9];  // подгоняемые добавки к ln mu, ln lambda: sum_{j<3,k<3} d[3j+k] s^k r^(j+1), r = rho/100
  EosCoef();       // значения по умолчанию
};

class AnalyticModel : public Model {
 public:
  AnalyticModel() {}
  explicit AnalyticModel(const EosCoef& c) : c_(c) {}
  // (поправочные члены сгруппированы: sum_i n_i th^(-t_i) r^(d_i) = sum_j r^(D_j) * A_j(T) – степени по T вычисляются
  //  по различным t один раз на вызов, степени по rho – умножением)
  bool load(const char* path, std::string* err = nullptr);   // текстовый файл коэффициентов
  bool save(const char* path) const;
  const EosCoef& coef() const { return c_; }
  EosCoef& coef() { dirty_ = true; return c_; }   // после изменений вызовите finalize() (иначе – при первом расчете)
  void finalize() const;                           // группировка членов по (t, d); не потокобезопасна при первом вызове

  const char* name() const override { return "analytic (Helmholtz-type EOS)"; }
  int thermo(double T, double P, unsigned need, Thermo& s, Hint* h) const override;
  int transport(double T, double P, double rho, unsigned which, double& lam, double& mu, Hint* h) const override;
  bool transport_needs_rho() const override { return true; }

  // Логарифм коэффициента летучести ln(phi) = g^r/(R T) – нужен построителю таблицы
  int lnphi(double T, double P, double& g) const;
  // Z(T, rho) для подгонки
  double Z_of(double T, double rho) const;
  double b_of(double T) const;  // b(T) в м3/кг
  // Плотностные поправки Энскога (для подгонки коэффициентов переноса)
  void enskog(double T, double rho, double& fmu, double& flam) const;

 private:
  EosCoef c_;
  struct Prep { double T, b, b1, b2; double A0[EOS_MAXC], A1[EOS_MAXC], A2[EOS_MAXC]; };
  mutable bool dirty_ = true;
  mutable int nt_ = 0, nd_ = 0;
  mutable double tv_[EOS_MAXC], dv_[EOS_MAXC], w_[EOS_MAXC][EOS_MAXC];
  struct Hd { double Phi, Pr, Prr, PT, PTT, PrT; };
  void prep(double T, Prep& p) const;
  void helm(const Prep& p, double rho, int lvl, Hd& d) const;
  int solve_rho(const Prep& p, double T, double P, double& rho) const;
  friend struct EosFitAccess;
};

}  // namespace he
