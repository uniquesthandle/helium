// he_ref.h - эталонная модель для построителя таблицы (и утилиты проверки)
#pragma once
#include <string>

namespace he {

struct RefState {
  double rho, cp, cv, U, H, lam, mu;       // СИ, нуль U и H – идеальный газ при 0 К
  double dUdT, drdT, drdU, drdH, drdP;     // (P=const), drdP при T=const
  double lnphi;                            // g^r/(R T) – логарифм коэффициента летучести
};

// Источник по умолчанию – аналитическая модель (he_eos); можно подгрузить файл ее коэффициентов.
bool ref_init(const char* coef_path, std::string* err);
bool ref_eval(double T, double P, RefState& out);
const char* ref_name();

}  // namespace he
