// he_api.h  -  общий программный интерфейс расчета свойств гелия-4 (версия 2)
//
// Диапазон: T = 50...3000 К, P = 10 Па ... 100 МПа (ниже 50 К – не поддерживается).
// Две реализации одного интерфейса:
//   * TableModel    (he_table.h) – таблица остаточного потенциала Гиббса, интерполяция C2;
//   * AnalyticModel (he_eos.h)   – компактное аналитическое уравнение состояния.
//
// Три вида вызова (до 13 величин + температура):
//   PT(P, T, mask, out)  – по давлению и температуре
//   PH(P, H, mask, out)  – по давлению и энтальпии     (T вычисляется)
//   PU(P, U, mask, out)  – по давлению и внутренней энергии (T вычисляется)
// Считаются только величины, указанные в маске (плюс то, без чего их не получить).
//
// Единицы СИ. Нуль U и H: идеальный газ при 0 К (U = 1.5 R T, H = 2.5 R T при P -> 0).
#pragma once
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

namespace he {

constexpr double T_MIN = 50.0, T_MAX = 3000.0;   // паспортный диапазон
constexpr double P_MIN = 10.0, P_MAX = 1.0e8;
// Жесткие границы: внутри них, но вне паспортного диапазона, считается экстраполяция
// с флагом; вне них возвращается HE_BAD_INPUT.
constexpr double T_HARD_LO = 20.0, T_HARD_HI = 6000.0;
constexpr double P_HARD_LO = 1.0e-3, P_HARD_HI = 1.0e9;
constexpr double R_SPEC = 8.314462618 / 4.002602e-3;  // Дж/(кг К)

// Маска запрашиваемых величин
enum Prop : unsigned {
  HP_T = 1u << 0,     // температура (для PT – входная, для PH/PU – выходная)
  HP_RHO = 1u << 1,   // плотность, кг/м3
  HP_CP = 1u << 2,    // изобарная теплоемкость, Дж/(кг К)
  HP_CV = 1u << 3,    // изохорная теплоемкость, Дж/(кг К)
  HP_U = 1u << 4,     // внутренняя энергия, Дж/кг
  HP_H = 1u << 5,     // энтальпия, Дж/кг
  HP_LAM = 1u << 6,   // теплопроводность, Вт/(м К)
  HP_MU = 1u << 7,    // динамическая вязкость, Па с
  HP_PR = 1u << 8,    // число Прандтля
  HP_DUDT = 1u << 9,  // dU/dT при P = const
  HP_DRDT = 1u << 10, // drho/dT при P = const
  HP_DRDU = 1u << 11, // drho/dU при P = const
  HP_DRDH = 1u << 12, // drho/dH при P = const
  HP_DRDP = 1u << 13, // drho/dP при T = const
  HP_ALL = (1u << 14) - 1
};

// Флаги экстраполяции (Result::flags)
enum Flag : unsigned { FL_LOW_T = 1, FL_HIGH_T = 2, FL_LOW_P = 4, FL_HIGH_P = 8 };

// Коды возврата
enum Err : int {
  HE_OK = 0,
  HE_BAD_INPUT = -1,    // NaN, неположительные значения или вне жестких границ
  HE_NO_SOLUTION = -2,  // для (P,H)/(P,U) нет температуры в допустимом интервале
  HE_MODEL_FAIL = -3    // не сошелся расчет плотности / ошибка модели
};

struct Result {
  double P, T, rho, cp, cv, U, H, lambda, mu, Pr;
  double dUdT, drhodT, drhodU, drhodH, drhodP;
  unsigned computed;  // маска фактически вычисленных величин (невычисленные = NaN)
  unsigned flags;     // FL_*: точка вне паспортного диапазона (значения – экстраполяция)
  int iters;          // число итераций по T (PH/PU), 0 для PT
};

// Подсказка между вызовами (потоково-локальная; один объект на поток/канал).
struct Hint {
  double T = 0.0;  // последняя найденная температура (начальное приближение для PH/PU)
  int cell = -1;   // последняя ячейка таблицы
};

// Термодинамическое ядро: то, что реализуют модели
enum Need : unsigned { N_RHO = 1, N_H = 2, N_CP = 4, N_DV = 8 };
struct Thermo {
  double rho = NAN;   // кг/м3
  double h = NAN;     // Дж/кг
  double cp = NAN;    // Дж/(кг К)
  double dvdT = NAN;  // (dv/dT)_P, м3/(кг К), v = 1/rho
  double dvdP = NAN;  // (dv/dP)_T, м3/(кг Па)
};

unsigned thermo_need(unsigned mask);  // какие части ядра нужны для маски
void assemble(double P, double T, unsigned mask, const Thermo& s, double lam, double mu, Result& r);

class Model {
 public:
  virtual ~Model() {}
  virtual const char* name() const = 0;

  int PT(double P, double T, unsigned mask, Result& r, Hint* h = nullptr) const;
  int PH(double P, double H, unsigned mask, Result& r, Hint* h = nullptr) const;
  int PU(double P, double U, unsigned mask, Result& r, Hint* h = nullptr) const;

  // Пакетные вызовы (для одномерных каналов): подсказка передается между точками.
  // status[i] – код возврата для точки i (может быть nullptr). Возвращает число ошибок.
  int PT_batch(size_t n, const double* P, const double* T, unsigned mask, Result* out, int* status = nullptr) const;
  int PH_batch(size_t n, const double* P, const double* H, unsigned mask, Result* out, int* status = nullptr) const;
  int PU_batch(size_t n, const double* P, const double* U, unsigned mask, Result* out, int* status = nullptr) const;

  // --- низкоуровневые функции моделей ---
  virtual int thermo(double T, double P, unsigned need, Thermo& s, Hint* h) const = 0;
  // which: HP_LAM/HP_MU/HP_PR; rho нужна только если transport_needs_rho()
  virtual int transport(double T, double P, double rho, unsigned which, double& lam, double& mu, Hint* h) const = 0;
  virtual bool transport_needs_rho() const { return false; }

 private:
  int invert(double P, double X, bool isU, unsigned mask, Result& r, Hint* h) const;
};

}  // namespace he
