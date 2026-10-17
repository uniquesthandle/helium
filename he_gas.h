// he_gas.h - интерфейс свойств газов для общего проекта.
//   ASPIG (Ng,P,T,...), ASPIGU(Ng,P,U,...), ASPIGH(Ng,P,H,...), AGHmax, AGUmax.
// Реализован только гелий-4 (Ng = 4); остальные газы - заглушки (err = 1, заполняются только ff[47..49]).
//
// Единицы: СИ (P - Па, T - К, U и H - Дж/кг). Нуль U и H как в he_api.h: идеальный газ при 0 К.
//
// ff[0..14] считаются только при nff[i] == 1 (nff == NULL - считать все 15); ff[47..49] - всегда.
//   ff[0] rho, ff[1] cp, ff[2] mu, ff[3] lambda, ff[4] U, ff[5] dU/dT (P=const, Дж/(кг К)),
//   ff[6] drho/dT (P), ff[7] drho/dP (T), ff[8] drho/dU (P), ff[9] Pr, ff[10] beta = -(1/rho)(drho/dT)_P,
//   ff[11] H, ff[12] cv, ff[13] drho/dH (P), ff[14] T,
//   ff[47] молярная масса (кг/моль; для кг/кмоль задайте -DHE_GAS_MOL_SCALE=1000.0),
//   ff[48] k = cp/cv (реальный газ), ff[49] R = Rm/Mol, Дж/(кг К).
// err = 0 - успех, 1 - ошибка (при ошибке ff[0..14] не изменяются).
// Точка вне паспортного диапазона (T 50-3000 К, P 10 Па - 100 МПа) в пределах жестких границ модели
// (T 20-6000 К, P 1 мПа - 1 ГПа) считается экстраполяцией с err = 0; признак - AGLastFlags()
// (биты he::FL_*), либо AGSetStrict(1) - тогда такие точки дают err = 1.
//
// Компоновка: по умолчанию C++ (имена с сигнатурой). Для extern "C": -DHE_GAS_C_LINKAGE.
#pragma once

#if defined(HE_GAS_C_LINKAGE) && defined(__cplusplus)
#define HE_GAS_API extern "C"
#else
#define HE_GAS_API
#endif

enum {
  AG_RHO = 0, AG_CP = 1, AG_MU = 2, AG_LAM = 3, AG_U = 4, AG_DUDT = 5, AG_DRDT = 6, AG_DRDP = 7,
  AG_DRDU = 8, AG_PR = 9, AG_BETA = 10, AG_H = 11, AG_CV = 12, AG_DRDH = 13, AG_T = 14,
  AG_NPROP = 15, AG_MOL = 47, AG_K = 48, AG_R = 49, AG_NFF = 50  // минимальный размер ff
};

HE_GAS_API void ASPIG (int Ng, double P, double T, double* ff, int* nff, int* err);
HE_GAS_API void ASPIGU(int Ng, double P, double U, double* ff, int* nff, int* err);
HE_GAS_API void ASPIGH(int Ng, double P, double H, double* ff, int* nff, int* err);
HE_GAS_API void AGHmax(int Ng, double* Hmin, double* Hmax);
HE_GAS_API void AGUmax(int Ng, double* Umin, double* Umax);

// --- сервис (необязательно) ---
// Задать файлы ДО первого расчета (NULL/"" - по умолчанию: переменные среды HE_TABLE_FILE / HE_COEF_FILE,
// затем he_table2.bin и data/he_coef_abinitio.txt). Возвращает 0, либо -1, если модель уже загружена.
HE_GAS_API int         ASPIGInit(const char* table_path, const char* coef_path);
HE_GAS_API int         AGLastFlags(void);        // флаги he::FL_* последнего расчета гелия в этом потоке
HE_GAS_API void        AGSetStrict(int strict);  // 1: вне паспортного диапазона err = 1
HE_GAS_API const char* AGLastError(void);        // причина, если модель гелия не загрузилась ("" - нет ошибки)
HE_GAS_API const char* AGModelInfo(void);        // источник модели гелия ("table: ...", "analytic: ...")

#ifdef __cplusplus
namespace he { class Model; }
const he::Model* AGGetHeModel();                  // модель гелия (nullptr, если не загрузилась); C++ only
#endif
