// he_calc.cpp - Программа 2 (версия 2): расчет свойств гелия-4 по (P,T), (P,H) или (P,U)
//
//   he_calc [--model table:ФАЙЛ | analytic[:коэф.txt]] [--mmap] [--mask список|all]
//           PT|PH|PU  P  X                 одна точка (P в Па; X = T в К, H или U в Дж/кг)
//   he_calc ... PT|PH|PU -f файл|-         пары "P X" из файла/stdin, вывод CSV
//
// Маска: список через запятую из T,rho,cp,cv,U,H,lambda,mu,Pr,dUdT,drhodT,drhodU,drhodH,drhodP или all.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

#include "he_eos.h"
#include "he_table.h"
#include "he_console.h"

using namespace he;

struct PropDef { const char* name; unsigned bit; const char* unit; };
static const PropDef kProps[] = {
    {"T", HP_T, "К"},           {"rho", HP_RHO, "кг/м3"},    {"cp", HP_CP, "Дж/(кг К)"},   {"cv", HP_CV, "Дж/(кг К)"},
    {"U", HP_U, "Дж/кг"},       {"H", HP_H, "Дж/кг"},        {"lambda", HP_LAM, "Вт/(м К)"}, {"mu", HP_MU, "Па с"},
    {"Pr", HP_PR, "-"},         {"dUdT", HP_DUDT, "Дж/(кг К)"}, {"drhodT", HP_DRDT, "кг/(м3 К)"},
    {"drhodU", HP_DRDU, "кг/м3 / (Дж/кг)"}, {"drhodH", HP_DRDH, "кг/м3 / (Дж/кг)"}, {"drhodP", HP_DRDP, "кг/(м3 Па)"}};

static double val(const Result& r, unsigned bit) {
  switch (bit) {
    case HP_T: return r.T; case HP_RHO: return r.rho; case HP_CP: return r.cp; case HP_CV: return r.cv;
    case HP_U: return r.U; case HP_H: return r.H; case HP_LAM: return r.lambda; case HP_MU: return r.mu;
    case HP_PR: return r.Pr; case HP_DUDT: return r.dUdT; case HP_DRDT: return r.drhodT; case HP_DRDU: return r.drhodU;
    case HP_DRDH: return r.drhodH; default: return r.drhodP;
  }
}

static bool parse_mask(const char* s, unsigned& m) {
  if (!std::strcmp(s, "all")) { m = HP_ALL; return true; }
  m = 0;
  std::string t(s);
  size_t p = 0;
  while (p <= t.size()) {
    size_t q = t.find(',', p);
    if (q == std::string::npos) q = t.size();
    const std::string w = t.substr(p, q - p);
    bool f = false;
    for (auto& d : kProps) if (w == d.name) { m |= d.bit; f = true; }
    if (!f) return false;
    p = q + 1;
  }
  return m != 0;
}

static void usage(const char* a) {
  std::fprintf(stderr,
               "Использование:\n  %s [--model table:ФАЙЛ|analytic[:коэф.txt]] [--mmap] [--mask список|all] PT|PH|PU P X\n"
               "  %s [...] PT|PH|PU -f файл|-      (пары \"P X\", вывод CSV)\n"
               "P – Па; X – T [К] для PT, H [Дж/кг] для PH, U [Дж/кг] для PU.\n", a, a);
}

int main(int argc, char** argv) {
  std::string model = "analytic";
  unsigned mask = HP_ALL, opts = 0;
  int k = 1;
  for (; k < argc && argv[k][0] == '-' && std::strcmp(argv[k], "-") != 0 && std::strcmp(argv[k], "-f") != 0; ++k) {
    if (!std::strcmp(argv[k], "--model") && k + 1 < argc) model = argv[++k];
    else if (!std::strcmp(argv[k], "--mask") && k + 1 < argc) {
      if (!parse_mask(argv[++k], mask)) { std::fprintf(stderr, "Ошибка: неверная маска\n"); return 1; }
    } else if (!std::strcmp(argv[k], "--mmap")) opts |= OPEN_MMAP;
    else { usage(argv[0]); return 1; }
  }
  if (k + 1 >= argc) { usage(argv[0]); return 1; }
  const std::string mode = argv[k++];
  if (mode != "PT" && mode != "PH" && mode != "PU") { usage(argv[0]); return 1; }

  std::unique_ptr<Model> m;
  std::string err;
  if (model.rfind("table:", 0) == 0) {
    auto t = TableModel::open(model.c_str() + 6, err, opts);
    if (!t) { std::fprintf(stderr, "Ошибка: %s\n", err.c_str()); return 2; }
    m = std::move(t);
  } else if (model.rfind("analytic", 0) == 0) {
    auto a = std::make_unique<AnalyticModel>();
    if (model.size() > 9 && model[8] == ':' && !a->load(model.c_str() + 9, &err)) { std::fprintf(stderr, "Ошибка: %s\n", err.c_str()); return 2; }
    m = std::move(a);
  } else { usage(argv[0]); return 1; }

  auto call = [&](double P, double X, Result& r, Hint* h) {
    return mode == "PT" ? m->PT(P, X, mask, r, h) : mode == "PH" ? m->PH(P, X, mask, r, h) : m->PU(P, X, mask, r, h);
  };
  const char* in = mode == "PT" ? "T" : (mode == "PH" ? "H" : "U");

  if (!std::strcmp(argv[k], "-f") || !std::strcmp(argv[k], "-")) {
    FILE* f = std::strcmp(argv[k], "-") == 0 ? stdin : (k + 1 < argc ? std::fopen(argv[k + 1], "r") : nullptr);
    if (!f) { std::fprintf(stderr, "Ошибка: не удалось открыть входной файл\n"); return 2; }
    std::printf("P,%s", mode == "PT" ? "T" : (mode == "PH" ? "H" : "U"));
    unsigned shown = mask | (mode == "PT" ? 0u : HP_T);
    for (auto& d : kProps) if ((shown & d.bit) && !(mode == "PT" && d.bit == HP_T) && !(mode == "PH" && d.bit == HP_H) && !(mode == "PU" && d.bit == HP_U)) std::printf(",%s", d.name);
    std::printf(",flags,status\n");
    char line[512];
    Hint h;
    int bad = 0;
    while (std::fgets(line, sizeof line, f)) {
      double P, X;
      if (line[0] == '#' || std::sscanf(line, "%lf %lf", &P, &X) != 2) continue;
      Result r;
      const int e = call(P, X, r, &h);
      std::printf("%.10g,%.10g", P, X);
      for (auto& d : kProps)
        if ((shown & d.bit) && !(mode == "PT" && d.bit == HP_T) && !(mode == "PH" && d.bit == HP_H) && !(mode == "PU" && d.bit == HP_U))
          std::printf(",%.10e", e ? NAN : val(r, d.bit));
      std::printf(",%u,%d\n", e ? 0u : r.flags, e);
      bad += e != 0;
    }
    if (f != stdin) std::fclose(f);
    return bad ? 1 : 0;
  }

  if (k + 1 >= argc) { usage(argv[0]); return 1; }
  char *e1, *e2;
  const double P = std::strtod(argv[k], &e1), X = std::strtod(argv[k + 1], &e2);
  if (*e1 || *e2) { std::fprintf(stderr, "Ошибка: P и %s должны быть числами\n", in); return 1; }
  Result r;
  const int e = call(P, X, r, nullptr);
  if (e) {
    std::fprintf(stderr, "Ошибка %d: %s\n", e, e == HE_BAD_INPUT ? "недопустимые входные данные (вне жестких границ)" : e == HE_NO_SOLUTION ? "для заданных (P, X) нет температуры в допустимом интервале" : "сбой модели");
    return 3;
  }
  std::printf("модель: %s | вызов %s | P = %.8g Па, %s = %.10g\n", m->name(), mode.c_str(), P, in, X);
  for (auto& d : kProps)
    if (r.computed & d.bit) std::printf("  %-7s = %.10e  %s\n", d.name, val(r, d.bit), d.unit);
  if (mode != "PT") std::printf("  итераций по T: %d\n", r.iters);
  if (r.flags) std::printf("  ВНИМАНИЕ: вне паспортного диапазона (экстраполяция), флаги = %u\n", r.flags);
  return 0;
}
