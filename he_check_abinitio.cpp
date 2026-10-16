// he_check_abinitio.cpp - проверка вязкости, теплопроводности (а также Pr и B) по ab initio-значениям
// для 4He в пределе нулевой плотности: E. Bich, R. Hellmann, E. Vogel, Mol. Phys. 105 (2007) 3035,
// таблица 3. Заявленная авторами неопределенность вязкости и теплопроводности – 0.02 % выше 15 К.
//
//   he_check_abinitio data/he4_abinitio_bich2007.csv [--model table:ФАЙЛ | analytic[:коэф.txt]]
//                     [--p 10] [--tol 0.01] [--csv-out файл] [--mmap]
//
// Модель вычисляется при давлении --p (по умолчанию 10 Па – нижняя граница диапазона): поправка на
// плотность при этом порядка 1e-7, т. е. результат сопоставим с пределом нулевой плотности.
// Сравниваются: mu, lambda, Pr = cp*mu/lambda (для идеального газа cp = 5R/2) и второй вириальный
// коэффициент B, вычисленный из сжимаемости модели при том же давлении: B = (Z - 1) R T / P.
// Критерий (код возврата 1 при невыполнении): максимальное отклонение mu, lambda и Pr в паспортном
// диапазоне 50-3000 К не более --tol (по умолчанию 1 %). Отклонение B выводится как справочное.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "he_eos.h"
#include "he_table.h"

using namespace he;

struct Row { double T, B, C, eta, lam; };

struct Stat {
  double maxd = 0, sum = 0, Tmax = 0;
  int n = 0;
  void add(double d, double T) { const double a = std::fabs(d); sum += a; ++n; if (a > maxd) { maxd = a; Tmax = T; } }
};

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "Использование: %s данные.csv [--model table:ФАЙЛ|analytic[:коэф.txt]] [--p 10] [--tol 0.01] [--csv-out файл] [--mmap]\n", argv[0]);
    return 1;
  }
  std::string data = argv[1], model = "analytic", csvout;
  double P = P_MIN, tol = 0.01;
  unsigned opts = 0;
  for (int i = 2; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--model" && i + 1 < argc) model = argv[++i];
    else if (a == "--p" && i + 1 < argc) P = std::atof(argv[++i]);
    else if (a == "--tol" && i + 1 < argc) tol = std::atof(argv[++i]);
    else if (a == "--csv-out" && i + 1 < argc) csvout = argv[++i];
    else if (a == "--mmap") opts |= OPEN_MMAP;
    else { std::fprintf(stderr, "неизвестный параметр %s\n", a.c_str()); return 1; }
  }

  std::vector<Row> rows;
  {
    FILE* f = std::fopen(data.c_str(), "r");
    if (!f) { std::fprintf(stderr, "Ошибка: не удалось открыть %s\n", data.c_str()); return 2; }
    char line[512];
    while (std::fgets(line, sizeof line, f)) {
      Row r;
      if (line[0] == '#' || !(line[0] == '-' || (line[0] >= '0' && line[0] <= '9'))) continue;
      if (std::sscanf(line, "%lf,%lf,%lf,%lf,%lf", &r.T, &r.B, &r.C, &r.eta, &r.lam) == 5) rows.push_back(r);
    }
    std::fclose(f);
  }
  if (rows.empty()) { std::fprintf(stderr, "Ошибка: в %s нет данных\n", data.c_str()); return 2; }

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
  } else { std::fprintf(stderr, "Ошибка: --model table:ФАЙЛ или analytic[:файл]\n"); return 1; }

  const double Ru = 8.314462618;
  FILE* fo = csvout.empty() ? nullptr : std::fopen(csvout.c_str(), "w");
  if (fo) std::fprintf(fo, "T,eta_ab,eta_model,eta_dev,lambda_ab,lambda_model,lambda_dev,Pr_ab,Pr_model,Pr_dev,B_ab,B_model,B_dev,flags\n");

  std::printf("Модель: %s | P = %g Па | данные: %s (%zu точек)\n\n", m->name(), P, data.c_str(), rows.size());
  std::printf("%8s | %9s %9s %7s | %9s %9s %7s | %7s %7s %7s | %8s %8s %7s\n", "T, K", "mu_ab", "mu_mod", "dev,%", "lam_ab", "lam_mod", "dev,%", "Pr_ab", "Pr_mod", "dev,%", "B_ab", "B_mod", "dev,%");
  std::printf("%8s | %9s %9s %7s | %9s %9s %7s | %7s %7s %7s | %8s %8s %7s\n", "", "мкПа·с", "мкПа·с", "", "мВт/м/К", "мВт/м/К", "", "", "", "", "см3/моль", "см3/моль", "");
  Stat sm[2], sl[2], sp[2], sb[2];  // [0] – паспортный диапазон 50-3000 К, [1] – вне его (экстраполяция)
  int skipped = 0;
  for (const Row& r : rows) {
    Result o;
    const int e = m->PT(P, r.T, HP_RHO | HP_MU | HP_LAM | HP_PR, o);
    if (e != HE_OK) { ++skipped; continue; }
    const double mu_ab = r.eta * 1e-6, lam_ab = r.lam * 1e-3;
    const double Pr_ab = 2.5 * R_SPEC * mu_ab / lam_ab;
    const double Z = P / (o.rho * R_SPEC * r.T);
    const double Bm = (Z - 1.0) * Ru * r.T / P * 1e6;
    const double dm = o.mu / mu_ab - 1, dl = o.lambda / lam_ab - 1, dp = o.Pr / Pr_ab - 1, db = Bm / r.B - 1;
    const int g = (r.T >= T_MIN && r.T <= T_MAX) ? 0 : 1;
    sm[g].add(dm, r.T); sl[g].add(dl, r.T); sp[g].add(dp, r.T);
    if (std::fabs(r.B) > 1.0) sb[g].add(db, r.T);  // у нуля B относительное отклонение не имеет смысла
    std::printf("%8.2f | %9.4f %9.4f %+7.3f | %9.3f %9.3f %+7.3f | %7.4f %7.4f %+7.3f | %8.3f %8.3f %+7.2f%s\n", r.T, r.eta, o.mu * 1e6, dm * 100, r.lam, o.lambda * 1e3, dl * 100,
                Pr_ab, o.Pr, dp * 100, r.B, Bm, db * 100, g ? "  *" : "");
    if (fo) std::fprintf(fo, "%g,%.6g,%.8g,%.6e,%.6g,%.8g,%.6e,%.6g,%.8g,%.6e,%.6g,%.8g,%.6e,%u\n", r.T, r.eta, o.mu * 1e6, dm, r.lam, o.lambda * 1e3, dl, Pr_ab, o.Pr, dp, r.B, Bm, db, o.flags);
  }
  if (fo) std::fclose(fo);
  std::printf("\n* – вне паспортного диапазона 50-3000 К (экстраполяция модели). Не вычислено точек: %d (за жесткими границами).\n\n", skipped);

  auto line = [](const char* name, const Stat& s) {
    if (s.n) std::printf("  %-10s max |dev| = %7.3f %%  (T = %7.2f К), среднее |dev| = %7.3f %%, точек %d\n", name, s.maxd * 100, s.Tmax, s.sum / s.n * 100, s.n);
  };
  for (int g = 0; g < 2; ++g) {
    std::printf(g == 0 ? "Паспортный диапазон 50-3000 К:\n" : "Вне паспортного диапазона (20-50 К и 3000-6000 К):\n");
    line("вязкость", sm[g]); line("теплопровод.", sl[g]); line("Pr", sp[g]); line("B (справ.)", sb[g]);
  }
  const double worst = std::fmax(sm[0].maxd, std::fmax(sl[0].maxd, sp[0].maxd));
  std::printf("\nКритерий: max |dev| вязкости, теплопроводности и Pr в 50-3000 К <= %.3g %%: %s (максимум %.3f %%)\n", tol * 100, worst <= tol ? "ВЫПОЛНЕН" : "НЕ ВЫПОЛНЕН", worst * 100);
  std::printf("Справка: авторы статьи заявляют неопределенность значений 0.02 %% выше 15 К.\n");
  if (sb[0].n && sb[0].maxd > 0.01)
    std::printf("ПРЕДУПРЕЖДЕНИЕ: второй вириальный коэффициент B модели отличается от ab initio более чем на 1 %% (до %.2f %% при T = %.0f К); это влияет на ρ, c_p и U при высоких давлениях.\n", sb[0].maxd * 100, sb[0].Tmax);
  return worst <= tol ? 0 : 1;
}
