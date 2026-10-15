// he_table.h - Вариант A: таблица остаточного потенциала Гиббса с интерполяцией класса C2.
//
// Хранится не набор свойств, а ОДНА функция  G(x, y) = g^r(T,P)/(R T) = ln(phi),  x = ln T, y = ln P
// (остаточная часть энергии Гиббса; идеально-газовая часть – аналитическая). Все свойства
// получаются дифференцированием интерполянта, поэтому тождества термодинамики выполняются точно:
//   Z = 1 + G_y,    h = R T (5/2 - G_x),    cp = R (5/2 - G_x - G_xx),
//   (dv/dT)_P = (R/P)(1 + G_y + G_xy),    (dv/dP)_T = (R T/P^2)(G_yy - 1 - G_y).
// Интерполянт – биквинтический эрмитов полином в каждой ячейке (36 коэффициентов),
// данные в узлах: G и производные до порядка (2,2) => G принадлежит C2 внутри согласованных областей.
// Сетка – квадродерево (локальное измельчение). На стыках ячеек разных уровней непрерывность
// выполняется с точностью до ошибки интерполяции (скачок контролируется построителем).
// Переносные свойства: ln(mu), ln(lambda) – бикубический эрмитов полином (C1) в тех же ячейках.
#pragma once
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "he_api.h"

namespace he {

struct LeafRec {          // запись листа в файле (560 байт)
  uint32_t level, i, j, pad;
  double c[36];           // G(u,v) = sum c[6*a+b] u^a v^b,  u,v in [0,1]
  double tr[2][16];       // ln lambda, ln mu: sum tr[k][4*a+b] u^a v^b
};
static_assert(sizeof(LeafRec) == 560, "LeafRec size");

struct TableHeader {      // заголовок файла (160 байт)
  char magic[8];          // "HE4TBL02"
  uint32_t endian;        // 0x01020304 – порядок байтов
  uint32_t version;       // 2
  uint32_t header_size, leaf_size;
  uint32_t rx, ry;        // число корневых ячеек по ln T и ln P
  uint32_t max_level;     // максимальный уровень листьев
  uint32_t nleaf;
  double Tmin, Tmax, Pmin, Pmax;
  double tol;             // допуск построения
  double verified_err;    // макс. относительная ошибка по контролю построителя
  uint64_t checksum;      // FNV-1a по массиву листьев
  char ref_name[64];
};
static_assert(sizeof(TableHeader) == 160, "TableHeader size");

struct LeafView { double x0, y0, ihx, ihy; const double* c; const double* tr; };

// Вычисление в листе (используется и построителем). (xc,yc) внутри листа; dx, dy – смещение
// исходной точки от (xc,yc) за границу области (экстраполяция), в координатах ln T, ln P.
void leaf_thermo(const LeafView& lv, double xc, double yc, double dx, double dy, double T, double P,
                 unsigned need, Thermo& s);
void leaf_transport(const LeafView& lv, double xc, double yc, double dx, double dy, unsigned which,
                    double& lam, double& mu);
uint64_t table_checksum(const LeafRec* rec, size_t n);

enum OpenOpt : unsigned { OPEN_MMAP = 1, OPEN_MLOCK = 2, OPEN_NOCHECK = 4 };

class TableModel : public Model {
 public:
  ~TableModel();
  static std::unique_ptr<TableModel> open(const char* path, std::string& err, unsigned opts = 0);
  static std::unique_ptr<TableModel> from_leaves(TableHeader hdr, std::vector<LeafRec>&& leaves);
  bool save(const char* path) const;

  const char* name() const override { return "table (C2 potential, quadtree)"; }
  int thermo(double T, double P, unsigned need, Thermo& s, Hint* h) const override;
  int transport(double T, double P, double rho, unsigned which, double& lam, double& mu, Hint* h) const override;

  const TableHeader& header() const { return hdr_; }
  size_t n_leaves() const { return n_; }
  const LeafRec& rec(size_t k) const { return rec_[k]; }
  LeafView view(int id) const;
  int locate(double x, double y, Hint* h) const;
  double X0() const { return X0_; }
  double X1() const { return X1_; }
  double Y0() const { return Y0_; }
  double Y1() const { return Y1_; }
  size_t bytes() const { return sizeof(TableHeader) + n_ * sizeof(LeafRec); }

 private:
  TableModel() {}
  void setup();
  struct Geom { double x0, y0, ihx, ihy; };
  TableHeader hdr_{};
  std::vector<LeafRec> own_;
  const LeafRec* rec_ = nullptr;
  size_t n_ = 0;
  void* map_ = nullptr;
  size_t maplen_ = 0;
  std::vector<Geom> geom_;
  std::vector<uint32_t> idx_;
  int Lm_ = 0;
  uint32_t nxg_ = 0, nyg_ = 0;
  double X0_, X1_, Y0_, Y1_, sx_, sy_;
};

}  // namespace he
