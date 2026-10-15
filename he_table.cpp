// he_table.cpp - Вариант A: вычисление свойств по таблице потенциала G = ln(phi)
#include "he_table.h"

#include <cstdio>
#include <cstring>

#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#define HE_HAVE_MMAP 1
#endif

namespace he {

enum : unsigned { W_GX = 1, W_GY = 2, W_GXX = 4, W_GXY = 8, W_GYY = 16, W_ALL = 31 };

static inline unsigned wants(unsigned need) {
  unsigned w = 0;
  if (need & N_RHO) w |= W_GY;
  if (need & N_H) w |= W_GX;
  if (need & N_CP) w |= W_GX | W_GXX;
  if (need & N_DV) w |= W_GY | W_GXY | W_GYY;
  return w;
}

// g = {G_x, G_y, G_xx, G_xy, G_yy}; значения G не нужны – только производные
static inline void eval_G(const double* c, double u, double v, double ihx, double ihy, unsigned w, double* g) {
  double U0[6], U1[6], U2[6], V0[6], V1[6], V2[6];
  U0[0] = V0[0] = 1.0;
  for (int k = 1; k < 6; ++k) { U0[k] = U0[k - 1] * u; V0[k] = V0[k - 1] * v; }
  U1[0] = V1[0] = 0.0;
  U2[0] = U2[1] = V2[0] = V2[1] = 0.0;
  for (int k = 1; k < 6; ++k) { U1[k] = k * U0[k - 1]; V1[k] = k * V0[k - 1]; }
  for (int k = 2; k < 6; ++k) { U2[k] = k * (k - 1) * U0[k - 2]; V2[k] = k * (k - 1) * V0[k - 2]; }
  const bool a = (w & (W_GX | W_GXX)) != 0, b = (w & (W_GY | W_GXY)) != 0, d = (w & W_GYY) != 0;
  double gx = 0, gy = 0, gxx = 0, gxy = 0, gyy = 0;
  for (int i = 0; i < 6; ++i) {
    const double* ci = c + 6 * i;
    double r0 = 0, r1 = 0, r2 = 0;
    if (a) for (int j = 0; j < 6; ++j) r0 += ci[j] * V0[j];
    if (b) for (int j = 0; j < 6; ++j) r1 += ci[j] * V1[j];
    if (d) for (int j = 0; j < 6; ++j) r2 += ci[j] * V2[j];
    gx += r0 * U1[i];
    gxx += r0 * U2[i];
    gy += r1 * U0[i];
    gxy += r1 * U1[i];
    gyy += r2 * U0[i];
  }
  g[0] = gx * ihx;
  g[1] = gy * ihy;
  g[2] = gxx * ihx * ihx;
  g[3] = gxy * ihx * ihy;
  g[4] = gyy * ihy * ihy;
}

void leaf_thermo(const LeafView& lv, double xc, double yc, double dx, double dy, double T, double P,
                 unsigned need, Thermo& s) {
  const double u = (xc - lv.x0) * lv.ihx, v = (yc - lv.y0) * lv.ihy;
  const bool ext = (dx != 0.0 || dy != 0.0);
  double g[5];
  eval_G(lv.c, u, v, lv.ihx, lv.ihy, ext ? W_ALL : wants(need), g);
  if (ext) {
    // экстраполяция: Taylor по T (ln T); по P вниз – закон G ~ P (вириальный предел), вверх – Taylor
    if (dx != 0.0) {
      g[0] += g[2] * dx;
      g[1] += g[3] * dx;
    }
    if (dy < 0.0) {
      const double f = std::exp(dy);
      for (int k = 0; k < 5; ++k) g[k] *= f;
    } else if (dy > 0.0) {
      g[1] += g[4] * dy;
      g[0] += g[3] * dy;
    }
  }
  const double R = R_SPEC, RT = R * T, opg = 1.0 + g[1];
  if (need & N_RHO) s.rho = P / (RT * opg);
  if (need & N_H) s.h = RT * (2.5 - g[0]);
  if (need & N_CP) s.cp = R * (2.5 - g[0] - g[2]);
  if (need & N_DV) {
    s.dvdT = (R / P) * (opg + g[3]);
    s.dvdP = (RT / (P * P)) * (g[4] - opg);
  }
}

void leaf_transport(const LeafView& lv, double xc, double yc, double dx, double dy, unsigned which, double& lam,
                    double& mu) {
  const double u = (xc - lv.x0) * lv.ihx, v = (yc - lv.y0) * lv.ihy;
  double U[4] = {1, u, u * u, u * u * u}, V[4] = {1, v, v * v, v * v * v};
  const bool ext = (dx != 0.0 || dy != 0.0);
  double dU[4] = {0, 1, 2 * u, 3 * u * u}, dV[4] = {0, 1, 2 * v, 3 * v * v};
  for (int k = 0; k < 2; ++k) {
    const bool want = (k == 0) ? (which & (HP_LAM | HP_PR)) != 0 : (which & (HP_MU | HP_PR)) != 0;
    if (!want) continue;
    const double* t = lv.tr + 16 * k;
    double f = 0, fx = 0, fy = 0;
    for (int a = 0; a < 4; ++a) {
      double r = 0, rd = 0;
      for (int b = 0; b < 4; ++b) { r += t[4 * a + b] * V[b]; if (ext) rd += t[4 * a + b] * dV[b]; }
      f += r * U[a];
      if (ext) { fx += r * dU[a]; fy += rd * U[a]; }
    }
    if (ext) f += fx * lv.ihx * dx + fy * lv.ihy * dy;
    (k == 0 ? lam : mu) = std::exp(f);
  }
}

uint64_t table_checksum(const LeafRec* rec, size_t n) {
  const unsigned char* p = reinterpret_cast<const unsigned char*>(rec);
  const size_t len = n * sizeof(LeafRec);
  uint64_t h = 1469598103934665603ULL;
  for (size_t k = 0; k < len; ++k) { h ^= p[k]; h *= 1099511628211ULL; }
  return h;
}

TableModel::~TableModel() {
#ifdef HE_HAVE_MMAP
  if (map_) munmap(map_, maplen_);
#endif
}

void TableModel::setup() {
  X0_ = std::log(hdr_.Tmin); X1_ = std::log(hdr_.Tmax);
  Y0_ = std::log(hdr_.Pmin); Y1_ = std::log(hdr_.Pmax);
  Lm_ = (int)hdr_.max_level;
  nxg_ = hdr_.rx << Lm_;
  nyg_ = hdr_.ry << Lm_;
  sx_ = nxg_ / (X1_ - X0_);
  sy_ = nyg_ / (Y1_ - Y0_);
  idx_.assign((size_t)nxg_ * nyg_, 0);
  geom_.resize(n_);
  for (size_t k = 0; k < n_; ++k) {
    const LeafRec& r = rec_[k];
    const double hx = (X1_ - X0_) / (double)(hdr_.rx << r.level), hy = (Y1_ - Y0_) / (double)(hdr_.ry << r.level);
    geom_[k] = {X0_ + r.i * hx, Y0_ + r.j * hy, 1.0 / hx, 1.0 / hy};
    const unsigned sh = Lm_ - r.level;
    for (uint32_t a = r.i << sh; a < (r.i + 1) << sh; ++a)
      for (uint32_t b = r.j << sh; b < (r.j + 1) << sh; ++b) idx_[(size_t)a * nyg_ + b] = (uint32_t)k;
  }
}

LeafView TableModel::view(int id) const {
  const Geom& g = geom_[id];
  return {g.x0, g.y0, g.ihx, g.ihy, rec_[id].c, &rec_[id].tr[0][0]};
}

int TableModel::locate(double x, double y, Hint* h) const {
  if (h && h->cell >= 0 && (size_t)h->cell < n_) {
    const Geom& g = geom_[h->cell];
    const double u = (x - g.x0) * g.ihx, v = (y - g.y0) * g.ihy;
    if (u >= 0.0 && u < 1.0 && v >= 0.0 && v < 1.0) return h->cell;
  }
  int ix = (int)((x - X0_) * sx_), iy = (int)((y - Y0_) * sy_);
  ix = ix < 0 ? 0 : (ix >= (int)nxg_ ? (int)nxg_ - 1 : ix);
  iy = iy < 0 ? 0 : (iy >= (int)nyg_ ? (int)nyg_ - 1 : iy);
  return (int)idx_[(size_t)ix * nyg_ + iy];
}

int TableModel::thermo(double T, double P, unsigned need, Thermo& s, Hint* h) const {
  const double xr = std::log(T), yr = std::log(P);
  const double xc = xr < X0_ ? X0_ : (xr > X1_ ? X1_ : xr);
  const double yc = yr < Y0_ ? Y0_ : (yr > Y1_ ? Y1_ : yr);
  const int id = locate(xc, yc, h);
  if (h) h->cell = id;
  const LeafView lv = view(id);
  leaf_thermo(lv, xc, yc, xr - xc, yr - yc, T, P, need, s);
  return HE_OK;
}

int TableModel::transport(double T, double P, double, unsigned which, double& lam, double& mu, Hint* h) const {
  const double xr = std::log(T), yr = std::log(P);
  const double xc = xr < X0_ ? X0_ : (xr > X1_ ? X1_ : xr);
  const double yc = yr < Y0_ ? Y0_ : (yr > Y1_ ? Y1_ : yr);
  const int id = locate(xc, yc, h);
  if (h) h->cell = id;
  leaf_transport(view(id), xc, yc, xr - xc, yr - yc, which, lam, mu);
  return HE_OK;
}

std::unique_ptr<TableModel> TableModel::from_leaves(TableHeader hdr, std::vector<LeafRec>&& leaves) {
  std::unique_ptr<TableModel> m(new TableModel());
  m->own_ = std::move(leaves);
  m->rec_ = m->own_.data();
  m->n_ = m->own_.size();
  uint32_t ml = 0;
  for (const LeafRec& r : m->own_) ml = r.level > ml ? r.level : ml;
  hdr.max_level = ml;
  hdr.nleaf = (uint32_t)m->n_;
  hdr.checksum = table_checksum(m->rec_, m->n_);
  m->hdr_ = hdr;
  m->setup();
  return m;
}

bool TableModel::save(const char* path) const {
  FILE* f = std::fopen(path, "wb");
  if (!f) return false;
  const bool ok = std::fwrite(&hdr_, sizeof hdr_, 1, f) == 1 && std::fwrite(rec_, sizeof(LeafRec), n_, f) == n_;
  return std::fclose(f) == 0 && ok;
}

std::unique_ptr<TableModel> TableModel::open(const char* path, std::string& err, unsigned opts) {
  std::unique_ptr<TableModel> m(new TableModel());
  FILE* f = std::fopen(path, "rb");
  if (!f) { err = "не удалось открыть файл таблицы"; return nullptr; }
  TableHeader h;
  if (std::fread(&h, sizeof h, 1, f) != 1) { std::fclose(f); err = "файл слишком короткий"; return nullptr; }
  if (std::memcmp(h.magic, "HE4TBL02", 8) != 0) { std::fclose(f); err = "неверная сигнатура (нужна таблица версии 2)"; return nullptr; }
  if (h.endian != 0x01020304u) { std::fclose(f); err = "порядок байтов файла не совпадает с порядком байтов машины"; return nullptr; }
  if (h.version != 2 || h.header_size != sizeof(TableHeader) || h.leaf_size != sizeof(LeafRec) || h.nleaf == 0 ||
      h.nleaf > 50000000u || h.max_level > 20 || h.rx == 0 || h.ry == 0) {
    std::fclose(f); err = "несовместимый формат файла"; return nullptr;
  }
  m->hdr_ = h;
  m->n_ = h.nleaf;
  bool mapped = false;
#ifdef HE_HAVE_MMAP
  if (opts & OPEN_MMAP) {
    const int fd = fileno(f);
    struct stat st;
    if (fstat(fd, &st) == 0 && (size_t)st.st_size == sizeof(TableHeader) + h.nleaf * sizeof(LeafRec)) {
      void* p = mmap(nullptr, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
      if (p != MAP_FAILED) {
        m->map_ = p;
        m->maplen_ = (size_t)st.st_size;
        m->rec_ = reinterpret_cast<const LeafRec*>(static_cast<const char*>(p) + sizeof(TableHeader));
        if (opts & OPEN_MLOCK) mlock(p, m->maplen_);  // ошибка (лимиты ОС) не критична
        mapped = true;
      }
    }
  }
#endif
  if (!mapped) {
    m->own_.resize(h.nleaf);
    if (std::fread(m->own_.data(), sizeof(LeafRec), h.nleaf, f) != h.nleaf) {
      std::fclose(f); err = "файл усечен"; return nullptr;
    }
    m->rec_ = m->own_.data();
  }
  std::fclose(f);
  if (!(opts & OPEN_NOCHECK) && table_checksum(m->rec_, m->n_) != h.checksum) {
    err = "контрольная сумма не совпадает – файл поврежден";
    return nullptr;
  }
  for (size_t k = 0; k < m->n_; ++k)
    if (m->rec_[k].level > h.max_level || m->rec_[k].i >= (h.rx << m->rec_[k].level) ||
        m->rec_[k].j >= (h.ry << m->rec_[k].level)) { err = "некорректные данные листа"; return nullptr; }
  m->setup();
  return m;
}

}  // namespace he
