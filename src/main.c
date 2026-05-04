/*
    ██████╗ ███╗   ██╗██╗   ██╗██████╗  ██████╗
    ██╔══██╗████╗  ██║╚██╗ ██╔╝██╔══██╗██╔════╝
    ██████╔╝██╔██╗ ██║ ╚████╔╝ ██████╔╝██║
    ██╔══██╗██║╚██╗██║  ╚██╔╝  ██╔═══╝ ██║
    ██████╔╝██║ ╚████║   ██║   ██║     ╚██████╗
    ╚═════╝ ╚═╝  ╚═══╝   ╚═╝   ╚═╝      ╚═════╝

        B N Y P C  -  Ballistics Engine

    Autor: Emanoel Rodrigues Costa Albuquerque

    Projeto: Calculadora Balística G1 (RK4)

    Destaques:
      - Integração numérica RK4 de alta precisão
      - Modelo G1 com interpolação
      - Simulação de vento lateral e frontal
      - Correção por densidade do ar e Mach

      Testado Com:
      Compilação em: Clang 22.1.4
      MSVC(Visual Studio 2026)

      Flags recomendadas:
      -O2 -march=native -ffast-math

    "Código limpo, trajetória precisa."
*/

#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GRAVITY 9.80665
#define RHO0 1.2250
#define SPEED_SOUND_0 340.29
#define PI 3.14159265358979323846

typedef struct {
  double mach;
  double cd;
} DragPoint;

static const DragPoint G1_TABLE[] = {
    {0.00, 0.2629}, {0.05, 0.2558},  {0.10, 0.2487}, {0.15, 0.2413},
    {0.20, 0.2335}, {0.25, 0.2258},  {0.30, 0.2184}, {0.35, 0.2113},
    {0.40, 0.2046}, {0.45, 0.1983},  {0.50, 0.1956}, {0.55, 0.1944},
    {0.60, 0.1935}, {0.65, 0.1930},  {0.70, 0.1928}, {0.725, 0.1928},
    {0.75, 0.1932}, {0.775, 0.1943}, {0.80, 0.1965}, {0.825, 0.1990},
    {0.85, 0.2028}, {0.875, 0.2076}, {0.90, 0.2138}, {0.925, 0.2224},
    {0.95, 0.2339}, {0.975, 0.2484}, {1.00, 0.2659}, {1.025, 0.2873},
    {1.05, 0.3118}, {1.10, 0.3607},  {1.15, 0.4043}, {1.20, 0.4418},
    {1.25, 0.4734}, {1.30, 0.4998},  {1.35, 0.5217}, {1.40, 0.5398},
    {1.45, 0.5546}, {1.50, 0.5664},  {1.55, 0.5756}, {1.60, 0.5828},
    {1.65, 0.5884}, {1.70, 0.5928},  {1.75, 0.5963}, {1.80, 0.5992},
    {1.85, 0.6016}, {1.90, 0.6037},  {1.95, 0.6057}, {2.00, 0.6074},
    {2.05, 0.6089}, {2.10, 0.6101},  {2.20, 0.6120}, {2.30, 0.6132},
    {2.40, 0.6140}, {2.50, 0.6145},  {2.60, 0.6148}, {2.80, 0.6150},
    {3.00, 0.6150}, {3.50, 0.6150},  {4.00, 0.6149}, {4.50, 0.6148},
    {5.00, 0.6130},
};
static const int G1_TABLE_SIZE = sizeof(G1_TABLE) / sizeof(G1_TABLE[0]);

typedef struct {
  double muzzle_vel;
  double bc_g1;
  double mass_g;
  double diameter_mm;
  double zero_range_m;
  double wind_speed;
  double wind_angle;
  double altitude_m;
} BulletParams;

typedef struct {
  double x;
  double y;
  double vx;
  double vy;
  double t;
  double wz;
  double vwz;
} State;

static double g1_cd(double mach) {
  if (mach <= G1_TABLE[0].mach)
    return G1_TABLE[0].cd;
  if (mach >= G1_TABLE[G1_TABLE_SIZE - 1].mach)
    return G1_TABLE[G1_TABLE_SIZE - 1].cd;

  for (int i = 0; i < G1_TABLE_SIZE - 1; i++) {
    if (mach >= G1_TABLE[i].mach && mach <= G1_TABLE[i + 1].mach) {
      double t =
          (mach - G1_TABLE[i].mach) / (G1_TABLE[i + 1].mach - G1_TABLE[i].mach);
      return G1_TABLE[i].cd + t * (G1_TABLE[i + 1].cd - G1_TABLE[i].cd);
    }
  }
  return G1_TABLE[G1_TABLE_SIZE - 1].cd;
}

static double density_ratio(double altitude_m) {

  double T = 288.15 - 0.0065 * altitude_m;
  double T0 = 288.15;
  return pow(T / T0, 4.256);
}

static double speed_of_sound(double altitude_m) {
  double T = 288.15 - 0.0065 * altitude_m;
  return 331.3 * sqrt(T / 273.15);
}

static double drag_accel(double v, double bc_g1, double altitude_m) {
  double dratio = density_ratio(altitude_m);
  double vs = speed_of_sound(altitude_m);
  double mach = v / vs;
  double cd = g1_cd(mach);
  double bc_si = bc_g1 * 703.0697;

  return (RHO0 * dratio * v * v * cd) / (2.0 * bc_si);
}

typedef struct {
  double dx, dy, dvx, dvy, dt, dwz, dvwz;
} Deriv;

static Deriv derivatives(const State *s, const BulletParams *p,
                         double crosswind_ms) {
  Deriv d;
  double v = sqrt(s->vx * s->vx + s->vy * s->vy);
  double a = (v > 1e-9) ? drag_accel(v, p->bc_g1, p->altitude_m) : 0.0;

  d.dx = s->vx;
  d.dy = s->vy;
  d.dvx = -(a * s->vx / v);
  d.dvy = -(a * s->vy / v) - GRAVITY;
  d.dt = 1.0;

  double rel_vwz = s->vwz - crosswind_ms;
  double v_tot = sqrt(v * v + rel_vwz * rel_vwz);
  double a_lat = (v_tot > 1e-9) ? drag_accel(v_tot, p->bc_g1, p->altitude_m) *
                                      (rel_vwz / v_tot)
                                : 0.0;
  d.dwz = s->vwz;
  d.dvwz = -a_lat;

  return d;
}

static State rk4_step(const State *s, const BulletParams *p, double crosswind,
                      double h) {
  Deriv k1 = derivatives(s, p, crosswind);

  State s2 = {s->x + 0.5 * h * k1.dx,
              s->y + 0.5 * h * k1.dy,
              s->vx + 0.5 * h * k1.dvx,
              s->vy + 0.5 * h * k1.dvy,
              s->t + 0.5 * h,
              s->wz + 0.5 * h * k1.dwz,
              s->vwz + 0.5 * h * k1.dvwz};
  Deriv k2 = derivatives(&s2, p, crosswind);

  State s3 = {s->x + 0.5 * h * k2.dx,
              s->y + 0.5 * h * k2.dy,
              s->vx + 0.5 * h * k2.dvx,
              s->vy + 0.5 * h * k2.dvy,
              s->t + 0.5 * h,
              s->wz + 0.5 * h * k2.dwz,
              s->vwz + 0.5 * h * k2.dvwz};
  Deriv k3 = derivatives(&s3, p, crosswind);

  State s4 = {s->x + h * k3.dx,    s->y + h * k3.dy, s->vx + h * k3.dvx,
              s->vy + h * k3.dvy,  s->t + h,         s->wz + h * k3.dwz,
              s->vwz + h * k3.dvwz};
  Deriv k4 = derivatives(&s4, p, crosswind);

  State next;
  next.x = s->x + (h / 6.0) * (k1.dx + 2 * k2.dx + 2 * k3.dx + k4.dx);
  next.y = s->y + (h / 6.0) * (k1.dy + 2 * k2.dy + 2 * k3.dy + k4.dy);
  next.vx = s->vx + (h / 6.0) * (k1.dvx + 2 * k2.dvx + 2 * k3.dvx + k4.dvx);
  next.vy = s->vy + (h / 6.0) * (k1.dvy + 2 * k2.dvy + 2 * k3.dvy + k4.dvy);
  next.t = s->t + h;
  next.wz = s->wz + (h / 6.0) * (k1.dwz + 2 * k2.dwz + 2 * k3.dwz + k4.dwz);
  next.vwz =
      s->vwz + (h / 6.0) * (k1.dvwz + 2 * k2.dvwz + 2 * k3.dvwz + k4.dvwz);
  return next;
}

static double y_at_range(double elev_rad, const BulletParams *p,
                         double target_range, double step) {
  double crosswind = p->wind_speed * sin(p->wind_angle * PI / 180.0);
  State s = {.x = 0,
             .y = 0,
             .vx = p->muzzle_vel * cos(elev_rad),
             .vy = p->muzzle_vel * sin(elev_rad),
             .t = 0,
             .wz = 0,
             .vwz = 0};
  while (s.x < target_range && s.y > -500.0) {
    s = rk4_step(&s, p, crosswind, step);
  }
  return s.y;
}

static double find_zero_elevation(const BulletParams *p) {
  double lo = -0.05, hi = 0.10;
  double step = 0.001;

  for (int i = 0; i < 40; i++) {
    double mid = 0.5 * (lo + hi);
    double y = y_at_range(mid, p, p->zero_range_m, step);
    if (y > 0)
      hi = mid;
    else
      lo = mid;
  }
  return 0.5 * (lo + hi);
}

static double read_double(const char *prompt, double default_val) {
  char buf[64];
  printf("  %s [%.4g]: ", prompt, default_val);
  fflush(stdout);
  if (fgets(buf, sizeof(buf), stdin) && buf[0] != '\n') {
    double v;
    if (sscanf(buf, "%lf", &v) == 1)
      return v;
  }
  return default_val;
}

static void read_params(BulletParams *p) {
  printf("\n+--------------------------------------+\n");
  printf("|    CALCULADORA BALISTICA  (G1/RK4)   |\n");
  printf("+--------------------------------------+\n");
  printf("\n-- Dados do Projetil -----------------\n");
  p->muzzle_vel = read_double("Velocidade de boca (m/s)", 900.0);
  p->bc_g1 = read_double("Coef. Balistico G1 (lb/in2)", 0.475);
  p->mass_g = read_double("Massa do projetil (gramas)", 10.9);
  p->diameter_mm = read_double("Diametro do projetil (mm)", 7.62);

  printf("\n-- Configuracao de Tiro --------------\n");
  p->zero_range_m = read_double("Zeragem (metros)", 100.0);
  p->altitude_m = read_double("Altitude do atirador (metros)", 0.0);

  printf("\n-- Vento -----------------------------\n");
  p->wind_speed = read_double("Velocidade do vento (m/s)", 0.0);
  p->wind_angle = read_double("Angulo do vento (90=cruzado, 0=frontal)", 90.0);
}

static void print_table(const BulletParams *p, double elev_rad) {
  double crosswind = p->wind_speed * sin(p->wind_angle * PI / 180.0);
  double headwind = p->wind_speed * cos(p->wind_angle * PI / 180.0);
  double mass_kg = p->mass_g / 1000.0;
  double step = 0.001;
  double table_step = 50.0;

  State s = {.x = 0,
             .y = 0,
             .vx = (p->muzzle_vel + headwind) * cos(elev_rad),
             .vy = (p->muzzle_vel + headwind) * sin(elev_rad),
             .t = 0,
             .wz = 0,
             .vwz = 0};

  double elev_mrad = elev_rad * 1000.0;
  double elev_moa = elev_rad * (180.0 / PI) * 60.0;

  printf("\n-- Solucao de Tiro -----------------------------------\n");
  printf("  Elevacao de boca   : %+.3f mrad  (%+.2f MOA)\n", elev_mrad,
         elev_moa);
  printf("  Zeragem            : %.0f m\n", p->zero_range_m);
  printf("  Vento cruzado      : %.1f m/s\n", crosswind);
  printf("  Altitude           : %.0f m (densidade: %.4f)\n", p->altitude_m,
         density_ratio(p->altitude_m));
  printf("  Vel. de boca efet. : %.1f m/s\n", p->muzzle_vel);
  printf("  Vel. do som local  : %.1f m/s\n", speed_of_sound(p->altitude_m));

  printf("\n%s\n", "+------+----------+----------+----------+---------+--------"
                   "-+----------+");
  printf("|%-6s|%-10s|%-10s|%-10s|%-9s|%-9s|%-10s|\n", " Alc.", "  Drop",
         " Drop", "  Vento", "  Vel.", "Mach", "  Energia");
  printf("|%-6s|%-10s|%-10s|%-10s|%-9s|%-9s|%-10s|\n", "  (m)", "  (cm)",
         " (MOA)", "  (cm)", "  (m/s)", "     ", "   (J)");
  printf("%s\n", "+------+----------+----------+----------+---------+---------+"
                 "----------+");

  double next_print = 0.0;

  while (s.x <= 1500.0 && s.y > -200.0) {
    if (s.x >= next_print - 1e-6) {
      double v = sqrt(s.vx * s.vx + s.vy * s.vy);
      double mach = v / speed_of_sound(p->altitude_m);
      double energy = 0.5 * mass_kg * v * v;

      double drop_cm = s.y * 100.0;

      double drop_moa =
          (s.x > 1e-3) ? atan2(-s.y, s.x) * (180.0 / PI) * 60.0 : 0.0;
      double wind_cm = s.wz * 100.0;

      printf("|%5.0f |%+9.1f |%+9.2f |%+9.1f |%8.1f |%8.3f |%9.1f |\n", s.x,
             drop_cm, drop_moa, wind_cm, v, mach, energy);

      next_print += table_step;
    }
    s = rk4_step(&s, p, crosswind, step);
  }

  printf("%s\n", "+------+----------+----------+----------+---------+---------+"
                 "----------+");
  printf("\n  Legenda:\n");
  printf("  Drop  = queda em relacao a linha de visada (negativo = abaixo)\n");
  printf("  MOA   = Minutes of Angle (1 MOA ~ 2.91 cm/100m)\n");
  printf("  Vento = deriva lateral (positivo = para a direita)\n\n");
}

int main(void) {
  BulletParams p;
  read_params(&p);

  printf("\n  Calculando trajetoria");
  fflush(stdout);

  double elev = find_zero_elevation(&p);
  printf(" - OK\n");

  print_table(&p, elev);

  printf("Pressione ENTER para sair...\n");
  getchar();
  return 0;
}
