# BNYPC Ballistic Engine

High-performance ballistic trajectory simulator written in C using RK4 numerical integration and G1 drag model.

---

## 🚀 Features

- RK4 numerical integration (4th order Runge-Kutta)
- G1 drag table interpolation
- Mach-dependent drag model
- Wind drift simulation (crosswind + headwind)
- Atmospheric density correction
- Full muzzle velocity trajectory solver

---

## ⚙️ Build

### Clang (recommended)
```bash
clang src/main.c -o ballistic -O2 -march=native -lm

### GCC :
gcc src/main.c -o ballistic -O2 -lm 