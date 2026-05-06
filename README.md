# 2D Ising Model — Monte Carlo Simulation

A C++ implementation of the **2D Ising model** on a 50×50 square lattice with periodic boundary conditions, simulated via the **Metropolis–Hastings Monte Carlo algorithm**. The program sweeps a range of temperatures, tracks magnetization across MC steps, and compares the equilibrium magnetization to **Onsager's analytic solution** for the spontaneous magnetization of the infinite 2D Ising model.

## What it does

1. Initializes an `L × L` lattice of spins (`±1`) at random.
2. For each temperature in `T = 0.1, 0.2, ..., 2.0` (in units of `J / k_B`):
   - Runs 100 Metropolis MC sweeps over the lattice.
   - At each site, attempts a spin flip and accepts with probability `min(1, exp(-ΔE / k_B T))`.
   - Records the average magnetization after every sweep.
   - Computes a time-averaged magnetization over the post-equilibration sweeps (last 50).
3. Evaluates Onsager's closed-form magnetization
   `m(T) = [1 − sinh⁻⁴(2J / k_B T)]^(1/8)`
   for the same temperature grid.
4. Writes everything to plain-text files for plotting in your tool of choice (gnuplot, Python, MATLAB, etc.).

## Physics notes

- Energy per spin: `E = −J · sᵢⱼ · (sum of 4 nearest neighbors)`, with the standard ½ factor to avoid double-counting bond energies.
- Units: `J = 1`, `k_B = 1`. Temperature is therefore in units of `J / k_B`.
- The exact 2D Ising critical temperature is `T_c = 2J / [k_B · ln(1 + √2)] ≈ 2.269`. The simulation sweeps from deep in the ordered phase (`T = 0.1`) up to just below `T_c`, so a sharp transition won't appear within this range — extending past `T = 2.5` will show it.
- Periodic boundary conditions are applied via modular arithmetic on lattice indices.

## Files generated

- `magnetization_vs_steps_T<temperature>.txt` — one file per temperature, two columns: MC step index and instantaneous magnetization.
- `onsager_magnetization_vs_temperature.txt` — temperature vs. Onsager analytic magnetization.

Each file is tab-separated and ready for `plot` or `loadtxt`.

## Build and run

No external dependencies — only the C++ standard library.

```bash
g++ -O2 -std=c++17 -o ising ising.cpp
./ising
```

Output files appear in the working directory.

## Quick plotting (Python)

```python
import numpy as np, matplotlib.pyplot as plt

T, m_onsager = np.loadtxt("onsager_magnetization_vs_temperature.txt", unpack=True)
plt.plot(T, m_onsager, "k--", label="Onsager (analytic)")
plt.xlabel("Temperature  (J / k_B)")
plt.ylabel("Magnetization per spin")
plt.legend(); plt.grid(); plt.show()
```

## Caveats and possible improvements

- **RNG quality.** `rand()` seeded with `time(NULL)` is fine for a class project but biased and low-period. Switching to `<random>` (`std::mt19937` with `std::uniform_real_distribution`) is a small change that improves statistical reliability.
- **Reseeding.** `srand` is currently called inside `createGrid()` only; subsequent random numbers come from the same seed. With a `std::mt19937` engine made global, this becomes cleaner.
- **Equilibration.** 50 sweeps of equilibration on a 50×50 lattice is short, especially near criticality where correlation times diverge. For temperatures near `T_c`, increase `MC_steps` and `equilibrium_steps` significantly.
- **Temperature range.** Sweeps stop at `T = 2.0`, just below `T_c ≈ 2.269`. To see the phase transition clearly, extend to `T = 3.5` or so.
- **Time-averaged magnetization is computed but not written out.** Adding a file `time_avg_magnetization_vs_temperature.txt` would let you overlay the simulated equilibrium magnetization on the Onsager curve directly.
- **`calculateDeltaEnergy` copies the entire grid** to flip one spin. The standard fast form is `ΔE = 2 · J · sᵢⱼ · (sum of 4 NN spins)`, which avoids the copy and is several orders of magnitude faster.

## Parameters

Edit at the top of the file or in `main()`:

| Constant | Value | Meaning |
|---|---|---|
| `L` | 50 | Lattice side length |
| `J` | 1 | Coupling constant |
| `kB` | 1 | Boltzmann constant |
| `MC_steps` | 100 | Total MC sweeps per temperature |
| `equilibrium_steps` | 50 | Sweeps discarded before time-averaging |
| `temperatures` | 0.1 → 2.0 | Temperature grid (units of `J / k_B`) |

