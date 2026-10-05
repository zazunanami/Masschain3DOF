# Report Errata

Corrections to [Report.pdf](Report.pdf), checked against the model in this repository. The PDF's source document is not part of the repository, so the PDF itself is unchanged.

## How the model values were obtained

- Chain: `m = 1 kg` and `k = 1000 N/m` everywhere unless the row says otherwise; total time 5 s, `dt = 0.001 s`; peak values are `max |x|` over the run.
- Undamped studies use the modal solver with the pulse converted to `Δv = F·Δt/m` (Auto mode with damping off). The mass studies use an impulse of 1 N·s on m1 (0.5 J at `m1 = 1 kg`); all other studies use 10 N·s (`F = 1000 N`, `Δt = 0.01 s`, 50 J).
- Damped studies use the RK4 solver with the force pulse (`F = 1000 N`, `Δt = 0.01 s`) on m1, the same damping value on every dashpot or mass, and `ε = 0.001 m/s`.
- "Upper bound" is the largest `|x|` the undamped response can ever reach: `Σ_j |Φ_ij| · |q̇_j(0)| / ω_j`.

The solver fixes made in this repository do not change these scenarios: for the 0.01 s pulse, damped and Coulomb runs differ from the original code by less than 3e-8 m. The original code showed a single-sample 1.6 J notch in the energy curve where the pulse ends, which no longer appears.

## Numerical claims that do not match the model

| Section | Report | Model | Assessment |
|---|---|---|---|
| m3 study, `m3 = 0.5 kg` (Table 3, text) | Tip moves most, `x3 ≈ 2.8 cm`, `x3 > x1` | `x1 = 2.57`, `x2 = 2.48`, `x3 = 2.21 cm`; upper bound for x3 is 2.24 cm | Incorrect: the tip moves least, and 2.8 cm is unreachable |
| m3 study, `m3 = 2.0 kg` | `x3 ≈ 1.6 cm` | `x3 = 1.86 cm` | About 14 % low |
| m1 study, `m1 = 1.0 kg` | Uniform, `≈ 2.4 cm`, `x1 ≈ x2` | `x1 = 2.55`, `x2 = 2.07`, `x3 = 2.49 cm` | x2 is about 19 % below x1 |
| k0 study (Table 4) | `x1 ≈ 0.28 / 0.20 / 0.15 m` for `k0 = 500 / 1000 / 2000 N/m` | `0.316 / 0.255 / 0.198 m` | Incorrect, 11–24 % low |
| Reference state S0 | `x1 ≈ 0.20 m` (Table 4) vs `0.24 m` (Table 14) | `0.255 m` | The report contradicts itself; the Table 4 value is wrong |
| k2 study, `k2 = 500 N/m` | m3 swings most, `x3 ≈ 0.26 m` | `x3 = 0.214 m`; x2 has the peak at `0.272 m` | Incorrect |
| k3 study (Table 7) | `x3 ≈ 0.26 / 0.16 m` for `k3 = 500 / 2000 N/m` | `0.280 / 0.195 m` | About 7 % and 18 % low |
| Viscous, `c = 0.5 N·s/m` (Table 15) | More than 10 % energy left at 5 s | 6.3 % of peak energy | Incorrect |
| Viscous, `c = 2.0 N·s/m` (Table 15) | Fewer than 5 cycles in 0–5 s | About 16 cycles above 5 % of the peak (20 zero-crossing cycles) | Incorrect |
| Viscous, all cases | `c = 2` described as a "shock absorber" that can "freeze" the system | Modal damping ratios are 0.6–1.5 %, 1.2–2.9 % and 2.4–5.8 % for `c = 0.5, 1, 2` | All three cases are lightly damped; the wording overstates them |

Claims that do match the model (within reading accuracy): m1 = 0.5 / 2.0 kg and m2 = 2.0 kg peaks, the excitation-location peaks (0.255 / 0.292 / 0.255 m), the input energies (`I²/2m`), the Coulomb energies at 5 s (27.4 / 12.0 / 0 J), the Coulomb stop time for `Fc = 2 N` (4.8 s), and the near-complete viscous decay for `c = 2` by 2.5 s (1.4 % left).

## Physics and reasoning corrections

- **Coulomb decay shape.** The report says energy drops linearly with time ("perfect triangle"). Under Coulomb friction the amplitude decays linearly, so energy (∝ amplitude²) decays along a parabola. In the model, `√E(t)` fits a straight line with R² = 1.000, while `E(t)` gives R² = 0.976 for `Fc = 2 N`.
- **Stiction.** The `tanh(v/ε)` friction model has no true sticking phase; masses near rest creep slowly instead of locking.
- **Impulse-duration study.** "Energy ∝ Δt² verifies the accuracy of the time-integration solver" is not a valid check. In Auto or Modal mode the pulse is converted to an initial velocity, so `E = I²/2m` holds by construction and the solver is never tested. The RK4 force-pulse solver gives peak energies of 12.45 / 49.2 / 187.1 J for `Δt = 0.005 / 0.01 / 0.02 s`. A 0.02 s pulse is not short compared with the fastest period (0.107 s), so response and energy are no longer exactly proportional to the impulse.

## Calculation and text errors in the fabrication section

- **q1, Case 1.** The method states `Δt ≈ 1.5 T = 0.6667 s`, but the frequency is computed as `1/0.6667 = 1.50 Hz`. With `T = 0.6667 / 1.5 = 0.444 s`, the result is `f ≈ 2.25 Hz` and `ω ≈ 14.1 rad/s`, not 1.50 Hz and 9.42 rad/s.
- **q1, Case 2.** The text says "This motion is happening in 0.523809 second", but 0.523809 s is `T`; the observed interval is `Δt = 0.785714 s`.
- **Case 2 heading.** It first says the impulse is applied on M2, then on M1.
- **Measured vs. model frequencies.** With the measured stiffnesses (55.6 / 88.3 / 93.8 / 79.0 N/m) and masses, the model gives 3.13 / 6.18 / 8.54 Hz for Case 1 and 4.23 / 23.6 / 25.8 Hz for Case 2. The measured 1.5–2.25 Hz are all below the lowest predicted mode. The report never compares the two, so the conclusion that the observations are "consistent with the expected behavior" is not supported quantitatively.

## Notation and labelling

- The layout line `Wall - k1 - m1 - k2 - m2 - k3 - m3 - k4 - Wall` numbers the springs k1–k4, while the equations, tables and the app use k0–k3.
- In the m2 damper derivative, "Term 3 (c3)" should read "Term 3 (c2)".
- The heading "3-DOF IMPULSE RESPONSE: MIDDLE MASS SENSITIVITY (m3)" names m3, which is the tip mass.
- The subject line "Analysis of Wall Stiffness / Root Constraint" is reused for the k1 and k2 studies; both are coupling springs, not wall springs.
- State labels S9 and S10 in Table 5 repeat labels already used for `m3 = 2.0 kg` (Table 3) and `k0 = 500 N/m` (Table 4).
- The figure number 43 is used twice (`c = 1` and `c = 2 N·s/m`).
