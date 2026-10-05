# Report Corrections

[Report.pdf](Report.pdf) was corrected in place in October 2026. Every claim below was checked against the model in this repository before the PDF was changed. The last page of the report ("Sanitization and Technical Correction Notes") summarises the same corrections.

## How the model values were obtained

- Chain: `m = 1 kg` and `k = 1000 N/m` everywhere unless a row says otherwise; total time 5 s, `dt = 0.001 s`; peak values are `max |x|` over the run.
- Undamped studies use the modal solver with the pulse converted to `Δv = F·Δt/m` (Auto mode with damping off). The mass studies use an impulse of 1 N·s on m1 (0.5 J at `m1 = 1 kg`); all other studies use 10 N·s (`F = 1000 N`, `Δt = 0.01 s`, 50 J).
- Damped studies use the RK4 solver with the force pulse (`F = 1000 N`, `Δt = 0.01 s`) on m1, the same damping value on every dashpot or mass, and `ε = 0.001 m/s`.
- "Upper bound" is the largest `|x|` the undamped response can ever reach: `Σ_j |Φ_ij| · |q̇_j(0)| / ω_j`.

The figures (simulation screenshots) are unchanged. The solver fixes made in this repository do not change the scenarios they show: with the 0.01 s pulse, damped and Coulomb runs differ from the original code by less than 3e-8 m. The only visible difference is that the original code drew a single-sample, 1.6 J notch in the energy curve at the end of the pulse (t = 0.01 s).

## Numerical claims

| Page | Report said | Model | Corrected to |
|---|---|---|---|
| 16–17 (m3 = 0.5 kg, Table 3) | Tip moves most, `x3 ≈ 2.8 cm`, `x3 > x1`, "whipping effect", use a light tip to amplify sensor signals | `x1 = 2.57`, `x2 = 2.48`, `x3 = 2.21 cm`; upper bound for x3 is 2.24 cm | "Tip-Reduced (x3 < x1)": the tip has the smallest peak (≈2.2 cm vs ≈2.6 cm at x1); a light tip does not amplify tip motion |
| 17 (m3 = 2.0 kg) | `x3 ≈ 1.6 cm`; x1 and x2 both grow | `x3 = 1.86 cm`; x2 grows from 2.07 to 2.46 cm, x1 stays at 2.55 cm | ≈1.9 cm; only the middle mass grows (≈2.1 → ≈2.5 cm) |
| 10, 11, 13 (reference state, Tables 1–2) | Uniform `≈ 2.4 cm`, `x1 ≈ x2` | `x1 = 2.55`, `x2 = 2.07`, `x3 = 2.49 cm` | `x1 ≈ x3`; peaks ≈2.1–2.6 cm |
| 19–20 (k0, Table 4) | `x1 ≈ 0.28 / 0.20 / 0.15 m` for `k0 = 500 / 1000 / 2000 N/m` | `0.316 / 0.255 / 0.198 m` | 0.32 / 0.26 / 0.20 m (text: 32 cm, 20 cm) |
| 22 (k1, Table 5) | `k1 = 2000 N/m`: x1 "High (Rebound)" | x1 falls from 0.255 m to 0.214 m | "Low (Coupled)" |
| 25–26 (k2, Table 6) | `k2 = 500 N/m`: m3 acts like a "flail", `x3 ≈ 0.26 m`, out of phase with m2 | `x3 = 0.214 m` (lower than 0.249 m); x2 has the peak, 0.272 m; corr(x2, x3) = +0.39 | Mass 3 is only weakly driven (≈0.21 m); Mass 2 has the largest peak (≈0.27 m); Table 6 "Low / Delayed" |
| 28–29 (k3, Table 7) | `x3 ≈ 0.26 / 0.16 m` for `k3 = 500 / 2000 N/m`; `k3 = 2000` raises x1 | `0.280 / 0.195 m`; x1 = 0.254 m vs 0.255 m | 0.28 / 0.20 m; x1 stays practically the same |
| 49–50 (excitation location, Table 14) | `0.24 / 0.29 / 0.25 m`; tip excitation gives a different modal participation and peak; place actuators at the root | `0.255 / 0.292 / 0.255 m`; the equal chain is symmetric, so tip excitation is the mirror image of root excitation | 0.26 / 0.29 / 0.26 m; mirror-image statement; actuators near either end |
| 52 (viscous, Table 15) | Oscillations: ">15", "moderate", "<5 cycles"; more than 10 % and ~0 % energy left at 5 s | About 20 cycles in 0–5 s for every `c`; 6.3 % (`c = 0.5`) and 1.4 % (`c = 1`) energy left | "≈20 cycles" with slow/fast decay; "≈6 %" and "≈1 %" |

Claims that already matched the model were left unchanged: the m1 = 0.5 / 2.0 kg and m2 = 2.0 kg peaks, the input energies (`I²/2m`), the Coulomb energies at 5 s (27.4 / 12.0 / 0 J), the Coulomb stop time for `Fc = 2 N` (4.8 s), and the near-complete viscous decay for `c = 2` by 2.5 s (1.4 % left).

## Physics and reasoning

- **Coulomb decay shape (page 55).** The report said energy drops linearly with time ("perfect triangle"). Under Coulomb friction the amplitude decays linearly, so energy (∝ amplitude²) follows a parabola. In the model, `√E(t)` fits a straight line with R² = 1.000, while `E(t)` gives R² = 0.976 for `Fc = 2 N`. The text now says this, and "energy line crashing straight down" became "energy curve falling steeply".
- **Stiction (page 55).** A note was added that the smoothed `tanh` friction model only approximates sticking; masses creep instead of locking.
- **Impulse-duration validation (page 47).** "Energy ∝ Δt² verifies the time-integration solver" was replaced by the correct statement. In Auto or Modal mode the pulse is converted to an initial velocity, so `E = I²/2m` holds by construction and only that conversion is checked. The RK4 force-pulse solver gives 12.45 / 49.2 / 187.1 J for `Δt = 0.005 / 0.01 / 0.02 s`.

## Calculations in the fabrication section

- **q1, Case 1 (page 62).** The method states `Δt = 1.5 T = 0.6667 s`, so `T = 0.444 s`. The frequency is therefore 2.25 Hz and ω = 14.14 rad/s, not 1.50 Hz and 9.42 rad/s.
- **q1, Case 2 (page 64).** The observed interval is 0.785714 s; 0.523809 s is `T`.
- **Case 2 heading (page 63).** The impulse is applied on M2, as the case title says.
- **Case 2 q1 (page 64).** After the Case 1 correction, Case 2 q1 (1.91 Hz) is lower than Case 1 (2.25 Hz), which is opposite to the mass–stiffness trend. The sentence now says so.
- **Measured vs. model frequencies (notes page).** With the measured stiffnesses (55.6 / 88.3 / 93.8 / 79.0 N/m) and masses, the model gives 3.13 / 6.18 / 8.54 Hz for Case 1 and 4.23 / 23.6 / 25.8 Hz for Case 2. The measured 1.5–2.25 Hz are all below the lowest predicted mode, so the agreement is qualitative only.

## Equations, notation and labels

- **Page 5.** The determinant expansion had a duplicated `(−k1)` factor and an unbalanced bracket in `(k1)²[[…`; both were removed.
- **Page 4.** The layout line now numbers the springs k0–k3, as the equations and the app do.
- **Page 7.** "Term 3 (c3)" became "Term 3 (c2)".
- **Page 15.** The heading now reads "TIP MASS SENSITIVITY (m3)", since m3 is the tip mass.
- **Pages 21, 24, 27.** The k1 and k2 subject lines say "Coupling Stiffness"; the k3 subject line says "Wall Stiffness / End Constraint".
- **Pages 22–23.** The k1 states were renamed S36/S37, because S9 and S10 were already used for `m3 = 2.0 kg` and `k0 = 500 N/m`.
- **Pages 51, 53, 54.** Figures were renumbered 45–48; numbers 43 and 44 had been used twice.

## How the PDF was edited

The text was replaced in place using the fonts already embedded in the report (Cambria Math, Calibri), subset and re-embedded under the font's editable-embedding permission. Word's simulated bold (fill plus stroke) and its script-size glyph variants were reproduced. Unchanged characters on an edited line were redrawn at their original positions. Text extraction and search work on every edited line, and all 34 unedited pages render pixel-identically to the previous version.
