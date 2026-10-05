# MassChain3DOF

MassChain3DOF is an educational Tkinter application for simulating a three-degree-of-freedom mass-spring chain. It includes modal analysis, undamped modal superposition, RK4 time-domain integration, optional viscous damping, optional smoothed Coulomb friction, and Matplotlib-based displacement/energy plots.

## Technologies

- Python 3.10+
- Tkinter
- NumPy
- Matplotlib
- Pytest for lightweight verification tests

## Features

- Builds 3x3 mass, stiffness, and damping matrices for a wall-connected three-mass chain.
- Solves the generalized modal problem and reports natural frequencies and mass-normalized mode shapes.
- Simulates undamped free vibration using modal superposition.
- Simulates damped/nonlinear time response using a fourth-order Runge-Kutta integrator.
- Supports a rectangular force pulse or equivalent initial-velocity interpretation.
- Shows displacement, mechanical energy, and a simple animated mass-chain view.

## Setup

Using `uv`:

```bash
uv venv
uv pip install -p .venv/bin/python -e ".[dev]"
```

Using standard Python tooling:

```bash
python -m venv .venv
.venv/bin/python -m pip install -e ".[dev]"
```

On Windows PowerShell, use:

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -e ".[dev]"
```

Tkinter must be available in the Python installation. On some Linux distributions it is provided by a separate system package such as `python3-tk`.

## Run

```bash
python -m masschain3dof
```

Or run the application module directly:

```bash
python src/masschain3dof/app.py
```

## Test

```bash
python -m compileall src tests
python -m pytest -q
```

The included tests verify the equal-chain mass/stiffness matrices, natural frequencies, mass normalization, undamped modal energy conservation, energy reduction under viscous damping, RK4 agreement with the modal solution on non-divisible time grids, exact impulse transfer of the force pulse, agreement with the closed-form rectangular-pulse response, and the RK4 stability estimate.

## Project Structure

```text
src/masschain3dof/
  app.py        # Tkinter GUI and simulation logic
  __main__.py   # `python -m masschain3dof` entry point
  __init__.py

tests/
  test_modal_core.py
  test_time_domain.py
```

## Known Limitations

This project is an educational numerical simulation. It is not certified for professional, safety-critical, or production engineering use. Results should be independently verified before any real-world use.

The Coulomb friction model uses a smooth `tanh(v/eps)` approximation to avoid a discontinuity at zero velocity. Near zero velocity this model is very stiff (slope `Fc/eps`), so small `eps` or large `Fc` values require a small time step; the app warns when the selected time step is outside the RK4 stability limit of the linearized model and refuses to plot diverged (non-finite) results. The displayed mechanical energy is kinetic plus spring potential energy; dissipated energy is not separately accumulated.

The time grid always ends at the requested total time. If the total time is not an integer multiple of the time step, the step is shortened slightly so that every stored sample matches its plotted time. The rectangular force pulse acts on `0 <= t < Δt`; the RK4 step containing the pulse end is split there, so the applied impulse is exactly `F·Δt`.

## Reports and Documents

A sanitized public copy of the project report is included at [docs/report/Report.pdf](docs/report/Report.pdf). Personal identifiers and course-submission metadata were removed from it; its numerical figures were produced by an earlier version of the code.

## License

Copyright (c) 2026 Zazu Nanami

The source code in this repository is licensed under the Apache License 2.0. See [LICENSE](LICENSE).
