import numpy as np
import pytest

from masschain3dof.app import (
    force_pulse,
    linearized_state_eigenvalues,
    make_mass_matrix,
    make_stiffness_matrix,
    make_time_grid,
    rk4_amplification,
    run_modal_sim,
    run_time_sim,
    solve_modal,
    uses_time_domain,
)


EQUAL_MASSES = np.array([1.0, 1.0, 1.0])
EQUAL_SPRINGS = np.array([1000.0, 1000.0, 1000.0, 1000.0])


def simulate(masses, springs, v0, F_imp, imp_dt, t_end, dt, dashpots=None, fc_vec=None, frictionEps=0.001):
    return run_time_sim(
        masses=masses,
        springs=springs,
        dashpots=np.zeros(4) if dashpots is None else dashpots,
        x0=np.zeros(3),
        v0=v0,
        imp_idx=0,
        F_imp=F_imp,
        imp_dt=imp_dt,
        fc_vec=np.zeros(3) if fc_vec is None else fc_vec,
        frictionEps=frictionEps,
        t_end=t_end,
        dt=dt,
    )


def exact_pulse_response(masses, springs, force, duration, t):
    """Closed-form undamped response to a rectangular pulse on mass 1 (step at 0 minus step at duration)."""
    omega, modes = solve_modal(make_mass_matrix(masses), make_stiffness_matrix(springs))
    modal_force = modes[0, :] * force
    q = np.zeros((len(t), 3))
    for j, w in enumerate(omega):
        static = modal_force[j] / w**2
        during = static * (1.0 - np.cos(w * t))
        after = static * (np.cos(w * (t - duration)) - np.cos(w * t))
        q[:, j] = np.where(t <= duration, during, after)
    return q @ modes.T


@pytest.mark.parametrize(
    "t_end, dt",
    [(1.0, 0.001), (1.0, 0.003), (0.7, 0.002), (0.3, 0.1), (1.0, 1e-5), (0.05, 0.2)],
)
def test_time_grid_ends_at_t_end_with_uniform_step_not_above_dt(t_end, dt):
    t = make_time_grid(t_end, dt)
    steps = np.diff(t)

    assert t[0] == 0.0
    assert t[-1] == pytest.approx(t_end, abs=1e-15)
    assert np.allclose(steps, steps[0], rtol=1e-9)
    assert steps[0] <= dt * (1.0 + 1e-9)


@pytest.mark.parametrize("t_end, dt", [(1.0, 0.003), (0.7, 0.002)])
def test_undamped_rk4_matches_modal_solution_on_returned_time_axis(t_end, dt):
    v0 = np.array([1.0, 0.0, 0.0])

    t_rk4, x_rk4, *_ = simulate(EQUAL_MASSES, EQUAL_SPRINGS, v0, 0.0, 0.0, t_end, dt)
    t_modal, x_modal, *_ = run_modal_sim(EQUAL_MASSES, EQUAL_SPRINGS, np.zeros(3), v0, t_end, dt)

    assert np.array_equal(t_rk4, t_modal)
    assert np.max(np.abs(x_rk4 - x_modal)) < 1e-5


@pytest.mark.parametrize(
    "imp_dt, dt",
    [(0.01, 0.001), (0.0005, 0.001), (0.0105, 0.001), (0.01, 0.003), (0.2, 0.05)],
)
def test_force_pulse_transfers_exact_impulse_to_free_mass(imp_dt, dt):
    masses = np.array([2.0, 1.0, 1.0])
    force = 1000.0

    _, _, v_hist, *_ = simulate(masses, np.zeros(4), np.zeros(3), force, imp_dt, 0.3, dt)

    assert v_hist[-1, 0] == pytest.approx(force * imp_dt / masses[0], rel=1e-12)
    assert np.all(v_hist[:, 1:] == 0.0)


def test_force_pulse_response_matches_closed_form_when_pulse_end_is_off_grid():
    force = 1000.0
    duration = 0.0105

    t, x_hist, *_ = simulate(EQUAL_MASSES, EQUAL_SPRINGS, np.zeros(3), force, duration, 0.5, 0.001)

    x_exact = exact_pulse_response(EQUAL_MASSES, EQUAL_SPRINGS, force, duration, t)
    assert np.max(np.abs(x_hist - x_exact)) < 1e-5 * np.max(np.abs(x_exact))


def test_force_pulse_is_active_on_half_open_interval():
    assert force_pulse(0.0, 1, 5.0, 0.01)[1] == 5.0
    assert force_pulse(0.00999, 1, 5.0, 0.01)[1] == 5.0
    assert force_pulse(0.01, 1, 5.0, 0.01)[1] == 0.0
    assert np.all(force_pulse(0.0, 1, 5.0, 0.0) == 0.0)


def rk4_growth(dt, dashpots=np.zeros(4), fc_vec=np.zeros(3), frictionEps=0.001):
    lam = linearized_state_eigenvalues(EQUAL_MASSES, EQUAL_SPRINGS, dashpots, fc_vec, frictionEps)
    return np.max(np.abs(rk4_amplification(lam * dt)))


def test_rk4_stability_estimate_matches_actual_divergence_for_heavy_viscous_damping():
    dashpots = np.full(4, 1000.0)
    v0 = np.array([1.0, 0.0, 0.0])

    assert rk4_growth(0.001, dashpots=dashpots) > 1.0
    with np.errstate(over="ignore", invalid="ignore"):
        _, _, _, energy_unstable, *_ = simulate(EQUAL_MASSES, EQUAL_SPRINGS, v0, 0.0, 0.0, 1.0, 0.001, dashpots=dashpots)
    assert not np.all(np.isfinite(energy_unstable)) or energy_unstable[-1] > energy_unstable[0]

    assert rk4_growth(0.0005, dashpots=dashpots) <= 1.0
    _, _, _, energy_stable, *_ = simulate(EQUAL_MASSES, EQUAL_SPRINGS, v0, 0.0, 0.0, 1.0, 0.0005, dashpots=dashpots)
    assert np.all(np.isfinite(energy_stable))
    assert energy_stable[-1] < energy_stable[0]


def test_rk4_stability_estimate_accounts_for_smoothed_coulomb_friction():
    assert rk4_growth(0.001) <= 1.0
    assert rk4_growth(0.001, fc_vec=np.full(3, 1.0)) <= 1.0
    assert rk4_growth(0.001, fc_vec=np.full(3, 20.0)) > 1.0
    assert rk4_growth(0.001, fc_vec=np.full(3, 20.0), frictionEps=0.01) <= 1.0


def test_solver_selection():
    assert uses_time_domain("Modal (undamped)", True, True) is False
    assert uses_time_domain("Time-domain (damped/nonlinear)", False, False) is True
    assert uses_time_domain("Auto", False, False) is False
    assert uses_time_domain("Auto", True, False) is True
    assert uses_time_domain("Auto", False, True) is True


@pytest.mark.parametrize("dt", [0.001, 0.01, 0.04])
def test_rk4_stability_estimate_has_no_false_alarm_for_free_free_chain(dt):
    lam = linearized_state_eigenvalues(EQUAL_MASSES, np.array([0.0, 1000.0, 1000.0, 0.0]), np.zeros(4), np.zeros(3), 0.001)

    # Fastest mode is about 54.8 rad/s, so every dt here is inside the RK4 stability region.
    assert np.max(np.abs(lam)) * dt < 2.5
    assert np.max(np.abs(rk4_amplification(lam * dt))) <= 1.0 + 1e-12
