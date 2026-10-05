import numpy as np
import pytest

from masschain3dof.app import (
    make_mass_matrix,
    make_stiffness_matrix,
    mass_normalize,
    run_modal_sim,
    run_time_sim,
    solve_modal,
)


def test_mass_and_stiffness_matrices_for_equal_chain():
    masses = np.array([1.0, 1.0, 1.0])
    springs = np.array([1000.0, 1000.0, 1000.0, 1000.0])

    mass_matrix = make_mass_matrix(masses)
    stiffness_matrix = make_stiffness_matrix(springs)

    assert np.allclose(mass_matrix, np.diag([1.0, 1.0, 1.0]))
    assert np.allclose(
        stiffness_matrix,
        np.array(
            [
                [2000.0, -1000.0, 0.0],
                [-1000.0, 2000.0, -1000.0],
                [0.0, -1000.0, 2000.0],
            ]
        ),
    )


def test_modal_solution_is_mass_normalized_and_energy_conserving():
    masses = np.array([1.0, 1.0, 1.0])
    springs = np.array([1000.0, 1000.0, 1000.0, 1000.0])
    x0 = np.zeros(3)
    v0 = np.array([1.0, 0.0, 0.0])

    mass_matrix = make_mass_matrix(masses)
    stiffness_matrix = make_stiffness_matrix(springs)
    omega, modes = solve_modal(mass_matrix, stiffness_matrix)
    modal_matrix = mass_normalize(mass_matrix, modes)

    expected_omega = np.array([24.20302538, 44.72135955, 58.43127213])
    assert np.allclose(omega, expected_omega, rtol=1e-7, atol=1e-7)
    assert np.allclose(modal_matrix.T @ mass_matrix @ modal_matrix, np.eye(3), atol=1e-12)

    _, _, _, energy, *_ = run_modal_sim(masses, springs, x0, v0, 1.0, 0.001)
    relative_drift = (energy.max() - energy.min()) / max(abs(energy[0]), 1e-12)
    assert relative_drift < 1e-10


def test_time_domain_viscous_damping_reduces_mechanical_energy():
    masses = np.array([1.0, 1.0, 1.0])
    springs = np.array([1000.0, 1000.0, 1000.0, 1000.0])
    dashpots = np.array([1.0, 1.0, 1.0, 1.0])
    x0 = np.zeros(3)
    v0 = np.array([1.0, 0.0, 0.0])

    _, _, _, energy, *_ = run_time_sim(
        masses=masses,
        springs=springs,
        dashpots=dashpots,
        x0=x0,
        v0=v0,
        imp_idx=0,
        F_imp=0.0,
        imp_dt=0.0,
        fc_vec=np.zeros(3),
        frictionEps=0.001,
        t_end=1.0,
        dt=0.001,
    )

    assert energy[-1] < energy[0]


def test_stiff_free_free_chain_is_accepted_and_rigid_body_mode_is_exactly_zero():
    # Without wall springs the chain has a rigid-body mode; at this stiffness the eigen-solver
    # round-off (~1e-4) used to trip the absolute negative-eigenvalue tolerance.
    mass_matrix = make_mass_matrix(np.array([1e-3, 1e-3, 1e-3]))
    stiffness_matrix = make_stiffness_matrix(np.array([0.0, 1e9, 1e9, 0.0]))

    omega, _ = solve_modal(mass_matrix, stiffness_matrix)

    assert omega[0] == 0.0
    assert np.allclose(omega[1:], [1e6, np.sqrt(3.0) * 1e6], rtol=1e-9)


def test_soft_free_free_chain_reports_zero_frequency_instead_of_round_off():
    omega, _ = solve_modal(make_mass_matrix(np.ones(3)), make_stiffness_matrix(np.array([0.0, 1000.0, 1000.0, 0.0])))

    assert omega[0] == 0.0
    assert np.allclose(omega[1:], [np.sqrt(1000.0), np.sqrt(3000.0)], rtol=1e-12)


def test_negative_stiffness_is_still_rejected():
    with pytest.raises(ValueError):
        solve_modal(make_mass_matrix(np.ones(3)), make_stiffness_matrix(np.array([-1500.0, 1000.0, 1000.0, 1000.0])))


def test_free_free_modal_response_moves_centre_of_mass_at_constant_velocity():
    masses = np.array([1.0, 2.0, 3.0])
    springs = np.array([0.0, 1000.0, 500.0, 0.0])
    v0 = np.array([1.0, 0.0, 0.0])

    t, x_hist, _, energy, *_ = run_modal_sim(masses, springs, np.zeros(3), v0, 2.0, 0.001)

    centre_of_mass = x_hist @ masses / masses.sum()
    assert np.allclose(centre_of_mass, t * (masses @ v0) / masses.sum(), atol=1e-12)
    assert np.allclose(energy, energy[0], rtol=1e-10)
