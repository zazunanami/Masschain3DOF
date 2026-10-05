"""GUI regression tests; skipped when Tk cannot open a display (e.g. headless CI)."""

import numpy as np
import pytest

tk = pytest.importorskip("tkinter")

import masschain3dof.app as app_module  # noqa: E402
from tkinter import messagebox  # noqa: E402


class FakeClock:
    """Deterministic replacement for the time module used by the animation clock."""

    def __init__(self):
        self.now = 1000.0

    def perf_counter(self):
        return self.now


@pytest.fixture
def app(monkeypatch):
    dialogs = []
    for name in ("showinfo", "showwarning", "showerror"):
        monkeypatch.setattr(messagebox, name, lambda title, text, _kind=name: dialogs.append((_kind, title, text)))

    clock = FakeClock()
    monkeypatch.setattr(app_module, "time", clock)

    try:
        application = app_module.MassChainModalApp()
    except tk.TclError as exc:
        pytest.skip(f"Tk display not available: {exc}")
    application.withdraw()
    application.dialogs = dialogs
    application.clock = clock
    yield application
    application.destroy()


def set_entry(entry, value):
    entry.configure(state="normal")
    entry.delete(0, tk.END)
    entry.insert(0, value)


def dialog_kinds(app):
    return [kind for kind, _title, _text in app.dialogs]


def pending_animation_frames(app):
    pending = []
    for after_id in app.tk.splitlist(app.tk.call("after", "info")):
        script = app.tk.splitlist(app.tk.call("after", "info", after_id))[0]
        if "schedule_animation_step" in str(script):
            pending.append(after_id)
    return pending


def advance_animation(app, seconds):
    """Let `seconds` of wall-clock time pass and run the pending animation frame."""
    if app.animAfterId is not None:
        app.after_cancel(app.animAfterId)
        app.animAfterId = None
    app.clock.now += seconds
    app.schedule_animation_step()


def assert_shows_playback_time(app, expected_time, dt):
    """The playback clock is at expected_time and the frame is the last sample at or before it."""
    assert app.animSimTime == pytest.approx(expected_time, abs=1e-9)
    shown_time = app.lastTimeArray[app.animFrameIndex]
    assert expected_time - dt - 1e-9 <= shown_time <= expected_time + 1e-9


def test_animation_plays_in_real_time_independent_of_time_step(app):
    set_entry(app.entryTotalTime, "1.0")
    for dt, samples in ((0.001, 1001), (0.0001, 10001)):
        set_entry(app.entryTimeStep, str(dt))
        app.on_run_simulation()
        app.on_animate()
        assert app.animFrameIndex == 0

        advance_animation(app, 0.5)
        assert_shows_playback_time(app, 0.5, dt)

        app.animSpeedVar.set(2.0)
        advance_animation(app, 0.2)
        assert_shows_playback_time(app, 0.9, dt)

        advance_animation(app, 0.2)
        assert app.animFrameIndex == samples - 1
        assert app.animAfterId is None  # stops after the last sample

        app.animSpeedVar.set(1.0)


def test_restart_rewinds_and_rearms_a_finished_animation(app):
    set_entry(app.entryTotalTime, "1.0")
    app.on_run_simulation()
    app.on_animate()
    advance_animation(app, 5.0)
    assert app.animAfterId is None

    app.on_anim_restart()

    assert app.animFrameIndex == 0
    assert len(pending_animation_frames(app)) == 1


def test_second_animate_click_keeps_a_single_update_loop(app):
    set_entry(app.entryTotalTime, "0.2")
    app.on_run_simulation()
    app.on_animate()
    first_window = app.animWindow

    app.on_animate()

    assert len(pending_animation_frames(app)) == 1
    assert not first_window.winfo_exists()


def test_new_run_stops_a_running_animation(app):
    set_entry(app.entryTotalTime, "1.0")
    app.on_run_simulation()
    app.on_animate()
    advance_animation(app, 0.9)

    set_entry(app.entryTotalTime, "0.5")
    app.on_run_simulation()

    assert app.animWindow is None
    assert pending_animation_frames(app) == []


def test_any_entry_edit_closes_animation_and_disables_animate(app):
    app.on_run_simulation()
    app.on_animate()

    set_entry(app.entryM2, "2.0")  # programmatic edit, like a mouse paste: no key events

    assert str(app.buttonAnimate.cget("state")) == "disabled"
    assert app.animWindow is None


@pytest.mark.parametrize("entry_name, value", [("entryTotalTime", "nan"), ("entryTotalTime", "inf"), ("entryImpulseForce", "nan")])
def test_non_finite_inputs_are_rejected(app, entry_name, value):
    set_entry(getattr(app, entry_name), value)

    app.on_run_simulation()

    assert dialog_kinds(app) == ["showerror"]


@pytest.mark.parametrize("t_end, dt", [("10000", "0.001"), ("1e300", "1e-300")])
def test_step_limit_rejects_oversized_runs(app, t_end, dt):
    set_entry(app.entryTotalTime, t_end)
    set_entry(app.entryTimeStep, dt)

    app.on_run_simulation()

    assert dialog_kinds(app) == ["showerror"]
    assert "Too many time steps" in app.dialogs[0][2]


def test_disabled_damping_inputs_are_not_validated_in_modal_mode(app):
    app.enableViscousVar.set(True)
    set_entry(app.entryC0, "-1")
    app.simMethodVar.set("Modal (undamped)")

    app.on_run_simulation()

    assert "showerror" not in dialog_kinds(app)
    assert app.lastTimeArray is not None


def test_modal_properties_do_not_warn_about_the_time_step(app):
    set_entry(app.entryTimeStep, "0.01")

    app.on_compute_modes()

    assert dialog_kinds(app) == ["showinfo"]


def test_diverged_run_is_reported_and_not_plotted(app):
    set_entry(app.entryTotalTime, "1.0")
    app.on_run_simulation()
    app.simMethodVar.set("Time-domain (damped/nonlinear)")
    app.enableViscousVar.set(True)
    for entry in (app.entryC0, app.entryC1, app.entryC2, app.entryC3):
        set_entry(entry, "1000")
    app.dialogs.clear()

    with np.errstate(over="ignore", invalid="ignore"):
        app.on_run_simulation()

    assert dialog_kinds(app) == ["showwarning", "showerror"]
    assert "RK4" in app.dialogs[0][2]
    assert app.lastTimeArray is None
    assert str(app.buttonAnimate.cget("state")) == "disabled"
    assert len(app.axDisp.lines) == 0 and len(app.axEnergy.lines) == 0


def test_mass_boxes_stay_narrower_than_the_closest_spacing(app):
    for entry, value in zip((app.entryL0, app.entryL1, app.entryL2, app.entryL3), ("0.02", "0.03", "0.03", "0.02")):
        set_entry(entry, value)
    app.on_run_simulation()

    app.on_animate()

    assert app.massWidth < 0.02


def test_energy_axis_labels_stay_readable_for_small_energies(app):
    set_entry(app.entryImpulseForce, "0.1")  # 0.001 N*s -> 5e-7 J

    app.on_run_simulation()
    app.canvas.draw()

    low, high = app.axEnergy.get_ylim()
    labels = [label.get_text() for label in app.axEnergy.get_yticklabels() if low <= label.get_position()[1] <= high]
    assert len(labels) >= 3
    assert len(set(labels)) == len(labels)
