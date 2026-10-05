import tkinter as tk
from tkinter import ttk, messagebox
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from matplotlib.patches import Rectangle


# ------------------------------------------------------------
# 3-DOF mass-chain modal and time-domain simulation GUI.
# The module is intentionally kept as a single-file Tkinter app so it
# can be run directly or packaged by users who need a standalone build.
# ------------------------------------------------------------
DEBUG = False  # Set to True for diagnostic console output.

# Eigenvalues smaller than this fraction of the largest one are round-off (double precision ~1e-16).
EIGEN_ROUNDOFF_REL_TOL = 1e-12

# Hard limit on time steps per run: RK4 needs ~40 us per step in pure Python and every
# stored sample costs memory, so larger runs would freeze the GUI for minutes.
MAX_TIME_STEPS = 2_000_000

# Animation redraw period in milliseconds (~33 frames per second).
ANIM_FRAME_MS = 30


# ------------ MODEL & MODAL ANALYSIS ------------

def make_mass_matrix(masses):
    """Mass matrix for the 3-DOF chain."""
    m1, m2, m3 = masses  # unpacking is a tiny bit friendlier to read

    if DEBUG:
        print("make_mass_matrix() ->", masses)

  
    M = np.diag([float(m1), float(m2), float(m3)])
    return M


def make_stiffness_matrix(springs):
    """Stiffness matrix K for the 3-DOF chain."""
    k0, k1, k2, k3 = springs

    # diagonal terms (each mass sees two springs, except the ends)
    k11 = k0 + k1
    k22 = k1 + k2
    k33 = k2 + k3

    K = np.array(
        [
            [k11, -k1, 0.0],
            [-k1, k22, -k2],
            [0.0, -k2, k33],
        ],
        dtype=float,
    )
    return K


def solve_modal(M, K, negativeEigenTol=-1e-10):
    """Solve the generalized eigenproblem K @ phi = omega**2 * M @ phi."""
    m_diag = np.diag(M).astype(float)

    # This should be caught by GUI validation
    if np.any(m_diag <= 0.0):
        raise ValueError("All masses must be > 0 for the modal solve.")

    with np.errstate(divide="raise", invalid="raise"):
        Minv_sqrt = np.diag(1.0 / np.sqrt(m_diag))

    A = Minv_sqrt @ K @ Minv_sqrt
    lam, U = np.linalg.eigh(A)  # real + sorted ascending (handy)

    # Eigen-solver round-off scales with the largest eigenvalue, so "tiny" is judged
    # relative to it; a purely absolute tolerance rejects valid stiff free-free chains.
    roundOff = EIGEN_ROUNDOFF_REL_TOL * float(np.max(np.abs(lam)))

    # Numerical guard: allow tiny negative eigenvalues from round-off; reject invalid setups.
    min_lam = float(np.min(lam))
    if min_lam < min(negativeEigenTol, -roundOff):
        raise ValueError(
            "Modal solve failed: negative eigenvalue detected. "
            "Double-check masses > 0, springs >= 0, and the topology."
        )

    # Eigenvalues within round-off of zero are rigid-body modes (e.g. no wall springs):
    # set them exactly to zero so they report omega = 0 instead of round-off noise.
    lam = np.where(lam <= roundOff, 0.0, lam)
    omega = np.sqrt(lam)

    Phi = Minv_sqrt @ U  # already mass-normal-ish due to the transformation
    return omega, Phi


def mass_normalize(M, modes):
    """Mass-normalize modes so that Phi.T @ M @ Phi = I (mostly a safety pass)."""
    Phi = modes.astype(float).copy()
    n_modes = Phi.shape[1]

    # A loop is totally fine here (only 3 modes). Also easier to read/debug.
    for i in range(n_modes):
        mode_vec = Phi[:, i]
        m_norm = float(np.sqrt(mode_vec.T @ (M @ mode_vec)))
        if m_norm == 0.0:
            # Unexpected zero mass norm; leave this mode unchanged.
            continue
        Phi[:, i] = mode_vec / m_norm

    return Phi


# ------------ TIME GRID ------------

def time_grid_intervals(t_end, dt):
    """Number of uniform steps on [0, t_end] such that each step is <= dt."""
    ratio = float(t_end) / float(dt)
    n_intervals = int(round(ratio))

    # Accept ratios that are integers up to round-off (e.g. 0.7 / 0.002 = 349.99999999999994);
    # otherwise round up so the step is slightly shortened instead of overshooting dt.
    if n_intervals < 1 or abs(ratio - n_intervals) > 1e-9 * ratio:
        n_intervals = max(1, int(np.ceil(ratio)))

    return n_intervals


def make_time_grid(t_end, dt):
    """Uniform sample times on [0, t_end].

    When t_end is not an integer multiple of dt, the step is shortened to
    t_end / ceil(t_end / dt), so the grid always ends at t_end and the solver
    step equals the spacing of the returned samples.
    """
    return np.linspace(0.0, float(t_end), time_grid_intervals(t_end, dt) + 1)


# ------------ MODAL TIME RESPONSE ------------

def run_modal_sim(masses, springs, x0, v0, t_end, dt):
    """Undamped free vibration via modal superposition: x(t) = Phi @ q(t)."""
    M = make_mass_matrix(masses)
    K = make_stiffness_matrix(springs)

    omega, modes = solve_modal(M, K)
    Phi = mass_normalize(M, modes)

    t = make_time_grid(t_end, dt)
    n_steps = len(t)

    # modal initial conditions (mass-normalized projection)
    q0 = Phi.T @ (M @ x0)
    qdot0 = Phi.T @ (M @ v0)

    # Use an explicit loop for readability; only three modes are present.
    q_hist = np.zeros((n_steps, 3), dtype=float)
    qdot_hist = np.zeros((n_steps, 3), dtype=float)

    for j in range(3):
        w = float(omega[j])

        if w == 0.0:
            # zero-ish frequency: treat as a drift (doesn't really happen in the usual setup)
            q_hist[:, j] = q0[j] + qdot0[j] * t
            qdot_hist[:, j] = qdot0[j]
            continue

        cos_wt = np.cos(w * t)
        sin_wt = np.sin(w * t)

        a = float(q0[j])
        b = float(qdot0[j]) / w

        q_hist[:, j] = a * cos_wt + b * sin_wt
        qdot_hist[:, j] = -a * w * sin_wt + b * w * cos_wt

    # back to physical space
    x_hist = q_hist @ Phi.T
    v_hist = qdot_hist @ Phi.T

    # quick energy sanity check (should be constant for undamped)
    E = np.zeros(n_steps, dtype=float)
    for i in range(n_steps):
        x = x_hist[i, :]
        v = v_hist[i, :]
        T = 0.5 * v.T @ (M @ v)
        V = 0.5 * x.T @ (K @ x)
        E[i] = T + V

    return (t, x_hist, v_hist, E,
            M, K, omega, Phi, q0, qdot0)


# ------------ DAMPING & TIME-DOMAIN SIMULATION ------------

def make_damping_matrix(dashpots):
    """Viscous damping matrix C for the same chain topology as the springs."""
    c0, c1, c2, c3 = dashpots

    c11 = c0 + c1
    c22 = c1 + c2
    c33 = c2 + c3

    C = np.array(
        [
            [c11, -c1, 0.0],
            [-c1, c22, -c2],
            [0.0, -c2, c33],
        ],
        dtype=float,
    )
    return C


def force_pulse(t_now, dofIndex, F_imp, imp_dt):
    """Simple rectangular force pulse: F for 0 <= t < imp_dt on a single DOF."""
    f = np.zeros(3, dtype=float)

    # Use explicit condition structure for readability.
    if F_imp != 0.0 and imp_dt > 0.0:
        if 0.0 <= t_now < imp_dt:
            f[dofIndex] = F_imp

    return f


def friction_smoothing_eps(smoothingEps):
    """Positive, finite smoothing velocity for the tanh friction model."""
    eps = float(smoothingEps) if smoothingEps is not None else 0.0
    if not np.isfinite(eps) or eps <= 0.0:
        eps = 1e-6  # fallback so we never divide by zero (or by NaN)
    return eps


def coulomb_friction(v, fc_vec, smoothingEps):
    """Regularized Coulomb friction: F = Fc * tanh(v / eps)."""
    eps = friction_smoothing_eps(smoothingEps)

    return fc_vec * np.tanh(v / eps)


def run_time_sim(
    masses, springs, dashpots,
    x0, v0,
    imp_idx, F_imp, imp_dt,
    fc_vec, frictionEps,
    t_end, dt
):
    """Time-domain simulation with viscous + Coulomb damping using RK4.

    Solves M x'' + C x' + K x = f_ext(t) - Fc * tanh(x' / eps).
    The rectangular force pulse is held constant over each RK4 (sub)step and the
    step containing the pulse end is split there, so the applied impulse is
    exactly F_imp * imp_dt however imp_dt aligns with the time grid.
    """
    M = make_mass_matrix(masses)
    K = make_stiffness_matrix(springs)
    C = make_damping_matrix(dashpots)

    Minv = np.linalg.inv(M)

    t = make_time_grid(t_end, dt)
    n_steps = len(t)

    x_hist = np.zeros((n_steps, 3), dtype=float)
    v_hist = np.zeros((n_steps, 3), dtype=float)
    E = np.zeros(n_steps, dtype=float)

    state = np.zeros(6, dtype=float)
    state[0:3] = x0
    state[3:6] = v0

    # End of the force pulse; 0.0 means no pulse, so no step is ever split.
    pulseEnd = float(imp_dt) if (F_imp != 0.0 and imp_dt > 0.0) else 0.0

    def dstate_dt(sVec, fExt):
        x = sVec[0:3]
        v = sVec[3:6]

        if DEBUG:
            # this gets spammy fast, so leave DEBUG off unless you're chasing a bug 
            # print('x=', x, 'v=', v)
            pass

        fCoulomb = coulomb_friction(v, fc_vec, frictionEps)

        acc = Minv @ (-C @ v - K @ x - fCoulomb + fExt)

        dst = np.zeros_like(sVec)
        dst[0:3] = v
        dst[3:6] = acc
        return dst

    def rk4_step(sVec, h, fExt):
        k1 = dstate_dt(sVec, fExt)
        k2 = dstate_dt(sVec + 0.5 * h * k1, fExt)
        k3 = dstate_dt(sVec + 0.5 * h * k2, fExt)
        k4 = dstate_dt(sVec + h * k3, fExt)

        return sVec + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    for i in range(n_steps):
        x = state[0:3]
        v = state[3:6]

        x_hist[i, :] = x
        v_hist[i, :] = v

        T = 0.5 * v.T @ (M @ v)
        V = 0.5 * x.T @ (K @ x)
        E[i] = T + V

        if i == n_steps - 1:
            break

        t_i = t[i]
        t_next = t[i + 1]

        # Split the step at the pulse end so the force discontinuity lies on a sub-step boundary.
        if t_i < pulseEnd < t_next:
            subSteps = ((t_i, pulseEnd), (pulseEnd, t_next))
        else:
            subSteps = ((t_i, t_next),)

        for t_a, t_b in subSteps:
            # The pulse is constant inside a sub-step, so sampling its midpoint is unambiguous.
            fExt = force_pulse(0.5 * (t_a + t_b), imp_idx, F_imp, imp_dt)
            state = rk4_step(state, t_b - t_a, fExt)

    return t, x_hist, v_hist, E, M, K, C


# ------------ RK4 STABILITY & SOLVER SELECTION ------------

# The RK4 stability region contains the left half-disc |z| <= ~2.6; 2.5 leaves a small margin.
RK4_SAFE_RADIUS = 2.5


def rk4_amplification(z):
    """RK4 stability function R(z) = 1 + z + z^2/2 + z^3/6 + z^4/24 (|R| > 1 means growth)."""
    z = np.asarray(z, dtype=complex)
    return 1.0 + z + z**2 / 2.0 + z**3 / 6.0 + z**4 / 24.0


def linearized_state_eigenvalues(masses, springs, dashpots, fc_vec, frictionEps):
    """Eigenvalues of the first-order system matrix [[0, I], [-M^-1 K, -M^-1 C_eff]].

    Coulomb friction is linearized at v = 0, where tanh(v / eps) is stiffest,
    so it adds Fc / eps to the diagonal of the damping matrix. The chain is
    passive (K and C_eff positive semidefinite), so any positive real part is
    eigen-solver round-off (e.g. a free-free rigid-body mode) and is clipped to 0.
    """
    M = make_mass_matrix(masses)
    K = make_stiffness_matrix(springs)
    C = make_damping_matrix(dashpots)

    eps = friction_smoothing_eps(frictionEps)
    C_eff = C + np.diag(np.asarray(fc_vec, dtype=float) / eps)

    Minv = np.linalg.inv(M)
    A = np.block(
        [
            [np.zeros((3, 3)), np.eye(3)],
            [-Minv @ K, -Minv @ C_eff],
        ]
    )
    lam = np.linalg.eigvals(A)
    return np.minimum(lam.real, 0.0) + 1j * lam.imag


def uses_time_domain(simMethod, viscousActive, coulombActive):
    """Resolve the solver: Modal -> exact modal, Time-domain -> RK4, Auto -> RK4 only if damping is active."""
    simMethodLower = str(simMethod).strip().lower()
    if simMethodLower.startswith("modal"):
        return False
    if simMethodLower.startswith("time"):
        return True
    return bool(viscousActive or coulombActive)


# ------------ GUI APPLICATION ------------

class MassChainModalApp(tk.Tk):
    """GUI for 3-DOF mass-on-rail chain with modal analysis."""

    def __init__(self):
        super().__init__()

        self.title("3-DOF Mass-on-Rail Chain (Modal + Damped Simulation)")

        # --- Start maximized. Fallback: fill screen geometry. ---
        try:
            self.state("zoomed")
        except tk.TclError:
            try:
                self.attributes("-zoomed", True)
            except tk.TclError:
                screenWidth = self.winfo_screenwidth()
                screenHeight = self.winfo_screenheight()
                self.geometry(f"{screenWidth}x{screenHeight}+0+0")

        self.lastTimeArray = None
        self.lastDisplacementArray = None
        self.lastEnergyArray = None
        self.animEqPositions = None

        self.animWindow = None
        self.animCanvas = None
        self.animAxis = None
        self.animDispAxis = None
        self.animEnergyAxis = None
        self.animLine = None
        self.animDispLines = None
        self.animEnergyLine = None
        self.animFrameIndex = 0
        self.animMaxIndex = 0
        self.animAfterId = None
        # Playback clock: simulated time on screen and the wall-clock time it was last advanced.
        self.animSimTime = 0.0
        self.animLastWallTime = 0.0

        self.animSpeedVar = tk.DoubleVar(value=1.0)
        self.animXZoomVar = tk.DoubleVar(value=1.0)
        self.animXSegmentsVar = tk.DoubleVar(value=6.0)

        self.animBaseXLeft = 0.0
        self.animBaseXRight = 1.0
        self.animLeftWall = 0.0
        self.animRightWall = 1.0
        self.animMaxDisp = 0.01

        # Visual objects for animation
        self.massRects = []
        self.springLines = []
        self.massWidth = 0.1
        self.massHeight = 0.05

        self._init_styles()

        self.create_layout()
        self.apply_ui_state()

    def _init_styles(self):
        """Initialize ttk styles for enabled/disabled labels."""
        self.style = ttk.Style(self)
        # Use conservative colors; platform themes vary. Just used first theme
        self.style.configure("Disabled.TLabel", foreground="gray60")

    def set_widget_state(self, widget, isEnabled):
        """Enable/disable a ttk widget (Entry/Combobox/Button/Checkbutton)."""
        try:
            # Keep Combobox in readonly mode when enabled so the text cannot be edited.
            from tkinter import ttk as _ttk_mod
            if isinstance(widget, _ttk_mod.Combobox):
                widget.configure(state="readonly" if isEnabled else "disabled")
            else:
                widget.configure(state="normal" if isEnabled else "disabled")
        except tk.TclError:
            try:
                from tkinter import ttk as _ttk_mod
                if isinstance(widget, _ttk_mod.Combobox):
                    widget.config(state="readonly" if isEnabled else "disabled")
                else:
                    widget.config(state="normal" if isEnabled else "disabled")
            except (tk.TclError, AttributeError, ValueError) as exc:
                # Some Tk/Matplotlib ops fail on certain backends/platforms; safe to ignore.
                if DEBUG:
                    print('[debug] ignored UI/backend exception:', repr(exc))
    def set_label_enabled(self, labelWidget, isEnabled):
        """Visually grey/ungrey a label to reflect disabled state."""
        if labelWidget is None:
            return
        try:
            labelWidget.configure(style="TLabel" if isEnabled else "Disabled.TLabel")
        except (tk.TclError, AttributeError, ValueError) as exc:
            # Some Tk/Matplotlib ops fail on certain backends/platforms; safe to ignore.
            if DEBUG:
                print('[debug] ignored UI/backend exception:', repr(exc))
    def apply_ui_state(self):
        """Grey-disable inputs that are not applicable to the selected simulation method."""
        simMethod = str(self.simMethodVar.get()).strip().lower()
        isModal = simMethod.startswith("modal")
        isTimeDomain = simMethod.startswith("time")
        isAuto = simMethod.startswith("auto")

        # Damping is never used in pure modal mode
        dampingControlsEnabled = (not isModal)

        self.set_widget_state(self.comboSimMethod, True)

        # Damping checkboxes themselves
        self.set_widget_state(self.checkViscous, dampingControlsEnabled)
        self.set_widget_state(self.checkCoulomb, dampingControlsEnabled)

        viscousOn = bool(self.enableViscousVar.get()) and dampingControlsEnabled
        coulombOn = bool(self.enableCoulombVar.get()) and dampingControlsEnabled

        # Grey labels: title/checks follow dampingControlsEnabled, sub-labels follow their group enable
        self.set_label_enabled(self.labelDampingTitle, dampingControlsEnabled)

        for labelWidget in [self.labelC0, self.labelC1, self.labelC2, self.labelC3]:
            self.set_label_enabled(labelWidget, viscousOn)
        for labelWidget in [self.labelFc1, self.labelFc2, self.labelFc3, self.labelFrictionEps]:
            self.set_label_enabled(labelWidget, coulombOn)

        # Viscous damping entry controls
        for entryWidget in [self.entryC0, self.entryC1, self.entryC2, self.entryC3]:
            self.set_widget_state(entryWidget, viscousOn)

        # Coulomb damping entry controls
        for entryWidget in [self.entryFc1, self.entryFc2, self.entryFc3, self.entryFrictionEps]:
            self.set_widget_state(entryWidget, coulombOn)

        # Equivalent Δv toggle: only meaningful when time-domain may be used (Time-domain or Auto)
        self.set_widget_state(self.checkUseEquivalentDeltaV, dampingControlsEnabled)

        # Keep impulse fields editable (they are used by either force pulse or Δv conversion)
        for entryWidget in [self.entryImpulseVelocity, self.entryImpulseForce, self.entryImpulseDuration]:
            self.set_widget_state(entryWidget, True)

        # Keep compute/simulate always enabled
        self.set_widget_state(self.buttonModes, True)
        self.set_widget_state(self.buttonSimulate, True)

        # Short hint about the method
        if isAuto:
            self.labelInfo.config(text="Simulation Method = Auto: uses Modal if damping is OFF, else uses Time-domain.")
        elif isModal:
            self.labelInfo.config(text="Simulation Method = Modal: undamped free vibration (force pulse converted to Δv).")
        elif isTimeDomain:
            self.labelInfo.config(text="Simulation Method = Time-domain: RK4 with viscous + Coulomb (nonlinear) damping.")

    def create_layout(self):
        """Root: 2-column grid: left (scrollable) controls, right plots."""
        self.grid_columnconfigure(0, weight=0, minsize=380)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # --- Scrollable control panel (left) ---
        controlOuterFrame = ttk.Frame(self)
        controlOuterFrame.grid(row=0, column=0, sticky="nsew")
        controlOuterFrame.grid_rowconfigure(0, weight=1)
        controlOuterFrame.grid_columnconfigure(0, weight=1)

        controlCanvas = tk.Canvas(controlOuterFrame, borderwidth=0, highlightthickness=0)
        controlScrollbar = ttk.Scrollbar(controlOuterFrame, orient="vertical", command=controlCanvas.yview)
        controlCanvas.configure(yscrollcommand=controlScrollbar.set)

        controlCanvas.grid(row=0, column=0, sticky="nsew")
        controlScrollbar.grid(row=0, column=1, sticky="ns")

        controlInnerFrame = ttk.Frame(controlCanvas, padding=10)
        controlCanvasWindow = controlCanvas.create_window((0, 0), window=controlInnerFrame, anchor="nw")

        def on_inner_configure(_event):
            # Update scroll region to match inner frame
            controlCanvas.configure(scrollregion=controlCanvas.bbox("all"))

        def on_canvas_configure(event):
            # Make inner frame track the canvas width
            controlCanvas.itemconfigure(controlCanvasWindow, width=event.width)

        controlInnerFrame.bind("<Configure>", on_inner_configure)
        controlCanvas.bind("<Configure>", on_canvas_configure)

        # Mouse wheel scrolling only while pointer is over the control panel or you cant scroll freely
        wheelRemainder = [0.0]

        def on_mousewheel_windows(event):
            # event.delta is 120 per notch on Windows; precision touchpads send smaller deltas,
            # so accumulate fractions instead of truncating every event to zero units.
            wheelRemainder[0] -= event.delta / 120.0
            units = int(wheelRemainder[0])
            if units != 0:
                wheelRemainder[0] -= units
                controlCanvas.yview_scroll(units, "units")

        def on_mousewheel_macos(event):
            # On macOS, event.delta is smaller and has opposite sign sometimes
            controlCanvas.yview_scroll(int(-1 * event.delta), "units")

        def on_mousewheel_linux_up(_event):
            controlCanvas.yview_scroll(-3, "units")

        def on_mousewheel_linux_down(_event):
            controlCanvas.yview_scroll(3, "units")

        def bind_mousewheel(_event):
            # Bind all so the wheel works even if focus is inside an Entry
            windowSystem = self.tk.call("tk", "windowingsystem")
            if windowSystem == "aqua":
                self.bind_all("<MouseWheel>", on_mousewheel_macos)
                # macOS reports horizontal swipes as Shift-MouseWheel; the panel only scrolls vertically.
                self.bind_all("<Shift-MouseWheel>", lambda _event: None)
            else:
                self.bind_all("<MouseWheel>", on_mousewheel_windows)
                self.bind_all("<Shift-MouseWheel>", on_mousewheel_windows)

            self.bind_all("<Button-4>", on_mousewheel_linux_up)
            self.bind_all("<Button-5>", on_mousewheel_linux_down)

        def unbind_mousewheel(_event):
            self.unbind_all("<MouseWheel>")
            self.unbind_all("<Shift-MouseWheel>")
            self.unbind_all("<Button-4>")
            self.unbind_all("<Button-5>")

        controlCanvas.bind("<Enter>", bind_mousewheel)
        controlCanvas.bind("<Leave>", unbind_mousewheel)

        self.controlCanvas = controlCanvas
        self.controlFrame = controlInnerFrame

        plotFrame = ttk.Frame(self, padding=10)
        plotFrame.grid(row=0, column=1, sticky="nsew")

        self.create_controls(controlInnerFrame)
        self.create_plots(plotFrame)

    def create_controls(self, frame):
        """Left control panel, split into entriesFrame + infoFrame."""

        entriesFrame = ttk.Frame(frame)
        entriesFrame.grid(row=0, column=0, columnspan=2, sticky="nw")
        entriesFrame.grid_columnconfigure(0, weight=0, minsize=140)
        entriesFrame.grid_columnconfigure(1, weight=0, minsize=80)

        infoFrame = ttk.Frame(frame)
        infoFrame.grid(row=1, column=0, columnspan=2, sticky="we", pady=(5, 0))

        # Masses
        ttk.Label(entriesFrame, text="Masses (kg)").grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 5))

        ttk.Label(entriesFrame, text="m1:").grid(row=1, column=0, sticky="e")
        self.entryM1 = ttk.Entry(entriesFrame, width=10)
        self.entryM1.grid(row=1, column=1, sticky="w")
        self.entryM1.insert(0, "1.0")

        ttk.Label(entriesFrame, text="m2:").grid(row=2, column=0, sticky="e")
        self.entryM2 = ttk.Entry(entriesFrame, width=10)
        self.entryM2.grid(row=2, column=1, sticky="w")
        self.entryM2.insert(0, "1.0")

        ttk.Label(entriesFrame, text="m3:").grid(row=3, column=0, sticky="e")
        self.entryM3 = ttk.Entry(entriesFrame, width=10)
        self.entryM3.grid(row=3, column=1, sticky="w")
        self.entryM3.insert(0, "1.0")

        # Springs
        ttk.Label(entriesFrame, text="Springs k (N/m)").grid(row=4, column=0, columnspan=2, sticky="w", pady=(10, 5))

        ttk.Label(entriesFrame, text="k0 (wall-m1):").grid(row=5, column=0, sticky="e")
        self.entryK0 = ttk.Entry(entriesFrame, width=10)
        self.entryK0.grid(row=5, column=1, sticky="w")
        self.entryK0.insert(0, "1000.0")

        ttk.Label(entriesFrame, text="k1 (m1-m2):").grid(row=6, column=0, sticky="e")
        self.entryK1 = ttk.Entry(entriesFrame, width=10)
        self.entryK1.grid(row=6, column=1, sticky="w")
        self.entryK1.insert(0, "1000.0")

        ttk.Label(entriesFrame, text="k2 (m2-m3):").grid(row=7, column=0, sticky="e")
        self.entryK2 = ttk.Entry(entriesFrame, width=10)
        self.entryK2.grid(row=7, column=1, sticky="w")
        self.entryK2.insert(0, "1000.0")

        ttk.Label(entriesFrame, text="k3 (m3-wall):").grid(row=8, column=0, sticky="e")
        self.entryK3 = ttk.Entry(entriesFrame, width=10)
        self.entryK3.grid(row=8, column=1, sticky="w")
        self.entryK3.insert(0, "1000.0")

        # Lengths
        ttk.Label(entriesFrame, text="Spring lengths (m)").grid(row=9, column=0, columnspan=2, sticky="w", pady=(10, 5))

        ttk.Label(entriesFrame, text="L0 (wall-m1):").grid(row=10, column=0, sticky="e")
        self.entryL0 = ttk.Entry(entriesFrame, width=10)
        self.entryL0.grid(row=10, column=1, sticky="w")
        self.entryL0.insert(0, "0.5")

        ttk.Label(entriesFrame, text="L1 (m1-m2):").grid(row=11, column=0, sticky="e")
        self.entryL1 = ttk.Entry(entriesFrame, width=10)
        self.entryL1.grid(row=11, column=1, sticky="w")
        self.entryL1.insert(0, "1.0")

        ttk.Label(entriesFrame, text="L2 (m2-m3):").grid(row=12, column=0, sticky="e")
        self.entryL2 = ttk.Entry(entriesFrame, width=10)
        self.entryL2.grid(row=12, column=1, sticky="w")
        self.entryL2.insert(0, "1.0")

        ttk.Label(entriesFrame, text="L3 (m3-wall):").grid(row=13, column=0, sticky="e")
        self.entryL3 = ttk.Entry(entriesFrame, width=10)
        self.entryL3.grid(row=13, column=1, sticky="w")
        self.entryL3.insert(0, "0.5")

        # Impulse
        ttk.Label(entriesFrame, text="Impulse Settings").grid(row=14, column=0, columnspan=2, sticky="w", pady=(10, 5))

        ttk.Label(entriesFrame, text="Initial velocity (m/s):").grid(row=15, column=0, sticky="e")
        self.entryImpulseVelocity = ttk.Entry(entriesFrame, width=10)
        self.entryImpulseVelocity.grid(row=15, column=1, sticky="w")
        self.entryImpulseVelocity.insert(0, "0.0")

        ttk.Label(entriesFrame, text="Impulse force F (N):").grid(row=16, column=0, sticky="e")
        self.entryImpulseForce = ttk.Entry(entriesFrame, width=10)
        self.entryImpulseForce.grid(row=16, column=1, sticky="w")
        self.entryImpulseForce.insert(0, "1000.0")

        ttk.Label(entriesFrame, text="Impulse duration Δt (s):").grid(row=17, column=0, sticky="e")
        self.entryImpulseDuration = ttk.Entry(entriesFrame, width=10)
        self.entryImpulseDuration.grid(row=17, column=1, sticky="w")
        self.entryImpulseDuration.insert(0, "0.01")

        ttk.Label(entriesFrame, text="Impulse mass:").grid(row=18, column=0, sticky="e")
        self.impulseMassVar = tk.StringVar(value="1")
        ttk.Radiobutton(
            entriesFrame, text="Mass 1", variable=self.impulseMassVar,
            value="1", command=self.on_parameters_changed
        ).grid(row=18, column=1, sticky="w")
        ttk.Radiobutton(
            entriesFrame, text="Mass 2", variable=self.impulseMassVar,
            value="2", command=self.on_parameters_changed
        ).grid(row=19, column=1, sticky="w")
        ttk.Radiobutton(
            entriesFrame, text="Mass 3", variable=self.impulseMassVar,
            value="3", command=self.on_parameters_changed
        ).grid(row=20, column=1, sticky="w")

        # One-toggle switch for time-domain impulse interpretation
        self.useEquivalentDeltaVVar = tk.BooleanVar(value=False)
        self.checkUseEquivalentDeltaV = ttk.Checkbutton(
            entriesFrame,
            text="Time-domain: use modal-equivalent Δv (ignore force pulse)",
            variable=self.useEquivalentDeltaVVar,
            command=self.on_parameters_changed
        )
        self.checkUseEquivalentDeltaV.grid(row=21, column=0, columnspan=2, sticky="w", pady=(6, 0))

        # Simulation method
        ttk.Label(entriesFrame, text="Simulation Method").grid(row=22, column=0, columnspan=2, sticky="w", pady=(10, 5))
        self.simMethodVar = tk.StringVar(value="Auto")
        self.comboSimMethod = ttk.Combobox(
            entriesFrame,
            textvariable=self.simMethodVar,
            values=["Auto", "Modal (undamped)", "Time-domain (damped/nonlinear)"],
            state="readonly",
            width=26
        )
        self.comboSimMethod.grid(row=23, column=0, columnspan=2, sticky="we")
        self.comboSimMethod.bind("<<ComboboxSelected>>", self.on_parameters_changed)

        self.labelMethodHelp = ttk.Label(
            entriesFrame,
            text=(
                "Auto: modal if damping OFF, else time-domain.\n"
                "Modal: undamped free vibration (force pulse -> Δv).\n"
                "Time-domain: RK4 with viscous/Coulomb; default uses force pulse.\n"
                "Time-domain toggle: convert F·Δt to Δv and use it as initial velocity."
            ),
            wraplength=320,
            justify="left"
        )
        self.labelMethodHelp.grid(row=24, column=0, columnspan=2, sticky="w", pady=(5, 5))

        # Damping
        self.labelDampingTitle = ttk.Label(entriesFrame, text="Damping (Optional)")
        self.labelDampingTitle.grid(row=25, column=0, columnspan=2, sticky="w", pady=(10, 5))

        self.enableViscousVar = tk.BooleanVar(value=False)
        self.checkViscous = ttk.Checkbutton(
            entriesFrame,
            text="Enable viscous damping (dashpots c0..c3)",
            variable=self.enableViscousVar,
            command=self.on_parameters_changed
        )
        self.checkViscous.grid(row=26, column=0, columnspan=2, sticky="w")

        self.labelC0 = ttk.Label(entriesFrame, text="c0 (wall-m1) (N·s/m):")
        self.labelC0.grid(row=27, column=0, sticky="e")
        self.entryC0 = ttk.Entry(entriesFrame, width=10)
        self.entryC0.grid(row=27, column=1, sticky="w")
        self.entryC0.insert(0, "0.0")

        self.labelC1 = ttk.Label(entriesFrame, text="c1 (m1-m2) (N·s/m):")
        self.labelC1.grid(row=28, column=0, sticky="e")
        self.entryC1 = ttk.Entry(entriesFrame, width=10)
        self.entryC1.grid(row=28, column=1, sticky="w")
        self.entryC1.insert(0, "0.0")

        self.labelC2 = ttk.Label(entriesFrame, text="c2 (m2-m3) (N·s/m):")
        self.labelC2.grid(row=29, column=0, sticky="e")
        self.entryC2 = ttk.Entry(entriesFrame, width=10)
        self.entryC2.grid(row=29, column=1, sticky="w")
        self.entryC2.insert(0, "0.0")

        self.labelC3 = ttk.Label(entriesFrame, text="c3 (m3-wall) (N·s/m):")
        self.labelC3.grid(row=30, column=0, sticky="e")
        self.entryC3 = ttk.Entry(entriesFrame, width=10)
        self.entryC3.grid(row=30, column=1, sticky="w")
        self.entryC3.insert(0, "0.0")

        self.enableCoulombVar = tk.BooleanVar(value=False)
        self.checkCoulomb = ttk.Checkbutton(
            entriesFrame,
            text="Enable Coulomb damping (rail friction Fc1..Fc3)",
            variable=self.enableCoulombVar,
            command=self.on_parameters_changed
        )
        self.checkCoulomb.grid(row=31, column=0, columnspan=2, sticky="w")

        self.labelFc1 = ttk.Label(entriesFrame, text="Fc1 (N):")
        self.labelFc1.grid(row=32, column=0, sticky="e")
        self.entryFc1 = ttk.Entry(entriesFrame, width=10)
        self.entryFc1.grid(row=32, column=1, sticky="w")
        self.entryFc1.insert(0, "0.0")

        self.labelFc2 = ttk.Label(entriesFrame, text="Fc2 (N):")
        self.labelFc2.grid(row=33, column=0, sticky="e")
        self.entryFc2 = ttk.Entry(entriesFrame, width=10)
        self.entryFc2.grid(row=33, column=1, sticky="w")
        self.entryFc2.insert(0, "0.0")

        self.labelFc3 = ttk.Label(entriesFrame, text="Fc3 (N):")
        self.labelFc3.grid(row=34, column=0, sticky="e")
        self.entryFc3 = ttk.Entry(entriesFrame, width=10)
        self.entryFc3.grid(row=34, column=1, sticky="w")
        self.entryFc3.insert(0, "0.0")

        self.labelFrictionEps = ttk.Label(entriesFrame, text="Friction smoothing ε (m/s):")
        self.labelFrictionEps.grid(row=35, column=0, sticky="e")
        self.entryFrictionEps = ttk.Entry(entriesFrame, width=10)
        self.entryFrictionEps.grid(row=35, column=1, sticky="w")
        self.entryFrictionEps.insert(0, "0.001")

        # Time
        ttk.Label(entriesFrame, text="Time Settings").grid(row=36, column=0, columnspan=2, sticky="w", pady=(10, 5))

        ttk.Label(entriesFrame, text="Total time (s):").grid(row=37, column=0, sticky="e")
        self.entryTotalTime = ttk.Entry(entriesFrame, width=10)
        self.entryTotalTime.grid(row=37, column=1, sticky="w")
        self.entryTotalTime.insert(0, "5.0")

        ttk.Label(entriesFrame, text="Time step (s):").grid(row=38, column=0, sticky="e")
        self.entryTimeStep = ttk.Entry(entriesFrame, width=10)
        self.entryTimeStep.grid(row=38, column=1, sticky="w")
        self.entryTimeStep.insert(0, "0.001")

        # Buttons
        self.buttonModes = ttk.Button(entriesFrame, text="Compute Modal Properties",
                                      command=self.on_compute_modes)
        self.buttonModes.grid(row=39, column=0, columnspan=2, pady=(15, 5), sticky="we")

        self.buttonSimulate = ttk.Button(entriesFrame, text="Run Simulation",
                                         command=self.on_run_simulation)
        self.buttonSimulate.grid(row=40, column=0, columnspan=2, pady=(5, 5), sticky="we")

        self.buttonAnimate = ttk.Button(entriesFrame, text="Animate Mass Motion",
                                        command=self.on_animate,
                                        state="disabled")
        self.buttonAnimate.grid(row=41, column=0, columnspan=2, pady=(5, 5), sticky="we")

        # Info labels
        self.labelInfo = ttk.Label(infoFrame, text="", foreground="blue",
                                   wraplength=320, anchor="w", justify="left")
        self.labelInfo.pack(side=tk.TOP, anchor="w")

        self.labelSimTime = ttk.Label(infoFrame, text="", foreground="darkgreen",
                                      wraplength=320, anchor="w", justify="left")
        self.labelSimTime.pack(side=tk.TOP, anchor="w")

        self.bind_parameter_change_events(entriesFrame)

    def bind_parameter_change_events(self, _entriesFrame):
        """Any change on inputs closes animation and disables Animate button."""
        entries = [
            self.entryM1, self.entryM2, self.entryM3,
            self.entryK0, self.entryK1, self.entryK2, self.entryK3,
            self.entryL0, self.entryL1, self.entryL2, self.entryL3,
            self.entryImpulseVelocity, self.entryImpulseForce, self.entryImpulseDuration,
            self.entryC0, self.entryC1, self.entryC2, self.entryC3,
            self.entryFc1, self.entryFc2, self.entryFc3, self.entryFrictionEps,
            self.entryTotalTime, self.entryTimeStep
        ]
        # A textvariable trace sees every edit (typing, paste, cut, middle-click paste),
        # while <KeyRelease> misses mouse edits and would leave a stale animation enabled.
        self.entryVars = []
        for entryWidget in entries:
            entryVar = tk.StringVar(self, value=entryWidget.get())
            entryWidget.configure(textvariable=entryVar)
            entryVar.trace_add("write", lambda *_args: self.on_parameters_changed())
            self.entryVars.append(entryVar)

        self.comboSimMethod.bind("<<ComboboxSelected>>", self.on_parameters_changed)

        # Track checkbox / variable changes
        for variable in [
            self.enableViscousVar, self.enableCoulombVar,
            self.simMethodVar, self.useEquivalentDeltaVVar
        ]:
            try:
                variable.trace_add("write", lambda *_args: self.on_parameters_changed())
            except (tk.TclError, AttributeError, ValueError) as exc:
                # Some Tk/Matplotlib ops fail on certain backends/platforms; safe to ignore.
                if DEBUG:
                    print('[debug] ignored UI/backend exception:', repr(exc))
    def on_parameters_changed(self, event=None):
        self.buttonAnimate.config(state="disabled")
        self.apply_ui_state()
        self.stop_animation()

    def stop_animation(self):
        """Cancel the pending animation frame and close the animation window."""
        if self.animAfterId is not None:
            try:
                self.after_cancel(self.animAfterId)
            except (tk.TclError, ValueError) as exc:
                if DEBUG:
                    print('[debug] ignored UI/backend exception:', repr(exc))
            self.animAfterId = None

        if self.animWindow is not None and self.animWindow.winfo_exists():
            self.animWindow.destroy()
        self.animWindow = None

    def create_plots(self, frame):
        """Right side: matplotlib figure with displacement & energy."""
        self.figure = plt.Figure(figsize=(7, 6), dpi=100)
        self.axDisp = self.figure.add_subplot(211)
        self.axEnergy = self.figure.add_subplot(212, sharex=self.axDisp)

        self.axDisp.set_ylabel("Displacement (m)")
        self.axEnergy.set_ylabel("Mechanical Energy (J)")
        self.axEnergy.set_xlabel("Time (s)")

        # Clean tick labels for energy axis: no offset text, precision follows the tick spacing
        # (a fixed "%.4f" printed every tick of a small energy as 0.0000).
        try:
            self.axEnergy.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
            self.axEnergy.yaxis.set_major_locator(MaxNLocator(nbins=6))
        except (tk.TclError, AttributeError, ValueError) as exc:
            # Some Tk/Matplotlib ops fail on certain backends/platforms; safe to ignore.
            if DEBUG:
                print('[debug] ignored UI/backend exception:', repr(exc))
        self.canvas = FigureCanvasTkAgg(self.figure, master=frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # ------------ INPUT / SIMULATION ------------

    def _warn_time_step(self, masses, springs, t_end, dt, rk4Damping=None):
        """Warn about extreme dt choices and RK4 instability.

        rk4Damping is (dashpots, frictionForces, frictionEps) when the RK4 solver
        will run, or None when the exact modal solution is used.
        """
        warnings = []

        n_steps = time_grid_intervals(t_end, dt)
        if n_steps > 1_000_000:
            warnings.append(f"Very large step count ({n_steps:,}). This may be slow or run out of memory.")

        if dt < 1e-6:
            warnings.append("Time step is extremely small (< 1e-6 s). Simulation may be very slow.")

        try:
            M = make_mass_matrix(masses)
            K = make_stiffness_matrix(springs)
            omega, _ = solve_modal(M, K)
            wMax = float(np.max(omega))
            if wMax > 0.0:
                smallestPeriod = 2.0 * np.pi / wMax
                dtMaxRecommended = smallestPeriod / 40.0  # ~40 points per smallest period
                dtMinRecommended = smallestPeriod / 5000.0  # below this is usually overkill

                if dt > dtMaxRecommended:
                    warnings.append(
                        f"Time step may be too large for accuracy. "
                        f"Fastest period ≈ {smallestPeriod:.4g} s, recommended dt ≤ {dtMaxRecommended:.4g} s."
                    )
                if dt < dtMinRecommended:
                    warnings.append(
                        f"Time step may be unnecessarily small. "
                        f"Fastest period ≈ {smallestPeriod:.4g} s, dt < {dtMinRecommended:.4g} s is likely overkill."
                    )
        except (np.linalg.LinAlgError, ValueError):
            # If modal solve fails here, parameter validation will handle it later.
            pass

        if rk4Damping is not None:
            try:
                dashpotsUsed, frictionUsed, frictionEps = rk4Damping
                lam = linearized_state_eigenvalues(masses, springs, dashpotsUsed, frictionUsed, frictionEps)
                growth = float(np.max(np.abs(rk4_amplification(lam * dt))))
                if growth > 1.0 + 1e-9:
                    dtStable = RK4_SAFE_RADIUS / float(np.max(np.abs(lam)))
                    hint = " or increase the friction smoothing ε" if np.any(frictionUsed > 0.0) else ""
                    warnings.append(
                        f"Time step is too large for stable RK4 integration with the current stiffness/damping; "
                        f"results will be wrong or diverge. Use dt ≤ {dtStable:.3g} s{hint}."
                    )
            except (np.linalg.LinAlgError, ValueError):
                pass

        if warnings:
            messagebox.showwarning("Time Step Warning", "\n".join(warnings))

    def get_parameters_from_gui(self, checkTimeStep=True):
        """Read numeric parameters from GUI entries.

        checkTimeStep enables the step-count limit and the time-step warnings; it is
        off for the modal-properties dialog, which does not use the time settings.
        """
        try:
            m1 = float(self.entryM1.get())
            m2 = float(self.entryM2.get())
            m3 = float(self.entryM3.get())

            k0 = float(self.entryK0.get())
            k1 = float(self.entryK1.get())
            k2 = float(self.entryK2.get())
            k3 = float(self.entryK3.get())

            L0 = float(self.entryL0.get())
            L1 = float(self.entryL1.get())
            L2 = float(self.entryL2.get())
            L3 = float(self.entryL3.get())

            impulseVelocity = float(self.entryImpulseVelocity.get())
            F_imp = float(self.entryImpulseForce.get())
            imp_dt = float(self.entryImpulseDuration.get())

            simMethod = str(self.simMethodVar.get())

            enableViscous = bool(self.enableViscousVar.get())
            c0 = float(self.entryC0.get())
            c1 = float(self.entryC1.get())
            c2 = float(self.entryC2.get())
            c3 = float(self.entryC3.get())

            enableCoulomb = bool(self.enableCoulombVar.get())
            Fc1 = float(self.entryFc1.get())
            Fc2 = float(self.entryFc2.get())
            Fc3 = float(self.entryFc3.get())
            frictionEps = float(self.entryFrictionEps.get())

            t_end = float(self.entryTotalTime.get())
            dt = float(self.entryTimeStep.get())
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numeric values.")
            return None

        # --- Input validation (hard errors) ---
        masses = np.array([m1, m2, m3], dtype=float)
        springs = np.array([k0, k1, k2, k3], dtype=float)
        lengths = np.array([L0, L1, L2, L3], dtype=float)
        dashpots = np.array([c0, c1, c2, c3], dtype=float)
        frictionForces = np.array([Fc1, Fc2, Fc3], dtype=float)

        # Damping inputs only matter when the time-domain (RK4) solver can be selected.
        dampingApplies = not simMethod.strip().lower().startswith("modal")
        viscousUsed = enableViscous and dampingApplies
        coulombUsed = enableCoulomb and dampingApplies

        scalarInputs = np.array([impulseVelocity, F_imp, imp_dt, t_end, dt], dtype=float)
        finiteInputs = (
            np.all(np.isfinite(masses)) and np.all(np.isfinite(springs)) and np.all(np.isfinite(lengths))
            and np.all(np.isfinite(scalarInputs))
            and (not viscousUsed or np.all(np.isfinite(dashpots)))
            and (not coulombUsed or (np.all(np.isfinite(frictionForces)) and np.isfinite(frictionEps)))
        )
        if not finiteInputs:
            messagebox.showerror("Input Error", "Parameters must be finite numbers (no NaN/Inf).")
            return None

        if np.any(masses <= 0.0):
            messagebox.showerror("Input Error", "Masses must be > 0.")
            return None

        if np.any(springs < 0.0):
            messagebox.showerror("Input Error", "Springs must be >= 0.")
            return None

        if np.any(lengths <= 0.0):
            messagebox.showerror("Input Error", "Spring lengths must be > 0.")
            return None

        if imp_dt < 0.0:
            messagebox.showerror("Input Error", "Impulse duration must be >= 0.")
            return None

        if t_end <= 0.0 or dt <= 0.0:
            messagebox.showerror("Input Error", "Total time and time step must be positive.")
            return None

        if viscousUsed and np.any(dashpots < 0.0):
            messagebox.showerror("Input Error", "Dashpots (c0..c3) must be >= 0.")
            return None

        if coulombUsed and np.any(frictionForces < 0.0):
            messagebox.showerror("Input Error", "Coulomb friction Fc values must be >= 0.")
            return None

        if coulombUsed and frictionEps <= 0.0:
            messagebox.showerror("Input Error", "Friction smoothing ε must be > 0.")
            return None

        # Validate modal problem early so we don't silently mask bad setups
        try:
            M = make_mass_matrix(masses)
            K = make_stiffness_matrix(springs)
            _w, _phi = solve_modal(M, K)
        except (np.linalg.LinAlgError, ValueError, TypeError) as ex:
            messagebox.showerror("Model Error", str(ex))
            return None

        imp_idx = int(self.impulseMassVar.get()) - 1
        if imp_idx < 0 or imp_idx > 2:
            messagebox.showerror("Input Error", "Impulse mass selection is invalid.")
            return None

        if checkTimeStep:
            # --- Hard step limit (also catches an overflowing t_end / dt) ---
            stepRatio = t_end / dt
            if not np.isfinite(stepRatio) or time_grid_intervals(t_end, dt) > MAX_TIME_STEPS:
                messagebox.showerror(
                    "Input Error",
                    f"Too many time steps (total time / time step ≈ {stepRatio:.3g}); "
                    f"the limit is {MAX_TIME_STEPS:,}. Increase the time step or reduce the total time."
                )
                return None

            # --- Soft warnings (dt size, step count, RK4 stability) ---
            dashpotsUsed = dashpots if viscousUsed else np.zeros(4, dtype=float)
            frictionUsed = frictionForces if coulombUsed else np.zeros(3, dtype=float)
            timeDomain = uses_time_domain(simMethod, np.any(dashpotsUsed > 0.0), np.any(frictionUsed > 0.0))
            rk4Damping = (dashpotsUsed, frictionUsed, frictionEps) if timeDomain else None
            self._warn_time_step(masses, springs, t_end, dt, rk4Damping)

        return (masses, springs, lengths,
                imp_idx, impulseVelocity, F_imp, imp_dt,
                t_end, dt,
                simMethod, enableViscous, dashpots,
                enableCoulomb, frictionForces, frictionEps)

    def on_compute_modes(self):
        params = self.get_parameters_from_gui(checkTimeStep=False)
        if params is None:
            return

        masses, springs = params[0], params[1]

        M = make_mass_matrix(masses)
        K = make_stiffness_matrix(springs)

        omega, modes = solve_modal(M, K)
        Phi = mass_normalize(M, modes)

        lines = []
        lines.append("Mass matrix M:")
        lines.append(str(M))
        lines.append("")
        lines.append("Stiffness matrix K:")
        lines.append(str(K))
        lines.append("")
        lines.append("Natural frequencies (rad/s, Hz):")
        for i, w in enumerate(omega, start=1):
            f = w / (2.0 * np.pi)
            lines.append(f"  Mode {i}: ω{i} = {w:.3f} rad/s, f{i} = {f:.3f} Hz")
        lines.append("")
        lines.append("Mass-normalized mode shapes (rows x1,x2,x3; columns modes):")
        for i in range(Phi.shape[0]):
            row_values = "  ".join(f"{Phi[i, j]: .3f}" for j in range(Phi.shape[1]))
            lines.append(f"x{i+1}:  {row_values}")

        messagebox.showinfo("Modal Properties", "\n".join(lines))
        self.labelInfo.config(text="Modal properties computed (M, K, ω, Φ).")

    def on_run_simulation(self):
        params = self.get_parameters_from_gui()
        if params is None:
            return

        # A new run replaces the arrays an open animation is reading (different length -> IndexError).
        self.stop_animation()

        (masses, springs, lengths,
         imp_idx, impulseVelocity, F_imp, imp_dt,
         t_end, dt,
         simMethod, enableViscous, dashpots,
         enableCoulomb, frictionForces, frictionEps) = params

        x0 = np.zeros(3, dtype=float)
        v0 = np.zeros(3, dtype=float)
        v0[imp_idx] = impulseVelocity

        viscousActive = bool(enableViscous) and np.any(np.abs(dashpots) > 0.0)
        coulombActive = bool(enableCoulomb) and np.any(np.abs(frictionForces) > 0.0)

        useTimeDomain = uses_time_domain(simMethod, viscousActive, coulombActive)

        useEquivalentDeltaV = bool(self.useEquivalentDeltaVVar.get())

        # Run
        t0 = time.perf_counter()

        if useTimeDomain:
            dashpotsUsed = dashpots if enableViscous else np.zeros(4, dtype=float)
            frictionUsed = frictionForces if enableCoulomb else np.zeros(3, dtype=float)

            # Two impulse interpretations for time-domain:
            # - Force pulse: apply F over Δt
            # - Equivalent Δv: convert F·Δt to Δv and use it as initial velocity with no external pulse.
            effectiveVelocity = float(impulseVelocity)
            impulseForceUsed = float(F_imp)
            impulseDurationUsed = float(imp_dt)

            if useEquivalentDeltaV and abs(F_imp) > 0.0 and imp_dt > 0.0:
                massValue = float(masses[imp_idx])
                effectiveVelocity = float(impulseVelocity) + (float(F_imp) * float(imp_dt)) / massValue

                v0 = np.zeros(3, dtype=float)
                v0[imp_idx] = effectiveVelocity

                impulseForceUsed = 0.0
                impulseDurationUsed = 0.0

            (t, x_hist, v_hist, E,
             M, K, C) = run_time_sim(
                masses=masses,
                springs=springs,
                dashpots=dashpotsUsed,
                x0=x0,
                v0=v0,
                imp_idx=imp_idx,
                F_imp=impulseForceUsed,
                imp_dt=impulseDurationUsed,
                fc_vec=frictionUsed,
                frictionEps=frictionEps,
                t_end=t_end,
                dt=dt
            )

        else:
            # Modal model has no external force term, so we convert force impulse to an equivalent initial velocity.
            baseVelocity = float(impulseVelocity)
            extraVelocity = 0.0
            if abs(F_imp) > 0.0 and imp_dt > 0.0:
                massValue = float(masses[imp_idx])
                extraVelocity = (float(F_imp) * float(imp_dt)) / massValue

            effectiveVelocity = baseVelocity + extraVelocity
            v0 = np.zeros(3, dtype=float)
            v0[imp_idx] = effectiveVelocity

            (t, x_hist, v_hist, E,
             M, K, omega,
             Phi, q0, qdot0) = run_modal_sim(
                masses, springs, x0, v0, t_end, dt
            )

        t1 = time.perf_counter()
        elapsed = t1 - t0

        if not (np.all(np.isfinite(x_hist)) and np.all(np.isfinite(E))):
            # RK4 blew up (dt beyond the stability limit); never plot or animate NaN/Inf data.
            self.lastTimeArray = None
            self.lastDisplacementArray = None
            self.lastEnergyArray = None
            self.buttonAnimate.config(state="disabled")
            self.labelInfo.config(text="Simulation diverged (non-finite values). Reduce the time step.")
            self.labelSimTime.config(text="")

            # Clear the previous run's curves so they are not mistaken for this run's result.
            self.axDisp.cla()
            self.axEnergy.cla()
            self.axDisp.set_ylabel("Displacement (m)")
            self.axEnergy.set_ylabel("Mechanical Energy (J)")
            self.axEnergy.set_xlabel("Time (s)")
            self.canvas.draw()

            messagebox.showerror(
                "Simulation Diverged",
                "The time-domain integration produced non-finite values.\n"
                "Reduce the time step (or increase the friction smoothing ε) and run again."
            )
            return

        # Geometry for animation
        L0, L1, L2, L3 = lengths
        eq1 = L0
        eq2 = L0 + L1
        eq3 = L0 + L1 + L2
        self.animEqPositions = np.array([eq1, eq2, eq3], dtype=float)

        self.lastTimeArray = t
        self.lastDisplacementArray = x_hist
        self.lastEnergyArray = E

        # X-domain and walls
        self.animLeftWall = 0.0
        self.animRightWall = L0 + L1 + L2 + L3

        maxDisp = float(np.max(np.abs(self.lastDisplacementArray)))
        if maxDisp <= 0:
            maxDisp = 0.01
        self.animMaxDisp = maxDisp

        margin = 0.05 * (self.animRightWall - self.animLeftWall) + maxDisp
        self.animBaseXLeft = self.animLeftWall - margin
        self.animBaseXRight = self.animRightWall + margin

        # Enable animation
        self.buttonAnimate.config(state="normal")

        # Plots
        self.axDisp.cla()
        self.axEnergy.cla()

        self.axDisp.plot(t, x_hist[:, 0], label="x1")
        self.axDisp.plot(t, x_hist[:, 1], label="x2")
        self.axDisp.plot(t, x_hist[:, 2], label="x3")
        self.axDisp.set_ylabel("Displacement (m)")
        self.axDisp.legend(loc="upper right")
        self.axDisp.grid(True)
        self.axEnergy.plot(t, E)
        self.axEnergy.set_ylabel("Mechanical Energy (J)")
        self.axEnergy.set_xlabel("Time (s)")

        # Auto-scale, but prevent the y-range from becoming absurdly tiny for near-constant energy.
        try:
            self.axEnergy.relim()
            self.axEnergy.autoscale_view(scalex=True, scaley=True)

            yMin, yMax = self.axEnergy.get_ylim()
            yRange = float(yMax - yMin)
            energyMax = float(np.max(E)) if len(E) > 0 else 0.0

            # Minimum visible range: 2% of max energy, so round-off ripple of a conserved energy is
            # not magnified. No absolute floor: a fixed 1e-4 J floor flattened small-energy curves.
            minRange = 0.02 * abs(energyMax) if energyMax != 0.0 else 1.0
            if yRange < minRange:
                yMid = 0.5 * (yMin + yMax)
                self.axEnergy.set_ylim(yMid - 0.5 * minRange, yMid + 0.5 * minRange)
        except (tk.TclError, AttributeError, ValueError) as exc:
            # Some Tk/Matplotlib ops fail on certain backends/platforms; safe to ignore.
            if DEBUG:
                print('[debug] ignored UI/backend exception:', repr(exc))
        # Show clean tick labels: no offset text, precision follows the tick spacing.
        try:
            self.axEnergy.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
            self.axEnergy.yaxis.set_major_locator(MaxNLocator(nbins=6))
        except (tk.TclError, AttributeError, ValueError) as exc:
            # Some Tk/Matplotlib ops fail on certain backends/platforms; safe to ignore.
            if DEBUG:
                print('[debug] ignored UI/backend exception:', repr(exc))
        self.axEnergy.grid(True)

        self.canvas.draw()

        # Name the solver, not the damping: Time-domain mode also runs with damping switched off.
        methodText = "Time-domain (RK4)" if useTimeDomain else "Modal (undamped)"
        if useTimeDomain and useEquivalentDeltaV:
            impulseText = "Impulse: modal-equivalent Δv (no force pulse applied)."
        elif useTimeDomain:
            impulseText = "Impulse: rectangular force pulse applied."
        else:
            impulseText = "Impulse: force pulse converted to Δv for modal model."

        self.labelInfo.config(
            text=f"{methodText}. {impulseText}\nExcited mass: {imp_idx+1}. Initial velocity used: {effectiveVelocity:.4f} m/s"
        )
        self.labelSimTime.config(
            text=f"Simulated 0-{t_end:.4g} s in {elapsed:.3f} s  (steps: {len(t) - 1})"
        )

        # Pop-up details
        lines = []
        lines.append(f"Method: {methodText}")
        lines.append(f"Impulse mass: {imp_idx+1}")
        lines.append(f"Initial velocity used on excited mass: {effectiveVelocity:.4f} m/s")

        if useTimeDomain:
            if useEquivalentDeltaV:
                lines.append("Impulse interpretation: Equivalent Δv = (F·Δt)/m (no external pulse).")
                lines.append(f"Computed from: F={F_imp:.3g} N, Δt={imp_dt:.3g} s")
            else:
                lines.append(f"Force pulse: F={F_imp:.3g} N, Δt={imp_dt:.3g} s")

            if viscousActive:
                lines.append(f"Viscous damping (c0..c3) = {dashpots}")
            else:
                lines.append("Viscous damping: OFF")
            if coulombActive:
                lines.append(f"Coulomb friction (Fc1..Fc3) = {frictionForces}, ε={frictionEps:.4g} m/s")
            else:
                lines.append("Coulomb friction: OFF")

            messagebox.showinfo("Simulation Setup (Time-Domain)", "\n".join(lines))
        else:
            lines.append("\nModal initial conditions:")
            lines.append(f"q(0)    = {q0}")
            lines.append(f"qdot(0) = {qdot0}")
            lines.append("\nNote: Modal model treats force impulse as initial momentum: Δv = (F·Δt)/m.")
            messagebox.showinfo("Simulation Setup (Modal)", "\n".join(lines))

    # ------------ ANIMATION HELPERS ------------

    
    def on_xzoom_changed(self, value):
        try:
            zoom = float(value)
        except (TypeError, ValueError):
            try:
                zoom = float(self.animXZoomVar.get())
            except (TypeError, ValueError, tk.TclError):
                zoom = 1.0

        if zoom < 0.5:
            zoom = 0.5
        if zoom > 3.0:
            zoom = 3.0

        self.animXZoomVar.set(zoom)

        if getattr(self, "animXZoomValueLabel", None) is not None:
            try:
                self.animXZoomValueLabel.config(text=f"{zoom:.2f}x")
            except (tk.TclError, AttributeError, ValueError) as exc:
                # Some Tk/Matplotlib ops fail on certain platforms; safe to ignore.
                if DEBUG:
                    print('[debug] ignored UI/backend exception:', repr(exc))
        # Re-apply axis limits with the new zoom
        try:
            self.update_anim_axes_limits()
        except (tk.TclError, AttributeError, ValueError) as exc:
            # Some Tk/Matplotlib ops fail on certain platforms; safe to ignore.
            if DEBUG:
                print('[debug] ignored UI/backend exception:', repr(exc))
    def update_anim_axes_limits(self, *args):
        """Update x/y limits in animation (do NOT change ticks here)."""
        if self.animAxis is None:
            return

        zoom = self.animXZoomVar.get()
        if zoom <= 0.1:
            zoom = 0.1

        xCenter = (self.animBaseXLeft + self.animBaseXRight) / 2.0
        halfWidth = (self.animBaseXRight - self.animBaseXLeft) / 2.0 / zoom
        self.animAxis.set_xlim(xCenter - halfWidth, xCenter + halfWidth)

        # Vertical limit: based on mass height and displacement
        maxLevel = max(self.animMaxDisp, self.massHeight)
        yMargin = 0.2 * maxLevel
        self.animAxis.set_ylim(-maxLevel - yMargin, maxLevel + yMargin)

    
    def on_speed_changed(self, value):
        try:
            speed = float(value)
        except (TypeError, ValueError):
            try:
                speed = float(self.animSpeedVar.get())
            except (TypeError, ValueError, tk.TclError):
                speed = 1.0

        if speed < 0.2:
            speed = 0.2
        if speed > 5.0:
            speed = 5.0

        # Update variable and label
        self.animSpeedVar.set(speed)
        if getattr(self, "animSpeedValueLabel", None) is not None:
            try:
                self.animSpeedValueLabel.config(text=f"{speed:.2f}x")
            except (tk.TclError, AttributeError, ValueError) as exc:
                # Some Tk/Matplotlib ops fail on certain backends/platforms; safe to ignore.
                if DEBUG:
                    print('[debug] ignored UI/backend exception:', repr(exc))
    def on_tick_segments_changed(self, *args):
        """Change only the number of divisions between leftWall and rightWall."""
        if self.animWindow is None or not self.animWindow.winfo_exists():
            return
        if self.animAxis is None or self.animCanvas is None:
            return

        value = self.animXSegmentsVar.get()
        nSegments = int(round(value))
        if nSegments < 1:
            nSegments = 1

        # Update numeric label if present
        if getattr(self, "animSegmentsValueLabel", None) is not None:
            try:
                self.animSegmentsValueLabel.config(text=str(nSegments))
            except (tk.TclError, AttributeError, ValueError) as exc:
                # Some Tk/Matplotlib ops fail on certain backends/platforms; safe to ignore.
                if DEBUG:
                    print('[debug] ignored UI/backend exception:', repr(exc))
        if self.animRightWall > self.animLeftWall:
            ticks = np.linspace(self.animLeftWall,
                                self.animRightWall,
                                nSegments + 1)
            self.animAxis.set_xticks(ticks)
            try:
                self.animCanvas.draw_idle()
            except (tk.TclError, AttributeError, ValueError) as exc:
                # Some Tk/Matplotlib ops fail on certain backends/platforms; safe to ignore.
                if DEBUG:
                    print('[debug] ignored UI/backend exception:', repr(exc))
    def build_spring_shape(self, x_start, x_end, n_coils=6, y_center=0.0, amplitude=None):
        """Return x,y points for a simple zigzag spring between two x-positions."""
        if amplitude is None:
            amplitude = 0.4 * self.massHeight

        length = x_end - x_start
        if length <= 0:
            # Degenerate: just a straight line
            return np.array([x_start, x_end]), np.array([y_center, y_center])

        # 2*n_coils + 1 points so start and end at center line
        n_points = 2 * n_coils + 1
        x_vals = np.linspace(x_start, x_end, n_points)

        y_vals = np.zeros(n_points)
        # pattern: 0, +a, -a, +a, -a, ..., 0
        for i in range(1, n_points - 1):
            if i % 2 == 1:
                y_vals[i] = y_center + amplitude
            else:
                y_vals[i] = y_center - amplitude

        return x_vals, y_vals

    # ------------ ANIMATION MAIN ------------

    def on_animate(self):
        if self.lastTimeArray is None or self.lastDisplacementArray is None:
            messagebox.showwarning("No Data", "Run a simulation before starting the animation.")
            return
        if self.lastEnergyArray is None:
            messagebox.showwarning("No Data", "Run a simulation before starting the animation.")
            return
        if self.animEqPositions is None:
            messagebox.showwarning("No Geometry", "Equilibrium positions are not set.")
            return

        # Close a previous animation window and stop its update loop, so only one loop ever runs.
        self.stop_animation()

        self.animWindow = tk.Toplevel(self)
        self.animWindow.title("Mass Motion Animation")

        # Make the animation window visible and prefer maximized geometry for backend stability.
        try:
            self.animWindow.state("zoomed")
        except tk.TclError:
            try:
                self.animWindow.attributes("-zoomed", True)
            except tk.TclError:
                screenWidth = self.winfo_screenwidth()
                screenHeight = self.winfo_screenheight()
                self.animWindow.geometry(f"{int(screenWidth*0.95)}x{int(screenHeight*0.90)}+20+20")

        self.animWindow.lift()
        try:
            self.animWindow.focus_force()
        except (tk.TclError, AttributeError, ValueError) as exc:
        
            if DEBUG:
                print('[debug] ignored UI/backend exception:', repr(exc))
        fig = plt.Figure(figsize=(10, 4.5), dpi=100)
        gridSpec = fig.add_gridspec(
            2, 2,
            width_ratios=[2.2, 1.3],
            height_ratios=[1.0, 1.0],
            wspace=0.35,
            hspace=0.35
        )

        self.animAxis = fig.add_subplot(gridSpec[:, 0])
        self.animDispAxis = fig.add_subplot(gridSpec[0, 1])
        self.animEnergyAxis = fig.add_subplot(gridSpec[1, 1], sharex=self.animDispAxis)

        # Box size follows the chain length but stays narrower than the closest spacing, so the
        # masses never look overlapped at equilibrium (a 1 m floor made short chains overlap).
        totalLength = self.animRightWall - self.animLeftWall
        nodePositions = np.concatenate(([self.animLeftWall], self.animEqPositions, [self.animRightWall]))
        minSpacing = float(np.min(np.diff(nodePositions)))
        self.massWidth = min(0.08 * totalLength, 0.5 * minSpacing)
        self.massHeight = 0.05 * totalLength

        # Walls
        self.animAxis.axvline(x=self.animLeftWall, linestyle="--")
        self.animAxis.axvline(x=self.animRightWall, linestyle="--")

        # Mass rectangles (3 pieces)
        self.massRects = []
        for _ in range(3):
            rect = Rectangle((0.0, -self.massHeight / 2.0),
                             self.massWidth, self.massHeight,
                             fill=False)
            self.animAxis.add_patch(rect)
            self.massRects.append(rect)

        # Spring lines (4 pieces)
        self.springLines = []
        for _ in range(4):
            line, = self.animAxis.plot([], [], "-", linewidth=1.5)
            self.springLines.append(line)

        # Adjust vertical scale so boxes are visible
        self.animMaxDisp = max(self.animMaxDisp, self.massHeight)

        # Initial zoom and segments
        self.animXZoomVar.set(1.0)
        self.animXSegmentsVar.set(6.0)

        self.update_anim_axes_limits()
        self.on_tick_segments_changed()

        self.animAxis.set_xlabel("Position along rail")
        self.animAxis.set_yticks([])
        self.animAxis.grid(True)

        # Runtime plots (right side)
        self.animDispAxis.set_ylabel("Displacement (m)")
        self.animEnergyAxis.set_ylabel("Mechanical Energy (J)")
        self.animEnergyAxis.set_xlabel("Time (s)")

        # Clean tick labels for energy axis.
        try:
            self.animEnergyAxis.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
            self.animEnergyAxis.yaxis.set_major_locator(MaxNLocator(nbins=6))
        except (tk.TclError, AttributeError, ValueError) as exc:
            # Some Tk/Matplotlib ops fail on certain backends/platforms; safe to ignore.
            if DEBUG:
                print('[debug] ignored UI/backend exception:', repr(exc))
        self.animDispAxis.grid(True)
        self.animEnergyAxis.grid(True)

        self.animDispLines = []
        for labelText in ["x1", "x2", "x3"]:
            line, = self.animDispAxis.plot([], [], label=labelText)
            self.animDispLines.append(line)
        self.animDispAxis.legend(loc="upper right")

        self.animEnergyLine, = self.animEnergyAxis.plot([], [])

        tMin = float(self.lastTimeArray[0])
        tMax = float(self.lastTimeArray[-1])
        if tMax <= tMin:
            tMax = tMin + 1.0
        self.animDispAxis.set_xlim(tMin, tMax)

        dispMax = float(np.max(np.abs(self.lastDisplacementArray)))
        if dispMax <= 0.0:
            dispMax = 1e-6
        self.animDispAxis.set_ylim(-1.1 * dispMax, 1.1 * dispMax)

        energyMax = float(np.max(self.lastEnergyArray))
        if energyMax <= 0.0:
            energyMax = 1.0
        self.animEnergyAxis.set_ylim(0.0, 1.1 * energyMax)

        self.animCanvas = FigureCanvasTkAgg(fig, master=self.animWindow)
        self.animCanvas.draw()
        self.animCanvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Controls: speed + x zoom + segment slider + restart
        ctrlFrame = ttk.Frame(self.animWindow, padding=5)
        ctrlFrame.pack(side=tk.BOTTOM, fill=tk.X)

        # Speed row
        speedFrame = ttk.Frame(ctrlFrame)
        speedFrame.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(speedFrame, text="Animation speed (0.2x - 5x real time)").pack(side=tk.LEFT)
        speedScale = ttk.Scale(
            speedFrame,
            from_=0.2,
            to=5.0,
            orient="horizontal",
            variable=self.animSpeedVar,
            command=self.on_speed_changed
        )
        speedScale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(10, 10))
        speedScale.set(self.animSpeedVar.get())

        # Current speed label (e.g. "1.00x")
        self.animSpeedValueLabel = ttk.Label(speedFrame, text=f"{self.animSpeedVar.get():.2f}x")
        self.animSpeedValueLabel.pack(side=tk.RIGHT, padx=(4, 0))

        restartBtn = ttk.Button(speedFrame, text="Restart animation",
                                command=self.on_anim_restart)
        restartBtn.pack(side=tk.RIGHT)

        # X-zoom row
        zoomFrame = ttk.Frame(ctrlFrame)
        zoomFrame.pack(side=tk.TOP, fill=tk.X, pady=(4, 0))
        ttk.Label(zoomFrame, text="X zoom (0.5x - 3x)").pack(side=tk.LEFT)
        xZoomScale = ttk.Scale(
            zoomFrame,
            from_=0.5,
            to=3.0,
            orient="horizontal",
            variable=self.animXZoomVar,
            command=self.on_xzoom_changed
        )
        xZoomScale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(10, 10))
        xZoomScale.set(self.animXZoomVar.get())

        # Current X-zoom label, e.g. "1.00x"
        self.animXZoomValueLabel = ttk.Label(zoomFrame, text=f"{self.animXZoomVar.get():.2f}x")
        self.animXZoomValueLabel.pack(side=tk.RIGHT, padx=(4, 0))

        # X segments row (only xticks, limits stay fixed)
        segFrame = ttk.Frame(ctrlFrame)
        segFrame.pack(side=tk.TOP, fill=tk.X, pady=(4, 0))
        ttk.Label(segFrame, text="X axis segments (pieces)").pack(side=tk.LEFT)
        segScale = ttk.Scale(
            segFrame,
            from_=4.0,
            to=20.0,
            orient="horizontal",
            variable=self.animXSegmentsVar,
            command=self.on_tick_segments_changed
        )
        segScale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(10, 10))
        segScale.set(self.animXSegmentsVar.get())

        # Current segment count label
        self.animSegmentsValueLabel = ttk.Label(segFrame, text=str(int(round(self.animXSegmentsVar.get()))))
        self.animSegmentsValueLabel.pack(side=tk.RIGHT, padx=(4, 0))

        self.animFrameIndex = 0
        self.animMaxIndex = len(self.lastTimeArray)

        # Start the playback clock only now, after the (slow) window and figure setup.
        self.animSimTime = float(self.lastTimeArray[0])
        self.animLastWallTime = time.perf_counter()
        self.schedule_animation_step()

    def on_anim_restart(self):
        if self.lastTimeArray is None or self.animWindow is None or not self.animWindow.winfo_exists():
            return

        # Rewind the playback clock; re-arm the update loop if it already stopped at the last frame.
        self.animSimTime = float(self.lastTimeArray[0])
        self.animLastWallTime = time.perf_counter()
        if self.animAfterId is None:
            self.schedule_animation_step()

    def schedule_animation_step(self):
        """Draw the frame for the current playback time and schedule the next one.

        Playback follows the wall clock: at 1x one second of simulated time takes one
        second, independent of dt and of how long a redraw takes (slow redraws skip samples).
        """
        self.animAfterId = None
        if self.animWindow is None or not self.animWindow.winfo_exists():
            return

        speed = self.animSpeedVar.get()
        if speed < 0.2:
            speed = 0.2
        if speed > 5.0:
            speed = 5.0

        # Advance the playback clock by the real time elapsed since the previous frame.
        now = time.perf_counter()
        self.animSimTime += (now - self.animLastWallTime) * speed
        self.animLastWallTime = now

        # Show the last sample at or before the playback time.
        frameIndex = int(np.searchsorted(self.lastTimeArray, self.animSimTime, side="right")) - 1
        self.animFrameIndex = min(max(frameIndex, 0), self.animMaxIndex - 1)

        disp = self.lastDisplacementArray[self.animFrameIndex, :]
        x_pos = self.animEqPositions + disp  # 3 masses

        # Update mass rectangles
        for i, rect in enumerate(self.massRects):
            cx = x_pos[i]
            rect.set_xy((cx - self.massWidth / 2.0, -self.massHeight / 2.0))

        # Update springs: 4 segments (wall-m1, m1-m2, m2-m3, m3-wall)
        segment_pairs = [
            (self.animLeftWall, x_pos[0]),
            (x_pos[0], x_pos[1]),
            (x_pos[1], x_pos[2]),
            (x_pos[2], self.animRightWall)
        ]
        for i, (xa, xb) in enumerate(segment_pairs):
            xs, ys = self.build_spring_shape(xa, xb)
            self.springLines[i].set_data(xs, ys)

        # Update runtime plots (show history up to current frame)
        if self.animDispLines is not None and self.animEnergyLine is not None:
            idx = self.animFrameIndex + 1
            tSlice = self.lastTimeArray[:idx]
            self.animDispLines[0].set_data(tSlice, self.lastDisplacementArray[:idx, 0])
            self.animDispLines[1].set_data(tSlice, self.lastDisplacementArray[:idx, 1])
            self.animDispLines[2].set_data(tSlice, self.lastDisplacementArray[:idx, 2])
            self.animEnergyLine.set_data(tSlice, self.lastEnergyArray[:idx])

        self.update_anim_axes_limits()
        try:
            self.animCanvas.draw_idle()
        except (tk.TclError, RuntimeError) as exc:
            # Avoid a hard crash if the backend fails during aggressive window resizing.
            if DEBUG:
                print('[debug] draw_idle failed:', repr(exc))

        if self.animFrameIndex >= self.animMaxIndex - 1:
            return  # last sample shown; "Restart animation" re-arms the loop

        self.animAfterId = self.after(ANIM_FRAME_MS, self.schedule_animation_step)


if __name__ == "__main__":
    app = MassChainModalApp()
    app.mainloop()