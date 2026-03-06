import numpy as np
from numpy import sin, cos
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# -----------------------------
# Triple Pendulum (planar) simulation + animation
# Generalized coordinates: q = [theta1, theta2, theta3]
# State vector y = [theta1, theta2, theta3, omega1, omega2, omega3]
# Uses Lagrange equations with:
#   M(q) * qdd + C(q, qd) + G(q) = 0
# where:
#   M(q) = sum_k m_k * J_k^T J_k
#   G(q) = dV/dq  (with V = sum_k m_k g y_k)
#   C computed via Christoffel symbols from M(q)
# -----------------------------

def triple_pendulum_rhs(t, y, params):
    """
    Returns dy/dt for triple pendulum.
    """
    th1, th2, th3, w1, w2, w3 = y
    q = np.array([th1, th2, th3], dtype=float)
    qd = np.array([w1, w2, w3], dtype=float)

    l1, l2, l3 = params["l1"], params["l2"], params["l3"]
    m1, m2, m3 = params["m1"], params["m2"], params["m3"]
    g = params["g"]
    damp = params["damp"]

    # Precompute sines/cosines for speed
    s1, c1 = sin(th1), cos(th1)
    s2, c2 = sin(th2), cos(th2)
    s3, c3 = sin(th3), cos(th3)

    # --- Jacobians of bob positions wrt q ---
    # Position:
    # x1 = l1 sin(th1);            y1 = -l1 cos(th1)
    # x2 = x1 + l2 sin(th2);       y2 = y1 - l2 cos(th2)
    # x3 = x2 + l3 sin(th3);       y3 = y2 - l3 cos(th3)

    # J1 = d[x1, y1]/d[th1, th2, th3]
    J1 = np.array([
        [l1 * c1,        0.0,       0.0],
        [l1 * s1,        0.0,       0.0]
    ], dtype=float)

    # J2 = d[x2, y2]/d[th1, th2, th3]
    J2 = np.array([
        [l1 * c1,  l2 * c2,  0.0],
        [l1 * s1,  l2 * s2,  0.0]
    ], dtype=float)

    # J3 = d[x3, y3]/d[th1, th2, th3]
    J3 = np.array([
        [l1 * c1,  l2 * c2,  l3 * c3],
        [l1 * s1,  l2 * s2,  l3 * s3]
    ], dtype=float)

    # --- Mass matrix M(q) = sum_k m_k * Jk^T Jk ---
    M = (m1 * (J1.T @ J1) +
         m2 * (J2.T @ J2) +
         m3 * (J3.T @ J3))

    # --- Potential energy gradient G(q) = dV/dq ---
    # y1 = -l1 cos th1
    # y2 = -l1 cos th1 - l2 cos th2
    # y3 = -l1 cos th1 - l2 cos th2 - l3 cos th3
    # V = sum m_k g y_k
    # dV/dth1 = (m1+m2+m3)*g*l1*sin(th1)
    # dV/dth2 = (m2+m3)*g*l2*sin(th2)
    # dV/dth3 = (m3)*g*l3*sin(th3)
    G = np.array([
        (m1 + m2 + m3) * g * l1 * s1,
        (m2 + m3)      * g * l2 * s2,
        (m3)           * g * l3 * s3
    ], dtype=float)

    # --- Compute Coriolis/centrifugal term C(q,qd) using Christoffel symbols ---
    # C_i = sum_{j,k} Gamma_{i,jk} qd_j qd_k
    # Gamma_{i,jk} = 0.5*(dM_{ij}/dq_k + dM_{ik}/dq_j - dM_{jk}/dq_i)

    # Numerical partial derivatives of M wrt each angle (central difference)
    eps = 1e-7
    dM = np.zeros((3, 3, 3), dtype=float)  # dM[:,:,k] = dM/dq_k
    for k in range(3):
        dq = np.zeros(3)
        dq[k] = eps
        qp = q + dq
        qm = q - dq

        def mass_matrix(q_local):
            t1, t2, t3 = q_local
            s1l, c1l = sin(t1), cos(t1)
            s2l, c2l = sin(t2), cos(t2)
            s3l, c3l = sin(t3), cos(t3)

            J1l = np.array([[l1*c1l, 0.0, 0.0],
                            [l1*s1l, 0.0, 0.0]])
            J2l = np.array([[l1*c1l, l2*c2l, 0.0],
                            [l1*s1l, l2*s2l, 0.0]])
            J3l = np.array([[l1*c1l, l2*c2l, l3*c3l],
                            [l1*s1l, l2*s2l, l3*s3l]])

            return (m1*(J1l.T@J1l) + m2*(J2l.T@J2l) + m3*(J3l.T@J3l))

        Mp = mass_matrix(qp)
        Mm = mass_matrix(qm)
        dM[:, :, k] = (Mp - Mm) / (2.0 * eps)

    # Build C vector
    C = np.zeros(3, dtype=float)
    for i in range(3):
        s = 0.0
        for j in range(3):
            for k in range(3):
                Gamma = 0.5 * (dM[i, j, k] + dM[i, k, j] - dM[j, k, i])
                s += Gamma * qd[j] * qd[k]
        C[i] = s

    # Damping torque (simple viscous)
    D = damp * qd

    # Solve: M qdd + C + G + D = 0  => qdd = -M^{-1} (C + G + D)
    rhs = -(C + G + D)
    qdd = np.linalg.solve(M, rhs)

    return np.array([w1, w2, w3, qdd[0], qdd[1], qdd[2]], dtype=float)

def simulate(params, y0, t_max=20.0, dt=1/240):
    """
    Integrate the equations of motion.
    """
    t_eval = np.arange(0.0, t_max + dt, dt)

    sol = solve_ivp(
        fun=lambda t, y: triple_pendulum_rhs(t, y, params),
        t_span=(0.0, t_max),
        y0=y0,
        t_eval=t_eval,
        method="DOP853",   # good for smooth problems; try "RK45" if needed
        rtol=1e-9,
        atol=1e-9,
        max_step=dt
    )
    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")
    return sol.t, sol.y

def angles_to_xy(params, th1, th2, th3):
    """
    Convert angles to cartesian positions for all bobs.
    """
    l1, l2, l3 = params["l1"], params["l2"], params["l3"]

    x1 = l1 * np.sin(th1)
    y1 = -l1 * np.cos(th1)

    x2 = x1 + l2 * np.sin(th2)
    y2 = y1 - l2 * np.cos(th2)

    x3 = x2 + l3 * np.sin(th3)
    y3 = y2 - l3 * np.cos(th3)

    return x1, y1, x2, y2, x3, y3

def animate_triple_pendulum(params, t, Y, trail=True, trail_length=600, fps=60):
    th1, th2, th3 = Y[0], Y[1], Y[2]

    x1, y1, x2, y2, x3, y3 = angles_to_xy(params, th1, th2, th3)

    L = params["l1"] + params["l2"] + params["l3"]
    pad = 0.15 * L
    lim = L + pad

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.grid(True, alpha=0.25)
    ax.set_title("Triple Pendulum")

    # Rod line
    line, = ax.plot([], [], lw=2)

    # Bobs
    bob1, = ax.plot([], [], marker="o", markersize=8)
    bob2, = ax.plot([], [], marker="o", markersize=8)
    bob3, = ax.plot([], [], marker="o", markersize=10)

    # Trail for third bob
    if trail:
        trail_line, = ax.plot([], [], lw=1, alpha=0.7)
        trail_x, trail_y = [], []
    else:
        trail_line = None

    # Time text
    time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

    # Determine frame step for desired fps
    dt = t[1] - t[0]
    step = max(1, int(round((1.0 / fps) / dt)))

    def init():
        line.set_data([], [])
        bob1.set_data([], [])
        bob2.set_data([], [])
        bob3.set_data([], [])
        if trail_line is not None:
            trail_line.set_data([], [])
        time_text.set_text("")
        return (line, bob1, bob2, bob3, time_text) if trail_line is None else (line, bob1, bob2, bob3, trail_line, time_text)

    def update(frame_idx):
        i = frame_idx * step
        if i >= len(t):
            i = len(t) - 1

        xs = [0.0, x1[i], x2[i], x3[i]]
        ys = [0.0, y1[i], y2[i], y3[i]]

        line.set_data(xs, ys)
        bob1.set_data([x1[i]], [y1[i]])
        bob2.set_data([x2[i]], [y2[i]])
        bob3.set_data([x3[i]], [y3[i]])

        if trail_line is not None:
            trail_x.append(x3[i])
            trail_y.append(y3[i])
            if len(trail_x) > trail_length:
                trail_x[:] = trail_x[-trail_length:]
                trail_y[:] = trail_y[-trail_length:]
            trail_line.set_data(trail_x, trail_y)

        time_text.set_text(f"t = {t[i]:.2f} s")
        return (line, bob1, bob2, bob3, time_text) if trail_line is None else (line, bob1, bob2, bob3, trail_line, time_text)

    nframes = int(np.ceil(len(t) / step))
    ani = FuncAnimation(fig, update, frames=nframes, init_func=init, blit=True, interval=1000/fps)

    plt.show()
    return ani

def main():
    # --- Parameters ---
    params = {
        "l1": 1.0, "l2": 1.0, "l3": 1.0,
        "m1": 1.0, "m2": 1.0, "m3": 1.0,
        "g": 9.81,
        "damp": 0.02,   # increase if it gets too wild numerically
    }

    # --- Initial state ---
    # Angles in radians. Try slightly different angles to see chaos explode.
    th1_0 = np.deg2rad(120.0)
    th2_0 = np.deg2rad(-10.0)
    th3_0 = np.deg2rad(40.0)

    w1_0 = 0.0
    w2_0 = 0.0
    w3_0 = 0.0

    y0 = np.array([th1_0, th2_0, th3_0, w1_0, w2_0, w3_0], dtype=float)

    # --- Simulate ---
    t, Y = simulate(params, y0, t_max=25.0, dt=1/300)

    # --- Animate ---
    animate_triple_pendulum(params, t, Y, trail=True, trail_length=800, fps=60)

if __name__ == "__main__":
    main()