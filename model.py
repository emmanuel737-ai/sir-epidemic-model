import numpy as np
from scipy.integrate import solve_ivp


def sir_rhs(t, y, beta, gamma, N):
    """
    Right-hand side of the SIR system.

    t     : current time (required by the solver; unused here)
    y     : [S, I, R] — the state at this instant
    beta  : transmission rate
    gamma : recovery rate (1 / infectious period)
    N     : total population

    Returns [dS/dt, dI/dt, dR/dt]
    """
    S, I, R = y

    new_infections = beta * S * I / N
    new_recoveries = gamma * I

    dS = -new_infections
    dI = new_infections - new_recoveries
    dR = new_recoveries

    return [dS, dI, dR]


def simulate_sir(beta, gamma, N, I0, days):
    """
    Integrate the SIR model forward in time.

    I0   : number infected at t = 0
    days : how many days to simulate

    Returns (t, S, I, R) as arrays.
    """
    S0 = N - I0
    R0_init = 0
    y0 = [S0, I0, R0_init]

    t_eval = np.linspace(0, days, days * 10)

    sol = solve_ivp(
        fun=sir_rhs,
        t_span=(0, days),
        y0=y0,
        args=(beta, gamma, N),
        t_eval=t_eval,
        method="RK45",
    )

    return sol.t, sol.y[0], sol.y[1], sol.y[2]


def sir_rhs_intervention(t, y, beta_before, beta_after, t_intervene, gamma, N):
    """
    SIR with a step change in transmission rate at t_intervene.

    This is where the `t` argument finally matters: beta is no longer
    a constant but a function of time.
    """
    S, I, R = y

    beta = beta_before if t < t_intervene else beta_after

    new_infections = beta * S * I / N
    new_recoveries = gamma * I

    return [-new_infections,
            new_infections - new_recoveries,
            new_recoveries]


def simulate_intervention(beta_before, beta_after, t_intervene,
                          gamma, N, I0, days):
    """
    Integrate the SIR model with an intervention part-way through.

    max_step is capped because adaptive solvers take large steps through
    smooth regions and can step straight over the discontinuity in beta.
    """
    y0 = [N - I0, I0, 0]
    t_eval = np.linspace(0, days, days * 10)

    sol = solve_ivp(
        fun=sir_rhs_intervention,
        t_span=(0, days),
        y0=y0,
        args=(beta_before, beta_after, t_intervene, gamma, N),
        t_eval=t_eval,
        method="RK45",
        max_step=0.1,
    )

    return sol.t, sol.y[0], sol.y[1], sol.y[2]


if __name__ == "__main__":
    t, S, I, R = simulate_sir(beta=1.5, gamma=0.3, N=763, I0=1, days=30)
    print("Peak infected:", I.max())
    print("Day of peak:", t[I.argmax()])
    print("Never infected:", S[-1])