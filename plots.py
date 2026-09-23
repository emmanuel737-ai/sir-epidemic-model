import numpy as np
import matplotlib.pyplot as plt

from model import simulate_sir, simulate_intervention
from fitting import load_data, GAMMA, N, I0


def plot_sir(beta, gamma, N, I0, days, outfile="figures/sir_curves.png"):
    """Plot S, I and R over time for one parameter set."""
    t, S, I, R = simulate_sir(beta, gamma, N, I0, days)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(t, S, label="Susceptible")
    ax.plot(t, I, label="Infectious")
    ax.plot(t, R, label="Recovered")

    ax.set_xlabel("Days")
    ax.set_ylabel("People")
    ax.set_title(f"SIR model (R0 = {beta / gamma:.1f})")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_r0_sweep(gamma, N, I0, days, r0_values,
                  outfile="figures/r0_sweep.png"):
    """Plot the infectious curve across a range of R0 values."""
    fig, ax = plt.subplots(figsize=(8, 5))

    for r0 in r0_values:
        beta = r0 * gamma
        t, S, I, R = simulate_sir(beta, gamma, N, I0, days)
        ax.plot(t, I, label=f"R0 = {r0}")

    ax.set_xlabel("Days")
    ax.set_ylabel("Infectious")
    ax.set_title("Effect of R0 on outbreak size")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_fit(beta, outfile="figures/fit.png"):
    """Plot the fitted model against the observed 1978 outbreak data."""
    days, observed = load_data()

    t, S, I, R = simulate_sir(beta, GAMMA, N, I0, days=16)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(days, observed, color="black", zorder=3,
               label="Observed (1978 outbreak)")
    ax.plot(t, I, color="crimson",
            label=f"SIR fit (R0 = {beta / GAMMA:.2f})")

    ax.set_xlabel("Day")
    ax.set_ylabel("Boys confined to bed")
    ax.set_title("SIR model fitted to boarding school influenza data")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_intervention(beta, gamma, N, I0, days, t_intervene, reduction,
                      outfile="figures/intervention.png"):
    """Compare an unmitigated outbreak with one where beta drops."""
    t0, S0, I_base, R_base = simulate_sir(beta, gamma, N, I0, days)

    beta_after = beta * (1 - reduction)
    t1, S1, I1, R1 = simulate_intervention(
        beta, beta_after, t_intervene, gamma, N, I0, days)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(t0, I_base, label="No intervention")
    ax.plot(t1, I1, label=f"{int(reduction * 100)}% drop on day {t_intervene}")
    ax.axvline(t_intervene, color="grey", linestyle="--", alpha=0.6)

    ax.set_xlabel("Day")
    ax.set_ylabel("Infectious")
    ax.set_title(f"Transmission cut by {int(reduction * 100)}% "
                 f"on day {t_intervene}")
    ax.set_ylim(0, 320)
    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)

    print(f"--- intervention on day {t_intervene} ---")
    print(f"Peak without intervention: {I_base.max():.0f}")
    print(f"Peak with intervention:    {I1.max():.0f}")
    print(f"Total infected, no action: {R_base[-1]:.0f}")
    print(f"Total infected, with:      {R1[-1]:.0f}")
    print(f"Saved {outfile}")


if __name__ == "__main__":
    plot_sir(beta=1.5, gamma=0.3, N=763, I0=1, days=30)

    plot_r0_sweep(gamma=0.3, N=763, I0=1, days=120,
                  r0_values=[0.8, 1.0, 1.5, 2.5, 5.0])

    plot_fit(beta=1.673)

    plot_intervention(beta=1.673, gamma=0.45, N=763, I0=1, days=30,
                      t_intervene=4, reduction=0.6,
                      outfile="figures/intervention_day4.png")

    plot_intervention(beta=1.673, gamma=0.45, N=763, I0=1, days=30,
                      t_intervene=8, reduction=0.6,
                      outfile="figures/intervention_day8.png")