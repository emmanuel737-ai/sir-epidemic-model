import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from model import simulate_sir

N = 763
I0 = 1
GAMMA = 0.45      # fixed: ~2.2 day infectious period for influenza


def load_data(path="data/boarding_school_flu.csv"):
    df = pd.read_csv(path)
    return df["day"].values, df["infected"].values


def model_at_days(beta, days):
    """Run the model, return predicted I at the observed days."""
    t, S, I, R = simulate_sir(beta, GAMMA, N, I0, days=int(days.max()) + 1)
    return np.interp(days, t, I)


def residuals(params, days, observed):
    """Difference between model prediction and observation."""
    beta = params[0]
    predicted = model_at_days(beta, days)

    return predicted - observed


if __name__ == "__main__":
    days, observed = load_data()

    result = least_squares(residuals, x0=[2.0], args=(days, observed))

    beta_hat = result.x[0]
    print("Estimated beta:", round(beta_hat, 3))
    print("Fixed gamma:", GAMMA)
    print("Implied R0:", round(beta_hat / GAMMA, 2))