import numpy as np
from numba import njit
from typing import Callable

@njit(cache=True)
def explicit_euler(x, y, f, n, h):
    for i in range(n):
        yi = y[i, :]
        y[i + 1, :] = yi + h * f(x[i + 1], yi)
    return y

@njit(cache=True)
def vec_norm(v):
    return np.sqrt(np.sum(v ** 2))

@njit(cache=True)
def implicit_euler(x, y, f, n, h):
    for i in range(n):
        yi = y[i, :]
        # Fixed-point iteration to solve implicit Euler method
        yn = yi + h * f(x[i + 1], yi)
        for _ in range(10):
            yn_ = yi + h * f(x[i + 1], yn)
            if vec_norm(yn_ - yn) < 1e-6:
                break
            yn = yn_
        y[i + 1, :] = yn_
    return y

@njit(cache=True)
def modified_midpoint(x, y, f, n, h):
    y[1, :] = y[0, :] + h * f(x[0], y[0, :])
    for i in range(1, n):
        yi = y[i, :]
        yp = y[i - 1, :]
        y[i + 1, :] = yp + 2 * h * f(x[i], yi)
    return y

@njit(cache=True)
def runge_kutta_3(x, y, f, n, h):
    for i in range(n):
        yi = y[i, :]
        k1 = f(x[i], yi)
        k2 = f(x[i] + 0.5 * h, yi + 0.5 * k1 * h)
        k3 = f(x[i] + 0.75 * h, yi + 0.75 * k2 * h)
        y[i + 1, :] = yi + h * (2 * k1 + 3 * k2 + 4 * k3) / 9.0
    return y

@njit(cache=True)
def runge_kutta_4(x, y, f, n, h):
    for i in range(n):
        yi = y[i, :]
        k1 = f(x[i], yi)
        k2 = f(x[i] + 0.5 * h, yi + 0.5 * k1 * h)
        k3 = f(x[i] + 0.5 * h, yi + 0.5 * k2 * h)
        k4 = f(x[i] + h, yi + k3 * h)
        y[i + 1, :] = yi + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return y

@njit(cache=True)
def apply(method: str, f: Callable, a: float, b: float, n: int, y0: np.ndarray):
    h = (b - a) / n
    x = np.arange(a, b + h, h)
    y = np.zeros((len(x), len(y0)))
    y[0, :] = y0

    if method == "Explicit Euler":
        return x, explicit_euler(x, y, f, n, h)
    elif method == "Implicit Euler":
        return x, implicit_euler(x, y, f, n, h)
    elif method == "Modified Midpoint":
        return x, modified_midpoint(x, y, f, n, h)
    elif method == "Runge-Kutta 3":
        return x, runge_kutta_3(x, y, f, n, h)
    elif method == "Runge-Kutta 4":
        return x, runge_kutta_4(x, y, f, n, h)
    else:
        raise ValueError(f"Unknown method: {method}")
