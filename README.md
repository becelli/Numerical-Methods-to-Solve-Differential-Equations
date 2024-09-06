# Predator-Prey Model Simulation

This project simulates the classic predator-prey model using various numerical methods. The results are visualized with both 2D and 3D plots, and different numerical methods are benchmarked for performance in terms of execution time and accuracy.

## Project Structure

```
__pycache__/
.gitignore
metrics.ipynb
outputs/
requirements.txt
solver.py
```

- **[`metrics.ipynb`](metrics.ipynb)**: Jupyter notebook that runs simulations, generates plots, and benchmarks numerical methods.
- **[`solver.py`](solver.py)**: Python module that implements various numerical methods for solving differential equations.
- **[`requirements.txt`](requirements.txt)**: Contains the list of dependencies required for the project.
- **[`outputs/`](outputs/)**: Directory where output plots and LaTeX files are stored.

## Requirements

The project depends on the following Python packages:

- `numpy`
- `matplotlib`
- `numba`

Install the required packages with:

```sh
pip install -r requirements.txt
```

## Usage

1. **Run Simulations**: Open [`metrics.ipynb`](metrics.ipynb) in Jupyter Notebook or Jupyter Lab and execute the cells to run simulations and generate plots.

2. **Numerical Methods**: The following methods are implemented in [`solver.py`](solver.py):
   - **Explicit Euler**
   - **Implicit Euler**
   - **Modified Midpoint**
   - **Runge-Kutta 3**
   - **Runge-Kutta 4**

3. **Benchmarking**: The notebook benchmarks the execution time and accuracy of each numerical method. Benchmark results are saved in the [`outputs/`](outputs/) directory.

## Methods Overview

### [`solver.py`](solver.py)

- **`explicit_euler(x, y, f, n, h)`**: Implements the Explicit Euler method for solving ordinary differential equations (ODEs). This method approximates the solution at each step using the current values and the derivative function `f`. The formula used is:

    $$y_{i+1} = y_i + h \cdot f(x_i, y_i)$$

    where $h$ is the step size, $x_i$ is the current value of the independent variable, and $y_i$ is the current value of the dependent variable.

- **`implicit_euler(x, y, f, n, h)`**: Implements the Implicit Euler method, a variant of the Euler method that uses the future value of the dependent variable to compute the next step. This method is more stable for stiff equations. The formula used is:

    $$y_{i+1} = y_i + h \cdot f(x_{i+1}, y_{i+1})$$

    This requires solving an implicit equation at each step, often using iterative solvers.

- **`modified_midpoint(x, y, f, n, h)`**: Implements the Modified Midpoint method, which is a second-order Runge-Kutta method. It improves on the basic midpoint method by incorporating a correction step. The method uses:

    $$k_1 = h \cdot f(x_i, y_i)$$

    $$k_2 = h \cdot f\left(x_i + \frac{h}{2}, y_i + \frac{k_1}{2}\right)$$

    $$y_{i+1} = y_i + k_2$$

    This approach provides better accuracy by averaging the slopes.

- **`runge_kutta_3(x, y, f, n, h)`**: Implements the third-order Runge-Kutta method (RK3), which offers higher accuracy than Euler methods by using three intermediate steps to estimate the next value. The method involves:

    $$k_1 = h \cdot f(x_i, y_i)$$

    $$k_2 = h \cdot f\left(x_i + \frac{h}{2}, y_i + \frac{k_1}{2}\right)$$

    $$k_3 = h \cdot f\left(x_i + h, y_i - k_1 + 2k_2\right)$$

    $$y_{i+1} = y_i + \frac{1}{6}(k_1 + 4k_2 + k_3)$$

- **`runge_kutta_4(x, y, f, n, h)`**: Implements the fourth-order Runge-Kutta method (RK4), known for its accuracy and stability. It calculates the next value by combining four intermediate estimates:

    $$k_1 = h \cdot f(x_i, y_i)$$

    $$k_2 = h \cdot f\left(x_i + \frac{h}{2}, y_i + \frac{k_1}{2}\right)$$

    $$k_3 = h \cdot f\left(x_i + \frac{h}{2}, y_i + \frac{k_2}{2}\right)$$

    $$k_4 = h \cdot f(x_i + h, y_i + k_3)$$

    $$y_{i+1} = y_i + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

    This method is widely used for solving ODEs due to its balance of accuracy and computational efficiency. It is implemented in various scientific computing libraries, such as SciPy.