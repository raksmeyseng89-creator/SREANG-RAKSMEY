# Root Finding Toolbox

A Python numerical root-finding toolbox that implements classical algorithms through a single, unified class interface.

---

## Project Description

`RootFindingProblem` solves equations of the form **f(x) = 0** using seven classical numerical methods. All algorithms are implemented from scratch using only the Python standard library and `cmath` — no external root-solver libraries are used.

---

## Implemented Methods

| Method | Description |
|---|---|
| **Bisection** | Repeatedly halves a bracket [a, b] where f changes sign |
| **Fixed-Point Iteration** | Iterates x = g(x) until convergence |
| **Newton–Raphson** | Uses f(x) and f′(x) for quadratic convergence |
| **Secant** | Two-point finite-difference approximation of Newton's method |
| **False Position (Regula Falsi)** | Bracketed method using a linear interpolant |
| **Steffensen** | Aitken-accelerated fixed-point iteration |
| **Horner's Method** | Efficient polynomial evaluation (used inside any solver) |
| **Muller's Method** | Three-point method that can find complex roots |

---

## How Each Algorithm Works

### Bisection
Given a bracket [a, b] with f(a)·f(b) < 0, the midpoint m = (a+b)/2 is computed. The sub-interval whose endpoints still have opposite signs is kept, halving the bracket each step. Guaranteed to converge.

### Fixed-Point Iteration
Rewrites f(x) = 0 as x = g(x) and iterates x_{n+1} = g(x_n). Converges when |g′(x)| < 1 near the root.

### Newton–Raphson
Uses the tangent line: x_{n+1} = x_n − f(x_n)/f′(x_n). Converges quadratically near a simple root. Requires the derivative f′.

### Secant Method
Approximates the derivative with a finite difference: x_{n+1} = x_n − f(x_n)·(x_n − x_{n-1})/(f(x_n) − f(x_{n-1})). No derivative needed.

### False Position (Regula Falsi)
Like bisection but uses a linear interpolant instead of the midpoint: c = b − f(b)·(b−a)/(f(b)−f(a)). Maintains a bracket at all times.

### Steffensen's Method
Applies Aitken's Δ² acceleration to fixed-point iteration, upgrading it to quadratic convergence without needing a derivative.

### Horner's Method
Evaluates a polynomial p(x) = a_n·x^n + … + a_0 in O(n) multiplications by the nested form p(x) = (…((a_n·x + a_{n-1})·x + a_{n-2})·x … + a_0). Can be plugged in as f inside any other solver.

### Muller's Method
Fits a parabola through three points (x_0, f(x_0)), (x_1, f(x_1)), (x_2, f(x_2)) and uses the quadratic formula to find the next approximation. Uses `cmath.sqrt` to handle complex discriminants, making it capable of finding complex roots.

---

## File Structure

```
root-finding-project/
├── root_finding.py   # RootFindingProblem class (all algorithms)
├── examples.py       # Runnable demonstrations of every method
└── README.md         # This file
```

---

## How to Run the Examples

No external libraries are required. Use Python 3.6+.

```bash
python examples.py
```

---

## Quick Code Example

```python
import math
from root_finding import RootFindingProblem

# --- Bisection ---
f = lambda x: x**3 - x - 2
p = RootFindingProblem(f=f)
print(p.solve("bisection", a=1, b=2))        # ≈ 1.5213797

# --- Newton's Method ---
df = lambda x: 3*x**2 - 1
p2 = RootFindingProblem(f=f, df=df)
print(p2.solve("newton", x0=1.5))            # ≈ 1.5213797

# --- Secant Method ---
print(p.solve("secant", x0=1, x1=2))         # ≈ 1.5213797

# --- False Position ---
print(p.solve("false_position", a=1, b=2))   # ≈ 1.5213797

# --- Fixed-Point Iteration ---
g = lambda x: (x + 2)**(1/3)
p3 = RootFindingProblem(g=g)
print(p3.solve("fixed_point", x0=1.5))       # ≈ 1.5213797

# --- Horner's Method (polynomial evaluation) ---
# Evaluates 2x^4 - 3x^3 + x^2 - 5x + 7 at x = 2
p4 = RootFindingProblem()
print(p4.solve("horner", coeffs=[2, -3, 1, -5, 7], x=2))  # 9

# --- Muller's Method (complex roots of x^2 + 1) ---
f_complex = lambda x: x**2 + 1
p5 = RootFindingProblem(f=f_complex)
print(p5.solve("muller", x0=0, x1=1, x2=2))  # ≈ 1j
```

---

## Constructor Parameters

```python
RootFindingProblem(f=None, df=None, g=None)
```

| Parameter | Description | Required by |
|---|---|---|
| `f` | Function f(x) whose root is sought | bisection, newton, secant, false_position, muller |
| `df` | Derivative f′(x) | newton |
| `g` | Fixed-point function g(x) | fixed_point, steffensen |

## Common Keyword Arguments

| Argument | Description | Default |
|---|---|---|
| `tol` | Convergence tolerance | `1e-10` |
| `max_iter` | Maximum number of iterations | `1000` |

---

## Error Handling

The class raises clear exceptions for:
- Invalid bracket (f(a) and f(b) same sign) in `bisection` / `false_position`
- Missing `df` in Newton's method
- Missing `g` in fixed-point / Steffensen methods
- Division by zero during iterations
- Non-convergence after `max_iter` steps
