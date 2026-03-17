"""
examples.py – Demonstrates all methods in RootFindingProblem.
 
Run with:
    python examples.py
"""
 
import cmath
import math
from Root_Finding import RootFindingProblem
 
 
def separator(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print('='*60)
 
 
# ──────────────────────────────────────────────────────────────
# 1. Bisection Method
#    f(x) = x^3 - x - 2   root near x ≈ 1.5214
# ──────────────────────────────────────────────────────────────
separator("1. Bisection Method  –  f(x) = x³ - x - 2")
 
f1 = lambda x: x**3 - x - 2
p1 = RootFindingProblem(f=f1)
root = p1.solve("bisection", a=1, b=2)
print(f"  Root ≈ {root:.10f}")
print(f"  Verification f({root:.10f}) = {f1(root):.2e}")
 
 
# ──────────────────────────────────────────────────────────────
# 2. Fixed-Point Iteration
#    f(x) = x^3 - x - 2   →   g(x) = (x + 2)^(1/3)
# ──────────────────────────────────────────────────────────────
separator("2. Fixed-Point Iteration  –  g(x) = (x + 2)^(1/3)")
 
g2 = lambda x: (x + 2) ** (1/3)
p2 = RootFindingProblem(f=f1, g=g2)
root = p2.solve("fixed_point", x0=1.5)
print(f"  Root ≈ {root:.10f}")
print(f"  Verification f({root:.10f}) = {f1(root):.2e}")
 
 
# ──────────────────────────────────────────────────────────────
# 3. Newton's Method
#    f(x) = cos(x) - x    root near x ≈ 0.7391 (Dottie number)
# ──────────────────────────────────────────────────────────────
separator("3. Newton's Method  –  f(x) = cos(x) - x")
 
f3  = lambda x: math.cos(x) - x
df3 = lambda x: -math.sin(x) - 1
p3  = RootFindingProblem(f=f3, df=df3)
root = p3.solve("newton", x0=0.5)
print(f"  Root ≈ {root:.10f}")
print(f"  Verification f({root:.10f}) = {f3(root):.2e}")
 
 
# ──────────────────────────────────────────────────────────────
# 4. Secant Method
#    f(x) = e^x - 3x      roots near x ≈ 0.6190 and x ≈ 1.5121
# ──────────────────────────────────────────────────────────────
separator("4. Secant Method  –  f(x) = eˣ - 3x")
 
f4 = lambda x: math.exp(x) - 3*x
p4 = RootFindingProblem(f=f4)
root = p4.solve("secant", x0=0, x1=1)
print(f"  Root ≈ {root:.10f}")
print(f"  Verification f({root:.10f}) = {f4(root):.2e}")
 
 
# ──────────────────────────────────────────────────────────────
# 5. Method of False Position (Regula Falsi)
#    f(x) = x^2 - 4       roots at x = ±2
# ──────────────────────────────────────────────────────────────
separator("5. False Position  –  f(x) = x² - 4")
 
f5 = lambda x: x**2 - 4
p5 = RootFindingProblem(f=f5)
root = p5.solve("false_position", a=0, b=3)
print(f"  Root ≈ {root:.10f}")
print(f"  Verification f({root:.10f}) = {f5(root):.2e}")
 
 
# ──────────────────────────────────────────────────────────────
# 6. Horner's Method
#    p(x) = 2x^4 - 3x^3 + x^2 - 5x + 7   evaluated at x = 2
# ──────────────────────────────────────────────────────────────
separator("6. Horner's Method  –  p(x) = 2x⁴ - 3x³ + x² - 5x + 7  at x = 2")
 
coeffs = [2, -3, 1, -5, 7]   # highest degree first
p6 = RootFindingProblem()
val = p6.solve("horner", coeffs=coeffs, x=2)
# Manual verification: 2*16 - 3*8 + 4 - 10 + 7 = 32 - 24 + 4 - 10 + 7 = 9
print(f"  p(2) = {val}  (expected 9)")
 
# Using Horner inside a solver: find a root of p(x) = x^3 - 6x^2 + 11x - 6
# roots are x = 1, 2, 3
poly_coeffs = [1, -6, 11, -6]
f6_poly  = lambda x: p6.solve("horner", coeffs=poly_coeffs, x=x)
df6_poly = lambda x: p6.solve("horner", coeffs=[3, -12, 11], x=x)   # derivative
p6b = RootFindingProblem(f=f6_poly, df=df6_poly)
root = p6b.solve("newton", x0=3.5)
print(f"  Root of x³-6x²+11x-6 near x=3.5  →  {root:.10f}  (expected 3.0)")
 
 
# ──────────────────────────────────────────────────────────────
# 7. Muller's Method – complex roots
#    p(x) = x^2 + 1   roots: x = ±i
# ──────────────────────────────────────────────────────────────
separator("7. Muller's Method  –  f(x) = x² + 1  (complex roots ±i)")
 
f7 = lambda x: x**2 + 1
p7 = RootFindingProblem(f=f7)
root = p7.solve("muller", x0=0, x1=1, x2=2)
print(f"  Root ≈ {root}")
print(f"  Verification |f(root)| = {abs(f7(root)):.2e}")
 
# Another example: x^3 - 1  (complex cube roots of unity)
separator("7b. Muller's Method  –  f(x) = x³ - 1  (complex cube roots)")
f7b = lambda x: x**3 - 1
p7b = RootFindingProblem(f=f7b)
root2 = p7b.solve("muller", x0=0.5+0.5j, x1=0.8+0.6j, x2=0.3+0.9j)
print(f"  Root ≈ {root2}")
print(f"  Verification |f(root)| = {abs(f7b(root2)):.2e}")
 
 
# ──────────────────────────────────────────────────────────────
# Bonus: Steffensen's Method
# ──────────────────────────────────────────────────────────────
separator("Bonus – Steffensen's Method  –  g(x) = (x + 2)^(1/3)")
 
p8 = RootFindingProblem(g=g2)
root = p8.solve("steffensen", x0=1.5)
print(f"  Root ≈ {root:.10f}")
print(f"  Verification f({root:.10f}) = {f1(root):.2e}")
 
 
print("\n" + "="*60)
print("  All examples completed successfully.")
print("="*60 + "\n")