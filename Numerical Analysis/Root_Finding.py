import cmath
 
 
class RootFindingProblem:
    """
    A numerical root-finding toolbox implementing classical algorithms.
 
    Parameters
    ----------
    f  : callable, optional  – function whose root is sought f(x) = 0
    df : callable, optional  – derivative f'(x)  (Newton's method)
    g  : callable, optional  – fixed-point function g(x)  (fixed-point / Steffensen)
    """
 
    def __init__(self, f=None, df=None, g=None):
        self.f = f
        self.df = df
        self.g = g
 
    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
 
    def solve(self, method: str, **kwargs):
        """
        Unified solver interface.
 
        Parameters
        ----------
        method : str
            One of 'bisection', 'fixed_point', 'newton', 'secant',
            'false_position', 'steffensen', 'horner', 'muller'.
        **kwargs
            Method-specific keyword arguments (see each private method).
 
        Returns
        -------
        float or complex
            Approximate root.
        """
        dispatch = {
            "bisection":      self._bisection,
            "fixed_point":    self._fixed_point,
            "newton":         self._newton,
            "secant":         self._secant,
            "false_position": self._false_position,
            "steffensen":     self._steffensen,
            "horner":         self._horner,
            "muller":         self._muller,
        }
        method = method.lower().replace("-", "_")
        if method not in dispatch:
            raise ValueError(
                f"Unknown method '{method}'. Choose from: {list(dispatch.keys())}"
            )
        return dispatch[method](**kwargs)
 
    # ------------------------------------------------------------------
    # Private methods
    # ------------------------------------------------------------------
 
    def _bisection(self, a, b, tol=1e-10, max_iter=1000):
        """
        Bisection method.
 
        Parameters
        ----------
        a, b     : float  – bracket endpoints with f(a)*f(b) < 0
        tol      : float  – convergence tolerance
        max_iter : int    – maximum iterations
        """
        if self.f is None:
            raise ValueError("f(x) must be provided for bisection.")
        fa, fb = self.f(a), self.f(b)
        if fa * fb > 0:
            raise ValueError(
                "Invalid interval: f(a) and f(b) must have opposite signs."
            )
        for i in range(max_iter):
            mid = (a + b) / 2.0
            fm = self.f(mid)
            if abs(fm) < tol or (b - a) / 2.0 < tol:
                return mid
            if fa * fm < 0:
                b, fb = mid, fm
            else:
                a, fa = mid, fm
        raise RuntimeError(f"Bisection did not converge after {max_iter} iterations.")
 
    def _fixed_point(self, x0, tol=1e-10, max_iter=1000):
        """
        Fixed-point iteration  x_{n+1} = g(x_n).
 
        Parameters
        ----------
        x0       : float  – initial guess
        tol      : float  – convergence tolerance
        max_iter : int    – maximum iterations
        """
        if self.g is None:
            raise ValueError("g(x) must be provided for fixed-point iteration.")
        x = x0
        for i in range(max_iter):
            x_new = self.g(x)
            if abs(x_new - x) < tol:
                return x_new
            x = x_new
        raise RuntimeError(
            f"Fixed-point iteration did not converge after {max_iter} iterations."
        )
 
    def _newton(self, x0, tol=1e-10, max_iter=1000):
        """
        Newton–Raphson method.
 
        Parameters
        ----------
        x0       : float  – initial guess
        tol      : float  – convergence tolerance
        max_iter : int    – maximum iterations
        """
        if self.f is None:
            raise ValueError("f(x) must be provided for Newton's method.")
        if self.df is None:
            raise ValueError("df(x) (derivative) must be provided for Newton's method.")
        x = x0
        for i in range(max_iter):
            fx = self.f(x)
            dfx = self.df(x)
            if dfx == 0:
                raise ZeroDivisionError(
                    f"Derivative is zero at x = {x}. Newton's method failed."
                )
            x_new = x - fx / dfx
            if abs(x_new - x) < tol:
                return x_new
            x = x_new
        raise RuntimeError(
            f"Newton's method did not converge after {max_iter} iterations."
        )
 
    def _secant(self, x0, x1, tol=1e-10, max_iter=1000):
        """
        Secant method.
 
        Parameters
        ----------
        x0, x1   : float  – two initial approximations
        tol      : float  – convergence tolerance
        max_iter : int    – maximum iterations
        """
        if self.f is None:
            raise ValueError("f(x) must be provided for the secant method.")
        for i in range(max_iter):
            f0, f1 = self.f(x0), self.f(x1)
            denom = f1 - f0
            if denom == 0:
                raise ZeroDivisionError(
                    "Division by zero in secant method (f(x1) == f(x0))."
                )
            x2 = x1 - f1 * (x1 - x0) / denom
            if abs(x2 - x1) < tol:
                return x2
            x0, x1 = x1, x2
        raise RuntimeError(
            f"Secant method did not converge after {max_iter} iterations."
        )
 
    def _false_position(self, a, b, tol=1e-10, max_iter=1000):
        """
        Method of false position (Regula Falsi).
 
        Parameters
        ----------
        a, b     : float  – bracket endpoints with f(a)*f(b) < 0
        tol      : float  – convergence tolerance
        max_iter : int    – maximum iterations
        """
        if self.f is None:
            raise ValueError("f(x) must be provided for false position.")
        fa, fb = self.f(a), self.f(b)
        if fa * fb > 0:
            raise ValueError(
                "Invalid interval: f(a) and f(b) must have opposite signs."
            )
        for i in range(max_iter):
            denom = fb - fa
            if denom == 0:
                raise ZeroDivisionError(
                    "Division by zero in false position (f(b) == f(a))."
                )
            c = b - fb * (b - a) / denom
            fc = self.f(c)
            if abs(fc) < tol:
                return c
            if fa * fc < 0:
                b, fb = c, fc
            else:
                a, fa = c, fc
        raise RuntimeError(
            f"False position did not converge after {max_iter} iterations."
        )
 
    def _steffensen(self, x0, tol=1e-10, max_iter=1000):
        """
        Steffensen's method (accelerated fixed-point using Aitken Δ² process).
 
        Parameters
        ----------
        x0       : float  – initial guess
        tol      : float  – convergence tolerance
        max_iter : int    – maximum iterations
        """
        if self.g is None:
            raise ValueError("g(x) must be provided for Steffensen's method.")
        x = x0
        for i in range(max_iter):
            gx  = self.g(x)
            ggx = self.g(gx)
            denom = ggx - 2 * gx + x
            if denom == 0:
                raise ZeroDivisionError(
                    "Division by zero in Steffensen's method."
                )
            x_new = x - (gx - x) ** 2 / denom
            if abs(x_new - x) < tol:
                return x_new
            x = x_new
        raise RuntimeError(
            f"Steffensen's method did not converge after {max_iter} iterations."
        )
 
    def _horner(self, coeffs, x):
        """
        Horner's method for efficient polynomial evaluation.
 
        Evaluates  p(x) = coeffs[0]*x^n + coeffs[1]*x^(n-1) + ... + coeffs[n]
 
        Parameters
        ----------
        coeffs : list of float  – polynomial coefficients, highest degree first
        x      : float or complex – evaluation point
 
        Returns
        -------
        float or complex – p(x)
        """
        if not coeffs:
            raise ValueError("Coefficient list must not be empty.")
        result = coeffs[0]
        for c in coeffs[1:]:
            result = result * x + c
        return result
 
    def _muller(self, x0, x1, x2, tol=1e-10, max_iter=1000):
        """
        Muller's method – finds real or complex roots.
 
        Parameters
        ----------
        x0, x1, x2 : float or complex  – three distinct initial approximations
        tol        : float             – convergence tolerance
        max_iter   : int               – maximum iterations
        """
        if self.f is None:
            raise ValueError("f(x) must be provided for Muller's method.")
 
        x0, x1, x2 = complex(x0), complex(x1), complex(x2)
 
        for i in range(max_iter):
            f0, f1, f2 = self.f(x0), self.f(x1), self.f(x2)
 
            h0 = x1 - x0
            h1 = x2 - x1
            delta0 = (f1 - f0) / h0 if h0 != 0 else 0
            delta1 = (f2 - f1) / h1 if h1 != 0 else 0
 
            denom_a = h0 + h1
            if denom_a == 0:
                raise ZeroDivisionError("Division by zero computing 'a' in Muller's method.")
 
            a = (delta1 - delta0) / denom_a
            b = a * h1 + delta1
            c = f2
 
            discriminant = b * b - 4 * a * c
            sqrt_disc = cmath.sqrt(discriminant)
 
            denom_plus  = b + sqrt_disc
            denom_minus = b - sqrt_disc
            denom = denom_plus if abs(denom_plus) >= abs(denom_minus) else denom_minus
 
            if denom == 0:
                raise ZeroDivisionError("Division by zero in Muller's method (denominator).")
 
            dx = -2 * c / denom
            x3 = x2 + dx
 
            if abs(dx) < tol:
                # Return real if imaginary part is negligible
                if abs(x3.imag) < tol:
                    return x3.real
                return x3
 
            x0, x1, x2 = x1, x2, x3
 
        raise RuntimeError(
            f"Muller's method did not converge after {max_iter} iterations."
        )