"""Polynômes sur F_p et factorisation de Berlekamp."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, List, Optional, Tuple


@dataclass(frozen=True)
class Poly:
    """Polynôme à coefficients dans Z/pZ, coeffs[i] = coefficient de X^i."""

    coeffs: Tuple[int, ...]
    p: int

    def __post_init__(self) -> None:
        if self.p < 2:
            raise ValueError("p doit être >= 2")
        normalized = tuple(c % self.p for c in self.coeffs)
        object.__setattr__(self, "coeffs", _trim(normalized))

    @classmethod
    def zero(cls, p: int) -> Poly:
        return cls((0,), p)

    @classmethod
    def one(cls, p: int) -> Poly:
        return cls((1,), p)

    @classmethod
    def monomial(cls, p: int, degree: int, coeff: int = 1) -> Poly:
        coeffs = [0] * (degree + 1)
        coeffs[degree] = coeff % p
        return cls(tuple(coeffs), p)

    @classmethod
    def x(cls, p: int) -> Poly:
        return cls.monomial(p, 1)

    def degree(self) -> int:
        if len(self.coeffs) == 1 and self.coeffs[0] == 0:
            return -1
        return len(self.coeffs) - 1

    def leading_coeff(self) -> int:
        d = self.degree()
        return 0 if d < 0 else self.coeffs[d] % self.p

    def is_zero(self) -> bool:
        return self.degree() < 0

    def is_one(self) -> bool:
        return self.degree() == 0 and self.coeffs[0] == 1

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Poly):
            return NotImplemented
        return self.p == other.p and self.coeffs == other.coeffs

    def __hash__(self) -> int:
        return hash((self.p, self.coeffs))

    def __neg__(self) -> Poly:
        return Poly(tuple((-c) % self.p for c in self.coeffs), self.p)

    def __add__(self, other: Poly) -> Poly:
        _check_same_field(self, other)
        n = max(len(self.coeffs), len(other.coeffs))
        coeffs = [
            (self.coeffs[i] if i < len(self.coeffs) else 0)
            + (other.coeffs[i] if i < len(other.coeffs) else 0)
            for i in range(n)
        ]
        return Poly(tuple(coeffs), self.p)

    def __sub__(self, other: Poly) -> Poly:
        return self + (-other)

    def __mul__(self, other: Poly) -> Poly:
        _check_same_field(self, other)
        if self.is_zero() or other.is_zero():
            return Poly.zero(self.p)
        result = [0] * (len(self.coeffs) + len(other.coeffs) - 1)
        for i, a in enumerate(self.coeffs):
            if a == 0:
                continue
            for j, b in enumerate(other.coeffs):
                if b == 0:
                    continue
                result[i + j] = (result[i + j] + a * b) % self.p
        return Poly(tuple(result), self.p)

    def scalar_mul(self, scalar: int) -> Poly:
        s = scalar % self.p
        if s == 0:
            return Poly.zero(self.p)
        return Poly(tuple((s * c) % self.p for c in self.coeffs), self.p)

    def divmod(self, divisor: Poly) -> Tuple[Poly, Poly]:
        _check_same_field(self, divisor)
        if divisor.is_zero():
            raise ZeroDivisionError("division par zéro")
        p = self.p
        remainder = Poly(self.coeffs, p)
        divisor_leading = divisor.leading_coeff()
        inv_leading = pow(divisor_leading, -1, p)
        quotient_coeffs: List[int] = []

        while remainder.degree() >= divisor.degree() and not remainder.is_zero():
            deg_diff = remainder.degree() - divisor.degree()
            factor = (remainder.leading_coeff() * inv_leading) % p
            quotient_coeffs.append((deg_diff, factor))
            shift = divisor.scalar_mul(factor)
            shift = _shift_poly(shift, deg_diff)
            remainder = remainder - shift

        quotient = Poly.zero(p)
        for deg, coeff in reversed(quotient_coeffs):
            quotient = quotient + Poly.monomial(p, deg, coeff)
        return quotient, remainder

    def __mod__(self, other: Poly) -> Poly:
        return self.divmod(other)[1]

    def powmod(self, exponent: int, modulus: Poly) -> Poly:
        if exponent < 0:
            raise ValueError("exposant négatif interdit dans powmod")
        result = Poly.one(self.p)
        base = self % modulus
        e = exponent
        while e > 0:
            if e & 1:
                result = (result * base) % modulus
            base = (base * base) % modulus
            e >>= 1
        return result

    def gcd(self, other: Poly) -> Poly:
        _check_same_field(self, other)
        a, b = self, other
        while not b.is_zero():
            a, b = b, a % b
        leading = a.leading_coeff()
        if leading == 0:
            return a
        inv = pow(leading, -1, self.p)
        return a.scalar_mul(inv)

    def derivative(self) -> Poly:
        if self.degree() <= 0:
            return Poly.zero(self.p)
        coeffs = [(i * self.coeffs[i]) % self.p for i in range(1, len(self.coeffs))]
        return Poly(tuple(coeffs), self.p)

    def evaluate(self, x: int) -> int:
        result = 0
        power = 1
        for coeff in self.coeffs:
            result = (result + coeff * power) % self.p
            power = (power * x) % self.p
        return result

    def __repr__(self) -> str:
        if self.is_zero():
            return "0"
        terms = []
        for i, c in enumerate(self.coeffs):
            if c % self.p == 0:
                continue
            if i == 0:
                terms.append(str(c % self.p))
            elif i == 1:
                terms.append("X" if c % self.p == 1 else f"{c % self.p}*X")
            else:
                terms.append("X^%d" % i if c % self.p == 1 else f"{c % self.p}*X^{i}")
        return " + ".join(reversed(terms))


def _trim(coeffs: Tuple[int, ...]) -> Tuple[int, ...]:
    i = len(coeffs)
    while i > 1 and coeffs[i - 1] == 0:
        i -= 1
    return coeffs[:i] if i > 0 else (0,)


def _check_same_field(a: Poly, b: Poly) -> None:
    if a.p != b.p:
        raise ValueError("polynômes sur des corps différents")


def _shift_poly(poly: Poly, degree: int) -> Poly:
    if degree == 0:
        return poly
    return Poly(tuple([0] * degree + list(poly.coeffs)), poly.p)


def poly_add(a: Poly, b: Poly) -> Poly:
    return a + b


def poly_sub(a: Poly, b: Poly) -> Poly:
    return a - b


def poly_mul(a: Poly, b: Poly) -> Poly:
    return a * b


def poly_mod(a: Poly, b: Poly) -> Poly:
    return a % b


def poly_powmod(base: Poly, exp: int, mod: Poly) -> Poly:
    return base.powmod(exp, mod)


def is_irreducible(f: Poly) -> bool:
    """Test d'irréductibilité via gcd(f, X^{p^i} - X) pour i = 1..deg(f)/2."""
    if f.is_zero():
        return False
    deg = f.degree()
    if deg <= 0:
        return False
    if deg == 1:
        return True
    p = f.p
    x = Poly.x(p)
    x_p_i = x
    half = deg // 2
    for _ in range(1, half + 1):
        x_p_i = x_p_i.powmod(p, f)
        diff = x_p_i - x
        if f.gcd(diff).degree() > 0:
            return False
    return True


def _matrix_nullspace(mod_p: int, matrix: List[List[int]]) -> List[List[int]]:
    """Base du noyau de (M - I) sur F_p, lignes = vecteurs de coefficients."""
    n = len(matrix)
    aug = [[(matrix[i][j] - (1 if i == j else 0)) % mod_p for j in range(n)] for i in range(n)]
    pivot_cols: List[int] = []
    row = 0
    for col in range(n):
        pivot = None
        for r in range(row, n):
            if aug[r][col] % mod_p != 0:
                pivot = r
                break
        if pivot is None:
            continue
        aug[row], aug[pivot] = aug[pivot], aug[row]
        inv = pow(aug[row][col], -1, mod_p)
        aug[row] = [(v * inv) % mod_p for v in aug[row]]
        for r in range(n):
            if r != row and aug[r][col] != 0:
                factor = aug[r][col]
                aug[r] = [(aug[r][c] - factor * aug[row][c]) % mod_p for c in range(n)]
        pivot_cols.append(col)
        row += 1
        if row == n:
            break

    free_cols = [c for c in range(n) if c not in pivot_cols]
    basis: List[List[int]] = []
    for free in free_cols:
        vec = [0] * n
        vec[free] = 1
        for i, pc in enumerate(pivot_cols):
            if i < row:
                vec[pc] = (-aug[i][free]) % mod_p
        basis.append(vec)
    return basis


def _poly_from_nullvector(vec: List[int], p: int) -> Poly:
    return Poly(tuple(vec), p)


def berlekamp_factor(f: Poly) -> List[Poly]:
    """
    Factorisation de Berlekamp d'un polynôme carré-libre sur F_p.
    Retourne des facteurs irréductibles moniques.
    """
    if f.is_zero():
        raise ValueError("polynôme nul")
    p = f.p
    deg = f.degree()
    if deg <= 1:
        leading = f.leading_coeff()
        inv = pow(leading, -1, p)
        return [f.scalar_mul(inv)]

    # Rendre f monique
    inv_leading = pow(f.leading_coeff(), -1, p)
    f = f.scalar_mul(inv_leading)

    # Construire la matrice Q : ligne i = coeffs de X^{p*i} mod f
    n = deg
    q_rows: List[List[int]] = []
    x = Poly.x(p)
    power = Poly.one(p)
    for _ in range(n):
        coeffs = list(power.coeffs) + [0] * (n - len(power.coeffs))
        q_rows.append([coeffs[i] % p for i in range(n)])
        power = power.powmod(p, f)

    nullspace = _matrix_nullspace(p, q_rows)
    if not nullspace:
        return [f]

    factors = [f]
    for vec in nullspace:
        h = _poly_from_nullvector(vec, p)
        if h.is_zero():
            continue
        new_factors: List[Poly] = []
        for factor in factors:
            if factor.degree() <= 1:
                new_factors.append(factor)
                continue
            split = _split_with_h(factor, h)
            new_factors.extend(split)
        factors = new_factors
        if all(g.degree() <= 1 for g in factors):
            break

    monic_factors = []
    for g in factors:
        if g.degree() <= 0:
            continue
        inv = pow(g.leading_coeff(), -1, p)
        monic_factors.append(g.scalar_mul(inv))
    return monic_factors


def _split_with_h(f: Poly, h: Poly) -> List[Poly]:
    """Sépare f en utilisant gcd(f, h - c)."""
    p = f.p
    result: List[Poly] = [f]
    for c in range(p):
        next_result: List[Poly] = []
        hc = h - Poly((c % p,), p)
        for factor in result:
            if factor.degree() <= 1:
                next_result.append(factor)
                continue
            g = factor.gcd(hc)
            if g.degree() == 0 or g.degree() == factor.degree():
                next_result.append(factor)
            else:
                next_result.append(g)
                q, _ = factor.divmod(g)
                next_result.append(q)
        result = next_result
    return [fac for fac in result if fac.degree() > 0]


def is_square_free(f: Poly) -> bool:
    g = f.gcd(f.derivative())
    return g.degree() == 0


def factor(f: Poly) -> List[Poly]:
    """Factorise f en irréductibles moniques (via Berlekamp si carré-libre)."""
    if f.is_zero():
        raise ValueError("polynôme nul")
    p = f.p
    inv_leading = pow(f.leading_coeff(), -1, p)
    f = f.scalar_mul(inv_leading)

    if f.degree() <= 1:
        return [f]

    if not is_square_free(f):
        # Extraire les multiplicités via gcd(f, f')
        g = f.gcd(f.derivative())
        if not g.is_one():
            q, _ = f.divmod(g)
            return factor(g) + factor(q)

    return berlekamp_factor(f)


def prime_factors(n: int) -> List[int]:
    """Facteurs premiers distincts de n."""
    if n < 2:
        return []
    factors: List[int] = []
    d = 2
    m = n
    while d * d <= m:
        if m % d == 0:
            factors.append(d)
            while m % d == 0:
                m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        factors.append(m)
    return factors


def iter_monic_polynomials(p: int, degree: int) -> Iterator[Poly]:
    """Énumère les polynômes moniques de degré `degree` sur F_p."""
    if degree < 1:
        return
    inner_count = p ** degree
    for idx in range(inner_count):
        inner = []
        x = idx
        for _ in range(degree):
            inner.append(x % p)
            x //= p
        coeffs = inner + [1]
        yield Poly(tuple(coeffs), p)


def find_primitive_polynomial(p: int, k: int) -> Poly:
    """
    Trouve un polynôme primitif de degré k sur F_p.
    Irréductible et une racine alpha engendre F_{p^k}^×.
    """
    if k < 1:
        raise ValueError("k doit être >= 1")
    if k == 1:
        return Poly.monomial(p, 1, 1)

    order = p ** k - 1
    prime_divs = prime_factors(order)

    for candidate in iter_monic_polynomials(p, k):
        if not is_irreducible(candidate):
            continue
        alpha = Poly.x(p)
        is_primitive = True
        for r in prime_divs:
            if alpha.powmod(order // r, candidate).is_one():
                is_primitive = False
                break
        if is_primitive:
            return candidate

    raise RuntimeError(f"aucun polynôme primitif trouvé pour F_{p}^{k}")
