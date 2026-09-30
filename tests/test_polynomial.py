"""Tests pour l'arithmétique polynomiale et Berlekamp."""

import pytest

from dobble.polynomial import (
    Poly,
    berlekamp_factor,
    factor,
    find_primitive_polynomial,
    is_irreducible,
    is_square_free,
)


class TestPolyArithmetic:
    def test_add_sub_mul_mod(self) -> None:
        p = 5
        a = Poly((1, 2, 3), p)  # 1 + 2X + 3X^2
        b = Poly((4, 0, 1), p)  # 4 + X^2
        assert (a + b).coeffs == (0, 2, 4)
        assert (a - b).coeffs == (2, 2, 2)
        assert (a * b).degree() == 4

        f = Poly((2, 0, 1), p)  # 2 + X^2
        g = Poly((1, 1), p)  # 1 + X
        q, r = (Poly((3, 2, 1), p)).divmod(f)
        assert r.degree() < f.degree()

    def test_powmod(self) -> None:
        p = 7
        x = Poly.x(p)
        f = Poly((3, 0, 1), p)  # 3 + X^2
        assert x.powmod(2, f).coeffs == (0, 0, 1) or x.powmod(2, f).degree() <= 1
        # X^2 mod (3+X^2) = -3 = 4 mod 7
        assert x.powmod(2, f).coeffs == (4,)


class TestIrreducibility:
    @pytest.mark.parametrize(
        "coeffs, p, expected",
        [
            ((1, 1), 2, True),  # X + 1
            ((1, 1, 1), 2, True),  # X^2 + X + 1 over F_2
            ((1, 0, 1), 3, True),  # X^2 + 1 over F_3
            ((1, 0, 1), 2, False),  # X^2 + 1 = (X+1)^2 over F_2
            ((2, 2, 1), 3, True),  # X^2 + 2X + 2, irréductible et primitif sur F_3
        ],
    )
    def test_is_irreducible(self, coeffs: tuple, p: int, expected: bool) -> None:
        f = Poly(coeffs, p)
        assert is_irreducible(f) == expected


class TestBerlekamp:
    def test_factor_x2_minus_1_over_f3(self) -> None:
        # X^2 - 1 = (X-1)(X+1) over F_3
        f = Poly((2, 0, 1), 3)  # 2 + X^2 = X^2 - 1
        factors = factor(f)
        assert len(factors) == 2
        assert all(g.degree() == 1 for g in factors)

    def test_irreducible_stays_single_factor(self) -> None:
        f = Poly((1, 1, 1), 2)  # X^2 + X + 1
        factors = berlekamp_factor(f)
        assert len(factors) == 1
        assert factors[0].degree() == 2


class TestPrimitivePolynomial:
    @pytest.mark.parametrize(
        "p, k, expected_coeffs",
        [
            (2, 1, (0, 1)),
            (2, 2, (1, 1, 1)),
            (2, 3, (1, 1, 0, 1)),
            (3, 2, (2, 1, 1)),  # premier primitif trouvé : X^2 + X + 2 sur F_3
            (2, 4, (1, 1, 0, 0, 1)),
        ],
    )
    def test_known_primitive_polynomials(
        self, p: int, k: int, expected_coeffs: tuple
    ) -> None:
        f = find_primitive_polynomial(p, k)
        assert f.coeffs == expected_coeffs
        assert is_irreducible(f)
