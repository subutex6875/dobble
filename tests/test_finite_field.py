"""Tests pour les corps finis F_{p^k} (étapes 2 et 3)."""

import pytest

from dobble.finite_field import FiniteField
from dobble.polynomial import Poly, is_irreducible


def _make_field(p: int, k: int, primitive_coeffs: tuple | None = None) -> FiniteField:
    if primitive_coeffs is not None:
        poly = Poly(primitive_coeffs, p)
        return FiniteField(p, k, primitive_poly=poly)
    return FiniteField(p, k)


class TestFieldSize:
    @pytest.mark.parametrize(
        "p, k, q",
        [
            (2, 1, 2),
            (2, 2, 4),
            (2, 3, 8),
            (3, 2, 9),
            (2, 4, 16),
        ],
    )
    def test_element_count(self, p: int, k: int, q: int) -> None:
        field = _make_field(p, k)
        elements = list(field.enumerate_elements())
        assert len(elements) == q
        assert field.q == q


class TestKnownFields:
    @pytest.mark.parametrize(
        "p, k, primitive_coeffs",
        [
            (2, 1, (0, 1)),
            (2, 2, (1, 1, 1)),
            (2, 3, (1, 1, 0, 1)),
            (3, 2, (2, 2, 1)),
            (2, 4, (1, 1, 0, 0, 1)),
        ],
    )
    def test_irreducible_and_primitive(
        self, p: int, k: int, primitive_coeffs: tuple
    ) -> None:
        field = _make_field(p, k, primitive_coeffs)
        assert field.modulus.coeffs == primitive_coeffs
        assert is_irreducible(field.modulus)
        assert field.is_irreducible_modulus()
        assert field.is_primitive_modulus()

    @pytest.mark.parametrize(
        "p, k, primitive_coeffs",
        [
            (2, 1, (0, 1)),
            (2, 2, (1, 1, 1)),
            (2, 3, (1, 1, 0, 1)),
            (3, 2, (2, 2, 1)),
            (2, 4, (1, 1, 0, 0, 1)),
        ],
    )
    def test_field_axioms(self, p: int, k: int, primitive_coeffs: tuple) -> None:
        field = _make_field(p, k, primitive_coeffs)
        assert field.verify_field_axioms()


class TestRepresentation:
    def test_f4_coefficients(self) -> None:
        field = _make_field(2, 2, (1, 1, 1))
        # Énumération : 0, 1, alpha, alpha+1
        elems = list(field.enumerate_elements())
        coeffs = [field.to_coeffs(e) for e in elems]
        assert len(coeffs) == 4
        assert coeffs[0] == [0, 0]
        assert coeffs[1] == [1, 0]

        # alpha = X mod f = deuxième élément base : (0, 1)
        alpha = field.from_coeffs([0, 1])
        assert field.mul(alpha, alpha) == field.from_coeffs([1, 1])  # alpha^2 = alpha + 1

    def test_from_int_roundtrip(self) -> None:
        field = _make_field(3, 2, (2, 2, 1))
        for i in range(field.q):
            e = field.from_int(i)
            assert field.to_int(e) == i
            assert field.to_coeffs(e) == [(i // (3**j)) % 3 for j in range(2)]

    def test_f9_multiplication_table(self) -> None:
        field = _make_field(3, 2, (2, 2, 1))
        a = field.from_coeffs([0, 1])  # alpha
        # alpha^2 = -2*alpha - 2 = alpha + 1 sur F_3
        assert field.pow(a, 2) == field.from_coeffs([1, 1])
        assert field.pow(a, 8) == field.one


class TestOperations:
    def test_f8_arithmetic(self) -> None:
        field = _make_field(2, 3, (1, 1, 0, 1))
        a = field.from_coeffs([0, 0, 1])
        b = field.from_coeffs([1, 1, 0])
        assert field.add(a, b) == field.from_coeffs([1, 1, 1])
        assert field.sub(b, a) == field.from_coeffs([1, 1, 1])
        c = field.mul(a, b)
        inv_c = field.inv(c)
        assert field.mul(c, inv_c) == field.one

    def test_auto_primitive_polynomial(self) -> None:
        for p, k in [(2, 2), (2, 3), (3, 2), (2, 4)]:
            field = FiniteField(p, k)
            assert field.is_irreducible_modulus()
            assert field.is_primitive_modulus()
            assert field.verify_field_axioms()
