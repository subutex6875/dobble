"""Tests pour l'étape 1 : recherche du corps minimal."""

import pytest

from dobble.field_finder import FieldParams, find_minimal_field, is_prime, minimal_exponent


@pytest.mark.parametrize(
    "n, expected",
    [
        (1, FieldParams(p=2, k=1, q=2)),
        (2, FieldParams(p=2, k=1, q=2)),
        (3, FieldParams(p=3, k=1, q=3)),
        (4, FieldParams(p=2, k=2, q=4)),
        (5, FieldParams(p=5, k=1, q=5)),
        (7, FieldParams(p=7, k=1, q=7)),
        (8, FieldParams(p=2, k=3, q=8)),
        (9, FieldParams(p=3, k=2, q=9)),
        (10, FieldParams(p=11, k=1, q=11)),
        (11, FieldParams(p=11, k=1, q=11)),
        (16, FieldParams(p=2, k=4, q=16)),
        (17, FieldParams(p=17, k=1, q=17)),
    ],
)
def test_find_minimal_field(n: int, expected: FieldParams) -> None:
    result = find_minimal_field(n)
    assert result == expected


def test_minimal_exponent() -> None:
    assert minimal_exponent(2, 1) == 1
    assert minimal_exponent(2, 2) == 1
    assert minimal_exponent(2, 3) == 2
    assert minimal_exponent(2, 8) == 3
    assert minimal_exponent(3, 10) == 3


def test_is_prime() -> None:
    assert is_prime(2)
    assert is_prime(97)
    assert not is_prime(1)
    assert not is_prime(9)
    assert not is_prime(100)
