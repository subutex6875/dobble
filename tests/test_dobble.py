"""Tests pour PG(2, q) et la génération Dobble (étape 4)."""

import pytest

from dobble.dobble_game import generate_dobble, generate_dobble_from_field, verify_dobble_invariants
from dobble.finite_field import FiniteField
from dobble.projective_plane import ProjectivePlane


@pytest.mark.parametrize(
    "q, expected_symbols, expected_per_card",
    [
        (2, 7, 3),
        (3, 13, 4),
        (4, 21, 5),
        (5, 31, 6),
        (7, 57, 8),
    ],
)
def test_projective_plane_sizes(q: int, expected_symbols: int, expected_per_card: int) -> None:
    field = _field_for_q(q)
    plane = ProjectivePlane(field)
    assert plane.num_points == expected_symbols
    assert len(plane.lines) == expected_symbols
    assert plane.points_per_line == expected_per_card
    assert plane.verify_plane_axioms()


@pytest.mark.parametrize("q", [2, 3, 4, 5, 7])
def test_generate_dobble_from_field(q: int) -> None:
    field = _field_for_q(q)
    game = generate_dobble_from_field(field)
    assert game.q == q
    assert game.symbols_count == q * q + q + 1
    assert game.symbols_per_card == q + 1
    assert len(game.cards) == game.symbols_count
    assert game.verify()


@pytest.mark.parametrize(
    "n, expected_q",
    [
        (1, 2),
        (7, 7),
        (8, 8),
        (57, 57),
        (58, 61),
    ],
)
def test_generate_dobble_minimal(n: int, expected_q: int) -> None:
    game = generate_dobble(n)
    assert game.q == expected_q
    assert game.symbols_count >= n
    assert game.verify()


def test_classic_dobble_57() -> None:
    """Jeu commercial Dobble : q = 7, 57 cartes, 8 symboles par carte."""
    game = generate_dobble(57)
    assert game.q == 7
    assert game.symbols_count == 57
    assert game.symbols_per_card == 8
    assert len(game.cards) == 57
    assert game.verify()


def test_card_intersection_unique() -> None:
    game = generate_dobble(7)
    cards = game.card_symbol_sets()
    for i in range(len(cards)):
        for j in range(i + 1, len(cards)):
            assert len(cards[i] & cards[j]) == 1


def test_symbol_pair_unique_card() -> None:
    game = generate_dobble(4)
    cards = game.card_symbol_sets()
    pairs: set[tuple[int, int]] = set()
    for card in cards:
        syms = sorted(card)
        for i in range(len(syms)):
            for j in range(i + 1, len(syms)):
                pair = (syms[i], syms[j])
                assert pair not in pairs
                pairs.add(pair)
    n = game.symbols_count
    assert len(pairs) == n * (n - 1) // 2


def _field_for_q(q: int) -> FiniteField:
    """Retourne un corps de taille q pour les petits q des tests."""
    known = {
        2: (2, 1),
        3: (3, 1),
        4: (2, 2),
        5: (5, 1),
        7: (7, 1),
        8: (2, 3),
        9: (3, 2),
    }
    if q not in known:
        raise ValueError(f"q={q} non supporté dans les tests")
    p, k = known[q]
    return FiniteField(p, k)

# test_classic_dobble_57()

# Jeu complet pour au moins 57 cartes (Dobble classique : q=7)
game = generate_dobble(3)
print(game.q, game.symbols_count, game.symbols_per_card)  # 7, 57, 8
print(game.verify())  # True
for card in game.cards:
    print(card)

"""
# Corps fini sous-jacent
params = find_minimal_field(10)  # p=11, k=1, q=11
field = FiniteField(params.p, params.k)
for element in field.enumerate_elements():
    print(field.to_coeffs(element))  # [a0, a1, ...]
"""
