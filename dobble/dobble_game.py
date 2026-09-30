"""Étape 4 : génération du jeu Dobble complet via PG(2, q)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from dobble.field_finder import FieldParams, find_minimal_field
from dobble.finite_field import FiniteField
from dobble.projective_plane import ProjectivePlane


@dataclass(frozen=True)
class DobbleGame:
    """
    Jeu Dobble complet construit sur PG(2, q).

    symbols : q² + q + 1 symboles (points)
    cards   : q² + q + 1 cartes (droites), chacune avec q + 1 symboles
    """

    params: FieldParams
    symbols_count: int
    symbols_per_card: int
    cards: Tuple[Tuple[int, ...], ...]

    @property
    def p(self) -> int:
        return self.params.p

    @property
    def k(self) -> int:
        return self.params.k

    @property
    def q(self) -> int:
        return self.params.q

    def verify(self) -> bool:
        """Vérifie les quatre invariants Dobble."""
        return verify_dobble_invariants(self.cards, self.symbols_count)

    def card_symbol_sets(self) -> List[frozenset[int]]:
        return [frozenset(card) for card in self.cards]

    def __repr__(self) -> str:
        return (
            f"DobbleGame(q={self.q}, symbols={self.symbols_count}, "
            f"cards={len(self.cards)}, per_card={self.symbols_per_card})"
        )


def verify_dobble_invariants(
    cards: Tuple[Tuple[int, ...], ...] | List[frozenset[int]],
    num_symbols: int,
) -> bool:
    """
    Vérifie les invariants d'un jeu Dobble complet :

    1. chaque carte contient exactement q + 1 symboles ;
    2. chaque symbole apparaît exactement q + 1 fois ;
    3. deux cartes distinctes partagent exactement un symbole ;
    4. deux symboles distincts apparaissent ensemble sur exactement une carte.
    """
    card_sets: List[frozenset[int]] = [frozenset(c) for c in cards]
    n_cards = len(card_sets)

    if n_cards == 0 or num_symbols == 0:
        return False

    symbols_per_card = len(card_sets[0])
    q = symbols_per_card - 1
    expected_symbols = q * q + q + 1

    if num_symbols != expected_symbols or n_cards != expected_symbols:
        return False

    # 1. Taille uniforme des cartes
    if any(len(card) != symbols_per_card for card in card_sets):
        return False

    # 2. Fréquence des symboles
    symbol_count: Dict[int, int] = {}
    for card in card_sets:
        for sym in card:
            if sym < 0 or sym >= num_symbols:
                return False
            symbol_count[sym] = symbol_count.get(sym, 0) + 1
    if len(symbol_count) != num_symbols:
        return False
    if any(count != symbols_per_card for count in symbol_count.values()):
        return False

    # 3. Intersection de deux cartes
    for i in range(n_cards):
        for j in range(i + 1, n_cards):
            if len(card_sets[i] & card_sets[j]) != 1:
                return False

    # 4. Deux symboles distincts : une carte commune
    pair_card: Dict[Tuple[int, int], int] = {}
    for card_idx, card in enumerate(card_sets):
        syms = sorted(card)
        for i in range(len(syms)):
            for j in range(i + 1, len(syms)):
                pair = (syms[i], syms[j])
                if pair in pair_card:
                    return False
                pair_card[pair] = card_idx

    expected_pairs = num_symbols * (num_symbols - 1) // 2
    if len(pair_card) != expected_pairs:
        return False

    return True


def generate_dobble(
    n: int,
    field: Optional[FiniteField] = None,
) -> DobbleGame:
    """
    Génère un jeu Dobble complet pour au moins n symboles/cartes.

    Choisit le plus petit F_q avec q >= n, construit PG(2, q), retourne le jeu.
    """
    params = find_minimal_field(n)
    ff = field if field is not None else FiniteField(params.p, params.k)
    if ff.q != params.q:
        raise ValueError("corps fourni incompatible avec les paramètres minimaux")

    plane = ProjectivePlane(ff)
    if not plane.verify_plane_axioms():
        raise RuntimeError("PG(2, q) invalide")

    cards = tuple(tuple(sorted(card)) for card in plane.cards())
    game = DobbleGame(
        params=params,
        symbols_count=plane.num_points,
        symbols_per_card=plane.points_per_line,
        cards=cards,
    )
    if not game.verify():
        raise RuntimeError("invariants Dobble non satisfaits")
    return game


def generate_dobble_from_field(field: FiniteField) -> DobbleGame:
    """Génère un Dobble directement à partir d'un corps fini F_q."""
    params = FieldParams(p=field.p, k=field.k, q=field.q)
    return generate_dobble(field.q, field=field)
