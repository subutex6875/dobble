"""Générateur de jeu Dobble via corps finis et géométrie projective."""

from dobble.dobble_game import DobbleGame, generate_dobble, generate_dobble_from_field
from dobble.field_finder import FieldParams, find_minimal_field
from dobble.finite_field import FieldElement, FiniteField
from dobble.polynomial import Poly, find_primitive_polynomial
from dobble.projective_plane import ProjectivePlane

__all__ = [
    "DobbleGame",
    "FieldElement",
    "FieldParams",
    "FiniteField",
    "Poly",
    "ProjectivePlane",
    "find_minimal_field",
    "find_primitive_polynomial",
    "generate_dobble",
    "generate_dobble_from_field",
]
