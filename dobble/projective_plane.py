"""Plan projectif PG(2, q) sur le corps fini F_q."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, List, Sequence, Tuple

from dobble.finite_field import FieldElement, FiniteField

ProjectiveTriple = Tuple[FieldElement, FieldElement, FieldElement]


@dataclass(frozen=True)
class ProjectivePoint:
    """Point de PG(2, q) : classe projective [x : y : z] ≠ [0 : 0 : 0]."""

    coords: ProjectiveTriple


@dataclass(frozen=True)
class ProjectiveLine:
    """
    Droite projective de PG(2, q).

    Représentée par un vecteur normal (a, b, c) : la droite est l'ensemble des
    points p tels que a·x + b·y + c·z = 0.
    """

    normal: ProjectiveTriple


class ProjectivePlane:
    """
    Plan projectif PG(2, q).

    q = p^k : nombre d'éléments du corps sous-jacent.
    Dimension projective : 2 (ne pas confondre avec k, degré de l'extension).
    """

    def __init__(self, field: FiniteField) -> None:
        self.field = field
        self.q = field.q
        self.p = field.p
        self.k = field.k
        self.num_points = self.q ** 2 + self.q + 1
        self.points_per_line = self.q + 1

        self._points = self._build_points()
        self._lines = self._build_lines()
        self._point_index = {self._point_key(p): i for i, p in enumerate(self._points)}
        self._line_index = {self._line_key(l): i for i, l in enumerate(self._lines)}

        self._line_points: List[frozenset[int]] = []
        for line in self._lines:
            indices = frozenset(
                self._point_index[self._point_key(p)]
                for p in self._points
                if self._incident(line, p)
            )
            self._line_points.append(indices)

        self._point_lines: List[frozenset[int]] = [frozenset() for _ in self._points]
        for line_idx, point_indices in enumerate(self._line_points):
            for point_idx in point_indices:
                self._point_lines[point_idx] = frozenset(
                    self._point_lines[point_idx] | {line_idx}
                )

    @staticmethod
    def _point_key(point: ProjectivePoint) -> Tuple[int, int, int]:
        x, y, z = point.coords
        return (x.coeffs, y.coeffs, z.coeffs)  # type: ignore[return-value]

    @staticmethod
    def _line_key(line: ProjectiveLine) -> Tuple[int, int, int]:
        a, b, c = line.normal
        return (a.coeffs, b.coeffs, c.coeffs)  # type: ignore[return-value]

    def _normalize(self, x: FieldElement, y: FieldElement, z: FieldElement) -> ProjectiveTriple:
        """Représentant canonique unique de [x : y : z]."""
        f = self.field
        if not f.is_zero(x):
            inv = f.inv(x)
            return (f.one, f.mul(y, inv), f.mul(z, inv))
        if not f.is_zero(y):
            inv = f.inv(y)
            return (f.zero, f.one, f.mul(z, inv))
        return (f.zero, f.zero, f.one)

    def _build_points(self) -> List[ProjectivePoint]:
        f = self.field
        points: List[ProjectivePoint] = []
        for a in f.enumerate_elements():
            for b in f.enumerate_elements():
                points.append(ProjectivePoint(self._normalize(f.one, a, b)))
        for c in f.enumerate_elements():
            points.append(ProjectivePoint(self._normalize(f.zero, f.one, c)))
        points.append(ProjectivePoint(self._normalize(f.zero, f.zero, f.one)))
        return points

    def _build_lines(self) -> List[ProjectiveLine]:
        f = self.field
        lines: List[ProjectiveLine] = []
        for a in f.enumerate_elements():
            for b in f.enumerate_elements():
                lines.append(ProjectiveLine(self._normalize(f.one, a, b)))
        for c in f.enumerate_elements():
            lines.append(ProjectiveLine(self._normalize(f.zero, f.one, c)))
        lines.append(ProjectiveLine(self._normalize(f.zero, f.zero, f.one)))
        return lines

    def _dot(self, normal: ProjectiveTriple, point: ProjectiveTriple) -> FieldElement:
        f = self.field
        a, b, c = normal
        x, y, z = point
        return f.add(f.add(f.mul(a, x), f.mul(b, y)), f.mul(c, z))

    def _incident(self, line: ProjectiveLine, point: ProjectivePoint) -> bool:
        return self.field.is_zero(self._dot(line.normal, point.coords))

    @property
    def points(self) -> Sequence[ProjectivePoint]:
        return self._points

    @property
    def lines(self) -> Sequence[ProjectiveLine]:
        return self._lines

    def point_coords(self, index: int) -> List[int]:
        """Coordonnées du point sous forme de triples d'entiers dans F_p."""
        p = self._points[index]
        f = self.field
        return [f.to_coeffs(c) for c in p.coords]

    def line_normal(self, index: int) -> List[int]:
        """Vecteur normal de la droite sous forme de triples d'entiers dans F_p."""
        l = self._lines[index]
        f = self.field
        return [f.to_coeffs(c) for c in l.normal]

    def points_on_line(self, line_index: int) -> frozenset[int]:
        return self._line_points[line_index]

    def lines_through_point(self, point_index: int) -> frozenset[int]:
        return self._point_lines[point_index]

    def cards(self) -> List[frozenset[int]]:
        """Cartes Dobble : chaque droite est une carte, ses points sont les symboles."""
        return list(self._line_points)

    def symbols(self) -> range:
        """Indices des symboles (points)."""
        return range(self.num_points)

    def verify_plane_axioms(self) -> bool:
        """Vérifie les axiomes du plan projectif fini PG(2, q)."""
        n = self.num_points
        q1 = self.points_per_line

        if len(self._points) != n or len(self._lines) != n:
            return False

        # Chaque droite contient exactement q + 1 points
        for pts in self._line_points:
            if len(pts) != q1:
                return False

        # Chaque point appartient à exactement q + 1 droites
        for lines in self._point_lines:
            if len(lines) != q1:
                return False

        # Deux droites distinctes : intersection d'un seul point
        for i in range(n):
            for j in range(i + 1, n):
                if len(self._line_points[i] & self._line_points[j]) != 1:
                    return False

        # Deux points distincts : une seule droite commune
        for i in range(n):
            for j in range(i + 1, n):
                if len(self._point_lines[i] & self._point_lines[j]) != 1:
                    return False

        return True

    def __repr__(self) -> str:
        return f"ProjectivePlane(PG(2, {self.q}), p={self.p}, k={self.k})"
