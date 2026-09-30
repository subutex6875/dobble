"""Étape 2-3 : corps fini F_{p^k} avec représentation polynomiale."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, List, Optional, Tuple

from dobble.polynomial import Poly, find_primitive_polynomial, is_irreducible, prime_factors


@dataclass(frozen=True)
class FieldElement:
    """Élément de F_{p^k} : vecteur (a_0, ..., a_{k-1}) avec a_i dans F_p."""

    coeffs: Tuple[int, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "coeffs", tuple(self.coeffs))


class FiniteField:
    """
    Corps fini F_{p^k}.

    p : caractéristique (premier)
    k : degré de l'extension
    q : p^k, nombre d'éléments
    """

    def __init__(
        self,
        p: int,
        k: int,
        primitive_poly: Optional[Poly] = None,
    ) -> None:
        if not isinstance(p, int) or p < 2:
            raise ValueError("p doit être un entier >= 2")
        if not isinstance(k, int) or k < 1:
            raise ValueError("k doit être un entier >= 1")

        self.p = p
        self.k = k
        self.q = p ** k

        if k == 1:
            self._modulus = Poly.monomial(p, 1, 1)
        elif primitive_poly is not None:
            if primitive_poly.p != p or primitive_poly.degree() != k:
                raise ValueError("polynôme primitif incompatible")
            inv = pow(primitive_poly.leading_coeff(), -1, p)
            self._modulus = primitive_poly.scalar_mul(inv)
        else:
            self._modulus = find_primitive_polynomial(p, k)

        self._zero = FieldElement(tuple([0] * k))
        self._one = FieldElement(tuple([1] + [0] * (k - 1)))

    @property
    def modulus(self) -> Poly:
        """Polynôme primitif f tel que F_{p^k} ≅ F_p[X]/(f)."""
        return self._modulus

    @property
    def zero(self) -> FieldElement:
        return self._zero

    @property
    def one(self) -> FieldElement:
        return self._one

    def from_coeffs(self, coeffs: List[int] | Tuple[int, ...]) -> FieldElement:
        """Construit un élément à partir de (a_0, ..., a_{k-1})."""
        normalized = [(c % self.p) for c in coeffs]
        if len(normalized) < self.k:
            normalized.extend([0] * (self.k - len(normalized)))
        elif len(normalized) > self.k:
            raise ValueError(f"attendu {self.k} coefficients, reçu {len(normalized)}")
        return FieldElement(tuple(normalized))

    def to_coeffs(self, element: FieldElement) -> List[int]:
        """Retourne [a_0, ..., a_{k-1}] pour l'élément donné."""
        return [c % self.p for c in element.coeffs]

    def from_int(self, n: int) -> FieldElement:
        """Encode n in [0, q) en base p sur k digits (a_0 + a_1 p + ...)."""
        n = n % self.q
        digits = []
        for _ in range(self.k):
            digits.append(n % self.p)
            n //= self.p
        return FieldElement(tuple(digits))

    def to_int(self, element: FieldElement) -> int:
        """Décodage base p de l'élément."""
        result = 0
        power = 1
        for c in element.coeffs:
            result += (c % self.p) * power
            power *= self.p
        return result

    def enumerate_elements(self) -> Iterator[FieldElement]:
        """Énumère les q éléments du corps."""
        for i in range(self.q):
            yield self.from_int(i)

    def add(self, a: FieldElement, b: FieldElement) -> FieldElement:
        coeffs = [(a.coeffs[i] + b.coeffs[i]) % self.p for i in range(self.k)]
        return FieldElement(tuple(coeffs))

    def sub(self, a: FieldElement, b: FieldElement) -> FieldElement:
        coeffs = [(a.coeffs[i] - b.coeffs[i]) % self.p for i in range(self.k)]
        return FieldElement(tuple(coeffs))

    def neg(self, a: FieldElement) -> FieldElement:
        return FieldElement(tuple((-c) % self.p for c in a.coeffs))

    def _to_poly(self, a: FieldElement) -> Poly:
        return Poly(a.coeffs, self.p)

    def _from_poly(self, poly: Poly) -> FieldElement:
        coeffs = list(poly.coeffs)
        if len(coeffs) < self.k:
            coeffs.extend([0] * (self.k - len(coeffs)))
        return FieldElement(tuple(coeffs[: self.k]))

    def mul(self, a: FieldElement, b: FieldElement) -> FieldElement:
        pa = self._to_poly(a)
        pb = self._to_poly(b)
        product = (pa * pb) % self._modulus
        return self._from_poly(product)

    def eq(self, a: FieldElement, b: FieldElement) -> bool:
        return a.coeffs == b.coeffs

    def is_zero(self, a: FieldElement) -> bool:
        return all(c % self.p == 0 for c in a.coeffs)

    def pow(self, a: FieldElement, exponent: int) -> FieldElement:
        if exponent < 0:
            a = self.inv(a)
            exponent = -exponent
        if exponent == 0:
            return self.one
        pa = self._to_poly(a)
        result = pa.powmod(exponent, self._modulus)
        return self._from_poly(result)

    def inv(self, a: FieldElement) -> FieldElement:
        if self.is_zero(a):
            raise ZeroDivisionError("inverse de zéro")
        # a^{-1} = a^{q-2} (Fermat dans le corps fini)
        return self.pow(a, self.q - 2)

    def div(self, a: FieldElement, b: FieldElement) -> FieldElement:
        return self.mul(a, self.inv(b))

    def is_irreducible_modulus(self) -> bool:
        return is_irreducible(self._modulus)

    def is_primitive_modulus(self) -> bool:
        """Vérifie que X mod f engendre le groupe multiplicatif."""
        if self.k == 1:
            return True
        order = self.q - 1
        alpha = Poly.x(self.p)
        for r in prime_factors(order):
            if alpha.powmod(order // r, self._modulus).is_one():
                return False
        return True

    def verify_field_axioms(self) -> bool:
        """Vérification exhaustive pour petits corps."""
        elements = list(self.enumerate_elements())
        if len(elements) != self.q:
            return False

        # Unicité
        ints = {self.to_int(e) for e in elements}
        if len(ints) != self.q:
            return False

        # Additif : commutativité, neutre, opposé
        for a in elements:
            if not self.eq(self.add(a, self.zero), a):
                return False
            if not self.eq(self.add(self.zero, a), a):
                return False
            found_neg = False
            for b in elements:
                if self.eq(self.add(a, b), self.zero):
                    found_neg = True
                    break
            if not found_neg:
                return False

        # Multiplicatif : neutre, inverse pour non-nuls
        for a in elements:
            if not self.eq(self.mul(a, self.one), a):
                return False
            if self.is_zero(a):
                continue
            found_inv = False
            for b in elements:
                if self.eq(self.mul(a, b), self.one):
                    found_inv = True
                    break
            if not found_inv:
                return False

        # Distributivité (échantillon exhaustif si q petit)
        for a in elements:
            for b in elements:
                for c in elements:
                    left = self.mul(a, self.add(b, c))
                    right = self.add(self.mul(a, b), self.mul(a, c))
                    if not self.eq(left, right):
                        return False

        return True

    def __repr__(self) -> str:
        return f"FiniteField(F_{self.p}^{self.k}), q={self.q}, f={self._modulus})"
