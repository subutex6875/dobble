"""Étape 1 : trouver le plus petit corps fini F_q avec q >= N."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, List, Tuple


@dataclass(frozen=True)
class FieldParams:
    """
    Paramètres du corps fini minimal.

    p : caractéristique (nombre premier)
    k : degré de l'extension F_{p^k} / F_p
    q : nombre d'éléments, q = p^k
    """

    p: int
    k: int
    q: int


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def sieve_primes(limit: int) -> List[int]:
    """Crible d'Ératosthène jusqu'à limit (inclus)."""
    if limit < 2:
        return []
    is_composite = [False] * (limit + 1)
    primes: List[int] = []
    for n in range(2, limit + 1):
        if not is_composite[n]:
            primes.append(n)
            if n * n <= limit:
                step = n
                start = n * n
                is_composite[start : limit + 1 : step] = [True] * len(
                    range(start, limit + 1, step)
                )
    return primes


def iter_primes() -> Iterator[int]:
    """Génère les nombres premiers sans limite a priori (trial division)."""
    yield 2
    n = 3
    primes: List[int] = [2]
    while True:
        if is_prime(n):
            primes.append(n)
            yield n
        n += 2


def minimal_exponent(p: int, n: int) -> int:
    """Plus petit k >= 1 tel que p^k >= n."""
    if n <= 1:
        return 1
    k = 1
    q = p
    while q < n:
        k += 1
        q *= p
    return k


def find_minimal_field(n: int) -> FieldParams:
    """
    Trouve le plus petit corps fini F_q avec q >= n.

    Minimise q = p^k sous la contrainte q >= n, avec p premier et k >= 1.
    """
    if n <= 1:
        return FieldParams(p=2, k=1, q=2)

    best: FieldParams | None = None

    for p in iter_primes():
        k = minimal_exponent(p, n)
        q = p ** k
        candidate = FieldParams(p=p, k=k, q=q)

        if best is None or candidate.q < best.q:
            best = candidate

        # Pour p > best.q, k=1 donne q = p > best.q : impossible d'améliorer.
        if best is not None and p > best.q:
            break

    assert best is not None
    return best
