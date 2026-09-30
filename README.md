# Dobble — générateur via corps finis

Construction d'un jeu Dobble complet à partir du plan projectif `PG(2, q)` sur le corps fini `F_q = F_{p^k}`.

## Notations

| Symbole | Signification |
|---------|---------------|
| `p` | Caractéristique (nombre premier) |
| `k` | Degré de l'extension `F_{p^k} / F_p` |
| `q = p^k` | Nombre d'éléments du corps |
| `PG(2, q)` | Plan projectif (57 cartes pour `q = 7`) |

**Attention :** le `k` de l'extension n'est pas la dimension du plan projectif (toujours 2).

## Structure

```
dobble/
  polynomial.py       # Polynômes sur F_p, Berlekamp, polynômes primitifs
  field_finder.py     # Étape 1 : plus petit F_q avec q >= N
  finite_field.py     # Étapes 2-3 : corps F_{p^k} et opérations
  projective_plane.py # Étape 4a : plan projectif PG(2, q)
  dobble_game.py      # Étape 4b : génération et vérification Dobble
tests/
```

## Installation et tests

```bash
pip install -e ".[dev]"
pytest
```

## Usage

```python
from dobble import generate_dobble, find_minimal_field, FiniteField

# Jeu complet pour au moins 57 cartes (Dobble classique : q=7)
game = generate_dobble(7)
print(game.q, game.symbols_count, game.symbols_per_card)  # 7, 57, 8
print(game.verify())  # True

# Corps fini sous-jacent
params = find_minimal_field(10)  # p=11, k=1, q=11
field = FiniteField(params.p, params.k)
for element in field.enumerate_elements():
    print(field.to_coeffs(element))  # [a0, a1, ...]
```

## Invariants vérifiés automatiquement

1. Chaque carte contient exactement `q + 1` symboles
2. Chaque symbole apparaît sur exactement `q + 1` cartes
3. Deux cartes distinctes partagent exactement un symbole
4. Deux symboles distincts apparaissent ensemble sur exactement une carte
