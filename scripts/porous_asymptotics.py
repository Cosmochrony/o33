#!/usr/bin/env python3
"""Exact limiting law for porous depths in the O33 nine-centre filtration.

Archived reproducibility script for the exact limit law of Section ``Asymptotic
Distribution of Porous Depths'' of the O33 paper (labels thm:porous-convergence and
thm:porous-law).  It certifies the 240-triangle rational polygon decomposition, the count
law (3/20, 7/24, 9/40, 5/24, 1/8), the piecewise affine intensity, the total mass 28/15,
and the conditional mean 491/10080.  Pure standard library (fractions, itertools,
dataclasses); every step is exact rational arithmetic; runtime is well under one second.

For uniformly sampled normalised parameters (u,v) in F_q^2, away from the negligible
collision locus, write the four positive centre distances as the sorted values of

    ||u||, ||v||, ||u+v||, ||u-v||.

The four paired cyclic spacings are the successive gaps from zero.  A spacing d gives a
porous depth exactly when d is even, at n=(d-4)/2.  The script partitions the continuous
two-torus into exact rational polygons, retains the two-scale parity weights, and computes
the piecewise-linear density of n/q in the q -> infinity limit.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations, product


Point = tuple[F, F]
Polygon = list[Point]


@dataclass(frozen=True)
class Affine:
    x: F
    y: F
    c: F

    def __call__(self, point: Point) -> F:
        return self.x * point[0] + self.y * point[1] + self.c

    def __sub__(self, other: "Affine") -> "Affine":
        return Affine(self.x - other.x, self.y - other.y, self.c - other.c)


@dataclass(frozen=True)
class Branch:
    distance: Affine
    inequalities: tuple[Affine, ...]
    parity_vector: tuple[int, int]
    parity_constant: int


def clip(poly: Polygon, inequality: Affine) -> Polygon:
    """Clip a convex polygon by inequality(point) >= 0 using exact arithmetic."""
    if not poly:
        return []
    result: Polygon = []
    for left, right in zip(poly, poly[1:] + poly[:1]):
        fl = inequality(left)
        fr = inequality(right)
        in_left = fl >= 0
        in_right = fr >= 0
        if in_left:
            result.append(left)
        if in_left != in_right:
            ratio = fl / (fl - fr)
            result.append(
                (
                    left[0] + ratio * (right[0] - left[0]),
                    left[1] + ratio * (right[1] - left[1]),
                )
            )
    cleaned: Polygon = []
    for point in result:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def twice_area(poly: Polygon) -> F:
    return abs(sum(
        left[0] * right[1] - left[1] * right[0]
        for left, right in zip(poly, poly[1:] + poly[:1])
    ))


def form_branches(ax: int, ay: int, k_values: range) -> list[Branch]:
    branches = []
    parity_vector = (ax % 2, ay % 2)
    raw = Affine(F(ax), F(ay), F(0))
    for k in k_values:
        centred = Affine(raw.x, raw.y, -F(k))
        for sign in (-1, 1):
            distance = Affine(sign * centred.x, sign * centred.y, sign * centred.c)
            branches.append(
                Branch(
                    distance=distance,
                    inequalities=(distance, Affine(-distance.x, -distance.y, F(1, 2) - distance.c)),
                    parity_vector=parity_vector,
                    parity_constant=k % 2,
                )
            )
    return branches


FORMS = (
    form_branches(1, 0, range(0, 2)),
    form_branches(0, 1, range(0, 2)),
    form_branches(1, 1, range(0, 3)),
    form_branches(1, -1, range(-1, 2)),
)


def parity_weight(left: Branch | None, right: Branch) -> F:
    """Probability that the first distance or a distance difference is even."""
    if left is None:
        return F(1, 2)
    vector = (
        left.parity_vector[0] ^ right.parity_vector[0],
        left.parity_vector[1] ^ right.parity_vector[1],
    )
    if vector != (0, 0):
        return F(1, 2)
    return F(int(left.parity_constant == right.parity_constant))


def cells() -> list[tuple[Polygon, tuple[Branch, ...]]]:
    square: Polygon = [(F(0), F(0)), (F(1), F(0)), (F(1), F(1)), (F(0), F(1))]
    result = []
    for branches in product(*FORMS):
        base = square
        for branch in branches:
            for inequality in branch.inequalities:
                base = clip(base, inequality)
                if len(base) < 3 or twice_area(base) == 0:
                    break
            if len(base) < 3 or twice_area(base) == 0:
                break
        if len(base) < 3 or twice_area(base) == 0:
            continue
        for order in permutations(range(4)):
            poly = base
            ordered = tuple(branches[index] for index in order)
            for left, right in zip(ordered, ordered[1:]):
                poly = clip(poly, right.distance - left.distance)
                if len(poly) < 3 or twice_area(poly) == 0:
                    break
            if len(poly) >= 3 and twice_area(poly) > 0:
                result.append((poly, ordered))
    return result


def triangle_contributions(
    partition: list[tuple[Polygon, tuple[Branch, ...]]],
) -> list[tuple[F, F, F, F, F]]:
    """Return exact (area, parity weight, minimum, middle, maximum) records."""
    result = []
    for poly, ordered in partition:
        for j in range(1, len(poly) - 1):
            triangle = (poly[0], poly[j], poly[j + 1])
            area = twice_area(list(triangle)) / 2
            previous: Branch | None = None
            for branch in ordered:
                gap = branch.distance if previous is None else branch.distance - previous.distance
                values = sorted(gap(point) for point in triangle)
                weight = parity_weight(previous, branch)
                if weight and values[-1] > values[0]:
                    result.append((area, weight, values[0], values[1], values[2]))
                previous = branch
    return result


def parity_value(branch: Branch, bits: tuple[int, int]) -> int:
    return (
        branch.parity_vector[0] * bits[0]
        + branch.parity_vector[1] * bits[1]
        + branch.parity_constant
    ) % 2


def limiting_count_law(
    partition: list[tuple[Polygon, tuple[Branch, ...]]],
) -> list[F]:
    """Return the limiting probabilities of exactly 0,...,4 porous depths."""
    probabilities = [F(0) for _ in range(5)]
    for poly, ordered in partition:
        area = twice_area(poly) / 2
        for bits in product((0, 1), repeat=2):
            parities = [parity_value(branch, bits) for branch in ordered]
            count = int(parities[0] == 0)
            count += sum(left == right for left, right in zip(parities, parities[1:]))
            probabilities[count] += area / 4
    assert sum(probabilities) == 1
    return probabilities


def density_g(term: tuple[F, F, F, F, F], t: F) -> F:
    area, weight, a, b, c = term
    if t < a or t > c or a == c:
        return F(0)
    if a == b:
        return weight * 2 * area * (c - t) / (c - a) ** 2
    if b == c:
        return weight * 2 * area * (t - a) / (c - a) ** 2
    if t <= b:
        return weight * 2 * area * (t - a) / ((b - a) * (c - a))
    return weight * 2 * area * (c - t) / ((c - a) * (c - b))


def piecewise_density(terms: list[tuple[F, F, F, F, F]]) -> list[tuple[F, F, F, F]]:
    """Return exact segments (left, right, slope, intercept) for y=n/q=g/2."""
    knots = sorted({F(0), F(1, 4)} | {value / 2 for term in terms for value in term[2:]})
    segments = []
    for left, right in zip(knots, knots[1:]):
        if left == right:
            continue
        y1 = (2 * left + right) / 3
        y2 = (left + 2 * right) / 3
        f1 = 2 * sum(density_g(term, 2 * y1) for term in terms)
        f2 = 2 * sum(density_g(term, 2 * y2) for term in terms)
        slope = (f2 - f1) / (y2 - y1)
        intercept = f1 - slope * y1
        if segments and segments[-1][2:] == (slope, intercept):
            segments[-1] = (segments[-1][0], right, slope, intercept)
        else:
            segments.append((left, right, slope, intercept))
    return segments


def integrate_segment(segment: tuple[F, F, F, F], moment: int = 0) -> F:
    left, right, slope, intercept = segment
    if moment == 0:
        return slope * (right**2 - left**2) / 2 + intercept * (right - left)
    if moment == 1:
        return slope * (right**3 - left**3) / 3 + intercept * (right**2 - left**2) / 2
    raise ValueError(moment)


def main() -> None:
    partition = cells()
    area_sum = sum(twice_area(poly) / 2 for poly, _ in partition)
    assert area_sum == 1, area_sum
    terms = triangle_contributions(partition)
    segments = piecewise_density(terms)
    count_law = limiting_count_law(partition)
    mass = sum(integrate_segment(segment) for segment in segments)
    first_moment = sum(integrate_segment(segment, 1) for segment in segments)
    print(f"Triangle contributions: {len(terms)}")
    print("Limiting intensity density for y=n/q:")
    for left, right, slope, intercept in segments:
        print(f"  [{left}, {right}]: f(y) = ({slope}) y + ({intercept})")
    print(f"Total mass (mean porous-depth count): {mass} = {float(mass):.12f}")
    print("Limiting count law P(N=k):")
    for count, probability in enumerate(count_law):
        print(f"  k={count}: {probability} = {float(probability):.12f}")
    print(
        "Mean normalised depth conditional on porosity: "
        f"{first_moment / mass} = {float(first_moment / mass):.12f}"
    )
    assert mass == sum(F(index) * probability for index, probability in enumerate(count_law))
    assert mass == F(28, 15)


if __name__ == "__main__":
    main()
