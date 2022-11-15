# SymbolicWedderburn.jl
[![CI](https://github.com/kalmarek/SymbolicWedderburn.jl/workflows/CI/badge.svg?branch=master)](https://github.com/kalmarek/SymbolicWedderburn.jl/actions)
[![codecov](https://codecov.io/gh/kalmarek/SymbolicWedderburn.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kalmarek/SymbolicWedderburn.jl)

Amazing package providing symbolic but explicit
* decomposition of group representations into (semi-)simple representations and
* the Wedderburn decompositions for endomorphisms of those representations.

We work with
* a (linear) **orthogonal** actions of a **finite group** `G` on **finite dimensional** vector space `V` over `K = ℝ` or `K = ℂ` (i.e. a `KG`-module `V`) and
* linear, `G`-equivariant maps (equivariant endomorphisms) `f : V → V`.

These objects are of primary inportance in the study of the representation theory for finite groups, but also naturally arise from (non)commutative polynomial optimization with group symmetry. The aim of the package is to facilitate such uses.

## A bit of theory
By Maschke's theorem `V` can be decomposed uniquely `V ≅ V₁ ⊕ ⋯ ⊕ Vᵣ` into isotypic/semisimple subspaces `Vᵢ` and each of `Vᵢ ≅ mᵢWᵢ` is (in a non-canonical fashion) isomorphic to a direct sum of `mᵢ` copies of irreducible/simple subspaces `Wᵢ`. By (symbolic) computation in (the group algebra) `KG`, `SymbolicWedderburn` is capable of producing the exact isomorphism `V ≅ V₁ ⊕ ⋯ ⊕ Vᵣ` in the form of a collection of projections `πᵢ : V → Vᵢ` either in the group algebra (lazy, unevaluated form), or in terms of projection matrices, when a basis for `V` is explicitly given.

The isomorphism produces a decomposition of `End_G(V)` (the set of linear `G`-equivariant self-maps of `V`) in the sense of Artin-Wedderburn theorem, i.e. the projections `πᵢ` block-diagonalize `f ≅ f₁ ⊕ ⋯ ⊕ fᵣ` where `fᵢ : Vᵢ → Vᵢ`. In terms of matrices if `f` is given by `n×n`-matrix, then we can rewrite it as a block diagonal matrix with blocks of sizes `nᵢ×nᵢ` where `Σᵢ nᵢ = n = dim V` and each `nᵢ = mᵢ · dim Wᵢ`.

### Semi-definite constraints
For example, if a basis for a semi-definite constraint admits an action of a finite group, then the semi-definite matrix of a **invariant solution** is such an equivariant endomorphism.
In particular if `P` is a positive semidefinite constraint, when searching for an `G`-**invariant** solution we may replace
> `0 ⪯ P[1:n, 1:n]`

by a sequence of constraints

> `0 ⪯ P[1:nᵢ, 1:nᵢ]` for `i = 1…r`,

greatly reducing the computational complexity (The size of psd constraint is reduced from `n²` to `Σᵢ nᵢ²`). Such replacement can be justified if e.g. the objective is symmetric and the set of linear constraints follows a similar group-symmetric structure.

If we are only interested in the _feasibility_ of an optimization problem, then such replacement is always justified (i.e. an _invariant solution_ is a _honest solution_ which might not attain the same objective).

### Rank one projection to simple components

Sometimes even stronger reduction is possible when the acting group `G` is _sufficiently complicated_ and we have a _minimal projection system_ at our disposal (The package tries to compute such system by a heuristic algorithm. If the approach fails please open an issue!). In such case we can (often) find subsequent projections `Vᵢ → K^{mᵢ}` (depending only on the multiplicity of the irreducible, **not** on its dimension!). This leads to an equivalent formulation for the psd constraint with `nᵢ = mᵢ` further reducing its size.

Moreover in the case of symmetric optimization problems it's possible to use the symmetry to reduce the number of linear constraints (since in that case only one constraint **per orbit** is needed). `SymbolicWedderburn` facilitates also this simplification.

## Example

In [_Aut(𝔽₅) has property (T)_](https://arxiv.org/abs/1712.07167) we use the trick above to successfully simplify and solve a large semidefinite problem coming from sum-of-squares optimization.

The original problem had one (symmetric) psd constraint of size `4641×4641` and `11_154_301` linear constraints. By exploiting its (admittedly -- pretty large) symmetry group (of order `3840`) we can reduce this problem to `20` (symmetric) psd constraints of sizes
```
[56  38  34  32  27  27  23  23  22  22  18  17  9  8  6  2  1  1  1  1]
```
which correspond to (the simple) `Wᵢ` blocks above. In particular, the number of variables in psd constraints was reduced from `10_771_761` to just `5_707`.

Moreover, the symmetry group has just `7 229` orbits (when acting on the subspace of linear constraints), so the symmetrized problem has equal number of (a bit denser) linear constraints.

The symmetrized problem is solvable in ~20 minutes on an average office laptop (with `16GB` of RAM).

For more examples you may have a look at [dihedral action example](https://github.com/kalmarek/SymbolicWedderburn.jl/blob/master/examples/ex_robinson_form.jl), or different [sum of squares formulations](https://github.com/kalmarek/SymbolicWedderburn.jl/blob/master/examples/sos_problem.jl).
# Related software

## Sum of Squares optimization
This package is used by [SumOfSquares](https://github.com/jump-dev/SumOfSquares.jl) to perform exactly this reduction, for an example use see its [documentation](https://jump.dev/SumOfSquares.jl/latest/generated/Symmetry/dihedral_symmetry_of_the_robinson_form/).

The software for sum of (hermitian) squares computations in a non-commutative setting (group algebra of a infinite group) using `SymbolicWedderburn` is my project [`PropertyT.jl`](https://github.com/kalmarek/PropertyT.jl/) (unregistered). There we used the sum of squares optimization to prove Property (T) for special automorphisms group of the free group. It's a cool result, [check it out!](https://annals.math.princeton.edu/2021/193-2/p03).

## Other symbolic decompositions
The main aim of `GAP` package [`Wedderga`](https://www.gap-system.org/Manuals/pkg/wedderga/doc/chap0.html) is to

> compute the simple components of the Wedderburn decomposition of semisimple group algebras of finite groups over finite fields and over subfields of finite cyclotomic extensions of the rationals.

The focus is thus on symbolic computations and identifying _isomorphism type_ of the simple components.
`SymbolicWedderburn` makes no efforts to compute the types or defining fields,
it's primary goal is to compute symbolic/numerical Wedderburn-Artin isomorphism in a form usable for (polynomial) optimization. `Wedderga` also contains much more sophisticated methods for computing _a complete set of orthogonal primitive idempotents_ (i.e. a minimal projection system) through Shoda pairs.
In principle those idempotents could be computed using [`Oscar`](https://github.com/oscar-system/Oscar.jl) and used in `SymbolicWedderburn`.

# Citing this package
If you happen to use `SymbolicWedderburn` please cite either of
* M. Kaluba, P.W. Nowak and N. Ozawa *$Aut(F₅)$ has property (T)* [1712.07167](https://arxiv.org/abs/1712.07167), and
* M. Kaluba, D. Kielak and P.W. Nowak *On property (T) for $Aut(Fₙ)$ and $SLₙ(Z)$* [1812.03456](https://arxiv.org/abs/1812.03456).

(Follow the arxiv link for proper link to the journal.)
