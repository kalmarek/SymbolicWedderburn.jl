# SymbolicWedderburn.jl
[![CI](https://github.com/kalmarek/SymbolicWedderburn.jl/workflows/CI/badge.svg?branch=master)](https://github.com/kalmarek/SymbolicWedderburn.jl/actions)
[![codecov](https://codecov.io/gh/kalmarek/SymbolicWedderburn.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kalmarek/SymbolicWedderburn.jl)

Amazing package providing symbolic Wedderburn decompositions for group *-algebras.
Group *-algebras are closely related to
* linear actions of **finite groups** on **finite dimensional** `k`-vector space `V` (i.e. a `k(G)`-module structure on `V`), which can be expressed in the language of
* representations `ρ : G → GL(V)` of a (finite) group `G` on a (finite dimensional) `k`-vector space `V`.

These objects naturally arise from (non)commutative polynomial optimization with symmetry and the aim of the package is to facilitate such uses.

By Masche's theorem, `V ≅ m₁V₁ ⊕ ⋯ ⊕ mᵣVᵣ ≅ W₁ ⊕ ⋯ ⊕ Wᵣ` can be decomposed uniquely into isotypic/semisimple components `Wᵢ` and nonuniquely into irreducible/simple components `Vᵢ` associated with a group `G`. By (symbolic) computation in `k(G)`, `SymbolicWedderburn.jl` is capable of producing exact
* isomorphism `V ≅ W₁ ⊕ ⋯ ⊕ Wᵣ`, or
* projection `V → V₁ ⊕ ⋯ ⊕ Vᵣ`.

The isomorphism gives a decomposition of `End_G(V)` in the sense of Wedderburn-Artin theorem. The projection is enough to reconstruct any matrix `P ∈ End_G(V)'` in the commutant of the endomorphism algebra.

In case of problems of positive semidefinite optimization which enjoy group symmetry this decomposition can be used to reduce an **invariant** psd constraint

> `0 ⪯ P[1:N, 1:N]`

(where `N = dim(V)` is the dimension of `V`), to a sequence of constraints
> `0 ⪯ P[1:nᵢ, 1:nᵢ]` where `i = 1...r`,

greatly reducing the computational complexity.
When computing decomposition into isotypic/semisimple components `Wᵢ` we have `nᵢ = mᵢ·dim(Vᵢ)`,
thus the size of psd constraint is reduced from `N²` to `Σᵢ nᵢ²` where `N = Σᵢ nᵢ`.
Sometimes even stronger reduction is possible (when the acting group `G` is _sufficiently complicated_ and a heuristic algorithm finding a minimal projection system is successful).
In such case each `nᵢ = mᵢ`.


This package is used by [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl) to perform exactly this reduction, for an example use see its [documentation](https://jump.dev/SumOfSquares.jl/latest/generated/Symmetry/dihedral_symmetry_of_the_robinson_form/).

#### Related software
The main aim of `GAP` package [`Wedderga`](https://www.gap-system.org/Manuals/pkg/wedderga/doc/chap0.html) is to

> compute the simple components of the Wedderburn decomposition of semisimple group algebras of finite groups over finite fields and over subfields of finite cyclotomic extensions of the rationals.

The focus is thus on symbolic computations and identifying _isomorphism type_ of the simple components.
`SymbolicWedderburn.jl` makes no efforts to compute the types or defining fields,
it's primary goal is to compute symbolic/numerical Wedderburn-Artin isomorphism in a form usable for (polynomial) optimization. `Wedderga` also contains much more sophisticated methods for computing _a complete set of orthogonal primitive idempotents_ (i.e. a minimal projection system) through Shoda pairs.
In principle those idempotents could be computed using `Oscar.jl` and used in `SymbolicWedderburn.jl`.
