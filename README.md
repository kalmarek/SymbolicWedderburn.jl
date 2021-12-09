# SymbolicWedderburn.jl
[![CI](https://github.com/kalmarek/SymbolicWedderburn.jl/workflows/CI/badge.svg?branch=master)](https://github.com/kalmarek/SymbolicWedderburn.jl/actions)
[![codecov](https://codecov.io/gh/kalmarek/SymbolicWedderburn.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kalmarek/SymbolicWedderburn.jl)

Amazing package providing symbolic Wedderburn decompositions for group *-algebras.

Group *-algebras are closely related to
* Representations ρ : G → GL(V) of a group G on a k-vector space V, or
* k(G)-modules structure on V,
which naturally arise from (non)commutative polynomial optimization with symmetry.

By Masche's theorem, `V ≅ m₁V₁ ⊕ ⋯ ⊕ mᵣVᵣ ≅ W₁ ⊕ ⋯ ⊕ Wᵣ` is decomposed uniquely into isotypic components `Wᵢ` and nonuniquely into irreducible/simple components `Vᵢ` associated with a group G. By (symbolic) computation in k(G), SymbolicWedderburn is capable of producing exact isomorphism `V ≅ W₁ ⊕ ⋯ ⊕ Wᵣ`, and (in general) numerical isomorphism `V ≅ m₁V₁ ⊕ ⋯ ⊕ mᵣVᵣ`. These isomorphisms give a decomposition of `End(V)` in the sense of Wedderburn-Artin theorem.

Important! The GAP package Wedderga does decomposition of the group ring k(G)
(for finite group G, and k is a finite field of characteristic coprime to |G| or
an abelian number field). SymbolicWedderburn extends these functionalities and give decomposition even for `V`.
