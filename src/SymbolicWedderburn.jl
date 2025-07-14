module SymbolicWedderburn

using LinearAlgebra
using SparseArrays
using Primes

using Cyclotomics
using GroupsCore
import AbstractPermutations as AP
import AbstractPermutations: degree
import PermutationGroups as PG
using StarAlgebras

### This is technically pirating from StarAlgebras;
# a stopgap solution until StarAlgebras is updated to v0.3 (?)
function Base.getindex(
    b::StarAlgebras.Basis{P},
    p::P) where P <: AP.AbstractPermutation
    return b.rbasis[p]
end

export symmetry_adapted_basis, WedderburnDecomposition
export basis,
    # degree, # too common name to export
    direct_summands,
    invariant_vectors,
    issimple,
    multiplicity

include("Characters/Characters.jl")
using .Characters
import .Characters: row_echelon_form!
import .Characters.FiniteFields

include("ext_homomorphisms.jl")
include("actions.jl")
include("group_action_error.jl")
include("action_characters.jl")
include("matrix_projections.jl")
include("image_basis.jl")
include("minimal_projections.jl")
include("direct_summands.jl")
include("sa_basis.jl")
include("wedderburn_decomposition.jl")

include("precompile.jl")

end # module
