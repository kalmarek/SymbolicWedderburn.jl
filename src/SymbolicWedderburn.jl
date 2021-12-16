module SymbolicWedderburn

using LinearAlgebra
using SparseArrays
using Primes

using Cyclotomics
using GroupsCore
using PermutationGroups
using StarAlgebras

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

include("util.jl")
include("actions.jl")
include("action_characters.jl")
include("matrix_projections.jl")
include("minimal_projections.jl")
include("direct_summands.jl")
include("sa_basis.jl")
include("wedderburn_decomposition.jl")

end # module
