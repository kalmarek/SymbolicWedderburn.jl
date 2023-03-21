# A template for new action definitions to use with SymbolicWedderburn
# Intended to be expanded by the user

# Step 1: define your group, or use existing group
# must be conformant with GroupsCore

using GroupsCore

struct NewGroup <: GroupsCore.Groups end
struct NewGroupElement <: GroupsCore.GroupsElement end

# see examples:
# dihedral.jl
# cyclic.jl
# PermutationGroups

# Step 2: define your object to be acted on, or use existing object

# for objects in abstract algebra, use the parent-element design pattern
struct NewObject end
struct NewObjectElement end

# for other objects, just define the object type
struct NewObject end

# see examples: 
# (i) abstract algebra objects: Monoids, StarAlgebras
# (ii) other objects: AbstractVector, or DynamicPolynomials

# Step 3: define action
# Notes:
# - a group can act on an object in many ways, i.e. defining different actions
# - a group can act on different objects, both mathematically and programmatically

import SymbolicWedderburn: action

struct NewAction <: Action end

function coeff_type(::NewAction) end
function action(a::NewAction, g::NewGroupElement, ω::NewObjectElement)
    return ω^g::NewObjectElement
end

# if encounter error messages, may need to define
function induce(a::NewAction, hom::ExtensionHomomorphism, g::NewGroupElement) end
function decompose(x, hom::InducedActionHomomorphism) end

# see examples: 
# actions.jl -> implements By(Signed)Permutations, ByLinearTransformations
# action_polynomials.jl -> implements VariablePermutation
# ex_C2_linear.jl
# ex_robinson_form.jl


# Additional guidelines: please consult the built-in actions defined
# in SymbolicWedderburn and subtype properly to avoid re-implementing
# existing methods (suboptimally)

using SymbolicWedderburn, TypeTree
tt(SymbolicWedderburn.Action)
"""
4-element Vector{String}:
 "SymbolicWedderburn.Action\n"
 " ├─ SymbolicWedderburn.ByLinearTransformation\n"
 " │   └─ SymbolicWedderburn.BySignedPermutations\n"
 " └─ SymbolicWedderburn.ByPermutations\n"
"""