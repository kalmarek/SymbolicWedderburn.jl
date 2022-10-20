# A template for new action definitions to use with SymbolicWedderburn
# may be converted to a documentation

# Q: action convention left/right in SymbolicWedderburn? does this matter?
# the syntax action(type, groupelement, object) suggest a LEFT action i.e.
# groupelement⋅object, but the docstring of ByPermutations indicates RIGHT action

# We follow the design pattern ... in abstract algebra ...
# Explain parent object, element...
# Q: can this be found somewhere, to be precise, and also copypaste here?

# Step 1: define your group, or use existing group
# must be conformant with GroupsCore

using GroupsCore

struct NewGroup <: GroupsCore.Groups end
struct NewGroupElement <: GroupsCore.GroupsElement end

# see examples: dihedral.jl, cyclic.jl, PermutationGroups

# Step 2: define your object to be acted on, or use existing object

struct NewObject end
struct NewObjectElement end

# see examples: 
# Julia built-in AbstractVector
# DynamicPolynomials
# StarAlgebras

# Step 3: define action
# Notes:
# - a group can act on an object in many ways, i.e. defining different actions
# - a group can act on different objects, both mathematically and programmatically

import SymbolicWedderburn: action # use this to extend method

struct NewAction <: Action end # vs abstract type NewAction <: Action end ?
# Q: which is better for performance, dispatch? generally better design style?

function action(a::NewAction, g::NewGroupElement, o::NewObjectElement) end
# Q: would you like the user to extend SymbolicWedderburn.action?

function coeff_type(::NewAction) end
function induce(a::NewAction, hom::ExtensionHomomorphism, g::NewGroupElement) end
# Q: other required methods? decompose?


# see examples: 
# actions.jl -> implements ByPermutations, BySignedPermutations, ByLinearTransformations
# action_polynomials.jl
# https://github.com/kalmarek/SymbolicWedderburn.jl/blob/master/examples/ex_C2_linear.jl
# https://github.com/kalmarek/SymbolicWedderburn.jl/blob/5d88715b73ea1b8133e7580d398077c99c0c7d8b/examples/ex_robinson_form.jl#L39
# 
# Additional guidelines:

1. Action type hierarchy:

using SymbolicWedderburn, TypeTree
tt(SymbolicWedderburn.Action)

"""
4-element Vector{String}:
 "SymbolicWedderburn.Action\n"
 " ├─ SymbolicWedderburn.ByLinearTransformation\n"
 " │   └─ SymbolicWedderburn.BySignedPermutations\n"
 " └─ SymbolicWedderburn.ByPermutations\n"
"""

2. Cached actions:

3. Disambiguation?