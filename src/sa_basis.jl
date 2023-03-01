function affordable_real(
    irreducible_characters,
    multiplicities=fill(1, length(irreducible_characters)),
)
    irr_real = similar(irreducible_characters, 0)
    mls_real = similar(multiplicities, 0)
    for (i, χ) in pairs(irreducible_characters)
        ι = Characters.frobenius_schur(χ)
        if abs(ι) == 1 # real or quaternionic
            @debug "real/quaternionic:" χ
            push!(irr_real, χ)
            push!(mls_real, multiplicities[i])
        else # complex one...
            cχ = conj(χ)
            k = findfirst(==(cχ), irreducible_characters)
            @assert k !== nothing
            @debug "complex" χ conj(χ)=irreducible_characters[k]
            if k > i # ... we haven't already observed a conjugate of
                @assert multiplicities[i] == multiplicities[k]
                push!(irr_real, χ + cχ)
                push!(mls_real, multiplicities[i])
            end
        end
    end

    return irr_real, mls_real
end

function symmetry_adapted_basis(
    G::PermutationGroups.AbstractPermutationGroup,
    S::Type = Rational{Int};
    semisimple::Bool = false,
)
    tbl = CharacterTable(S, G)
    return symmetry_adapted_basis(eltype(tbl), tbl, semisimple = semisimple)
end

"""
    symmetry_adapted_basis([T::Type,] G::AbstractPermutationGroup[, S=Rational{Int};
        semisimple=false])
Compute a decomposition of `V=𝕂ⁿ` (with `n = degree(G)`) into
subspaces invariant under the natural permutation action of `G`.

Consider the decomposition of `V` into irreducible (simple) subspaces
    `V ≅ m₁V₁ ⊕ ⋯ ⊕ mᵣVᵣ`
and write `Wᵢ = mᵢVᵢ` (the `mᵢ`-fold direct sum of `Vᵢ`). The decomposition is
returned as a vector of `DirectSummand{T}`s (blocks) corresponding
to the distinct irreducible characters of `G` (i.e. types of action, here from `1`
to `r`). Each block contains a basis for a `G`-invariant subspace in `V` (`Vᵢ` or `Wᵢ`). The blocks are guaranteed to be orthogonal to **each other**, however vectors within a single block may *not* be orthogonal.

If `T<:LinearAlgebra.BlasFloat` BLAS routines will be used to orthogonalize
vectors within each block.

!!! note:
Each returned block is invariant (as a subspace) under the action of `G`, which
means that the action may still e.g. permute (row) vectors, but only *within* each block.

Arguments:
* `T` controls the type of coefficients of the returned basis. Unless specified,
the coefficients will be computed exactly in the field of cyclotomic numbers.
If you know that the group has rational characters only (which happens e.g. for
the full symmetric groups) You may specify `Rational{Int}` here. For a group
with complex characters specifying `T<:Real` will result in the computation of the realified basis.
* `S` controls the types of `Cyclotomic`s used in the computation of
character table. Exact type are preferred. For larger groups `G` (or if
overflow occurs during the computation of characters) specifying
`Rational{BigInt}` might be necessary.
* `semisimple::Bool` controls the nature of the the returned basis.
  - `semisimple=true`: the returned basis consists of orthogonal blocks
  (`DirectSummand`s) which define an isomorphism
  `V ≅ Wₖ ⊕ ⋯ ⊕ Wₖ`.
  associated to isotypical components `Wₖ`, which are (in general) semi-simple.
  I.e. each direct summand `ds` may further decopose into `G`-invariant
  subspaces `Wₖ ≅ mₖVₖ`, all of the same _type_. Multiplicity `mₖ` can be
  obtained by calling `multiplicity(@ref)`.
  - `semisimple=false`: (the default) In addition to finding blocks `Wₖ`, an
  effort to find _minimal projection system_ is made, i.e all, some (or none!)
  of the returned blocks corresponds to a **projection** `V → Wₖ ≅ mₖVₖ → Vₖ` for a single irreducible subspace `Vₖ`. This means that some blocks can not be further decomposed into nontrivial `G`-invariant subspaces.
"""
symmetry_adapted_basis(
    T::Type,
    G::PermutationGroups.AbstractPermutationGroup,
    S::Type = Rational{Int};
    semisimple::Bool = false,
) = symmetry_adapted_basis(T, CharacterTable(S, G), semisimple = semisimple)

function symmetry_adapted_basis(
    T::Type,
    tbl::CharacterTable;
    semisimple::Bool = false,
)
    irr, multips = _constituents_decomposition(
        action_character(T, conjugacy_classes(tbl), tbl),
        tbl,
    )

    if T <: Real
        irr, multips = affordable_real(irr, multips)
    end

    if semisimple || all(isone ∘ degree, irr)
        return _symmetry_adapted_basis(T, irr, multips)
    else
        RG = _group_algebra(parent(tbl))
        return _symmetry_adapted_basis(T, irr, multips, RG)
    end
end

"""
    symmetry_adapted_basis([T::Type,] G::Group, action, basis[, S=Rational{Int}];
        semisimple=false)
Compute a decomposition of `V=⟨basis⟩` into subspaces invariant under the
given action of `G`.

It is assumed that `G` acts on a subset of basis and the action needs to be
extended to the whole `basis`. If `G` is a permutation group already acting on
the whole `basis`, a call to `symmetry_adapted_basis(G)` is preferred.

* For inducing the action `basis` needs to be indexable and iterable
(e.g. in the form of an `AbstractVector`).

"""
function symmetry_adapted_basis(
    G::Group,
    action::Action,
    basis,
    S::Type = Rational{Int};
    semisimple=false,
)
    tbl = CharacterTable(S, G)
    ehom = CachedExtensionHomomorphism(parent(tbl), action, basis, precompute=true)
    return symmetry_adapted_basis(eltype(tbl), tbl, ehom, semisimple=semisimple)
end

function symmetry_adapted_basis(
    T::Type,
    G::Group,
    action::Action,
    basis,
    S::Type = Rational{Int};
    semisimple=false,
)
    tbl = CharacterTable(S, G)
    ehom = CachedExtensionHomomorphism(parent(tbl), action, basis, precompute=true)
    return symmetry_adapted_basis(T, tbl, ehom, semisimple=semisimple)
end

function symmetry_adapted_basis(
    T::Type,
    tbl::CharacterTable,
    ehom::InducedActionHomomorphism;
    semisimple=false,
)
    ψ = action_character(ehom, tbl)

    irr, multips = _constituents_decomposition(ψ, tbl)
    if T <: Real
        irr, multips = affordable_real(irr, multips)
        @debug "Decomposition into real character spaces:
        degrees:        $(join([lpad(d, 6) for d in degree.(irr)], ""))
        multiplicities: $(join([lpad(m, 6) for m in multips], ""))"

    end

    if semisimple || all(isone ∘ degree, irr)
        return _symmetry_adapted_basis(T, irr, multips, ehom)
    else
        RG = _group_algebra(parent(tbl))
        return _symmetry_adapted_basis(T, irr, multips, RG, ehom)
    end
end

function _constituents_decomposition(ψ::Character, tbl::CharacterTable)
    irr = irreducible_characters(tbl)
    degrees = degree.(irr)
    multiplicities = constituents(ψ)

    @debug "Decomposition into character spaces:
    degrees:        $(join([lpad(d, 6) for d in degrees], ""))
    multiplicities: $(join([lpad(m, 6) for m in multiplicities], ""))"

    @assert dot(multiplicities, degrees) == degree(ψ)
    "Something went wrong: characters do not constitute a complete basis for action:
    $(dot(multiplicities, degrees)) ≠ $(degree(ψ))"

    present_irreps = [i for (i, m) in enumerate(multiplicities) if m ≠ 0]
    return irr[present_irreps], multiplicities[present_irreps]
end

function _symmetry_adapted_basis(
    T::Type,
    irr::AbstractVector{<:Character},
    multiplicities::AbstractVector{<:Integer},
    hom=nothing
)
    res = map(zip(irr, multiplicities)) do (µ, m)
        Threads.@spawn begin
            µT = eltype(µ) == T ? µ : Character{T}(µ)
            image = hom === nothing ? image_basis(µT) : image_basis(hom, µT)
            simple = size(image, 1) == m
            deg = degree(µ)
            @assert size(image, 1) == (simple ? m : m*deg) "incompatible projection dimension: $(size(image, 1)) ≠ $(simple ? m : m*deg)"
            if deg == 1
                @assert simple "Central projection associated to character is not simple unless its degree == 1"
            end
            DirectSummand(image, m, deg, simple)
        end
    end
    return fetch.(res)
end

function _symmetry_adapted_basis(
    T::Type,
    irr::AbstractVector{<:Character},
    multiplicities::AbstractVector{<:Integer},
    RG::StarAlgebra{<:Group},
    hom=nothing,
)
    mps, simples = minimal_projection_system(irr, RG)
    degrees = degree.(irr)
    res = map(zip(mps, multiplicities, degrees, simples)) do (µ, m, deg, simple)
        Threads.@spawn begin
            µT = eltype(µ) == T ? µ : AlgebraElement{T}(µ)
            image = hom === nothing ? image_basis(µT) : image_basis(hom, µT)
            @assert size(image, 1) == (simple ? m : m*deg) "incompatible projection dimension: $(size(image, 1)) ≠ $(simple ? m : m*deg)"
            DirectSummand(image, m, deg, simple)
        end
    end
    direct_summands = fetch.(res)

    for (χ, ds) in zip(irr, direct_summands)
        if issimple(ds) && (d = size(ds, 1)) != (e = multiplicity(ds)*sum(constituents(χ).>0))
            throw("The dimension of the projection doesn't match with simple summand multiplicity: $d ≠ $e")
        end
    end

    return direct_summands
end
