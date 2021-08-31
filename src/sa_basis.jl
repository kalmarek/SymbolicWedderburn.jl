function image_basis(χ::Character)
    mpr = isirreducible(χ) ? matrix_projection_irr(χ) : matrix_projection(χ)
    image, pivots = image_basis!(mpr)
    return image[1:length(pivots), :]
end

function image_basis(hom::InducedActionHomomorphism, χ::Character)
    mpr = isirreducible(χ) ? matrix_projection_irr(hom, χ) : matrix_projection(hom, χ)
    image, pivots = image_basis!(mpr)
    return image[1:length(pivots), :]
end

function image_basis(α::AlgebraElement)
    image, pivots = image_basis!(matrix_projection(α))
    return image[1:length(pivots), :]
end

function image_basis(hom::InducedActionHomomorphism, α::AlgebraElement)
    image, pivots = image_basis!(matrix_projection(hom, α))
    return image[1:length(pivots), :]
end

struct DirectSummand{T, M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    basis::M
    multiplicity::Int
    degree::Int
    simple::Bool
end

StarAlgebras.basis(ds::DirectSummand) = ds.basis
issimple(ds::DirectSummand) = ds.simple
degree(ds::DirectSummand) = ds.degree
multiplicity(ds::DirectSummand) = ds.multiplicity

Base.size(ds::DirectSummand) = size(ds.basis)
Base.@propagate_inbounds Base.getindex(ds::DirectSummand, i...) =
    basis(ds)[i...]

function affordable_real(
    irreducible_characters,
    multiplicities=fill(1, length(irreducible_characters)),
)
    irr_real = similar(irreducible_characters, 0)
    mls_real = similar(multiplicities, 0)
    for (i, χ) in pairs(irreducible_characters)
        ι = frobenius_schur(χ)
        if abs(ι) == 1 # real or quaternionic
            @debug "real/quaternionic:" χ
            push!(irr_real, χ)
            push!(mls_real, multiplicities[i])
        else # complex one...
            cχ = conj(χ)
            k = findfirst(==(cχ), irreducible_characters)
            @assert !isnothing(k)
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

"""
    symmetry_adapted_basis([T::Type,] G::AbstractPermutationGroup[, S=Rational{Int};
        semisimple=false])
Compute a basis for the linear space `ℝⁿ` (where `n = degree(G)`) which is
invariant under the symmetry of `G`.

The coefficients of the invariant basis are returned in (orthogonal) blocks
corresponding to irreducible characters of `G`.

Arguments:
* `S` controls the types of `Cyclotomic`s used in the computation of
character table. Exact type are preferred. For larger groups `G` `Rational{BigInt}`
might be necessary.
* `T` controls the type of coefficients of the returned basis.
* `semisimple`: if set to `false` (the default) an effort to find minimal
projection system is made, so that the blocks give isomorphism to simple
summands. Otherwise semisimple decomposition is computed.

!!! Note:
Each returned block (a `DirectSummand`) is invariant under the action of `G`,
which means that the action may still e.g. permute (row) vectors , but only
*within* each block. The blocks are guaranteed to be orthogonal.
If `T<:LinearAlgebra.BlasFloat` BLAS routines will be used to orthogonalize
vectors within each block.
"""
function symmetry_adapted_basis(
    G::AbstractPermutationGroup,
    S::Type = Rational{Int};
    semisimple::Bool = false,
)
    tbl = CharacterTable(S, G)
    return symmetry_adapted_basis(eltype(tbl), tbl, semisimple = semisimple)
end

symmetry_adapted_basis(
    T::Type,
    G::AbstractPermutationGroup,
    S::Type = Rational{Int};
    semisimple::Bool = false,
) = symmetry_adapted_basis(T, CharacterTable(S, G), semisimple = semisimple)

function symmetry_adapted_basis(
    T::Type,
    tbl::CharacterTable;
    semisimple::Bool = false,
)
    irr, multips = _constituents_decomposition(
        action_character(conjugacy_classes(tbl), tbl),
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
Compute a decomposition of basis into (semi)simple subspaces which are invariant under
the action of `G`.

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
    return symmetry_adapted_basis(eltype(tbl), tbl, action, basis, semisimple=semisimple)
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
    return symmetry_adapted_basis(T, tbl, action, basis, semisimple=semisimple)
end

function symmetry_adapted_basis(
    T::Type,
    tbl::CharacterTable,
    action::Action,
    basis;
    semisimple=false,
)
    ehom = CachedExtensionHomomorphism(parent(tbl), action, basis, precompute=true)
    ψ = action_character(ehom, tbl)

    irr, multips = _constituents_decomposition(ψ, tbl)
    if T <: Real
        irr, multips = affordable_real(irr, multips)
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

function _group_algebra(G::Group)
    @assert isfinite(G)
    b = StarAlgebras.Basis{UInt16}(vec(collect(G)))
    RG = if order(Int, G) <= (typemax(UInt16)>>2)
        StarAlgebra(G, b, (length(b), length(b)), precompute=true)
        # cache is about ~ 1Gb
    else
        StarAlgebra(G, b)
    end
    return RG
end

function _symmetry_adapted_basis(
    T::Type,
    irr::AbstractVector{<:Character},
    multiplicities::AbstractVector{<:Integer},
    hom=nothing
)
    res = map(zip(irr, multiplicities)) do (µ, m)
        @spawn_compat begin
            µT = eltype(µ) == T ? µ : Character{T}(µ)
            image = isnothing(hom) ? image_basis(µT) : image_basis(hom, µT)
            simple = size(image, 1) == m
            deg = degree(µ)
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
        @spawn_compat begin
            µT = eltype(µ) == T ? µ : AlgebraElement{T}(µ)
            image = isnothing(hom) ? image_basis(µT) : image_basis(hom, µT)
            if simple
                @assert size(image, 1) == m "The dimension of the projection doesn't match with simple summand multiplicity"
            end
            DirectSummand(image, m, deg, simple)
        end
    end
    return fetch.(res)
end
