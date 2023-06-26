struct DirectSummand{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    basis::M
    multiplicity::Int
    degree::Int

    function DirectSummand(
        basis::M,
        multiplicity::Integer,
        degree::Integer,
    ) where {M<:AbstractMatrix}
        @assert size(basis, 1) < size(basis, 2)
        # @assert 1 ≤ degree ≤ size(basis, 1)
        let (pr_rank, r) = divrem(size(basis, 1), multiplicity)
            @assert 1 ≤ pr_rank ≤ degree
            @assert r == 0
        end
        return new{eltype(M),M}(basis, multiplicity, degree)
    end
end

image_basis(ds::DirectSummand) = ds.basis
PermutationGroups.degree(ds::DirectSummand) = ds.degree
multiplicity(ds::DirectSummand) = ds.multiplicity
pr_rank(ds::DirectSummand) = div(size(image_basis(ds), 1), multiplicity(ds))
issimple(ds::DirectSummand) = isone(pr_rank(ds))

Base.size(ds::DirectSummand) = size(ds.basis)
Base.@propagate_inbounds function Base.getindex(ds::DirectSummand, i...)
    return image_basis(ds)[i...]
end

function SparseArrays.sparse(ds::DirectSummand)
    sp = sparse(image_basis(ds))
    return DirectSummand(sp, multiplicity(ds), degree(ds), issimple(ds))
end

function SparseArrays.droptol!(
    ds::DirectSummand{T,<:AbstractSparseArray},
    tol,
) where {T}
    droptol!(image_basis(ds), tol)
    return ds
end

function simple_projection(
    ds::SymbolicWedderburn.DirectSummand,
    G::Group,
    hom::SymbolicWedderburn.InducedActionHomomorphism;
    atol=eps(eltype(ds))*max(size(ds)...)
)
    issimple(ds) && return ds
    @assert degree(ds) > 1
    @assert rem(size(ds, 1), degree(ds)) == 0
    @assert multiplicity(ds) == size(ds, 1)÷degree(ds)

    # ds is a projection onto subspace of dimension m·d
    # where m=multiplicity(ds), d = degreee(ds)
    # here we separate the image into m orthogonal (simple) subspaces
    # each of dimension d

    Uπ = SymbolicWedderburn.image_basis(ds)
    T = eltype(ds)
    @assert T == Float64
    subspaces = Vector{typeof(Uπ)}()
    qrs = LinearAlgebra.QRCompactWY{T, Matrix{T}}[]
    sizehint!(subspaces, multiplicity(ds))
    sizehint!(qrs, multiplicity(ds))

    @assert isapprox(Uπ*Uπ', I, atol=atol)

    for (i, b) in enumerate(eachrow(Uπ))
        b_orth = if isempty(subspaces)
            b
        else
            b_orth = b - sum(_orth_proj(V, b, QR) for (V, QR) in zip(subspaces, qrs))
        end
        norm(b_orth, 2) < atol && continue

        # we found a new subspace!
        # zero_small!(b_orth, atol=atol)
        b_orth ./= norm(b_orth, 2)
        # seach subspace is of dimension degree(ds)
        orbit_subspace = Matrix{T}(undef, length(b_orth), degree(ds))
        orbit_subspace[:, 1] .= b_orth

        dim = 1
        # cached qrfact to be reused when b_orth is already in the orbit_subspace
        qrfact = qr(@view orbit_subspace[:, 1:dim])
        for g in G
            isone(g) && continue
            k = action(hom, g, b_orth)
            @assert isapprox(norm(k, 2), 1.0, atol=atol) "Norm of translated vector should be ≈ 1.0, got $(norm(k, 2))"

            # check if we've seen this one
            for c in eachcol(orbit_subspace)
                isapprox(k, c, atol=atol) && continue
            end
            # or maybe it's already in the span of others?
            residual = k - _orth_proj(@view(orbit_subspace[:, 1:dim]), k, qrfact)
            (rnorm = norm(residual, 2)) < atol && continue
            if !isempty(subspaces)
                before_norm = norm(residual, 2)
                residual = residual - sum(_orth_proj(V, residual, QR) for (V, QR) in zip(subspaces, qrs))
                @info "norms" before_norm after_norm=norm(residual, 2)
                (rnorm = norm(residual, 2)) < atol && continue
            end

            # ok, it's not so we found a new vector in the subspace!
            dim += 1
            orbit_subspace[:, dim] .= residual ./ norm(residual, 2)
            qrfact = qr(@view orbit_subspace[:, 1:dim])

            if dim == degree(ds)
                @info orbit_subspace' * orbit_subspace
                @assert isapprox(orbit_subspace' * orbit_subspace, I, atol=atol)
                push!(subspaces, orbit_subspace)
                push!(qrs, qrfact)
                break
            end
        end
        # check that we exited via break, i.e. found all vectors in the subspace
        @assert dim == degree(ds)
    end
    # and that we found exactly the right number of subspaces
    @assert length(subspaces) == multiplicity(ds) "Expected $(multiplicity(ds)) subspaces, found $(length(subspaces))"

    # as the projection we take the average of all vectors in the subspace
    # new_Uπ = typeof(Uπ)(transpose(reduce(hcat, sum.(subspaces, dims=2))))
    # new_Uπ ./= degree(ds)
    # TODO: math says the the first one would be equally fine
    new_Uπ = typeof(Uπ)(transpose(reduce(
        hcat,
        (@view(s[:, 1]) for s in subspaces)
    )))

    return SymbolicWedderburn.DirectSummand(
        new_Uπ,
        degree(ds),
        multiplicity(ds),
        true
    )
end
