struct DirectSummand{T, M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    basis::M
    multiplicity::Int
    degree::Int
    simple::Bool
end

image_basis(ds::DirectSummand) = ds.basis
issimple(ds::DirectSummand) = ds.simple
PermutationGroups.degree(ds::DirectSummand) = ds.degree
multiplicity(ds::DirectSummand) = ds.multiplicity

Base.size(ds::DirectSummand) = size(ds.basis)
Base.@propagate_inbounds Base.getindex(ds::DirectSummand, i...) =
    image_basis(ds)[i...]

function SparseArrays.sparse(ds::DirectSummand)
    sp = sparse(image_basis(ds))
    return DirectSummand(sp, multiplicity(ds), degree(ds), issimple(ds))
end

function SparseArrays.droptol!(
    ds::DirectSummand{T, <:AbstractSparseArray},
    tol
) where T
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
