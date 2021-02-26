Base.promote_rule(
    ::Type{Character{T,Cl}},
    ::Type{VirtualCharacter{S,Cl}},
) where {T,S,Cl} = VirtualCharacter{promote_type(T, S),Cl}

for f in (:+, :-)
    @eval begin
        function Base.$f(χ::ClassFunction, ψ::ClassFunction)
            @assert χ.inv_of == ψ.inv_of
            @assert conjugacy_classes(χ) === conjugacy_classes(ψ)
            return VirtualCharacter(
                $f.(values(χ), values(ψ)),
                χ.inv_of,
                conjugacy_classes(χ),
            )
        end
    end
end

for C in (:Character, :VirtualCharacter)
    @eval begin
        Base.:*(χ::$C, c::Number) =
            VirtualCharacter(c .* values(χ), χ.inv_of, conjugacy_classes(χ))

        Base.:*(c::Number, χ::$C) =
            VirtualCharacter(c .* values(χ), χ.inv_of, conjugacy_classes(χ))

        Base.:/(χ::$C, c::Number) =
            VirtualCharacter(values(χ) ./ c, χ.inv_of, conjugacy_classes(χ))
    end
end
