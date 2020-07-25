AbstractAlgebra.exponent(G::AbstractAlgebra.Group) = AbstractAlgebra.exponent(AbstractAlgebra.conjugacy_classes(G))
AbstractAlgebra.exponent(cclasses::AbstractVector) = lcm(order.(first.(cclasses)))
dixon_prime(G::AbstractAlgebra.Group) = dixon_prime(order(G), AbstractAlgebra.exponent(G))

function dixon_prime(cclasses::AbstractVector)
    ordG = sum(length, cclasses)
    m = AbstractAlgebra.exponent(cclasses)
    return dixon_prime(ordG, m)
end

function dixon_prime(ordG::Integer, exponent::Integer)
    p = 2floor(Int, sqrt(ordG))
    while true
        p = nextprime(p+1)
        isone(p % exponent) && break # we need -1 to be in the field
    end
    return p
end

