using SpecialPolynomials, Polynomials

export R_nl

laguerre_poly(n, α, x) = basis(Laguerre{α}, n)(x)

function R_nl(n, l, r; Z=1)
    nr = n - l - 1
    a = 1 / Z 
    r0 = 2 * r / (n * a)

    # normalization coeffecient
    C = sqrt((2 * Z / n)^3 * factorial(nr) / (2 * n * factorial(n + l)))

    return C * r0^l * exp(-r0 / 2) * laguerre_poly(nr, 2 * l + 1, r0)
end