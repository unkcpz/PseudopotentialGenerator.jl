using Roots
using SpecialFunctions
using Polynomials: fit, coeffs, Polynomial, derivative
using ForwardDiff

# nbes for the number of bessel functions used 
function rrkj(aewfc::Vector{Float64}, l::Int, rc::Float64, rgrid::Vector{Float64}, dr::Vector{Float64}, nbessel::Int=3)
    # find the rc in rgrid
    rc, ic = find_rc_ic(rgrid, rc)
    ae_norm = compute_ae_norm(aewfc, dr, ic)
    ae_deriv = compute_ae_deriv(aewfc, rgrid, ic)   # calculate to 3 order derivative

    # using bisect to find q (s) that [rc*jl(qi*rc)]'/(rc*jl(qi*rc)) = phi'(rc)/phi(rc)
    println(ae_deriv)
    ld_rbessel(q) = logder_rbessel(l, q, rc) - ae_deriv[2] / ae_deriv[1]
    qs = find_qs(ld_rbessel, nbessel, 0.0, 20.0)
    println("qs=", qs)
    println("Estimate cutoff= ", 0.5 * qs[end]^2)
    
    # find rrkj coefficients
    x0 = 0.0
    residual(x) = compute_rrkj_residual(x, qs, ae_norm, ae_deriv, rc, rgrid, dr, l)
    D(f) = x -> ForwardDiff.derivative(f, float(x))
    c2 = find_zero((residual, D(residual)), x0, Roots.Newton())

    # coefficient found
    c = solve_rrkj_linear(c2, qs, ae_deriv, rc, l)
    println("rrkj coefficient: ", c)

    pswfc .= aewfc # > rc the wfc is the same
    for i in 1:ic
        pswfc[i] = rrkj_func(rgrid[i], l, c, qs)
    end

    d2pswfc .= pswfc    # TODO 2nd derivative of pswfc

    pswfc, d2pswfc
end

qbess(l, q, r) = sphericalbesselj_fast(l, q * r)
qbessp(l, q, r) = sphericalbesseljp_fast(l, q * r)
function qbesspp(l, q, r) 
    x = q * r
    res = (l * (l+1) - x * x) * sphericalbesselj_fast(l,x) - 2 * x * sphericalbesseljp_fast(l,x) / r^2

    res
end
rqbess(l, q, r) = r * qbess(l, q, r)
rqbessp(l, q, r) = qbess(l, q, r) + r * qbessp(l, q, r)
rqbesspp(l, q, r) = 2 * qbessp(l, q, r) + r * qbesspp(l, q, r)
logder_rbessel(l, q, r) = rqbessp(l, q, r) ./ rqbess(l, q, r)

function find_rc_ic(rgrid, rc)
    ic = 1
    for i in eachindex(rgrid)
        if rgrid[i] <= rc < rgrid[i+1]
            rc = rgrid[i]
            ic = i
        end
    end 

    rc, ic
end

function compute_rrkj_residual(c2, qs, ae_norm, ae_deriv, rc, rgrid, dr, l)
    c = solve_rrkj_linear(c2, qs, ae_deriv, rc, l)

    _, ic = find_rc_ic(rgrid, rc)
    r = rgrid[1:ic]
    ps_norm = sum(rrkj_func(r, l, c, qs).^2 .* dr[1:ic])

    res = ps_norm - ae_norm
    println("residual", res)
    res
end

function solve_rrkj_linear(c2, qs, ae_deriv, rc, l)
    # rhs side
    b = zeros(Float64, 2)
    b[1] = ae_deriv[1]
    b[2] = ae_deriv[3]

    # lhs side
    A = zeros(Float64, 2, 2)
    for i in 1:3
        rq = rqbess(l, qs[i], rc)
        rqp = rqbesspp(l, qs[i], rc)
        if i < 3
            A[1, i] = rq
            A[2, i] = rqp
        else
            b[1] -= c2 * rq
            b[2] -= c2 * rq
        end
    end

    x = A \ b

    append!(x, c2)
end

function rrkj_func(r::Float64, l, c, qs)
    res = 0.0
    # println("c is: ", c)
    # println("qs is: ", qs)
    for i in eachindex(qs)
        res += c[i] * rqbesspp(l, qs[i], r)
    end

    res
end

function rrkj_func(rs::Vector{Float64}, l, c, qs) 
    res = zeros(Float64, length(rs))
    for (i, r) in enumerate(rs)
        res[i] = rrkj_func(r, l, c, qs)
    end

    res
end

function compute_ae_norm(aewfc, dr, ic)
    # density sum up
    ds = 0.0
    # println(aewfc)
    for i in 1:ic
        ds += aewfc[i] * aewfc[i] * dr[i]
    end

    ds
end

# fitting a polynamial and compute derivative
function compute_ae_deriv(aewfc, rgrid, ic)
    # if ic < 11 throw 
    rs = rgrid[ic-10:ic+10]
    fs = aewfc[ic-10:ic+10]
    poly = fit(rs, fs, 6) |> p -> round.(coeffs(p), digits=5) |> Polynomial
    polyp = derivative(poly)
    polypp = derivative(polyp)

    c = aewfc[ic]
    deriv = [c, polyp(c), polypp(c)]

    deriv
end

function find_qs(f, nbessel, q_min, q_max)
    qs = zeros(Float64, nbessel)
    q = q_min
    intv = 0.05 # the step to find root
    for i in 1:nbessel
        while q < q_max
            try
                qs[i] = find_zero(f, (q, q+intv), Bisection())
            catch exc
                if isa(exc, ArgumentError)
                    q += intv
                    continue
                end
            end

            q += intv
            
            break
        end
    end
    
    qs
    # qs = find_zeros(f, q_min, q_max)

    # qs[1:nbessel]
end

"""
    Copy from DFTK spherical_bessels.jl
    sphericalbesselj_fast(l::Integer, x::Number)
Returns the spherical Bessel function of the first kind j_l(x). Consistent with
https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions and with
`SpecialFunctions.sphericalbesselj`. Specialized for integer `0 <= l <= 5`.
"""
@fastmath function sphericalbesselj_fast(l::Integer, x::T)::T where {T}
    if l == 0
        iszero(x) && return one(T)
        return sin(x) / x
    end

    iszero(x) && return zero(T)

    l == 1 && return (sin(x) - cos(x) * x) / x^2
    l == 2 && return (sin(x) * (3 - x^2) + cos(x) * (-3x)) / x^3
    l == 3 && return (sin(x) * (15 - 6x^2) + cos(x) * (x^3 - 15x)) / x^4
    l == 4 && return (sin(x) * (105 - 45x^2 + x^4) + cos(x) * (10x^3 - 105x)) / x^5
    l == 5 && return (sin(x) * (945 - 420x^2 + 15x^4) + cos(x) * (-945x + 105x^3 - x^5)) / x^6
    error("The case l = $l is not implemented")
end

# d/dx(sin(x)/x) = (x cos(x) - sin(x))/x^2
# d/dx((sin(x) - cos(x) x)/x^2) = ((x^2 - 2) sin(x) + 2 x cos(x))/x^3
# d/dx((sin(x) (3 - x^2) + cos(x) (-3 x))/x^3) = ((4 x^2 - 9) sin(x) - x (x^2 - 9) cos(x))/x^4
# d/dx((sin(x) (15 - 6 x^2) + cos(x) (x^3 - 15 x))/x^4) = ((60 x - 7 x^3) cos(x) - (x^4 - 27 x^2 + 60) sin(x))/x^5
# d/dx((sin(x) (105 - 45 x^2 + x^4) + cos(x) (10 x^3 - 105 x))/x^5) = ((-11 x^4 + 240 x^2 - 525) sin(x) + x (x^4 - 65 x^2 + 525) cos(x))/x^6
# d/dx((sin(x) (945 - 420 x^2 + 15 x^4) + cos(x) (-945 x + 105 x^3 - x^5))/x^6) = (x (16 x^4 - 735 x^2 + 5670) cos(x) + (x^6 - 135 x^4 + 2625 x^2 - 5670) sin(x))/x^7

@fastmath function sphericalbesseljp_fast(l::Integer, x::T)::T where {T}
    if l == 0
        iszero(x) && return zero(T)
        return (x * cos(x) - sin(x)) / x^2
    end

    iszero(x) && return zero(T)

    l == 1 && return ((x^2 - 2) * sin(x) + 2x * cos(x)) / x^3
    l == 2 && return ((4x^2 - 9) * sin(x) - x * (x^2 - 9) * cos(x)) / x^4
    l == 3 && return ((60x - 7x^3) * cos(x) - (x^4 - 27x^2 + 60) * sin(x)) / x^5
    l == 4 && return ((-11x^4 + 240x^2 - 525) * sin(x) + x * (x^4 - 65x^2 + 525) * cos(x)) / x^6
    l == 5 && return (x * (16x^4 - 735x^2 + 5670) * cos(x) + (x^6 - 135x^4 + 2625x^2 - 5670) * sin(x)) / x^7
    error("The case l = $l is not implemented")
end