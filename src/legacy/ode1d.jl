function adams(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = trapz1(rp[1:4], y[1:4])
    for i in 4:N-1
        s += adams_interp_outward(g, i)
    end
    s
end

"""
    adams_extrapolation_outward(y, i)

Extrapolate the Adams method outward from the center of the interval.
"""
function adams_extrapolation_outward(y::Vector{Float64}, i::Int)::Float64
    r = +(55 * y[i] - 59 * y[i-1] + 37 * y[i-2] - 9 * y[i-3]) / 24
    r
end

"""
    adams_interp_outward(y, i)

Interpolate the Adams method outward from the center of the interval.
"""
function adams_interp_outward(y::Vector{Float64}, i::Int)::Float64
    r = +(9 * y[i+1] + 19 * y[i] - 5 * y[i-1] + y[i-2]) / 24
    r
end

"""
    adams_interp_outward_implicit(y, i)

Interpolate the Adams method outward from the center of the interval (implicit version).
"""
function adams_interp_outward_implicit(y::Vector{Float64}, i::Int)::Float64
    r = +(19 * y[i] - 5 * y[i-1] + y[i-2]) / 24
    r
end

"""
    adams_interp_inward_implicit(y, i)

Interpolate the Adams method inward from the center of the interval (implicit version).
"""
function adams_interp_inward_implicit(y::Vector{Float64}, i::Int)::Float64
    r = -(19 * y[i] - 5 * y[i+1] + y[i+2]) / 24
    r
end

"""
    rk4_integrate(f1, f2, r, y0, max_val)

Integrate the system of ODEs using the Runge-Kutta 4th order method.
Integrates the following set of ODEs:

    y1' = f1(y1, y2, r)
    y2' = f2(y1, y2, r)
"""
function rk4_integrate(f1::Function, f2::Function, r::Vector{Float64}, y0::Vector{Float64}, max_val::Float64)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    N = length(r)
    y1 = zeros(N)
    y2 = zeros(N)
    y = y0
    y1[1] = y0[1]
    y2[1] = y0[2]
    imax = N

    # Initialize the derivatives
    dym = zeros(2)
    dyt = zeros(2)
    yt = zeros(2)
    dydx = zeros(2)

    for i in 2:N
        # Step size
        h = r[i] - r[i-1]

        dydx[1] = f1(y[1], y[2], r[i-1])
        dydx[2] = f2(y[1], y[2], r[i-1])
        @. yt = y + h / 2 * dydx

        dyt[1] = f1(yt[1], yt[2], r[i-1] + h / 2)
        dyt[2] = f2(yt[1], yt[2], r[i-1] + h / 2)
        @. yt = y + h / 2 * dyt

        dym[1] = f1(yt[1], yt[2], r[i-1] + h / 2)
        dym[2] = f2(yt[1], yt[2], r[i-1] + h / 2)
        @. yt = y + h * dym
        @. dym += dyt

        dyt[1] = f1(yt[1], yt[2], r[i])
        dyt[2] = f2(yt[1], yt[2], r[i])
        @. y += h / 6 * (dydx + dyt + 2 * dym)

        y1[i] = y[1]
        y2[i] = y[2]

        if y1[i] > max_val
            imax = i
            return y1, y2, imax
        end
    end

    y1, y2, imax
end

"""
    rk4_integrate(r, y0, C1, C2, C1mid, C2mid, max_val)

Integrate the system of ODEs using the Runge-Kutta 4th order method.
Integrates the following set of ODEs:

    y1' = y2
    y2' = C1 * y1 + C2 * y2
"""
function rk4_integrate(r::Vector{Float64}, y0::Vector{Float64}, C1::Vector{Float64}, C2::Vector{Float64}, C1mid::Vector{Float64}, C2mid::Vector{Float64}, max_val::Float64)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    N = length(r)
    y1 = zeros(N)
    y2 = zeros(N)
    y = y0
    y1[1] = y0[1]
    y2[1] = y0[2]
    imax = N

    # Initialize the derivatives
    dym = zeros(2)
    dyt = zeros(2)
    yt = zeros(2)
    dydx = zeros(2)

    for i in 2:N
        # Step size
        h = r[i] - r[i-1]

        dydx[1] =                            y[2]
        dydx[2] = C1[i-1] * y[1] + C2[i-1] * y[2]
        @. yt = y + h / 2 * dydx

        dyt[1] =                                   yt[2]
        dyt[2] = C1mid[i-1] * yt[1] + C2mid[i-1] * yt[2]
        @. yt = y + h / 2 * dyt

        dym[1] =                                   yt[2]
        dym[2] = C1mid[i-1] * yt[1] + C2mid[i-1] * yt[2]
        @. yt = y + h * dym
        @. dym += dyt

        dyt[1] =                         yt[2]
        dyt[2] = C1[i] * yt[1] + C2[i] * yt[2]
        @. y += h / 6 * (dydx + dyt + 2 * dym)

        y1[i] = y[1]
        y2[i] = y[2]

        if y1[i] > max_val
            imax = i
            return y1, y2, imax
        end
    end

    y1, y2, imax
end

"""
    rk4_integrate_poisson(r, y0, C1, C2, C1mid, C2mid, max_val)

Integrate the system of ODEs using the Runge-Kutta 4th order method.
Integrates the following set of ODEs:

    y1' = y2
    y2' = C1 + C2 * y2
"""
function rk4_integrate_poisson(r::Vector{Float64}, y0::Vector{Float64}, C1::Vector{Float64}, C2::Vector{Float64}, C1mid::Vector{Float64}, C2mid::Vector{Float64}, max_val::Float64)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    N = length(r)
    y1 = zeros(N)
    y2 = zeros(N)
    y = y0
    y1[1] = y0[1]
    y2[1] = y0[2]
    imax = N

    # Initialize the derivatives
    dym = zeros(2)
    dyt = zeros(2)
    yt = zeros(2)
    dydx = zeros(2)

    for i in 2:N
        # Step size
        h = r[i] - r[i-1]

        dydx[1] =                            y[2]
        dydx[2] = C1[i-1] + C2[i-1] * y[2]
        @. yt = y + h / 2 * dydx

        dyt[1] =                                   yt[2]
        dyt[2] = C1mid[i-1] + C2mid[i-1] * yt[2]
        @. yt = y + h / 2 * dyt

        dym[1] =                                   yt[2]
        dym[2] = C1mid[i-1] + C2mid[i-1] * yt[2]
        @. yt = y + h * dym
        @. dym += dyt

        dyt[1] =                         yt[2]
        dyt[2] = C1[i] + C2[i] * yt[2]
        @. y += h / 6 * (dydx + dyt + 2 * dym)

        y1[i] = y[1]
        y2[i] = y[2]

        if y1[i] > max_val
            imax = i
            return y1, y2, imax
        end
    end

    y1, y2, imax
end
