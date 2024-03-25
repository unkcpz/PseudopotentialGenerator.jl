"""
    integrate(y, x; method::Symbol)

Integrate the function y(x) over x using the specified method (default is F_Trapezoidal7).
The method must be one of the following (f_ prefix means the implementation is in Fortran):
    * `trapz1`
    * `trapz3`
    * `trapz5`
    * `trapz7`
    * `simpson`
    * `adams`
"""
function integrate(y, rp; method::Symbol=:trapz7)::Float64
    # TODO: attach Fortran method here

    N = length(rp)
    if N < 8
        throw(ArgumentError("Length of integrated f must be at least 8"))
    end

    if method == :trapz1
        s = trapz1(y, rp)
    elseif method == :trapz3
        s = trapz3(y, rp)
    elseif method == :trapz5
        s = trapz5(y, rp)
    elseif method == :trapz7
        s = trapz7(y, rp)
    elseif method == :simpson
        s = simpson(y, rp)
    elseif method == :adams
        s = adams(y, rp)
    else
        throw(ArgumentError("Invalid method: $method"))
    end

    s
end

function trapz1(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = (g[1] + g[N]) / 2
    s += sum(g[2:N-1])
    s
end

function trapz3(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = (
        9 * (g[1] + g[N]) + 
        28 * (g[2] + g[N-1]) + 
        23 * (g[3] + g[N-2])
    ) / 24
    s += sum(g[4:N-3])
    s
end

function trapz5(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = (
        475 * (g[1] + g[N]) +
        1902 * (g[2] + g[N-1]) +
        1104 * (g[3] + g[N-2]) +
        1586 * (g[4] + g[N-3]) +
        1413 * (g[5] + g[N-4])
    ) / 1440
    s += sum(g[6:N-5])
    s
end

function trapz7(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = (
        36799 * (g[1] + g[N]) +
        176648 * (g[2] + g[N-1]) +
        54851 * (g[3] + g[N-2]) +
        177984 * (g[4] + g[N-3]) +
        89437 * (g[5] + g[N-4]) +
        130936 * (g[6] + g[N-5]) +
        119585 * (g[7] + g[N-6])
    ) / 120960
    s += sum(g[8:N-7])
    s
end

function simpson(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = 0
    for i in 2:2:N-1
        s += g[i-1] + 4 * g[i] + g[i+1]
    end
    s /= 3
    if N % 2 == 0
        # If N is even, add the last slice
        s += (5 * g[N] + 8 * g[N-1] - g[N-2]) / 12
    end
    s
end

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
    adams_interp_outward(y, i)

Interpolate the Adams method outward from the center of the interval.
"""
function adams_interp_outward(y::Vector{Float64}, i::Int)::Float64
    r = +(9 * y[i+1] + 19 * y[i] - 5 * y[i-1] + y[i-2]) / 24
    r
end

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
        #println("dftatom dydx: ", dydx)

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
