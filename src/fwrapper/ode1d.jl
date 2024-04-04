"""
    integrate(y, x; method::Symbol)

Integrate the function y(x) over x using the specified method (default is F_Trapezoidal7).
The method must be one of the following (f_ prefix means the implementation is in Fortran):
    * `f_trapz1`
    * `f_trapz3`
    * `f_trapz5`
    * `f_trapz7`
    * `f_simpson`
    * `f_adams`
"""
function integrate(y, rp; method::Symbol=:f_trapz7)::Float64
    N = length(rp)
    if N < 8
        throw(ArgumentError("Length of integrated f must be at least 8"))
    end

    if method == :f_trapz1
        s = @ccall libDFTATOM.integrate_trapz_1(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Cint})::Float64
    elseif method == :f_trapz3
        s = @ccall libDFTATOM.integrate_trapz_3(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Cint})::Float64
    elseif method == :f_trapz5
        s = @ccall libDFTATOM.integrate_trapz_5(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Cint})::Float64
    elseif method == :f_trapz7
        s = @ccall libDFTATOM.integrate_trapz_7(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Cint})::Float64
    elseif method == :f_simpson
        s = @ccall libDFTATOM.integrate_simpson(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Cint})::Float64
    elseif method == :f_adams
        s = @ccall libDFTATOM.integrate_adams(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Cint})::Float64
    else
        throw(ArgumentError("Invalid method: $method"))
    end

    s
end

"""
    rk4_integrate(r, y0, C1, C2, C1mid, C2mid, max_val)

Integrate the 2nd-order ode using the Runge-Kutta 4th order method (convert to 2 1st-order odes).
"""
function rk4_integrate(r::Vector{Float64}, y0::Vector{Float64}, C1::Vector{Float64}, C2::Vector{Float64}, C1mid::Vector{Float64}, C2mid::Vector{Float64}, max_val::Float64)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    N = length(r)
    y1 = zeros(Float64, N)
    y2 = zeros(Float64, N)
    imax = @ccall libDFTATOM.rk4_integrate(r::Ptr{Float64}, y0::Ptr{Float64}, C1::Ptr{Float64}, C2::Ptr{Float64}, C1mid::Ptr{Float64}, C2mid::Ptr{Float64}, max_val::Ref{Float64}, y1::Ptr{Float64}, y2::Ptr{Float64}, N::Ref{Cint})::Int64

    y1, y2, imax
end

function midpoints(x::Vector{Float64}, r::Vector{Float64})::Vector{Float64}
    N = length(x)
    x_mid = zeros(Float64, N-1)
    @ccall libDFTATOM.get_midpoints(r::Ref{Float64}, x::Ref{Float64}, N::Ref{Cint}, x_mid::Ptr{Float64})::Cvoid
    x_mid
end