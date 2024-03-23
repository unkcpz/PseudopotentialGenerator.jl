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
        s = @ccall libDFTATOM.integrate_trapz_1(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Int32})::Float64
    elseif method == :f_trapz3
        s = @ccall libDFTATOM.integrate_trapz_3(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Int32})::Float64
    elseif method == :f_trapz5
        s = @ccall libDFTATOM.integrate_trapz_5(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Int32})::Float64
    elseif method == :f_trapz7
        s = @ccall libDFTATOM.integrate_trapz_7(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Int32})::Float64
    elseif method == :f_simpson
        s = @ccall libDFTATOM.integrate_simpson(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Int32})::Float64
    elseif method == :f_adams
        s = @ccall libDFTATOM.integrate_adams(rp::Ptr{Float64}, y::Ptr{Float64}, N::Ref{Int32})::Float64
    else
        throw(ArgumentError("Invalid method: $method"))
    end

    s
end
