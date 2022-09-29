"""
Function format from paper: 
Berghe, G. Vanden, V. Fack, and H. E. De Meyer. Journal of computational and applied mathematics 28 (1989): 391-401.
"""

"""
V = u0 / (1+t) + u1 * t / (1+t)^2, 
where,
t = exp((r - r0) / a)
"""
function woods_saxon(rgrid::Vector{Float64}, u0::Float64, u1::Float64, r0::Float64, a::Float64)
    t = @. exp((rgrid - r0) / a)
    v = @. u0 / (1+t) + u1 * t / (1+t)^2 

    v
end

function hulthen(rgrid::Vector{Float64}, Z::Int, α::Float64)
    v = @. (-2 * Z / α) * (exp(-rgrid) / (1 - exp(-rgrid)))

    v
end