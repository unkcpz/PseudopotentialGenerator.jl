
function tf(Z::Int, r::Float64)
    b = (0.69395656/Z) ^ (1.0/3.0)
    x = r / b
    xs = √x

    t = Z / (1.0 + xs * (0.02747 - x*(0.1486 - 0.007298*x)) + x * (1.243 + x*(0.2302 + 0.006944*x)))

    if t < 1.0
        t = 1.0
    end

    -t / r
end

function tf(Z::Int, rgrid::Vector{Float64})
    pot = zeros(Float64, length(rgrid))
    for i in 1:length(rgrid)
        pot[i] = tf(Z, rgrid[i])
    end

    pot
end