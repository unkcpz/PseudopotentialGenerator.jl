struct Mesh
    r_min::Float64
    r_max::Float64
    a::Float64
    N::Int64
    r::Vector{Float64}
    rp::Vector{Float64}
    function Mesh(r_min::Float64, r_max::Float64, a::Float64, N::Int64)
        r = zeros(Float64, N+1)
        rp = zeros(Float64, N+1)

        if a < 0
            throw(ArgumentError("require a > 0"))
        elseif (abs(a - 1) < eps(typeof(a)))
            # Uniform grid if a = 1
            r = range(r_min, r_max, length=N+1)
            # rp = 0 for all i
        else
            # Exponential grid
            if (N > 1)
                beta = log(a) / (N - 1)
                alpha = (r_max - r_min) / (exp(beta * N) - 1)
                for i = 1:N+1
                    r[i] = r_min + alpha * (exp(beta * (i - 1)) - 1)
                    rp[i] = alpha * beta * exp(beta * (i - 1))
                end
            elseif N == 1
                r = [r_min, r_max]
            else
                throw(ArgumentError("require N > 0"))
            end
        end

        new(r_min, r_max, a, N, r, rp)
    end
end


function mesh_exp_deriv2(mesh::Mesh)::Vector{Float64}
    rpp = mesh_exp_deriv2(mesh.r_min, mesh.r_max, mesh.a, mesh.N)
    rpp
end


function mesh_exp_deriv2(r_min::Float64, r_max::Float64, a::Float64, N::Int64)::Vector{Float64}
    rpp = zeros(Float64, N+1)
    if a < 0
        throw(ArgumentError("require a > 0"))
    elseif (abs(a - 1) < eps(typeof(a)))
        # Uniform grid if a = 1, rpp = 0 for all i
        nothing
    else
        # Exponential grid
        if (N > 1)
            beta = log(a) / (N - 1)
            alpha = (r_max - r_min) / (exp(beta * N) - 1)
            for i = 1:N+1
                rpp[i] = alpha * beta^2 * exp(beta * (i - 1))
            end
        elseif N == 1
            rpp = [0, 0]
        else
            throw(ArgumentError("require N > 0"))
        end
    end

    rpp
end