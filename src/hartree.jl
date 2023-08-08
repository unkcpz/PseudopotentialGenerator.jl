export hartree

# solve Possion equation and for hartree energy


hartree(R::RGrid, ρ::Vector{Float64}, q::Float64) = hartree_oncv(R::RGrid, ρ::Vector{Float64}, q::Float64)

function hartree_oncv(R::RGrid, ρ::Vector{Float64}, q::Float64)
    # interpolation inward integral
    aii(yy::Vector{Float64}, jj::Int64) = -4.166666666667e-2 * (9.0*yy[jj-1] + 19.0 * yy[jj] - 5.0 * yy[jj+1] + yy[jj+2])
    
    rvp = zeros(Float64, R.N)
    rv = zeros(Float64, R.N)

    al = 0.01 * log(R.r[101] / R.r[1])
    println(al)
    for i in 1:R.N
        rvp[i] = ρ[i] * al * R.r[i] ^ 3
    end

    rv[end-2:end] .= q
    for i in R.N-2:-1:2
        #rv[i-1] = rv[i] + integrate(R, rvp, R.r[1], R.r[i])
        rv[i-1] = rv[i] + aii(rvp, i)
    end

    for i in 1:R.N
        rvp[i] = ρ[i] * al * R.r[i] ^ 2
    end

    tv = 0.0
    for i in R.N-2:-1:2
        #tv += integrate(R, rvp, R.r[1], R.r[i])
        tv += aii(rvp, i)
        rv[i-1] = rv[i-1] - R.r[i-1] * tv
    end
    Y₀ = rv 

    # calculate the hartree potential
    ϕ = Y₀ ./ R.r

    ϕ
end

function hartree_pe(R::RGrid, ρ::Vector{Float64}, q::Float64)
    # calculate the hartree screening function Y₀
    Y₀ = zeros(Float64, R.N)
    f = ρ .* R.r .^ 2
    for i in 1:R.N
        r = R.r[i]
        rmin = R.r[1]
        rmax = R.r[end]
        Y₀[i] = integrate(R, f, rmin, r) + r * integrate(R, f ./ R.r, r, rmax)
    end

    # Force potential asymptotic to q
    if q != 0 && Y₀[end] !=0
        Y₀ = Y₀ / Y₀[end] * q
    end

    # calculate the hartree potential
    ϕ = Y₀ ./ R.r

    ϕ
end