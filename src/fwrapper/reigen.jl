function solve_radial_eigenproblem(n, l, Z, V, mesh; tol=1e-9, max_iter=200, E_window=[-8000.0, 0], E_ini = -1000.0, rel = 0, perturb = false)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Bool}
    N = length(mesh.r)
    P = zeros(Float64, N)
    Q = zeros(Float64, N)
    Vs = V.(mesh.r)
    c = SPEED_OF_LIGHT
    Emin, Emax = E_window

    E = Ref{Float64}(0.0)
    converged = Ref{Int32}(0)

    @ccall libDFTATOM.solve_radial_eigenproblem(n::Ref{Int32}, l::Ref{Int32}, E_ini::Ref{Float64}, tol::Ref{Float64}, max_iter::Ref{Int32}, mesh.r::Ptr{Float64}, mesh.rp::Ptr{Float64}, Vs::Ptr{Float64}, Z::Ref{Int32}, c::Ref{Float64}, rel::Ref{Int32}, perturb::Ref{Cint}, Emin::Ref{Float64}, Emax::Ref{Float64}, converged::Ref{Int32}, E::Ref{Float64}, Q::Ptr{Float64}, P::Ptr{Float64}, N::Ref{Int32})::Cvoid

    if converged[] != 0
        is_converged = false
    else
        is_converged = true
    end

    E[], Q, P, is_converged
end