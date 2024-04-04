function solve_radial_eigenproblem(n, l, Z, V, mesh; tol=1e-9, max_iter=200, E_window=[-8000.0, 0], E_ini = -1000.0, rel = 0, perturb = false)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Bool}
    N = length(mesh.r)
    P = zeros(Float64, N)
    Q = zeros(Float64, N)
    Vs = V.(mesh.r)
    c = SPEED_OF_LIGHT
    Emin, Emax = E_window

    E = Ref{Float64}(0.0)
    converged = Ref{Cint}(0)

    @ccall libDFTATOM.solve_radial_eigenproblem(n::Ref{Cint}, l::Ref{Cint}, E_ini::Ref{Float64}, tol::Ref{Float64}, max_iter::Ref{Cint}, mesh.r::Ptr{Float64}, mesh.rp::Ptr{Float64}, Vs::Ptr{Float64}, Z::Ref{Cint}, c::Ref{Float64}, rel::Ref{Cint}, perturb::Ref{Cint}, Emin::Ref{Float64}, Emax::Ref{Float64}, converged::Ref{Cint}, E::Ref{Float64}, Q::Ptr{Float64}, P::Ptr{Float64}, N::Ref{Cint})::Cvoid

    if converged[] != 0
        is_converged = false
    else
        is_converged = true
    end

    E[], Q, P, is_converged
end