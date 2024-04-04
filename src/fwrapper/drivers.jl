function get_atom_orb(Z::Int64)::Int64
    n = @ccall libDFTATOM.__drivers_MOD_get_atom_orb(Z::Ref{Cint})::Int64
    n
end

# wrap
#subroutine atom_lda(Z, r_min, r_max, a, N, n_orb, &
#    no, lo, fo, ks_energies, E_tot, R, Rp, V_tot, density, orbitals, &
#    reigen_eps, reigen_max_iter, mixing_eps, mixing_alpha, &
#    mixing_max_iter, perturb) bind(c)
function atom_lda(Z, mesh; reigen_atol=1e-6, mixing_atol=1e-6, mixing_alpha=0.5, mixing_max_iter=100, reigen_max_iter=50, perturb=true)
    r_min = mesh.r_min
    r_max = mesh.r_max
    r = mesh.r
    rp = mesh.rp
    a = mesh.a
    N = length(r)

    n_orbs = get_atom_orb(Z)
    no = zeros(Cint, n_orbs)
    lo = zeros(Cint, n_orbs)
    fo = zeros(Float64, n_orbs)
    ks_energies = zeros(Float64, n_orbs)
    E_tot = Ref{Float64}(0.0)
    r = mesh.r
    rp = mesh.rp
    V_tot = zeros(Float64, N)
    ρ = zeros(Float64, N)
    orbitals = Array{Float64}(undef, N, n_orbs)
    

    @ccall libDFTATOM.atom_lda(Z::Ref{Cint}, r_min::Ref{Float64}, r_max::Ref{Float64}, a::Ref{Float64}, N::Ref{Cint}, n_orbs::Ref{Cint}, no::Ptr{Cint}, lo::Ptr{Cint}, fo::Ptr{Float64}, ks_energies::Ptr{Float64}, E_tot::Ptr{Float64}, r::Ptr{Float64}, rp::Ptr{Float64}, V_tot::Ptr{Float64}, ρ::Ptr{Float64}, orbitals::Ptr{Float64}, reigen_atol::Ref{Float64}, reigen_max_iter::Ref{Cint}, mixing_atol::Ref{Float64}, mixing_alpha::Ref{Float64}, mixing_max_iter::Ref{Cint}, perturb::Ref{Cint})::Cvoid
    info = (; E_tot=E_tot[], energies=ks_energies, ρ=ρ, orbitals=orbitals, V_tot=V_tot)
    info
end

