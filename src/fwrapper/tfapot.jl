# Thomas Fermi potential implemented from ONCVPSP src/tfapot.f90
# ! generalized Thomas-Fermi atomic potential
# (Copy from the comment from ONCVPSP)
# !...to an article of N. H. March ( "The Thomas-Fermi Approximation in 
# ! Quantum Mechanics", Adv. Phys. 6, 1 (1957)). He has the formula,
# ! but it is not the result of his work. The original publication is: 
# !     R. Latter, Phys. Rev. 99, 510 (1955).
# ! He says that it''s an analytic fit to an improved calculation of the 
# ! potential distribution of a Thomas-Fermi atom without exchange first 
# ! performed by Miranda (C. Miranda, Mem. Acc. Italia 5, 285 (1934)).
# !                                 Alexander Seidl, TU Munich
export tf
function tf(Z::Int, r::Float64)
    zz = convert(Float64, Z)
    pot = ccall((:tfapot_, LIBONCV_PATH), Float64, (Ref{Float64}, Ref{Float64}), r, zz)
    pot
end