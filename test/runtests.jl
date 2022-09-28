using Test
using PspGen: woods_saxon, compute_rgrid, sol_orb

# lschfb test with woods_saxon 
u0 = -50.0
a = 0.6
r0 = 7.0
u1 = - u0 / a
Z = 1
rgrid, dr = compute_rgrid(Z)
v_ws = woods_saxon(rgrid, u0, u1, r0, a)

n = 12
l = 1
sol_orb(n, l, rgrid, v_ws, Z, true)



