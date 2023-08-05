module PspGen

include("potential.jl")
include("ae.jl")
include("pseudolize.jl")
include("hydrogen_like.jl")
include("radial_grids.jl")

# oncvpsp wrapper
module ONCVPSP
    export ONCVPSP

    include("oncvpsp_wrapper/oncvpsp.jl")
end

end
