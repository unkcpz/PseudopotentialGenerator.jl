module PspGen

include("potential.jl")
include("ae.jl")
include("pseudolize.jl")
include("hydrogen_like.jl")

# oncvpsp wrapper
module ONCVPSP
    export ONCVPSP

    include("oncvpsp_wrapper/oncvpsp.jl")
end

end
