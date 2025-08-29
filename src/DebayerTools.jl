module DebayerTools

import Base.maximum
using ColorTypes: colorview
using ColorVectorSpace
using Statistics: mean
using LinearAlgebra: norm
using ImageCore

include("debayer.jl")
include("utils.jl")

end
