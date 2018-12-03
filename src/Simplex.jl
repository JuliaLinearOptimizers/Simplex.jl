module Simplex
using Base.SparseArrays.halfperm!, Base.SparseArrays.increment!, Base.permute!!, Base.ipermute!!
include("simplexauxiliar.jl")
include("simplexluup.jl")
include("simplexinv.jl")

end
