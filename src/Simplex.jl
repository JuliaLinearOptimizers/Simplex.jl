module Simplex

using LinearAlgebra, SparseArrays
import SparseArrays.halfperm!, Base.permute!!
include("simplexauxiliar.jl")
include("simplexluup.jl")
include("simplexinv.jl")

end
