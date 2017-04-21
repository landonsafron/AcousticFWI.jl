module AcousticFWI

# Julia modules
using Seismic

# Inversion
export FWI
include("Inversion/FWI.jl")

# Modeling
export Helmholtz
include("Modeling/Helmholtz.jl")

# Operators
export Attenuation,
Laplacian,
MassMatrix,
Restriction,
Sensitivity,
Source
include("Operators/Attenuation.jl")
include("Operators/Laplacian.jl")
include("Operators/MassMatrix.jl")
include("Operators/Restriction.jl")
include("Operators/Sensitivity.jl")
include("Operators/Source.jl")

end