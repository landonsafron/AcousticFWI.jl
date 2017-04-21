module AcousticFWI

# Julia modules
using Seismic

# Inversion
export FWI
include("Inversion/FWI.jl")

# Modeling
export Helmholtz,
Source
include("Modeling/Helmholtz.jl")
include("Modeling/Source.jl")

# Operators
export Attenuation,
Laplacian,
MassMatrix,
Restriction,
SensitivityMultiShot
include("Operators/Attenuation.jl")
include("Operators/Laplacian.jl")
include("Operators/MassMatrix.jl")
include("Operators/Restriction.jl")
include("Operators/SensitivityMultiShot.jl")

end