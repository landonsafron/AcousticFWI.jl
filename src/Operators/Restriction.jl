function Restriction(nz::Int,nx::Int,ext::Int,igz::Array{Int,1},igx::Array{Int,1})

# This function creates the reciver restriction/sampling operator/matrix. 
#
# INPUTS:     nz        - Number of grid points in the z-direction
#             nx        - Number of grid points in the x-direction
#             igz       - z-indices for receiver locations
#             igx       - x-indices for receiver locations
#
# OUTPUTS:    R         - Restriction operator

    if length(igx)!=length(igz)
        error("The lengths of 'igx' and 'igz' must be the same.")
    end

    I = zeros(Int,length(igz))
    J = zeros(Int,length(igz))
    V = zeros(Int,length(igz))
    for i = 1:length(igz)
        I[i] = i
        J[i] = (ext+igx[i]-1)*nz + ext + igz[i]
        V[i] = 1
    end

    R = sparse(I,J,V,length(igz),nz*nx)

    return R

end