function Laplacian{T<:AbstractFloat}(nz::Int,nx::Int,dz::T,dx::T)

# This function computes the 2D Laplace operator that acts on a vectorized
# version of a 2D field. The vectorization is assumed to be done using column
# major order.
#
# INPUTS:     nz        - Number of grid points in the z-direction
#             nx        - Number of grid points in the x-direction
#             dz        - Grid spacing in z-direction
#             dx        - Grid spacing in x-direction
#
# OUTPUTS:    L         - 2D Laplace operator

    Dzz = (1/dz^2)*spdiagm((ones(nz-1),-2*ones(nz),ones(nz-1)),(-1,0,1))
    Dxx = (1/dx^2)*spdiagm((ones(nx-1),-2*ones(nx),ones(nx-1)),(-1,0,1))

    Mzz = kron(speye(nx),Dzz)
    Mxx = kron(Dxx,speye(nz))

    L = Mzz + Mxx

    return L

end