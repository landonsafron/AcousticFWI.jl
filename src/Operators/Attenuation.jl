function Attenuation{T<:AbstractFloat}(nz::Int,nx::Int,ext::Int,atten_max::T)

# This function builds a complex-valued diagonal matrix, where the elements
# of the main diagonal contain the attenuation factor for each point in the
# in the modeling region AND the surrounding boundary region. It is assumed 
# that no attenuation occurs in the modeling regionl. It is also assumed 
# that all 4 sides of the modeling region are surrounded by an attenuating
# boundary material to prevent boundary reflections. The imaginary component
# in the attenuation factor is zero in the modeling region and increases
# parabolically outward in the boundary region.
#
# INPUTS:     nz        - Number of grid points in z-direction INCLUDING the absorbing boundary region
#             nx        - Number of grid points in x-direction INCLUDING the absorbing boundary region
#             ext       - Thickness (number of grid points) of attenuating boundary region
#             atten_max - Maximum complex amplitude in the attenuating boundary layer
#
# OUTPUTS:    A         - Diagonal matrix containing the attenuation factors

    atten = (ext:-1:1).^2/(ext.^2)*atten_max
    a_boundary1 = kron(ones(nz,1),atten')
    a_boundary2 = kron(atten,ones(1,nx))
    a = zeros(T,nz,nx)
    a[:,1:ext] += a_boundary1
    a[:,nx-ext+1:nx] += flipdim(a_boundary1,2)
    a[1:ext,:] += a_boundary2
    a[nz-ext+1:nz,:] += flipdim(a_boundary2,1)

    A = 1 - im*a

    return spdiagm((A[:]),(0))

end