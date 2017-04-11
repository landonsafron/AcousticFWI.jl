function FWI{T<:AbstractFloat}(vp0::Array{T,2},d::Array{T,2},wav::Array{T,1},isz::Array{Int,1},isx::Array{Int,1},igz::Array{Int,1},igx::Array{Int,1},ot::Array{T,1},fmin::T,fmax::T,nf::Int,dz::T,dx::T,dt::T,ext::Int=50,atten_max::T=2.,alpha=1,maxiter=5)

# This function performs frequency-domain acoustic full waveform inversion
# using the Gauss-Newton method. A muti-scale approach is adopted by
# treating each frequency as a new iteration of Gauss-Newton. The data
# is assumed to be sorted into a single super-shot. The receivers are
# assumed to be located in the same position for each shot within the
# super-shot.
#
# INPUTS:     vp0       - Initial guess for the velocity model
#             d         - Data (time domain) organized into a single super-shot
#             wav       - Source wavelet (time domain)
#             isz       - z-indices for source locations
#             isx       - x-indices for source locations
#             igz       - z-indices for receiver locations
#             igx       - x-indices for receiver locations
#             ot        - Activation time for each of the sources
#             fmin      - Minimum frequency to model
#             fmax      - Maximum frequency to model
#             nf        - Number of frequency bins; 'nf' must be greater than or equal to 'length(wav)'
#             dz        - Grid spacing in z-direction
#             dx        - Grid spacing in x-direction
#             dt        - Sampling interval (time step)
#             ext       - Thickness (number of grid points) of absorbing boundary region
#             atten_max - Maximum complex amplitude in the attenuating boundary layer
#             alpha     - Step size in gradient descent method; how far to step along the gradient
#             maxiter   - Maximum number of iterations of Gauss-Newton
#
# OUTPUTS:    vp        - Estimated velocity model

    nz = size(vp0,1) + 2*ext
    nx = size(vp0,2) + 2*ext

    L = Laplacian(nz,nx,dz,dx)
    A = Attenuation(nz,nx,ext,atten_max)
    M = SlownessSquared(vp0,nz,nx,ext)
    R = Restriction(nz,nx,ext,igz,igx)

    D = fft([d ; zeros(nf-size(d,1),size(d,2))])
    WAV = fft([wav;zeros(nf-length(wav))])
    fs = 1/dt
    df = fs/nf
    faxis = fftshift(-fs/2:df:fs/2-df)
    waxis = 2*pi*faxis

    _,iwmin = findmin(abs(fmin-faxis))
    _,iwmax = findmin(abs(fmax-faxis))

    iter = 1
    while iter<=maxiter #&& (NOT CONVERGED)
        for iw = iwmin:iwmax
            w = waxis[iw]
            s = Source(isz,isx,ot,WAV,waxis,w,nz,nx,ext)
            MassMatrix = w^2*M*A
            H = L + MassMatrix
            U = H\s
            r = R*U - D[iw,:]
            grad = -real(MassMatrix'*(conj(H)\(R'*r)))
            M = M + alpha*spdiagm(grad)
        end
        iter += 1
    end

    vp = sqrt(1./M)

    return vp

end
