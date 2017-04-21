function Helmholtz{T<:AbstractFloat}(isz::Array{Int,1},isx::Array{Int,1},ot::Array{T,1},vp::Array{T,2},wav::Array{T,1},fmin::T,fmax::T,nf::Int,nt::Int,dz::T,dx::T,dt::T,ext::Int=10,atten_max::T=1.,flag=2)

# This function solves the acoustic Helmholtz equation for a prescribed
# source wavelet and range of frequencies. A complex-numbered attenuating
# boundary is used to damp the boundary reflections. This code is intended
# for simultaneous sources.
#
# INPUTS:     isz       - z-index for each source location
#             isx       - x-index for each source location 
#             ot        - Activation time for each source
#             vp        - Acoustic velocity model
#             wav       - Source wavelet (time domain)
#             fmin      - Minimum frequency to model
#             fmax      - Maximum frequency to model
#             nf        - Number of frequency bins; 'nf' must be greater than or equal to 'length(wav)' and 'nt'
#             nt        - Number of time steps to perform modeling for
#             dz        - Grid spacing in z-direction
#             dx        - Grid spacing in x-direction
#             dt        - Sampling interval (time step)
#             ext       - Thickness (number of grid points) of absorbing boundary region
#             atten_max - Maximum complex amplitude in the attenuating boundary layer
#             flag      - Output frequency-domain data if flag==1, Output time-domain data if flag==2
#
# OUTPUTS:    u         - Solution to Helmhotz equation; stored as a 3D volume with dimensions (w,z,x) or (t,z,x)

    nz = size(vp,1) + 2*ext
    nx = size(vp,2) + 2*ext

    L = Laplacian(nz,nx,dz,dx)
    M = MassMatrix(vp,nz,nx,ext) .* Attenuation(nz,nx,ext,atten_max)

    WAV = fft([wav;zeros(nf-length(wav))])
    fs = 1/dt
    df = fs/nf
    faxis = fftshift(-fs/2:df:fs/2-df)
    waxis = 2*pi*faxis

    _,iwmin = findmin(abs(fmin-faxis))
    _,iwmax = findmin(abs(fmax-faxis))
    U = complex(zeros(nf,nz,nx))
    for iw = iwmin:iwmax
        w = waxis[iw]
        H = L + w^2*M
        s = Source(isz,isx,ot,WAV,waxis,w,nz,nx,ext)
        U_tmp = H\s
        U_tmp = reshape(U_tmp,nz,nx)
        U[iw,:,:] = U_tmp
        U[nf-iw+2,:,:] = conj(U_tmp)
    end

    if flag==1
        return U
    elseif flag==2 
        return real(ifft(U,1))[1:nt,ext+1:end-ext,ext+1:end-ext]
    end

end





function Helmholtz{T<:AbstractFloat}(isz::Int,isx::Int,ot::T,vp::Array{T,2},wav::Array{T,1},fmin::T,fmax::T,nf::Int,nt::Int,dz::T,dx::T,dt::T,ext::Int=10,atten_max::T=1.;flag=2)

# This function solves the acoustic Helmholtz equation for a prescribed
# source wavelet and range of frequencies. A complex-numbered attenuating
# boundary is used to damp the boundary reflections. This code is intended
# for a single source.
#
# INPUTS:     isz       - z-index for source location
#             isx       - x-index for source location 
#             ot        - Activation time for the source
#             vp        - Acoustic velocity model
#             wav       - Source wavelet (time domain)
#             fmin      - Minimum frequency to model
#             fmax      - Maximum frequency to model
#             nf        - Number of frequency bins; 'nf' must be greater than or equal to 'length(wav)' and 'nt'
#             nt        - Number of time steps to perform modeling for
#             dz        - Grid spacing in z-direction
#             dx        - Grid spacing in x-direction
#             dt        - Sampling interval (time step)
#             ext       - Thickness (number of grid points) of absorbing boundary region
#             atten_max - Maximum complex amplitude in the attenuating boundary layer
#             flag      - Output frequency-domain data if flag==1, Output time-domain data if flag==2
#
# OUTPUTS:    U         - Solution to Helmhotz equation; stored as a 3D volume with dimensions (w,z,x) or (t,z,x)
#
# Contributors: Landon Safron
# Institution: University of Alberta
# Date: March 1, 2017

    nz = size(vp,1) + 2*ext
    nx = size(vp,2) + 2*ext

    L = Laplacian(nz,nx,dz,dx)
    M = MassMatrix(vp,nz,nx,ext) .* Attenuation(nz,nx,ext,atten_max)

    WAV = fft([wav;zeros(nf-length(wav))])
    fs = 1/dt
    df = fs/nf
    faxis = fftshift(-fs/2:df:fs/2-df)
    waxis = 2*pi*faxis

    _,iwmin = findmin(abs(fmin-faxis))
    _,iwmax = findmin(abs(fmax-faxis))
    U = complex(zeros(nf,nz,nx))
    for iw = iwmin:iwmax
        w = waxis[iw]
        H = L + w^2*M
        s = Source(isz,isx,ot,WAV,waxis,w,nz,nx,ext)
        U_tmp = H\s
        U_tmp = reshape(U_tmp,nz,nx)
        U[iw,:,:] = U_tmp
        U[nf-iw+2,:,:] = conj(U_tmp)
    end

    if flag==1
        return U[:,ext+1:end-ext,ext+1:end-ext]
    elseif flag==2 
        return real(ifft(U,1))[1:nt,ext+1:end-ext,ext+1:end-ext]
    end

end



