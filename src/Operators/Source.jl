function Source{T1<:AbstractFloat,T2<:Complex}(isz::Array{Int,1},isx::Array{Int,1},ot::Array{T1,1},WAV::Array{T2,1},waxis::Array{T1,1},w::T1,nz::Int,nx::Int,ext::Int)

# This function computes the source field at the given frequency. The source
# field is vectorized using column major order. This code is intended for
# simultaneous sources.
# 
# INPUTS:     isz       - z-index for the source locations
#             isx       - x-index for the source locations
#             ot        - Activation time for the sources
#             WAV       - Source wavelet (frequency domain)
#             waxis     - Angular frequency axis; contains the angular frequency for each entry in 'WAV'
#             w         - Angular frequency to perform modeling for
#             nz        - Number of grid points in the z-direction
#             nx        - Number of grid points in the x-direction
#             ext       - Thickness (number of grid points) of absorbing boundary region
#
# OUTPUTS:    src       - Source field (vector)

    _,iw = findmin(abs(waxis-w))
    src = complex(zeros(nz*nx))
    for isrc = 1:length(isz)
        k = (isx[isrc]+ext-1)*nz + isz[isrc] + ext
        coeff = -WAV[iw]*exp(-im*w*ot[isrc])
        src[k] = coeff
    end

    return src

end





function Source{T1<:AbstractFloat,T2<:Complex}(isz::Int,isx::Int,ot::T1,WAV::Array{T2,1},waxis::Array{T1,1},w::T1,nz::Int,nx::Int,ext::Int)

# This function computes the source field at the given frequency. The source
# field is vectorized using column major order. This code is intended for only
# a single source.
# 
# INPUTS:     isz       - z-index for the source location
#             isx       - x-index for the source location
#             ot        - Activation time for the source
#             WAV       - Source wavelet (frequency domain)
#             waxis     - Angular frequency axis; contains the angular frequency for each entry in 'WAV'
#             w         - Angular frequency to perform modeling for
#             nz        - Number of grid points in the z-direction
#             nx        - Number of grid points in the x-direction
#             ext       - Thickness (number of grid points) of absorbing boundary region
#
# OUTPUTS:    src       - Source field (vector)

    _,iw = findmin(abs(waxis-w))
    k = (isx+ext-1)*nz + isz + ext
    coeff = -WAV[iw]*exp(-im*w*ot)

    src = complex(zeros(nz*nx))
    src[k] = coeff

    return src

end
