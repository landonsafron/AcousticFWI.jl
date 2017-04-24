function SensitivityMultiShot{T1<:Complex,T2<:Int,T3<:AbstractFloat}(IN::Array{T1,1},adj::Bool;R::SparseMatrixCSC{T2,T2}=sparse(ones(T2,1,1)),H::Base.SparseArrays.UMFPACK.UmfpackLU{T1,T2}=lufact(sparse(ones(T1,1,1))),U::Array{T1,2}=im*ones(T1,1,1),A::SparseMatrixCSC{T1,T2}=sparse(ones(T1,1,1)),w::T3=1.)

# This function is the forward and adjoint sensitivity/jacobian operator
# for multiple shots.
#
# INPUTS:     IN        - Data residuals for adjoint; Model perturbation for forward
#             adj       - true for adjoint; false for forward
#             R         - Restriction operator
#             H         - LU factorization of the Helmholtz operator
#             U         - Predicted forward-propagated wavefield for each shot (stored as a matrix with dimensions nz*nx,ns)
#             A         - Diagonal matrix containing the complex-valued attenuation coefficients
#             w         - Angular frequency
#
# OUTPUTS:    OUT       - Model perturbation for adjoint; Data residuals for forward

    ng = size(R,1)
    ns = size(U,2)
    n = size(A,2)

    if adj==false
        OUT = zeros(T1,ng*ns)
        for ishot = 1:ns
            OUT[(ishot-1)*ng+1:ishot*ng] = Sensitivity1Shot(IN,adj,R=R,H=H,U=U[:,ishot],A=A,w=w)
        end
    end

    if adj==true
        OUT = zeros(T1,n)
        for ishot = 1:ns
            OUT += Sensitivity1Shot(IN[(ishot-1)*ng+1:ishot*ng],adj,R=R,H=H,U=U[:,ishot],A=A,w=w)
        end
    end

    return OUT

end





function Sensitivity1Shot{T1<:Complex,T2<:Int,T3<:AbstractFloat}(IN::Array{T1,1},adj::Bool;R::SparseMatrixCSC{T2,T2}=sparse(ones(T2,1,1)),H::Base.SparseArrays.UMFPACK.UmfpackLU{T1,T2}=lufact(sparse(ones(T1,1,1))),U::Array{T1,1}=im*ones(T1,1),A::SparseMatrixCSC{T1,T2}=sparse(ones(T1,1,1)),w::T3=1.)

# This function is the forward and adjoint sensitivity/jacobian operator
# for a single shot.
#
# INPUTS:     IN        - Data residual for adjoint; Model perturbation for forward
#             adj       - true for adjoint; false for forward
#             R         - Restriction operator
#             H         - LU factorization of the Helmholtz operator
#             U         - Predicted forward-propagated wavefield
#             A         - Diagonal matrix containing the complex-valued attenuation coefficients
#             w         - Angular frequency
#
# OUTPUTS:    OUT       - Model perturbation for adjoint; Data residual for forward

    if adj==false
        OUT = -w^2*R*(H\(spdiagm(U)*A*IN))
    end

    if adj==true
        OUT = -w^2*A'*spdiagm(U)'*(H'\(R'*IN))
    end

    return OUT

end
