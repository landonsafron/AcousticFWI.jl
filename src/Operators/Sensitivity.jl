
    ####### NEED TO TAKE THE REAL PART OF OUT???????????????


    #### HOW TO INPUT SPARSE MATRICES INTO A FUNCTION WITH PREDEFINED TYPES???????


function Sensitivity{T<:Complex}()

# This function is the forward and adjoint sensitivity/jacobian operator
# for multiple shots.
#
# INPUTS:     IN        - Data residuals for adjoint; Model perturbation for forward
#             adj       - true for adjoint; false for forward
#             R         - Restriction operator
#             H         - LU factorization of the Helmholtz operator
#             U         - Predicted forward-propagated wavefield for each shot (stored as a matrix with dimensions nz*nx,ns)
#             A         - Diagonal matrix containing the complex-valued attenuation coefficients
#
# OUTPUTS:    OUT       - Model perturbation for adjoint; Data residuals for forward

    ng = size(R,1)
    ns = size(U,2)
    n = size(A,2)

    if adj==false
        OUT = zeros(T,ng*ns)
        for ishot = 1:ns
            OUT[(ishot-1)*ng+1:ishot*ng] = Sensitivity(IN,adj,R=R,H=H,U=U[:,ishot],A=A)
        end
    end

    if adj==true
        OUT = zeros(T,n)
        for ishot = 1:ns
            OUT += Sensitivity(IN[(ishot-1)*ng+1:ishot*ng],adj,R=R,H=H,U=U[:,ishot],A=A)
        end
    end

    return OUT

end




function Sensitivity{T<:Complex}(IN::Array{T,1},adj::Bool;R::Array{Int,2},H,U::Array{T,1},A::Array{T,1})

# This function is the forward and adjoint sensitivity/jacobian operator
# for a single shot.
#
# INPUTS:     IN        - Data residual for adjoint; Model perturbation for forward
#             adj       - true for adjoint; false for forward
#             R         - Restriction operator
#             H         - LU factorization of the Helmholtz operator
#             U         - Predicted forward-propagated wavefield
#             A         - Diagonal matrix containing the complex-valued attenuation coefficients
#
# OUTPUTS:    OUT       - Model perturbation for adjoint; Data residual for forward

    if adj==false
        OUT = -w^2*R*H\(spdiagm(U)*A*IN)
    end

    if adj==true
        OUT = -w^2*A'*spdiagm(U)'*(H'\(R'*IN))
    end

    return OUT

end
