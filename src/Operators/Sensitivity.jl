function Sensitivity()

# This function is the forward and adjoint sensitivity/jacobian operator
# for multiple shots.

    for ishot = 1:ns
        SOMETHING = Sensitivity()
    end

end




function Sensitivity()

# This function is the forward and adjoint sensitivity/jacobian operator
# for a single shot.

# S = real(R*dudm) = real( -w^2 * R * inv(H) * diag(U) * diag(A) )
# S' = real( -w^2 * diag(A)' * diag(U)' * inv(H)' * R' )

    # grad = -real(w^2*A'*spdiagm(U)'*(H'\(R'*r)))

end