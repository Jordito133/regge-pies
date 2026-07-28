# defects.jl
#
# 

# sign time- or space-like depending on Q
function hinge_signature_from_Q(
    q::Real;
    null_tol::Real=DEFAULT_NULL_TOL,
)
    # Q > 0 is spacelike, Q < 0 is timelike.
    if q > null_tol
        return :spacelike
    end
    if q < -null_tol
        return :timelike
    end
    # otherwise its zero
    error("Null hinge")
end

