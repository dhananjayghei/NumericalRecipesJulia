#=
## --------------------------- Golden Section Search --------------------------- ##
Given a function FUNC, and given a bracketing triplet of abscissas ax, bx, cx
(such that bx is between ax and cx, and FUNC(bx) is less than both FUNC(ax) and FUNC(cx)),
this routine performs a golden section search for the minimum, isolating it to a
fractional precision of about tol. The abscissa of the minimum is returned as xmin,
and the minimum function value is returned as golden, the returned function value.
=#

function golden(ax::Float64, bx::Float64, cx::Float64, FUNC::Function; tol=√eps(Float64))
    R = .61803399 ## Golden Ratios
    C = 1. - R
    x0 = ax
    x3 = cx
    if abs(cx-bx) > abs(bx-ax)
        x1 = bx
        x2 = bx + C*(cx-bx)
    else
        x2 = bx
        x1 = bx - C*(bx-ax)
    end
    f1 = FUNC(x1)
    f2 = FUNC(x2)
    while (abs(x3-x0) > tol*(abs(x1)+abs(x2)))
        if (f2 < f1)
            x0 = x1
            x1 = x2
            x2 = R*x1 + C*x3
            f1 = f2
            f2 = FUNC(x2)
        else
            x3 = x2
            x2 = x1
            x1 = R*x1 + C*x0
            f2 = f1
            f1 = FUNC(x1)
        end
    end
    if f1 < f2
        golden = f1
        xmin = x1
    else
        golden = f2
        xmin = x2
    end
    return xmin, golden
end

## This is the golden section search algorithm using the brackets derived from the mnbrak function.
