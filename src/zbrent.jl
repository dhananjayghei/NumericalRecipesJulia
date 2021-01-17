#=
Using Brent’s method, find the root of a function `func` known to lie between `x1` and `x2`.
The root, returned as zbrent, will be refined until its accuracy is `tol`.
=#

function zbrent(x1::Float64, x2::Float64, FUNC::Function; itmax=200, tol=√eps(Float64))
    a, b, c, d, e, min1, min2 = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    a = x1
    b = x2
    fa = FUNC(a)
    fb = FUNC(b)
    # xmin = Float64[]
    if (fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)
        error("roots must be bracketed for zbrent")
    end
    c = b
    fc = fb
    for _ in 1:1:itmax
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
            c = a
            fc = fa
            d = b-a
            e = d
        end
        if (abs(fc) < abs(fb))
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end
        tol1 = 2.0*eps(Float64)*abs(b) + 0.5*tol      ## Convergence check
        xm = 0.5*(c-b)
        if (abs(xm) <= tol1) || (fb == 0.0)
            xmin = b
            return xmin
        end
        if (abs(e) >= tol1) && (abs(fa) > abs(fb))
            s = fb/fa
            if (a == c)
                p = 2.0*xm*s
                q = 1.0-s
            else
                q = fa/fc
                r = fb/fc
                p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
                q = (q-1.0)*(r-1.0)*(s-1.0)
            end
            if (p > 0.0)
                q = -q
            end
            p = abs(p)
            min1 = 3.0*xm*q-abs(tol1*q)
            min2 = abs(e*q)
            if (2.0*p < (min1 < min2 ? min1 : min2))
                e = d                               ## Accept interpolation.
                d = p/q
            else
                d = xm                              ## Interpolation failed, use bisection.
                e = d
            end
        else
            d = xm                                  ## Bounds decreasing too slowly, use bisection
            e = d
        end
        a = b
        fa = fb
        if (abs(d) > tol1)
            b = b+d
        else
            b = b+copysign(tol1, xm)
        end
        fb = FUNC(b)
        xmin = b
    end
    println("Maximum number of iterations reached.")
    return xmin
end


## Example and comparizon with the Brent Root Finding method of Roots.jl
# lxa, lxb = 1.e-8, 1.
# a, b, c, fc = mnbrak(lxa, lxb, x -> fngridtrans(x, true))
# @time br1 = zbrent(a, c, x -> fngridtrans(x, true))
# isapprox(fngridtrans(br1, true), 0.0, atol=√eps(Float64))
# using Roots
# @time br2 = find_zero(x -> fngridtrans(x, true), (a, c), Roots.Brent())
# isapprox(br1, br2)
