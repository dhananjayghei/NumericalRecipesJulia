#=
## ------------------- Bracket the minimum of a univariate function ------------------- ##
Given a function FUNC(X), and given distinct initial points AX and BX, this routine
searches in the downhill direction (defined by the function as evaluated at the
initial points) and returns new points AX, BX, CX which bracket a minimum of the function.
Also returned are the function values at the three points, FA, FB and FC.
=#

function mnbrak(ax::Float64, bx::Float64, FUNC::Function)
    ## Parameters
    gold = 1.618034
    glimit = 100
    tiny = 1.e-20
    # The first parameter is the default ratio by which successive intervals are magnified;
    # the second is the maximum magnification allowed for a parabolic-fit step.
    fa = FUNC(ax)
    fb = FUNC(bx)
    if (fb > fa)
        dum = ax
        ax = bx
        bx = dum
        dum = fb
        fb = fa
        fa = dum
    end
    cx = bx + gold*(bx-ax)
    fc = FUNC(cx)
##    iter = 1
    while (fb >= fc) ## & iter <= 25
        r = (bx-ax)*(fb-fc)                    ## Compute u by parabolic  extrapolation from a b c.
        q = (bx-cx)*(fb-fa)                    ## TINY is used to prevent any possible division by zero.
        u = bx-((bx-cx)*q-(bx-ax)*r)/(2. *copysign(maximum([abs(q-r),tiny]),q-r))
        ulim = bx+glimit*(cx-bx)
        if (bx-u)*(u-cx) > 0                   ## Parabolic u is between b and c :  try  it
            fu = FUNC(u)
            if fu < fc                         ## Got a minimum between b and c.
                ax = bx
                fa = fb
                bx = u
                fb = fu
            elseif fu > fb                     ## Got a minimum between a and u.
                cx = u
                fc = fu
            end
            u = cx + gold*(cx-bx)              ## Parabolic fit was no use. Use default magnification.
            fu = FUNC(u)
        elseif (cx-u)*(u-ulim) > 0             ## Parabolic fit is between c and its allowed limit.
            fu = FUNC(u)
            if fu < fc
                bx = cx
                cx = u
                u = cx + gold*(cx-bx)
                fb = fc
                fc = fu
                fu = FUNC(u)
            end
        elseif (u-ulim)*(ulim-cx) >= 0        ## Limit  parabolic u to maximum allowed value.
            u = ulim
            fu = FUNC(u)
        else                                  ## Reject parabolic u, use default magnification.
            u = cx + gold*(cx-bx)
            fu = FUNC(u)
        end
        ax=bx                                 ## Eliminate  oldest  point  and  continue.
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
##        iter = iter + 1
    end
    return ax, bx, cx, fc
end
