#=
Using the Newton-Raphson method, find the root of a function known to lie in the interval
[x1, x2]. The root `rtnewt` will be refined until its accuracy is known within ±xacc.
`funcd` is a user-supplied routine that returns both the function value and the first derivative of the
function at the point x.

Reference: Page 365, Chatper 9, Press et al (1992) Numerical Recipes in C
=#

function rtnewt(funcd::Function, x1::Float64, x2::Float64; xacc=√eps(Float64), jmax=20)
    df, dx, f, rtnewt = Float64[], Float64[], Float64[], Float64[]
    rtnewt = 0.5*(x1+x2)                                       ## Initial guess
    for j in 1:jmax
        f, df= funcd(rtnewt)
        dx = f/df
        rtnewt = rtnewt - dx
        if ((x1-rtnewt)*(rtnewt-x2) < 0.0)
            error("Jumped out of brackets in rtnewt.")
        end
        if (abs(dx) < xacc)                                    ## Convergence
            return rtnewt
        end
    end
    println("Maximum number of iterations reached.")
end
