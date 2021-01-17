#=
Using the secant method, find the root of a function `FUNC` thought to lie between x1 and x2 .
The root, returned as rtsec, is refined until its accuracy is ±xacc .
=#
function rtsec(FUNC::Function, x1::Float64, x2::Float64, xacc::Float64)
    maxit = 30
    fl = FUNC(x1)
    f = FUNC(x2)
    if (abs(fl) < abs(f))                          ## Pick the bound with the smaller function value as the most recent guess.
        rtsec = x1
        xl = x2
        swap = fl
        fl = f
        f = swap
    else
        xl = x1
        rtsec = x2
    end
    for j in 1:maxit                               ## Secant loop
		dx=(xl-rtsec)*f/(f-fl)                     ## Increment with respect to latest value
		xl=rtsec
		fl=f
		rtsec=rtsec+dx
		f=FUNC(rtsec)
		if ((abs(dx) < xacc) || f == 0.0)          ## Convergence
			return rtsec
		end
	end
	println("rtsec: exceeded maximum iterations.\n")
end
