# File contains implementation of the NN based solver.
using BSON: @load

include("initcond.jl")
include("fluxes.jl")
include("hotimefluxes.jl")
include("param.jl")
include("enosolver.jl")

using Plots
using Flux
using Zygote
using DifferentialEquations
pyplot()

include("nets.jl")
@load wdatafn weights
Flux.loadparams!(model, weights)

function getalpha!(alp, u)
	N, K = size(u)
	for i=snapw+1:N-snapw
		inp = vec(u[i-snapw:i+snapw-1, :])
		press = zeros(2*snapw)
		for l = 1:2*snapw
			press[l] = p(u[i-snapw+l, :])
		end
		inp = vcat(inp, press)
		alp[i, :] = max.(0.0, min.(1.0, model(inp)))
		alp[i, :] .= maximum(alp[i, :])
	end
	alp[1:snapw, :] .= 1.0
	alp[N-snapw:end, :] .= 1.0
end

function flux(alp, u, k, lambda)
	N, K = size(u)
	width=2
	if k < width || k > N-width
		return hEulerLLF(u[k, :], u[k+1, :])
	else
		return eulerHOconsflux(u[k-width+1:k+width, :], width, alp[k])
	end
end

function EfluxDD(alp, u, k, lambda)
        width=order+1
        width=6
        N, K = size(u)
        if k < width || k > N-width
                return HEulerLLF(u[k, :], u[k+1, :])
        else
                return alp[k, 1]*HEulerLLF2(u[k-order+1:k+order, :], lambda) + (1.0 - alp[k, 1])*FEulerST43(u[k-width+1:k+width, :], lambda)
        end
end

function FuncTRI(ul, um, ur, lambda, Func)
        unext = um + lambda*(hEulerCons(ul, um) - hEulerCons(um, ur))
        return Func(unext)
end

function FuncLLF(ul, um, ur, lambda, Func)
        unext = um + lambda*(hEulerLLF(ul, um) - hEulerLLF(um, ur))
        return Func(unext)
end

function FuncGT(ul, um, ur, lambda, alpl, alpr, Func)
        flux(alp, a, b) = alp.*hEulerLLF(a, b) + (1.0 .-alp).*hEulerCons(a, b)
        unext = um + lambda*(flux(alpl, ul, um) - flux(alpr, um, ur))
        return Func(unext)
end


function FuncLMRL(u, lambda, Func, alp)
        p = 2
        unext = u[p+1, :] + lambda*(eulerHOconsflux(u, p, alp) - fEuler(u[p+1, :]))
        return Func(unext)
end
function FuncLMRR(u, lambda, Func, alp)
        p = 2
        unext = u[p, :] + lambda*(fEuler(u[p, :]) - eulerHOconsflux(u, p, alp))
        return Func(unext)
end


rhofunc(u) = -u[1]
pfunc(u) = -p(u)

function LowAlpBoundpp(uar, lambda, func)
        width = 2
        N, K = size(uar)
        alps = zeros(N)
        pl, pr, sl, sr = zeros(N), zeros(N), zeros(N), zeros(N)
        for i=1:N-1
                if i < width || i > N-width
                # Low order
                        pr[i] = FuncTRI(uar[i, :], uar[i, :], uar[i+1, :], 2*lambda, func)
                        pl[i] = FuncTRI(uar[i, :], uar[i+1, :], uar[i+1, :], 2*lambda, func)
                        sr[i] = FuncLLF(uar[i,:], uar[i, :], uar[i+1, :], 2*lambda, func)
                        sl[i] = FuncLLF(uar[i, :], uar[i+1, :], uar[i+1, :], 2*lambda, func)

#               High Order
                else
                        pr[i] = FuncLMRR(uar[i-1:i+2, :], 2*lambda, func, 0.0)
                        pl[i] = FuncLMRL(uar[i-1:i+2, :], 2*lambda, func, 0.0)
                        sr[i] = FuncLMRR(uar[i-1:i+2, :], 2*lambda, func, 1.0)
                        sl[i] = FuncLMRL(uar[i-1:i+2, :], 2*lambda, func, 1.0)

                end
                sr[i] = FuncLLF(uar[i,:], uar[i, :], uar[i+1, :], 2*lambda, func)
                sl[i] = FuncLLF(uar[i, :], uar[i+1, :], uar[i+1, :], 2*lambda, func)
        end
        for i=1:N-1
                alpl = eps^2
                alpr = eps^2
                if pl[i]  > -eps^2
                        alpl = abs(pl[i])/(abs(pl[i])-sl[i])
                end
                if pr[i]  > -eps^2
                        alpr = abs(pr[i])/(abs(pr[i])-sr[i])
                end
                alps[i] = min(1.0, max(alpl, alpr))

        end
        return alps
end


function CheckProdDD(uar, lambda)
        N, K = size(uar)
        pc =  zeros(N)
	alp = zeros(N, 1)
	getalpha!(alp, uar)
        for i=2:N-1
                unew = uar[i, :] + lambda*(flux(alp, uar, i-1, lambda) - flux(alp, uar, i, lambda))
                pc[i] = PUEuler(unew) - PUEuler(uar[i, :]) + lambda*(EfluxDD(alp, uar, i, lambda) - EfluxDD(alp, uar, i-1,lambda))
        end
        return alp, pc
end

function EulerDD!(du, u, p, t)
	N, K = size(u)
	dx = p[1]
	dt = p[2]
	lambda = dt/dx
	alp = zeros(N, 1)
	getalpha!(alp, u)
        alpp = LowAlpBoundpp(u, lambda, pfunc)
        alprho = LowAlpBoundpp(u, lambda, rhofunc)
        alp = max.(alp, alprho, alpp)

	for i=2:N-1
		du[i, :] = (flux(alp, u, i-1, lambda) - flux(alp, u, i, lambda))/dx               
	end
	du[1, :] .= 0.0
	du[N, :] .= 0.0
end

function EulerDDsolid!(du, u, par, t)
        N, K = size(u)
        dx = par[1]
        dt = par[2]
        lambda = dt/dx
           # enforce reflecting boundary
        pl = p(u[2, :])
        vl = u[2, 2] / u[2, 1]
        u[1, :] = cons_var(u[2, 1], -vl, pl)

        pr = p(u[N-1, :])
        vr = u[N-1, 2] / u[N-1, 1]
        u[N, :] = cons_var(u[N-1, 1], -vr, pr)

        alp = zeros(N, 1)
        getalpha!(alp, u)
	alpp = LowAlpBoundpp(u, lambda, pfunc)
        alprho = LowAlpBoundpp(u, lambda, rhofunc)
        alp = max.(alp, alprho, alpp)
        
	for i=2:N-1
                du[i, :] = (flux(alp, u, i-1, lambda) - flux(alp, u, i, lambda))/dx
        end
        du[1, :] .= 0.0
        du[N, :] .= 0.0
end


function NT(u0, tmax; TG = false, lcfl = CFL, tinteg = SSPRK33())
	if TG == false
                u0ar = u0(xgridtest)
                locdx = dxtest
                locdt = lcfl*dxtest
        else
                u0ar = u0(TG)
                locdx = TG[2]-TG[1]
                locdt = lcfl*locdx
        end
	tspan = (0, tmax)
	prob = ODEProblem(EulerDD!, u0ar, tspan, [locdx, locdt])
	sol = solve(prob, tinteg, dt=locdt, progress=true)
	return sol
end


function NTsolid(u0, tmax; TG = false, lcfl = CFL)
        if TG == false
                u0ar = u0(xgridtest)
                locdx = dxtest
                locdt = lcfl*dxtest
        else
                u0ar = u0(TG)
                locdx = TG[2]-TG[1]
                locdt = lcfl*locdx
        end
        tspan = (0, tmax)
        prob = ODEProblem(EulerDDsolid!, u0ar, tspan, [locdx, locdt])
        sol = solve(prob, SSPRK33(), dt=locdt, progress=true)
        return sol
end
