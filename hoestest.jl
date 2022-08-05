# File contains implementation of the PA operator based scheme and all needed functions to perform tests.

using DifferentialEquations
using Plots
pyplot()

include("param.jl")
include("hotimefluxes.jl")
include("hoes.jl")

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


function hoesalpha(u)
	N, K = size(u)
	par = zeros(N)
	no, np = zeros(3), 0
	for i=1:N
		par[i] = p(u[i, :])
		no = no + abs.(u[i, :])./N
		np = np + abs(par[i])./N
	end
	h1, hn1 = hoes(u[:, 1])
	h2, hn2 = hoes(u[:, 2])
	h3, hn3 = hoes(u[:, 3])
	hp, hnp = hoes(par)
	fac = 10.0
	return min.(1, max.(0, fac*max.(h1/(no[1]+hn1+eps), h2/(no[2]+hn2+eps), h3/(no[3] + hn3+eps), hp/(np+hnp+eps))))
end

function fluxh(alp, u, k, lambda)
	width=order+1
	width=2
        N, K = size(u)
        if k < width || k > N-width
                return hEulerLLF(u[k, :], u[k+1, :])
        else
                return eulerHOconsflux(u[k-width+1:k+width, :], width, alp[k])

        end
end


function CheckProdHoes(uar, lambda)
        N, K = size(uar)
        pc =  zeros(N)
        alp = hoesalpha(uar)
	alp[3:end-2] = max.(0.33.*alp[1:end-4], 0.66.*alp[2:end-3], alp[3:end-2], 0.66.*alp[4:end-1], 0.33.*alp[5:end])
        for i=2:N-1
                unew = uar[i, :] + lambda*(fluxh(alp, uar, i-1, lambda) - fluxh(alp, uar, i, lambda))
                pc[i] = PUEuler(unew) - PUEuler(uar[i, :]) + lambda*(Efluxh(alp, uar, i, lambda) - Efluxh(alp, uar, i-1,lambda))
        end
        return alp, pc
end

function EulerHoes!(du, u, p, t)
        N, K = size(u)
        dx = p[1]
        dt = p[2]
        lambda = dt/dx
        alp = hoesalpha(u)
        alpp = LowAlpBoundpp(u, lambda, pfunc)
        alprho = LowAlpBoundpp(u, lambda, rhofunc)
        alp = max.(alp, alprho, alpp)

	alp[3:end-2] = max.(0.33.*alp[1:end-4], 0.66.*alp[2:end-3], alp[3:end-2], 0.66.*alp[4:end-1], 0.33.*alp[5:end])
	for i=2:N-1
                du[i, :] = (fluxh(alp, u, i-1, lambda) - fluxh(alp, u, i, lambda))/dx
        end
        du[1, :] .= 0.0
        du[N, :] .= 0.0
end

function EulerHoessolid!(du, u, para, t)
        N, K = size(u)
        dx = para[1]
        dt = para[2]
        lambda = dt/dx
        # enforce reflecting boundary
        pl = p(u[2, :])
        vl = u[2, 2] / u[2, 1]
        u[1, :] = cons_var(u[2, 1], -vl, pl)
        
        pr = p(u[N-1, :])
        vr = u[N-1, 2] / u[N-1, 1]
        u[N, :] = cons_var(u[N-1, 1], -vr, pr)


        alp = hoesalpha(u)
        alp[3:end-2] = max.(0.33.*alp[1:end-4], 0.66.*alp[2:end-3], alp[3:end-2], 0.66.*alp[4:end-1], 0.33.*alp[5:end]) #.+ (50/N)^2
        for i=2:N-1
                du[i, :] = (fluxh(alp, u, i-1, lambda) - fluxh(alp, u, i, lambda))/dx
        end
        du[1, :] .= 0.0
        du[N, :] .= 0.0
end


function hoesalpha(u)
        N, K = size(u)
        par = zeros(N)
        no, np = zeros(3), 0
        for i=1:N
                par[i] = p(u[i, :])
                no = no + abs.(u[i, :])./N
                np = np + abs(par[i])./N
        end
        h1, hn1 = hoes(u[:, 1])
        h2, hn2 = hoes(u[:, 2])
        h3, hn3 = hoes(u[:, 3])
        hp, hnp = hoes(par)
        fac = 10.0
        return min.(1, max.(0, fac*max.(h1/(no[1]+hn1+eps), h2/(no[2]+hn2+eps), h3/(no[3] + hn3+eps), hp/(np+hnp+eps))))
end


function NThoesVG(u0, tmax, TG, CFL)
        tspan = (0, tmax)
	locdx = TG[2]-TG[1]
        locdt = CFL*locdx
        u0ar = u0(TG)
        prob = ODEProblem(EulerHoes!, u0ar, tspan, [locdx, locdt])
        sol = solve(prob, SSPRK104(), dt=locdt, progess=true)
        return sol
end

function NThoes(u0, tmax; locfl=CFL)
        tspan = (0, tmax)
        u0ar = (u0(xgridtest .- dxtest./2) .+ u0(xgridtest .+ dxtest./2)) ./ 2
        locdx = dxtest
        locdt = locfl*dxtest
        prob = ODEProblem(EulerHoes!, u0ar, tspan, [locdx, locdt])
        sol = solve(prob, SSPRK33(), dt=locdt, progress=true)
        return sol
end

function NThoessolid(u0, tmax; locfl=CFL)
  	tspan = (0, tmax)
        u0ar = (u0(xgridtest .- dxtest./2) .+ u0(xgridtest .+ dxtest./2)) ./ 2
        locdx = dxtest
        locdt = locfl*dxtest
        prob = ODEProblem(EulerHoessolid!, u0ar, tspan, [locdx, locdt])
       	sol = solve(prob, SSPRK33(), dt=locdt, progress=true)
        return sol
end



