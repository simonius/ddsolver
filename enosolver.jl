# File contains the eno2/MUSCL solvers used to calculate reference solutions

using DifferentialEquations

include("eno.jl")
include("fluxes.jl")
include("param.jl")

function EulerRef(du::Matrix{Float64}, u::Matrix{Float64}, p::Vector{Float64}, t::Float64)
	dx = p[1]
	N = size(u)[1]
	ul, ur = MUSCLMat(u)
	for i = 1:N
		du[i, :] = (hEulerLLF(ur[mod1(i-1, N), :], ul[mod1(i, N), :]) - hEulerLLF(ur[mod1(i, N), :], ul[mod1(i+1, N), :]))/dx
	end
end

function EulerRefNP(du::Matrix{Float64}, u::Matrix{Float64}, p::Vector{Float64}, t::Float64)
        dx = p[1]
        N::Int64 = size(u)[1]
        ul, ur = MUSCLMat(u)
	for i = 3:N-2
                du[i, :] = (hEulerLLF(ur[i-1, :], ul[i, :]) - hEulerLLF(ur[i, :], ul[i+1, :]))/dx
	end
	du[1:2, :] .= 0
        du[N-1:N, :] .= 0
end

function EulerRefNPsolid(du::Matrix{Float64}, u::Matrix{Float64}, par::Vector{Float64}, t::Float64)
        dx = par[1]
        N::Int64 = size(u)[1]
         # enforce reflecting boundary
        pl = p(u[3, :])
        vl = u[3, 2] / u[3, 1]
        u[1, :] = cons_var(u[3, 1], -vl, pl)
        u[2, :] = u[1, :]

        pr = p(u[N-2, :])
        vr = u[N-2, 2] / u[N-2, 1]
        u[N, :] = cons_var(u[N-2, 1], -vr, pr)
        u[N-1, :] = u[N, :]

        ul, ur = MUSCLMat(u)

        for i = 3:N-2
                du[i, :] = (hEulerLLF(ur[i-1, :], ul[i, :]) - hEulerLLF(ur[i, :], ul[i+1, :]))/dx
        end
        du[1:2, :] .= 0
        du[N-1:N, :] .= 0
end


function EulerRefNPE(du::Matrix{Float64}, u::Matrix{Float64}, p::Vector{Float64}, t::Float64)
        dx = p[1]
        N::Int64 = size(u)[1]
        ul, ur = ENO2Mat(u)

        for i = 2:N-1
                du[i, :] = (hEulerLLF(ur[i-1, :], ul[i, :]) - hEulerLLF(ur[i, :], ul[i+1, :]))/dx
        end
        du[1, :] .= 0
        du[N, :] .= 0
end


function CalcRef(u0f=u0, xg=xgridcont)
	u0ar = u0f(xg)
	p = zeros(1)
	p[1] = dxcont
	prob = ODEProblem(EulerRef, u0ar, t, p)
	sol = solve(prob, SSPRK33(), dt=dtcont, saveat=4*dt)
	return sol
end

function CalcRefNP(;u0f = shuosh6, tmax=1.8, xg = xgridcont, locfl=CFL)
	u0ar = u0f(xg)
	p = zeros(1)
	p[1] = xg[2]-xg[1]
	locdt = locfl*p[1]
        prob = ODEProblem(EulerRefNP, u0ar, tmax, p)
        sol = solve(prob, SSPRK22(), dt=locdt)
	return sol
end

function CalcRefNPsolid(;u0f = shuosh6, tmax=1.8, xg = xgridcont, locfl=CFL)
        u0ar = u0f(xg)
        p = zeros(1)
        p[1] = xg[2]-xg[1]
        locdt = locfl*p[1]
        prob = ODEProblem(EulerRefNPsolid, u0ar, tmax, p)
        sol = solve(prob, SSPRK22(), dt=locdt)
        return sol
end


function CalcRefNPE(;u0f = shuosh6, tmax=1.8, xg = xgridcont)
        u0ar = u0f(xg)
        p = zeros(1)
        p[1] = xg[2]-xg[1]
        locdt = 3.0*CFL*p[1]
        prob = ODEProblem(EulerRefNPE, u0ar, tmax, p)
        sol = solve(prob, SSPRK33(), dt=locdt)
        return sol
end

