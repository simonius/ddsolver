using LinearAlgebra

# File contains the Flux funcitons and Riemann solvers for the 1D Euler equations
# u = (rho, rho u, E)

gamma = 1.4

@inline function p(u)
        p = (gamma-1)*(u[3]-0.5*u[2]^2/u[1])
        return p
end

@inline function a(u, p)
        a = sqrt(abs(gamma*p / u[1]))
        return a
end

@inline function fEuler(u, p)
        fEuler1 = u[2]
        fEuler2 = u[2]^2/u[1] + p
        fEuler3 = u[2]/u[1]*(u[3]+p)
        return [fEuler1, fEuler2, fEuler3]
end

@inline function fEuler(u)
        return fEuler(u, p(u))
end

@inline function hEulerLF(ul, ur, lambda)
        hEulerLF = zeros(3)
        pl = p(ul)
        pr = p(ur)
        hEulerLF = 0.5*(fEuler(ul, pl) + fEuler(ur, pr) - (ur-ul)*lambda)
        return hEulerLF
end

@inline function hEulerLLF(ul, ur)
	pl = p(ul)
        pr = p(ur)
        amax = max(a(ul, pl), a(ur, pr))
        lambound=  max(abs(ul[2]/ul[1]), abs(ur[2]/ur[1])) + amax
       	return hEulerLF(ul, ur, lambound) 
end

@inline meanv(ul, ur) =  0.5*(ul + ur)
# out of Roe 2009 Appendix A and B
@inline function logmean(ul, ur)
        eps = 10^(-2)
        zeta = ul/ur
        f = (zeta - 1)/(zeta + 1)
        u = f^2
        if u < eps
                F = 1 + u/3 + u^2/ 5 + u^3/7
        else
                F = log(abs(zeta))/(2*f)
        end
        return (ul + ur)/(2*F)
end

# out of ARBITRARILY HIGH-ORDER ACCURATE ENTROPY STABLEESSENTIALLY NONOSCILLATORY SCHEMES FOR SYSTEMS OFCONSERVATION LAWS
# By Fjordhlm, Mishra and Tadmor
@inline function hEulerCons(ul, ur)
        hEulerCons = zeros(3)
        pl, pr = p(ul), p(ur)
        zl = sqrt(abs(ul[1]/pl))*[1.0, ul[2] / ul[1], pl]
        zr = sqrt(abs(ur[1]/pr))*[1.0, ur[2] / ur[1], pr]

        rhomean =meanv(zl[1], zr[1])*logmean(zl[3], zr[3])
        umean = meanv(zl[2], zr[2])/meanv(zl[1], zr[1])
        p1mean = meanv(zl[3], zr[3])/meanv(zl[1], zr[1])
        p2mean = (gamma + 1)/gamma*logmean(zl[3], zr[3])/logmean(zl[1], zr[1])/2 + (gamma-1)/(2*gamma)*meanv(zl[3], zr[3])/meanv(zl[1], zr[1])
        amean = sqrt(gamma*p2mean/rhomean)
        Hmean = amean^2/(gamma-1) + umean^2/2
        hEulerCons[1] = rhomean*umean
        hEulerCons[2] = rhomean*umean^2 + p1mean
        hEulerCons[3] = rhomean*umean*Hmean
        return hEulerCons
end


# The physical Euler Entropy see Harten (1983) On the symmetric form of systems of conservation laws
# scaled with gamma-1 as in Fjordholm (2012)
function PUEuler(u)
	S = log(p(u)) - gamma*log(u[1])
        return -u[1]*S/(gamma-1)
end

function PFEuler(u)
        S = log(p(u)*u[1]^(-gamma))
        return -u[2]*S/(gamma-1)
end


function GetAlphaVec(k)
        b = zeros(k)
        b[1] = 1.0
        A = zeros(k,k)
        for r=1:k
                A[1, r] = r #fuck i say r, they say 2r
        end
        for s=2:k
                for r = 1:k
                        A[s, r] = r^(2*s-1)
                end
        end
        alpha = A \ b
        return alpha
end


HOcoefs = [GetAlphaVec(1), GetAlphaVec(2), GetAlphaVec(3), GetAlphaVec(4), GetAlphaVec(5), GetAlphaVec(6), GetAlphaVec(7), GetAlphaVec(8)]

@inline function eulerHOconsflux(u, order, alpha = 0.0)
	fnum = hEulerCons
        k = Int(order)

        h = zeros(3)
        for r = 1:k #    p = k, r = i
                hi = zeros(3)
                for s = 0:r-1
                        hi = hi .+ fnum(u[k-s, :],u[k-s+r, :])
                end
		h = h .+ hi.*((HOcoefs[k])[r])
        end
        ul, ur = u[k, :], u[k+1, :]
        h = h + alpha*HOcoefs[k][1] * (hEulerLLF(ul, ur) - hEulerCons(ul, ur))
        return h
end

# gives the euler entropy variables 
# u = [\rho,  u\rho, E]
function eulerEvar(cu)
        rho = cu[1]
        u = cu[2]/cu[1]
        E = cu[3]
        p = (gamma-1)*(E - 0.5*rho*u^2)
        s = log(p) - gamma*log(rho)
        v1 = (gamma-s)/(gamma-1) - (rho*u^2) / (2*p)
        v2 = rho*u/p
        v3 = - rho / p
        return [v1, v2, v3]
end

function eulerEpot(cu)
	return cu[2]
end

# From Tadmor 1987 (4.10b)
function HeulerTRI(ul, ur)
	vl, vr = eulerEvar(ul), eulerEvar(ur)
	g = hEulerCons(ul, ur)
	phil, phir = eulerEpot(ul), eulerEpot(ur)
	return 0.5*(dot(vl + vr, g) - phil - phir)
end

function HEulerLF(ul, ur, lambda)
        HEulerLF = zeros(3)
        pl = p(ul)
        pr = p(ur)
        HEulerLF = 0.5*(PFEuler(ul) + PFEuler(ur) - (PUEuler(ur)-PUEuler(ul))*lambda)
        return HEulerLF
end

function HEulerLLF(ul, ur)
        pl = p(ul)
        pr = p(ur)
        amax = max(a(ul, pl), a(ur, pr))
        lambound=  max(abs(ul[2]/ul[1]), abs(ur[2]/ur[1])) + amax
        return HEulerLF(ul, ur, lambound)
end

function eulerHOconsEflux(u, order, alpha=0.0)
        fnum = HeulerTRI
        k = Int(order)
        falpha(k, r) = (HOcoefs[k])[r]
        h = 0
        for r = 1:k #    p = k, r = i
                hi = 0
                for s = 0:r-1
                        hi = hi .+ fnum(u[k-s, :],u[k-s+r, :])
                end
                h = h .+ hi.*falpha(k, r)
        end
        ul, ur = u[k, :], u[k+1, :]
        h = h + alpha*HOcoefs[k][1] * (HEulerLLF(ul, ur) - fnum(ul, ur))
        return h
end

