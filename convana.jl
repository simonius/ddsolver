# File produces plots of a convergence analysis

include("enosolver.jl")
include("pptest.jl")
include("etest3.jl")
include("tests.jl")
include("hoestest.jl")

nstart = 5
nend = 8
tend = 1.8

Narray = 2 .^collect(1:nend)
HE = zeros(nend)
DDE = zeros(nend)
DEE = zeros(nend)

Nref = 2^13
dxref = (x[2] - x[1]) / Nref
Gref = collect(x[1]:dxref:x[2])[1:Nref] .+ 0.5.*dxref


ilDD = []
ilHS = []
ilTG = []
cflfac = 0.1*Narray[1]
for k=nstart:nend
	N = Narray[k]
	dx = (x[2] - x[1]) / N
	TG = collect(x[1]:dx:x[2])[1:N] .+ 0.5.*dx

	solHOES = NThoesVG(smootheuler, tend,TG, 0.05)
	soldd = NT(smootheuler, tend, TG = TG, lcfl = 0.05, tinteg=SSPRK104())
	solde = NTdeHO22(smootheuler, tend, TG = TG, lcfl=0.025)
        sref = zeros(N, 3)
        projref(x) = uex(x, tend)
        for l = 1:N
                sref[l, :] = projref(TG[l])
        end
	HE[k] =  norm(sref - solHOES(tend), 1)/N
	DEE[k] = norm(sref - solde(tend), 1)/N
	DDE[k] = norm(sref - soldd(tend), 1)/N

	plot(TG, solHOES(tend)[:, 1], label = "HLFT")
	plot!(TG, soldd(tend)[:, 1], label = "DD")
	plot!(TG, solde(tend)[:,1], label = "DE")
	savefig("pics/convana"*string(N)*".pdf")
	push!(ilDD, soldd)
	push!(ilHS, solHOES)
	push!(ilTG, TG)
end

# Measurements
scatter(Narray[nstart:nend], HE[nstart:nend], xlabel = "points", ylabel = "Error", label = "PALFT", color=:black, markershape=:x, xscale=:log2, yscale=:log10)

scatter!(Narray[nstart:nend], DDE[nstart:nend], label = "DDLFT", color=:black, markershape=:+)
scatter!(Narray[nstart:nend], DEE[nstart:nend], label = "DELFT", color=:black, markershape=:o)
Nl = [Narray[nstart], Narray[nend]]
nfac = Nl[2]/Nl[1]
HE1, HEend = HE[nstart], HE[nend]
# Order lines
plot!(Nl, [HE1, HE1/nfac], color=:black, linestyle=:dot, label="First Order")
plot!(Nl, [HE1, HE1/nfac^2], color=:black, linestyle=:dash, label = "Second Order")
plot!(Nl, [HE1, HE1/nfac^3], color=:black, linestyle=:dashdot, label = "Third Order")
plot!(Nl, [HE1, HE1/nfac^4], color=:black, linestyle=:solid, label = "Fourth Order", legend=:bottomleft)
#plot!(Nl, [HE1, HE1/nfac^5], color=:black, linestyle=:dot, label = "Fift Order")
#plot!(Nl, [HE1, HE1/nfac^6], color=:black, linestyle=:dash, label = "Sixt Order", legend=:bottomleft)


savefig("pics/convana.pdf")

