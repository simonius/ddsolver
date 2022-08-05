using LaTeXStrings
include("enosolver.jl")
include("etest3.jl")

### ShuOsher 6 Pictures for discrete entropy stable scheme

solDE = NTdeHO22(shuosh6a, 2.0)
refsol = CalcRefNP(u0f=shuosh6a, tmax=2.0)

du = zeros(Ntest, 3)

alp, pl, pr, sl, sr =  EulerDEHO22!(du, solDE(2.0), (dxtest, dttest, false), 0.0) 
plot(xgridcont, refsol(2.0)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
plot!(xgridtest, alp, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, solDE(2.0)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "DELFT at "*string(Ntest))
savefig("pics/shuosh6DE"*string(Ntest)*"rho.pdf")

plot(xgridcont, refsol(2.0)[:, 2] ./ refsol(2.0)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont), legend=:topleft)
plot!(xgridtest, alp, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, solDE(2.0)[:, 2] ./ solDE(2.0)[:, 1], markershape=:x, color=:black, ylabel = L"v", label = "DELFT at "*string(Ntest))
savefig("pics/shuosh6DE"*string(Ntest)*"v.pdf")


pf(i) = p(refsol(2.0)[i, :])
plot(xgridcont, pf.(collect(1:Ncont)), color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
plot!(xgridtest, alp, color=:black, ls=:dot, label = L"\alpha")
pf(i) = p(solDE(2.0)[i, :])
scatter!(xgridtest, pf.(collect(1:Ntest)), markershape=:x, color=:black, ylabel = L"p", label = "DELFT at "*string(Ntest))
savefig("pics/shuosh6DE"*string(Ntest)*"p.pdf")


#### ShuOsher 6b Pictures

solDE = NTdeHO22(shuosh6b, 1.3)
#println("Warn: deactivated reference calc")
refsol = CalcRefNP(u0f=shuosh6b)

alp, pl, pr, sl, sr =  EulerDEHO22!(du, solDE(1.2), (dxtest, dttest, false), 0.0)

plot(xgridcont, refsol(1.3)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont), legend=:topleft)
plot!(xgridtest, alp, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, solDE(1.3)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "DELFT at "*string(Ntest))
savefig("pics/shuosh6bDE"*string(Ntest)*".pdf")

#### ShuOsher 8

solDE = NTdeHO22(shuosh8, 1.8)
#println("Warn: deactivated reference calc")
refsol = CalcRefNP(u0f=shuosh8)

alp, pl, pr, sl, sr =  EulerDEHO22!(du, solDE(1.8), (dxtest, dttest, false), 0.0)

plot(xgridcont, refsol(1.8)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
plot!(xgridtest, alp, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, solDE(1.8)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "DELFT at "*string(Ntest))
savefig("pics/shuosh8DE"*string(Ntest)*".pdf")


#### Toro 123

solDE = NTdeHO22(toro123, 1.2)
refsol = CalcRefNP(u0f=toro123)

alp, pl, pr, sl, sr =  EulerDEHO22!(du, solDE(1.2), (dxtest, dttest, false), 0.0)

plot(xgridcont, refsol(1.2)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
plot!(xgridtest, alp, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, solDE(1.2)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "DELFT at "*string(Ntest))
savefig("pics/toro123DE"*string(Ntest)*".pdf")

pf(i) = p(refsol(1.2)[i, :])
plot(xgridcont, pf.(collect(1:Ncont)), color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
pf(i) = p(solDE(1.2)[i, :])
scatter!(xgridtest, pf.(collect(1:Ntest)), markershape=:x, color=:black, ylabel = L"p", label = "DELFT at "*string(Ntest))
savefig("pics/toro123DE"*string(Ntest)*"p.pdf")

####################### Woodward Colella Hoes


solDE = NTdeHOrefl22(shuosh7, 0.38, cfl = 0.01)
refsol = CalcRefNPsolid(u0f=shuosh7, tmax = 0.38, locfl=0.01)

alp, pl, pr, sl, sr =  EulerDEHO22!(du, solDE(0.38), (dxtest, dxtest*0.01, true), 0.0)

plot(xgridcont, refsol(0.38)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))

plot!(xgridtest, alp, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, solDE(0.38)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "DELFT at "*string(Ntest))
savefig("pics/shuosh7DE"*string(Ntest)*"rho.pdf")

plot(xgridcont, refsol(0.38)[:, 2]./refsol(0.38)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont), legend=:topleft)

plot!(xgridtest, alp, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, solDE(0.38)[:, 2] ./solDE(0.38)[:, 1], markershape=:x, color=:black, ylabel = L"v", label = "DELFT at "*string(Ntest))
savefig("pics/shuosh7DE"*string(Ntest)*"v.pdf")

pf(i) = p(refsol(0.38)[i, :])
plot(xgridcont, pf.(collect(1:Ncont)), color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))

plot!(xgridtest, alp, color=:black, ls=:dot, label = L"\alpha")

pf(i) = p(solDE(0.38)[i, :])
scatter!(xgridtest, pf.(collect(1:Ntest)), markershape=:x, color=:black, ylabel = L"p", label = "DELFT at "*string(Ntest))
savefig("pics/shuosh7DE"*string(Ntest)*"p.pdf")


