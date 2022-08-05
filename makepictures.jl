using LaTeXStrings
include("enosolver.jl")

###################### ShuOsher 6 Pictures for HOES steered blending
include("hoestest.jl")
suHOES = NThoes(shuosh6a, 2.0)
#println("Warn: deactivated reference calc")
refsol = CalcRefNP(u0f=shuosh6a, tmax=2.0)


plot(xgridcont, refsol(2.0)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
plot!(xgridtest, hoesalpha(suHOES(2.0)), color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, suHOES(2.0)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "PALFT at "*string(Ntest))
savefig("pics/shuosh6hoes"*string(Ntest)*"rho.pdf")

plot(xgridcont, refsol(2.0)[:, 2] ./ refsol(2.0)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont), legend=:topleft)
plot!(xgridtest, hoesalpha(suHOES(2.0)), color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, suHOES(2.0)[:, 2] ./ suHOES(2.0)[:, 1], markershape=:x, color=:black, ylabel = L"v", label = "PALFT at "*string(Ntest))
savefig("pics/shuosh6hoes"*string(Ntest)*"v.pdf")


pf(i) = p(refsol(2.0)[i, :])
plot(xgridcont, pf.(collect(1:Ncont)), color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
plot!(xgridtest, hoesalpha(suHOES(2.0)), color=:black, ls=:dot, label = L"\alpha")
pf(i) = p(suHOES(2.0)[i, :])
scatter!(xgridtest, pf.(collect(1:Ntest)), markershape=:x, color=:black, ylabel = L"p", label = "PALFT at "*string(Ntest))
savefig("pics/shuosh6hoes"*string(Ntest)*"p.pdf")

#####6b
suHOES = NThoes(shuosh6b, 1.3)
#println("Warn: deactivated reference calc")
refsol = CalcRefNP(u0f=shuosh6b)

plot(xgridcont, refsol(1.3)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont), legend=:topleft)
plot!(xgridtest, hoesalpha(suHOES(1.3)), color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, suHOES(1.3)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "PALFT at "*string(Ntest))
savefig("pics/shuosh6bhoes"*string(Ntest)*".pdf")




####################### ShuOsher 8 Pictures for HOES steered blending
include("hoestest.jl")
suHOES = NThoes(shuosh8, 1.8)
#println("Warn: deactivated reference calc")
refsol = CalcRefNP(u0f=shuosh8)

plot(xgridcont, refsol(1.8)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
plot!(xgridtest, hoesalpha(suHOES(1.8)), color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, suHOES(1.8)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "PALFT at "*string(Ntest))
savefig("pics/shuosh8hoes"*string(Ntest)*".pdf")



######################### ShuOsher6 Pictures for data driven blending
alps = zeros(Ntest, 1)
include("tests.jl")

sudd = NT(shuosh6a, 2.0)
refsol = CalcRefNP(u0f=shuosh6a, tmax = 2.0)



plot(xgridcont, refsol(2.0)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
getalpha!(alps, sudd(2.0))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, sudd(2.0)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "DDLFT at "*string(Ntest))
savefig("pics/shuosh6dd"*string(Ntest)*"rho.pdf")

plot(xgridcont, refsol(2.0)[:, 2]./refsol(2.0)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont), legend=:topleft)
getalpha!(alps, sudd(2.0))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, sudd(2.0)[:, 2] ./sudd(2.0)[:, 1], markershape=:x, color=:black, ylabel = L"v", label = "DDLFT at "*string(Ntest))
savefig("pics/shuosh6dd"*string(Ntest)*"v.pdf")

pf(i) = p(refsol(2.0)[i, :])
plot(xgridcont, pf.(collect(1:Ncont)), color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
getalpha!(alps, sudd(2.0))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")

pf(i) = p(sudd(2.0)[i, :])
scatter!(xgridtest, pf.(collect(1:Ntest)), markershape=:x, color=:black, ylabel = L"p", label = "DDLFT at "*string(Ntest))
savefig("pics/shuosh6dd"*string(Ntest)*"p.pdf")

####6b
sudd = NT(shuosh6b, 1.3)
refsol = CalcRefNP(u0f=shuosh6b)
plot(xgridcont, refsol(1.3)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont), legend=:topleft)
getalpha!(alps, sudd(1.3))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, sudd(1.3)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "DDLFT at "*string(Ntest))
savefig("pics/shuosh6bdd"*string(Ntest)*".pdf")


######################## ShuOsher8 Pictures for data driven blending
include("tests.jl")

sudd = NT(shuosh8, 1.8)
refsol = CalcRefNP(u0f=shuosh8)



getalpha!(alps, sudd(1.8))
plot(xgridcont, refsol(1.8)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, sudd(1.8)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "DDLFT at "*string(Ntest))
savefig("pics/shuosh8dd"*string(Ntest)*".pdf")

####################### Totoro test 123 DD
sudd = NT(toro123, 1.2)
refsol = CalcRefNP(u0f=toro123)

getalpha!(alps, sudd(1.2))
alpp = LowAlpBound(sudd(1.2), CFL, pfunc)
alprho = LowAlpBound(sudd(1.2), CFL, rhofunc)
alps = max.(alps, alprho, alpp)

plot(xgridcont, refsol(1.2)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, sudd(1.2)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "DDLFT at "*string(Ntest))
savefig("pics/toro123dd"*string(Ntest)*".pdf")

pf(i) = p(refsol(1.2)[i, :])
plot(xgridcont, pf.(collect(1:Ncont)), color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
pf(i) = p(sudd(1.2)[i, :])
scatter!(xgridtest, pf.(collect(1:Ntest)), markershape=:x, color=:black, ylabel = L"p", label = "DDLFT at "*string(Ntest))
savefig("pics/toro123dd"*string(Ntest)*"p.pdf")



####################### Totoro test 123 Hoes
suHOES = NThoes(toro123, 1.2)
refsol = CalcRefNP(u0f=toro123)

plot(xgridcont, refsol(1.2)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
plot!(xgridtest, hoesalpha(suHOES(1.2)), color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, suHOES(1.2)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "PALFT at "*string(Ntest))
savefig("pics/toro123hoes"*string(Ntest)*".pdf")

pf(i) = p(refsol(1.2)[i, :])
plot(xgridcont, pf.(collect(1:Ncont)), color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
pf(i) = p(suHOES(1.2)[i, :])
scatter!(xgridtest, pf.(collect(1:Ntest)), markershape=:x, color=:black, ylabel = L"p", label = "PALFT at "*string(Ntest))
savefig("pics/toro123hoes"*string(Ntest)*"p.pdf")

###################### Woodward Colella DD
include("tests.jl")

sudd = NTsolid(shuosh7, 0.38, lcfl = 0.01)
refsol = CalcRefNPsolid(u0f=shuosh7, tmax = 0.38, locfl=0.01)


plot(xgridcont, refsol(0.38)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
getalpha!(alps, sudd(0.38))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, sudd(0.38)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "DDLFT at "*string(Ntest))
savefig("pics/shuosh7dd"*string(Ntest)*"rho.pdf")

plot(xgridcont, refsol(0.38)[:, 2]./refsol(0.38)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont), legend=:topleft)
getalpha!(alps, sudd(0.38))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, sudd(0.38)[:, 2] ./sudd(0.38)[:, 1], markershape=:x, color=:black, ylabel = L"v", label = "DDLFT at "*string(Ntest))
savefig("pics/shuosh7dd"*string(Ntest)*"v.pdf")

pf(i) = p(refsol(0.38)[i, :])
plot(xgridcont, pf.(collect(1:Ncont)), color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
getalpha!(alps, sudd(0.38))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")

pf(i) = p(sudd(0.38)[i, :])
scatter!(xgridtest, pf.(collect(1:Ntest)), markershape=:x, color=:black, ylabel = L"p", label = "DDLFT at "*string(Ntest))
savefig("pics/shuosh7dd"*string(Ntest)*"p.pdf")

####################### Woodward Colella Hoes

include("tests.jl")

suhoes = NThoessolid(shuosh7, 0.38, locfl = 0.01)
refsol = CalcRefNPsolid(u0f=shuosh7, tmax = 0.38, locfl=0.01)


plot(xgridcont, refsol(0.38)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
getalpha!(alps, suhoes(0.38))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, suhoes(0.38)[:, 1], markershape=:x, color=:black, ylabel = L"\rho", label = "PALFT at "*string(Ntest))
savefig("pics/shuosh7hoes"*string(Ntest)*"rho.pdf")

plot(xgridcont, refsol(0.38)[:, 2]./refsol(0.38)[:, 1], color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont), legend=:topleft)
getalpha!(alps, suhoes(0.38))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")
scatter!(xgridtest, suhoes(0.38)[:, 2] ./suhoes(0.38)[:, 1], markershape=:x, color=:black, ylabel = L"v", label = "PALFT at "*string(Ntest))
savefig("pics/shuosh7hoes"*string(Ntest)*"v.pdf")

pf(i) = p(refsol(0.38)[i, :])
plot(xgridcont, pf.(collect(1:Ncont)), color=:black, xlabel = L"x", label = "ENO2 at "*string(Ncont))
getalpha!(alps, suhoes(0.38))
plot!(xgridtest, alps, color=:black, ls=:dot, label = L"\alpha")

pf(i) = p(suhoes(0.38)[i, :])
scatter!(xgridtest, pf.(collect(1:Ntest)), markershape=:x, color=:black, ylabel = L"p", label = "PALFT at "*string(Ntest))
savefig("pics/shuosh7hoes"*string(Ntest)*"p.pdf")




