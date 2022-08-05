# File contains routines for training. 
# Consult also rpt2 which calls these routines

include("param.jl")

using HDF5

using Flux: train!
using Flux: Optimise.Optimiser
using BSON: @save
using Plots
pyplot()
include("nets.jl")
include("fluxes.jl")

tdatafn = "hdf5test.h5"


rxdata = h5read(tdatafn, "u")
rydata = h5read(tdatafn, "alpha")

function fam(u, alp)
        R, L, K, C = size(rxdata)

        xdata = zeros(Float32, 8*snapw,K, L, R)
        ydata = zeros(Float32, 1, K, L, R)
        # adding pressures
        for r = 1:R #tests
        for l=1:L #timesteps
                for k=1:K #koordinates
                        for i=1:snapw
	                        xdata[i, k, l, r] = rxdata[r,l, mod1(k-snapw+i, K), 1]
                                xdata[snapw+i, k, l, r] = rxdata[r, l, mod1(k+i, K), 1]
                                xdata[2*snapw+i,k, l, r] = rxdata[r, l, mod1(k-snapw+i, K), 2]
                                xdata[3*snapw+i, k, l, r] = rxdata[r, l, mod1(k+i, K), 2]
                                xdata[4*snapw+i, k, l, r] = rxdata[r, l, mod1(k-snapw+i, K), 3]
                                xdata[5*snapw+i, k, l, r] = rxdata[r, l, mod1(k+i, K), 3]

                        end
	                for s=1:2*snapw
		                xdata[6*snapw+s, k, l, r] = p([xdata[s, k, l, r], xdata[2*snapw+s,k,  l, r], xdata[4*snapw+s, k, l, r]])
	                end
                
                        ydata[1, k, l, r] = rydata[r, l, k, 1]
                end
        end
        end

        return reshape(xdata, (8*snapw, L*R*K)), reshape(ydata, (1, L*R*K))

end

xdata, ydata = fam(rxdata, rydata)
K, L = size(xdata)

ndata = 1024*1024

rindexes = mod1.(rand(Int, 1024), L)
xdataTest = xdata[:, rindexes]
ydataTest = ydata[:, rindexes]

rindexes = mod1.(rand(Int, ndata), L)
xdataTrain = xdata[:, rindexes]
ydataTrain = ydata[:, rindexes]


data = [(xdata, ydata)] 

sqnorm(x) = sum(abs2, x)
loss(x, y) = Flux.Losses.mse(model(x), y)

function trainnet!(model, loss, xdata, ydata, fac = 1, runs = 100, bs=16, opt=ADAM())
	K, L = size(xdata)
        lossar = zeros(runs)
	for epoch=1:runs
        	println("Epoch ", epoch)
        	for k=1:floor(Int, fac*L/bs)
			indexes = mod1.(rand(Int, bs), L)
			data = [(xdata[:, indexes], ydata[:, indexes])]
                	train!(loss, Flux.params(model), data, opt)
		end
		lossar[epoch] = loss(xdataTest, ydataTest)
       		println("Loss ", lossar[epoch])
	end

	weights = Flux.params(model);
	@save wdatafn weights

	plot(lossar, xlabel="epochs", ylabel="loss")
	return lossar
end

function diagnosenet(model, xdata, ydata, ind)               # some plots to test convergence
        display(scatter(ydata[ind, :], model(xdata)[ind, :], xlabel="reference", ylabel = "model"))
        println("loss is: ", loss(xdata, ydata))
end

