# File contains the training procedure for our NN

include("flick.jl")

lossar = []
nepochs = 20
lossar = vcat(lossar, trainnet!(model, loss, xdata, ydata, 0.1, nepochs, 32, ADAM()))
lossar = vcat(lossar, trainnet!(model, loss, xdata, ydata, 0.3, nepochs, 256, ADAM()))
lossar = vcat(lossar, trainnet!(model, loss, xdata, ydata, 1.0, nepochs, 1024, ADAM()))
lossar = vcat(lossar, trainnet!(model, loss, xdata, ydata, 1.0, nepochs, 4096, ADAM()))
lossar = vcat(lossar, trainnet!(model, loss, xdata, ydata, 1.0, nepochs, 4096*8, ADAM()))
lossar = vcat(lossar, trainnet!(model, loss, xdata, ydata, 1.0, nepochs, 4096*8, ADAM(0.0001)))
lossar = vcat(lossar, trainnet!(model, loss, xdata, ydata, 1.0, nepochs, 4096*8, ADAM(0.00001)))
lossar = vcat(lossar, trainnet!(model, loss, xdata, ydata, 1.0, nepochs, 4096*8, ADAM(0.000001)))
plot(lossar, label = "training curve", xlabel="epochs", ylabel="loss")
savefig("trainloss.pdf")

