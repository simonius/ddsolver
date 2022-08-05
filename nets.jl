# File contains the definition of our Network
# The entire programm uses the Flux ML library for networks and training.

using Flux
s=80
pdrop = 0.2
act = elu

id(x) = x
actend = id

layer1 = Dense(4*2*snapw, s, act)                
layer1b = Dropout(pdrop)                                             
layer2 = Dense(s, s, act)
layer2b = Dropout(pdrop)			
layer3 = Dense(s, s, act)
layer3b = Dropout(pdrop)
layer4 = Dense(s, s, act)
layer4b = Dropout(pdrop)
layer5 = Dense(s, s, act)
layer5b = Dropout(pdrop)
layer6 = Dense(s,1, actend)                  

model = Chain(layer1, layer1b, layer2, layer2b, layer3, layer3b, layer4, layer4b, layer5, layer5b, layer6)

