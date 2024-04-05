##
include("CPS.jl")
using Plots
using .CPS

## continuous signals
x = -2:0.001:2

##
plot(x, CPS.ramp_wave.(x))
##
plot(x, CPS.sawtooth_wave.(x))
##
plot(x, CPS.triangular_wave.(x))
##
plot(x, CPS.square_wave.(x))
##
plot(x, CPS.pulse_wave.(x, 0.75))



## discrete signals
n = -100:1:100

##
plot(n, CPS.kronecker.(n))
##
plot(n, CPS.heaviside.(n))