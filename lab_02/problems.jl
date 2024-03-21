##
begin
    using Plots
    using Statistics
    using LinearAlgebra
    using Random
    using Waveforms
end

## problem 2.1
sampling_rate = 1000
sample_step = 1 / sampling_rate
first_sample_time = 0.25
sample_count = 256

A = 2
f₀ = 25
ϕ = π / 4

t = range(;
    start=first_sample_time,
    step=sample_step,
    length=sample_count
)
signal = A * cos.(2π * f₀ * t .+ ϕ)

plot(t, signal)

## problem 2.2
sampling_rate = 2048
sample_step = 1 / sampling_rate
first_sample_time = 5
last_sample_time = 10

A = 0.25
f₀ = π / 2
ϕ = π

t = range(;
    start=first_sample_time,
    stop=last_sample_time,
    step=sample_step
)
signal = A * exp.(im * (2π * f₀ * t .+ ϕ))

p1 = plot(signal)
plot(
    plot(t, real(signal)),
    plot(t, imag(signal)),
    layout=(2, 1),
    title=["real" "imaginary"],
    label=["re" "im"]
)

## problem 2.3
function white_noise(n, power)
    noise = sqrt(power) * randn(n)
    return noise
end
signal = white_noise(1000, 0.25)
var(signal)

## problem 2.4
function complex_white_noise(n, power)
    noise = sqrt(power) * (randn(ComplexF64, n))
    return noise
end
signal = complex_white_noise(1000, 3)
var(signal)

## problem 2.5
function cw_rectangular(T, t)
    impulse_value = 1 / T


end

## problem 2.9
function ramp_wave(x)
    output = 2 * rem(x, 1, RoundNearest)
    return output
end

x = 0:0.001:2
plot(x, ramp_wave.(x))

## problem 2.10
function sawtooth_wave(x)
    output = -2 * rem(x, 1, RoundNearest)
    return output
end

x = 0:0.001:2
plot(x, sawtooth_wave.(x))

## problem 2.11
function triangle_wave(x)
    #FIXME: conflicting requirements
end

x = 0:0.001:2
plot(x, triangle_wave.(x))

## problem 2.12
function square_wave(x)
    ifelse(rem(x, 1) < 0.5, 1, -1)
end

trianglewave1(x)

x = 0:0.001:2
plot(x, square_wave.(x))

## problem 2.13
function pulse_wave(x, ρ)
    ifelse(rem(x, 1) < ρ, 1, 0)
end

x = 0:0.001:2
plot(x, pulse_wave.(x, 0.2))