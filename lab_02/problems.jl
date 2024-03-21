##
begin
    using Plots
    using Statistics
    using Random
end

## problem 2.1
sampling_rate = 1000
sample_step = 1 / sampling_rate
first_sample_time = 0.25
sample_count = 256

A = 2
f₀ = 25
ϕ = π / 4

t = range(; start=first_sample_time, step=sample_step, length=sample_count)
signal = A * cos.(2π * f₀ * t .+ ϕ)

plot(signal)

## problem 2.2
sampling_rate = 2048
sample_step = 1 / sampling_rate
first_sample_time = 5
last_sample_time = 10

A = 0.25
f₀ = π / 2
ϕ = π

t = range(; start=first_sample_time, stop=last_sample_time, step=sample_step)
signal = A * exp.(im * (2π * f₀ * t .+ ϕ))
p1 = plot(signal)
plot(plot(real(signal)), plot(imag(signal)), layout=(2, 1), title=["real" "imaginary"], label=["re" "im"])

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


end