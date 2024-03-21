##
begin
    using CairoMakie
    using Plots
    using Statistics
    using Random
end

## problem 2.1
begin
    sampling_rate = 1000
    sample_step = 1 / sampling_rate
    amplitude = 2
    frequency = 25
    phase_shift = π / 4
    first_sample_time = 0.25
    sample_count = 256
    t = range(; start=first_sample_time, step=sample_step, length=sample_count)
    signal = amplitude * cos.((2 * π * frequency * t) .+ phase_shift)
end

## problem 2.2
begin
    sampling_rate = 2048
    sample_step = 1 / sampling_rate
    amplitude = 0.25
    frequency = π / 2
    phase_shift = π
    first_sample_time = 5
    last_sample_time = 10
    t = range(; start=first_sample_time, stop=last_sample_time, step=sample_step)
    signal = amplitude * exp.((im * 2 * π * frequency * t) .+ phase_shift)
end

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