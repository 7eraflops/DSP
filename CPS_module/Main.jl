##
begin
    include("CPS.jl")
    using Plots
    using FFTW
    using .CPS
    using BenchmarkTools
end

## continuous signals
t = -2:0.001:2

## Impulses

## problem 2.5
plot(t, CPS.cw_rectangular.(t))
## problem 2.6
plot(t, CPS.cw_triangle.(t))
## problem 2.7
plot(t, CPS.cw_literka_M.(t))
## problem 2.8
plot(t, CPS.cw_literka_U.(t))

## infinite bandwidth waves

## problem 2.9
plot(t, CPS.ramp_wave.(t))
## problem 2.10
plot(t, CPS.sawtooth_wave.(t))
## problem 2.11
plot(t, CPS.triangular_wave.(t))
## problem 2.12
plot(t, CPS.square_wave.(t))
## problem 2.13
plot(t, CPS.pulse_wave.(t, 0.75))
## problem 2.14
plot(t, CPS.impulse_repeater(cos, -π / 4, π / 4).(t))

## finite bandwidth waves 

## problem 2.15
plot(t, CPS.ramp_wave_bl.(t; A=1.0, T=1.0, band=20.0))
## problem 2.16
plot(t, CPS.sawtooth_wave_bl.(t; A=1.0, T=1.0, band=20.0))
## problem 2.17
plot(t, CPS.triangular_wave_bl.(t; A=1.0, T=1.0, band=20.0))
## problem 2.18
plot(t, CPS.square_wave_bl.(t; A=1.0, T=1.0, band=20.0))
## problem 2.19
plot(t, CPS.pulse_wave_bl.(t; ρ=0.75, A=1.0, T=1.0, band=20.0,))
## problem 2.20
plot(t, CPS.impulse_repeater_bl(cos, 0, π / 2, 20).(t))
## problem 2.21
N = length(t)
fs = 1 / 0.001
freqs = [n * fs / N for n in 0:N-1]
signal = CPS.rand_signal_bl(100, 300).(t)
plot(t, signal)
plot(freqs, (2 / N) .* abs.(fft(signal)))

## discrete signals
n = -100:1:100

## problem 2.22
plot(n, CPS.kronecker.(n))
## problem 2.23
plot(n, CPS.heaviside.(n))

## windows
N = 100

## problem 2.24
plot(CPS.rect.(N))
## problem 2.25
plot(CPS.triang.(N))
## problem 2.26
plot(CPS.hanning.(N))
## problem 2.27
plot(CPS.hamming.(N))
## problem 2.28
plot(CPS.blackman.(N))

## signal parameters
signal = randn(1024)

## problem 3.1
CPS.mean(signal)
## problem 3.2
CPS.peak2peak(signal)
## problem 3.3
CPS.energy(signal)
## problem 3.4
CPS.power(signal)
## problem 3.5
CPS.rms(signal)
## problem 3.6
plot(CPS.running_mean(signal, 10))
## problem 3.7
plot(CPS.running_energy(signal, 10))
## problem 3.6
plot(CPS.running_power(signal, 10))

## sampling
## problem 4.5
f = t -> 5cos(t) - 2sin(3t) + cos(10t)
m = -2:0.2:2
s = f.(m)
interpolated_signal = CPS.interpolate(m, s)
t = -2:0.001:2
plot(t, [f.(t), interpolated_signal.(t)]) # y1 - original signal y2 - interpolated signal
plot(t, f.(t) - interpolated_signal.(t)) # error

## quantization
f = t -> 5cos(2π * t) - 2sin(2π * 3t) + cos(2π * 10t)
t = -2:0.001:2
signal = f.(t)
N = 16 # N-bit ADC
L = LinRange(minimum(signal), maximum(signal), 2^N)

## problem 5.1
signal_quantized = CPS.quantize(L).(signal)
noise = signal - signal_quantized
plot(t, [signal, signal_quantized, noise])

## problem 5.2
CPS.SQNR(N)

## problem 5.3
CPS.SNR(CPS.power(signal), CPS.power(noise))

## discrete fourier transform
f = t -> sin(2π * t) + (3 * im) * cos(2π * 5t)
fs = 100
t = -2:(1/fs):2
signal = f.(t)
end_freq = 20
freqs = -end_freq:1:end_freq

## problem 6.3
result = CPS.dtft.(freqs; signal, fs)
plot(t, real.(signal), imag.(signal))
plot(freqs, real.(result), imag.(result))

## problem 6.4-6.5
signal = [0, 1, 0, im, im, 0, 1, 0]
signal_dft = CPS.dft(signal)
signal_idft = CPS.idft(signal_dft)

## problem 6.6-6.7
signal = [0, 1, 0, 1, 1, 0, 1, 0]
signal_rdft = CPS.rdft(signal)
signal_irdft = CPS.irdft(result_r, 8)

## problem 6.9
signal = [0, 1, 0, 1, 1, 0, 1, 0]
h1 = CPS.fft_radix2_dit_r(signal)
h1 ≈ fft(signal)

@benchmark CPS.fft_radix2_dit_r(x) setup = (x = rand(2^20) .|> ComplexF64) seconds = 10
@benchmark CPS.fft_algorithm(x) setup = (x = rand(2^20) .|> ComplexF64) seconds = 10
@benchmark CPS.fft_algorithm_optimized(x) setup = (x = rand(2^20) .|> ComplexF64) seconds = 10