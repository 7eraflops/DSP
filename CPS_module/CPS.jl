# module CPS

using LinearAlgebra

author = Dict{Symbol,String}(
    :index => "417096",
    :name => "Michał Charaszkiewicz",
    :email => "mchar@student.agh.edu.pl",
    :group => "0",
)

# Sygnały ciągłe
cw_rectangular(t::Real; T=1.0)::Real = abs(t) < T / 2 ? 1 : (abs(t) == T / 2 ? 0.5 : 0)
cw_triangle(t::Real; T=1.0)::Real = abs(t) < T ? 1 - abs(t) : 0
cw_literka_M(t::Real; T=1.0)::Real = missing
cw_literka_U(t::Real; T=1.0)::Real = missing

ramp_wave(t::Real)::Real = 2 * rem(t, 1, RoundNearest)
sawtooth_wave(t::Real)::Real = -2 * rem(t, 1, RoundNearest)
triangular_wave(t::Real)::Real = ifelse(mod(t + 1 / 4, 1.0) < 1 / 2, 4mod(t + 1 / 4, 1.0) - 1, -4mod(t + 1 / 4, 1.0) + 3)
square_wave(t::Real)::Real = ifelse(rem(t, 1) < 0.5, 1, -1)
pulse_wave(t::Real, ρ::Real)::Real = ifelse(rem(t, 1) < ρ, 1, 0)
impulse_repeater(g::Function, t1::Real, t2::Real)::Function = missing

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    temp = 0
    n = 1
    while (arg = 2π * n * (1 / T)) < band
        temp += (-1)^n * sin.(arg * t) / n
        n += 1
    end
    signal += -2A / π * temp
    return signal
end

function sawtooth_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    n = 1
    while (arg = 2π * n * (1 / T)) < band
        signal += (-1)^n * sin.(arg * t) / n
        n += 1
    end
    signal *= 2A / π
    return signal
end

function triangular_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    n = 1
    while (arg = 2π * n * (1 / T)) < band
        signal += (-1)^((n - 1) / 2) * sin.(arg * t) / n^2
        n += 2
    end
    signal *= (8A / π^2)
    return signal
end

function square_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    n = 1
    while (arg = 2π * (2n - 1) * (1 / T)) < band
        signal += sin.(arg * t) / (2n - 1)
        n += 1
    end
    signal *= 4 * A / π
    return signal
end

function pulse_wave_bl(t; ρ=0.2, A=1.0, T=1.0, band=20.0)
    missing
end

function impulse_repeater_bl(g::Function, t0::Real, t1::Real, band::Real)::Function
    missing
end

function rand_signal_bl(f1::Real, f2::Real)::Function
    missing
end

# Sygnały dyskretne
kronecker(n::Integer)::Real = ifelse(n == 0, 1, 0)
heaviside(n::Integer)::Real = ifelse(n < 0, 0, 1)

# Okna
rect(N::Integer)::AbstractVector{<:Real} = ones(N)
triang(N::Integer)::AbstractVector{<:Real} = [1 - (2abs(n - ((N - 1) / 2))) / (N - 1) for n = 0:N-1]
hanning(N::Integer)::AbstractVector{<:Real} = [0.5(1 - cos(2π * n / (N - 1))) for n = 0:N-1]
hamming(N::Integer)::AbstractVector{<:Real} = [0.54 - 0.46cos(2π * n / (N - 1)) for n = 0:N-1]
blackman(N::Integer)::AbstractVector{<:Real} = [0.42 - 0.5cos(2π * n / (N - 1)) + 0.08cos(4π * n / (N - 1)) for n = 0:N-1]

# Parametry sygnałów
mean(x::AbstractVector)::Number = sum(x) / length(x)
peak2peak(x::AbstractVector)::Real = abs(max(x) - min(x))
energy(x::AbstractVector)::Real = sum(abs2, x)
power(x::AbstractVector)::Real = sum(abs2, x) / length(x)
rms(x::AbstractVector)::Real = sqrt(sum(abs2, x) / length(x))

function running_mean(x::AbstractVector, M::Integer)::Vector
    result::AbstractVector = zeros(length(x))
    for k in 1:length(x)
        n₁ = k - M < 1 ? 1 : k - M
        n₂ = k + M > lastindex(x) ? lastindex(x) : k + M
        result[k] = (1 / (n₂ - n₁ + 1)) * sum(x[n₁:n₂])
    end
    return result
end

function running_energy(x::AbstractVector, M::Integer)::Vector
    result::AbstractVector = zeros(length(x))
    for k in 1:length(x)
        n₁ = k - M < 1 ? 1 : k - M
        n₂ = k + M > lastindex(x) ? lastindex(x) : k + M
        result[k] = sum(abs2, x[n₁:n₂])
    end
    return result
end

function running_power(x::AbstractVector, M::Integer)::Vector
    result::AbstractVector = zeros(length(x))
    for k in 1:length(x)
        n₁ = k - M < 1 ? 1 : k - M
        n₂ = k + M > lastindex(x) ? lastindex(x) : k + M
        result[k] = (1 / (n₂ - n₁ + 1)) * sum(abs2, x[n₁:n₂])
    end
    return result
end

# Próbkowanie
function interpolate(
    m::AbstractVector,
    s::AbstractVector,
    kernel::Function=sinc
)::Real
    missing
end

# Kwantyzacja
quantize(L::AbstractVector)::Function = x -> L[argmin(abs.(-L .+ x))]
SQNR(N::Integer)::Real = 6.02N
SNR(Psignal, Pnoise)::Real = 10 * log10(Psignal / Pnoise)

# Obliczanie DFT
function dtft(f::Real; signal::AbstractVector, fs::Real)
    missing
end

function dft(x::AbstractVector)::Vector
    missing
end

function idft(X::AbstractVector)::Vector
    missing
end

function rdft(x::AbstractVector)::Vector
    missing
end

function irdft(X::AbstractVector, N::Integer)::Vector
    missing
end

function fft_radix2_dit_r(x::AbstractVector)::Vector
    missing
end

function ifft_radix2_dif_r(X::AbstractVector)::Vector
    missing
end

function fft(x::AbstractVector)::Vector
    dft(x) # Może da rade lepiej?
end

function ifft(X::AbstractVector)::Vector
    idft(X) # Może da rade lepiej?
end

fftfreq(N::Integer, fs::Real)::Vector = missing
rfftfreq(N::Integer, fs::Real)::Vector = missing
amplitude_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = missing
power_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = missing
psd(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector = missing

function periodogram(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector
    missing
end

function stft(x::AbstractVector, w::AbstractVector, L::Integer)::Matrix
    missing
end

function istft(X::AbstractMatrix{<:Complex}, w::AbstractVector{<:Real}, L::Integer)::AbstractVector{<:Real}
    missing
end