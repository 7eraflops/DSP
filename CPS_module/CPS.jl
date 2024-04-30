module CPS

using LinearAlgebra
using QuadGK
using OffsetArrays

author = Dict{Symbol,String}(
    :index => "417096",
    :name => "Michał Charaszkiewicz",
    :email => "mchar@student.agh.edu.pl",
    :group => "1",
)

# Sygnały ciągłe
cw_rectangular(t::Real; T=1.0)::Real = abs(t) < T / 2 ? 1 : (abs(t) == T / 2 ? 0.5 : 0)
cw_triangle(t::Real; T=1.0)::Real = abs(t) < T ? 1 - abs(t) : 0
cw_literka_M(t::Real; T=1.0)::Real = abs(t) < T ? (t < 0 ? -t + 1 : t + 1) : 0
cw_literka_U(t::Real; T=1.0)::Real = abs(t) < T ? t^2 : 0

ramp_wave(t::Real)::Real = 2 * rem(t, 1, RoundNearest)
sawtooth_wave(t::Real)::Real = -2 * rem(t, 1, RoundNearest)
triangular_wave(t::Real)::Real = ifelse(mod(t + 1 / 4, 1.0) < 1 / 2, 4mod(t + 1 / 4, 1.0) - 1, -4mod(t + 1 / 4, 1.0) + 3)
square_wave(t::Real)::Real = ifelse(mod(t, 1) < 0.5, 1, -1)
pulse_wave(t::Real, ρ::Real)::Real = ifelse(mod(t, 1) < ρ, 1, 0)
impulse_repeater(g::Function, t1::Real, t2::Real)::Function = x -> g(mod(x - t1, t2 - t1) + t1)

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    temp = 0
    n = 1
    while (arg = 2π * n * (1 / T)) < band * 2π
        temp += (-1)^n * sin.(arg * t) / n
        n += 1
    end
    signal += -2A / π * temp
    return signal
end

function sawtooth_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    n = 1
    while (arg = 2π * n * (1 / T)) < band * 2π
        signal += (-1)^n * sin.(arg * t) / n
        n += 1
    end
    signal *= 2A / π
    return signal
end

function triangular_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    n = 1
    while (arg = 2π * n * (1 / T)) < band * 2π
        signal += (-1)^((n - 1) / 2) * sin.(arg * t) / n^2
        n += 2
    end
    signal *= (8A / π^2)
    return signal
end

function square_wave_bl(t; A=1.0, T=1.0, band=20.0)
    signal = 0
    n = 1
    while (arg = 2π * (2n - 1) * (1 / T)) < band * 2π
        signal += sin.(arg * t) / (2n - 1)
        n += 1
    end
    signal *= 4 * A / π
    return signal
end

function pulse_wave_bl(t; ρ=0.2, A=1.0, T=1.0, band=20.0)
    signal = (sawtooth_wave_bl.(t .- (T / 2); A, T, band) - sawtooth_wave_bl.(t .- ((T / 2) + ρ); A, T, band)) .+ (2 * A * ρ)
    return signal
end

function impulse_repeater_bl(g::Function, t1::Real, t2::Real, band::Real)::Function
    T::Float64 = t2 - t1
    ω₀::Float64 = (2π / T)
    n_terms::Integer = div(band * 2π, ω₀)

    a0 = 1 / T * quadgk(g, t1, t2)[1]
    an_coeffs = zeros(Float64, n_terms)
    bn_coeffs = zeros(Float64, n_terms)

    for n in 1:n_terms
        an_coeffs[n] = 2 / T * quadgk(t -> g(t) * cos(ω₀ * n * t), t1, t2)[1]
        bn_coeffs[n] = 2 / T * quadgk(t -> g(t) * sin(ω₀ * n * t), t1, t2)[1]
    end

    function fourier_series_output(t::Float64)
        result = a0 / 2
        for n in 1:n_terms
            result += an_coeffs[n] * cos(ω₀ * n * t) + bn_coeffs[n] * sin(ω₀ * n * t)
        end
        return result
    end

    return fourier_series_output
end

# TODO: naive approach, might change later
function rand_signal_bl(f1::Real, f2::Real)::Function
    f = f1 .+ rand(1000) .* (f2 - f1)
    ϕ = -π .+ rand(1000) * 2π
    A = randn(1000)
    A = A ./ sqrt(0.5 * sum(A .^ 2))
    return t -> sum(A .* sin.(2π * f .* t .+ ϕ))
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
peak2peak(x::AbstractVector)::Real = abs(maximum(x) - minimum(x))
energy(x::AbstractVector)::Real = sum(abs2, x)
power(x::AbstractVector)::Real = energy(x) / length(x)
rms(x::AbstractVector)::Real = sqrt(power(x))

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
)::Function
    return x -> begin
        sum = 0.0
        Δt = m[2] - m[1]
        for i in eachindex(m)
            sum += s[i] * kernel((x - m[i]) / Δt)
        end
        return sum
    end
end

# Kwantyzacja
quantize(L::AbstractVector)::Function = x -> L[argmin(abs.(-L .+ x))]
SQNR(N::Integer)::Real = 1.76 + 6.02 * N # 6.02N [dB] also correct
SNR(Psignal, Pnoise)::Real = 10 * log10(Psignal / Pnoise)

# Obliczanie DFT
function dtft(freq::Real; signal::AbstractVector, fs::Real)
    dtft_val::ComplexF64 = 0.0
    for n in eachindex(signal)
        dtft_val += signal[n] * cispi(-2 * freq * n / fs)
    end
    return dtft_val
end

function dft(x::AbstractVector)::Vector
    N = length(x)
    ζ = OffsetArray(
        [cispi(-2 * n / N) for n in 0:(N-1)],
        0:(N-1)
    )
    [
        sum((
            x[n+1] * ζ[(n*f)%N]
            for n in 0:(N-1)
        ))
        for f in 0:(N-1)
    ]
end

function idft(X::AbstractVector)::Vector
    N = length(X)
    ζ = OffsetArray(
        [cispi(2 * n / N) for n in 0:(N-1)],
        0:(N-1)
    )
    [
        (1 / N) * sum((
            X[n+1] * ζ[(n*f)%N]
            for n in 0:(N-1)
        ))
        for f in 0:(N-1)
    ]
end

function rdft(x::AbstractVector)::Vector
    N = length(x)
    ζ = OffsetArray(
        [cispi(-2 * n / N) for n in 0:(N-1)],
        0:(N-1)
    )
    [
        sum((
            x[n+1] * ζ[(n*f)%N]
            for n in 0:(N-1)
        ))
        for f in 0:(N÷2)
    ]
end

# TODO: ask Mr. Woźniak about the N in the argument list
function irdft(X::AbstractVector, N::Integer)::Vector
    S = length(X)
    X₁ = [n <= S ? X[n] : conj(X[2S-n]) for n in 1:N]
    idft(X₁)
end

# TODO: still not complete
function fft_radix2_dit_r(x::AbstractVector)::Vector
    N = length(x)
    if !ispow2(N)
        error("Input vector length must be a power of 2")
    end
    x = Complex{Float64}.(x)
    for e = 1:floor(Int, log2(N))
        L = 2^e
        M = 2^(e - 1)
        Wi = 1
        W = cispi(-2 / L)
        for m = 1:M
            for g = m:L:N
                d = g + M
                @inbounds T = x[d] * Wi
                x[d] = x[g] - T
                x[g] = x[g] + T
            end
            Wi *= W
        end
    end
    return x
end

function fft_algorithm(x::AbstractVector)::Vector
    N = length(x)
    x = Complex{Float64}.(x)
    for e = 1:floor(Int, log2(N))
        L = 2^e
        M = 2^(e - 1)
        Wi = 1
        W = cispi(-2 / L)
        for m = 1:M
            for g = m:L:N
                d = g + M
                T = x[d] * Wi
                x[d] = x[g] - T
                x[g] = x[g] + T
            end
            Wi = Wi * W
        end
    end
    return x
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

function conv(f::Vector, g::Vector)::Vector
    missing
end

function fast_conv(f::Vector, g::Vector)::Vector
    missing
end

function overlap_add(x::Vector, h::Vector, L::Integer)::Vector
    missing
end

function overlap_save(x::Vector, h::Vector, L::Integer)::Vector
    missing
end

function lti_filter(b::Vector, a::Vector, x::Vector)::Vector
    missing
end
end