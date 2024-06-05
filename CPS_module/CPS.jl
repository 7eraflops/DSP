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

# TODO: QuadGK -> own quadrature integral implementation
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
function interpolate(m::AbstractVector, s::AbstractVector, kernel::Function=sinc)::Function
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

function irdft(X::AbstractVector, N::Integer)::Vector
    S = length(X)
    X₁ = [n <= S ? X[n] : conj(X[2S-n+(N % 2 == 0 ? 0 : 1)]) for n in 1:N]
    real.(idft(X₁))
end

function fft_radix2_dit_r!(x::Vector{Complex{T}}) where {T}
    N = length(x)
    if N == 0 || (N & (N - 1)) != 0
        throw(ArgumentError("Length of input must be a power of 2"))
    end

    bits = Int(log2(N))
    function bitreverse(n, bits)
        reversed = 0
        for i in 1:bits
            reversed <<= 1
            reversed |= (n & 1)
            n >>= 1
        end
        return reversed
    end

    for i in 1:N
        j = bitreverse(i - 1, bits) + 1
        if i < j
            x[i], x[j] = x[j], x[i]
        end
    end

    m = 2
    while m <= N
        half_m = div(m, 2)
        w_m = cispi(-2 / m)
        for k in 1:m:N
            w = one(Complex{T})
            for j in 0:half_m-1
                t = w * x[k+j+half_m]
                u = x[k+j]
                x[k+j] = u + t
                x[k+j+half_m] = u - t
                w *= w_m
            end
        end
        m *= 2
    end

    return x
end

function fft_radix2_dit_r(x::Vector{Complex{T}}) where {T}
    x_copy = copy(x)
    return fft_radix2_dit_r!(x_copy)
end

function fft_radix2_dit_r(x::Vector{T}) where {T}
    x_complex = Complex{promote_type(T, Float64)}[complex(xi, 0.0) for xi in x]
    return fft_radix2_dit_r!(x_complex)
end

function ifft_radix2_dif_r!(X::Vector{Complex{T}}) where {T}
    N = length(X)
    if N == 0 || (N & (N - 1)) != 0
        throw(ArgumentError("Length of input must be a power of 2"))
    end

    m = N
    while m >= 2
        half_m = div(m, 2)
        w_m = cispi(2 / m)
        for k in 1:m:N
            w = one(Complex{T})
            for j in 0:half_m-1
                t = X[k+j]
                u = X[k+j+half_m]
                X[k+j] = t + u
                X[k+j+half_m] = (t - u) * w
                w *= w_m
            end
        end
        m = half_m
    end

    bits = Int(log2(N))
    function bitreverse(n, bits)
        reversed = 0
        for i in 1:bits
            reversed <<= 1
            reversed |= (n & 1)
            n >>= 1
        end
        return reversed
    end

    for i in 1:N
        j = bitreverse(i - 1, bits) + 1
        if i < j
            X[i], X[j] = X[j], X[i]
        end
    end

    scale = inv(T(N))
    for i in 1:N
        X[i] *= scale
    end

    return X
end

function ifft_radix2_dif_r(X::Vector{Complex{T}}) where {T}
    X_copy = copy(X)
    return ifft_radix2_dif_r!(X_copy)
end

function fft(x::AbstractVector)::Vector
    dft(x) # Może da rade lepiej?
end

function ifft(X::AbstractVector)::Vector
    idft(X) # Może da rade lepiej?
end

fftfreq(N::Integer, fs::Real)::Vector = [n * N / fs for n in 0:(N-1)]
rfftfreq(N::Integer, fs::Real)::Vector = [n * N / fs for n in 0:(N÷2)]
amplitude_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = abs.(fft(x .* w)) / (length(x) * mean(w))
power_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = amplitude_spectrum(x, w) .^ 2
psd(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector = abs2.(fft(x .* w)) / (sum(abs2, w) * fs)

function periodogram(x::AbstractVector, w::AbstractVector=rect(length(x)), L::Integer=0, fs::Real=1.0)::Vector
    missing
end

function stft(x::AbstractVector, w::AbstractVector, L::Integer)::Matrix
    missing
end

function istft(X::AbstractMatrix{<:Complex}, w::AbstractVector{<:Real}, L::Integer)::AbstractVector{<:Real}
    missing
end

function conv(f::AbstractVector, g::AbstractVector)::Vector
    n = length(f)
    m = length(g)
    y = zeros(eltype(f), n + m - 1)
    for i in 1:n
        for j in 1:m
            y[i+j-1] += f[i] * g[j]
        end
    end
    return y
end

function fast_conv(f::Vector, g::Vector)::Vector
    N = length(f) + length(g) - 1

    f_padded = vcat(f, zeros(N - length(f)))
    g_padded = vcat(g, zeros(N - length(g)))

    F = fft(f_padded)
    G = fft(g_padded)
    Y = F .* G
    y = real(ifft(Y))

    return y
end

# FIXME: fix reading array out of bounds
function overlap_add(x::Vector, h::Vector, L::Integer)::Vector
    M = length(h)
    N = L + M - 1

    padded_h = vcat(h, zeros(N - M))
    H = fft(padded_h)

    y = zeros(eltype(x), length(x) + M - 1)

    for k in 1:L:length(x)
        xk = x[k:min(k + L - 1, end)]
        padded_xk = vcat(xk, zeros(N - length(xk)))
        Xk = fft(padded_xk)
        Yk = ifft(H .* Xk)
        y[k:k+N-1] += real(Yk)
    end

    return y[1:length(x)+M-1]
end

function overlap_save(x::Vector, h::Vector, L::Integer)::Vector
    M = length(h)
    N = L + M - 1

    padded_h = vcat(h, zeros(N - M))
    H = fft(padded_h)

    y = []
    padded_x = vcat(zeros(M - 1), x, zeros(N - 1))

    for k in 1:L:(length(padded_x)-N+1)
        xk = padded_x[k:k+N-1]
        Xk = fft(xk)
        Yk = ifft(H .* Xk)
        y = vcat(y, real(Yk[M:end]))
    end

    return y[1:(length(x)+M-1)]
end

function lti_filter(b::Vector, a::Vector, x::Vector)::Vector
    N = length(x)
    M = length(b) - 1
    K = length(a) - 1
    y = zeros(Float64, N)

    for n in 1:N
        for k in 0:M
            if n - k > 0
                y[n] += b[k+1] * x_padded[n-k]
            end
        end
        for k in 1:K
            if n - k > 0
                y[n] -= a[k+1] * y[n-k]
            end
        end
    end
    return y
end

function filtfilt(b::Vector, a::Vector, x::Vector)::Vector
    missing
end

function lti_amp(f::Real, b::Vector, a::Vector)::Real
    M = length(b)
    K = length(a)
    num = sum(b[m+1] * cispi(-2 * f * m) for m in 0:M-1)
    denom = sum(a[k+1] * cispi(-2 * f * m) for k in 0:K-1)
    H_f = num / denom
    return abs(H_f)
end

function lti_phase(f::Real, b::Vector, a::Vector)::Real
    M = length(b)
    K = length(a)
    num = sum(b[m+1] * cispi(-2 * f * m) for m in 0:M-1)
    denom = sum(a[k+1] * cispi(-2 * f * m) for k in 0:K-1)
    H_f = num / denom
    return angle(H_f)
end

function firwin_lp_I(order::Integer, F0::Float64)::Vector
    return [2F0 * sinc(2F0 * n) for n in -order/2:order/2]
end

function firwin_hp_I(order::Integer, F0::Float64)::Vector
    return [kronecker(Int(n)) - 2F0 * sinc(2F0 * n) for n in -order/2:order/2]
end

function firwin_bp_I(order::Integer, F1::Float64, F2::Float64)::Vector
    return [2F2 * sinc(2F2 * n) - 2F1 * sinc(2F1 * n) for n in -order/2:order/2]
end

function firwin_bs_I(order::Integer, F1::Float64, F2::Float64)::Vector
    return [kronecker(Int(n)) - (2F2 * sinc(2F2 * n) - 2F1 * sinc(2F1 * n)) for n in -order/2:order/2]
end

function firwin_lp_II(N::Integer, F0::Float64)::Vector
    missing
end

function firwin_bp_II(N::Integer, F1::Float64, F2::Float64)::Vector
    missing
end

function resample(x::Vector, M::Int, N::Int)
    missing
end

end