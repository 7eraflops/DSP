begin
    include("../CPS_module/CPS.jl")
    using Plots
    using FFTW
    using .CPS
    using BenchmarkTools
    using Profile
end

function normalize_frequency_to_pi(freq::AbstractVector{<:Real}, fs::Real)
    nyquist_freq = fs / 2
    normalized_freq = 2π * (freq .- nyquist_freq) / fs
    return normalized_freq
end

function shift_array(arr::Vector{T}) where {T}
    len = length(arr)
    mid_index = (len + 1) ÷ 2  # Middle index

    if len % 2 == 0
        # For even-length arrays, split into two equal halves
        first_half = arr[1:mid_index]
        second_half = arr[mid_index+1:end]
    else
        # For odd-length arrays, include the middle element in the first half
        first_half = arr[1:mid_index]
        second_half = arr[mid_index+1:end]
    end

    # Reconstruct the array with the desired shift
    shifted_array = [arr[mid_index]; second_half; first_half[1:end-1]]

    return shifted_array
end

x = n -> cospi(n / 3) + 2sinpi(n / 4)

signal = n -> (cospi(n / 4) + 0.001cospi(21 * n / 64))

n = 0:1:63

y = signal.(n)

fs = length(y)

# plot(y, title="Time-domain Signal", xlabel="Sample", ylabel="Amplitude")

# Calculate the amplitude spectrum
spectrum = shift_array(CPS.amplitude_spectrum(y))

# Calculate the frequency bins
freq = CPS.fftfreq(length(y), fs)

# Normalize frequencies to range from -pi to pi
normalized_freq = normalize_frequency_to_pi(freq, fs)

plot(
    normalized_freq,
    abs.(shift_array(CPS.amplitude_spectrum(x.(n)))),
    title="Frequency Spectrum",
    xlabel="Frequency (radians/sample)",
    ylabel="Amplitude",
    xticks=(-π:π/4:π, ["-π", "-3π/4", "-π/2", "-π/4", "0", "π/4", "π/2", "3π/4", "π"])
)

# Plot the frequency spectrum
plot(
    normalized_freq,
    abs.(spectrum),
    title="Frequency Spectrum",
    xlabel="Frequency (radians/sample)",
    ylabel="Amplitude",
    xticks=(-π:π/4:π, ["-π", "-3π/4", "-π/2", "-π/4", "0", "π/4", "π/2", "3π/4", "π"])
)

X = CPS.dft([20, 5])
x = CPS.idft(X)

p = CPS.fftfreq(2000, 1000)
p[1001]