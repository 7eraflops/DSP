begin
    using Random
    using Profile
    using BenchmarkTools
end
x = rand(Complex{Float64}, 2^30)
sizeof(x) / 1024^3

function optimized_fft(x::AbstractVector)::Vector
    N = length(x)
    x = Complex{Float64}.(x)
    W = cispi(-2 / N)  # Precompute W
    Wi = 1.0

    for e = 1:floor(Int, log2(N))
        L = 1 << e  # Use bitwise shift for powers of 2
        M = L >> 1  # Equivalent to 2^(e - 1)

        for m = 1:M
            for g = m:L:N
                d = g + M
                T = x[d] * Wi
                x[d] = x[g] - T
                x[g] = x[g] + T
            end
            Wi *= W
        end
    end

    return x
end



fft_algorithm(x)

Profile.clear()
for i in 1:10
    @profile fft_algorithm(x)
end
@profview fft_algorithm(x)

@benchmark optimized_fft(x) setup = (x = rand(2^20) .|> ComplexF64) seconds = 20
