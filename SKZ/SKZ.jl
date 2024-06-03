begin
    using Plots
    using LinearAlgebra
end

#= Zadanie 1: 
#* correct result
Oblicz wartość skuteczną dyskretnego sygnału x ∈ R^54. Dyskretny sygnał x powstał w wyniku 
pobrania N = 54 próbek z ciągłego sygnału y(t) = 1.8 * g(4.8 * t - 0.5) z szybkością fp = 385.71 
próbke na sekundę. Pierwsza próbka x1 = y(t1) została pobrana w chwili t1 = 5.39. Funkcja g(t) 
zwraca wartości sygnału fali piłokształtnej o opadającym zboczu i następujących parametrach: 
amplituda 1, okres 1 sekunda, składowa stała 0, g(0) = 0, oraz dg/dt |t=0 = -1.
=#
begin
    function rozwiazanie_1(;
        fp::Float64=385.71,
        t1::Float64=5.39,
        N::Int=54,
    )
        g = t -> -2 * rem(t, 1, RoundNearest)
        t = range(; start=t1, step=(1 / fp), length=N)
        y = 1.8 * g.(4.8 .* t .- 0.5)
        rms = sqrt(sum(abs2, y) / length(y))
        return rms
        1.1370311872575378
    end
    rms = rozwiazanie_1()
end
#= Zadanie 2:
#* correct solution
Z ciągłego sygnału f(t) ∈ R o paśmie ograniczonym od dołu i od góry przez częstotliwość |B| < 1/2Δm,
zostało pobrane 69 próbek w równych odstępach czasu Δm = m_{n+1} - m_n. wartości sygnału oraz momenty
w których zostały pobrane kolejne próbki, znajdują się odpowiednio w wektorze s ∈ R^69 oraz w wektorze
m ∈ R^69, gdzie s_n = f(m_n). Na podstawie wektorów m oraz s, znajdź sygnał g(t), będący rekonstrukcją
sygnału f(t) otrzymaną z wykorzystaniem wzoru interpolacyjnego Whittakera-Shannona. Jako rozwiazanie
podaj sumę wartości sygnału g(t) dla momentów t ∈ R^4, to znaczy ∑ n=1 -> 4 g(t_n)
=#
begin
    function rozwiazanie_2(;
        m::Vector{Float64}=[2.0, 2.0089, 2.0178, 2.0267, 2.0356, 2.0445, 2.0534, 2.0623, 2.0712, 2.0801, 2.089, 2.0979, 2.1068, 2.1157, 2.1246, 2.1335, 2.1424, 2.1513, 2.1602, 2.1691, 2.178, 2.1869, 2.1958, 2.2047, 2.2136, 2.2225, 2.2314, 2.2403, 2.2492, 2.2581, 2.267, 2.2759, 2.2848, 2.2937, 2.3026, 2.3115, 2.3204, 2.3293, 2.3382, 2.3471, 2.356, 2.3649, 2.3738, 2.3827, 2.3916, 2.4005, 2.4094, 2.4183, 2.4272, 2.4361, 2.445, 2.4539, 2.4628, 2.4717, 2.4806, 2.4895, 2.4984, 2.5073, 2.5162, 2.5251, 2.534, 2.5429, 2.5518, 2.5607, 2.5696, 2.5785, 2.5874, 2.5963, 2.6052],
        s::Vector{Float64}=[0.2096, 0.078, 0.5424, 0.1371, 0.8907, 0.8571, 0.9054, 0.8251, 0.8257, 0.0806, 0.0627, 0.4738, 0.3664, 0.3949, 0.2519, 0.8641, 0.0585, 0.2472, 0.7465, 0.6552, 0.4091, 0.0134, 0.3141, 0.1118, 0.1763, 0.9125, 0.2146, 0.0039, 0.0613, 0.1577, 0.7172, 0.8305, 0.3476, 0.1374, 0.9451, 0.6567, 0.17, 0.8519, 0.0879, 0.0532, 0.5134, 0.099, 0.9124, 0.6079, 0.8783, 0.8155, 0.4724, 0.2542, 0.1515, 0.4467, 0.6567, 0.8456, 0.8264, 0.6765, 0.1342, 0.5729, 0.193, 0.5653, 0.3538, 0.5603, 0.0119, 0.7747, 0.8549, 0.4851, 0.3145, 0.4572, 0.952, 0.0376, 0.4361],
        t::Vector{Float64}=[2.00445, 2.50107, 2.16732, 2.46547],
    )
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
        g = interpolate(m, s)
        solution = sum(g.(t))
        return solution
        1.861809777968849
    end
    solution = rozwiazanie_2()
end
#= Zadanie 3:
#* correct solution
Dany jest idealny równomierny 5-bitowy kwantyzator q(x), którego najmniejszy poziom kwantyzacji ma
wartość a = -0.88, natomiast największy poziom kwantyzacji ma wartość b = 1.1. Oblicz sygnał błęd
kwantyzacji tego przetwornika dla dyskretnego sygnału x ∈ R^78. Jako rozwiazanie podaj energię sygnału
błędu kwantyzacji
=#
begin
    function rozwiazanie_3(;
        a::Float64=-0.88,
        b::Float64=1.1,
        x::Vector{Float64}=[0.7, 0.67743, 0.65485, 0.63228, 0.60971, 0.58713, 0.56456, 0.54199, 0.51941, 0.49684, 0.47427, 0.45169, 0.42912, 0.40655, 0.38397, 0.3614, 0.33883, 0.31625, 0.29368, 0.27111, 0.24853, 0.22596, 0.20339, 0.18081, 0.15824, 0.13567, 0.11309, 0.09052, 0.06795, 0.04537, 0.0228, 0.00023, -0.02235, -0.04492, -0.06749, -0.09007, -0.11264, -0.13521, -0.15779, -0.18036, -0.20293, -0.22551, -0.24808, -0.27065, -0.29323, -0.3158, -0.33837, -0.36095, -0.38352, -0.40609, -0.42867, -0.45124, -0.47381, -0.49639, -0.51896, -0.54153, -0.56411, -0.58668, -0.60926, -0.63183, -0.6544, -0.67698, -0.69955, -0.72212, -0.7447, -0.76727, -0.78984, -0.81242, -0.83499, -0.85756, -0.88014, 1.09729, 1.07472, 1.05214, 1.02957, 1.007, 0.98442, 0.96185],
    )
        L = range(; start=a, stop=b, length=2^5)
        quantize(L::AbstractVector)::Function = x -> L[argmin(abs.(-L .+ x))]
        q = quantize(L)
        x_quantized = q.(x)
        quantization_error = x .- x_quantized
        energy = sum(abs2, quantization_error)
        return energy
        0.025959334412591073
    end
    energy = rozwiazanie_3()
end
#= Zadanie 4:
#* correct solution
DFT, x ∈ C^26, f_p = 546, suma faz dla f
=#
begin
    function rozwiazanie_4(;
        fp::Int=546,
        x::Vector{ComplexF64}=ComplexF64[0.72-0.06im, -1.18-0.02im, -0.81+0.53im, 0.51+0.6im, -0.33-1.04im, -0.31-0.44im, 0.24+0.15im, 0.74+0.34im, -1.41-0.11im, 0.14+0.95im, 0.23-0.83im, -1.32+0.49im, 0.16+0.54im, -0.6-0.28im, -0.3-0.52im, 0.18-0.4im, -0.69+0.58im, -0.95-0.84im, -0.13+0.37im, 0.58+2.19im, 0.61+0.27im, -0.9-0.2im, -1.48-0.01im, 0.55+0.2im, 0.28+0.0im, 0.32-0.76im],
        f::Vector{Int}=[105, -189, 273, -105],
    )
        function dft(x::AbstractVector)::Vector
            N = length(x)
            ζ = [cispi(-2 * n / N) for n in 0:(N-1)]
            [
                sum((
                    x[n+1] * ζ[(n*f)%N+1]
                    for n in 0:(N-1)
                ))
                for f in 0:(N-1)
            ]
        end
        function freq_to_index(f, N, fp)
            k = f * N / fp
            return Int(round(mod(k, N))) + 1
        end
        X = dft(x)
        phases = [angle(X[freq_to_index(freq, length(x), fp)]) for freq in f]
        phase_sum = sum(phases)
        return phase_sum
        -1.1645946595139476
    end
    phase_sum = rozwiazanie_4()
end
#= Zadanie 5:
#* correct solution
Dyskretny sygnał x ∈ R^69 został przetworzony przez dyskretny nierekurencyjny układ liniowy 
o odpowiedź impulsowej h ∈ R^18. Znajdź dyskretny sygnał y[n] będący sygnałem wyjściowym 
z tego układu. Jako rozwiązanie podaj energię otrzymanego sygnału.
=#
begin
    function rozwiazanie_5(;
        x::Vector{Float64}=[2.7, -3.51, 3.34, -3.94, 3.77, 4.02, 0.0, 1.62, -2.47, -2.06, 1.13, 0.87, -2.5, -0.96, -1.35, -1.47, -0.5, -3.13, 1.89, -2.1, 1.73, 3.28, 2.66, -0.77, 2.01, 3.66, -1.21, -4.0, 1.6, 3.35, 3.43, 1.72, -1.23, 4.95, -0.57, -3.3, -0.32, -2.68, -3.26, -4.23, -4.65, 2.56, -4.8, 1.99, 1.36, 4.06, -2.78, -1.04, -2.44, 4.84, 2.5, -1.48, 2.86, 1.6, 0.71, -1.9, -3.22, 0.15, -4.23, -0.3, 2.21, 0.32, -4.88, -1.67, 3.85, -4.37, 0.05, 3.44, -2.71],
        h::Vector{Float64}=[-4.81, 4.18, -1.56, -3.35, 3.63, 2.73, 2.5, -1.57, -4.17, -2.95, 2.05, 2.5, -4.25, -1.95, -2.86, -3.12, -4.87, -0.76],
    )
        n = length(x)
        m = length(h)
        y = zeros(eltype(x), n + m - 1)
        for i in 1:n
            for j in 1:m
                y[i+j-1] += x[i] * h[j]
            end
        end
        energy = sum(abs2, y)
        return energy
        87875.41206419999
    end
    energy = rozwiazanie_5()

end

function rozwiazanie_6(;
    b::Vector{Float64}=[0.5910613303216834, -6.112781006955527, 29.875443473694183, -91.01885078406374, 192.2281779559463, -296.2040756929196, 341.2822493432759, -296.2040756929196, 192.22817795594628, -91.01885078406374, 29.875443473694183, -6.112781006955528, 0.5910613303216835],
    a::Vector{Float64}=[1.0, -9.606819706803597, 43.657203099722054, -123.83997365628642, 243.91699531133975, -351.1790246082447, 378.8724249931937, -308.64912149362937, 188.52967696239602, -84.28273922211159, 26.212348155433993, -5.102166378720933, 0.4714208299517704],
    x::Vector{Float64}=[-0.09, -0.05, -0.19, -0.9, -0.48, 0.78, -0.84, -0.79, -0.17, -0.64, 0.88, -0.83, 0.79, -0.84, -0.82, 0.54, -0.51, 0.0, -0.67, -0.21, -0.36, -0.94, 0.37, -0.86, -0.53, -0.9, -0.62, 0.22, 1.0, 0.62, 0.14, 0.22, -0.04, 0.26, 0.71, 0.94, 0.36, 0.95, 0.8, -0.24, 0.46, -0.03, -0.82, -0.2, -0.96],
    L::Int=85,
)
    -0.04716978912749659
    missing
end

#= Zadanie 8
#* correct solution
dyskretny system liniowy, stacjonarny, wyznacz stabilność (1 - tak/0 - na granicy/-1 - nie)
a - mianownik, b - licznik
=#
begin
    function rozwiazanie_8(;
        b::Vector{Float64}=[0.25420348460046177, -1.2710174230023088, 2.5420348460046176, -2.5420348460046176, 1.2710174230023088, -0.25420348460046177],
        a::Vector{Float64}=[1.0, -4.62624610600513, 8.01514803904025, -6.118940600685512, 2.2495829771575107, -0.32498311549368936],
    )
        a = a / a[1]
        H = Matrix(I, length(a) - 2, length(a) - 2)
        Z = transpose(zeros(length(a) - 2))
        H = vcat(Z, H)
        H = hcat(H, -1 * reverse(a[2:end]))
        roots = eigvals(H)
        lengths = abs.(roots)
        for i in eachindex(lengths)
            if lengths[i] > 1
                return -1
            end
        end
        for i in eachindex(lengths)
            if lengths[i] == 1
                return 0
            end
        end
        return 1
        -1
    end
    stab = rozwiazanie_8()
end