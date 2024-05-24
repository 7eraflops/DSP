begin
    using Plots
end

#= Zadanie 1: 
* correct result
Oblicz wartość skuteczną dyskretnego sygnału x ∈ R^54. Dyskretny sygnał x powstał w wyniku 
pobrania N = 54 próbek z ciągłego sygnału y(t) = 1.8 * g(4.8 * t - 0.5) z szybkością fp = 385.71 
próbke na sekundę. Pierwsza próbka x1 = y(t1) została pobrana w chwili t1 = 5.39. Funkcja g(t) 
zwraca wartości sygnału fali piłokształtnej o opadającym zboczu i następujących parametrach: 
amplituda 1, okres 1 sekunda, składowa stała 0, g(0) = 0, oraz dg/dt |t=0 = -1.
=#
function rozwiazanie(;
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
rms = rozwiazanie()

#= Zadanie 2:
* correct solution
Oblicz energię dyskretnego sygnału x ∈ R^58. Dyskretny sygnał x powstał w wyniku 
pobrania N = 58 próbek z ciągłego sygnału y(t) = 1.5 * g(2.8 * t - 4.1) z szybkością 
fp = 360.44 próbke na sekundę. Pierwsza próbka x1 = y(t1) została pobrana w chwili 
t1 = -2.7. Funkcja g(t) zwraca wartości sygnału bipolarnej fali prostokątnej o 
następujących parametrach: amplituda 1, okres 1 sekunda, składowa stała 0, g(0) = 0, 
oraz g(t) = 1 dla t ∈ (0, 1/2).
=#
function rozwiazanie(;
    fp::Float64=360.44,
    t1::Float64=-2.7,
    N::Int=58,
)
    g = t -> ifelse(mod(t, 1) < 0.5, 1, -1)
    t = t = range(; start=t1, step=(1 / fp), length=N)
    y = 1.5 * g.(2.8 .* t .- 4.1)
    energy = sum(abs2, y)
    return energy
    130.5
end
energy = rozwiazanie()

#= Zadanie 3:
* correct solution
Dyskretny sygnał x ∈ R^69 został przetworzony przez dyskretny nierekurencyjny układ liniowy 
o odpowiedź impulsowej h ∈ R^18. Znajdź dyskretny sygnał y[n] będący sygnałem wyjściowym 
z tego układu. Jako rozwiązanie podaj energię otrzymanego sygnału.
=#
function rozwiazanie(;
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
energy = rozwiazanie()

#= Zadanie 4
* correct solution
Oblicz moc dyskretnego sygnału x ∈ R^259. Dyskretny sygnał x powstał w wyniku 
pobrania N = 259 próbek z ciągłego sygnału y(t) = 3.7 * g(3.1 * t - 5.0) z szybkością 
fp = 124.42 próbke na sekundę. Pierwsza próbka x1 = y(t1) została pobrana w chwili 
t1 = -4.16. Funkcja g(t) zwraca wartości sygnału fali piłokształtnej o opadającym 
zboczu i następujących parametrach: amplituda 1, okres 1 sekunda, składowa stała 0,
g(0) = 0, oraz dg/dt |t=0 = -1.
=#
function rozwiazanie(;
    fp::Float64=124.42,
    t1::Float64=-4.16,
    N::Int=259,
)
    g = t -> -2 * rem(t, 1, RoundNearest)
    t = range(; start=t1, step=(1 / fp), length=N)
    y = 3.7 * g.(3.1 .* t .- 5.0)
    power = sum(abs2, y) / length(y)
    return power
    4.681103148861615
end
power = rozwiazanie()

#= Zadanie 5
* correct solution
Oblicz moc dyskretnego sygnału x ∈ R^377. Dyskretny sygnał x powstał w wyniku 
pobrania N = 377 próbek z ciągłego sygnału y(t) = 3.0 * g(0.2 * t - 4.2) z szybkością 
fp = 325.1 próbke na sekundę. Pierwsza próbka x1 = y(t1) została pobrana w chwili 
t1 = 4.39. Funkcja g(t) zwraca wartości sygnału fali piłokształtnej o opadającym 
zboczu i następujących parametrach: amplituda 1, okres 1 sekunda, składowa stała 0,
g(0) = 0, oraz dg/dt |t=0 = -1.
=#
function rozwiazanie(;
    fp::Float64=325.1,
    t1::Float64=4.39,
    N::Int=377,
)
    g = t -> -2 * rem(t, 1, RoundNearest)
    t = range(; start=t1, step=(1 / fp), length=N)
    y = 3.0 * g.(0.2 .* t .- 4.2)
    rms = sqrt(sum(abs2, y) / length(y))
    return rms
    1.3016002840023688
end
rms = rozwiazanie()

#= Zadanie 6:
* correct solution
Dyskretny sygnał x ∈ R^54 został przetworzony przez dyskretny nierekurencyjny układ liniowy 
o odpowiedź impulsowej h ∈ R^22. Znajdź dyskretny sygnał y[n] będący sygnałem wyjściowym 
z tego układu. Jako rozwiązanie podaj moc otrzymanego sygnału.
=#
function rozwiazanie(;
    x::Vector{Float64}=[1.47, 2.35, 2.26, 4.42, -4.96, 3.11, 1.56, 3.21, 1.62, -1.66, 0.98, -1.16, 2.09, 1.92, 1.12, -1.75, 4.83, -2.07, 2.14, -4.91, 4.01, -0.79, 0.11, -0.22, 0.86, 4.81, 3.3, 4.83, 4.46, -3.47, -4.75, 0.07, -0.93, 4.98, -0.32, -3.15, 3.13, -1.81, -2.42, 0.52, -1.77, -2.59, 0.06, 0.19, 1.87, -3.85, -0.1, 3.92, -3.99, 1.63, 4.18, -2.82, -0.72, 4.45],
    h::Vector{Float64}=[-0.25, 2.09, 0.89, -4.39, 3.75, 0.57, 4.49, -4.74, 3.63, -4.37, -3.6, -4.32, 4.27, 1.45, 0.06, 3.52, 2.98, -1.54, 3.74, 0.43, 0.47, -2.71],
)
    n = length(x)
    m = length(h)
    y = zeros(eltype(x), n + m - 1)
    for i in 1:n
        for j in 1:m
            y[i+j-1] += x[i] * h[j]
        end
    end
    power = sum(abs2, y) / length(y)
    return power
    1219.5582110730666
end
power = rozwiazanie()
∑
#= Zadanie 7:
* correct solution
Z ciągłego sygnału f(t) ∈ R o paśmie ograniczonym od dołu i od góry przez częstotliwość |B| < 1/2Δm,
zostało pobrane 69 próbek w równych odstępach czasu Δm = m_{n+1} - m_n. wartości sygnału oraz momenty
w których zostały pobrane kolejne próbki, znajdują się odpowiednio w wektorze s ∈ R^69 oraz w wektorze
m ∈ R^69, gdzie s_n = f(m_n). Na podstawie wektorów m oraz s, znajdź sygnał g(t), będący rekonstrukcją
sygnału f(t) otrzymaną z wykorzystaniem wzoru interpolacyjnego Whittakera-Shannona. Jako rozwiazanie
podaj sumę wartości sygnału g(t) dla momentów t ∈ R^4, to znaczy ∑ n=1 -> 4 g(t_n)
=#
function rozwiazanie(;
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
solution = rozwiazanie()

#= Zadanie 8:
* correct solution
Dyskretny sygnał x ∈ R^70 został przetworzony przez dyskretny nierekurencyjny układ liniowy 
o odpowiedź impulsowej h ∈ R^18. Znajdź dyskretny sygnał y[n] będący sygnałem wyjściowym 
z tego układu. Jako rozwiązanie podaj energię otrzymanego sygnału.
=#
function rozwiazanie(;
    x::Vector{Float64}=[4.19, 3.72, -1.2, -3.52, 2.13, -4.31, -2.54, -0.11, -0.7, 1.35, -3.43, 0.42, 3.77, -2.07, -2.58, 4.89, 1.55, 2.45, 4.76, -0.25, 4.78, 0.41, 4.1, 0.04, -2.11, 3.61, 0.69, 2.72, 4.48, -4.07, 1.09, -3.33, -2.06, -1.8, -2.19, -3.92, -0.66, 3.21, 0.16, 4.1, 2.15, 1.76, 2.58, 4.89, -1.51, -3.53, -2.75, 0.35, 1.54, -3.81, -2.02, 2.03, 4.4, 2.54, -4.59, -0.72, 3.7, -4.43, -0.13, -3.32, -4.69, 3.92, -2.94, -1.67, -2.92, 4.89, 2.33, -2.87, -1.45, -4.99],
    h::Vector{Float64}=[2.19, -3.19, 4.47, -3.65, -2.36, -0.93, 4.74, 0.87, 3.42, 2.0, 1.47, -2.75, -4.74, 0.37, 0.98, 0.95, 4.42, -0.91],
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
    68175.53795130999
end
energy = rozwiazanie()

#= Zadanie 9:
* correct solution
Dany jest idealny równomierny 5-bitowy kwantyzator q(x), którego najmniejszy poziom kwantyzacji ma
wartość a = -0.88, natomiast największy poziom kwantyzacji ma wartość b = 1.1. Oblicz sygnał błęd
kwantyzacji tego przetwornika dla dyskretnego sygnału x ∈ R^78. Jako rozwiazanie podaj energię sygnału
błędu kwantyzacji
=#
function rozwiazanie(;
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
energy = rozwiazanie()