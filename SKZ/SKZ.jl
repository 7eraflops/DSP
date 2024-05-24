begin
    using Plots
end

#= Zadanie 1: 
?# check result
Oblicz wartość skuteczną dyskretnego sygnału x ∈ R^184. Dyskretny sygnał x powstał w wyniku 
pobrania N = 184 próbek z ciągłego sygnału y(t) = 2.6 * g(2.9 * t - 4.8) z szybkością fp = 154.69 
próbke na sekundę. Pierwsza próbka x1 = y(t1) została pobrana w chwili t1 = -7.91. Funkcja g(t) 
zwraca wartości sygnału fali piłokształtnej o opadającym zboczu i następujących parametrach: 
amplituda 1, okres 1 sekunda, składowa stała 0, g(0) = 0, oraz dg/dt |t=0 = -1.
=#
function rozwiazanie(;
    fp::Float64=154.69,
    t1::Float64=-7.91,
    N::Int=184,
)
    g = t -> -2 * rem(t, 1, RoundNearest)
    t = range(; start=t1, step=(1 / fp), length=N)
    y = 2.6 * g.(2.9 .* t .- 4.8)
    rms = sqrt(sum(abs2, y) / length(y))
    return rms
    1.5816485480484588
end
rms = rozwiazanie()

#= Zadanie 2:
*# correct solution
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
*# correct solution
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