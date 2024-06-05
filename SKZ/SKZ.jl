begin
    using LinearAlgebra
end

# funkcje używne w zadaniach
begin
    # parametry sygnałów
    mean(x::AbstractVector)::Number = sum(x) / length(x)
    energy(x::AbstractVector)::Real = sum(abs2, x)
    power(x::AbstractVector)::Real = energy(x) / length(x)
    rms(x::AbstractVector)::Real = sqrt(power(x))
    # sygnały
    ramp_wave(t::Real)::Real = 2 * rem(t, 1, RoundNearest)
    sawtooth_wave(t::Real)::Real = -2 * rem(t, 1, RoundNearest)
    triangular_wave(t::Real)::Real = ifelse(mod(t + 1 / 4, 1.0) < 1 / 2, 4mod(t + 1 / 4, 1.0) - 1, -4mod(t + 1 / 4, 1.0) + 3)
    square_wave(t::Real)::Real = ifelse(mod(t, 1) < 0.5, 1, -1)
    # interploacja
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
    # dft i zamiana częstotliwości na index dft
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
    # splot
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
    # odpowiedź systemu LTI
    function lti_filter(b::Vector, a::Vector, x::Vector)::Vector
        N = length(x)
        M = length(b) - 1
        K = length(a) - 1
        y = zeros(Float64, N)

        for n in 1:N
            for k in 0:M
                if n - k > 0
                    y[n] += b[k+1] * x[n-k]
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
    # wzmocnienie/przesunięcie fazowye systemu LTI
    function lti_amp(f::Real, b::Vector, a::Vector)::Real
        M = length(b)
        K = length(a)
        num = sum(b[m+1] * cispi(-2 * f * m) for m in 0:M-1)
        denom = sum(a[k+1] * cispi(-2 * f * k) for k in 0:K-1)
        H_f = num / denom
        return abs(H_f)
    end
    function lti_phase(f::Real, b::Vector, a::Vector)::Real
        M = length(b)
        K = length(a)
        num = sum(b[m+1] * cispi(-2 * f * m) for m in 0:M-1)
        denom = sum(a[k+1] * cispi(-2 * f * k) for k in 0:K-1)
        H_f = num / denom
        return angle(H_f)
    end
    # obliczanie piwrwiastków wielomianu przy użyciu companion matrix (https://www.wikiwand.com/en/Companion_matrix)
    function roots(a::AbstractVector)::Vector
        H = Matrix(I, length(a) - 2, length(a) - 2)
        Z = zeros(length(a) - 2)
        H = vcat(Z', H)
        H = hcat(H, -1 * reverse(a[2:end]))
        return eigvals(H)
    end

    function poly_from_roots(r::AbstractVector)
        p = [1.0 + 0im]
        for i in eachindex(r)
            p = conv(p, [1, -r[i]])
        end
        return p
    end
    # filtry
    kronecker(n::Integer)::Real = ifelse(n == 0, 1, 0)
    firwin_lp_I(order::Integer, F0::Float64)::Vector = [2F0 * sinc(2F0 * n) for n in -order/2:order/2]
    firwin_hp_I(order::Integer, F0::Float64)::Vector = [kronecker(Int(n)) - 2F0 * sinc(2F0 * n) for n in -order/2:order/2]
    firwin_bp_I(order::Integer, F1::Float64, F2::Float64)::Vector = [2F2 * sinc(2F2 * n) - 2F1 * sinc(2F1 * n) for n in -order/2:order/2]
    firwin_bs_I(order::Integer, F1::Float64, F2::Float64)::Vector = [kronecker(Int(n)) - (2F2 * sinc(2F2 * n) - 2F1 * sinc(2F1 * n)) for n in -order/2:order/2]
    # okna do filtrów
    triang(M::Integer)::AbstractVector{<:Real} = [1 - abs(n) / (M + 1) for n = -M:M]
    hanning(M::Integer)::AbstractVector{<:Real} = [0.5(1 + cos(2π * n / (2M + 1))) for n = -M:M]
    hamming(M::Integer)::AbstractVector{<:Real} = [0.54 + 0.46cos(2π * n / (2M + 1)) for n = -M:M]
    blackman(M::Integer)::AbstractVector{<:Real} = [0.42 + 0.5cos(2π * n / (2M + 1)) + 0.08cos(4π * n / (2M + 1)) for n = -M:M]
end

#= Zadanie 1: 
#* correct solution
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
        g = ramp_wave # ramp/sawtooth/triangular/square
        t = range(; start=t1, step=(1 / fp), length=N)
        y = 1.8 * g.(4.8 .* t .- 0.5)
        return rms(y) # energy/power/rms
        1.1370311872575378
    end
    out_1 = rozwiazanie_1()
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
        g = interpolate(m, s)
        return sum(g.(t))
        1.861809777968849
    end
    out_2 = rozwiazanie_2()
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
        N = 5 # N-bitowy kwantyzator
        L = range(; start=a, stop=b, length=2^N)
        q = quantize(L)
        x_quantized = q.(x)
        quantization_error = x .- x_quantized
        return energy(quantization_error) # energy/power/rms
        0.025959334412591073
    end
    out_3 = rozwiazanie_3()
end

#= Zadanie 4:
#* correct solution
DFT, A lub ϕ
=#
# wariant A
begin
    function rozwiazanie_4_1(;
        fp::Int=1620,
        x::Vector{ComplexF64}=ComplexF64[0.36-0.96im, -0.12-0.24im, -0.71-0.15im, -0.41+0.04im, 0.41-0.01im, -0.42+0.46im, -1.07+0.16im, -0.2-0.5im, 0.91+1.18im, -0.04-0.11im, -0.95-0.29im, 0.06-0.89im, -0.56+0.46im, 0.68+0.04im, -0.67-1.77im, 0.33-0.24im, 1.08+0.44im, -1.14+0.16im, 0.17-0.18im, -0.15-0.3im, 0.85+0.3im, 1.04+1.48im, -0.37-0.03im, 0.3-0.79im, 0.21-0.37im, -0.29-1.15im, -0.44-0.6im, 0.67+0.53im, 1.17+0.93im, -0.09+0.82im, -0.7+0.27im, -1.03+0.99im, 1.1+0.58im, -0.36-0.49im, -1.18+0.14im, -0.03-0.33im, 1.1-0.51im, -0.23+0.36im, -0.49+0.01im, 0.32+0.77im, 0.22+0.58im, 0.25-0.72im, -0.3+0.04im, 0.85+1.03im, 0.43+0.92im],
        f::Vector{Int}=[-36, 252, 360, 396],
    )
        X = dft(x)
        A = [abs(X[freq_to_index(freq, length(x), fp)]) / length(X) for freq in f]
        return sum(A)
        0.5240689812893122
    end
    out_4 = rozwiazanie_4_1()
end

# wariant ϕ 
begin
    function rozwiazanie_4_2(;
        fp::Int=546,
        x::Vector{ComplexF64}=ComplexF64[0.72-0.06im, -1.18-0.02im, -0.81+0.53im, 0.51+0.6im, -0.33-1.04im, -0.31-0.44im, 0.24+0.15im, 0.74+0.34im, -1.41-0.11im, 0.14+0.95im, 0.23-0.83im, -1.32+0.49im, 0.16+0.54im, -0.6-0.28im, -0.3-0.52im, 0.18-0.4im, -0.69+0.58im, -0.95-0.84im, -0.13+0.37im, 0.58+2.19im, 0.61+0.27im, -0.9-0.2im, -1.48-0.01im, 0.55+0.2im, 0.28+0.0im, 0.32-0.76im],
        f::Vector{Int}=[105, -189, 273, -105],
    )
        X = dft(x)
        phases = [angle(X[freq_to_index(freq, length(x), fp)]) for freq in f] # angle -> abs dla wariantu z sumą amplitud
        return sum(phases)
        -1.1645946595139476
    end
    out_4 = rozwiazanie_4_2()
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
        y = conv(x, h)
        return energy(y) # energy/power/rms
        87875.41206419999
    end
    out_5 = rozwiazanie_5()
end

#= Zadanie 6:
#* correct solution
system LTI, odpowiedź na impuls x, policz energię/moc/rms dla pierwszych L próbek
=#
begin
    function rozwiazanie_6(;
        b::Vector{Float64}=[0.003170506666592921, 0.015852533332964606, 0.03170506666592921, 0.03170506666592921, 0.015852533332964606, 0.003170506666592921],
        a::Vector{Float64}=[1.0, -2.905844673274107, 4.319099303556442, -3.7648421381727513, 1.9147028507947308, -0.4616591295733407],
        x::Vector{Float64}=[-0.45, -0.17, -0.75, 0.0, 0.86, 0.28, -0.23, 0.57, 0.5, 0.62, 0.34, -0.17, -0.86, 0.27, 0.8, 0.0, -0.63, -0.19, -0.72, 1.0, 0.69, 0.62, -0.78, -0.01, 0.77, -0.14, 0.89, 0.55, -0.43, 0.58, -0.91, -0.45, 0.27, -0.08, -0.47, -0.96, -0.53, -0.58, -0.92, 0.0, -0.08, 0.69, 0.54, -0.26, -0.11, -0.5, 0.56, 0.1, 0.57, -0.49],
        L::Int=98,
    )
        N = length(x)
        x_padded = vcat(x, zeros(Float64, L - N))
        y = lti_filter(b, a, x_padded)

        return energy(y) # energy/power/rms
        4.975588325718084
    end
    out_6 = rozwiazanie_6()
end

#= Zadanie 7:
#* correct solution
wzmocnienie/przesunięcie fazowe systemu LTI
=#
# wariant a,b
begin
    function rozwiazanie_7_1(;
        b::Vector{Float64}=[0.005474558453857019, -0.021161090277260635, 0.05061387091651935, -0.07935032965393377, 0.09320417435199661, -0.07935032965393378, 0.05061387091651936, -0.021161090277260642, 0.0054745584538570214],
        a::Vector{Float64}=[1.0, -3.9521388331525866, 9.56265700381189, -14.849104406012085, 16.97637773213497, -13.784716624391761, 8.240829684422375, -3.1615178439273466, 0.7426218157666638],
        F::Vector{Float64}=[0.17, 0.23, 0.28, 0.48],
    )
        ϕ = [lti_phase(f, b, a) for f in F] # lti_amp/phase
        return mean(ϕ)
        -0.2671329372585373
    end
    out_7 = rozwiazanie_7_1()
end

# wariant zz,pp,k
begin
    function rozwiazanie_7_2(;
        zz::Vector{ComplexF64}=ComplexF64[0.9562125977134925-0.2926729710342487im, 0.6431481557756844+0.7657417643842708im, 0.9562125977134925+0.2926729710342487im, 0.6431481557756844-0.7657417643842708im, 0.9494540691139096-0.3139059901356441im, 0.6842289643189376+0.7292672516896902im, 0.9494540691139096+0.3139059901356441im, 0.6842289643189376-0.7292672516896902im, 0.9154123917349581-0.4025172705090848im, 0.8016754572775892+0.5977595345946131im, 0.9154123917349581+0.4025172705090848im, 0.8016754572775892-0.5977595345946131im],
        pp::Vector{ComplexF64}=ComplexF64[0.6064983645211166+0.787230370569453im, 0.9585153738231527-0.277256566212053im, 0.6064983645211166-0.787230370569453im, 0.9585153738231527+0.277256566212053im, 0.5363569616798629+0.7946801478013233im, 0.9537353148226411-0.25407497984321314im, 0.5363569616798629-0.7946801478013233im, 0.9537353148226411+0.25407497984321314im, 0.12166088338468982+0.6657516938860896im, 0.9290620344841986-0.14298335418965522im, 0.12166088338468982-0.6657516938860896im, 0.9290620344841986+0.14298335418965522im],
        k::Float64=0.28119079368088185,
        F::Vector{Float64}=[0.14, 0.17, 0.2],
    )
        a = poly_from_roots(pp)
        b = poly_from_roots(zz) .* k
        ϕ = [lti_amp(f, b, a) for f in F] # lti_amp/phase
        return mean(ϕ)
        0.5191880641409167
    end
    out_7 = rozwiazanie_7_2()
end

#= Zadanie 8:
#* correct solution
dyskretny system liniowy, stacjonarny, wyznacz stabilność (1 - tak/0 - na granicy/-1 - nie)
a - mianownik, b - licznik
=#
# wariant a,b
begin
    function rozwiazanie_8_1(;
        b::Vector{Float64}=[0.20496999142745434, -1.0248499571372718, 2.0496999142745436, -2.0496999142745436, 1.0248499571372718, -0.20496999142745434],
        a::Vector{Float64}=[1.0, -2.04328699669025, 2.205686901762555, -1.0925050400038385, 0.2795267139181644, 0.06196592669626961],
    )
        p = roots(a)
        radii = abs.(p)
        if all(radii .< 1)
            return 1.0
        elseif all(radii .<= 1) && any(radii .== 1)
            return 0.0
        else
            return -1.0
        end
        1.0
        stable
    end
    out_8 = rozwiazanie_8_1()
end

# wariant z,p,k
begin
    function rozwiazanie_8_2(;
        z::Vector{ComplexF64}=ComplexF64[0.9632631043548601+0.2685594753283481im, 0.9632631043548601-0.2685594753283481im, 0.9632631043548601+0.2685594753283481im, 0.9632631043548601-0.2685594753283481im, 0.9632631043548601+0.2685594753283481im, 0.9632631043548601-0.2685594753283481im, 0.9632631043548601+0.2685594753283481im, 0.9632631043548601-0.2685594753283481im, 0.9632631043548601+0.2685594753283481im, 0.9632631043548601-0.2685594753283481im],
        p::Vector{ComplexF64}=ComplexF64[0.9374548821453359-0.34810679962028074im, 0.9614739634351532+0.20095615586824897im, 0.9093219192859163+0.33766013615801665im, 0.9614739634351532-0.20095615586824897im, 0.8826481683827513-0.28771647831629865im, 0.9251374712100688+0.20909477610205304im, 0.8826481683827513+0.28771647831629865im, 0.9251374712100688-0.20909477610205304im, 0.892983758401963-0.2380350661011338im, 0.892983758401963+0.2380350661011338im],
        k::Float64=0.7753165952469473,
    )
        radii = abs.(p)
        if all(radii .< 1)
            return 1.0
        elseif all(radii .<= 1) && any(radii .== 1)
            return 0.0
        else
            return -1.0
        end
        0.0
        semistable
    end
    out_8 = rozwiazanie_8_2()
end

#= Zadanie 9:
#* correct solution
odpowiedzi impulsowe filtrów
=#
begin
    function rozwiazanie_9(;
        order::Int=74,
        fp::Float64=118.0,
        f1::Float64=10.03,
        f2::Float64=44.25,
        z::Vector{Int}=[40, 22, 74],
    )
        h = firwin_bs_I(order, f1 / fp, f2 / fp) # lp/hp/bp/bs
        h = h .* hamming(Int(order / 2)) # triang/hamming/hanning/blackman
        h_z = [h[i] for i in z]
        return sum(h_z)
        0.30685976285261213
    end
    out_9 = rozwiazanie_9()
end