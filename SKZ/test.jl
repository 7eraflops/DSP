begin
    using LinearAlgebra
    using Plots
end

function interp(m::Vector, s::Vector)::Function
    x -> begin
        sum::Float64 = 0.0
        Δt = m[2] - m[1]
        for i in eachindex(m)
            sum += s[i] * sinc((x - m[i]) / Δt)
        end
        return sum
    end
end


function rozwiazanie_1(;
    m::Vector{Float64}=[-4.7, -4.693, -4.686, -4.679, -4.672, -4.665, -4.658, -4.651, -4.644, -4.637, -4.63, -4.623, -4.616, -4.609, -4.602, -4.595, -4.588, -4.581, -4.574, -4.567, -4.56, -4.553, -4.546, -4.539, -4.532, -4.525, -4.518, -4.511, -4.504, -4.497, -4.49, -4.483, -4.476, -4.469, -4.462, -4.455, -4.448, -4.441, -4.434, -4.427, -4.42, -4.413, -4.406, -4.399, -4.392, -4.385, -4.378, -4.371, -4.364, -4.357, -4.35, -4.343, -4.336, -4.329, -4.322, -4.315, -4.308, -4.301, -4.294, -4.287, -4.28, -4.273, -4.266, -4.259, -4.252, -4.245, -4.238, -4.231, -4.224, -4.217, -4.21, -4.203, -4.196, -4.189, -4.182, -4.175, -4.168, -4.161, -4.154, -4.147, -4.14, -4.133, -4.126, -4.119, -4.112, -4.105, -4.098, -4.091, -4.084, -4.077, -4.07, -4.063, -4.056, -4.049],
    s::Vector{Float64}=[0.1341, 0.0009, 0.5749, 0.0896, 0.2634, 0.087, 0.6688, 0.3917, 0.5232, 0.2585, 0.7158, 0.9403, 0.9996, 0.575, 0.8647, 0.3559, 0.291, 0.9353, 0.747, 0.9326, 0.9747, 0.3191, 0.6705, 0.7298, 0.6991, 0.1069, 0.3095, 0.2437, 0.521, 0.7105, 0.3574, 0.8307, 0.3258, 0.3897, 0.5269, 0.0307, 0.5127, 0.8353, 0.3738, 0.8465, 0.0423, 0.3549, 0.7833, 0.172, 0.7738, 0.9263, 0.3024, 0.1555, 0.9879, 0.4182, 0.5143, 0.4226, 0.9599, 0.8262, 0.1294, 0.9552, 0.8977, 0.8782, 0.7279, 0.7376, 0.8788, 0.5089, 0.8415, 0.9273, 0.7598, 0.4772, 0.4259, 0.2801, 0.2404, 0.8741, 0.0215, 0.2299, 0.4406, 0.8229, 0.6719, 0.2625, 0.5542, 0.3627, 0.6211, 0.7575, 0.1887, 0.5563, 0.8977, 0.0621, 0.3229, 0.6237, 0.4767, 0.8019, 0.959, 0.0907, 0.8299, 0.3788, 0.241, 0.6729],
    t::Vector{Float64}=[-4.3227, -4.203, -4.3059, -4.686, -4.4676, -4.3605, -4.1883, -4.2842],
)
    f = interp(m, s)
    return sum(f.(t))
end
r_1 = rozwiazanie_1()

triangle_wave(t::Real)::Real = 2 / π * asin(sin(2π * t))

function rozwiazanie_2(;
    fp::Float64=465.65,
    t1::Float64=-0.09,
    N::Int=976,
)
    t = range(start=t1, step=1 / fp, length=N)
    g = triangle_wave
    y = 3.7 .* g.(4.8 .* t .- 2.5)
    return (sum(y) / length(y))
end
r_2 = rozwiazanie_2()

function dft(x::AbstractVector)::Vector
    N = length(x)
    w = [cispi(-2 * n / N) for n in 0:(N-1)]
    [
        sum((
            x[n+1] * w[n*f%N+1]
            for n in 0:(N-1)
        ))
        for f in 0:(N-1)
    ]
end

function freq2idx(f, fp, N)
    k = f * N / fp
    return Int(round(mod(k, N))) + 1
end

function rozwiazanie_3(;
    fp::Int=1452,
    x::Vector{ComplexF64}=ComplexF64[0.14+0.69im, -0.02+0.03im, 0.44+1.39im, -1.05-0.06im, 0.81+0.29im, -0.44-0.51im, 0.25-0.81im, -0.18-0.7im, 0.39+0.54im, -0.82+0.8im, 0.02-0.5im, -0.51-0.96im, 0.66-0.41im, 0.02-0.83im, -0.78+0.63im, -0.4+1.17im, 0.2+0.14im, 0.45-0.57im, 0.27-1.13im, 0.57+0.02im, 0.97+0.52im, 0.5+0.97im, -0.36+0.64im, 1.16-0.72im, 0.67+0.33im, 0.4+0.31im, -0.08-0.0im, 0.18+0.72im, 0.36+0.52im, -0.19+1.05im, 0.07+0.51im, 2.1+0.15im, -0.18-1.92im, 0.62-0.75im, 0.37+0.61im, 0.35-0.65im, 0.54-0.44im, 0.83-0.36im, -0.95-0.43im, -0.05-0.32im, -0.14-0.24im, -0.49+0.1im, -1.42+0.01im, 0.09+1.01im],
    f::Vector{Int}=[-297, -165, -99, 33, 495, 693],
)
    N = length(x)
    X = dft(x)
    A = [abs(X[freq2idx(freq, fp, N)]) / N for freq in f]
    return sum(A)
end
r_3 = rozwiazanie_3()

function conv(x::Vector, h::Vector)::Vector
    N = length(x)
    M = length(h)
    y = zeros(eltype(x), N + M - 1)
    for n in 1:N
        for m in 1:M
            y[n+m-1] += x[n] * h[m]
        end
    end
    return y
end

function rozwiazanie_4(;
    x::Vector{Float64}=[-1.53, 4.0, -2.5, 3.88, -1.65, -3.99, 0.24, 0.57, 2.37, -3.78, 4.59, 0.62, -1.03, -1.57, 4.03, 3.13, 2.46, 4.64, -4.07, -0.51, 3.09, 2.98, 1.53, 4.22, -0.56, 0.69, 4.36, 3.06, 4.96, 4.96, 1.54, 4.16, -2.39, 1.04, 2.22, -1.78, 0.8, -1.57, -1.75, -3.23, -0.38, -1.25, 4.48, -0.84, -3.84, 1.08, 0.0, -0.26, -3.99, 4.67, 1.76, 4.55, -3.94, -3.68, -2.55, 4.14, 1.36, -1.78, -1.3, 2.28, -3.79, -2.55, 0.7, -2.2, -2.6, -1.91, 3.26, 0.31, 0.2, -2.78, -4.9, 2.27, -2.09],
    h::Vector{Float64}=[-0.54, 2.18, -4.31, -1.52, -3.67, 1.25, -4.63, -3.61, -1.87, 3.65, -3.3, -2.14, -3.97, -1.6, -3.55, 3.39, 0.58, 1.79, 3.65, 1.93, 4.99, -4.3, -2.56, -4.01],
)
    y = conv(x, h)
    return sqrt(sum(abs2, y) / length(y))
end
r_4 = rozwiazanie_4()

function lti_filter(a::Vector, b::Vector, x::Vector)::Vector
    M = length(b) - 1
    K = length(a) - 1
    N = length(x)
    y = zeros(Float64, N)
    for n in 1:N
        for m in 0:M
            if n - m > 0
                y[n] += b[m+1] * x[n-m]
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

function rozwiazanie_5(;
    b::Vector{Float64}=[0.6901301478594652, -5.477897747610723, 19.065794944226962, -38.00417187606593, 47.452289104354726, -38.00417187606593, 19.06579494422696, -5.477897747610722, 0.6901301478594652],
    a::Vector{Float64}=[1.0, -7.205562937574693, 22.78751128545781, -41.31727406834236, 46.98266631589755, -34.31224083896434, 15.717681986016956, -4.129061402471896, 0.47627970115523915],
    x::Vector{Float64}=[-0.81, -0.12, 0.45, 0.55, 0.13, -0.95, -0.07, -0.9, 0.7, 0.54, 0.07, -0.07, -0.22, 0.37, 0.06, -0.35, 0.13, 0.47, 0.9, -0.71, -0.2, -0.45, -0.81, 0.84, -0.09],
    L::Int=64,
)
    x = vcat(x, zeros(L - length(x)))
    y = lti_filter(a, b, x)
    return sum(y) / length(y)
end
r_5 = rozwiazanie_5()
plot(r_5)


kronecker(n::Integer)::Real = ifelse(n == 0, 1, 0)
function firwin_bs_I(order::Integer, F1::Float64, F2::Float64)::Vector
    return [kronecker(Int(n)) - (2F2 * sinc(2F2 * n) - 2F1 * sinc(2F1 * n)) for n in -order/2:order/2]
end
hamming(M::Integer)::AbstractVector{<:Real} = [0.54 + 0.46cos(2π * n / (2M + 1)) for n = -M:M]

function rozwiazanie_6(;
    order::Int=32,
    fp::Float64=162.0,
    f1::Float64=3.24,
    f2::Float64=12.96,
    z::Vector{Int}=[19, 13, 28, 32, 31],
)
    h = firwin_bs_I(order, f1 / fp, f2 / fp)
    h .*= hamming(Int(order / 2))
    return sum([h[i] for i in z])
end
r_6 = rozwiazanie_6()



function rozwiazanie_7(;
    z::Vector{ComplexF64}=ComplexF64[0.8007730344245322+0.5989679017597914im, 0.8007730344245322-0.5989679017597914im, 0.9627435618567148+0.2704160389167881im, 0.9627435618567148-0.2704160389167881im],
    p::Vector{ComplexF64}=ComplexF64[0.21093504073151423+0.6783578605946071im, 0.21093504073151423-0.6783578605946071im, 1.3442975761847227+1.4808322074643971im, 0.1952953343794459-0.21513065725923242im],
    k::Float64=0.20099697813825557,
)
    radii = abs.(p)
    return radii
    if all(radii .< 1)
        return 1
    else
        return -1
    end
end
r_7 = rozwiazanie_7()
scatter(r_7)


function lti_phase(x, a::Vector, b::Vector)
    
end


function rozwiazanie_8(;
    zz::Vector{ComplexF64}=ComplexF64[0.1788605244924674-0.98387443953905im, 0.1788605244924674+0.98387443953905im, -0.19787415949776652-0.9802274312643228im, -0.19787415949776652+0.9802274312643228im, -1.0+0.0im],
    pp::Vector{ComplexF64}=ComplexF64[0.4587979265746062-0.8300138634274059im, 0.4587979265746062+0.8300138634274059im, 0.5552092908985323-0.5701035718407141im, 0.5552092908985323+0.5701035718407141im, 0.6538846599640508+0.0im],
    k::Float64=0.02257972299894665,
    F::Vector{Float64}=[0.02, 0.24, 0.25, 0.3, 0.4, 0.5],
)
    missing
end