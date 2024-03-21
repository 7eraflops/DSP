# module CPS

using LinearAlgebra

author = Dict{Symbol,String}(
    :index => "417096",
    :name => "Michał Charaszkiewicz",
    :email => "mchar@student.agh.edu.pl",
    :group => "0",
)

# Sygnały ciągłe
cw_rectangular(t::Real; T=1.0)::Real = missing
cw_triangle(t::Real; T=1.0)::Real = missing
cw_literka_M(t::Real; T=1.0)::Real = missing
cw_literka_U(t::Real; T=1.0)::Real = missing

ramp_wave(t::Real)::Real = 2 * rem(t, 1, RoundNearest)
sawtooth_wave(t::Real)::Real = -2 * rem(t, 1, RoundNearest)
triangular_wave(t::Real)::Real = ifelse(modx < 1 / 2, 4modx - 1, -4modx + 3)
square_wave(t::Real)::Real = ifelse(rem(t, 1) < 0.5, 1, -1)
pulse_wave(t::Real, ρ::Real)::Real = ifelse(rem(t, 1) < ρ, 1, 0)
impulse_repeater(g::Function, t1::Real, t2::Real)::Function = missing

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function sawtooth_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function triangular_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function square_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function pulse_wave_bl(t; ρ=0.2, A=1.0, T=1.0, band=20.0)
    missing
end


function impuse_repeater_bl(g::Function, t0::Real, t1::Real, band::Real)::Function

    function rand_siganl_bl(f1::Real, f2::Real)::Function
        missing
    end

    function rand_siganl_bl(f1::Real, f2::Real)::Function
        missing
    end



    # Sygnały dyskretne
    kronecker(n::Integer)::Real = ifelse(n == 0, 1, 0)
    heaviside(n::Integer)::Real = ifelse(n < 0, 0, 1)

    # Okna
    rect(N::Integer)::AbstractVector{<:Real} = missing
    triang(N::Integer)::AbstractVector{<:Real} = missing
    hanning(N::Integer)::AbstractVector{<:Real} = missing
    hamming(N::Integer)::AbstractVector{<:Real} = missing
    blackman(N::Integer)::AbstractVector{<:Real} = missing

    # Parametry sygnałów
    mean(x::AbstractVector)::Number = missing
    peak2peak(x::AbstractVector)::Real = missing
    energy(x::AbstractVector)::Real = missing
    power(x::AbstractVector)::Real = missing
    rms(x::AbstractVector)::Real = missing

    function running_mean(x::AbstractVector, M::Integer)::Vector
        missing
    end

    function running_energy(x::AbstractVector, M::Integer)::Vector
        missing
    end

    function running_power(x::AbstractVector, M::Integer)::Vector
        missing
    end



    # Próbkowanie
    function interpolate(
        t::Real;
        m::AbstractVector,
        s::AbstractVector,
        kernel::Function=sinc
    )::Real
        missing
    end

    # Kwantyzacja
    quantize(x::Real; L::AbstractVector)::Real = missing
    SQNR(N::Integer)::Real = missing
    SNR(Psignal, Pnoise)::Real = missing


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

end
