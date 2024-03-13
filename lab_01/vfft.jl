using BenchmarkTools
using FFTW

bit_reverse(::Val{10}, num) = begin
    num = ((num & 0x3e0) >> 5) | ((num & 0x01f) << 5)
    num = ((num & 0x318) >> 3) | (num & 0x084) | ((num & 0x063) << 3)
    ((num & 0x252) >> 1) | (num & 0x084) | ((num & 0x129) << 1)
end

bit_reverse(::Val{64}, num) = bit_reverse(Val(32), (num & 0xffffffff00000000) >> 32) | (bit_reverse(Val(32), num & 0x00000000ffffff) << 32)
bit_reverse(::Val{32}, num) = bit_reverse(Val(16), (num & 0xffff0000) >> 16) | (bit_reverse(Val(16), num & 0x0000ffff) << 16)
bit_reverse(::Val{16}, num) = bit_reverse(Val(8), (num & 0xff00) >> 8) | (bit_reverse(Val(8), num & 0x00ff) << 8)
bit_reverse(::Val{8}, num) = bit_reverse(Val(4), (num & 0xf0) >> 4) | (bit_reverse(Val(4), num & 0x0f) << 4)
bit_reverse(::Val{4}, num) = bit_reverse(Val(2), (num & 0xc) >> 2) | (bit_reverse(Val(2), num & 0x3) << 2)
bit_reverse(::Val{3}, num) = ((num & 0x1) << 2) | ((num & 0x4) >> 2) | (num & 0x2)
bit_reverse(::Val{2}, num) = ((num & 0x2) >> 1) | ((num & 0x1) << 1)
bit_reverse(::Val{1}, num) = num

bit_reverse(::Val{9}, num) = begin
    num = ((num & 0x1e0) >> 5) | (num & 0x010) | ((num & 0x00f) << 5)
    num = ((num & 0x18c) >> 2) | (num & 0x010) | ((num & 0x063) << 2)
    ((num & 0x14a) >> 1) | (num & 0x010) | ((num & 0x0a5) << 1)
end

bit_reverse(::Val{31}, num) = begin
    bit_reverse(Val(15), num & 0x7fff0000 >> 16) | (num & 0x00008000) | (bit_reverse(Val(7), num & 0x00007fff) << 16)
end
bit_reverse(::Val{15}, num) = bit_reverse(Val(7), (num & 0x7f00) >> 8) | (num & 0x0080) | (bit_reverse(Val(7), num & 0x007f) << 8)
bit_reverse(::Val{7}, num) = bit_reverse(Val(3), (num & 0x70) >> 4) | (num & 0x08) | (bit_reverse(Val(3), num & 0x07) << 4)

function reverse_bit_order!(X, order)
    N = length(X)
    for i in 0:(N-1)
        j = bit_reverse(order, i)
        if i < j
            X[i+1], X[j+1] = X[j+1], X[i+1]
        end
    end
    X
end

function reverse_bit_order_double!(x, order)
    N = length(x)
    for i in 0:(N÷2-1)
        j = bit_reverse(order, i)
        if i < j
            # swap real part
            x[2*i+1], x[2*j+1] = x[2*j+1], x[2*i+1]
            # swap imaginary part
            x[2*i+2], x[2*j+2] = x[2*j+2], x[2*i+2]
        end
    end
    x
end

function my_fft(x)
    # Stop condition, the TF of an array of size 1 is this same array.
    if length(x) <= 1
        x
    else
        N = length(x)
        # Xᵒ contains the TF of odd terms and Xᵉ that of even terms.
        # The subtlety being that Julia's tablals start at 1 and not 0.
        Xᵒ = my_fft(x[2:2:end])
        Xᵉ = my_fft(x[1:2:end])
        factors = @. exp(-2im * π * (0:(N/2-1)) / N)
        [Xᵉ .+ factors .* Xᵒ; Xᵉ .- factors .* Xᵒ]
    end
end

function my_fft_2(x)
    N = length(x)
    order = Int(log2(N))
    @inbounds reverse_bit_order!(x, Val(order))
    n₁ = 0
    n₂ = 1
    for i = 1:order # i done the number of the column we are in.
        n₁ = n₂ # n₁ = 2ⁱ-¹
        n₂ *= 2 # n₂ = 2ⁱ

        step_angle = -2π / n₂
        angle = 0
        for j = 1:n₁ # j is the index in Xᵉ and Xᵒ
            factors = exp(im * angle) # z = exp(-2im*π*(j-1)/n₂)
            angle += step_angle # a = -2π*(j+1)/n₂

            # We combine the element j of each group of subarrays
            for k = j:n₂:N
                @inbounds x[k], x[k+n₁] = x[k] + factors * x[k+n₁], x[k] - factors * x[k+n₁]
            end
        end
    end
    x
end

function my_fft_3(x)
    N = length(x) ÷ 2
    order = Int(log2(N))
    @inbounds reverse_bit_order_double!(x, Val(order))

    n₁ = 0
    n₂ = 1
    for i = 1:order # i done the number of the column we are in.
        n₁ = n₂ # n₁ = 2ⁱ-¹
        n₂ *= 2 # n₂ = 2ⁱ

        step_angle = -2π / n₂
        angle = 0
        for j = 1:n₁ # j is the index in Xᵉ and Xᵒ
            re_factor = cos(angle)
            im_factor = sin(angle)
            angle += step_angle # a = -2π*j/n₂

            # We combine element j from each group of subarrays
            @inbounds for k = j:n₂:N
                re_xₑ = x[2*k-1]
                im_xₑ = x[2*k]
                re_xₒ = x[2*(k+n₁)-1]
                im_xₒ = x[2*(k+n₁)]
                x[2*k-1] = re_xₑ + re_factor * re_xₒ - im_factor * im_xₒ
                x[2*k] = im_xₑ + im_factor * re_xₒ + re_factor * im_xₒ
                x[2*(k+n₁)-1] = re_xₑ - re_factor * re_xₒ + im_factor * im_xₒ
                x[2*(k+n₁)] = im_xₑ - im_factor * re_xₒ - re_factor * im_xₒ
            end
        end
    end
    # We build the final version of the TF
    # N half the size of x
    # Special case n=0
    x[1] = x[1] + x[2]
    x[2] = 0

    step_angle = -π / N
    angle = step_angle
    @inbounds for n = 1:(N÷2)
        re_factor = cos(angle)
        im_factor = sin(angle)
        re_h = x[2*n+1]
        im_h = x[2*n+2]
        re_h_sym = x[2*(N-n)+1]
        im_h_sym = x[2*(N-n)+2]
        x[2*n+1] = 1 / 2 * (re_h + re_h_sym + im_h * re_factor + re_h * im_factor + im_h_sym * re_factor - re_h_sym * im_factor)
        x[2*n+2] = 1 / 2 * (im_h - im_h_sym - re_h * re_factor + im_h * im_factor + re_h_sym * re_factor + im_h_sym * im_factor)
        x[2*(N-n)+1] = 1 / 2 * (re_h_sym + re_h - im_h_sym * re_factor + re_h_sym * im_factor - im_h * re_factor - re_h * im_factor)
        x[2*(N-n)+2] = 1 / 2 * (im_h_sym - im_h + re_h_sym * re_factor + im_h_sym * im_factor - re_h * re_factor + im_h * im_factor)
        angle += step_angle
    end
    x
end

function my_fft_4(x)
    N = length(x) ÷ 2
    order = Int(log2(N))
    @inbounds reverse_bit_order_double!(x, Val(order))

    n₁ = 0
    n₂ = 1

    i = 1
    while i <= order # i done the number of the column we are in.
        n₁ = n₂ # n₁ = 2ⁱ-¹
        n₂ *= 2 # n₂ = 2ⁱ

        step_angle = -2π / n₂
        α = 2sin(step_angle / 2)^2
        β = sin(step_angle)
        cj = 1
        sj = 0
        j = 1
        while j <= n₁ # j is the index in Xᵉ and Xᵒ
            # We combine the element j from each group of subarrays
            k = j
            @inbounds while k <= N
                re_xₑ = x[2*k-1]
                im_xₑ = x[2*k]
                re_xₒ = x[2*(k+n₁)-1]
                im_xₒ = x[2*(k+n₁)]
                x[2*k-1] = re_xₑ + cj * re_xₒ - sj * im_xₒ
                x[2*k] = im_xₑ + sj * re_xₒ + cj * im_xₒ
                x[2*(k+n₁)-1] = re_xₑ - cj * re_xₒ + sj * im_xₒ
                x[2*(k+n₁)] = im_xₑ - sj * re_xₒ - cj * im_xₒ

                k += n₂
            end
            # We compute the next cosine and sine.
            cj, sj = cj - (α * cj + β * sj), sj - (α * sj - β * cj)
            j += 1
        end
        i += 1
    end
    # We build the final version of the TF
    # N half the size of x
    # Special case n=0
    x[1] = x[1] + x[2]
    x[2] = 0

    step_angle = -π / N
    α = 2sin(step_angle / 2)^2
    β = sin(step_angle)
    cj = 1
    sj = 0
    j = 1
    @inbounds while j <= (N ÷ 2)
        # We calculate the cosine and sine before the main calculation here to compensate for the first
        # step of the loop that was skipped.
        cj, sj = cj - (α * cj + β * sj), sj - (α * sj - β * cj)

        re_h = x[2*j+1]
        im_h = x[2*j+2]
        re_h_sym = x[2*(N-j)+1]
        im_h_sym = x[2*(N-j)+2]
        x[2*j+1] = 1 / 2 * (re_h + re_h_sym + im_h * cj + re_h * sj + im_h_sym * cj - re_h_sym * sj)
        x[2*j+2] = 1 / 2 * (im_h - im_h_sym - re_h * cj + im_h * sj + re_h_sym * cj + im_h_sym * sj)
        x[2*(N-j)+1] = 1 / 2 * (re_h_sym + re_h - im_h_sym * cj + re_h_sym * sj - im_h * cj - re_h * sj)
        x[2*(N-j)+2] = 1 / 2 * (im_h_sym - im_h + re_h_sym * cj + im_h_sym * sj - re_h * cj + im_h * sj)

        j += 1
    end
    x
end

@benchmark fft(a) setup = (a = rand(1024))
@benchmark fft!(a) setup = (a = rand(1024) .|> complex)
@benchmark rfft(a) setup = (a = rand(1024))

@benchmark my_fft(a) setup = (a = rand(1024))
@benchmark my_fft_2(a) setup = (a = rand(1024) .|> complex)
@benchmark my_fft_3(x) setup = (x = rand(1024))
@benchmark my_fft_4(x) setup = (x = rand(1048576))