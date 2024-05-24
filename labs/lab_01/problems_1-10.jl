begin
    using LinearAlgebra
    using Makie
end
## Problem 1.1

function factorial_r(n)
    if n == 1 || n == 0
        return 1
    end
    n * factorial_r(n - 1)
end

##

factorial_r(5)

## Problem 1.2

function factorial_i(n)
    x = 1
    for i in 1:n
        x *= i
    end
    return x
end

##

factorial_i(5)

## Problem 1.3

function is_even(n)
    if n % 2 == 0
        return true
    else
        return false
    end
end

##
is_even(4)
##
is_even(5)

## problem 1.4

function is_prime(n)
    if n <= 1
        return false
    end
    for i in 2:sqrt(n)
        if n % i == 0
            return false
        end
    end
    return true
end

##
is_prime(2)
##
is_prime(3)
##
is_prime(9)

## problem 1.5

function reverse_order(s)
    z = s[end:-1:1]
    println(z)
end

##
reverse_order("blue")

## problem 1.6

function is_pallindrome(s)
    for i in eachindex(s)
        if s[i] != s[end+1-i]
            return false
        end
    end
    return true
end

##
is_pallindrome("blue")
##
is_pallindrome("kayak")

## problem 1.7

function sierpinski_area(order)
    if order == 0
        return 1
    else
        area = 0.75^(order)
    end
end

##
sierpinski_area(0)
##
sierpinski_area(1)
##
sierpinski_area(2)
##
sierpinski_area(3)
##
sierpinski_area(4)

## problem 1.8

function newton_root(n, tol)
    x = n
    root = 0
    while true
        root = 0.5 * (x + (n / x))
        if abs(root - x) < tol
            break
        end
        x = root
    end
    return root
end

##
newton_root(4, 0.000001)
##
newton_root(18, 0.000001)

## problem 1.9

    

## problem 1.10

function is_circular_prime(i)
    if is_prime(i)
        if ndigits(i) < 1
            return true
        else
            dig = digits(i)
            for j in 1:ndigits(i)-1
                if !(is_prime(evalpoly(10, circshift(dig, j))))
                    return false
                end
            end
            return true
        end
    end
    return false
end

##
circular_primes = Int[]
for i in 1:1000000
    if is_circular_prime(i)
        push!(circular_primes, i)
    end
end

##
@show circular_primes
##
@show length(circular_primes)