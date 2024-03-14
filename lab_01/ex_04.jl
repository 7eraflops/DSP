##
function is_prime(n, i)
    if n == 0 || n == 1
        return false
    end
    if n == i
        return true
    end
    if n % i == 0
        return false
    end
    i += 1
    is_prime(n, i)
end

##
is_prime(2, 2)
##
is_prime(3, 2)
##
is_prime(9, 2)