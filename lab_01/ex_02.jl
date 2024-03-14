##
function factorial_i(n)
    x = 1
    for i in 1:n
        x *= i
    end
    return x
end
##
factorial_i(5)