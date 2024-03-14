function factorial_r(n)
    if n == 1 || n == 0
        return 1
    end
    n * factorial_r(n - 1)
end

factorial_r(5)