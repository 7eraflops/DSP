##
function newton_root(n, tol)
    x = n
    count = 0
    root = 0
    while true
        count += 1
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