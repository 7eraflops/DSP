
##
function sierpinski_area(order)
    if order == 0
        return 1
    else
        area = 0.75 ^ (order)
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