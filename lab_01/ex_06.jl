##
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
