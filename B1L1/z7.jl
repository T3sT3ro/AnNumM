function stringfloat(str)
    s = str[1] == '1' ? -1 : 1

    e = 0
    for i = 2:12
        e = e<<1 + (str[i] == '1' ? 1 : 0)
    end
    e -= 1023 #BIAS
    m = t = Float64(1)
    for i = 13:64
        t/=2; m += str[i] == '1' ? t : 0
    end

    s*m*2.0^e
end
