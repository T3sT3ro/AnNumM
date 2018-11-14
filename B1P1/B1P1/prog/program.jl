# original method used by Henry Briggs in the 17th century
function briggs1(z::Complex{BigFloat}, maxiter=54)
    for i in 1:maxiter
        z = sqrt(z)
    end
    (z-1.0)*2.0^(maxiter)
end

# improved method proposed in the article by Al-Mohy
function briggs2(z::Complex{BigFloat}, maxiter=54)
    k = maxiter
    if angle(z) >= π/2 # square root of complex rotates by -ϕ/2
            z = sqrt(z)
            k = k-1;
    end
    z0 = z-1 # numerator
    z = sqrt(z)
    r = 1 + z
    for j in 1:(k-1)
        z = sqrt(z)
        r = r*(1+z)
    end
    (z0/r)*2.0^(maxiter)
end