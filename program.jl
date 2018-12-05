# Maksymilian Polarczyk 300791
# Pracownia z Analizy Numerycznej(M)
# zadania P2.13

using Printf
using Plots;

# P2D to wykorzystywana przeze mnie forma przechowywania punktów 2D w postaci macierzy nx2 o pierwszej kolumnie x.

# Iloraz różnicowy -funkcja licząca kolejne współczynniki. Ostatni wyraz to szukany iloraz.
# Zakładam brak wielokrotnych punktów xs (1pierwszej kolumny P2D)
function _diff(P2D) 
    b, n = Array(P2D[:, 2]), length(P2D[:, 1])
    for  i = 2:n
        for j = n:-1:i 
            b[j] = (b[j] - b[j-1])/(P2D[j, 1] - P2D[j-i+1, 1])
        end
    end
    return b 
end

# Zwraca wartość ilorazu.
function diff(P2D)
    last(_diff(P2D))
end

function interp_newton(P2D)
    return function(x)
        b, n = _diff(P2D), length(P2D[:, 1])-1
        sum, p = b[1], 1
        for i in 1:n
            p *= (x - P2D[i, 1])
            sum += b[i+1]*p
        end
        sum
    end
end

