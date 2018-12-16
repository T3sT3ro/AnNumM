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
#   kolejne wartości ilorazu róznicowego są współczynnikami wielomianu interpolacyjnego
        b, n = _diff(P2D), length(P2D[:, 1])-1
        sum, p = b[1], 1
        for i in 1:n
            p *= (x - P2D[i, 1])
            sum += b[i+1]*p
        end
        sum
        
#   a, b, n = P2D[:, 1], _diff(P2D), length(P2D[:, 1])
#   polynomial = 0
#   for i in 1:n
#       polynomial += b[i]*polynomial*poly([a[i]])
#   end
#   polynomial(x)
    end
end

##################################################################################################################################

# zwraca kopię tablicy punktów o wielkości a*k+m, gdzie k należy do całkowitych,
# kopiując ostatni punkt odpowiednią ilość razy
function _P2D_fixlength(P2D, a, m)
    lacking = a-((length(P2D[:, 1]) - m)%a)
    if lacking == 0
        return Array(P2D)
    else
        _P2D = Array(P2D)
        for i in 1:lacking
            _P2D = vcat(_P2D, reshape(P2D[length(P2D[:, 1]), :], (1, :)))
        end
        return _P2D
    end
end

##################################################################################################################################

# https://www.cs.helsinki.fi/group/goa/mallinnus/curves/curves.html
# P = a₂t²+a₁t¹+a₀ gdzie a to wektory, t to parametr[0:1]
# kwadratowa krzywa sklejana
function interp_qubic_param_spline(P2D)
    # każdy kawałek przybliżany jest trzema punktami. Potrzeba 2k+1 punktów
    P2D_2k1 = _P2D_fixlength(P2D, 2, 1)
    # zwracana funkcja przyjmuje t z zakresu 
    # dla n punktów (x_1, x_2, ..., x_n) jest n-1 zakresów dla funkcji parametrycznej
    # t  z zakresu t:[1..n]
    return function(t)
        p = []
        if t == 1
            t = 0
            p = P2D[1:3, :]
        else
            p = P2D[2*(ceil(Int, (t-1)/2))-1:2*(ceil(Int, (t-1)/2))+1, :]
            t = (t - 1 - 2*((ceil(Int, (t-1)/2))-1))/2
        end
        if t == 0
            p[1, :]
        elseif t == 1
            p[3, :]
        else 
            a = [p[1, :], 0, 0]
            a[3] = ((p[2, :]-p[1, :])-t*(p[3, :]-p[1, :]))/(t*(t-1))
            a[2] = p[3, :]-p[1, :] - a[3]
            (a[3]*t + a[2])*t + a[1] # prosty schemat Hornera
        end
    end
end


##################################################################################################################################

#Wzór funkcji sklejanej
function _spline_equation(x,P2D,h,M)
    n = length(P2D[:, 1])-1
    #Wybor przedzialu dla wzoru
    k = 1
    for i in 1:n
        if x <= P2D[i+1, 1]
            k = i
            break
        end
    end 
    # zgodnie ze wzorem
    res =  (M[k] * (P2D[k+1, 1] - x)^3 + M[k+1] * ((x-P2D[k, 1])^3))/6
    res += (P2D[k, 2] - M[k]*(h[k+1]^2)/6) * (P2D[k+1, 1] - x) 
    res += (P2D[k+1, 2] - M[k+1] * (h[k+1]^2)/6) * (x-P2D[k, 1])
    return res * (1 / h[k+1]) 
end

# macierz drugich pochodnych funkcji sklejanej
function _spline_M_derivatives(P2D)
    n = length(P2D[:, 1])-1
    #############
    # d[k]: ilorazy różnicowe
    d=zeros(n+1)
    for k in 1 : (n-1)
        d[k+1] = 6 * diff(P2D[k:k+2, :])
    end
    #############
    h = zeros(n+1)
    for k in 1:n 
        h[k+1] = P2D[k+1, 1]-P2D[k, 1]
    end
    #############    
    lambda = zeros(n+1)
    for k in 0:(n-1) 
        lambda[k+1] = h[k+1] / (h[k+1]+h[k+1+1]) 
    end
    #############
    # M[k]: druga pochodna
    M, q, p, u=zeros(n+1),zeros(n+1),zeros(n+1), zeros(n+1)
    for k in 1:(n-1)
        p[k+1] = lambda[k+1] * q[k] + 2
        q[k+1] = (lambda[k+1] - 1) / p[k+1]
        u[k+1] = (d[k+1] - lambda[k+1] * u[k]) / p[k+1]
    end
    #############
    M[n]=u[n]
    for k in (n-2):-1:1
        M[k+1]=u[k+1] + q[k+1] *M[k+2]
    end
    return  (h,M)
end

# zwraca naturalną funkcję sklejanoą 3 stopnia 
function _interp_spline(P2D) 
    h, M = _spline_M_derivatives(P2D)
    return x -> _spline_equation(x, P2D, h, M)
end

##################################################################################################################################


# macierz pochodnych dla parametrycznych funkcji sklejanych, liczenie dwóch funkcji naraz
# ts: tablica wartości parametru t (funkcje fx(t) i fy(t) interpolowane w tym samym czasie, bo te same węzły)
function _spline_M_derivatives_2(P2D,ts)
    n = length(P2D[:, 1])-1
    #############
    # d[k]: ilorazy różnicowe x i y
    d = zeros((n+1,n+1))
    for k in 1 : (n-1)
        d[k+1, 1] = 6 * diff( hcat(ts[k:k+2], P2D[k:k+2, 1]))
        d[k+1, 2] = 6 * diff( hcat(ts[k:k+2], P2D[k:k+2, 2]))
    end
    #############
    h = zeros(n+1)
    for k in 1:n 
        h[k+1] = ts[k+1]-ts[k]
    end
    #############
    lambda = zeros(n+1)
    for k in 0:(n-1) 
        lambda[k+1] = h[k+1] / (h[k+1]+h[k+2]) 
    end
    #############
    # M[k]: druga pochodna
    q, p = zeros(n+1), zeros(n+1)
    u, M = zeros((n+1, n+1)), zeros((n+1, n+1))
    for k in 1:(n-1)
        p[k+1] = lambda[k+1] * q[k] + 2
        q[k+1] = (lambda[k+1] - 1) / p[k+1]
        u[k+1] = (d[k+1] - lambda[k+1] .* u[k]) / p[k+1]
    end
    #############
    M[n]=u[n]
    for k in (n-2):-1:1
        M[k+1]=u[k+1] + q[k+1] .* M[k+2]
    end
    return  (h,M)
end


# zwraca parę tablicę wartości dla funkcji parametrycznych interpolujących kształt zadany przez tablicę P2D
# przybliżanie za pomocą 5000 równoodległych punktów na krzywej parametrycznej.
function splain_param_interp_P2D(P2D)
    n = length(P2D[:, 1])
    fx = _interp_spline(hcat([1:n;],P2D[:, 1]))
    fy = _interp_spline(hcat([1:n;],P2D[:, 2]))
    points = [range(1, stop=n, length=5000);]
    return (map(fx,points), map(fy,points))
end

# zwraca tuple (fx, fy) wzorów parametrycznych funkcji sklejanych 
function _interpSplineParam(P2D,ts)
    h,M = _spline_M_derivatives_2(P2D,ts)
    return (x->_spline_equation(x, hcat(ts, P2D[:, 1]), h, M[:, 1]),
            y->_spline_equation(y, hcat(ts, P2D[:, 2]), h, M[:, 2]))
end

function make_param_curve(fx,fy,ts,xss, mytitle)
    xs, ys= map(fx,ts), map(fy,ts)
    SplineX,SplineY = _interpSplineParam(hcat(xs, ys),ts)

    # Krzywa parametryczna
    shape = begin
        plot(map(SplineX, xss), map(SplineY, xss), label = "sklejana", line=(0.8, [:blue]))
        plot!(map(fx,xss), map(fy,xss), label="parametryczna", line=(0.8, [:orange]), title="wykres")
    end

    # Funkcje błędow i ich wykres
    err = begin
        plot(xss, map(t->abs(SplineX(t) - fx(t)), xss), label = "|Sx(t) - fx(t)|", line=(0.8, [:blue]))
        plot!(xss, map(t->abs(SplineY(t) - fy(t)), xss), label = "|Sy(t) - fy(t)|", line=(0.8, [:orange]), title="błąd")
    end
    
    # norma jako max
    norm =  reduce(max, (map(t->sqrt((fx(t) - SplineX(t))^2 + (fy(t) - SplineY(t))^2), xss)))

    return (plot(shape, err, layout = (1, 2), size=(900, 300)),
            norm)
end

function plot_shape_interpolate(P2D)
    xs, ys = splain_param_interp_P2D(P2D)
    curve = begin
        scatter(markersize=2, markercolor=:blue, markerstrokealpha=0,P2D[:, 1], P2D[:, 2], label="dane")
        plot!(xs, ys, size=(900, 300), aspect_ratio=:equal, label="funkcja sklejana", title="funkcja")
    end
    plot(curve, layout=(1, 2))
end