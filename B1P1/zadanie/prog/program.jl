# Maksymilian Polarczyk 300791
# Pracownia z Analizy Numerycznej(M)
# zadania P1.9

using Printf
using Plots;

# procedury kończące się ...R zwracają wartość log(a)/2^k
# oryginalna metoda Henrego Briggsa użyta w 17. wieku
function briggs1R(z, maxiter=54)
    for i in 1:maxiter
        z = sqrt(z)
    end
    z-1.0
end
function briggs1(z, maxiter=54)    briggs1R(z, maxiter)*2.0^(maxiter) end

# werjsa udoskonalona zaproponowana przez H. Al-Mohy-ego
function briggs2R(z, maxiter=54)
    k = maxiter
    if angle(z) >= π/2 # sprowadzenie 'z' na prawą połowę płaszczyzny zespolonej
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
    z0/r
end
function briggs2(z, maxiter=54)    briggs2R(z, maxiter)*2.0^(maxiter) end

# funkcja zwracająca błąd względny dla liczb zespolonych zdefiniowany jako |z-exact|/|exact|
function complexMaxRelErr(z, exact)
    if abs(exact) == 0
        abs(z)
    else
        abs(z - exact)/abs(exact)
    end
end

function epsilons(;maxk::Integer=54)
    # generujemy #samples równo odległych wartości dla ℜz oraz ℑz w zakresie
    # (0, +ε] oraz [-ε, 0) (+{0} dla ℑz) równoodległych o ε[i]/samples
    samples = 50
    ε = Float64[10e-16, 10e-8, 10e1, 10e8, 10e16]
    relErrors = zeros(Float64, 2, length(ε), length(ε))
    # dla każdej kombinacji zakresów ε
    for εRe in 1:length(ε)
        for εIm in 1:length(ε)
            # dla każdej próbki
            for k in 1:2*samples+1
                for j in 1:2*samples+1
                    # uzupełniamy tablicę wartościami błędów względnych
                    ℜz, ℑz = ε[εRe]*(k-samples-1)/samples, ε[εIm]*(k-samples-1)/samples
                    z = Complex{Float64}(Float64(ℜz)+Float64(ℑz)*im)
                    # odrzuć wartości dla argumentów z poza dziedziny
                    if !(ℜz <= 0 && ℑz == 0)
                        z1 = briggs1(z, maxk)
                        z2 = briggs2(z, maxk)
                        relErrors[1, εRe, εIm] = max(relErrors[1, εRe, εIm], 
                                                 complexMaxRelErr(z1, Complex{Float64}(log(z))))
                        relErrors[2, εRe, εIm] = max(relErrors[2, εRe, εIm], 
                                                 complexMaxRelErr(z2, Complex{Float64}(log(z))))
                    end
                end
            end
        end
    end
    @printf("Maksymalne błędy względne dla briggs1 i briggs2 względem bibliotecznej funkcji log.\n\t")
    for i in 1:5 @printf("\tℑz ≈ %.2e", ε[i]) end; println()
    println("> briggs1():")
    for i in 1:5 @printf("ℜz ≈ %.2e|", ε[i]);
        for j in 1:5 @printf("\t[%+.2e]", relErrors[1, i, j]) end; println() end # max. błędy względne dla alg1
    println("> briggs2():")
    for i in 1:5 @printf("ℜz ≈ %.2e|", ε[i]);
        for j in 1:5 @printf("\t[%+.2e]", relErrors[2, i, j]) end; println() end # max. błędy względne dla alg2
end

function epsilonK(;ε=10e1, samples::Integer=200, maxk::Integer=50)
    # generujemy #samples równo odległych wartości dla ℜz oraz ℑz w zakresie
    # (0, +ε] oraz [-ε, 0) (+{0} dla ℑz) równoodległych o ε/samples
    relErrors = zeros(Float64, 2, maxk)
    # da każdego k
    for k in 1:maxk
        # dla każdej próbki
        for i in 1:2*samples+1
            for j in 1:2*samples+1
                # uzupełniamy tablicę wartościami błędów względnych
                ℜz, ℑz = ε*(i-samples-1)/samples, ε*(i-samples-1)/samples
                z = Complex{Float64}(Float64(ℜz)+Float64(ℑz)*im)
                # odrzuć wartości dla argumentów z poza dziedziny
                if !(ℜz <= 0 && ℑz == 0)
                    z1, z2 = briggs1(z, k), briggs2(z, k)
                    relErrors[1, k] = max(relErrors[1, k], complexMaxRelErr(z1, Complex{Float64}(log(z))))
                    relErrors[2, k] = max(relErrors[2, k], complexMaxRelErr(z2, Complex{Float64}(log(z))))
                end
            end
        end
    end
    plot(title = "Complex case: relative error for given epsilon", xlabel="k", ylabel="rel. err.", yscale=:log10, fmt=:png,
          label = "briggs1", collect(1:maxk), relErrors[1, :], linecolor=:red)
    plot!(label = "briggs2", collect(1:maxk), relErrors[2, :], linecolor=:blue, linestyle=:dash)
    xticks!(0:5:maxk)
end

function increasingK(;logging=false)
    # liczby rzeczywiste; sprawdzanie dokładności obu algorytmów przy różnej ilości iteracji na 20 próbkach i 256 bitowej precyzji
    samples, maxk = 50, 250
    ε = BigFloat[10e-16, 10e2, 10e16]
    errors = zeros(BigFloat, 2, length(ε), maxk)
    for range in 1:length(ε)
        for k in 1:maxk
            # policz maksymalny błąd względny dla wszystkich próbek
            for j in 2:samples # wyklucz przypadek j=1 ε=1, gdzie log = 0, więc bład wzgl. = NaN lub +-Inf
                _0 = log(Complex{BigFloat}(ε[range]*j))
                _1 = briggs1(Complex{BigFloat}(ε[range]*j), k)
                _2 = briggs2(Complex{BigFloat}(ε[range]*j), k)
                errors[1, range, k] = max(errors[1, range, k], abs((real(_1)-real(_0))/real(_0)))
                errors[2, range, k] = max(errors[2, range, k], abs((real(_2)-real(_0))/real(_0)))
            end
            if logging @printf("[r, k: %d %3d] [%.e] [%.e]\n",range, k, errors[1, range, k], errors[2, range, k]) end
        end
    end
    
    plot(title = "Real case: relative error comparison", xlabel="k", ylabel="rel. err.", yscale=:log10, fmt=:png,
          label = "briggs1 10e-16", collect(1:maxk), errors[1, 1, :], linecolor=:black, linewidth=2)
    plot!(label = "briggs1 10e+2",  collect(1:maxk), errors[1, 2, :], linecolor=:red)
    plot!(label = "briggs1 10e+16", collect(1:maxk), errors[1, 3, :], linecolor=:magenta)

    plot!(label = "briggs2 10e-16", collect(1:maxk), errors[2, 1, :], linecolor=:cyan, linestyle=:dash)
    plot!(label = "briggs2 10e+2",  collect(1:maxk), errors[2, 2, :], linecolor=:blue, linestyle=:dash)
    plot!(label = "briggs2 10e+16", collect(1:maxk), errors[2, 3, :], linecolor=:green, linestyle=:dash)
    xticks!(0:20:maxk)
end

function bigK(;samples=20, maxk=500)
    setprecision(512) do # na 256 precyzji nie otrzymamy lepszych wyników niż ~10e-75
        ε = BigFloat[10e-8, 1, 10e8]
        errors2 = zeros(BigFloat, length(ε), maxk)
        for range in 1:length(ε)
            for k in 1:maxk
                # policz maksymalny błąd względny dla wszystkich próbek
                for j in 2:samples # wyklucz przypadek j=1 ε=1, gdzie log = 0, więc bład wzgl. = NaN lub +-Inf
                    _0 = log(Complex{BigFloat}(ε[range]*j))
                    _2 = briggs2(Complex{BigFloat}(ε[range]*j), k)
                    errors2[range, k] = max(errors2[range, k], abs((real(_2)-real(_0))/real(_0)))
                end
            end
        end

        plot(title = "Real case: relative error for big k", xlabel="k", ylabel="rel. err.", fmt=:png,
          label = "briggs2 10e-8", collect(1:maxk), errors2[1, :], yscale=:log10)
        plot!(label = "briggs2 10e-1", collect(1:maxk), errors2[2, :])
        plot!(label = "briggs2 10e8", collect(1:maxk), errors2[3, :])
        xticks!(0:50:maxk)
    end
end
