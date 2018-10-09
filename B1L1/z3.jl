W = [   x -> x^3 - 6x^2 + 3x - .149;
        x -> ((x-6)*x + 3)*x - .149 ]

x = 4.71
η = -14.636489

using Printf
for tp in [Float16, Float32, Float64]
    for w in W
        val = w(tp(x))
        @printf("%s:    w(%f) = %.10f   err = %.20f\n", tp, x, val, abs((val-η)/η))
    end
end
