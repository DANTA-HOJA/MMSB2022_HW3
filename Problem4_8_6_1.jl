using DifferentialEquations
using ModelingToolkit
using Plots
using LinearAlgebra
Plots.gr(lw=2)

# Convenience functions
hill(x, k) = x / (x + k)
hill(x, k, n) = hill(x^n, k^n)

@parameters k_1 k_2 k_3 k_4 k_5 n
@variables t A(t) B(t)
D = Differential(t)

eqs = [ D(A) ~ k_1 * hill(1, B, n) - (k_3 + k_5)* A,
        D(B) ~ k_2 + k_5 * A - k_4 * B]
@named sys = ODESystem(eqs)

params = Dict(k_1=>20.0, k_2=>5.0, k_3=>5.0, k_4=>5.0, k_5=>2.0, n=>4)

u0s = (Dict(A=>0.0, B=>0.0), 
       Dict(A=>0.5, B=>0.6),
       Dict(A=>0.17, B=>1.1),
       Dict(A=>0.25, B=>1.9),
       Dict(A=>1.85, B=>1.70))

tend = 1.5

sols = map(u0s) do u0
    prob = ODEProblem(sys, u0, tend, params)
    sol = solve(prob)
end;


p3 = plot()
for sol in sols
    plot!(p3, sol, linealpha=0.5, legend = nothing)
end

plot!(p3, xlabel="Time", ylabel="Concentration", title="Fig. 4.3A (Time series)")


p4 = plot()
for sol in sols
    plot!(p4, sol, vars=(1, 2), linealpha=0.7, legend = nothing)
end

plot!(p4, aspect_ratio=:equal, title="Fig. 4.3B (Phase plot)", xlabel="[A]", ylabel="[B]", ylims=(0.0, 2.0), xlims=(0.0, 2.0), size=(600, 600))