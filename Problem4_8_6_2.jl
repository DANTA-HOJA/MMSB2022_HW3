using DifferentialEquations
using ModelingToolkit
using Plots
using LinearAlgebra
Plots.gr(lw=2)

# Convenience functions
hill(x, k) = x / (x + k)
hill(x, k, n) = hill(x^n, k^n)

@parameters k_1 k_2 k_3 A B
@variables t X(t) Y(t)
D = Differential(t)

eqs = [ D(X) ~ A*k_1*X - k_2*X*Y,
        D(Y) ~ k_2*X*Y - k_3*Y]
@named sys = ODESystem(eqs)

params = Dict(k_1=>1.0, k_2=>1.0, k_3=>1.0, A=>1.0, B=>1.0)

u0s = (Dict(X=>0.0, Y=>0.0),
       Dict(X=>0.1, Y=>0.3),
       Dict(X=>0.17, Y=>0.6),
       Dict(X=>0.25, Y=>1.3),
       Dict(X=>1.5, Y=>1.70))

tend = 10

sols = map(u0s) do u0
    prob = ODEProblem(sys, u0, tend, params)
    sol = solve(prob)
end;


# p3 = plot()
# for sol in sols
#     plot!(p3, sol, linealpha=0.5, legend = nothing)
# end

# plot!(p3, xlabel="Time", ylabel="Concentration", title="Fig. 4.3A (Time series)")


p4 = plot()
for sol in sols
    plot!(p4, sol, vars=(1, 2), linealpha=0.7, legend = nothing)
end

plot!(p4, aspect_ratio=:equal, title="Fig. 4.8.6-1 (Phase plot)", xlabel="[X]", ylabel="[Y]", ylims=(0.0, 2.0), xlims=(0.0, 2.0), size=(600, 600))