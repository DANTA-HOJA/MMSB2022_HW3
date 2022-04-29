using DifferentialEquations
using ModelingToolkit
using Plots
using LinearAlgebra
Plots.gr(lw=2)

# # Convenience functions
# hill(x, k) = x / (x + k)
# hill(x, k, n) = hill(x^n, k^n)

@parameters V k_2 k_3 k_4 k_5 K n
@variables t S1(t) S2(t)
D = Differential(t)

eqs = [ D(S1) ~ V - k_3*S1 + k_5*S2,
        D(S2) ~ k_2/(1+(S1/K)^n) - k_4*S2 - k_5*S2]
@named sys = ODESystem(eqs)

params = Dict(V=>1.0, k_2=>1.0, k_3=>1.0, k_4=>1.0, k_5=>1.0, K=>1.0, n=>1.0)

u0s = (Dict(S1=>0.0, S2=>0.0),
       Dict(S1=>0.1, S2=>0.3),
       Dict(S1=>0.17, S2=>0.6),
       Dict(S1=>0.25, S2=>1.3),
       Dict(S1=>1.5, S2=>1.70))

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

plot!(p4, aspect_ratio=:equal, title="Fig. 4.8.7 (Phase plot)", xlabel="[S1]", ylabel="[S2]", ylims=(0.0, 2.0), xlims=(0.0, 2.0), size=(600, 600))



# Nullclines

nullcline_S1(S1) = (params[k_3]*S1 - params[V])/params[k_5]
nullcline_S2(S1) = params[k_2]/((1+(S1/params[K])^params[n])*(params[k_4]+params[k_5]))

# rhs = [x.rhs for x in collect(eqs)]

# oop, iip = eval.(ModelingToolkit.build_function(rhs, [A, B], [k_1, k_2, k_3, k_4, k_5, n], t))

# oop([0, 0], [20, 5, 5, 5, 2, 4], 0.0)

# # function for vector field
# function ∂F(x, y, params; scale=20)
# 	du = oop([x, y], params, 0.0)
# 	return du ./ (norm(du)^0.5 * scale)
# end

# # Mesh points
# xx = [x for y in 0.0:0.1:2.0, x in 0.0:0.1:2.0];
# yy = [y for y in 0.0:0.1:2.0, x in 0.0:0.1:2.0];

# p1 = quiver(xx, yy, quiver=(x, y)->∂F(x, y, [20, 5, 5, 5, 2, 4]), line=(:lightblue))
# for sol in sols
#     plot!(p1, sol, vars=(1, 2), legend = nothing)
# end

# plot!(p1, aspect_ratio=:equal, title="Fig. 4.4A (Phase Plot with vector field)", 
#   xlabel="[A]", ylabel="[B]", xlim=(0.0, 2.0), ylim=(0.0, 2.0), size=(600, 600))

# [TextBook] Figure 4.5A（nullcline）
p487 = plot(aspect_ratio=:equal, title="Fig. 4.8.7, Phase plot with nullclines")

# Phase plots
for sol in sols
    plot!(p487, sol, vars=(2, 1), linealpha=0.7, lab=nothing)
end

# Parametric plotting for nullcline
plot!(p487, nullcline_S1, identity, 0.0, 2.0, label="S1 nullcline", line=(:black, :dot))
plot!(p487, nullcline_S2, identity, 0.0, 2.0, label="S2 nullcline", line=(:black, :dash))
plot!(p487, xlim=(0.0, 2.0), ylim=(0.0, 2.0), legend=:bottomright, size=(600, 600), xlabel="[S1]", ylabel="[S2]")