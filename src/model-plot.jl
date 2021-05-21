export simulation_plot

# function simulation_plot(
#     cell::Cell,
#     time::Array{Float64, 1},
#     toxin::Array{Float64, 1},
#     n_cells::Int64=0,
#     sde_alpha=0.5
# )

#     # Run simulalions. Control, normal and stochastic.
#     sol_c, sol_c_cell = simulate(Cell(cell; filamentation_rate=0.0), time, toxin)
#     sol, sol_cell = simulate(cell, time, toxin)
#     sols_sde = [simulate(cell, time, toxin; stochastic=true)[1] for _ in 1:n_cells];

#     # Establish defaults for plotting.
#     ll = @layout [A; B]
#     ll = Dict(:layout => grid(2, 1), :ga => 0.5, :leg => false)
#     ode = Dict(:lw => 2.0, :lc => :black)
#     sde = Dict(:lw => 1.0, :lc => "#676BA0", :alpha => sde_alpha)
#     bkg_control = Dict( :color => "#FF811A", :alpha=>0.5)
#     bkg_normal = Dict(:color => "#72D861", :alpha=>0.5)
#     lines = Dict(:ls => :dash, :lw=> 1.5)
#     points = Dict(:msc => :white, :ms => 5)

#     # Extract information for positions of geometries.
#     ## τs.
#     τ = Dict(
#         :sos => [sol_c_cell.τ_filamentation, :blue, L"\tau_{sos}"],
#         :kill_c => [sol_c_cell.τ_kill, "#0c1413", L"\tau_{kill}"],
#         :kill => [sol_cell.τ_kill, :red, L"\tau_{kill}"]
#     )

#     ## Axis.
#     x_max = sol.t[end] + 1
#     y_max = cell.kill_threshold * 1.06

#     # Initialize empty plot.
#     p = plot(; ll...)

#     # Geometric annotations part 1.
#     hline!(p, [cell.filamentation_threshold], c=τ[:sos][2], subplot=1; lines...)
#     hline!(p, [cell.kill_threshold], c=τ[:kill][2], subplot=1; lines...)
#     for i in 1:2

#         # Add shadow regions.
#         if !isinf(τ[:kill_c][1])
#             vspan!(p, [τ[:sos][1], τ[:kill_c][1]], subplot=i; bkg_control...)

#             max = isinf(τ[:kill][1]) ? x_max - 1 : τ[:kill][1]
#             vspan!(p, [τ[:kill_c][1], max], subplot=i; bkg_normal...)
#         end

#         # Add vertical lines.
#         for (k, (value, color, txt)) in τ

#             if !isinf(value)
#                 vline!(p, [value], color=color, subplot=i; lines...)
#             end
#         end

#     end

#     # Add simulation lines.
#     map(x->plot!(p, x; sde...), sols_sde)
#     plot!(p, sol_c, ls=:dash; ode...)
#     plot!(p, sol; ode...)

#     # Geometric annotations part 2 (points).
#     for i in 1:2
#         for (k, (value, color, txt)) in τ
#             if !isinf(value)
#                 v = k == :kill ? sol(value)[i] : sol_c(value)[i]
#                 scatter!(p, [value], [v], c=color, subplot=i; points...)
#                 annotate!(p, value, y_max, text(txt, color, 12), subplot=1)
#             end
#         end
#     end

#     # Text Annotations
#     xlabel!(p, "", subplot=1)
#     xlabel!(p, "Time (a.u.)", subplot=2)
#     ylabel!(p, "Internal toxin", subplot=1)
#     ylabel!(p, "Length", subplot=2)
#     xlims!(p, (0, x_max))

#     return p
# end;

function simulation_plot(
    cell::Cell,
    time::Array{Float64, 1},
    toxin::Array{Float64, 1},
    trajectories::Int64=0;
    kwargs...
)

    # Run simulalions. Control, normal and stochastic.
    sol_c, sol_c_cell = simulate_filamentation(Cell(cell; filamentation_rate=0.0), time, toxin, trajectories=0)
    sol, sol_cell = simulate_filamentation(cell, time, toxin, trajectories=0)
    sols_sde, _ = simulate_filamentation(cell, time, toxin, trajectories=trajectories; kwargs...)
    sols_sde = EnsembleSummary(sols_sde, time)


    # Establish defaults for plotting.
    ll = Dict(:layout => (2, 1), :ga => 0.4, :legend => false, :link => :x)
    ode = Dict(:lw => 2.0, :lc => :black)
    sde = Dict(:lw => 1.0, :fillcolor => "#676BA0", :linecolor => "#676BA0", :alpha => 1.0)
    bkg_control = Dict( :color => "#FF811A", :alpha=>0.5)
    bkg_normal = Dict(:color => "#72D861", :alpha=>0.5)
    lines = Dict(:ls => :dash, :lw=> 1.5)
    points = Dict(:msc => :white, :ms => 5)

    # Extract information for positions of geometries.
    ## τs.
    τ = Dict(
        :sos => [sol_c_cell.τ_filamentation, :blue, L"\tau_{sos}"],
        :kill_c => [sol_c_cell.τ_kill, "#0c1413", L"\tau_{kill}"],
        :kill => [sol_cell.τ_kill, :red, L"\tau_{kill}"]
    )

    ## Axis.
    x_max = sol.t[end] + 1
    y_max = cell.kill_threshold * 1.07
    
    # Initialize empty plot.
    p = plot(; ll...)

    # Geometric annotations part 1.
    hline!(p, [cell.filamentation_threshold], c=τ[:sos][2], subplot=1; lines...)
    hline!(p, [cell.kill_threshold], c=τ[:kill][2], subplot=1; lines...)
    for i in 1:2
    
        # Add shadow regions.
        if !isinf(τ[:kill_c][1])
            vspan!(p, [τ[:sos][1], τ[:kill_c][1]], subplot=i; bkg_control...)
    
            max = isinf(τ[:kill][1]) ? x_max - 1 : τ[:kill][1]
            vspan!(p, [τ[:kill_c][1], max], subplot=i; bkg_normal...)
        end
    
        # Add vertical lines.
        for (k, (value, color, txt)) in τ
    
            if !isinf(value)
                vline!(p, [value], color=color, subplot=i; lines...)
            end
        end
    end

    plot!(p, sols_sde; sde...)
    plot!(p, sol_c, ls=:dash; ode...)
    plot!(p, sol; ode...)

    # Geometric annotations part 2 (points).
    for i in 1:2
        for (k, (value, color, txt)) in τ
            if !isinf(value)
                v = k == :kill ? sol(value)[i] : sol_c(value)[i]
                scatter!(p, [value], [v], c=color, subplot=i; points...)
                annotate!(p, value, y_max, text(txt, color, 12), subplot=1)
            end
        end
    end

    # Text Annotations
    xlabel!(p, "", subplot=1)
    xlabel!(p, "Time (a.u.)", subplot=2)
    ylabel!(p, "Internal toxin", subplot=1)
    ylabel!(p, "Length", subplot=2)
    xlims!(p, (0, x_max))
    
    return p
end