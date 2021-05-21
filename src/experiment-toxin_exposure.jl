export single_toxin_exposure_experiment
export toxin_exposure_experiment
export plot_toxin_exposure_experiment

function toxin_exposure_experiment(
    cell::Cell,
    amount_toxins::Array{Float64,1},
    exposure_times::Array{Float64,1},
    n_points::Int64 = 100,
    stable::Bool = true,
    n_cells::Int64 = 1,
    control::Bool = false;
    kwargs...,
)

    experiments = control ? [:Control, :Normal] : [:Normal]
    dfs = []
    sort!(exposure_times)
    for amount_toxin ∈ amount_toxins, experiment ∈ experiments

        # Setout.
        cell_to_evaluate =
            experiment == :Control ? Cell(cell; filamentation_rate = 0.0) : cell
        max_τ = exposure_times[end]

        # Create toxin signal.
        time, toxin = toxin_entry(0.0, max_τ, n_points, amount_toxin, stable)

        # Stochastic or ordinary.
        if n_cells > 1

            states = repeat([0.0], length(exposure_times))
            for n_cell ∈ 1:n_cells
                sol, sol_cell = simulate(
                    cell_to_evaluate,
                    time,
                    toxin;
                    stochastic = true,
                    status = false,
                    kwargs...,
                )
                status = _check_cell_status(sol_cell, exposure_times)
                states += map(x -> x ∈ [:Normal, :Stressed] ? 1.0 : 0.0, status)
            end
            states ./= n_cells
        else
            sol, sol_cell = simulate(
                cell_to_evaluate,
                time,
                toxin;
                sthocastic = false,
                status = false,
                kwargs...,
            )
            states = _check_cell_status(sol_cell, exposure_times)
        end
        push!(
            dfs,
            DataFrame(
                exposure_time = exposure_times,
                amount_toxin = amount_toxin,
                state = states,
                experiment = experiment,
            ),
        )
    end

    return reduce(vcat, dfs)
end

function plot_toxin_exposure_experiment(df, width = 300, height = 300)

    state_type = eltype(df[!, :state])
    df[!, :experiment] = map(
        x -> x == :Control ? "Without filamentation" : "Filamentation",
        df[!, :experiment],
    )

    if state_type == Symbol
        return @vlplot(
            data = df,
            width = width,
            height = height,
            mark = :rect,
            x = {"exposure_time:q", title = "Exposure time", axis = {labelAngle = 0}},
            y = {"amount_toxin:q", sort = "descending", title = "Amount toxin"},
            column = {"experiment:n", title = "Experiment", sort = [:Control, :Normal]},
            color = {"state:n", title = "Cell state", sort = [:Normal, :Stressed, :Dead]}
        )
    elseif state_type == Float64
        return @vlplot(
            data = df,
            width = width,
            height = height,
            mark = :rect,
            x = {"exposure_time:o", title = "Exposure time", axis = {labelAngle = 0}},
            y = {"amount_toxin:o", sort = "descending", title = "Amount toxin"},
            column = {"experiment:n", title = "Experiment", sort = [:Control, :Normal]},
            color = {"state:q", title = "Proportion of\nSurvivors"}
        )
    end
    return nothing
end
