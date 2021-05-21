export antitoxin_experiment
export plot_antitoxin_experiment

function antitoxin_experiment(
    cell::Cell,
    time,
    toxin;
    n_cells = 100,
    σs = [0.0],
    seed = nothing,
    normalize::Bool = true,
    kwargs...,
)

    colnames = [:Normal, :Stressed, :Dead]
    df = DataFrame(time = time)
    for col in colnames
        df[!, col] .= 0
    end

    dfs = []
    antitoxins_df = []
    for σ ∈ σs

        tmp_df = deepcopy(df)
        tmp_df[!, :sigma] .= σ

        D = Normal(cell.antitoxin, σ)
        if seed ≠ nothing
            Random.seed!(seed)
        end

        antitoxins = map(x -> x < 0 ? 0 : x, rand(D, n_cells))

        push!(antitoxins_df, DataFrame(sigma = σ, antitoxin = antitoxins))
        for antitoxin ∈ antitoxins

            tmp_cell = Cell(cell; antitoxin = antitoxin)
            status = simulate(tmp_cell, time, toxin; status = true, kwargs...)[3]

            for (i, s) ∈ enumerate(status)
                tmp_df[i, s] += 1
            end
        end
        push!(dfs, tmp_df)
    end

    df = reduce(vcat, dfs)
    df = stack(df, Not([:time, :sigma]))
    df[!, :μ] .= cell.antitoxin

    if normalize
        df[!, :value] = df[!, :value] ./ n_cells
    end

    return df, reduce(vcat, antitoxins_df)
end

function plot_antitoxin_experiment(df, plot_heatmap::Bool = false)

    if plot_heatmap
        new_df = @linq df |>
            transform(value = :value / maximum(:value) * 100.0) |>
            where(:variable .== "Dead") |>
            transform(value = 100.0 .- :value)

        m = Matrix(unstack(new_df, :time, :sigma, :value)[!, Not(:time)])'
        p = heatmap(unique(new_df[!, :time]), unique(new_df[!, :sigma]), m)
        xlabel!(p, "Exposure time")
        ylabel!(p, "Population variability")
        title!(p, "Survival probability (%)")
        return p

    end

    title = string("μ = ", df[!, :μ][1], " with σ values of:")
    p = @vlplot(
        :area,
        data = df,
        x = {"time:q", title = "Time"},
        y = {"value:Q", stack = :normalize, title = "Proportion"},
        color =
            {"variable:n", title = "Cell state", sort = [:Normal, :Stressed, :Dead]},
        column = {"sigma:n", title = title},
    )

    return p
end

function plot_antitoxin_experiment(
    cell_df,
    antitoxin_df;
    width = nothing,
    height = nothing,
    maxbins = 10,
)

    # Proportion of states in population.
    p1 = @vlplot(
        width = width,
        height = height,
        data = cell_df,
        mark = :line,
        x = {"time:q", title = "Exposure time"},
        y = {"value:q", title = "Fraction population"},
        color = {
            "variable:n",
            sort = [:Normal, :Stressed, :Dead],
            legend = {title = "Cell state", orient = "bottom"},
        },
        column = {"sigma:n", title = ""}
    )

    # Distribution of antitoxins.
    p2 = @vlplot(
        width = width,
        height = height,
        data = antitoxin_df,
        mark = :bar,
        x = {"antitoxin:q", bin = {maxbins = maxbins}, title = "Amount antitoxin (binned)"},
        y = {aggregate = "count", type = "quantitative"},
        color = {"antitoxin:q", bin = true, legend = false, scale = {scheme = "blues"}},
        column = {
            "sigma:n",
            title = string("μ = ", cell_df[!, :μ][1], " with σ values of:"),
            legend = false,
        },
        #resolve = {scale = {y = "independent"}}
    )

    return vcat(p2, p1)
end
