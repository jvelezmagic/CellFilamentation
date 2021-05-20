export simulate
export simulate_filamentation

function simulate(
    cell::Cell,
    time::Array{Float64,1},
    toxin::Array{Float64,1};
    stochastic::Bool = false,
    status::Bool = false,
    kwargs...,
)
    # This helps not to modify the original input cell.
    cell = Cell(cell)

    # Set initial conditions, time-span and parameters needed for the model.
    u0 = @> [cell.internal_toxin, height(cell)] reshape((2, 1))
    tspan = (time[1], time[end])
    p = [cell, extrapolate(interpolate((time,), toxin, Gridded(Linear())), Flat())]

    # Definine the problem.
    cb = filamentation_callbacks!()
    prob =
        stochastic ?
        SDEProblem(filamentation_de, filamentation_noise, u0, tspan, p; callback = cb) :
        ODEProblem(filamentation_de, u0, tspan, p; callback = cb)

    # Solve the problem.
    sol = solve(prob; isoutofdomain = (y, p, t) -> any(x -> x < 0, y), kwargs...)

    if status
        return sol, cell, _check_cell_status(cell, time)
    else
        return sol, cell
    end
end

function simulate_filamentation(
    cell::Cell,
    time::Array{Float64,1},
    toxin::Array{Float64,1};
    trajectories::Int64 = 0,
    kwargs...,
)
    # This helps not to modify the original input cell.
    cell = Cell(cell)

    # Set initial conditions, time-span and parameters needed for the model.
    u0 = @> [cell.internal_toxin, height(cell)] reshape((2, 1))
    tspan = (time[1], time[end])
    p = [cell, extrapolate(interpolate((time,), toxin, Gridded(Linear())), Flat())]

    # Set parameters to control the global behaivour of resolution of the problem.
    cb = filamentation_callbacks!()
    isoutofdomain = (y, p, t) -> any(x -> x < 0, y)

    # Define and solve the problem; SDE or ODE.
    if trajectories > 0
        sol = @> begin
            SDEProblem(filamentation_de, filamentation_noise, u0, tspan, p; callback = cb)
            EnsembleProblem()
            solve(
                EnsembleThreads(),
                trajectories = trajectories,
                isoutofdomain = isoutofdomain,
                kwargs...,
            )
        end
    else
        sol = @> begin
            ODEProblem(filamentation_de, u0, tspan, p; callback = cb)
            solve(isoutofdomain = isoutofdomain, kwargs...)
        end
    end

    return sol, cell
end

function _check_cell_status(cell::Cell, τs)
    # Correct operation is guaranteed only within the evaluation range of the differential equation.
    @unpack τ_filamentation, τ_kill = cell
    return map(x -> x < τ_filamentation ? :Normal : x < τ_kill ? :Stressed : :Dead, τs)
end
