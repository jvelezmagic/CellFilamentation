export filamentation_callbacks!

function filamentation_callbacks!()

    function conditions(out, u, t, integrator)
        out[1] = u[1] - integrator.p[1].filamentation_threshold
        out[2] = u[1] - integrator.p[1].kill_threshold
        out[3] = u[2] - integrator.p[1].max_height
    end

    function affect!(integrator, idx)
        if idx == 1
            integrator.p[1].state = :Stressed
            if isinf(integrator.p[1].τ_filamentation)
                integrator.p[1].τ_filamentation = integrator.t
            end
        elseif idx == 2
            integrator.p[1].state = :Dead
            integrator.p[1].τ_kill = integrator.t
            terminate!(integrator)
        elseif idx == 3
            integrator.u[2] = integrator.p[1].max_height
        end
    end

    # Apply callbacks only when condition is found to be zero and upcrossing.
    cb = VectorContinuousCallback(conditions, affect!, 3; affect_neg! = nothing)

    return cb
end
