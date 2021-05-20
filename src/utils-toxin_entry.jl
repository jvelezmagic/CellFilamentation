export toxin_entry

function toxin_entry(start, stop, len, max_toxin = 1.0, stable = true)

    time = range(0.0, stop = stop, length = len) |> collect
    toxin = zeros(len)

    entry = @. start <= time <= stop

    if stable
        toxin[entry] .= 1 * max_toxin
    else
        toxin[entry] = range(0.0, stop = max_toxin, length = sum(entry)) |> collect
    end

    return (time, toxin)
end
