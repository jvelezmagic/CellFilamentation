export Cell
export height, surface_area, volume, surface_area_over_volume
export geometric_properties

@with_kw mutable struct Cell{T<:AbstractFloat}
    @deftype T
    side_length = 0.6 #@assert side_length > 0.0
    radius = 0.2
    @assert radius > 0.0
    internal_toxin = 0.0
    @assert internal_toxin ≥ 0.0
    toxin_diffusion_rate = 0.3
    antitoxin = 0.0
    antitoxin_efficacy = 0.0
    filamentation_rate = 0.1
    filamentation_threshold = 0.5
    kill_threshold = 1.0
    max_height = 2 * (side_length + 2 * radius)
    noise_int_toxin = 0.10
    noise_height = 0.0
    τ_filamentation = Inf
    τ_kill = Inf
    τ_delay = 0.0
    state::Symbol = :Normal
    @assert state ∈ (:Normal, :Stressed, :Dead)
end

height(cell::Cell) = cell.side_length + 2 * cell.radius
surface_area(cell::Cell) = 2π * cell.radius * height(cell)
volume(cell::Cell) = π * cell.radius^2 * (4 / 3 * cell.radius + cell.side_length)
surface_area_over_volume(cell::Cell) = surface_area(cell) / volume(cell)

function geometric_properties(cell::Cell)
    map(x -> x(cell), [height, surface_area, volume, surface_area_over_volume])
end
