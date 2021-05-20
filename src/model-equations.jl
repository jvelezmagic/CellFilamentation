export filamentation_de_sym
export filamentation_de
export filementation_noise

function filamentation_de(du, u, p, t)

    internal_toxin, height = u

    cell, interpolate_external_toxin = p

    external_toxin = interpolate_external_toxin(t)

    # Compare geometric properties of the initial cell
    # with respect to the current one.
    c_cell = Cell(cell; side_length = height - 2 * cell.radius)
    surfa_area_ratio = surface_area(cell) / surface_area(c_cell)
    volume_ratio = volume(cell) / volume(c_cell)

    toxin_volume = internal_toxin / volume_ratio
    toxin_surface_area = cell.toxin_diffusion_rate * surfa_area_ratio

    # Growth cell if :Stressed.
    dheight =
        cell.state == :Stressed &&
        t ≥ (cell.τ_filamentation + cell.τ_delay) &&
        height < cell.max_height ? cell.filamentation_rate * height : 0.0

    # Calculate internal toxin change.
    dinternal_toxin = (
        toxin_surface_area * (external_toxin - toxin_volume) -
        cell.antitoxin_efficacy * cell.antitoxin * internal_toxin
    )

    du .= dinternal_toxin, dheight
end

function filamentation_noise(du, u, p, t)
    du[1] = u[1] * p[1].noise_int_toxin
    du[2] = u[2] * p[1].noise_height
end

filamentation_de_sym = ODEFunction(filamentation_de, syms=[:internal_toxin, :length])