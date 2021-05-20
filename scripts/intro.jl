using DrWatson
@quickactivate :CellFilamentation

cell = Cell(max_height = 4.0);
t_ramp, toxin_ramp = toxin_entry(10, 100, 100, 5, false);

sim, out_cell = simulate_filamentation(cell, t_ramp, toxin_ramp, trajectories = 1000)
simm = EnsembleSummary(sim)

plot(simm, fillalpha = 0.2, fillcolor = "gray", linecolor = "black", layout = (2, 1))
