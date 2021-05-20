using DrWatson
@quickactivate :CellFilamentation

cell = Cell(max_height = 4.0);
t_ramp, toxin_ramp = toxin_entry(10, 100, 100, 5, false);

simulation, out_cell = simulate_filamentation(cell, t_ramp, toxin_ramp, trajectories = 100)
plot(simulation)