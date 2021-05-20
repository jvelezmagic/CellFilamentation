using DrWatson
@quickactivate :CellFilamentation

using Random

cell = Cell(max_height = 4.0);
t_ramp, toxin_ramp = toxin_entry(10, 100, 100, 5, false);

Random.seed!(1234)
simulation_plot(cell, t_ramp, toxin_ramp, 1000)