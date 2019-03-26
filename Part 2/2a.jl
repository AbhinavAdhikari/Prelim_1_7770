#Part 2a. This file uses the ICT1.json file provided by Dr. Varner. 
#Output is a plot that simulates the response of P1, P2 and P3 versus
# change in initiator concentration.

using GRNSimKit
using PyPlot

# time step -
time_step_size = 1.0*(1/60)   

path_to_model_file = "$(pwd())/ICT1.json"
color_P1 = "red"
color_P2 = "green"
color_P3 = "blue"

# Build a data dictionary from a model file -
ddd = build_discrete_dynamic_data_dictionary_from_model_file(time_step_size, path_to_model_file)

# Run the model to steady-state, before we do anything -
steady_state = GRNSteadyStateSolve(ddd)

# Run the model for 60 mins *before* we add inducer -
ddd[:initial_condition_array] = steady_state
(T0, X0) = GRNDiscreteDynamicSolve((0.0,1.0,time_step_size), ddd)

# Add inducer for 300 min after the extra 60 min run-
ddd[:initial_condition_array] = X0[end,:]
ddd[:initial_condition_array][7] = 10.0
tstart_1 = T0[end]
tstop_1 = tstart_1 + 5.0
(T1, X1) = GRNDiscreteDynamicSolve((tstart_1,tstop_1, time_step_size), ddd)

# Package -
T = [T0 ; T1]
X = [X0 ; X1]

# make a plot -
plot(T*(60),X[:,4],color_P1,linewidth=2, linestyle="--")
plot(T*(60),X[:,5],color_P2,linewidth=2, linestyle="--")
plot(T*(60),X[:,6],color_P3,linewidth=2, linestyle="--")

# axis -
xlabel("Time (min)", fontsize=16)
ylabel("Protein (nmol/gDW)", fontsize=16)
pygui(true)
