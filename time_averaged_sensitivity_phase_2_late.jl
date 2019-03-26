#This file is used to calculate the time averaged sensitivity matrix for late phase 2. time unit 331 to 350.

using GRNSimKit
using PyPlot

#initialize array to keep track of changes in P1,P2,P3 levels due to change in parameters
X_i_total = Array{}([])
# time step -
time_step_size = 1.0*(1/60)

#build a FOR loop to iterate over all the different changed parameter values. Goal is to get an array of
#P1, P2 and P3 for a time period of 20 mins for phase 1.
#Another goal is to plot the graphs vs parameter changes.


for i in 1:13
path_to_model_file = "$(pwd())/$i.json"
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

# make a plot to see how the parameters change the protein production -
plot(T*(60),X[:,4],color_P1,linewidth=2, linestyle="--")
plot(T*(60),X[:,5],color_P2,linewidth=2, linestyle="--")
plot(T*(60),X[:,6],color_P3,linewidth=2, linestyle="--")

# axis -
xlabel("Time (min)", fontsize=16)
ylabel("Protein (nmol/gDW)", fontsize=16)
pygui(true)

#create array to look at 20 minute changes of each protein for 1 parameter
X_i = X[331:350,4:6]
#array to store changes over all parameters
push!(X_i_total,X_i)
end

#Calculate P1, P2, P3 concentrations for original (unchanged) parameter values
#Everything the same as in part a except final X array is relabeled as X_original

using GRNSimKit

# time step -
time_step_size = 1.0*(1/60)

path_to_model_file = "$(pwd())/ICT1.json"

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
X_original = [X0 ; X1]


#Now, calculate the derivative dx/dp where x=P1, P2 or P3.
#dp array lists change in parameter values. Data in excel sheet.

dp = [
    0.048        ;   # transcription_saturation_constant
    90.928       ;   # translation_saturation_constant
    8.4          ;   # characteristic_initiation_time_transcription
    3            ;   # characteristic_initiation_time_translation
    0.132        ;   # cell_doubling_time
    5.6E-14      ;   # mass_of_single_cell
    0.14         ;   # fraction_of_water_per_cell
    230          ;   # copies_of_rnapII_per_cell
    9000         ;   # copies_of_ribosome_per_cell
    3.3          ;   # translation_elongation_rate
    12           ;   # transcription_elongation_rate
    200          ;   # characteristic_transcript_length
    66           ;   # characteristic_protein_length
];


#initialize array to keep track of derivatives
dProt_j_by_dp_i_time_331 = rand(13,3)
dProt_j_by_dp_i_time_341 = rand(13,3)
dProt_j_by_dp_i_time_350 = rand(13,3)

#For time value 331
for i in 1:13          # 13 parameter values
    for j in 1:3       # P1, P2, P3
        dProt_j_by_dp_i_time_331[i,j] = (X_i_total[i][1,j] - X_original[331,j+3])/(dp[i])
    end
end

#For time value 341
for i in 1:13          # 13 parameter values
    for j in 1:3       # P1, P2, P3
        dProt_j_by_dp_i_time_341[i,j] = (X_i_total[i][10,j] - X_original[341,j+3])/(dp[i])
    end
end

#For time value 350
for i in 1:13          # 13 parameter values
    for j in 1:3       # P1, P2, P3
        dProt_j_by_dp_i_time_350[i,j] = (X_i_total[i][20,j] - X_original[350,j+3])/(dp[i])
    end
end


#Now, to get the time averaged sensitivity coefficients, calculate mean of the three matrices above:
#Ideally, I would have used all 20 time values, but due to time constraints, I chose only 3 time values.

s_time_averaged_sensitivity_phase_2_late = (dProt_j_by_dp_i_time_331+dProt_j_by_dp_i_time_341+dProt_j_by_dp_i_time_350)/3
