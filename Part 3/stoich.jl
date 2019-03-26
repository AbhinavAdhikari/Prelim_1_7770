include("calculate_optimal_flux_distribution.jl")
include("tx_tl_rates.jl")

#define objective value array that stores v5 flux corresponding to each inducer concentration
v5_flux_array = similar(I)  # I defined in tx_tl_rates.jl 
for i in 1:length(I)
a=1; #stoichiometric matrix variable
n=1; #stoichiometric matrix variable
stoichiometric_matrix = Array{Float64}([
    #v1 to v6 and b1 to b9 as columns
   -1  1  0  0  0  0  0  0  0  0  0  0  0  0  0   ; # G
    1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0   ; # G*
   -1  1  0  0  0  0  0  0  0  0  0  0  0  0  0   ; # RNAP
    0  1 -1 -1  1  0  0  0  0  0  0  0  0  0  0   ; # mRNA
    0  0  0 -1  1  0  0  0  0  0  0  0  0  0  0   ; # rib
    0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0   ; # rib*
    0  0  0  0 -a  1  0  0  0  0  0  0  0  0  0   ; # AAtRNA
    0  0  0  0  a -1  0  0  0  0  0  0  0  0  0   ; # tRNA
    0  2n 0  0  2a 2  0  0  0  0  0  0  0  0 -1   ; # P_i
    0 -n  0  0  0  0  0  1  0  0  0  0  0  0  0   ; # NTP
    0  0  n  0  0  0  0  0  0 -1  0  0  0  0  0   ; # NMP
    0  0  0  0 -2a 0  0  0  0  0  0  0  1  0  0   ; # GTP
    0  0  0  0  2a 0  0  0  0  0  0  0  0 -1  0   ; # GDP
    0  0  0  0  1  0  0  0 -1  0  0  0  0  0  0   ; # protein
    0  0  0  0  0 -1  1  0  0  0  0  0  0  0  0   ; # AA
    0  0  0  0  0 -1  0  0  0  0  1  0  0  0  0   ; # ATP
    0  0  0  0  0  1  0  0  0  0  0 -1  0  0  0   ; # AMP
]);


# Setup default flux bounds array -
global(default_bounds_array)
default_bounds_array = Array{Float64}([

    0	100000	;	# v1
    r_x_hat[i]	r_x_hat[i]  ;	# v2
    0	100000	;	# v3
    0	100000	;	# v4
    0	r_l_hat[i]	    ;	# v5
    0	100000	;	# v6
    -100000	100000	;	# b1
    -100000	100000	;	# b2
    -100000	100000	;	# b3
    -100000	100000	;	# b4
    -100000	100000	;	# b5
    -100000	100000	;	# b6
    -100000	100000	;	# b7
    -100000	100000	;	# b8
    -100000	100000	;	# b9

]);

species_bounds_array = Array{Float64}([

     0	  0	;	# 1
     0	  0	;	# 2
     0	  0	;   # 3
     0	  0	;	# 4
     0	  0	;	# 5
     0	  0	;   # 6
     0	  0	;   # 7
     0	  0	;   # 8
     0	  0	;	# 9
     0	  0	;   # 10
     0	  0	;   # 11
     0	  0	;   # 12
     0	  0	;   # 13
     0	  0	;	# 14
     0	  0	;   # 15
     0	  0	;   # 16
     0	  0	;   # 17

]);

# Setup the objective coefficient array. Urea production = 1.0
objective_coefficient_array = Array{Float64}([

    0	;	# v1
    0	;	# v2
    0	;	# v3
    0	;	# v4
    1	;	# v5
    0   ;   # v6
    0	;	# b1
    0	;	# b2
    0	;	# b3
    0	;	# b4
    0	;	# b5
    0	;	# b6
    0	;	# b7
    0	;	# b8
    0	;	# b9

]);


#Run calculate_optimal_flux_distribution function
(objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(stoichiometric_matrix,default_bounds_array,species_bounds_array,objective_coefficient_array,min_flag= false)

global(v5_flux_array)
v5_flux_array[i]=objective_value
end

logarithm_I = similar(I)
for i in 1:length(I)
    logarithm_I[i] = log10(I[i])
end

using PyPlot
pygui(true)
plot(logarithm_I,v5_flux_array)
title("Protein levels vs Inducer concentration")
xlabel("Inducer concentration in mM")
ylabel("Protein levels in μM/hour")
