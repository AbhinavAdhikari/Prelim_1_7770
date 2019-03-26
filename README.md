# Prelim_1_7770: Abhinav Adhikari

## Part 1
`Part 1.pdf` file contains the writeup
## Part 2:

### Part 2a

2a.jl: This file uses the `ICT1.json` file provided by Dr. Varner. Output is a plot that simulates the response of P1, P2 and P3 versus change in initiator concentration.
	
### Part 2b
	
There are several files pertinent to this part of the problem. 

*.json files (* = 1 to 13) are edited json files. In these files, each parameter is updated to a value 20% larger than its original value. For example, in 1.json, parameter 1 is increased by 20% with all other parameters held constant. 

The list of parameters is in file `percentage_change.pdf`. 

The `ICIT1.json` file holds the original parameters and is also required.

Run files: `time_averaged_sensitivity_phase_1.jl`, `time_averaged_sensitivity_phase_2.jl`, `time_averaged_sensitivity_phase_2_late.jl` solve for the time-averaged sensitivity matrices for each of the phases. 

`time_averaged_sensitivity_phase_1.jl`: This file is used to calculate the time averaged sensitivity matrix for phase 1, i.e. from time unit 21 to 40- an interval of 20 mins.

`time_averaged_sensitivity_phase_2.jl`: This file is used to calculate the time averaged sensitivity matrix for early phase 1, i.e. from time unit 71 to 90- an interval of 20 mins.

`time_averaged_sensitivity_phase_2_late.jl`: This file is used to calculate the time averaged sensitivity matrix for early phase 1, i.e. from time unit 331 to 350- an interval of 20 mins.

### Part 2c

`SVD.jl`: This file performs SVD on the 3 time averaged sensitivity arrays calculated in the 3 run files above in part 2b. The rankings and explanation are also given in this file.

## Part 3:

### Part 3a

`stoich.jl`: This file has the stoichiometric matrix according to table 1 of Allen and Palsson. It is also the main run file. It solves for the FBA problem and plots the protein levels versus different inducer concentrations.

`tx_tl_rates.jl`: This file is used to calculate the rate of transcription and translation. These rates will act as bounds for the fluxes.

### Part 3b

`stoich.jl`: See above description

`calculate_optimal_flux_distribution`: required function file for FBA problem.

### Part 3c

`change_b*.jl`: this file changes the bounds individually for each of the * values (* = 1 to 9) to be a small value [-0.1, 0.1] keeping other b* constant.

`change exchange bounds.pdf`: contains graphs of how changing each bound affects the protein levels. Also has the explanation for the trends of the graphs.
	

