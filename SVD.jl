include("time_averaged_sensitivity_phase_1.jl")
include("time_averaged_sensitivity_phase_2.jl")
include("time_averaged_sensitivity_phase_2_late.jl")

using LinearAlgebra

#Perform SVD on the 3 time averaged sensitivity arrays

svd_1 = svd(s_time_averaged_sensitivity_phase_1)
svd_2 = svd(s_time_averaged_sensitivity_phase_2_early)
svd_3 = svd(s_time_averaged_sensitivity_phase_2_late)

svd_1.U
# We see that parameter 6 has the largest magnitude followed by p7, p5, p1.
svd_2.U
# We see that parameter 6 has the largest magnitude followed by p7, p1, p4.
svd_3.U
# We see that parameter 2 has the largest magnitude followed by p3, p8, p6.

#Overall results:
# 1st phase rankings: mass_of_single_cell, fraction_of_water_per_cell, 
#                   cell_doubling_time, transcription_saturation_constant

# Early 2nd phase rankings: mass_of_single_cell, fraction_of_water_per_cell,
#               transcription_saturation_constant, characteristic_initiation_time_translation

# Late 2nd phase rankings: translation_saturation_constant, characteristic_initiation_time_transcription
#               copies_of_rnapII_per_cell, mass_of_single_cell
