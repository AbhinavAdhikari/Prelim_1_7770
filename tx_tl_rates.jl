# This file is used to calculate the rate of transcription and translation. These rates will act as bounds for the fluxes in stoich.jl file.

#Given and calculated parameters
w1=0.26 #for control function
w2=300  #for control function
G_j=0.005 #gene concentration in microMolar.
n=1.5
K=0.30 #K for control function in microMolar
I=Array{Float64}([0.0001:0.0001:10...]) #inducer concentration in microMolar

#initiate arrays
F_I=similar(I) #F in control function equation
u_I=similar(I) #control function
r_x_hat=similar(I) #rate of transcription
r_l_hat=similar(I) #rate of translation
m_j=similar(I) #mRNA levels

#mRNA
k_e=(60/924)*3600 #mRNA elongation constant in 1/h. 60 nt/s divided by 924 nt gene length and converted to 1/h
R_xt=0.15 #RNAP concentration in microMolar
K_x=0.3 #Saturation constant transcription in microMolar
tau_x=2.7 #time constant transcription

#protein
k_e_p=(16.5/308)*3600 #protein elongation constant in 1/h. 16.5 aa/s divided by 308 aa peptide length and converted to 1/h
R_lt=1.6 #Ribosome concentration in microMolar
K_x_p=57.0 #Saturation constant translation in microMolar
tau_l=0.8 #time constant translation

#Write control function terms and rate of transcription
for i in 1:length(I)
    F_I[i] =(I[i]^n/(K^n+I[i]^n))
    u_I[i]=(w1+w2*F_I[i])/(1+w1+w2*F_I[i])
    r_x_hat[i]=(u_I[i]*(k_e*R_xt*G_j))/(tau_x*K_x + (tau_x + 1)*G_j)
    m_j[i]=(r_x_hat[i]/(8.35))   #from the equation dmRNA/dt = r_x_hat - (μ+k_d)*mRNA. At steady state, d/dt = 0
    r_l_hat[i]= (1*(k_e_p*R_lt*m_j[i]))/(tau_l*K_x_p + (tau_l + 1)*m_j[i])
end
