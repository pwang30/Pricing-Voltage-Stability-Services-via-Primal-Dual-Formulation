# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Pricing voltage stability
# 4.July.2025


#------------------------------------------------------------------------------
#--------Price voltage stability with restricted method------------------------
#------------------------------------------------------------------------------



import Pkg
using JuMP,Gurobi, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles,Plots,MAT,Optimization, OptimizationMetaheuristics
include("dataset_gene.jl")
include("admittance_matrix_calculation.jl") 
include("offline_trainning.jl")
include("calculate_Tao.jl")
include("calculate_delta_gamma.jl")
# SGs, buses:2,3,4,5,27,30    VSG, buses:1,     IBGs, buses:23,24

data_bus = CSV.read("IEEE30_Bus_Data.csv", DataFrame)
Bus_Pd = CSV.read("Bus_Pd_24h.csv", DataFrame)
Bus_Qd = CSV.read("Bus_Qd_24h.csv", DataFrame)

bus_connection = Vector{Vector{Int}}(undef, 30) # pre-allocate space for storing the connected buses with the given one
bus_connection[1] = [2, 3]
bus_connection[2] = [1, 4, 5, 6]
bus_connection[3] = [1, 4]
bus_connection[4] = [2, 3, 6, 12]
bus_connection[5] = [2, 7]
bus_connection[6] = [2, 4, 7, 8, 9, 10, 28]
bus_connection[7] = [5, 6]
bus_connection[8] = [6, 28]
bus_connection[9] = [6, 11, 10]
bus_connection[10] = [6, 9, 17, 20, 21, 22]
bus_connection[11] = [9]
bus_connection[12] = [4, 13, 14, 15, 16]
bus_connection[13] = [12]
bus_connection[14] = [12, 15]
bus_connection[15] = [12, 14, 18, 23]
bus_connection[16] = [12, 17]
bus_connection[17] = [10, 16]
bus_connection[18] = [15, 19]
bus_connection[19] = [18, 20]
bus_connection[20] = [10, 19]
bus_connection[21] = [10, 22]
bus_connection[22] = [10, 21, 24]
bus_connection[23] = [15, 24]
bus_connection[24] = [22, 23, 25]
bus_connection[25] = [24, 26, 27]
bus_connection[26] = [25]
bus_connection[27] = [25, 28, 29, 30]
bus_connection[28] = [6, 8, 27]
bus_connection[29] = [27, 30]
bus_connection[30] = [27, 29]



#-----------------------------------Representation and Approximation of VOLTAGE STABILITY Constraints  &  Visualization-----------------------------------
nᵥ= []              # The interval is evenly divided into nᵥ parts
for i in 0:0.001:1
    push!(nᵥ,i)
end 
IBG=[23,24]         # GFL.   The location (Bus) of IBGs.
Gᵥ=[1]              # FM.    The location (Bus) of Virtual Synchronous Generators (VSGs).
zᵟ_c, matrix_ω =dataset_gene(nᵥ,Gᵥ,IBG)                                                             # data set generation      +sum(K_c_m[2,:] .*matrix_ω[i,14:91])                
K_c_gc, K_c_gv, K_c_m, sh_1, sh_2 = offline_trainning(zᵟ_c, matrix_ω)  # offline_trainning



#-----------------------------------Define Parameters for Optimization-----------------------------------  
bus_num=30                              # number of buses
T=24                                    # number of periods

Pᴰ=Matrix(Bus_Pd[:,(2:25)]  )                  # system demand of active power at each bus
Qᴰ=Matrix(Bus_Qd[:,(2:25)]  )               # system demand of reactive power at each bus

Pᴰ[23,:]=Pᴰ[23,:] + Pᴰ[5,:].*2/5
Pᴰ[24,:]=Pᴰ[24,:] + Pᴰ[5,:].*2/5
Pᴰ[5,:]=Pᴰ[5,:]*1/5

Qᴰ[23,:]=Qᴰ[23,:] + Qᴰ[5,:]*2/5
Qᴰ[24,:]=Qᴰ[24,:] + Qᴰ[5,:]*2/5
Qᴰ[5,:]=Qᴰ[5,:]*1/5

S_lines=ones(30,30)*70                  # line capacity for apparent power

Pˢᴳₘₐₓ= [55, 55, 55, 55, 55, 55]
Pˢᴳₘᵢₙ= Pˢᴳₘₐₓ*0.3                      # ACtive min generation of SGs,  buses:2,3,4,5,27,30
Qˢᴳₘₐₓ= Pˢᴳₘₐₓ*0.6                          # REACtive max generation of SGs,  buses:2,3,4,5,27,30
Qˢᴳₘᵢₙ= -Pˢᴳₘₐₓ*0.6                      # REACtive min generation of SGs,  buses:2,3,4,5,27,30

Pⱽˢᴳₘₐₓ=[100]                           # ACtive max generation of VSGs, buses:1
Pⱽˢᴳₘᵢₙ=[0]                             # ACtive min generation of VSGs,  buses:1
Qⱽˢᴳₘₐₓ= Pⱽˢᴳₘₐₓ *0.6                         # REACtive max generation of SGs,  buses:2,3,4,5,27,30
Qⱽˢᴳₘᵢₙ= -Pⱽˢᴳₘₐₓ*0.6  

Pᴵᴮᴳₘₐₓ= [100, 100]                      # ACtive max generation of IBGs, buses:23,24
Pᴵᴮᴳₘᵢₙ=  [0, 0]                        # ACtive min generation of IBGs, buses:23,24

Qᴵᴮᴳₘₐₓ= Pᴵᴮᴳₘₐₓ*0.6
Qᴵᴮᴳₘᵢₙ= -Pᴵᴮᴳₘₐₓ*0.6


Sᵐᵃˣ_gc=Pˢᴳₘₐₓ*1.2
Sᵐᵃˣ_gv=Pⱽˢᴳₘₐₓ*1.2
Sᵐᵃˣ_c=Pᴵᴮᴳₘₐₓ*1.2

α_VSG=0.9
α_IBG=1
ratio=1

Kˢᵗ=[125, 125, 125, 125, 125, 125]              #  Startup cost of SGs                     SGs, buses:2,3,4,5,27,30
Kˢʰ=[50, 50, 50, 50, 50, 50]            #  Shutdown cost of SGs                    SGs, buses:2,3,4,5,27,30
Oᵐ=[10.47, 10.47, 10.47, 10.47, 10.47, 10.47]                   #  Marginal generation cost of SGs   in    SGs, buses:2,3,4,5,27,30
Oⁿˡ=Oᵐ*10 
