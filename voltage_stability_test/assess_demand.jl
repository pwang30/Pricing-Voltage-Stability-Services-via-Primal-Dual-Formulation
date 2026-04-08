# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Pricing voltage stability
# 4.July.2025

import Pkg
using JuMP,Gurobi, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles,Plots,MAT,Ipopt
include("dataset_gene.jl")
include("admittance_matrix_calculation.jl") 
include("offline_trainning.jl")
include("admittance_matrix.jl")
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

G_matrix, B_matrix = build_GB_matrices("Linespara.csv")  # build G and B matrices for AC power flow

model= Model()

Pˢᴳₘₐₓ=[120, 100, 75, 50, 45, 45]                               # ACtive max generation of SGs,  buses:2,3,4,5,27,30
Pˢᴳₘᵢₙ= Pˢᴳₘₐₓ*0.3                                              # ACtive min generation of SGs,  buses:2,3,4,5,27,30
Pⱽˢᴳₘₐₓ=[100]                                                   # ACtive max generation of VSGs,  buses:1
Pⱽˢᴳₘᵢₙ=[0]                                                     # ACtive min generation of VSGs,  buses:1


Qˢᴳₘₐₓ=   Pˢᴳₘₐₓ*0.5                                           # Reactive max generation of SGs,  buses:2,3,4,5,27,30
Qˢᴳₘᵢₙ=   -Pˢᴳₘₐₓ*0.5                                            # Reactive min generation of SGs,  buses:2,3,4,5,27,30
Qⱽˢᴳₘₐₓ=  Pⱽˢᴳₘₐₓ*0.5                                           # Reactive max generation of VSGs,  buses:1
Qⱽˢᴳₘᵢₙ=  -Pⱽˢᴳₘₐₓ*0.5                                           # Reactive min generation of VSGs,  buses:1

@variable(model, c[1:30,1:30])
@variable(model, s[1:30,1:30])
@variable(model,obj_P[1:30])
@variable(model,obj_Q[1:30])

for i in 1:bus_num
        @constraint(model, obj_P[i] ==  G_matrix[i,i]*c[i,i] + 
            sum( G_matrix[i,j] *c[i,j] - B_matrix[i,j] *s[i,j] for j in bus_connection[i] ) )   
        @constraint(model, obj_Q[i] == -B_matrix[i,i]*c[i,i] - 
            sum( B_matrix[i,j] *c[i,j] + G_matrix[i,j] *s[i,j] for j in bus_connection[i] ))    

        @constraint(model, c[data_bus[i,1],data_bus[i,1]] <= (data_bus[i,12])^2)    # c_{i,i} <=  V_{i,max}^2   eq. 50(c)
        @constraint(model, c[data_bus[i,1],data_bus[i,1]] >= (data_bus[i,13])^2)    # V_{i,min}^2  <= c_{i,i}   eq. 50(c)
        for j in bus_connection[i]
            @constraint(model, c[i,j] - c[j,i] == 0)                              # c_{i,j} == c_{j,i}   eq. 50(f)
            @constraint(model, s[i,j] + s[j,i] == 0)                              # s_{i,j} == -s_{j,i}   eq. 50(f)
            @constraint(model, c[i,j]^2 + s[i,j]^2 <= c[i,i] * c[j,j])        # c_{i,j}^2+s_{i,j}^2 <= c_{i,i}*c_{j,j}   eq. 50(g)       
        end    
end 


@objective(model, Max, sum(obj_Q)+ sum(obj_P))  # single-level objective function
#-------Solve and Output Results
set_optimizer(model , Gurobi.Optimizer)
optimize!(model)

obj_P=JuMP.value.(obj_P)
obj_Q=JuMP.value.(obj_Q)
plot(obj_P)
plot!(obj_Q)


