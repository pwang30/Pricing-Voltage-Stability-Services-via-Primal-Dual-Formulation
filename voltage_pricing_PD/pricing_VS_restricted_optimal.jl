# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Pricing voltage stability
# 4.July.2025


#------------------------------------------------------------------------------
#--------Price voltage stability with restricted method------------------------
#------------------------------------------------------------------------------



import Pkg
using JuMP,Gurobi, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles, Plots, DataStructures, Dualization
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
K_c_gc, K_c_gv, K_c_m, MAPE_z_23_1, MAPE_z_24_1 = offline_trainning(zᵟ_c, matrix_ω)  # offline_trainning



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
S_Base=50

α_VSG=0.9
α_IBG=1
ratio=1

Kˢᵗ=[125, 125, 125, 125, 125, 125]              #  Startup cost of SGs                     SGs, buses:2,3,4,5,27,30
Kˢʰ=[50, 50, 50, 50, 50, 50]            #  Shutdown cost of SGs                    SGs, buses:2,3,4,5,27,30
Oᵐ=[10.47, 10.47, 10.47, 10.47, 10.47, 10.47]                   #  Marginal generation cost of SGs   in    SGs, buses:2,3,4,5,27,30
Oⁿˡ=Oᵐ*10 

yˢᴳ²_0 = 1
yˢᴳ³_0 = 1
yˢᴳ⁴_0 = 0
yˢᴳ⁵_0 = 0
yˢᴳ²⁷_0 = 0
yˢᴳ³⁰_0 = 0


#-------------------------------------------------------
#-------------------Define model-------------------
#-------------------------------------------------------

#=
model= Model()                     

@variable(model, yˢᴳ²[1:T] >= 0)                         # status of SGs, buses:2,3,4,5,27,30.     ,Bin
@variable(model, yˢᴳ³[1:T] >= 0)            
@variable(model, yˢᴳ⁴[1:T] >= 0)           
@variable(model, yˢᴳ⁵[1:T] >= 0)           
@variable(model, yˢᴳ²⁷[1:T] >= 0)            
@variable(model, yˢᴳ³⁰[1:T] >= 0)  

@constraint(model, yˢᴳ² .<= 1) 
@constraint(model, yˢᴳ³ .<= 1) 
@constraint(model, yˢᴳ⁴ .<= 1) 
@constraint(model, yˢᴳ⁵ .<= 1) 
@constraint(model, yˢᴳ²⁷ .<= 1) 
@constraint(model, yˢᴳ³⁰ .<= 1) 

@variable(model, Cᵁ²[1:T]>=0)                           # startup costs for SGs                                          
@variable(model, Cᵁ³[1:T]>=0)                            
@variable(model, Cᵁ⁴[1:T]>=0)                                
@variable(model, Cᵁ⁵[1:T]>=0)                                
@variable(model, Cᵁ²⁷[1:T]>=0)                               
@variable(model, Cᵁ³⁰[1:T]>=0)                 

@variable(model, Pᴳ[1:30,1:T]) 
@variable(model, Qᴳ[1:30,1:T])  

for t in 2:T                                                # startup costs and shutdown costs occured at hour t
    @constraint(model, Cᵁ²[t]>=(yˢᴳ²[t]-yˢᴳ²[t-1])*Kˢᵗ[1])         
    @constraint(model, Cᵁ³[t]>=(yˢᴳ³[t]-yˢᴳ³[t-1])*Kˢᵗ[2]) 
    @constraint(model, Cᵁ⁴[t]>=(yˢᴳ⁴[t]-yˢᴳ⁴[t-1])*Kˢᵗ[3])
    @constraint(model, Cᵁ⁵[t]>=(yˢᴳ⁵[t]-yˢᴳ⁵[t-1])*Kˢᵗ[4])
    @constraint(model, Cᵁ²⁷[t]>=(yˢᴳ²⁷[t]-yˢᴳ²⁷[t-1])*Kˢᵗ[5])
    @constraint(model, Cᵁ³⁰[t]>=(yˢᴳ³⁰[t]-yˢᴳ³⁰[t-1])*Kˢᵗ[6])
end   



#--------------------------AC power flow constraints 

for t in 1:T                                            # bounds for the active and reactive output of each bus
    @constraint(model, Pᴳ[2,t] <= yˢᴳ²[t]*Pˢᴳₘₐₓ[1])    # bounds for active power of SGs: buses 2 3 4 5 27 30           
    @constraint(model, yˢᴳ²[t]*Pˢᴳₘᵢₙ[1] <= Pᴳ[2,t])       
    @constraint(model, Pᴳ[3,t] <= yˢᴳ³[t]*Pˢᴳₘₐₓ[2])       
    @constraint(model, yˢᴳ³[t]*Pˢᴳₘᵢₙ[2] <= Pᴳ[3,t])         
    @constraint(model, Pᴳ[4,t] <= yˢᴳ⁴[t]*Pˢᴳₘₐₓ[3])       
    @constraint(model, yˢᴳ⁴[t]*Pˢᴳₘᵢₙ[3] <= Pᴳ[4,t])    
    @constraint(model, Pᴳ[5,t] <= yˢᴳ⁵[t]*Pˢᴳₘₐₓ[4])       
    @constraint(model, yˢᴳ⁵[t]*Pˢᴳₘᵢₙ[4] <= Pᴳ[5,t])
    @constraint(model, Pᴳ[27,t] <= yˢᴳ²⁷[t]*Pˢᴳₘₐₓ[5])       
    @constraint(model, yˢᴳ²⁷[t]*Pˢᴳₘᵢₙ[5] <=Pᴳ[27,t])
    @constraint(model, Pᴳ[30,t] <= yˢᴳ³⁰[t]*Pˢᴳₘₐₓ[6])       
    @constraint(model, yˢᴳ³⁰[t]*Pˢᴳₘᵢₙ[6] <= Pᴳ[30,t])

    @constraint(model, Pᴳ[1,t] <= α_VSG*Pⱽˢᴳₘₐₓ[1])            # bounds for active power of VSGs and IBGs: buses 1 23 24           
    @constraint(model, Pⱽˢᴳₘᵢₙ[1] <= Pᴳ[1,t])       
    @constraint(model, Pᴳ[23,t] <= α_IBG*Pᴵᴮᴳₘₐₓ[1])       
    @constraint(model, Pᴵᴮᴳₘᵢₙ[1] <= Pᴳ[23,t])         
    @constraint(model, Pᴳ[24,t] <= α_IBG*Pᴵᴮᴳₘₐₓ[2])       
    @constraint(model, Pᴵᴮᴳₘᵢₙ[2] <= Pᴳ[24,t])   

    @constraint(model, Qᴳ[2,t] <= yˢᴳ²[t]*Qˢᴳₘₐₓ[1])    # bounds for REactive power of SGs: buses 2 3 4 5 27 30           
    @constraint(model, yˢᴳ²[t]*Qˢᴳₘᵢₙ[1] <= Qᴳ[2,t])       
    @constraint(model, Qᴳ[3,t] <= yˢᴳ³[t]*Qˢᴳₘₐₓ[2])       
    @constraint(model, yˢᴳ³[t]*Qˢᴳₘᵢₙ[2] <= Qᴳ[3,t])         
    @constraint(model, Qᴳ[4,t] <= yˢᴳ⁴[t]*Qˢᴳₘₐₓ[3])       
    @constraint(model, yˢᴳ⁴[t]*Qˢᴳₘᵢₙ[3] <= Qᴳ[4,t])    
    @constraint(model, Qᴳ[5,t] <= yˢᴳ⁵[t]*Qˢᴳₘₐₓ[4])       
    @constraint(model, yˢᴳ⁵[t]*Qˢᴳₘᵢₙ[4] <= Qᴳ[5,t])
    @constraint(model, Qᴳ[27,t] <= yˢᴳ²⁷[t]*Qˢᴳₘₐₓ[5])       
    @constraint(model, yˢᴳ²⁷[t]*Qˢᴳₘᵢₙ[5] <=Qᴳ[27,t])
    @constraint(model, Qᴳ[30,t] <= yˢᴳ³⁰[t]*Qˢᴳₘₐₓ[6])       
    @constraint(model, yˢᴳ³⁰[t]*Qˢᴳₘᵢₙ[6] <= Qᴳ[30,t])           

    @constraint(model, Qᴳ[1,t] <= Qⱽˢᴳₘₐₓ[1])       
    @constraint(model, Qⱽˢᴳₘᵢₙ[1] <= Qᴳ[1,t])

    @constraint(model, Qᴳ[23,t] <= ratio*Qᴵᴮᴳₘₐₓ[1])       
    @constraint(model, ratio*Qᴵᴮᴳₘᵢₙ[1] <= Qᴳ[23,t])
    @constraint(model, Qᴳ[24,t] <= ratio*Qᴵᴮᴳₘₐₓ[2])       
    @constraint(model, ratio*Qᴵᴮᴳₘᵢₙ[2] <= Qᴳ[24,t])

    @constraint(model, [ 1, Pᴳ[2,t] /(Sᵐᵃˣ_gc[1]), Qᴳ[2,t] /(Sᵐᵃˣ_gc[1]) ] in SecondOrderCone()  ) 
    @constraint(model, [ 1, Pᴳ[3,t] /(Sᵐᵃˣ_gc[2]), Qᴳ[3,t] /(Sᵐᵃˣ_gc[2]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[4,t] /(Sᵐᵃˣ_gc[3]), Qᴳ[4,t] /(Sᵐᵃˣ_gc[3]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[5,t] /(Sᵐᵃˣ_gc[4]), Qᴳ[5,t] /(Sᵐᵃˣ_gc[4]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[27,t] /(Sᵐᵃˣ_gc[5]), Qᴳ[27,t] /(Sᵐᵃˣ_gc[5]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[30,t] /(Sᵐᵃˣ_gc[6]), Qᴳ[30,t] /(Sᵐᵃˣ_gc[6]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[1,t] /(Sᵐᵃˣ_gv[1]), Qᴳ[1,t] /(Sᵐᵃˣ_gv[1]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[23,t] /(Sᵐᵃˣ_c[1]), Qᴳ[23,t] /(Sᵐᵃˣ_c[1]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[24,t] /(Sᵐᵃˣ_c[2]), Qᴳ[24,t] /(Sᵐᵃˣ_c[2]) ] in SecondOrderCone()  )

    gens_SGs = [2, 3, 4, 5, 27, 30]    # buses with SGs
    gens_VSGs = [1]                    # buses with VSGs
    gens_IBGs = [23, 24]               # buses with IBGs gens_all = union(gens_SGs, gens_VSGs, gens_IBGs)

    for i in 1:bus_num                 # no power is injected into the non-generation bus 
        if !(i in gens_SGs) && !(i in gens_VSGs) && !(i in gens_IBGs)
            @constraint(model, Pᴳ[i, t] == 0)
            @constraint(model, Qᴳ[i, t] == 0)
        end
    end
end

@variable(model, Pˡⁱⁿᵉ[1:30,1:30,1:T])
@variable(model, Qˡⁱⁿᵉ[1:30,1:30,1:T])



for i in 1:bus_num
    for t in 1:T

        @constraint(model, Pᴳ[i,t] - Pᴰ[i,t]*1.1 ==  sum( Pˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )   
        @constraint(model, Qᴳ[i,t] - Qᴰ[i,t] ==  sum( Qˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )   

        for j in bus_connection[i]
            @constraint(model, Pˡⁱⁿᵉ[i,j,t] + Pˡⁱⁿᵉ[j,i,t] == 0)
            @constraint(model, Qˡⁱⁿᵉ[i,j,t] + Qˡⁱⁿᵉ[j,i,t] == 0)
        end    
    end
end 



#--------------------------Voltage stability constraints--------------------------

@variable(model, Γ_23[1:T])
@variable(model, Γ_24[1:T])

@variable(model, ηₘ_1[1:15,1:T] >= 0)        # product of each pair of SG 
@constraint(model, ηₘ_1 .<= 1)
@variable(model, ηₘ_2[1:6,1:T])             # product of each pair of SG and VSG 


for t in 1:T
    # 15+6
    @constraint(model, ηₘ_1[1,t]>=yˢᴳ²[t]+yˢᴳ³[t]-1)  
    @constraint(model, ηₘ_1[1,t]<=yˢᴳ²[t])
    @constraint(model, ηₘ_1[1,t]<=yˢᴳ³[t])
    @constraint(model, ηₘ_1[2,t]>=yˢᴳ²[t]+yˢᴳ⁴[t]-1)  
    @constraint(model, ηₘ_1[2,t]<=yˢᴳ²[t])
    @constraint(model, ηₘ_1[2,t]<=yˢᴳ⁴[t])
    @constraint(model, ηₘ_1[3,t]>=yˢᴳ²[t]+yˢᴳ⁵[t]-1)  
    @constraint(model, ηₘ_1[3,t]<=yˢᴳ²[t])
    @constraint(model, ηₘ_1[3,t]<=yˢᴳ⁵[t])
    @constraint(model, ηₘ_1[4,t]>=yˢᴳ²[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model, ηₘ_1[4,t]<=yˢᴳ²[t])
    @constraint(model, ηₘ_1[4,t]<=yˢᴳ²⁷[t])
    @constraint(model, ηₘ_1[5,t]>=yˢᴳ²[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model, ηₘ_1[5,t]<=yˢᴳ²[t])
    @constraint(model, ηₘ_1[5,t]<=yˢᴳ³⁰[t])

    @constraint(model, ηₘ_1[6,t]>=yˢᴳ³[t]+yˢᴳ⁴[t]-1)  
    @constraint(model, ηₘ_1[6,t]<=yˢᴳ³[t])
    @constraint(model, ηₘ_1[6,t]<=yˢᴳ⁴[t])
    @constraint(model, ηₘ_1[7,t]>=yˢᴳ³[t]+yˢᴳ⁵[t]-1)  
    @constraint(model, ηₘ_1[7,t]<=yˢᴳ³[t])
    @constraint(model, ηₘ_1[7,t]<=yˢᴳ⁵[t])
    @constraint(model, ηₘ_1[8,t]>=yˢᴳ³[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model, ηₘ_1[8,t]<=yˢᴳ³[t])
    @constraint(model, ηₘ_1[8,t]<=yˢᴳ²⁷[t])
    @constraint(model, ηₘ_1[9,t]>=yˢᴳ³[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model, ηₘ_1[9,t]<=yˢᴳ³[t])
    @constraint(model, ηₘ_1[9,t]<=yˢᴳ³⁰[t])

    @constraint(model, ηₘ_1[10,t]>=yˢᴳ⁴[t]+yˢᴳ⁵[t]-1)  
    @constraint(model, ηₘ_1[10,t]<=yˢᴳ⁴[t])
    @constraint(model, ηₘ_1[10,t]<=yˢᴳ⁵[t])
    @constraint(model, ηₘ_1[11,t]>=yˢᴳ⁴[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model, ηₘ_1[11,t]<=yˢᴳ⁴[t])
    @constraint(model, ηₘ_1[11,t]<=yˢᴳ²⁷[t])
    @constraint(model, ηₘ_1[12,t]>=yˢᴳ⁴[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model, ηₘ_1[12,t]<=yˢᴳ⁴[t])
    @constraint(model, ηₘ_1[12,t]<=yˢᴳ³⁰[t])

    @constraint(model, ηₘ_1[13,t]>=yˢᴳ⁵[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model, ηₘ_1[13,t]<=yˢᴳ⁵[t])
    @constraint(model, ηₘ_1[13,t]<=yˢᴳ²⁷[t])
    @constraint(model, ηₘ_1[14,t]>=yˢᴳ⁵[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model, ηₘ_1[14,t]<=yˢᴳ⁵[t])
    @constraint(model, ηₘ_1[14,t]<=yˢᴳ³⁰[t])

    @constraint(model, ηₘ_1[15,t]>=yˢᴳ²⁷[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model, ηₘ_1[15,t]<=yˢᴳ²⁷[t])
    @constraint(model, ηₘ_1[15,t]<=yˢᴳ³⁰[t])

    @constraint(model, ηₘ_2[1,t]==yˢᴳ²[t]*α_VSG)
    @constraint(model, ηₘ_2[2,t]==yˢᴳ³[t]*α_VSG)
    @constraint(model, ηₘ_2[3,t]==yˢᴳ⁴[t]*α_VSG)
    @constraint(model, ηₘ_2[4,t]==yˢᴳ⁵[t]*α_VSG)
    @constraint(model, ηₘ_2[5,t]==yˢᴳ²⁷[t]*α_VSG)
    @constraint(model, ηₘ_2[6,t]==yˢᴳ³⁰[t]*α_VSG)

    @constraint(model, Γ_23[t]== 1/2 * ( ( K_c_gc[1,1] *yˢᴳ²[t] +K_c_gc[1,2] *yˢᴳ³[t] +K_c_gc[1,3] *yˢᴳ⁴[t] 
    +K_c_gc[1,4] *yˢᴳ⁵[t] +K_c_gc[1,5] *yˢᴳ²⁷[t] +K_c_gc[1,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[1,1]
    + sum(K_c_m[1,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[1,16:end] .*ηₘ_2[:,t]   ) ) )

    @constraint(model, Γ_24[t]== 1/2 * ( ( K_c_gc[2,1] *yˢᴳ²[t] +K_c_gc[2,2] *yˢᴳ³[t] +K_c_gc[2,3] *yˢᴳ⁴[t] 
    +K_c_gc[2,4] *yˢᴳ⁵[t] +K_c_gc[2,5] *yˢᴳ²⁷[t] +K_c_gc[2,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[2,1]
    + sum(K_c_m[2,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[2,16:end] .*ηₘ_2[:,t]   ) ) )

    @constraint(model, [Qᴳ[23,t] + Γ_23[t]*S_Base, Pᴳ[23,t], Qᴳ[23,t]] in SecondOrderCone()  )      
    @constraint(model, [Qᴳ[24,t] + Γ_24[t]*S_Base, Pᴳ[24,t], Qᴳ[24,t]] in SecondOrderCone()  )

end


#--------------------------------
#  Model solving
#--------------------------------

cost_onoff_Primal=sum(Cᵁ²)+sum(Cᵁ³)+sum(Cᵁ⁴)+sum(Cᵁ⁵)+sum(Cᵁ²⁷)+sum(Cᵁ³⁰)
cost_nl_Primal=sum(Oⁿˡ[1].*(yˢᴳ²))+sum(Oⁿˡ[2].*(yˢᴳ³))+sum(Oⁿˡ[3].*(yˢᴳ⁴))+sum(Oⁿˡ[4].*(yˢᴳ⁵))+sum(Oⁿˡ[5].*(yˢᴳ²⁷))+sum(Oⁿˡ[6].*(yˢᴳ³⁰))    
cost_gene_Primal=sum(Oᵐ[1].*Pᴳ[2,:])+ sum(Oᵐ[2].*Pᴳ[3,:])+sum(Oᵐ[3].*Pᴳ[4,:])+sum(Oᵐ[4].*Pᴳ[5,:])+sum(Oᵐ[5].*Pᴳ[27,:])+sum(Oᵐ[6].*Pᴳ[30,:])   
obj=cost_onoff_Primal +cost_nl_Primal +cost_gene_Primal 

@objective(model, Min, obj)  # single-level objective function
#-------Solve and Output Results
set_optimizer(model , Gurobi.Optimizer)
set_optimizer_attribute(model, "QCPDual", 1)
optimize!(model)

=#



model= Model()                     

@variable(model, yˢᴳ² >= 0)                         # status of SGs, buses:2,3,4,5,27,30.     ,Bin
@variable(model, yˢᴳ³ >= 0)            
@variable(model, yˢᴳ⁴ >= 0)           
@variable(model, yˢᴳ⁵ >= 0)           
@variable(model, yˢᴳ²⁷ >= 0)            
@variable(model, yˢᴳ³⁰ >= 0)  

@constraint(model, yˢᴳ² <= 1) 
@constraint(model, yˢᴳ³ <= 1) 
@constraint(model, yˢᴳ⁴ <= 1) 
@constraint(model, yˢᴳ⁵ <= 1) 
@constraint(model, yˢᴳ²⁷ <= 1) 
@constraint(model, yˢᴳ³⁰ <= 1) 

@variable(model, Cᵁ²>=0)                           # startup costs for SGs                                          
@variable(model, Cᵁ³>=0)                            
@variable(model, Cᵁ⁴>=0)                                
@variable(model, Cᵁ⁵>=0)                                
@variable(model, Cᵁ²⁷>=0)                               
@variable(model, Cᵁ³⁰>=0)                 

@variable(model, Pᴳ[1:30]) 
@variable(model, Qᴳ[1:30])  
                                               # startup costs and shutdown costs occured at hour t
@constraint(model, Cᵁ²>=(yˢᴳ²-yˢᴳ²_0)*Kˢᵗ[1])         
@constraint(model, Cᵁ³>=(yˢᴳ³-yˢᴳ³_0)*Kˢᵗ[2]) 
@constraint(model, Cᵁ⁴>=(yˢᴳ⁴-yˢᴳ⁴_0)*Kˢᵗ[3])
@constraint(model, Cᵁ⁵>=(yˢᴳ⁵-yˢᴳ⁵_0)*Kˢᵗ[4])
@constraint(model, Cᵁ²⁷>=(yˢᴳ²⁷-yˢᴳ²⁷_0)*Kˢᵗ[5])
@constraint(model, Cᵁ³⁰>=(yˢᴳ³⁰-yˢᴳ³⁰_0)*Kˢᵗ[6])


#--------------------------AC power flow constraints 
                                          # bounds for the active and reactive output of each bus
    @constraint(model, Pᴳ[2] <= yˢᴳ²*Pˢᴳₘₐₓ[1])    # bounds for active power of SGs: buses 2 3 4 5 27 30           
    @constraint(model, yˢᴳ²*Pˢᴳₘᵢₙ[1] <= Pᴳ[2])       
    @constraint(model, Pᴳ[3] <= yˢᴳ³*Pˢᴳₘₐₓ[2])       
    @constraint(model, yˢᴳ³*Pˢᴳₘᵢₙ[2] <= Pᴳ[3])         
    @constraint(model, Pᴳ[4] <= yˢᴳ⁴*Pˢᴳₘₐₓ[3])       
    @constraint(model, yˢᴳ⁴*Pˢᴳₘᵢₙ[3] <= Pᴳ[4])    
    @constraint(model, Pᴳ[5] <= yˢᴳ⁵*Pˢᴳₘₐₓ[4])       
    @constraint(model, yˢᴳ⁵*Pˢᴳₘᵢₙ[4] <= Pᴳ[5])
    @constraint(model, Pᴳ[27] <= yˢᴳ²⁷*Pˢᴳₘₐₓ[5])       
    @constraint(model, yˢᴳ²⁷*Pˢᴳₘᵢₙ[5] <=Pᴳ[27])
    @constraint(model, Pᴳ[30] <= yˢᴳ³⁰*Pˢᴳₘₐₓ[6])       
    @constraint(model, yˢᴳ³⁰*Pˢᴳₘᵢₙ[6] <= Pᴳ[30])

    @constraint(model, Pᴳ[1] <= α_VSG*Pⱽˢᴳₘₐₓ[1])            # bounds for active power of VSGs and IBGs: buses 1 23 24           
    @constraint(model, Pⱽˢᴳₘᵢₙ[1] <= Pᴳ[1])       
    @constraint(model, Pᴳ[23] <= α_IBG*Pᴵᴮᴳₘₐₓ[1])       
    @constraint(model, Pᴵᴮᴳₘᵢₙ[1] <= Pᴳ[23])         
    @constraint(model, Pᴳ[24] <= α_IBG*Pᴵᴮᴳₘₐₓ[2])       
    @constraint(model, Pᴵᴮᴳₘᵢₙ[2] <= Pᴳ[24])   

    @constraint(model, Qᴳ[2] <= yˢᴳ²*Qˢᴳₘₐₓ[1])    # bounds for REactive power of SGs: buses 2 3 4 5 27 30           
    @constraint(model, yˢᴳ²*Qˢᴳₘᵢₙ[1] <= Qᴳ[2])       
    @constraint(model, Qᴳ[3] <= yˢᴳ³*Qˢᴳₘₐₓ[2])       
    @constraint(model, yˢᴳ³*Qˢᴳₘᵢₙ[2] <= Qᴳ[3])         
    @constraint(model, Qᴳ[4] <= yˢᴳ⁴*Qˢᴳₘₐₓ[3])       
    @constraint(model, yˢᴳ⁴*Qˢᴳₘᵢₙ[3] <= Qᴳ[4])    
    @constraint(model, Qᴳ[5] <= yˢᴳ⁵*Qˢᴳₘₐₓ[4])       
    @constraint(model, yˢᴳ⁵*Qˢᴳₘᵢₙ[4] <= Qᴳ[5])
    @constraint(model, Qᴳ[27] <= yˢᴳ²⁷*Qˢᴳₘₐₓ[5])       
    @constraint(model, yˢᴳ²⁷*Qˢᴳₘᵢₙ[5] <=Qᴳ[27])
    @constraint(model, Qᴳ[30] <= yˢᴳ³⁰*Qˢᴳₘₐₓ[6])       
    @constraint(model, yˢᴳ³⁰*Qˢᴳₘᵢₙ[6] <= Qᴳ[30])           

    @constraint(model, Qᴳ[1] <= Qⱽˢᴳₘₐₓ[1])       
    @constraint(model, Qⱽˢᴳₘᵢₙ[1] <= Qᴳ[1])

    @constraint(model, Qᴳ[23] <= ratio*Qᴵᴮᴳₘₐₓ[1])       
    @constraint(model, ratio*Qᴵᴮᴳₘᵢₙ[1] <= Qᴳ[23])
    @constraint(model, Qᴳ[24] <= ratio*Qᴵᴮᴳₘₐₓ[2])       
    @constraint(model, ratio*Qᴵᴮᴳₘᵢₙ[2] <= Qᴳ[24])

    @constraint(model, [ 1, Pᴳ[2] /(Sᵐᵃˣ_gc[1]), Qᴳ[2] /(Sᵐᵃˣ_gc[1]) ] in SecondOrderCone()  ) 
    @constraint(model, [ 1, Pᴳ[3] /(Sᵐᵃˣ_gc[2]), Qᴳ[3] /(Sᵐᵃˣ_gc[2]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[4] /(Sᵐᵃˣ_gc[3]), Qᴳ[4] /(Sᵐᵃˣ_gc[3]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[5] /(Sᵐᵃˣ_gc[4]), Qᴳ[5] /(Sᵐᵃˣ_gc[4]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[27] /(Sᵐᵃˣ_gc[5]), Qᴳ[27] /(Sᵐᵃˣ_gc[5]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[30] /(Sᵐᵃˣ_gc[6]), Qᴳ[30] /(Sᵐᵃˣ_gc[6]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[1] /(Sᵐᵃˣ_gv[1]), Qᴳ[1] /(Sᵐᵃˣ_gv[1]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[23] /(Sᵐᵃˣ_c[1]), Qᴳ[23] /(Sᵐᵃˣ_c[1]) ] in SecondOrderCone()  )
    @constraint(model, [ 1, Pᴳ[24] /(Sᵐᵃˣ_c[2]), Qᴳ[24] /(Sᵐᵃˣ_c[2]) ] in SecondOrderCone()  )

    gens_SGs = [2, 3, 4, 5, 27, 30]    # buses with SGs
    gens_VSGs = [1]                    # buses with VSGs
    gens_IBGs = [23, 24]               # buses with IBGs gens_all = union(gens_SGs, gens_VSGs, gens_IBGs)

    for i in 1:bus_num                 # no power is injected into the non-generation bus 
        if !(i in gens_SGs) && !(i in gens_VSGs) && !(i in gens_IBGs)
            @constraint(model, Pᴳ[i] == 0)
            @constraint(model, Qᴳ[i] == 0)
        end
    end


@variable(model, Pˡⁱⁿᵉ[1:30,1:30])
@variable(model, Qˡⁱⁿᵉ[1:30,1:30])


for i in 1:bus_num

        @constraint(model, Pᴳ[i] - Pᴰ[i] ==  sum( Pˡⁱⁿᵉ[i,j] for j in bus_connection[i] ) )   #  *1.1
        @constraint(model, Qᴳ[i] - Qᴰ[i] ==  sum( Qˡⁱⁿᵉ[i,j] for j in bus_connection[i] ) )   

        for j in bus_connection[i]
            @constraint(model, Pˡⁱⁿᵉ[i,j] + Pˡⁱⁿᵉ[j,i] == 0)
            @constraint(model, Qˡⁱⁿᵉ[i,j] + Qˡⁱⁿᵉ[j,i] == 0)
        end    

end 



#--------------------------Voltage stability constraints--------------------------

@variable(model, Γ_23)
@variable(model, Γ_24)

@variable(model, ηₘ_1[1:15,1:T] >= 0)        # product of each pair of SG 
@constraint(model, ηₘ_1 .<= 1)
@variable(model, ηₘ_2[1:6,1:T])             # product of each pair of SG and VSG 


for t in 1:T
    # 15+6
    @constraint(model, ηₘ_1[1,t]>=yˢᴳ²+yˢᴳ³-1)  
    @constraint(model, ηₘ_1[1,t]<=yˢᴳ²)
    @constraint(model, ηₘ_1[1,t]<=yˢᴳ³)
    @constraint(model, ηₘ_1[2,t]>=yˢᴳ²+yˢᴳ⁴-1)  
    @constraint(model, ηₘ_1[2,t]<=yˢᴳ²)
    @constraint(model, ηₘ_1[2,t]<=yˢᴳ⁴)
    @constraint(model, ηₘ_1[3,t]>=yˢᴳ²+yˢᴳ⁵-1)  
    @constraint(model, ηₘ_1[3,t]<=yˢᴳ²)
    @constraint(model, ηₘ_1[3,t]<=yˢᴳ⁵)
    @constraint(model, ηₘ_1[4,t]>=yˢᴳ²+yˢᴳ²⁷-1)  
    @constraint(model, ηₘ_1[4,t]<=yˢᴳ²)
    @constraint(model, ηₘ_1[4,t]<=yˢᴳ²⁷)
    @constraint(model, ηₘ_1[5,t]>=yˢᴳ²+yˢᴳ³⁰-1)  
    @constraint(model, ηₘ_1[5,t]<=yˢᴳ²)
    @constraint(model, ηₘ_1[5,t]<=yˢᴳ³⁰)

    @constraint(model, ηₘ_1[6,t]>=yˢᴳ³+yˢᴳ⁴-1)  
    @constraint(model, ηₘ_1[6,t]<=yˢᴳ³)
    @constraint(model, ηₘ_1[6,t]<=yˢᴳ⁴)
    @constraint(model, ηₘ_1[7,t]>=yˢᴳ³+yˢᴳ⁵-1)  
    @constraint(model, ηₘ_1[7,t]<=yˢᴳ³)
    @constraint(model, ηₘ_1[7,t]<=yˢᴳ⁵)
    @constraint(model, ηₘ_1[8,t]>=yˢᴳ³+yˢᴳ²⁷-1)  
    @constraint(model, ηₘ_1[8,t]<=yˢᴳ³)
    @constraint(model, ηₘ_1[8,t]<=yˢᴳ²⁷)
    @constraint(model, ηₘ_1[9,t]>=yˢᴳ³+yˢᴳ³⁰-1)  
    @constraint(model, ηₘ_1[9,t]<=yˢᴳ³)
    @constraint(model, ηₘ_1[9,t]<=yˢᴳ³⁰)

    @constraint(model, ηₘ_1[10,t]>=yˢᴳ⁴+yˢᴳ⁵-1)  
    @constraint(model, ηₘ_1[10,t]<=yˢᴳ⁴)
    @constraint(model, ηₘ_1[10,t]<=yˢᴳ⁵)
    @constraint(model, ηₘ_1[11,t]>=yˢᴳ⁴+yˢᴳ²⁷-1)  
    @constraint(model, ηₘ_1[11,t]<=yˢᴳ⁴)
    @constraint(model, ηₘ_1[11,t]<=yˢᴳ²⁷)
    @constraint(model, ηₘ_1[12,t]>=yˢᴳ⁴+yˢᴳ³⁰-1)  
    @constraint(model, ηₘ_1[12,t]<=yˢᴳ⁴)
    @constraint(model, ηₘ_1[12,t]<=yˢᴳ³⁰)

    @constraint(model, ηₘ_1[13,t]>=yˢᴳ⁵+yˢᴳ²⁷-1)  
    @constraint(model, ηₘ_1[13,t]<=yˢᴳ⁵)
    @constraint(model, ηₘ_1[13,t]<=yˢᴳ²⁷)
    @constraint(model, ηₘ_1[14,t]>=yˢᴳ⁵+yˢᴳ³⁰-1)  
    @constraint(model, ηₘ_1[14,t]<=yˢᴳ⁵)
    @constraint(model, ηₘ_1[14,t]<=yˢᴳ³⁰)

    @constraint(model, ηₘ_1[15,t]>=yˢᴳ²⁷+yˢᴳ³⁰-1)  
    @constraint(model, ηₘ_1[15,t]<=yˢᴳ²⁷)
    @constraint(model, ηₘ_1[15,t]<=yˢᴳ³⁰)

    @constraint(model, ηₘ_2[1,t]==yˢᴳ²*α_VSG)
    @constraint(model, ηₘ_2[2,t]==yˢᴳ³*α_VSG)
    @constraint(model, ηₘ_2[3,t]==yˢᴳ⁴*α_VSG)
    @constraint(model, ηₘ_2[4,t]==yˢᴳ⁵*α_VSG)
    @constraint(model, ηₘ_2[5,t]==yˢᴳ²⁷*α_VSG)
    @constraint(model, ηₘ_2[6,t]==yˢᴳ³⁰*α_VSG)

    @constraint(model, Γ_23== 1/2 * ( ( K_c_gc[1,1] *yˢᴳ² +K_c_gc[1,2] *yˢᴳ³ +K_c_gc[1,3] *yˢᴳ⁴ 
    +K_c_gc[1,4] *yˢᴳ⁵ +K_c_gc[1,5] *yˢᴳ²⁷ +K_c_gc[1,6] *yˢᴳ³⁰ )   +α_VSG*K_c_gv[1,1]
    + sum(K_c_m[1,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[1,16:end] .*ηₘ_2[:,t]   ) ) )

    @constraint(model, Γ_24== 1/2 * ( ( K_c_gc[2,1] *yˢᴳ² +K_c_gc[2,2] *yˢᴳ³ +K_c_gc[2,3] *yˢᴳ⁴ 
    +K_c_gc[2,4] *yˢᴳ⁵ +K_c_gc[2,5] *yˢᴳ²⁷ +K_c_gc[2,6] *yˢᴳ³⁰ )   +α_VSG*K_c_gv[2,1]
    + sum(K_c_m[2,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[2,16:end] .*ηₘ_2[:,t]   ) ) )

    @constraint(model, [Qᴳ[23,t] + Γ_23*S_Base, Pᴳ[23,t], Qᴳ[23,t]] in SecondOrderCone()  )      
    @constraint(model, [Qᴳ[24,t] + Γ_24*S_Base, Pᴳ[24,t], Qᴳ[24,t]] in SecondOrderCone()  )

end


#--------------------------------
#  Model solving
#--------------------------------

cost_onoff_Primal=sum(Cᵁ²)+sum(Cᵁ³)+sum(Cᵁ⁴)+sum(Cᵁ⁵)+sum(Cᵁ²⁷)+sum(Cᵁ³⁰)
cost_nl_Primal=sum(Oⁿˡ[1].*(yˢᴳ²))+sum(Oⁿˡ[2].*(yˢᴳ³))+sum(Oⁿˡ[3].*(yˢᴳ⁴))+sum(Oⁿˡ[4].*(yˢᴳ⁵))+sum(Oⁿˡ[5].*(yˢᴳ²⁷))+sum(Oⁿˡ[6].*(yˢᴳ³⁰))    
cost_gene_Primal=sum(Oᵐ[1].*Pᴳ[2,:])+ sum(Oᵐ[2].*Pᴳ[3,:])+sum(Oᵐ[3].*Pᴳ[4,:])+sum(Oᵐ[4].*Pᴳ[5,:])+sum(Oᵐ[5].*Pᴳ[27,:])+sum(Oᵐ[6].*Pᴳ[30,:])   
obj=cost_onoff_Primal +cost_nl_Primal +cost_gene_Primal 

@objective(model, Min, obj)  # single-level objective function
#-------Solve and Output Results
set_optimizer(model , Gurobi.Optimizer)
set_optimizer_attribute(model, "QCPDual", 1)
optimize!(model)



#------------------------------
# Analyze results
#------------------------------
dual_model = dualize(model)  # Automatically create dual problem
set_optimizer(dual_model, Gurobi.Optimizer)
optimize!(dual_model) 




# =================================================================
# Print and save the DUAL MODEL structure (objective + constraints)
# Format is fully aligned with your original model printing code
# =================================================================

# Save full dual model to file
open("dual_model_summary.txt", "w") do io
    show(io, dual_model)
end

# Save detailed dual model (objective, variables, constraints)
open("dual_model_details.txt", "w") do io

    # --------------------------
    # Dual Objective Function
    # --------------------------
    println(io, "========================")
    println(io, "DUAL OBJECTIVE FUNCTION")
    println(io, "========================")
    println(io, objective_function(dual_model))
    println(io, "\n")

    # --------------------------
    # Dual Variables
    # --------------------------
    println(io, "========================")
    println(io, "DUAL VARIABLES")
    println(io, "========================")
    for var in all_variables(dual_model)
        println(io, var)
    end
    println(io, "\n")

    # --------------------------
    # Dual Constraints
    # --------------------------
    println(io, "========================")
    println(io, "DUAL CONSTRAINTS")
    println(io, "========================")

    # Print all constraint types in the dual model
    for (F, S) in list_of_constraint_types(dual_model)
        println(io, "\n[Constraint Type: $F in $S]")
        cons = all_constraints(dual_model, F, S)
        for con in cons
            println(io, "  ", con)
        end
    end

end

println("✅ Dual model saved to: dual_model_summary.txt + dual_model_details.txt")



yˢᴳ²_opt=JuMP.value.(yˢᴳ²)       # store optimal values of decision variables
yˢᴳ³_opt=JuMP.value.(yˢᴳ³)
yˢᴳ⁴_opt=JuMP.value.(yˢᴳ⁴)
yˢᴳ⁵_opt=JuMP.value.(yˢᴳ⁵)
yˢᴳ²⁷_opt=JuMP.value.(yˢᴳ²⁷)
yˢᴳ³⁰_opt=JuMP.value.(yˢᴳ³⁰)

zᵟ_23_opt=JuMP.value.(Γ_23)
zᵟ_24_opt=JuMP.value.(Γ_24)

UC_schedule= [yˢᴳ²_opt yˢᴳ³_opt yˢᴳ⁴_opt yˢᴳ⁵_opt yˢᴳ²⁷_opt yˢᴳ³⁰_opt]
Γ_23_1_opt = zeros(T)
Γ_24_1_opt = zeros(T)

for t in 1:T
    Γ_23_1_opt[t], Γ_24_1_opt[t] = calculate_Tao(α_VSG,Gᵥ,IBG,UC_schedule[t,:])  
end


plot(2*zᵟ_23_opt)
plot!(Γ_23_1_opt)



plot(2*zᵟ_24_opt)
plot!(Γ_24_1_opt)