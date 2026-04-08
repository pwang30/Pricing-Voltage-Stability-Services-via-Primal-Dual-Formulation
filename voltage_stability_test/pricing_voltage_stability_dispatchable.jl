import Pkg

using JuMP,Gurobi, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles,Plots,MAT,MathOptInterface,Dualization,HiGHS
include("dataset_gene.jl")
include("admittance_matrix_calculation.jl") 
include("offline_trainning.jl")
include("calculate_Tao.jl")
include("calculate_delta_gamma.jl")
include("calculate_ ZZ.jl")
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



Kˢᵗ=[125, 120, 115, 110, 105, 100]              #  Startup cost of SGs                     SGs, buses:2,3,4,5,27,30
Kˢʰ=[65, 60, 55, 50, 45, 40]            #  Shutdown cost of SGs                    SGs, buses:2,3,4,5,27,30
Oᵐ=[15.47, 14.47, 13.47, 12.47, 11.47, 10.47]                   #  Marginal generation cost of SGs   in    SGs, buses:2,3,4,5,27,30
Oⁿˡ=Oᵐ*10 

#-----------------
# Define model
#-----------------

model= Model()                     

@variable(model, yˢᴳ²[1:T]>=0)  # ,Bin  >=0                      # status of SGs, buses:2,3,4,5,27,30.    
@variable(model, yˢᴳ³[1:T]>=0)            
@variable(model, yˢᴳ⁴[1:T]>=0)           
@variable(model, yˢᴳ⁵[1:T]>=0)           
@variable(model, yˢᴳ²⁷[1:T]>=0)            
@variable(model, yˢᴳ³⁰[1:T]>=0)  

for t in 1:T  
    @constraint(model, yˢᴳ²[t]<=1)        
    @constraint(model, yˢᴳ³[t]<=1)        
    @constraint(model, yˢᴳ⁴[t]<=1)        
    @constraint(model, yˢᴳ⁵[t]<=1)        
    @constraint(model, yˢᴳ²⁷[t]<=1)        
    @constraint(model, yˢᴳ³⁰[t]<=1)        
end

@variable(model, Cᵁ²[1:T]>=0)                           # startup costs and shutdown costs for SGs          
@variable(model, Cᴰ²[1:T]>=0)                                  
@variable(model, Cᵁ³[1:T]>=0)               
@variable(model, Cᴰ³[1:T]>=0)                
@variable(model, Cᵁ⁴[1:T]>=0)                 
@variable(model, Cᴰ⁴[1:T]>=0)                 
@variable(model, Cᵁ⁵[1:T]>=0)                 
@variable(model, Cᴰ⁵[1:T]>=0)                 
@variable(model, Cᵁ²⁷[1:T]>=0)                 
@variable(model, Cᴰ²⁷[1:T]>=0)                
@variable(model, Cᵁ³⁰[1:T]>=0)                 
@variable(model, Cᴰ³⁰[1:T]>=0) 

@variable(model, Pᴳ[1:30,1:T]) 
@variable(model, Qᴳ[1:30,1:T])  

for t in 2:T                                                # startup costs and shutdown costs occured at hour t
    @constraint(model, Cᵁ²[t]>=(yˢᴳ²[t]-yˢᴳ²[t-1])*Kˢᵗ[1])        
    @constraint(model, Cᴰ²[t]>=(yˢᴳ²[t-1]-yˢᴳ²[t])*Kˢʰ[1])  
    @constraint(model, Cᵁ³[t]>=(yˢᴳ³[t]-yˢᴳ³[t-1])*Kˢᵗ[2]) 
    @constraint(model, Cᴰ³[t]>=(yˢᴳ³[t-1]-yˢᴳ³[t])*Kˢʰ[2])
    @constraint(model, Cᵁ⁴[t]>=(yˢᴳ⁴[t]-yˢᴳ⁴[t-1])*Kˢᵗ[3])
    @constraint(model, Cᴰ⁴[t]>=(yˢᴳ⁴[t-1]-yˢᴳ⁴[t])*Kˢʰ[3])
    @constraint(model, Cᵁ⁵[t]>=(yˢᴳ⁵[t]-yˢᴳ⁵[t-1])*Kˢᵗ[4])
    @constraint(model, Cᴰ⁵[t]>=(yˢᴳ⁵[t-1]-yˢᴳ⁵[t])*Kˢʰ[4])
    @constraint(model, Cᵁ²⁷[t]>=(yˢᴳ²⁷[t]-yˢᴳ²⁷[t-1])*Kˢᵗ[5])
    @constraint(model, Cᴰ²⁷[t]>=(yˢᴳ²⁷[t-1]-yˢᴳ²⁷[t])*Kˢʰ[5])
    @constraint(model, Cᵁ³⁰[t]>=(yˢᴳ³⁰[t]-yˢᴳ³⁰[t-1])*Kˢᵗ[6])
    @constraint(model, Cᴰ³⁰[t]>=(yˢᴳ³⁰[t-1]-yˢᴳ³⁰[t])*Kˢʰ[6])
end   



#---------------------------------------------------------------------------------------------------------------------------
#  AC power flow constraints
#---------------------------------------------------------------------------------------------------------------------------
#@variable(model, anci_var)
#@constraint(model, anci_var ==1)

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

@variable(model, Pˡⁱⁿᵉ_abs[1:30,1:30,1:T] >=0)
@variable(model, Qˡⁱⁿᵉ_abs[1:30,1:30,1:T] >=0)

Power_balance_P = [Dict{Int, ConstraintRef}() for _ in 1:bus_num]
Power_balance_Q = [Dict{Int, ConstraintRef}() for _ in 1:bus_num]

for i in 1:bus_num
    for t in 1:T

        Power_balance_P[i][t] = @constraint(model, Pᴳ[i,t] - Pᴰ[i,t] ==  sum( Pˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # p_i^g-p_i^d == G_{i,i}*c_{i,i}+...   eq. 50(a)
        Power_balance_Q[i][t] = @constraint(model, Qᴳ[i,t] - Qᴰ[i,t] ==  sum( Qˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # q_i^g-q_i^d == -B_{i,i}*c_{i,i}+...  eq. 50(b)

        for j in bus_connection[i]
            @constraint(model, Pˡⁱⁿᵉ[i,j,t] + Pˡⁱⁿᵉ[j,i,t] == 0)
            @constraint(model, Qˡⁱⁿᵉ[i,j,t] + Qˡⁱⁿᵉ[j,i,t] == 0)

            @constraint(model, Pˡⁱⁿᵉ[i,j,t] <= Pˡⁱⁿᵉ_abs[i,j,t])
            @constraint(model, -Pˡⁱⁿᵉ[i,j,t] <= Pˡⁱⁿᵉ_abs[i,j,t])
            @constraint(model, Qˡⁱⁿᵉ[i,j,t] <= Qˡⁱⁿᵉ_abs[i,j,t])
            @constraint(model, -Qˡⁱⁿᵉ[i,j,t] <= Qˡⁱⁿᵉ_abs[i,j,t])
            @constraint(model, Pˡⁱⁿᵉ_abs[i,j,t] <= S_lines[i,j])
            @constraint(model, Qˡⁱⁿᵉ_abs[i,j,t] <= S_lines[i,j])
            @constraint(model, Pˡⁱⁿᵉ_abs[i,j,t] + Qˡⁱⁿᵉ_abs[i,j,t] <= S_lines[i,j]*sqrt(2))

        end    
    end
end 




#--------------------------------
#  Voltage stability constraints
#--------------------------------

@variable(model, P_23_eq[1:T])              # equivalent P and Q of IBGs, buses:23,24    and Γ_c  
@variable(model, P_24_eq[1:T]) 
@variable(model, Q_23_eq[1:T])                 
@variable(model, Q_24_eq[1:T])
@variable(model, Γ_23[1:T])
@variable(model, Γ_24[1:T])

@variable(model, ηₘ_1[1:15,1:T]>=0)        # product of each pair of SG 
@variable(model, ηₘ_2[1:6,1:T])             # product of each pair of SG and VSG 
for i in 1:15, j in 1:T
   @constraint(model, ηₘ_1[i,j] <=1)        
end


Voltage_constraint_23=Dict()
Voltage_constraint_24=Dict()

@variable(model, var_Mc_g2_P23[1:T])
@variable(model, var_Mc_g3_P23[1:T])
@variable(model, var_Mc_g4_P23[1:T])
@variable(model, var_Mc_g5_P23[1:T])
@variable(model, var_Mc_g27_P23[1:T])
@variable(model, var_Mc_g30_P23[1:T])
@variable(model, var_Mc_g2_P24[1:T])
@variable(model, var_Mc_g3_P24[1:T])
@variable(model, var_Mc_g4_P24[1:T])
@variable(model, var_Mc_g5_P24[1:T])
@variable(model, var_Mc_g27_P24[1:T])
@variable(model, var_Mc_g30_P24[1:T])

@variable(model, var_Mc_g2_Q23[1:T])
@variable(model, var_Mc_g3_Q23[1:T])
@variable(model, var_Mc_g4_Q23[1:T])
@variable(model, var_Mc_g5_Q23[1:T])
@variable(model, var_Mc_g27_Q23[1:T])
@variable(model, var_Mc_g30_Q23[1:T])
@variable(model, var_Mc_g2_Q24[1:T])
@variable(model, var_Mc_g3_Q24[1:T])
@variable(model, var_Mc_g4_Q24[1:T])
@variable(model, var_Mc_g5_Q24[1:T])
@variable(model, var_Mc_g27_Q24[1:T])
@variable(model, var_Mc_g30_Q24[1:T])

@variable(model, var_Mc_ηₘ_1_P23[1:15,1:T])
@variable(model, var_Mc_ηₘ_1_P24[1:15,1:T])
@variable(model, var_Mc_ηₘ_1_Q23[1:15,1:T])
@variable(model, var_Mc_ηₘ_1_Q24[1:15,1:T])
@variable(model, var_Mc_ηₘ_2_P23[1:6,1:T])
@variable(model, var_Mc_ηₘ_2_P24[1:6,1:T])
@variable(model, var_Mc_ηₘ_2_Q23[1:6,1:T])
@variable(model, var_Mc_ηₘ_2_Q24[1:6,1:T])

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


    @constraint(model, var_Mc_g2_P23[t] <=yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g2_P23[t] <= Pᴳ[23,t] +yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g2_P23[t] >=Pᴳ[23,t]+yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g2_P23[t] >= yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g3_P23[t] <=yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g3_P23[t] <= Pᴳ[23,t] +yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g3_P23[t] >=Pᴳ[23,t]+yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g3_P23[t] >= yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g4_P23[t] <=yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g4_P23[t] <= Pᴳ[23,t] +yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g4_P23[t] >=Pᴳ[23,t]+yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g4_P23[t] >= yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g5_P23[t] <=yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g5_P23[t] <= Pᴳ[23,t] +yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g5_P23[t] >=Pᴳ[23,t]+yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])   
    @constraint(model, var_Mc_g5_P23[t] >= yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g27_P23[t] <=yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g27_P23[t] <= Pᴳ[23,t] +yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g27_P23[t] >=Pᴳ[23,t]+yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])   
    @constraint(model, var_Mc_g27_P23[t] >= yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g30_P23[t] <=yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g30_P23[t] <= Pᴳ[23,t] +yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g30_P23[t] >=Pᴳ[23,t]+yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])   
    @constraint(model, var_Mc_g30_P23[t] >= yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g2_P24[t] <=yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g2_P24[t] <= Pᴳ[24,t] +yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g2_P24[t] >=Pᴳ[24,t]+yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g2_P24[t] >= yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_g3_P24[t] <=yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g3_P24[t] <= Pᴳ[24,t] +yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g3_P24[t] >=Pᴳ[24,t]+yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g3_P24[t] >= yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_g4_P24[t] <=yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g4_P24[t] <= Pᴳ[24,t] +yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g4_P24[t] >=Pᴳ[24,t]+yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g4_P24[t] >= yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_g5_P24[t] <=yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g5_P24[t] <= Pᴳ[24,t] +yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g5_P24[t] >=Pᴳ[24,t]+yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])   
    @constraint(model, var_Mc_g5_P24[t] >= yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_g27_P24[t] <=yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g27_P24[t] <= Pᴳ[24,t] +yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g27_P24[t] >=Pᴳ[24,t]+yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])   
    @constraint(model, var_Mc_g27_P24[t] >= yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_g30_P24[t] <=yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g30_P24[t] <= Pᴳ[24,t] +yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g30_P24[t] >=Pᴳ[24,t]+yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])   
    @constraint(model, var_Mc_g30_P24[t] >= yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[2] )

    
    
    @constraint(model, var_Mc_g2_Q23[t] <= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g2_Q23[t] <= Qᴳ[23,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g2_Q23[t] >= Qᴳ[23,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g2_Q23[t] >= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g3_Q23[t] <= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g3_Q23[t] <= Qᴳ[23,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g3_Q23[t] >= Qᴳ[23,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g3_Q23[t] >= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g4_Q23[t] <= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g4_Q23[t] <= Qᴳ[23,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g4_Q23[t] >= Qᴳ[23,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g4_Q23[t] >= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g5_Q23[t] <= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g5_Q23[t] <= Qᴳ[23,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g5_Q23[t] >= Qᴳ[23,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])   
    @constraint(model, var_Mc_g5_Q23[t] >= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g27_Q23[t] <= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g27_Q23[t] <= Qᴳ[23,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g27_Q23[t] >= Qᴳ[23,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])   
    @constraint(model, var_Mc_g27_Q23[t] >= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g30_Q23[t] <= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model, var_Mc_g30_Q23[t] <= Qᴳ[23,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_g30_Q23[t] >= Qᴳ[23,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])   
    @constraint(model, var_Mc_g30_Q23[t] >= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_g2_Q24[t] <= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g2_Q24[t] <= Qᴳ[24,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g2_Q24[t] >= Qᴳ[24,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g2_Q24[t] >= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )   

    @constraint(model, var_Mc_g3_Q24[t] <= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g3_Q24[t] <= Qᴳ[24,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g3_Q24[t] >= Qᴳ[24,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g3_Q24[t] >= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_g4_Q24[t] <= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g4_Q24[t] <= Qᴳ[24,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g4_Q24[t] >= Qᴳ[24,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g4_Q24[t] >= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_g5_Q24[t] <= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g5_Q24[t] <= Qᴳ[24,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g5_Q24[t] >= Qᴳ[24,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])   
    @constraint(model, var_Mc_g5_Q24[t] >= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_g27_Q24[t] <= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g27_Q24[t] <= Qᴳ[24,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] ) 
    @constraint(model, var_Mc_g27_Q24[t] >= Qᴳ[24,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])   
    @constraint(model, var_Mc_g27_Q24[t] >= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_g30_Q24[t] <= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model, var_Mc_g30_Q24[t] <= Qᴳ[24,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_g30_Q24[t] >= Qᴳ[24,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])   
    @constraint(model, var_Mc_g30_Q24[t] >= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_ηₘ_1_P23[:,t] <= ηₘ_1[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model, var_Mc_ηₘ_1_P23[:,t] <= Pᴳ[23,t] .+ ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[1] .- Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_ηₘ_1_P23[:,t] >= Pᴳ[23,t] .+ ηₘ_1[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[1] .- α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model, var_Mc_ηₘ_1_P23[:,t] >= ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_ηₘ_2_P23[:,t] <= ηₘ_2[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model, var_Mc_ηₘ_2_P23[:,t] <= Pᴳ[23,t] .+ ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[1] .- Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_ηₘ_2_P23[:,t] >= Pᴳ[23,t] .+ ηₘ_2[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[1] .- α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model, var_Mc_ηₘ_2_P23[:,t] >= ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_ηₘ_1_P24[:,t] <= ηₘ_1[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model, var_Mc_ηₘ_1_P24[:,t] <= Pᴳ[24,t] .+ ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[2] .- Pᴵᴮᴳₘᵢₙ[2] )    
    @constraint(model, var_Mc_ηₘ_1_P24[:,t] >= Pᴳ[24,t] .+ ηₘ_1[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[2] .- α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model, var_Mc_ηₘ_1_P24[:,t] >= ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_ηₘ_2_P24[:,t] <= ηₘ_2[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model, var_Mc_ηₘ_2_P24[:,t] <= Pᴳ[24,t] .+ ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[2] .- Pᴵᴮᴳₘᵢₙ[2] )    
    @constraint(model, var_Mc_ηₘ_2_P24[:,t] >= Pᴳ[24,t] .+ ηₘ_2[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[2] .- α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model, var_Mc_ηₘ_2_P24[:,t] >= ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_ηₘ_1_Q23[:,t] <= ηₘ_1[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[1] )  
    @constraint(model, var_Mc_ηₘ_1_Q23[:,t] <= Qᴳ[23,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] .- ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_ηₘ_1_Q23[:,t] >= Qᴳ[23,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘₐₓ[1] .- ratio*Qᴵᴮᴳₘₐₓ[1] )
    @constraint(model, var_Mc_ηₘ_1_Q23[:,t] >= ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] )
    
    @constraint(model, var_Mc_ηₘ_2_Q23[:,t] <= ηₘ_2[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[1] )  
    @constraint(model, var_Mc_ηₘ_2_Q23[:,t] <= Qᴳ[23,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] .- ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model, var_Mc_ηₘ_2_Q23[:,t] >= Qᴳ[23,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘₐₓ[1] .- ratio*Qᴵᴮᴳₘₐₓ[1] )
    @constraint(model, var_Mc_ηₘ_2_Q23[:,t] >= ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model, var_Mc_ηₘ_1_Q24[:,t] <= ηₘ_1[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[2] )  
    @constraint(model, var_Mc_ηₘ_1_Q24[:,t] <= Qᴳ[24,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] .- ratio*Qᴵᴮᴳₘᵢₙ[2] )    
    @constraint(model, var_Mc_ηₘ_1_Q24[:,t] >= Qᴳ[24,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘₐₓ[2] .- ratio*Qᴵᴮᴳₘₐₓ[2] )
    @constraint(model, var_Mc_ηₘ_1_Q24[:,t] >= ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model, var_Mc_ηₘ_2_Q24[:,t] <= ηₘ_2[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[2] )  
    @constraint(model, var_Mc_ηₘ_2_Q24[:,t] <= Qᴳ[24,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] .- ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model, var_Mc_ηₘ_2_Q24[:,t] >= Qᴳ[24,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘₐₓ[2] .- ratio*Qᴵᴮᴳₘₐₓ[2] )
    @constraint(model, var_Mc_ηₘ_2_Q24[:,t] >= ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] )




    @constraint(model, P_23_eq[t]== Pᴳ[23,t]+  ( sh_1*Pᴳ[24,t]- ( K_c_gc[1,1] *var_Mc_g2_P24[t] +K_c_gc[1,2] *var_Mc_g3_P24[t] +K_c_gc[1,3] *var_Mc_g4_P24[t] 
    +K_c_gc[1,4] *var_Mc_g5_P24[t] +K_c_gc[1,5] *var_Mc_g27_P24[t] +K_c_gc[1,6] *var_Mc_g30_P24[t]    +α_VSG*K_c_gv[1,1]*Pᴳ[24,t] 
    + sum(K_c_m[1,1:15] .*var_Mc_ηₘ_1_P24[:,t]) +  sum(K_c_m[1,16:end] .*var_Mc_ηₘ_2_P24[:,t] ) ) )    )


    @constraint(model, P_24_eq[t]== Pᴳ[24,t]+  ( sh_2*Pᴳ[23,t]- ( K_c_gc[2,1] *var_Mc_g2_P23[t] +K_c_gc[2,2] *var_Mc_g3_P23[t] +K_c_gc[2,3] *var_Mc_g4_P23[t] 
    +K_c_gc[2,4] *var_Mc_g5_P23[t] +K_c_gc[2,5] *var_Mc_g27_P23[t] +K_c_gc[2,6] *var_Mc_g30_P23[t]    +α_VSG*K_c_gv[2,1]*Pᴳ[23,t]
    + sum(K_c_m[2,1:15] .*var_Mc_ηₘ_1_P23[:,t]) +  sum(K_c_m[2,16:end] .*var_Mc_ηₘ_2_P23[:,t] ) ) )  )

    @constraint(model, Q_23_eq[t]== Qᴳ[23,t]+  ( sh_1*Qᴳ[24,t]- ( K_c_gc[1,1] *var_Mc_g2_Q24[t] +K_c_gc[1,2] *var_Mc_g3_Q24[t] +K_c_gc[1,3] *var_Mc_g4_Q24[t] 
    +K_c_gc[1,4] *var_Mc_g5_Q24[t] +K_c_gc[1,5] *var_Mc_g27_Q24[t] +K_c_gc[1,6] *var_Mc_g30_Q24[t]    +α_VSG*K_c_gv[1,1]*Qᴳ[24,t] 
    + sum(K_c_m[1,1:15] .*var_Mc_ηₘ_1_Q24[:,t]) +  sum(K_c_m[1,16:end] .*var_Mc_ηₘ_2_Q24[:,t] ) ) )    )

    @constraint(model, Q_24_eq[t]== Qᴳ[24,t]+  ( sh_2*Qᴳ[23,t]- ( K_c_gc[2,1] *var_Mc_g2_Q23[t] +K_c_gc[2,2] *var_Mc_g3_Q23[t] +K_c_gc[2,3] *var_Mc_g4_Q23[t] 
    +K_c_gc[2,4] *var_Mc_g5_Q23[t] +K_c_gc[2,5] *var_Mc_g27_Q23[t] +K_c_gc[2,6] *var_Mc_g30_Q23[t]    +α_VSG*K_c_gv[2,1]*Pᴳ[23,t]
    + sum(K_c_m[2,1:15] .*var_Mc_ηₘ_1_Q23[:,t]) +  sum(K_c_m[2,16:end] .*var_Mc_ηₘ_2_Q23[:,t] ) ) )    )

    
    @constraint(model, Γ_23[t]== 1/2 * ( K_c_gc[3,1] *yˢᴳ²[t] +K_c_gc[3,2] *yˢᴳ³[t] +K_c_gc[3,3] *yˢᴳ⁴[t] 
    +K_c_gc[3,4] *yˢᴳ⁵[t] +K_c_gc[3,5] *yˢᴳ²⁷[t] +K_c_gc[3,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[3,16:end] .*ηₘ_2[:,t]   ) )

    @constraint(model, Γ_24[t]== 1/2 * ( K_c_gc[4,1] *yˢᴳ²[t] +K_c_gc[4,2] *yˢᴳ³[t] +K_c_gc[4,3] *yˢᴳ⁴[t] 
    +K_c_gc[4,4] *yˢᴳ⁵[t] +K_c_gc[4,5] *yˢᴳ²⁷[t] +K_c_gc[4,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[4,1]
    + sum(K_c_m[4,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[4,16:end] .*ηₘ_2[:,t]   ) )

    Voltage_constraint_23[t] = @constraint(model, [Q_23_eq[t] + Γ_23[t]*100, P_23_eq[t], Q_23_eq[t]] in SecondOrderCone()  )      
    Voltage_constraint_24[t] = @constraint(model, [Q_24_eq[t] + Γ_24[t]*100, P_24_eq[t], Q_24_eq[t]] in SecondOrderCone()  )

end



#--------------------------------
#  Model solving
#--------------------------------

cost_onoff_Primal=sum(Cᵁ²)+sum(Cᴰ²)+sum(Cᵁ³)+sum(Cᴰ³)+sum(Cᵁ⁴)+sum(Cᴰ⁴)+sum(Cᵁ⁵)+sum(Cᴰ⁵)+sum(Cᵁ²⁷)+sum(Cᴰ²⁷)+sum(Cᵁ³⁰)+sum(Cᴰ³⁰)  
cost_nl_Primal=sum(Oⁿˡ[1].*(yˢᴳ²))+sum(Oⁿˡ[2].*(yˢᴳ³))+sum(Oⁿˡ[3].*(yˢᴳ⁴))+sum(Oⁿˡ[4].*(yˢᴳ⁵))+sum(Oⁿˡ[5].*(yˢᴳ²⁷))+sum(Oⁿˡ[6].*(yˢᴳ³⁰))    
cost_gene_Primal=sum(Oᵐ[1].*Pᴳ[2,:])+ sum(Oᵐ[2].*Pᴳ[3,:])+sum(Oᵐ[3].*Pᴳ[4,:])+sum(Oᵐ[4].*Pᴳ[5,:])+sum(Oᵐ[5].*Pᴳ[27,:])+sum(Oᵐ[6].*Pᴳ[30,:])   
obj=cost_onoff_Primal +cost_nl_Primal +cost_gene_Primal 
@objective(model, Min, cost_onoff_Primal +cost_nl_Primal +cost_gene_Primal )  # single-level objective function - 1.06802712e+04

set_optimizer(model , Gurobi.Optimizer)
set_optimizer_attribute(model, "QCPDual", 1)
#set_optimizer_attribute(model, "Method", 1)
#set_optimizer_attribute(model, "BarHomogeneous", 1) 
optimize!(model) 
# optimal obj of primal model ==  1.347853305046e+04




#------------------------------
# Analyze results
#------------------------------
energy_prices_P=zeros(bus_num,T)
energy_prices_P = [dual(Power_balance_P[i][t]) for i in 1:bus_num, t in 1:T]

bar(energy_prices_P')

println(energy_prices_P[1,:])

#------------------------------
# Calculate energy profits of units
#------------------------------
Pᴳ=JuMP.value.(Pᴳ)

Cᵁ²=value.(Cᵁ²)
Cᴰ²=value.(Cᴰ²)
Cᵁ³=value.(Cᵁ³)
Cᴰ³=value.(Cᴰ³)
Cᵁ⁴=value.(Cᵁ⁴)
Cᴰ⁴=value.(Cᴰ⁴)
Cᵁ⁵=value.(Cᵁ⁵)
Cᴰ⁵=value.(Cᴰ⁵)
Cᵁ²⁷=value.(Cᵁ²⁷)
Cᴰ²⁷=value.(Cᴰ²⁷)
Cᵁ³⁰=value.(Cᵁ³⁰)
Cᴰ³⁰=value.(Cᴰ³⁰)

yˢᴳ²=JuMP.value.(yˢᴳ²)
yˢᴳ³=JuMP.value.(yˢᴳ³)
yˢᴳ⁴=JuMP.value.(yˢᴳ⁴)
yˢᴳ⁵=JuMP.value.(yˢᴳ⁵)
yˢᴳ²⁷=JuMP.value.(yˢᴳ²⁷)
yˢᴳ³⁰=JuMP.value.(yˢᴳ³⁰)

profit_VSG_1=sum(Pᴳ[1,:] .*energy_prices_P[1,:])
profitˢᴳ² = sum(Pᴳ[2,:] .*energy_prices_P[2,:]) -sum(Pᴳ[2,:] *Oᵐ[1]) -sum(Cᵁ²) -sum(Cᴰ²) -sum(yˢᴳ²*Oⁿˡ[1])
profitˢᴳ³ = sum(Pᴳ[3,:] .*energy_prices_P[3,:]) -sum(Pᴳ[3,:] *Oᵐ[2]) -sum(Cᵁ³) -sum(Cᴰ³) -sum(yˢᴳ³*Oⁿˡ[2])
profitˢᴳ⁴ = sum(Pᴳ[4,:] .*energy_prices_P[4,:]) -sum(Pᴳ[4,:] *Oᵐ[3]) -sum(Cᵁ⁴) -sum(Cᴰ⁴) -sum(yˢᴳ⁴*Oⁿˡ[3])
profitˢᴳ⁵ = sum(Pᴳ[5,:] .*energy_prices_P[5,:]) -sum(Pᴳ[5,:] *Oᵐ[4]) -sum(Cᵁ⁵) -sum(Cᴰ⁵) -sum(yˢᴳ⁵*Oⁿˡ[4])
profitˢᴳ²⁷ = sum(Pᴳ[27,:] .*energy_prices_P[27,:]) -sum(Pᴳ[27,:] *Oᵐ[5]) -sum(Cᵁ²⁷) -sum(Cᴰ²⁷) -sum(yˢᴳ²⁷*Oⁿˡ[5])
profitˢᴳ³⁰ = sum(Pᴳ[30,:] .*energy_prices_P[30,:]) -sum(Pᴳ[30,:] *Oᵐ[6]) -sum(Cᵁ³⁰) -sum(Cᴰ³⁰) -sum(yˢᴳ³⁰*Oⁿˡ[6])
profit_IBG_23=sum(Pᴳ[23,:] .*energy_prices_P[23,:])
profit_IBG_24=sum(Pᴳ[24,:] .*energy_prices_P[24,:])

plot(Pᴳ[2,:] .*energy_prices_P[2,:] .-Pᴳ[2,:] *Oᵐ[1] .-Cᵁ² .-Cᴰ² .-yˢᴳ²*Oⁿˡ[1])
plot!(Pᴳ[3,:] .*energy_prices_P[3,:] .-Pᴳ[3,:] *Oᵐ[2] .-Cᵁ³ .-Cᴰ³ .-yˢᴳ³*Oⁿˡ[2])
plot!(Pᴳ[4,:] .*energy_prices_P[4,:] .-Pᴳ[4,:] *Oᵐ[3] .-Cᵁ⁴ .-Cᴰ⁴ .-yˢᴳ⁴*Oⁿˡ[3])
plot!(Pᴳ[5,:] .*energy_prices_P[5,:] .-Pᴳ[5,:] *Oᵐ[4] .-Cᵁ⁵ .-Cᴰ⁵ .-yˢᴳ⁵*Oⁿˡ[4])
plot!(Pᴳ[27,:] .*energy_prices_P[27,:] .-Pᴳ[27,:] *Oᵐ[5] .-Cᵁ²⁷ .-Cᴰ²⁷ .-yˢᴳ²⁷*Oⁿˡ[5])
plot!(Pᴳ[30,:] .*energy_prices_P[30,:] .-Pᴳ[30,:] *Oᵐ[6] .-Cᵁ³⁰ .-Cᴰ³⁰ .-yˢᴳ³⁰*Oⁿˡ[6])








#------------------------------
# Calculate VS service revenues of units
#------------------------------
Pᴳ=JuMP.value.(Pᴳ)
Qᴳ=JuMP.value.(Qᴳ)

vs_price_23 = [dual(Voltage_constraint_23[t]) for t in 1:T]
vs_price_24 = [dual(Voltage_constraint_24[t]) for t in 1:T]

a=hcat(vs_price_23)
b=hcat(vs_price_24)

mu=zeros(T,2)   # IBGs at bus 23 and 24
lambda_1=zeros(T,2)
lambda_2=zeros(T,2)

for t in 1:T
    aa=a[t,:]
    bb=b[t,:]
    aaa = aa[1]
    bbb = bb[1]
    mu[t,1] = aaa[1]
    lambda_1[t,1] = aaa[2]
    lambda_2[t,1] = aaa[3]
    mu[t,2] = bbb[1]
    lambda_1[t,2] = bbb[2]
    lambda_2[t,2] = bbb[3]
end


plot(mu[:,1])
plot!(mu[:,2])

plot!(-lambda_1[:,1])
plot!(-lambda_1[:,2])

plot!(-lambda_2[:,1]+mu[:,1])
plot!(-lambda_2[:,2]+mu[:,2])


sum(mu[:,1])/T
sum(mu[:,2])/T
sum(-lambda_1[:,1])/T
sum(-lambda_1[:,2])/T
sum(-lambda_2[:,1]+mu[:,1])/T
sum(-lambda_2[:,2]+mu[:,2])/T

println(mu[:,2])

pppp=-lambda_2[:,1]

plot(-lambda_1[:,1])
plot!(mu[:,1])
plot!(-lambda_2[:,1]+mu[:,1])
plot!(-lambda_2[:,2]+mu[:,2])

#----revenue for IBGs
yˢᴳ²=JuMP.value.(yˢᴳ²)
yˢᴳ³=JuMP.value.(yˢᴳ³)
yˢᴳ⁴=JuMP.value.(yˢᴳ⁴)
yˢᴳ⁵=JuMP.value.(yˢᴳ⁵)
yˢᴳ²⁷=JuMP.value.(yˢᴳ²⁷)
yˢᴳ³⁰=JuMP.value.(yˢᴳ³⁰)

Pᴳ=JuMP.value.(Pᴳ)
Qᴳ=JuMP.value.(Qᴳ)

zᵟ_23_24_list=zeros(T)
zᵟ_24_23_list=zeros(T)

revenue_IBG_23_give_active=zeros(T)
revenue_IBG_24_give_active=zeros(T)
revenue_IBG_23_give_reactive=zeros(T)
revenue_IBG_24_give_reactive=zeros(T)

UC_schedule = zeros(6,T)

for t in 1:T
    UC_schedule[:,t]=[yˢᴳ²[t];yˢᴳ³[t];yˢᴳ⁴[t];yˢᴳ⁵[t];yˢᴳ²⁷[t];yˢᴳ³⁰[t]]
    zᵟ_23_24, zᵟ_24_23 =  calculate_ZZ(α_VSG,Gᵥ,IBG,UC_schedule[:,t]) 
    zᵟ_23_24_list[t]=zᵟ_23_24
    zᵟ_24_23_list[t]=zᵟ_24_23
    #P_23_eqivalent_list[t]=Pᴳ[23,t] + zᵟ_23_24_list[t] *Pᴳ[24,t]
    #P_24_eqivalent_list[t]=Pᴳ[24,t] + zᵟ_24_23_list[t] *Pᴳ[23,t]
    #Q_23_eqivalent_list[t]=Qᴳ[23,t] + zᵟ_23_24_list[t] *Qᴳ[24,t]
    #Q_24_eqivalent_list[t]=Qᴳ[24,t] + zᵟ_24_23_list[t] *Qᴳ[23,t]

    revenue_IBG_23_give_active[t]= (-lambda_1[t,1] * Pᴳ[23,t] )- (lambda_1[t,2] *zᵟ_24_23_list[t] * Pᴳ[23,t]  )
    revenue_IBG_24_give_active[t]= (-lambda_1[t,2] * Pᴳ[24,t] )-  (lambda_1[t,1] *zᵟ_23_24_list[t] * Pᴳ[24,t]  )

    revenue_IBG_23_give_reactive[t]=  (-lambda_2[t,1] + mu[t,1]  )* Qᴳ[23,t] - (lambda_2[t,2] - mu[t,2]  ) *zᵟ_24_23_list[t] * Qᴳ[23,t]  
    revenue_IBG_24_give_reactive[t]=  (-lambda_2[t,2] + mu[t,2]  )* Qᴳ[24,t] - (lambda_2[t,1] - mu[t,1]  ) *zᵟ_23_24_list[t] * Qᴳ[24,t]  
end

plot(revenue_IBG_23_give_active)
plot!(revenue_IBG_24_give_active)
plot!(revenue_IBG_23_give_reactive)
plot!(revenue_IBG_24_give_reactive)

println(revenue_IBG_23_give_reactive)

sum(revenue_IBG_23_give_active)+sum(revenue_IBG_23_give_reactive)

sum(revenue_IBG_24_give_active)+sum(revenue_IBG_24_give_reactive)




#----revenue for (V)SGs
VS_revenue_VSG_1=zeros(2,T)
VS_revenue_SGs_23=zeros(6,T)
VS_revenue_SGs_24=zeros(6,T)
UC_schedule = zeros(6,T)

for i in [1,2,3,4,5,6]
    for t in 1:T
        #α_VSG=0.9
        UC_schedule[:,t]=[yˢᴳ²[t];yˢᴳ³[t];yˢᴳ⁴[t];yˢᴳ⁵[t];yˢᴳ²⁷[t];yˢᴳ³⁰[t]]
        zᵟ_23_1, zᵟ_24_1 =  calculate_Tao(α_VSG,Gᵥ,IBG,UC_schedule[:,t]) 
        UC_schedule[i,t] = 0
        #α_VSG=0
        zᵟ_23_1_delta, zᵟ_24_1_delta = calculate_delta_gamma(α_VSG,Gᵥ,IBG,UC_schedule[:,t])
        VS_revenue_SGs_23[i,t] = mu[t,1]* (zᵟ_23_1 - zᵟ_23_1_delta) *100/2
        VS_revenue_SGs_24[i,t] = mu[t,2]* (zᵟ_24_1 - zᵟ_24_1_delta) *100/2
        #VS_revenue_VSG_1[1,t] =  mu[t,1]* (zᵟ_23_1 - zᵟ_23_1_delta) *100/2
        #VS_revenue_VSG_1[2,t] =  mu[t,2]* (zᵟ_24_1 - zᵟ_24_1_delta) *100/2
    end
end

sum(VS_revenue_SGs_23[1,:]) + sum(VS_revenue_SGs_24[1,:]) +profitˢᴳ²
sum(VS_revenue_SGs_23[2,:]) + sum(VS_revenue_SGs_24[2,:]) +profitˢᴳ³
sum(VS_revenue_SGs_23[3,:]) + sum(VS_revenue_SGs_24[3,:]) +profitˢᴳ⁴
sum(VS_revenue_SGs_23[4,:]) + sum(VS_revenue_SGs_24[4,:]) +profitˢᴳ⁵
sum(VS_revenue_SGs_23[5,:]) + sum(VS_revenue_SGs_24[5,:]) +profitˢᴳ²⁷
sum(VS_revenue_SGs_23[6,:]) + sum(VS_revenue_SGs_24[6,:]) +profitˢᴳ³⁰
sum(VS_revenue_VSG_1[1,:]) + sum(VS_revenue_VSG_1[2,:])


plot( VS_revenue_VSG_1[1,:]+ VS_revenue_VSG_1[2,:] )
plot!( VS_revenue_SGs_23[1,:]+ VS_revenue_SGs_24[1,:] )
plot!( VS_revenue_SGs_23[2,:]+ VS_revenue_SGs_24[2,:] )
plot!( VS_revenue_SGs_23[3,:]+ VS_revenue_SGs_24[3,:] )
plot!( VS_revenue_SGs_23[4,:]+ VS_revenue_SGs_24[4,:] )
plot!( VS_revenue_SGs_23[5,:]+ VS_revenue_SGs_24[5,:] )
plot!( VS_revenue_SGs_23[6,:]+ VS_revenue_SGs_24[6,:] )


println(  VS_revenue_VSG_1[1,:]+ VS_revenue_VSG_1[2,:])






#--------------------------------
#  Dualization   optimal obj of dual model ==1.06802707e+04
#--------------------------------
  
dual_model = dualize(model)

set_optimizer(dual_model , Gurobi.Optimizer)
optimize!(dual_model)

all_vars = all_variables(dual_model)




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







Kˢᵗ=[125, 125, 125, 125, 125, 125]              #  Startup cost of SGs                     SGs, buses:2,3,4,5,27,30
Kˢʰ=[50, 50, 50, 50, 50, 50]            #  Shutdown cost of SGs                    SGs, buses:2,3,4,5,27,30
Oᵐ=[10.47, 10.47, 10.47, 10.47, 10.47, 10.47]                   #  Marginal generation cost of SGs   in    SGs, buses:2,3,4,5,27,30
Oⁿˡ=Oᵐ*10 