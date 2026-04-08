import Pkg
using JuMP,Gurobi, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles,Plots,MAT
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


#------loop--------

ratio_list=[0, 0.2, 0.4, 0.6, 0.8, 1]
costs_ratio_voltage_support=zeros(6)
UC_list=zeros(6)
Tao_list_23=zeros(6)
Tao_list_24=zeros(6)
Tao_list_23_real=zeros(6)
Tao_list_24_real=zeros(6)
P_1_curtail_ave=zeros(6)
P_23_curtail_ave=zeros(6)
P_24_curtail_ave=zeros(6)


UC_schedule=zeros(6,24)

a=0
for ratio in ratio_list
a=a+1

#-----------------
# Define model
#-----------------

model= Model()                     

@variable(model, yˢᴳ²[1:T],Bin)                         # status of SGs, buses:2,3,4,5,27,30.    
@variable(model, yˢᴳ³[1:T],Bin)            
@variable(model, yˢᴳ⁴[1:T],Bin)           
@variable(model, yˢᴳ⁵[1:T],Bin)           
@variable(model, yˢᴳ²⁷[1:T],Bin)            
@variable(model, yˢᴳ³⁰[1:T],Bin)  

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

    @constraint(model, (Pᴳ[2,t] /(Sᵐᵃˣ_gc[1]))^2 +(Qᴳ[2,t] /(Sᵐᵃˣ_gc[1]))^2 <= 1)   # bounds for apparent power of SGs
    @constraint(model, (Pᴳ[3,t] /(Sᵐᵃˣ_gc[2]))^2 +(Qᴳ[3,t] /(Sᵐᵃˣ_gc[2]))^2 <= 1)
    @constraint(model, (Pᴳ[4,t] /(Sᵐᵃˣ_gc[3]))^2 +(Qᴳ[4,t] /(Sᵐᵃˣ_gc[3]))^2 <= 1)
    @constraint(model, (Pᴳ[5,t] /(Sᵐᵃˣ_gc[4]))^2 +(Qᴳ[5,t] /(Sᵐᵃˣ_gc[4]))^2 <= 1)
    @constraint(model, (Pᴳ[27,t] /(Sᵐᵃˣ_gc[5]))^2 +(Qᴳ[27,t] /(Sᵐᵃˣ_gc[5]))^2 <= 1)
    @constraint(model, (Pᴳ[30,t] /(Sᵐᵃˣ_gc[6]))^2 +(Qᴳ[30,t] /(Sᵐᵃˣ_gc[6]))^2 <= 1)

    @constraint(model, (Pᴳ[1,t] /(Sᵐᵃˣ_gv[1]))^2 +(Qᴳ[1,t] /(Sᵐᵃˣ_gv[1]))^2 <= 1)   # bounds for apparent power of VSGs

    @constraint(model, (Pᴳ[23,t] /(Sᵐᵃˣ_c[1]))^2 +(Qᴳ[23,t] /(Sᵐᵃˣ_c[1]))^2 <= 1)    # bounds for apparent power of IBGs 
    @constraint(model, (Pᴳ[24,t] /(Sᵐᵃˣ_c[2]))^2 +(Qᴳ[24,t] /(Sᵐᵃˣ_c[2]))^2 <= 1)

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

#@variable(model, c[1:30,1:30,1:T])
#@variable(model, s[1:30,1:30,1:T])
for i in 1:bus_num
    for t in 1:T
        #@constraint(model, Pᴳ[i,t] - Pᴰ[i,t] ==  G_matrix[i,i]*c[i,i,t] + 
            #sum( G_matrix[i,j] *c[i,j,t] - B_matrix[i,j] *s[i,j,t] for j in bus_connection[i] ) )    # p_i^g-p_i^d == G_{i,i}*c_{i,i}+...   eq. 50(a)
        #@constraint(model, Qᴳ[i,t] - Qᴰ[i,t] == -B_matrix[i,i]*c[i,i,t] - 
            #sum( B_matrix[i,j] *c[i,j,t] + G_matrix[i,j] *s[i,j,t] for j in bus_connection[i] ))    # q_i^g-q_i^d == -B_{i,i}*c_{i,i}+...  eq. 50(b)
        #@constraint(model, Pᴳ[i,t] - Pᴰ[i,t] ==  G_matrix[i,i]*c[i,i,t] + 
            #sum( Pˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # p_i^g-p_i^d == G_{i,i}*c_{i,i}+...   eq. 50(a)
        #@constraint(model, Qᴳ[i,t] - Qᴰ[i,t] == -B_matrix[i,i]*c[i,i,t] + 
            #sum( Qˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # q_i^g-q_i^d == -B_{i,i}*c_{i,i}+...  eq. 50(b)

        @constraint(model, Pᴳ[i,t] - Pᴰ[i,t] ==  sum( Pˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # p_i^g-p_i^d == G_{i,i}*c_{i,i}+...   eq. 50(a)
        @constraint(model, Qᴳ[i,t] - Qᴰ[i,t] ==  sum( Qˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # q_i^g-q_i^d == -B_{i,i}*c_{i,i}+...  eq. 50(b)

        #@constraint(model, c[data_bus[i,1],data_bus[i,1],t] <= (data_bus[i,12])^2)    # c_{i,i} <=  V_{i,max}^2   eq. 50(c)
        #@constraint(model, c[data_bus[i,1],data_bus[i,1],t] >= (data_bus[i,13])^2)    # V_{i,min}^2  <= c_{i,i}   eq. 50(c)
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

            #@constraint(model, c[i,j,t] - c[j,i,t] == 0)                              # c_{i,j} == c_{j,i}   eq. 50(f)
            #@constraint(model, s[i,j,t] + s[j,i,t] == 0)                              # s_{i,j} == -s_{j,i}   eq. 50(f)
            #@constraint(model, c[i,j,t]^2 + s[i,j,t]^2 <= c[i,i,t] * c[j,j,t])        # c_{i,j}^2+s_{i,j}^2 <= c_{i,i}*c_{j,j}   eq. 50(g)       
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

@variable(model, ηₘ_1[1:15,1:T],Bin)        # product of each pair of SG 
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

    @constraint(model, P_23_eq[t]== Pᴳ[23,t]+  ( sh_1- ( K_c_gc[1,1] *yˢᴳ²[t] +K_c_gc[1,2] *yˢᴳ³[t] +K_c_gc[1,3] *yˢᴳ⁴[t] 
    +K_c_gc[1,4] *yˢᴳ⁵[t] +K_c_gc[1,5] *yˢᴳ²⁷[t] +K_c_gc[1,6] *yˢᴳ³⁰[t]    +α_VSG*K_c_gv[1,1] 
    + sum(K_c_m[1,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[1,16:end] .*ηₘ_2[:,t] ) ) )    *Pᴳ[24,t])

    @constraint(model, P_24_eq[t]== Pᴳ[24,t]+  ( sh_2- ( K_c_gc[2,1] *yˢᴳ²[t] +K_c_gc[2,2] *yˢᴳ³[t] +K_c_gc[2,3] *yˢᴳ⁴[t] 
    +K_c_gc[2,4] *yˢᴳ⁵[t] +K_c_gc[2,5] *yˢᴳ²⁷[t] +K_c_gc[2,6] *yˢᴳ³⁰[t]    +α_VSG*K_c_gv[2,1]
    + sum(K_c_m[2,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[2,16:end] .*ηₘ_2[:,t] ) ) )    *Pᴳ[23,t])

    @constraint(model, Q_23_eq[t]== Qᴳ[23,t]+  ( sh_1- ( K_c_gc[1,1] *yˢᴳ²[t] +K_c_gc[1,2] *yˢᴳ³[t] +K_c_gc[1,3] *yˢᴳ⁴[t] 
    +K_c_gc[1,4] *yˢᴳ⁵[t] +K_c_gc[1,5] *yˢᴳ²⁷[t] +K_c_gc[1,6] *yˢᴳ³⁰[t]    +α_VSG*K_c_gv[1,1] 
    + sum(K_c_m[1,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[1,16:end] .*ηₘ_2[:,t] ) ) )    *Qᴳ[24,t])

    @constraint(model, Q_24_eq[t]== Qᴳ[24,t]+  ( sh_2- ( K_c_gc[2,1] *yˢᴳ²[t] +K_c_gc[2,2] *yˢᴳ³[t] +K_c_gc[2,3] *yˢᴳ⁴[t] 
    +K_c_gc[2,4] *yˢᴳ⁵[t] +K_c_gc[2,5] *yˢᴳ²⁷[t] +K_c_gc[2,6] *yˢᴳ³⁰[t]    +α_VSG*K_c_gv[2,1]
    + sum(K_c_m[2,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[2,16:end] .*ηₘ_2[:,t] ) ) )    *Qᴳ[23,t])

    @constraint(model, Γ_23[t]== 1/2 * ( K_c_gc[3,1] *yˢᴳ²[t] +K_c_gc[3,2] *yˢᴳ³[t] +K_c_gc[3,3] *yˢᴳ⁴[t] 
    +K_c_gc[3,4] *yˢᴳ⁵[t] +K_c_gc[3,5] *yˢᴳ²⁷[t] +K_c_gc[3,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[3,16:end] .*ηₘ_2[:,t]   ) )

    @constraint(model, Γ_24[t]== 1/2 * ( K_c_gc[4,1] *yˢᴳ²[t] +K_c_gc[4,2] *yˢᴳ³[t] +K_c_gc[4,3] *yˢᴳ⁴[t] 
    +K_c_gc[4,4] *yˢᴳ⁵[t] +K_c_gc[4,5] *yˢᴳ²⁷[t] +K_c_gc[4,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[4,1]
    + sum(K_c_m[4,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[4,16:end] .*ηₘ_2[:,t]   ) )

    @constraint(model, [Q_23_eq[t] + Γ_23[t]*100, P_23_eq[t], Q_23_eq[t]] in SecondOrderCone()  )      
    @constraint(model, [Q_24_eq[t] + Γ_24[t]*100, P_24_eq[t], Q_24_eq[t]] in SecondOrderCone()  )

end



#--------------------------------
#  Model solving
#--------------------------------

cost_onoff_Primal=sum(Cᵁ²)+sum(Cᴰ²)+sum(Cᵁ³)+sum(Cᴰ³)+sum(Cᵁ⁴)+sum(Cᴰ⁴)+sum(Cᵁ⁵)+sum(Cᴰ⁵)+sum(Cᵁ²⁷)+sum(Cᴰ²⁷)+sum(Cᵁ³⁰)+sum(Cᴰ³⁰)  
cost_nl_Primal=sum(Oⁿˡ[1].*(yˢᴳ²))+sum(Oⁿˡ[2].*(yˢᴳ³))+sum(Oⁿˡ[3].*(yˢᴳ⁴))+sum(Oⁿˡ[4].*(yˢᴳ⁵))+sum(Oⁿˡ[5].*(yˢᴳ²⁷))+sum(Oⁿˡ[6].*(yˢᴳ³⁰))    
cost_gene_Primal=sum(Oᵐ[1].*Pᴳ[2,:])+ sum(Oᵐ[2].*Pᴳ[3,:])+sum(Oᵐ[3].*Pᴳ[4,:])+sum(Oᵐ[4].*Pᴳ[5,:])+sum(Oᵐ[5].*Pᴳ[27,:])+sum(Oᵐ[6].*Pᴳ[30,:])   
obj=cost_onoff_Primal +cost_nl_Primal +cost_gene_Primal 


@objective(model, Min, obj)  # single-level objective function
#-------Solve and Output Results
set_optimizer(model , Gurobi.Optimizer)
set_optimizer_attribute(model, "QCPDual", 1)
optimize!(model)

obj=JuMP.value.(obj)

yˢᴳ²=JuMP.value.(yˢᴳ²)
yˢᴳ³=JuMP.value.(yˢᴳ³)
yˢᴳ⁴=JuMP.value.(yˢᴳ⁴)
yˢᴳ⁵=JuMP.value.(yˢᴳ⁵)
yˢᴳ²⁷=JuMP.value.(yˢᴳ²⁷)
yˢᴳ³⁰=JuMP.value.(yˢᴳ³⁰)

zᵟ_23_1_list=zeros(T)
zᵟ_24_1_list=zeros(T)
for t in 1:T
    UC_schedule[:,t]=[yˢᴳ²[t];yˢᴳ³[t];yˢᴳ⁴[t];yˢᴳ⁵[t];yˢᴳ²⁷[t];yˢᴳ³⁰[t]]
    zᵟ_23_1, zᵟ_24_1 =  calculate_Tao(α_VSG,Gᵥ,IBG,UC_schedule[:,t]) 
    zᵟ_23_1_list[t]=zᵟ_23_1
    zᵟ_24_1_list[t]=zᵟ_24_1
end

Γ_23=JuMP.value.(Γ_23)
Γ_24=JuMP.value.(Γ_24)

P_1_curtail_ave[a]= sum(Pⱽˢᴳₘₐₓ[1] .-JuMP.value.(Pᴳ[1,:]))/T
P_23_curtail_ave[a]= sum(Pᴵᴮᴳₘₐₓ[1] .-JuMP.value.(Pᴳ[23,:]))/T
P_24_curtail_ave[a]= sum(Pᴵᴮᴳₘₐₓ[2] .-JuMP.value.(Pᴳ[24,:]))/T

UC_list[a]=sum(yˢᴳ²) +sum(yˢᴳ³) +sum(yˢᴳ⁴) +sum(yˢᴳ⁵) +sum(yˢᴳ²⁷) +sum(yˢᴳ³⁰)
costs_ratio_voltage_support[a]=obj

Tao_list_23[a]=2*sum(Γ_23)/T
Tao_list_24[a]=2*sum(Γ_24)/T

Tao_list_23_real[a]=sum(zᵟ_23_1_list)/T
Tao_list_24_real[a]=sum(zᵟ_24_1_list)/T

end

Pᴳ=JuMP.value.(Pᴳ)
Qᴳ=JuMP.value.(Qᴳ)

plot(costs_ratio_voltage_support)
plot(Tao_list_23)
plot!(Tao_list_24)
plot(Tao_list_23_real)
plot!(Tao_list_24_real)

plot(UC_list)
plot(P_1_curtail_ave)
plot(P_23_curtail_ave)
plot!(P_24_curtail_ave)






matwrite("costs_full.mat", Dict("costs_full" => costs_full))

