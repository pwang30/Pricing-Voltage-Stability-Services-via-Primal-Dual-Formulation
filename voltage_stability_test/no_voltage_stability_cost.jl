# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Pricing voltage stability
# 4.July.2025

import Pkg
using JuMP,Gurobi, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles,Plots,MAT


#-----------------------------------Define Parameters for Optimization-----------------------------------  
      
Load_total=[18.42,17.95,18.29,18.51,18.13,17.88,19.46,21.97,23.17,23.87,
23.91,23.77,23.80,23.82,24.23,23.79,26.01,26.91,25.26,23.69,22.12,20.04,18.17,18.01]*10^3/50   # (MW)
T=length(Load_total)
#plot(Load_total)

costs_no_VS=zeros(1,13)
α=0.4
for i in 1:13
    

VSG₁=α*400*ones(T,1)      # VSG capacity
IBG₂₃=α*100*ones(T,1)     # IBG capacity
IBG₂₄=α*100*ones(T,1)     # IBG capacity


Pˢᴳₘₐₓ=[6.584, 5.760, 3.781, 3.335, 3.252, 2.880]*15            # Max generation of SGs                    SGs, buses:2,3,4,5,27,30
Pˢᴳₘᵢₙ=[3.292, 2.880, 1.512, 0.667, 0.650, 0.288]*15            #  Min generation of SGs                   SGs, buses:2,3,4,5,27,30
Rₘₐₓ=[1.317, 1.152, 1.512, 1.334, 1.951, 1.728]*15              #  Ramp limits of SGs                      SGs, buses:2,3,4,5,27,30
Kˢᵗ=[200, 125, 92.5, 72, 55, 31]*10                           #  Startup cost of SGs                     SGs, buses:2,3,4,5,27,30
Kˢʰ=[50, 28.5, 18.5, 14.4, 12, 10]*10                        #  Shutdown cost of SGs                    SGs, buses:2,3,4,5,27,30
Oᵐ=[6.20, 7.10, 10.47, 12.28, 13.53, 15.36]                        #  Marginal generation cost of SG 1  in    SGs, buses:2,3,4,5,27,30
Oⁿˡ=[17.431, 15.005, 13.755, 10.930, 9.900, 8.570]*10         #  No-load cost of SGs                     SGs, buses:2,3,4,5,27,30
P_g₀=[5.268 4.608 3.025 2.668 2.602 0]*15                       #  Initial generation (t=0) of SGs         SGs, buses:2,3,4,5,27,30
yˢᴳ₀=[1 1 1 1 1 0]  



#-----------------------------------Define Primal-Dual Model-----------------------------------
model= Model()
#-------Define Primal Variales
 
@variable(model, Pˢᴳ²[1:T])        # generation of SGs , buses:2,3,4,5,27,30.  
@variable(model, Pˢᴳ³[1:T])                    
@variable(model, Pˢᴳ⁴[1:T])               
@variable(model, Pˢᴳ⁵[1:T])               
@variable(model, Pˢᴳ²⁷[1:T])                
@variable(model, Pˢᴳ³⁰[1:T])                 

@variable(model, Pⱽˢᴳ¹[1:T]>=0)         # wind generation, buses:1, 23, 24 , dual variables: ζᵐⁱⁿₜ  
@variable(model, Pᴵᴮᴳ²³[1:T]>=0)              
@variable(model, Pᴵᴮᴳ²⁴[1:T]>=0)  

@variable(model, yˢᴳ²[1:T],Bin)      # status of SGs, buses:2,3,4,5,27,30.    
@variable(model, yˢᴳ³[1:T],Bin)            
@variable(model, yˢᴳ⁴[1:T],Bin)           
@variable(model, yˢᴳ⁵[1:T],Bin)           
@variable(model, yˢᴳ²⁷[1:T],Bin)            
@variable(model, yˢᴳ³⁰[1:T],Bin)                          


@variable(model, Cᵁ²[1:T]>=0)                 # startup costs and shutdown costs for SGs , dual variables: ρˢᵗₜ , ρˢʰₜ                
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




#-------Define Primal Constraints
#-----Voltage Stability Constraints

@constraint(model, Pˢᴳ² +Pˢᴳ³ +Pˢᴳ⁴ +Pˢᴳ⁵ +Pˢᴳ²⁷ +Pˢᴳ³⁰ +Pⱽˢᴳ¹ +Pᴵᴮᴳ²³ +Pᴵᴮᴳ²⁴ ==Load_total)     # power balance , dual variable: λᴱₜ

@constraint(model, Pˢᴳ².<=yˢᴳ²*Pˢᴳₘₐₓ[1])           # bounds for the output of SGs with UC , dual variables: μᵐⁱⁿₜ , μᵐᵃˣₜ
@constraint(model, yˢᴳ²*Pˢᴳₘᵢₙ[1].<=Pˢᴳ²)       
@constraint(model, Pˢᴳ³.<=yˢᴳ³*Pˢᴳₘₐₓ[2])       
@constraint(model, yˢᴳ³*Pˢᴳₘᵢₙ[2].<=Pˢᴳ³)         
@constraint(model, Pˢᴳ⁴.<=yˢᴳ⁴*Pˢᴳₘₐₓ[3])       
@constraint(model, yˢᴳ⁴*Pˢᴳₘᵢₙ[3].<=Pˢᴳ⁴)    
@constraint(model, Pˢᴳ⁵.<=yˢᴳ⁵*Pˢᴳₘₐₓ[4])       
@constraint(model, yˢᴳ⁵*Pˢᴳₘᵢₙ[4].<=Pˢᴳ⁵)
@constraint(model, Pˢᴳ²⁷.<=yˢᴳ²⁷*Pˢᴳₘₐₓ[5])       
@constraint(model, yˢᴳ²⁷*Pˢᴳₘᵢₙ[5].<=Pˢᴳ²⁷)
@constraint(model, Pˢᴳ³⁰.<=yˢᴳ³⁰*Pˢᴳₘₐₓ[6])       
@constraint(model, yˢᴳ³⁰*Pˢᴳₘᵢₙ[6].<=Pˢᴳ³⁰)


@constraint(model, Pˢᴳ²[1]-P_g₀[1]<=Rₘₐₓ[1])        # bounds for the ramp of SGs , dual variables: πʳᵈₜ , πʳᵘₜ
@constraint(model, -Rₘₐₓ[1]<=Pˢᴳ²[1]-P_g₀[1])  
@constraint(model, Pˢᴳ³[1]-P_g₀[2]<=Rₘₐₓ[2])
@constraint(model, -Rₘₐₓ[2]<=Pˢᴳ³[1]-P_g₀[2])
@constraint(model, Pˢᴳ⁴[1]-P_g₀[3]<=Rₘₐₓ[3])
@constraint(model, -Rₘₐₓ[3]<=Pˢᴳ⁴[1]-P_g₀[3])
@constraint(model, Pˢᴳ⁵[1]-P_g₀[4]<=Rₘₐₓ[4])
@constraint(model, -Rₘₐₓ[4]<=Pˢᴳ⁵[1]-P_g₀[4])
@constraint(model, Pˢᴳ²⁷[1]-P_g₀[5]<=Rₘₐₓ[5])
@constraint(model, -Rₘₐₓ[5]<=Pˢᴳ²⁷[1]-P_g₀[5])
@constraint(model, Pˢᴳ³⁰[1]-P_g₀[6]<=Rₘₐₓ[6])
@constraint(model, -Rₘₐₓ[6]<=Pˢᴳ³⁰[1]-P_g₀[6])

for t in 2:T                                                   
    @constraint(model, Pˢᴳ²[t]-Pˢᴳ²[t-1]<=Rₘₐₓ[1])        
    @constraint(model, -Rₘₐₓ[1]<=Pˢᴳ²[t]-Pˢᴳ²[t-1])  

    @constraint(model, Pˢᴳ³[t]-Pˢᴳ³[t-1]<=Rₘₐₓ[2])
    @constraint(model, -Rₘₐₓ[2]<=Pˢᴳ³[t]-Pˢᴳ³[t-1])  
    
    @constraint(model, Pˢᴳ⁴[t]-Pˢᴳ⁴[t-1]<=Rₘₐₓ[3])P
    @constraint(model, -Rₘₐₓ[3]<=Pˢᴳ⁴[t]-Pˢᴳ⁴[t-1])

    @constraint(model, Pˢᴳ⁵[t]-Pˢᴳ⁵[t-1]<=Rₘₐₓ[4])
    @constraint(model, -Rₘₐₓ[4]<=Pˢᴳ⁵[t]-Pˢᴳ⁵[t-1])

    @constraint(model, Pˢᴳ²⁷[t]-Pˢᴳ²⁷[t-1]<=Rₘₐₓ[5])
    @constraint(model, -Rₘₐₓ[5]<=Pˢᴳ²⁷[t]-Pˢᴳ²⁷[t-1])

    @constraint(model, Pˢᴳ³⁰[t]-Pˢᴳ³⁰[t-1]<=Rₘₐₓ[6])
    @constraint(model, -Rₘₐₓ[6]<=Pˢᴳ³⁰[t]-Pˢᴳ³⁰[t-1])
end

@constraint(model, Cᵁ²[1]>=(yˢᴳ²[1]-yˢᴳ₀[1])*Kˢᵗ[1])        # startup costs and shutdown costs for SGs , dual variables: σˢᵗₜ , σˢʰₜ
@constraint(model, Cᴰ²[1]>=(yˢᴳ₀[1]-yˢᴳ²[1])*Kˢʰ[1])  
@constraint(model, Cᵁ³[1]>=(yˢᴳ³[1]-yˢᴳ₀[2])*Kˢᵗ[2])
@constraint(model, Cᴰ³[1]>=(yˢᴳ₀[2]-yˢᴳ³[1])*Kˢʰ[2])
@constraint(model, Cᵁ⁴[1]>=(yˢᴳ⁴[1]-yˢᴳ₀[3])*Kˢᵗ[3])
@constraint(model, Cᴰ⁴[1]>=(yˢᴳ₀[3]-yˢᴳ⁴[1])*Kˢʰ[3])
@constraint(model, Cᵁ⁵[1]>=(yˢᴳ⁵[1]-yˢᴳ₀[4])*Kˢᵗ[4])
@constraint(model, Cᴰ⁵[1]>=(yˢᴳ₀[4]-yˢᴳ⁵[1])*Kˢʰ[4])
@constraint(model, Cᵁ²⁷[1]>=(yˢᴳ²⁷[1]-yˢᴳ₀[5])*Kˢᵗ[5])
@constraint(model, Cᴰ²⁷[1]>=(yˢᴳ₀[5]-yˢᴳ²⁷[1])*Kˢʰ[5])
@constraint(model, Cᵁ³⁰[1]>=(yˢᴳ³⁰[1]-yˢᴳ₀[6])*Kˢᵗ[6])
@constraint(model, Cᴰ³⁰[1]>=(yˢᴳ₀[6]-yˢᴳ³⁰[1])*Kˢʰ[6])

for t in 2:T
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
    
for t in 1:T
    @constraint(model, Pⱽˢᴳ¹[t] <= VSG₁[t])        # wind power limit  , dual variable: ζᵐᵃˣₜ
    @constraint(model, Pᴵᴮᴳ²³[t]<= IBG₂₃[t])       
    @constraint(model, Pᴵᴮᴳ²⁴[t]<= IBG₂₄[t])       
                
end




#-------Define Objective Functions 
#-Primal obj
cost_onoff_Primal=sum(Cᵁ²)+sum(Cᴰ²)+sum(Cᵁ³)+sum(Cᴰ³)+sum(Cᵁ⁴)+sum(Cᴰ⁴)+sum(Cᵁ⁵)+sum(Cᴰ⁵)+sum(Cᵁ²⁷)+sum(Cᴰ²⁷)+sum(Cᵁ³⁰)+sum(Cᴰ³⁰)  
cost_nl_Primal=sum(Oⁿˡ[1].*(yˢᴳ²))+sum(Oⁿˡ[2].*(yˢᴳ³))+sum(Oⁿˡ[3].*(yˢᴳ⁴))+sum(Oⁿˡ[4].*(yˢᴳ⁵))+sum(Oⁿˡ[5].*(yˢᴳ²⁷))+sum(Oⁿˡ[6].*(yˢᴳ³⁰))    
cost_gene_Primal=sum(Oᵐ[1].*Pˢᴳ²)+ sum(Oᵐ[2].*Pˢᴳ³)+sum(Oᵐ[3].*Pˢᴳ⁴)+sum(Oᵐ[4].*Pˢᴳ⁵)+sum(Oᵐ[5].*Pˢᴳ²⁷)+sum(Oᵐ[6].*Pˢᴳ³⁰)   

obj_Primal=cost_onoff_Primal +cost_nl_Primal +cost_gene_Primal 



@objective(model, Min, obj_Primal)  # single-level objective function
#-------Solve and Output Results
set_optimizer(model , Gurobi.Optimizer)
optimize!(model)

obj_Primal=objective_value(model)

costs_no_VS[i]=obj_Primal
α +=0.05

end

plot(costs_no_VS', label="No Voltage Stability Cost")
plot!(costs_full', label="Full Reactive Power Support")
plot!(costs', label="50% Reactive Power Support")