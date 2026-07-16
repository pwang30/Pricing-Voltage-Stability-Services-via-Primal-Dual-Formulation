# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Pricing voltage stability via Primal-Dual formulation
# 12.Apr.2026



function pricing_VS_dispatchable(K_c_gc, K_c_gv, K_c_m,
                        T,Pᴰ,Qᴰ,Pˢᴳₘₐₓ,Pˢᴳₘᵢₙ,Qˢᴳₘₐₓ,Qˢᴳₘᵢₙ,Pⱽˢᴳₘₐₓ,Pⱽˢᴳₘᵢₙ,Qⱽˢᴳₘₐₓ,Qⱽˢᴳₘᵢₙ,Pᴵᴮᴳₘₐₓ,Pᴵᴮᴳₘᵢₙ,Qᴵᴮᴳₘₐₓ,Qᴵᴮᴳₘᵢₙ,
                        Sᵐᵃˣ_gc,Sᵐᵃˣ_gv,Sᵐᵃˣ_c,S_Base,α_VSG,α_IBG,ratio,
                        cˢᵗ,cᵐ,cⁿˡ,
                        yˢᴳ²_0,yˢᴳ³_0,yˢᴳ⁴_0,yˢᴳ⁵_0,yˢᴳ²⁷_0,yˢᴳ³⁰_0)


    function get_dual_var(model, name)
        ref = variable_by_name(model, name)
        ref === nothing && error("none: $name")
        return value(ref)
    end

#-----------------
# Define model for pricing
#-----------------

model= Model()                     

@variable(model, yˢᴳ²[1:T]>=0)          
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

@variable(model, Cᵁ²[1:T] >=0)                                                
@variable(model, Cᵁ³[1:T] >=0)                            
@variable(model, Cᵁ⁴[1:T] >=0)                                
@variable(model, Cᵁ⁵[1:T] >=0)                                
@variable(model, Cᵁ²⁷[1:T] >=0)                               
@variable(model, Cᵁ³⁰[1:T] >=0)   

@variable(model, Pᴳ[1:9,1:T]) 
@variable(model, Qᴳ[1:9,1:T])   

@constraint(model, Cᵁ²[1]>=(yˢᴳ²[1]-yˢᴳ²_0)*cˢᵗ[1])    
@constraint(model, Cᵁ³[1]>=(yˢᴳ³[1]-yˢᴳ³_0)*cˢᵗ[2])    
@constraint(model, Cᵁ⁴[1]>=(yˢᴳ⁴[1]-yˢᴳ⁴_0)*cˢᵗ[3])     
@constraint(model, Cᵁ⁵[1]>=(yˢᴳ⁵[1]-yˢᴳ⁵_0)*cˢᵗ[4])    
@constraint(model, Cᵁ²⁷[1]>=(yˢᴳ²⁷[1]-yˢᴳ²⁷_0)*cˢᵗ[5])  
@constraint(model, Cᵁ³⁰[1]>=(yˢᴳ³⁰[1]-yˢᴳ³⁰_0)*cˢᵗ[6])  


for t in 2:T

    @constraint(model, Cᵁ²[t]>=(yˢᴳ²[t]-yˢᴳ²[t-1])*cˢᵗ[1])     
    @constraint(model, Cᵁ³[t]>=(yˢᴳ³[t]-yˢᴳ³[t-1])*cˢᵗ[2])    
    @constraint(model, Cᵁ⁴[t]>=(yˢᴳ⁴[t]-yˢᴳ⁴[t-1])*cˢᵗ[3])     
    @constraint(model, Cᵁ⁵[t]>=(yˢᴳ⁵[t]-yˢᴳ⁵[t-1])*cˢᵗ[4])     
    @constraint(model, Cᵁ²⁷[t]>=(yˢᴳ²⁷[t]-yˢᴳ²⁷[t-1])*cˢᵗ[5])  
    @constraint(model, Cᵁ³⁰[t]>=(yˢᴳ³⁰[t]-yˢᴳ³⁰[t-1])*cˢᵗ[6]) 

end 

for t in 1:T                                           
    @constraint(model, Pᴳ[1,t] <= yˢᴳ²[t]*Pˢᴳₘₐₓ[1])  
    @constraint(model, yˢᴳ²[t]*Pˢᴳₘᵢₙ[1] <= Pᴳ[1,t])     
    @constraint(model, Pᴳ[2,t] <= yˢᴳ³[t]*Pˢᴳₘₐₓ[2])  
    @constraint(model, yˢᴳ³[t]*Pˢᴳₘᵢₙ[2] <= Pᴳ[2,t])    
    @constraint(model, Pᴳ[3,t] <= yˢᴳ⁴[t]*Pˢᴳₘₐₓ[3])    
    @constraint(model, yˢᴳ⁴[t]*Pˢᴳₘᵢₙ[3] <= Pᴳ[3,t])     
    @constraint(model, Pᴳ[4,t] <= yˢᴳ⁵[t]*Pˢᴳₘₐₓ[4])     
    @constraint(model, yˢᴳ⁵[t]*Pˢᴳₘᵢₙ[4] <= Pᴳ[4,t])    
    @constraint(model, Pᴳ[5,t] <= yˢᴳ²⁷[t]*Pˢᴳₘₐₓ[5])  
    @constraint(model, yˢᴳ²⁷[t]*Pˢᴳₘᵢₙ[5] <=Pᴳ[5,t])   
    @constraint(model, Pᴳ[6,t] <= yˢᴳ³⁰[t]*Pˢᴳₘₐₓ[6])   
    @constraint(model, yˢᴳ³⁰[t]*Pˢᴳₘᵢₙ[6] <= Pᴳ[6,t])   

    @constraint(model, Qᴳ[1,t] <= yˢᴳ²[t]*Qˢᴳₘₐₓ[1])  
    @constraint(model, yˢᴳ²[t]*Qˢᴳₘᵢₙ[1] <= Qᴳ[1,t])  
    @constraint(model, Qᴳ[2,t] <= yˢᴳ³[t]*Qˢᴳₘₐₓ[2])     
    @constraint(model, yˢᴳ³[t]*Qˢᴳₘᵢₙ[2] <= Qᴳ[2,t])          
    @constraint(model, Qᴳ[3,t] <= yˢᴳ⁴[t]*Qˢᴳₘₐₓ[3])     
    @constraint(model, yˢᴳ⁴[t]*Qˢᴳₘᵢₙ[3] <= Qᴳ[3,t])   
    @constraint(model, Qᴳ[4,t] <= yˢᴳ⁵[t]*Qˢᴳₘₐₓ[4])  
    @constraint(model, yˢᴳ⁵[t]*Qˢᴳₘᵢₙ[4] <= Qᴳ[4,t])    
    @constraint(model, Qᴳ[5,t] <= yˢᴳ²⁷[t]*Qˢᴳₘₐₓ[5]) 
    @constraint(model, yˢᴳ²⁷[t]*Qˢᴳₘᵢₙ[5] <=Qᴳ[5,t])   
    @constraint(model, Qᴳ[6,t] <= yˢᴳ³⁰[t]*Qˢᴳₘₐₓ[6]) 
    @constraint(model, yˢᴳ³⁰[t]*Qˢᴳₘᵢₙ[6] <= Qᴳ[6,t])  

    @constraint(model, Pᴳ[7,t] <= α_VSG[t]*Pⱽˢᴳₘₐₓ[1])                
    @constraint(model, Pⱽˢᴳₘᵢₙ[1] <= Pᴳ[7,t])       
    @constraint(model, Pᴳ[8,t] <= α_IBG[1,t]*Pᴵᴮᴳₘₐₓ[1])   
    @constraint(model, Pᴵᴮᴳₘᵢₙ[1] <= Pᴳ[8,t])         
    @constraint(model, Pᴳ[9,t] <= α_IBG[2,t]*Pᴵᴮᴳₘₐₓ[2])    
    @constraint(model, Pᴵᴮᴳₘᵢₙ[2] <= Pᴳ[9,t])        

    @constraint(model, Qᴳ[7,t] <= Qⱽˢᴳₘₐₓ[1])   
    @constraint(model, Qⱽˢᴳₘᵢₙ[1] <= Qᴳ[7,t])    
    @constraint(model, Qᴳ[8,t] <= ratio*Qᴵᴮᴳₘₐₓ[1])    
    @constraint(model, Qᴵᴮᴳₘᵢₙ[1] <= Qᴳ[8,t])    
    @constraint(model, Qᴳ[9,t] <= ratio*Qᴵᴮᴳₘₐₓ[2])  
    @constraint(model, Qᴵᴮᴳₘᵢₙ[2] <= Qᴳ[9,t])    

    @constraint(model, [ Sᵐᵃˣ_gc[1], Pᴳ[1,t], Qᴳ[1,t] ] in SecondOrderCone()  )  
    @constraint(model, [ Sᵐᵃˣ_gc[2], Pᴳ[2,t], Qᴳ[2,t] ] in SecondOrderCone()  ) 
    @constraint(model, [ Sᵐᵃˣ_gc[3], Pᴳ[3,t], Qᴳ[3,t] ] in SecondOrderCone()  )   
    @constraint(model, [ Sᵐᵃˣ_gc[4], Pᴳ[4,t], Qᴳ[4,t] ] in SecondOrderCone()  )  
    @constraint(model, [ Sᵐᵃˣ_gc[5], Pᴳ[5,t], Qᴳ[5,t] ] in SecondOrderCone()  )  
    @constraint(model, [ Sᵐᵃˣ_gc[6], Pᴳ[6,t], Qᴳ[6,t] ] in SecondOrderCone()  )  
    @constraint(model, [ Sᵐᵃˣ_gv[1], Pᴳ[7,t], Qᴳ[7,t] ] in SecondOrderCone()  )    
    @constraint(model, [ Sᵐᵃˣ_c[1], Pᴳ[8,t], Qᴳ[8,t] ] in SecondOrderCone()  )    
    @constraint(model, [ Sᵐᵃˣ_c[2], Pᴳ[9,t], Qᴳ[9,t] ] in SecondOrderCone()  )    

    @constraint(model, sum(Qᴳ[:,t]) == Qᴰ[t] )    
    
end

@constraint(model, Power_balance_P[t in 1:T], sum(Pᴳ[:,t]) == Pᴰ[t])

#--------------------------Voltage stability constraints--------------------------

@variable(model, z_23[1:T])
@variable(model, z_24[1:T])

@variable(model, ηₘ_1[1:15,1:T] >= 0)       
@variable(model, ηₘ_2[1:6,1:T])     
for i in 1:15, j in 1:T
   @constraint(model, ηₘ_1[i,j] <=1)        
end

Voltage_constraint_23=Dict()
Voltage_constraint_24=Dict()

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

    @constraint(model, ηₘ_2[1,t]==yˢᴳ²[t]*α_VSG[t])
    @constraint(model, ηₘ_2[2,t]==yˢᴳ³[t]*α_VSG[t])
    @constraint(model, ηₘ_2[3,t]==yˢᴳ⁴[t]*α_VSG[t])
    @constraint(model, ηₘ_2[4,t]==yˢᴳ⁵[t]*α_VSG[t])
    @constraint(model, ηₘ_2[5,t]==yˢᴳ²⁷[t]*α_VSG[t])
    @constraint(model, ηₘ_2[6,t]==yˢᴳ³⁰[t]*α_VSG[t])

    @constraint(model, z_23[t] ==  ( ( K_c_gc[1, 1] *yˢᴳ²[t] +K_c_gc[1, 2] *yˢᴳ³[t] +K_c_gc[1, 3] *yˢᴳ⁴[t] 
    +K_c_gc[1, 4] *yˢᴳ⁵[t] +K_c_gc[1, 5] *yˢᴳ²⁷[t] +K_c_gc[1, 6] *yˢᴳ³⁰[t] )   +α_VSG[t]*K_c_gv[1, 1]
    + sum(K_c_m[1, 1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[1, 16:end] .*ηₘ_2[:,t]   ) ) )     

    @constraint(model, z_24[t] == ( ( K_c_gc[2, 1] *yˢᴳ²[t] +K_c_gc[2, 2] *yˢᴳ³[t] +K_c_gc[2, 3] *yˢᴳ⁴[t] 
    +K_c_gc[2, 4] *yˢᴳ⁵[t] +K_c_gc[2, 5] *yˢᴳ²⁷[t] +K_c_gc[2, 6] *yˢᴳ³⁰[t] )   +α_VSG[t]*K_c_gv[2, 1]
    + sum(K_c_m[2, 1:15] .*ηₘ_1[:,t] ) +  sum(K_c_m[2, 16:end] .*ηₘ_2[:,t]   ) ) )     

end

@constraint(model, Voltage_constraint_23[t in 1:T], [Qᴳ[8,t] + 1/2*z_23[t]*S_Base,  Pᴳ[8,t],  Qᴳ[8,t]] in SecondOrderCone())
@constraint(model, Voltage_constraint_24[t in 1:T], [Qᴳ[9,t] + 1/2*z_24[t]*S_Base,  Pᴳ[9,t],  Qᴳ[9,t]] in SecondOrderCone())

#--------------------------------
#  Model solving
#--------------------------------
@variable(model, obj_primal[1:T])
@variable(model, obj_primal_final)
for t in 1:T
    @constraint(model, obj_primal[t] == Cᵁ²[t]+Cᵁ³[t]+Cᵁ⁴[t]+Cᵁ⁵[t]+Cᵁ²⁷[t]+Cᵁ³⁰[t] 
                    + cⁿˡ[1]*yˢᴳ²[t] + cⁿˡ[2]*yˢᴳ³[t] + cⁿˡ[3]*yˢᴳ⁴[t] + cⁿˡ[4]*yˢᴳ⁵[t] + cⁿˡ[5]*yˢᴳ²⁷[t] + cⁿˡ[6]*yˢᴳ³⁰[t] 
                    + cᵐ[1]*Pᴳ[1,t] + cᵐ[2]*Pᴳ[2,t] + cᵐ[3]*Pᴳ[3,t] + cᵐ[4]*Pᴳ[4,t] + cᵐ[5]*Pᴳ[5,t] + cᵐ[6]*Pᴳ[6,t] )
end
@constraint(model, obj_primal_final == sum(obj_primal) )

@objective(model, Min,  obj_primal_final)  # single-level objective function

dual_model = dualize(model; dual_names = DualNames("dual_var_", "dual_con_"))
set_optimizer(dual_model, Gurobi.Optimizer)
optimize!(dual_model) 


#------------------------------
# Get prices
#------------------------------

mu=zeros(T,2)   
lambda_2=zeros(T,2)
energy_prices_P = zeros(T)  

for t in 1:T

    energy_prices_P[t] = get_dual_var(dual_model, "dual_var_Power_balance_P[$t]_1")   

    mu[t,1] = get_dual_var(dual_model, "dual_var_Voltage_constraint_23[$t]_1")  
    mu[t,2] = get_dual_var(dual_model, "dual_var_Voltage_constraint_24[$t]_1")   

    lambda_2[t,1] = get_dual_var(dual_model, "dual_var_Voltage_constraint_23[$t]_3")  
    lambda_2[t,2] = get_dual_var(dual_model, "dual_var_Voltage_constraint_24[$t]_3")   
end


#-----------------
# Dispatch units
#-----------------


model= Model()                     

@variable(model, yˢᴳ²[1:T] ,Bin)    
@variable(model, yˢᴳ³[1:T] ,Bin)            
@variable(model, yˢᴳ⁴[1:T] ,Bin)           
@variable(model, yˢᴳ⁵[1:T] ,Bin)           
@variable(model, yˢᴳ²⁷[1:T] ,Bin)            
@variable(model, yˢᴳ³⁰[1:T] ,Bin)  

@variable(model, Cᵁ²[1:T] >=0)                                                
@variable(model, Cᵁ³[1:T] >=0)                            
@variable(model, Cᵁ⁴[1:T] >=0)                                
@variable(model, Cᵁ⁵[1:T] >=0)                                
@variable(model, Cᵁ²⁷[1:T] >=0)                               
@variable(model, Cᵁ³⁰[1:T] >=0)                 

@variable(model, Pᴳ[1:9,1:T]) 
@variable(model, Qᴳ[1:9,1:T])  

@constraint(model, Cᵁ²[1]>=(yˢᴳ²[1]-yˢᴳ²_0)*cˢᵗ[1])     
@constraint(model, Cᵁ³[1]>=(yˢᴳ³[1]-yˢᴳ³_0)*cˢᵗ[2])     
@constraint(model, Cᵁ⁴[1]>=(yˢᴳ⁴[1]-yˢᴳ⁴_0)*cˢᵗ[3])     
@constraint(model, Cᵁ⁵[1]>=(yˢᴳ⁵[1]-yˢᴳ⁵_0)*cˢᵗ[4])     
@constraint(model, Cᵁ²⁷[1]>=(yˢᴳ²⁷[1]-yˢᴳ²⁷_0)*cˢᵗ[5])  
@constraint(model, Cᵁ³⁰[1]>=(yˢᴳ³⁰[1]-yˢᴳ³⁰_0)*cˢᵗ[6])  

for t in 2:T

    @constraint(model, Cᵁ²[t]>=(yˢᴳ²[t]-yˢᴳ²[t-1])*cˢᵗ[1])     
    @constraint(model, Cᵁ³[t]>=(yˢᴳ³[t]-yˢᴳ³[t-1])*cˢᵗ[2])     
    @constraint(model, Cᵁ⁴[t]>=(yˢᴳ⁴[t]-yˢᴳ⁴[t-1])*cˢᵗ[3])     
    @constraint(model, Cᵁ⁵[t]>=(yˢᴳ⁵[t]-yˢᴳ⁵[t-1])*cˢᵗ[4])     
    @constraint(model, Cᵁ²⁷[t]>=(yˢᴳ²⁷[t]-yˢᴳ²⁷[t-1])*cˢᵗ[5]) 
    @constraint(model, Cᵁ³⁰[t]>=(yˢᴳ³⁰[t]-yˢᴳ³⁰[t-1])*cˢᵗ[6])  

end

for t in 1:T

    @constraint(model, Pᴳ[1,t] <= yˢᴳ²[t]*Pˢᴳₘₐₓ[1])        
    @constraint(model, yˢᴳ²[t]*Pˢᴳₘᵢₙ[1] <= Pᴳ[1,t])    
    @constraint(model, Pᴳ[2,t] <= yˢᴳ³[t]*Pˢᴳₘₐₓ[2])    
    @constraint(model, yˢᴳ³[t]*Pˢᴳₘᵢₙ[2] <= Pᴳ[2,t])     
    @constraint(model, Pᴳ[3,t] <= yˢᴳ⁴[t]*Pˢᴳₘₐₓ[3])     
    @constraint(model, yˢᴳ⁴[t]*Pˢᴳₘᵢₙ[3] <= Pᴳ[3,t])     
    @constraint(model, Pᴳ[4,t] <= yˢᴳ⁵[t]*Pˢᴳₘₐₓ[4])     
    @constraint(model, yˢᴳ⁵[t]*Pˢᴳₘᵢₙ[4] <= Pᴳ[4,t])     
    @constraint(model, Pᴳ[5,t] <= yˢᴳ²⁷[t]*Pˢᴳₘₐₓ[5])   
    @constraint(model, yˢᴳ²⁷[t]*Pˢᴳₘᵢₙ[5] <=Pᴳ[5,t])    
    @constraint(model, Pᴳ[6,t] <= yˢᴳ³⁰[t]*Pˢᴳₘₐₓ[6])   
    @constraint(model, yˢᴳ³⁰[t]*Pˢᴳₘᵢₙ[6] <= Pᴳ[6,t])   

    @constraint(model, Qᴳ[1,t] <= yˢᴳ²[t]*Qˢᴳₘₐₓ[1])   
    @constraint(model, yˢᴳ²[t]*Qˢᴳₘᵢₙ[1] <= Qᴳ[1,t])    
    @constraint(model, Qᴳ[2,t] <= yˢᴳ³[t]*Qˢᴳₘₐₓ[2])   
    @constraint(model, yˢᴳ³[t]*Qˢᴳₘᵢₙ[2] <= Qᴳ[2,t])      
    @constraint(model, Qᴳ[3,t] <= yˢᴳ⁴[t]*Qˢᴳₘₐₓ[3])   
    @constraint(model, yˢᴳ⁴[t]*Qˢᴳₘᵢₙ[3] <= Qᴳ[3,t])   
    @constraint(model, Qᴳ[4,t] <= yˢᴳ⁵[t]*Qˢᴳₘₐₓ[4])   
    @constraint(model, yˢᴳ⁵[t]*Qˢᴳₘᵢₙ[4] <= Qᴳ[4,t])    
    @constraint(model, Qᴳ[5,t] <= yˢᴳ²⁷[t]*Qˢᴳₘₐₓ[5])  
    @constraint(model, yˢᴳ²⁷[t]*Qˢᴳₘᵢₙ[5] <=Qᴳ[5,t])  
    @constraint(model, Qᴳ[6,t] <= yˢᴳ³⁰[t]*Qˢᴳₘₐₓ[6])  
    @constraint(model, yˢᴳ³⁰[t]*Qˢᴳₘᵢₙ[6] <= Qᴳ[6,t])  

    @constraint(model, Pᴳ[7,t] <= α_VSG[t]*Pⱽˢᴳₘₐₓ[1])              
    @constraint(model, 0 <= Pᴳ[7,t])       
    @constraint(model, Pᴳ[8,t] <= α_IBG[1,t]*Pᴵᴮᴳₘₐₓ[1])     
    @constraint(model, 0 <= Pᴳ[8,t])         
    @constraint(model, Pᴳ[9,t] <= α_IBG[2,t]*Pᴵᴮᴳₘₐₓ[2])     
    @constraint(model, 0 <= Pᴳ[9,t])        

    @constraint(model, Qᴳ[7,t] <= Qⱽˢᴳₘₐₓ[1])    
    @constraint(model, Qⱽˢᴳₘᵢₙ[1] <= Qᴳ[7,t])    
    @constraint(model, Qᴳ[8,t] <= ratio*Qᴵᴮᴳₘₐₓ[1])    
    @constraint(model, Qᴵᴮᴳₘᵢₙ[1] <= Qᴳ[8,t])    
    @constraint(model, Qᴳ[9,t] <= ratio*Qᴵᴮᴳₘₐₓ[2])    
    @constraint(model, Qᴵᴮᴳₘᵢₙ[2] <= Qᴳ[9,t])    

    @constraint(model, [ Sᵐᵃˣ_gc[1], Pᴳ[1,t], Qᴳ[1,t] ] in SecondOrderCone()  )   
    @constraint(model, [ Sᵐᵃˣ_gc[2], Pᴳ[2,t], Qᴳ[2,t] ] in SecondOrderCone()  )   
    @constraint(model, [ Sᵐᵃˣ_gc[3], Pᴳ[3,t], Qᴳ[3,t] ] in SecondOrderCone()  )   
    @constraint(model, [ Sᵐᵃˣ_gc[4], Pᴳ[4,t], Qᴳ[4,t] ] in SecondOrderCone()  )   
    @constraint(model, [ Sᵐᵃˣ_gc[5], Pᴳ[5,t], Qᴳ[5,t] ] in SecondOrderCone()  )  
    @constraint(model, [ Sᵐᵃˣ_gc[6], Pᴳ[6,t], Qᴳ[6,t] ] in SecondOrderCone()  )  
    @constraint(model, [ Sᵐᵃˣ_gv[1], Pᴳ[7,t], Qᴳ[7,t] ] in SecondOrderCone()  )    
    @constraint(model, [ Sᵐᵃˣ_c[1], Pᴳ[8,t], Qᴳ[8,t] ] in SecondOrderCone()  )    
    @constraint(model, [ Sᵐᵃˣ_c[2], Pᴳ[9,t], Qᴳ[9,t] ] in SecondOrderCone()  )    

    @constraint(model, sum(Pᴳ[:,t]) == Pᴰ[t]  )    
    @constraint(model, sum(Qᴳ[:,t]) == Qᴰ[t] )   

end


#--------------------------Voltage stability constraints--------------------------

@variable(model, z_23[1:T])
@variable(model, z_24[1:T])

@variable(model, ηₘ_1[1:15,1:T] >= 0)       
@variable(model, ηₘ_2[1:6,1:T])            

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

    @constraint(model, ηₘ_2[1,t]==yˢᴳ²[t]*α_VSG[t])        
    @constraint(model, ηₘ_2[2,t]==yˢᴳ³[t]*α_VSG[t])        
    @constraint(model, ηₘ_2[3,t]==yˢᴳ⁴[t]*α_VSG[t])        
    @constraint(model, ηₘ_2[4,t]==yˢᴳ⁵[t]*α_VSG[t])        
    @constraint(model, ηₘ_2[5,t]==yˢᴳ²⁷[t]*α_VSG[t])       
    @constraint(model, ηₘ_2[6,t]==yˢᴳ³⁰[t]*α_VSG[t])       

    @constraint(model, z_23[t] ==  ( ( K_c_gc[1, 1] *yˢᴳ²[t] +K_c_gc[1, 2] *yˢᴳ³[t] +K_c_gc[1, 3] *yˢᴳ⁴[t] 
    +K_c_gc[1, 4] *yˢᴳ⁵[t] +K_c_gc[1, 5] *yˢᴳ²⁷[t] +K_c_gc[1, 6] *yˢᴳ³⁰[t] )   +α_VSG[t]*K_c_gv[1, 1]
    + sum(K_c_m[1, 1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[1, 16:end] .*ηₘ_2[:,t]   ) ) )     

    @constraint(model, z_24[t] == ( ( K_c_gc[2, 1] *yˢᴳ²[t] +K_c_gc[2, 2] *yˢᴳ³[t] +K_c_gc[2, 3] *yˢᴳ⁴[t] 
    +K_c_gc[2, 4] *yˢᴳ⁵[t] +K_c_gc[2, 5] *yˢᴳ²⁷[t] +K_c_gc[2, 6] *yˢᴳ³⁰[t] )   +α_VSG[t]*K_c_gv[2, 1]
    + sum(K_c_m[2, 1:15] .*ηₘ_1[:,t] ) +  sum(K_c_m[2, 16:end] .*ηₘ_2[:,t]   ) ) )     

    @constraint(model, [Qᴳ[8,t] + 1/2*z_23[t]*S_Base,  Pᴳ[8,t],  Qᴳ[8,t]] in SecondOrderCone()  )     
    @constraint(model, [Qᴳ[9,t] + 1/2*z_24[t]*S_Base,  Pᴳ[9,t],  Qᴳ[9,t]] in SecondOrderCone()  )     

end

#--------------------------------
#  Model solving
#--------------------------------
@variable(model, obj_primal[1:T])
@variable(model, obj_primal_final)
for t in 1:T
    @constraint(model, obj_primal[t] == Cᵁ²[t]+Cᵁ³[t]+Cᵁ⁴[t]+Cᵁ⁵[t]+Cᵁ²⁷[t]+Cᵁ³⁰[t] 
                    + cⁿˡ[1]*yˢᴳ²[t] + cⁿˡ[2]*yˢᴳ³[t] + cⁿˡ[3]*yˢᴳ⁴[t] + cⁿˡ[4]*yˢᴳ⁵[t] + cⁿˡ[5]*yˢᴳ²⁷[t] + cⁿˡ[6]*yˢᴳ³⁰[t] 
                    + cᵐ[1]*Pᴳ[1,t] + cᵐ[2]*Pᴳ[2,t] + cᵐ[3]*Pᴳ[3,t] + cᵐ[4]*Pᴳ[4,t] + cᵐ[5]*Pᴳ[5,t] + cᵐ[6]*Pᴳ[6,t] )
end
@constraint(model, obj_primal_final == sum(obj_primal) )

@objective(model, Min,  obj_primal_final)  # single-level objective function
#-------Solve and Output Results
set_optimizer(model,  Gurobi.Optimizer)   
optimize!(model)

Cᵁ²=value.(Cᵁ²)
Cᵁ³=value.(Cᵁ³)
Cᵁ⁴=value.(Cᵁ⁴)
Cᵁ⁵=value.(Cᵁ⁵)
Cᵁ²⁷=value.(Cᵁ²⁷)
Cᵁ³⁰=value.(Cᵁ³⁰)

yˢᴳ²=JuMP.value.(yˢᴳ²)
yˢᴳ³=JuMP.value.(yˢᴳ³)
yˢᴳ⁴=JuMP.value.(yˢᴳ⁴)
yˢᴳ⁵=JuMP.value.(yˢᴳ⁵)
yˢᴳ²⁷=JuMP.value.(yˢᴳ²⁷)
yˢᴳ³⁰=JuMP.value.(yˢᴳ³⁰)

Pᴳ=JuMP.value.(Pᴳ)
Qᴳ=JuMP.value.(Qᴳ)

ηₘ_1=JuMP.value.(ηₘ_1)
ηₘ_2=JuMP.value.(ηₘ_2)

#------------------------------
# Analyze results
#------------------------------
# Output optimized results of decision variables and dual variables for constraints


energy_profit_VSG_1, energy_profitˢᴳ², energy_profitˢᴳ³, energy_profitˢᴳ⁴, energy_profitˢᴳ⁵, energy_profitˢᴳ²⁷, energy_profitˢᴳ³⁰, energy_profit_IBG_23, energy_profit_IBG_24,
       revenue_IBG_23_give_reactive, revenue_IBG_24_give_reactive,
       VS_revenue_VSG_1, VS_revenue_SGs = settlement_dispatchable( Pᴳ,Qᴳ,mu,lambda_2,energy_prices_P, 
                    yˢᴳ²,yˢᴳ³,yˢᴳ⁴,yˢᴳ⁵,yˢᴳ²⁷,yˢᴳ³⁰,
                    Cᵁ²,Cᵁ³,Cᵁ⁴,Cᵁ⁵,Cᵁ²⁷,Cᵁ³⁰,
                    cᵐ,cⁿˡ,α_VSG,Gᵥ,IBG,S_Base,
                    K_c_gc, K_c_m, ηₘ_1, ηₘ_2)

energy_profit_SG = zeros(6,T)
energy_profit_SG[1,:] = energy_profitˢᴳ²
energy_profit_SG[2,:] = energy_profitˢᴳ³
energy_profit_SG[3,:] = energy_profitˢᴳ⁴
energy_profit_SG[4,:] = energy_profitˢᴳ⁵
energy_profit_SG[5,:] = energy_profitˢᴳ²⁷
energy_profit_SG[6,:] = energy_profitˢᴳ³⁰

return energy_prices_P, mu[:,1], mu[:,2], mu[:,1] .-lambda_2[:,1], mu[:,2] .-lambda_2[:,2], energy_profit_SG, VS_revenue_SGs, energy_profit_VSG_1, energy_profit_IBG_23, energy_profit_IBG_24, VS_revenue_VSG_1, revenue_IBG_23_give_reactive, revenue_IBG_24_give_reactive

end
