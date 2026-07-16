# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Pricing voltage stability via Primal-Dual formulation
# 12.Apr.2026



function pricing_VS_restricted(K_c_gc, K_c_gv, K_c_m,
                        T,Pᴰ,Qᴰ,Pˢᴳₘₐₓ,Pˢᴳₘᵢₙ,Qˢᴳₘₐₓ,Qˢᴳₘᵢₙ,Pⱽˢᴳₘₐₓ,Pⱽˢᴳₘᵢₙ,Qⱽˢᴳₘₐₓ,Qⱽˢᴳₘᵢₙ,Pᴵᴮᴳₘₐₓ,Pᴵᴮᴳₘᵢₙ,Qᴵᴮᴳₘₐₓ,Qᴵᴮᴳₘᵢₙ,
                        Sᵐᵃˣ_gc,Sᵐᵃˣ_gv,Sᵐᵃˣ_c,S_Base,α_VSG,α_IBG,ratio,
                        cˢᵗ,cᵐ,cⁿˡ,
                        yˢᴳ²_0,yˢᴳ³_0,yˢᴳ⁴_0,yˢᴳ⁵_0,yˢᴳ²⁷_0,yˢᴳ³⁰_0)

#------------------------------------------------------------------------------
#--------First, fix the optimal UC schedule------------------------
#------------------------------------------------------------------------------



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

yˢᴳ²_opt=JuMP.value.(yˢᴳ²)       # store optimal values of decision variables
yˢᴳ³_opt=JuMP.value.(yˢᴳ³)
yˢᴳ⁴_opt=JuMP.value.(yˢᴳ⁴)
yˢᴳ⁵_opt=JuMP.value.(yˢᴳ⁵)
yˢᴳ²⁷_opt=JuMP.value.(yˢᴳ²⁷)
yˢᴳ³⁰_opt=JuMP.value.(yˢᴳ³⁰)

#=
ηₘ_1 = JuMP.value.(ηₘ_1)
ηₘ_2 = JuMP.value.(ηₘ_2)

# Contribution to bus 23

contri_SG2_bus23 = ηₘ_2[1,:] * K_c_m[1,16] + K_c_gc[1, 1] *yˢᴳ²_opt  +  ηₘ_1[1:5,:]' * K_c_m[1,1:5] 
contri_SG3_bus23 = ηₘ_2[2,:] * K_c_m[1,17] +K_c_gc[1, 2] *yˢᴳ³_opt  +  ηₘ_1[6:9,:]' * K_c_m[1,6:9]   +  ηₘ_1[1,:] * K_c_m[1,1]
contri_SG4_bus23 = ηₘ_2[3,:] * K_c_m[1,18] +K_c_gc[1, 3] *yˢᴳ⁴_opt  +  ηₘ_1[10:12,:]' * K_c_m[1,10:12]   +  ηₘ_1[2,:] * K_c_m[1,2] +  ηₘ_1[6,:] * K_c_m[1,6]
contri_SG5_bus23 = ηₘ_2[4,:] * K_c_m[1,19] +K_c_gc[1, 4] *yˢᴳ⁵_opt  +  ηₘ_1[13:14,:]' * K_c_m[1,13:14]   +  ηₘ_1[3,:] * K_c_m[1,3] +  ηₘ_1[7,:] * K_c_m[1,7] +  ηₘ_1[10,:] * K_c_m[1,10]
contri_SG27_bus23 = ηₘ_2[5,:] * K_c_m[1,20] +K_c_gc[1, 5] *yˢᴳ²⁷_opt  +  ηₘ_1[4,:] * K_c_m[1,4] +  ηₘ_1[8,:] * K_c_m[1,8] +  ηₘ_1[11,:] * K_c_m[1,11] +  ηₘ_1[13,:] * K_c_m[1,13] +  ηₘ_1[15,:] * K_c_m[1,15]
contri_SG30_bus23 = ηₘ_2[6,:] * K_c_m[1,21] +K_c_gc[1, 6] *yˢᴳ³⁰_opt  +  ηₘ_1[5,:] * K_c_m[1,5] +  ηₘ_1[9,:] * K_c_m[1,9] +  ηₘ_1[12,:] * K_c_m[1,12] +  ηₘ_1[14,:] * K_c_m[1,14] +  ηₘ_1[15,:] * K_c_m[1,15]

plot(contri_SG2_bus23)
plot!(contri_SG3_bus23)
plot!(contri_SG4_bus23)
plot!(contri_SG5_bus23)
plot!(contri_SG27_bus23)
plot!(contri_SG30_bus23)

println(contri_VSG_bus23)




contri_SG2_bus24 = ηₘ_2[1,:] * K_c_m[2,16] + K_c_gc[2, 1] *yˢᴳ²_opt  +  ηₘ_1[1:5,:]' * K_c_m[2,1:5] 
contri_SG3_bus24 = ηₘ_2[2,:] * K_c_m[2,17] +K_c_gc[2, 2] *yˢᴳ³_opt  +  ηₘ_1[6:9,:]' * K_c_m[2,6:9]   +  ηₘ_1[1,:] * K_c_m[2,1]
contri_SG4_bus24 = ηₘ_2[3,:] * K_c_m[2,18] +K_c_gc[2, 3] *yˢᴳ⁴_opt  +  ηₘ_1[10:12,:]' * K_c_m[2,10:12]   +  ηₘ_1[2,:] * K_c_m[2,2] +  ηₘ_1[6,:] * K_c_m[2,6]
contri_SG5_bus24 = ηₘ_2[4,:] * K_c_m[2,19] +K_c_gc[2, 4] *yˢᴳ⁵_opt  +  ηₘ_1[13:14,:]' * K_c_m[2,13:14]   +  ηₘ_1[3,:] * K_c_m[2,3] +  ηₘ_1[7,:] * K_c_m[2,7] +  ηₘ_1[10,:] * K_c_m[2,10]
contri_SG27_bus24 = ηₘ_2[5,:] * K_c_m[2,20] +K_c_gc[2, 5] *yˢᴳ²⁷_opt  +  ηₘ_1[4,:] * K_c_m[2,4] +  ηₘ_1[8,:] * K_c_m[2,8] +  ηₘ_1[11,:] * K_c_m[2,11] +  ηₘ_1[13,:] * K_c_m[2,13] +  ηₘ_1[15,:] * K_c_m[2,15]
contri_SG30_bus24 = ηₘ_2[6,:] * K_c_m[2,21] +K_c_gc[2, 6] *yˢᴳ³⁰_opt  +  ηₘ_1[5,:] * K_c_m[2,5] +  ηₘ_1[9,:] * K_c_m[2,9] +  ηₘ_1[12,:] * K_c_m[2,12] +  ηₘ_1[14,:] * K_c_m[2,14] +  ηₘ_1[15,:] * K_c_m[2,15]



plot(contri_SG2_bus23)
plot!(contri_SG3_bus23)
plot!(contri_SG4_bus23)
plot!(contri_SG5_bus23)
plot!(contri_SG27_bus23)
plot!(contri_SG30_bus23)
plot!(contri_SG2_bus24)
plot!(contri_SG3_bus24)
plot!(contri_SG4_bus24)
plot!(contri_SG5_bus24)
plot!(contri_SG27_bus24)
plot!(contri_SG30_bus24)

contri_VSG_bus23 = α_VSG*K_c_gv[1, 1] + ηₘ_2' *  K_c_m[1, 16:end]
contri_VSG_bus24 = α_VSG*K_c_gv[2, 1] + ηₘ_2' *  K_c_m[2, 16:end]

plot(contri_VSG_bus23)
plot!(contri_VSG_bus24)

reactive_support_b23 = JuMP.value.(Qᴳ[8,:]) 
reactive_support_b24 = JuMP.value.(Qᴳ[9,:]) 
active_b23 = JuMP.value.(Pᴳ[8,:]) 
active_b24 = JuMP.value.(Pᴳ[9,:]) 
active_b1 = JuMP.value.(Pᴳ[7,:]) 

plot(active_b1)

plot(reactive_support_b23)
plot!(reactive_support_b24)
bar!(active_b23 )
bar!(active_b24)

sum(reactive_support_b23)/T
sum(reactive_support_b24)/T

z_23 = JuMP.value.(z_23)
z_24 = JuMP.value.(z_24)

println(Pᴰ)

=#

#------------------------------------------------------------------------------
#--------Second, calculate prices unde the attained UC schedule------------------------
#------------------------------------------------------------------------------

model= Model()                     

    @variable(model, yˢᴳ²[1:T])    
    @variable(model, yˢᴳ³[1:T])            
    @variable(model, yˢᴳ⁴[1:T])           
    @variable(model, yˢᴳ⁵[1:T])           
    @variable(model, yˢᴳ²⁷[1:T])            
    @variable(model, yˢᴳ³⁰[1:T])  

    price_yˢᴳ²=Dict()
    price_yˢᴳ³=Dict()
    price_yˢᴳ⁴=Dict()
    price_yˢᴳ⁵=Dict()
    price_yˢᴳ²⁷=Dict()
    price_yˢᴳ³⁰=Dict()

for t in 1:T
    price_yˢᴳ²[t] = @constraint(model, yˢᴳ²[t] ==yˢᴳ²_opt[t])      # assign the optimal UC schedule to the decision variables
    price_yˢᴳ³[t] = @constraint(model, yˢᴳ³[t] ==yˢᴳ³_opt[t])
    price_yˢᴳ⁴[t] = @constraint(model, yˢᴳ⁴[t] ==yˢᴳ⁴_opt[t])
    price_yˢᴳ⁵[t] = @constraint(model, yˢᴳ⁵[t] ==yˢᴳ⁵_opt[t])
    price_yˢᴳ²⁷[t] = @constraint(model, yˢᴳ²⁷[t] ==yˢᴳ²⁷_opt[t])
    price_yˢᴳ³⁰[t] = @constraint(model, yˢᴳ³⁰[t] ==yˢᴳ³⁰_opt[t])
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

Power_balance_P = Dict()

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

    Power_balance_P[t] = @constraint(model, sum(Pᴳ[:,t]) == Pᴰ[t]  )    
    @constraint(model, sum(Qᴳ[:,t]) == Qᴰ[t] )    

end


#--------------------------Voltage stability constraints--------------------------

@variable(model, z_23[1:T])
@variable(model, z_24[1:T])

@variable(model, ηₘ_1[1:15,1:T] )       
@variable(model, ηₘ_2[1:6,1:T])            

Voltage_constraint_23=Dict()
Voltage_constraint_24=Dict()

for t in 1:T
    # 15+6
    @constraint(model, ηₘ_1[1,t] == yˢᴳ²_opt[t] *yˢᴳ³_opt[t]  )
    @constraint(model, ηₘ_1[2,t] == yˢᴳ²_opt[t] *yˢᴳ⁴_opt[t] )
    @constraint(model, ηₘ_1[3,t] == yˢᴳ²_opt[t] *yˢᴳ⁵_opt[t] )
    @constraint(model, ηₘ_1[4,t] == yˢᴳ²_opt[t] *yˢᴳ²⁷_opt[t] )
    @constraint(model, ηₘ_1[5,t] == yˢᴳ²_opt[t] *yˢᴳ³⁰_opt[t] )
    @constraint(model, ηₘ_1[6,t] == yˢᴳ³_opt[t] *yˢᴳ⁴_opt[t] )
    @constraint(model, ηₘ_1[7,t] == yˢᴳ³_opt[t] *yˢᴳ⁵_opt[t] )
    @constraint(model, ηₘ_1[8,t] == yˢᴳ³_opt[t] *yˢᴳ²⁷_opt[t] )
    @constraint(model, ηₘ_1[9,t] == yˢᴳ³_opt[t] *yˢᴳ³⁰_opt[t] )
    @constraint(model, ηₘ_1[10,t] == yˢᴳ⁴_opt[t] *yˢᴳ⁵_opt[t] )
    @constraint(model, ηₘ_1[11,t] == yˢᴳ⁴_opt[t] *yˢᴳ²⁷_opt[t] )
    @constraint(model, ηₘ_1[12,t] == yˢᴳ⁴_opt[t] *yˢᴳ³⁰_opt[t] )
    @constraint(model, ηₘ_1[13,t] == yˢᴳ⁵_opt[t] *yˢᴳ²⁷_opt[t] )
    @constraint(model, ηₘ_1[14,t] == yˢᴳ⁵_opt[t] *yˢᴳ³⁰_opt[t] )
    @constraint(model, ηₘ_1[15,t] == yˢᴳ²⁷_opt[t] *yˢᴳ³⁰_opt[t] )

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

    Voltage_constraint_23[t] = @constraint(model, [Qᴳ[8,t] + 1/2*z_23[t]*S_Base,  Pᴳ[8,t],  Qᴳ[8,t]] in SecondOrderCone()  )     
    Voltage_constraint_24[t] = @constraint(model, [Qᴳ[9,t] + 1/2*z_24[t]*S_Base,  Pᴳ[9,t],  Qᴳ[9,t]] in SecondOrderCone()  )     

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
set_optimizer_attribute(model, "QCPDual", 1)
optimize!(model)




#------------------------------------------------------------------------------
#--------Analyze results------------------------
#------------------------------------------------------------------------------

UC_price_yˢᴳ²= [max(dual(price_yˢᴳ²[t]), 0.0) for t in 1:T]
UC_price_yˢᴳ³= [max(dual(price_yˢᴳ³[t]), 0.0) for t in 1:T]
UC_price_yˢᴳ⁴= [max(dual(price_yˢᴳ⁴[t]), 0.0) for t in 1:T]
UC_price_yˢᴳ⁵= [max(dual(price_yˢᴳ⁵[t]), 0.0) for t in 1:T]
UC_price_yˢᴳ²⁷= [max(dual(price_yˢᴳ²⁷[t]), 0.0) for t in 1:T]
UC_price_yˢᴳ³⁰= [max(dual(price_yˢᴳ³⁰[t]), 0.0) for t in 1:T]

Cᵁ²=value.(Cᵁ²)
Cᵁ³=value.(Cᵁ³)
Cᵁ⁴=value.(Cᵁ⁴)
Cᵁ⁵=value.(Cᵁ⁵)
Cᵁ²⁷=value.(Cᵁ²⁷)
Cᵁ³⁰=value.(Cᵁ³⁰)

Pᴳ=JuMP.value.(Pᴳ)
Qᴳ=JuMP.value.(Qᴳ)

energy_price = [dual(Power_balance_P[t]) for t in 1:T]

vs_price_23 = [dual(Voltage_constraint_23[t]) for t in 1:T]
vs_price_24 = [dual(Voltage_constraint_24[t]) for t in 1:T]

a=hcat(vs_price_23)
b=hcat(vs_price_24)

mu=zeros(T,2)   # IBGs at bus 23 and 24
lambda_2=zeros(T,2)

for t in 1:T
    aa=a[t,:]
    bb=b[t,:]
    aaa = aa[1]
    bbb = bb[1]
    mu[t,1] = aaa[1]
    lambda_2[t,1] = aaa[3]
    mu[t,2] = bbb[1]
    lambda_2[t,2] = bbb[3]
end

ηₘ_1 = JuMP.value.(ηₘ_1)
ηₘ_2 = JuMP.value.(ηₘ_2)


energy_profit_VSG_1, energy_profitˢᴳ², energy_profitˢᴳ³, energy_profitˢᴳ⁴, energy_profitˢᴳ⁵, energy_profitˢᴳ²⁷, energy_profitˢᴳ³⁰, energy_profit_IBG_23, energy_profit_IBG_24,
       pro_SG_no_vs,
       revenue_IBG_23_give_reactive, revenue_IBG_24_give_reactive,
       VS_revenue_VSG_1, VS_revenue_SGs = settlement_approx( Pᴳ,Qᴳ,mu,lambda_2,energy_price, 
                    UC_price_yˢᴳ²,UC_price_yˢᴳ³,UC_price_yˢᴳ⁴,UC_price_yˢᴳ⁵,UC_price_yˢᴳ²⁷,UC_price_yˢᴳ³⁰,
                    yˢᴳ²_opt,yˢᴳ³_opt,yˢᴳ⁴_opt,yˢᴳ⁵_opt,yˢᴳ²⁷_opt,yˢᴳ³⁰_opt,
                    Cᵁ²,Cᵁ³,Cᵁ⁴,Cᵁ⁵,Cᵁ²⁷,Cᵁ³⁰,
                    cᵐ,cⁿˡ,α_VSG,S_Base,
                    K_c_gc, K_c_m, ηₘ_1, ηₘ_2)

uplifts = zeros(6,T)
for t in 1:T
    uplifts[1,t] = UC_price_yˢᴳ²[t]*yˢᴳ²_opt[t] 
    uplifts[2,t] = UC_price_yˢᴳ³[t]*yˢᴳ³_opt[t]
    uplifts[3,t] = UC_price_yˢᴳ⁴[t]*yˢᴳ⁴_opt[t]
    uplifts[4,t] = UC_price_yˢᴳ⁵[t]*yˢᴳ⁵_opt[t]
    uplifts[5,t] = UC_price_yˢᴳ²⁷[t]*yˢᴳ²⁷_opt[t]
    uplifts[6,t] = UC_price_yˢᴳ³⁰[t]*yˢᴳ³⁰_opt[t]
end

energy_profit_SGs = [sum(energy_profitˢᴳ²), sum(energy_profitˢᴳ³), sum(energy_profitˢᴳ⁴), sum(energy_profitˢᴳ⁵), sum(energy_profitˢᴳ²⁷), sum(energy_profitˢᴳ³⁰)]

SCR_23 = JuMP.value.(z_23)
SCR_24 = JuMP.value.(z_24)
SCR = zeros(2)
SCR[1] = sum(SCR_23)/T
SCR[2] = sum(SCR_24)/T

wind_curtail_23 = zeros(T)
wind_curtail_24 = zeros(T)
for t in 1:T
    wind_curtail_23[t] = α_IBG[1,t]*Pᴵᴮᴳₘₐₓ[1] - Pᴳ[8,t]
    wind_curtail_24[t] = α_IBG[2,t]*Pᴵᴮᴳₘₐₓ[2] - Pᴳ[9,t]
end
wind_curtailment = zeros(2)
wind_curtailment[1] = sum(wind_curtail_23)/T
wind_curtailment[2] = sum(wind_curtail_24)/T

return energy_price, mu[:,1], mu[:,2], mu[:,1] .-lambda_2[:,1], mu[:,2] .-lambda_2[:,2], pro_SG_no_vs, VS_revenue_SGs, uplifts, energy_profit_VSG_1, energy_profit_IBG_23, energy_profit_IBG_24, VS_revenue_VSG_1, revenue_IBG_23_give_reactive, revenue_IBG_24_give_reactive, energy_profit_SGs,
        wind_curtailment, SCR

end
