# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Pricing voltage stability via Primal-Dual formulation
# 12.Apr.2026


#-------------------------------------------------------
#-------------------Define primal model-------------------
#-------------------------------------------------------

function pricing_VS_PD(K_c_gc, K_c_gv, K_c_m,
                        T,Pᴰ,Qᴰ,Pˢᴳₘₐₓ,Pˢᴳₘᵢₙ,Qˢᴳₘₐₓ,Qˢᴳₘᵢₙ,Pⱽˢᴳₘₐₓ,Pⱽˢᴳₘᵢₙ,Qⱽˢᴳₘₐₓ,Qⱽˢᴳₘᵢₙ,Pᴵᴮᴳₘₐₓ,Pᴵᴮᴳₘᵢₙ,Qᴵᴮᴳₘₐₓ,Qᴵᴮᴳₘᵢₙ,
                        Sᵐᵃˣ_gc,Sᵐᵃˣ_gv,Sᵐᵃˣ_c,S_Base,α_VSG,α_IBG,ratio,
                        cˢᵗ,cᵐ,cⁿˡ,
                        yˢᴳ²_0,yˢᴳ³_0,yˢᴳ⁴_0,yˢᴳ⁵_0,yˢᴳ²⁷_0,yˢᴳ³⁰_0,
                        LCOE_gv, LCOE_gf)

model= Model()                     

@variable(model, yˢᴳ²[1:T] ,Bin)    # ,Bin  >= 0
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

@constraint(model, Cᵁ²[1]>=(yˢᴳ²[1]-yˢᴳ²_0)*cˢᵗ[1])     #  σˢᴳ²_ˢᵗ
@constraint(model, Cᵁ³[1]>=(yˢᴳ³[1]-yˢᴳ³_0)*cˢᵗ[2])     #  σˢᴳ³_ˢᵗ
@constraint(model, Cᵁ⁴[1]>=(yˢᴳ⁴[1]-yˢᴳ⁴_0)*cˢᵗ[3])     #  σˢᴳ⁴_ˢᵗ
@constraint(model, Cᵁ⁵[1]>=(yˢᴳ⁵[1]-yˢᴳ⁵_0)*cˢᵗ[4])     #  σˢᴳ⁵_ˢᵗ
@constraint(model, Cᵁ²⁷[1]>=(yˢᴳ²⁷[1]-yˢᴳ²⁷_0)*cˢᵗ[5])  #  σˢᴳ²⁷_ˢᵗ
@constraint(model, Cᵁ³⁰[1]>=(yˢᴳ³⁰[1]-yˢᴳ³⁰_0)*cˢᵗ[6])  #  σˢᴳ³⁰_ˢᵗ

for t in 2:T

    @constraint(model, Cᵁ²[t]>=(yˢᴳ²[t]-yˢᴳ²[t-1])*cˢᵗ[1])     #  σˢᴳ²_ˢᵗ
    @constraint(model, Cᵁ³[t]>=(yˢᴳ³[t]-yˢᴳ³[t-1])*cˢᵗ[2])     #  σˢᴳ³_ˢᵗ
    @constraint(model, Cᵁ⁴[t]>=(yˢᴳ⁴[t]-yˢᴳ⁴[t-1])*cˢᵗ[3])     #  σˢᴳ⁴_ˢᵗ
    @constraint(model, Cᵁ⁵[t]>=(yˢᴳ⁵[t]-yˢᴳ⁵[t-1])*cˢᵗ[4])     #  σˢᴳ⁵_ˢᵗ
    @constraint(model, Cᵁ²⁷[t]>=(yˢᴳ²⁷[t]-yˢᴳ²⁷[t-1])*cˢᵗ[5])  #  σˢᴳ²⁷_ˢᵗ
    @constraint(model, Cᵁ³⁰[t]>=(yˢᴳ³⁰[t]-yˢᴳ³⁰[t-1])*cˢᵗ[6])  #  σˢᴳ³⁰_ˢᵗ

end

for t in 1:T

    @constraint(model, Pᴳ[1,t] <= yˢᴳ²[t]*Pˢᴳₘₐₓ[1])    #  τᵖ_ᵘᵇ_ˢᴳ²     
    @constraint(model, yˢᴳ²[t]*Pˢᴳₘᵢₙ[1] <= Pᴳ[1,t])     #  τᵖ_ˡᵇ_ˢᴳ² 
    @constraint(model, Pᴳ[2,t] <= yˢᴳ³[t]*Pˢᴳₘₐₓ[2])    #  τᵖ_ᵘᵇ_ˢᴳ³
    @constraint(model, yˢᴳ³[t]*Pˢᴳₘᵢₙ[2] <= Pᴳ[2,t])     #  τᵖ_ˡᵇ_ˢᴳ³ 
    @constraint(model, Pᴳ[3,t] <= yˢᴳ⁴[t]*Pˢᴳₘₐₓ[3])     #  τᵖ_ᵘᵇ_ˢᴳ⁴
    @constraint(model, yˢᴳ⁴[t]*Pˢᴳₘᵢₙ[3] <= Pᴳ[3,t])     #  τᵖ_ˡᵇ_ˢᴳ⁴
    @constraint(model, Pᴳ[4,t] <= yˢᴳ⁵[t]*Pˢᴳₘₐₓ[4])     #  τᵖ_ᵘᵇ_ˢᴳ⁵
    @constraint(model, yˢᴳ⁵[t]*Pˢᴳₘᵢₙ[4] <= Pᴳ[4,t])     #  τᵖ_ˡᵇ_ˢᴳ⁵
    @constraint(model, Pᴳ[5,t] <= yˢᴳ²⁷[t]*Pˢᴳₘₐₓ[5])   #  τᵖ_ᵘᵇ_ˢᴳ²⁷ 
    @constraint(model, yˢᴳ²⁷[t]*Pˢᴳₘᵢₙ[5] <=Pᴳ[5,t])    #  τᵖ_ˡᵇ_ˢᴳ²⁷
    @constraint(model, Pᴳ[6,t] <= yˢᴳ³⁰[t]*Pˢᴳₘₐₓ[6])   #  τᵖ_ᵘᵇ_ˢᴳ³⁰     
    @constraint(model, yˢᴳ³⁰[t]*Pˢᴳₘᵢₙ[6] <= Pᴳ[6,t])    #  τᵖ_ˡᵇ_ˢᴳ³⁰

    @constraint(model, Qᴳ[1,t] <= yˢᴳ²[t]*Qˢᴳₘₐₓ[1])   # τᵠ_ᵘᵇ_ˢᴳ²
    @constraint(model, yˢᴳ²[t]*Qˢᴳₘᵢₙ[1] <= Qᴳ[1,t])   # τᵠ_ˡᵇ_ˢᴳ²    
    @constraint(model, Qᴳ[2,t] <= yˢᴳ³[t]*Qˢᴳₘₐₓ[2])   # τᵠ_ᵘᵇ_ˢᴳ³      
    @constraint(model, yˢᴳ³[t]*Qˢᴳₘᵢₙ[2] <= Qᴳ[2,t])   # τᵠ_ˡᵇ_ˢᴳ³        
    @constraint(model, Qᴳ[3,t] <= yˢᴳ⁴[t]*Qˢᴳₘₐₓ[3])   # τᵠ_ᵘᵇ_ˢᴳ⁴       
    @constraint(model, yˢᴳ⁴[t]*Qˢᴳₘᵢₙ[3] <= Qᴳ[3,t])   # τᵠ_ˡᵇ_ˢᴳ⁴   
    @constraint(model, Qᴳ[4,t] <= yˢᴳ⁵[t]*Qˢᴳₘₐₓ[4])   # τᵠ_ᵘᵇ_ˢᴳ⁵       
    @constraint(model, yˢᴳ⁵[t]*Qˢᴳₘᵢₙ[4] <= Qᴳ[4,t])    # τᵠ_ˡᵇ_ˢᴳ⁵
    @constraint(model, Qᴳ[5,t] <= yˢᴳ²⁷[t]*Qˢᴳₘₐₓ[5])  # τᵠ_ᵘᵇ_ˢᴳ²⁷       
    @constraint(model, yˢᴳ²⁷[t]*Qˢᴳₘᵢₙ[5] <=Qᴳ[5,t])   # τᵠ_ˡᵇ_ˢᴳ²⁷
    @constraint(model, Qᴳ[6,t] <= yˢᴳ³⁰[t]*Qˢᴳₘₐₓ[6])  # τᵠ_ᵘᵇ_ˢᴳ³⁰      
    @constraint(model, yˢᴳ³⁰[t]*Qˢᴳₘᵢₙ[6] <= Qᴳ[6,t])  # τᵠ_ˡᵇ_ˢᴳ³⁰      

    @constraint(model, Pᴳ[7,t] <= α_VSG[t]*Pⱽˢᴳₘₐₓ[1])      # τᵖ_ᵘᵇ_VSG_1                
    @constraint(model, Pⱽˢᴳₘᵢₙ[1] <= Pᴳ[7,t])       
    @constraint(model, Pᴳ[8,t] <= α_IBG[1,t]*Pᴵᴮᴳₘₐₓ[1])     # τᵖ_ᵘᵇ_IBG_23   
    @constraint(model, Pᴵᴮᴳₘᵢₙ[1] <= Pᴳ[8,t])         
    @constraint(model, Pᴳ[9,t] <= α_IBG[2,t]*Pᴵᴮᴳₘₐₓ[2])     # τᵖ_ᵘᵇ_IBG_24  
    @constraint(model, Pᴵᴮᴳₘᵢₙ[2] <= Pᴳ[9,t])        

    @constraint(model, Qᴳ[7,t] <= Qⱽˢᴳₘₐₓ[1])    # τᵠ_ᵘᵇ_VSG_1   
    @constraint(model, Qⱽˢᴳₘᵢₙ[1] <= Qᴳ[7,t])    # τᵠ_ˡᵇ_VSG_1
    @constraint(model, Qᴳ[8,t] <= ratio*Qᴵᴮᴳₘₐₓ[1])    # τᵠ_ᵘᵇ_IBG_23   
    @constraint(model, Qᴵᴮᴳₘᵢₙ[1] <= Qᴳ[8,t])    # τᵠ_ˡᵇ_IBG_23
    @constraint(model, Qᴳ[9,t] <= ratio*Qᴵᴮᴳₘₐₓ[2])    # τᵠ_ᵘᵇ_IBG_24   
    @constraint(model, Qᴵᴮᴳₘᵢₙ[2] <= Qᴳ[9,t])    # τᵠ_ˡᵇ_IBG_24

    @constraint(model, [ Sᵐᵃˣ_gc[1], Pᴳ[1,t], Qᴳ[1,t] ] in SecondOrderCone()  )   #  υ_1_ˢᴳ²，  υ_2_ˢᴳ²，  υ_3_ˢᴳ²
    @constraint(model, [ Sᵐᵃˣ_gc[2], Pᴳ[2,t], Qᴳ[2,t] ] in SecondOrderCone()  )   #  υ_1_ˢᴳ³，  υ_2_ˢᴳ³，  υ_3_ˢᴳ³
    @constraint(model, [ Sᵐᵃˣ_gc[3], Pᴳ[3,t], Qᴳ[3,t] ] in SecondOrderCone()  )   #  υ_1_ˢᴳ⁴，  υ_2_ˢᴳ⁴，  υ_3_ˢᴳ⁴
    @constraint(model, [ Sᵐᵃˣ_gc[4], Pᴳ[4,t], Qᴳ[4,t] ] in SecondOrderCone()  )   #  υ_1_ˢᴳ⁵，  υ_2_ˢᴳ⁵，  υ_3_ˢᴳ⁵
    @constraint(model, [ Sᵐᵃˣ_gc[5], Pᴳ[5,t], Qᴳ[5,t] ] in SecondOrderCone()  )  #  υ_1_ˢᴳ²⁷，  υ_2_ˢᴳ²⁷，  υ_3_ˢᴳ²⁷
    @constraint(model, [ Sᵐᵃˣ_gc[6], Pᴳ[6,t], Qᴳ[6,t] ] in SecondOrderCone()  )  #  υ_1_ˢᴳ³⁰，  υ_2_ˢᴳ³⁰，  υ_3_ˢᴳ³⁰
    @constraint(model, [ Sᵐᵃˣ_gv[1], Pᴳ[7,t], Qᴳ[7,t] ] in SecondOrderCone()  )    #  υ_1_VSG_1，  υ_2_VSG_1，  υ_3_VSG_1
    @constraint(model, [ Sᵐᵃˣ_c[1], Pᴳ[8,t], Qᴳ[8,t] ] in SecondOrderCone()  )    #  υ_1_IBG_23，  υ_2_IBG_23，  υ_3_IBG_23
    @constraint(model, [ Sᵐᵃˣ_c[2], Pᴳ[9,t], Qᴳ[9,t] ] in SecondOrderCone()  )    #  υ_1_IBG_24，  υ_2_IBG_24，  υ_3_IBG_24

    @constraint(model, sum(Pᴳ[:,t]) == Pᴰ[t]  )    #   λᴱ
    @constraint(model, sum(Qᴳ[:,t]) == Qᴰ[t] )    #    ϕ

end


#--------------------------Voltage stability constraints--------------------------

@variable(model, z_23[1:T])
@variable(model, z_24[1:T])

@variable(model, ηₘ_1[1:15,1:T] >= 0)       
@variable(model, ηₘ_2[1:6,1:T])            

for t in 1:T
    # 15+6
    @constraint(model, ηₘ_1[1,t]>=yˢᴳ²[t]+yˢᴳ³[t]-1)    # γ_ᵐⁱⁿ[1]
    @constraint(model, ηₘ_1[1,t]<=yˢᴳ²[t])              # γ_ᵐᵃˣ[1,1]
    @constraint(model, ηₘ_1[1,t]<=yˢᴳ³[t])              # γ_ᵐᵃˣ[2,1]
    @constraint(model, ηₘ_1[2,t]>=yˢᴳ²[t]+yˢᴳ⁴[t]-1)    # γ_ᵐⁱⁿ[2]
    @constraint(model, ηₘ_1[2,t]<=yˢᴳ²[t])              # γ_ᵐᵃˣ[1,2] 
    @constraint(model, ηₘ_1[2,t]<=yˢᴳ⁴[t])              # γ_ᵐᵃˣ[2,2]
    @constraint(model, ηₘ_1[3,t]>=yˢᴳ²[t]+yˢᴳ⁵[t]-1)    # γ_ᵐⁱⁿ[3]
    @constraint(model, ηₘ_1[3,t]<=yˢᴳ²[t])              # γ_ᵐᵃˣ[1,3]
    @constraint(model, ηₘ_1[3,t]<=yˢᴳ⁵[t])              # γ_ᵐᵃˣ[2,3]
    @constraint(model, ηₘ_1[4,t]>=yˢᴳ²[t]+yˢᴳ²⁷[t]-1)   # γ_ᵐⁱⁿ[4]
    @constraint(model, ηₘ_1[4,t]<=yˢᴳ²[t])              # γ_ᵐᵃˣ[1,4]
    @constraint(model, ηₘ_1[4,t]<=yˢᴳ²⁷[t])             # γ_ᵐᵃˣ[2,4]
    @constraint(model, ηₘ_1[5,t]>=yˢᴳ²[t]+yˢᴳ³⁰[t]-1)   # γ_ᵐⁱⁿ[5]
    @constraint(model, ηₘ_1[5,t]<=yˢᴳ²[t])              # γ_ᵐᵃˣ[1,5]
    @constraint(model, ηₘ_1[5,t]<=yˢᴳ³⁰[t])             # γ_ᵐᵃˣ[2,5]

    @constraint(model, ηₘ_1[6,t]>=yˢᴳ³[t]+yˢᴳ⁴[t]-1)   # γ_ᵐⁱⁿ[6]
    @constraint(model, ηₘ_1[6,t]<=yˢᴳ³[t])             # γ_ᵐᵃˣ[1,6]
    @constraint(model, ηₘ_1[6,t]<=yˢᴳ⁴[t])             # γ_ᵐᵃˣ[2,6]
    @constraint(model, ηₘ_1[7,t]>=yˢᴳ³[t]+yˢᴳ⁵[t]-1)   # γ_ᵐⁱⁿ[7]
    @constraint(model, ηₘ_1[7,t]<=yˢᴳ³[t])             # γ_ᵐᵃˣ[1,7]
    @constraint(model, ηₘ_1[7,t]<=yˢᴳ⁵[t])             # γ_ᵐᵃˣ[2,7]
    @constraint(model, ηₘ_1[8,t]>=yˢᴳ³[t]+yˢᴳ²⁷[t]-1)  # γ_ᵐⁱⁿ[8] 
    @constraint(model, ηₘ_1[8,t]<=yˢᴳ³[t])             # γ_ᵐᵃˣ[1,8]
    @constraint(model, ηₘ_1[8,t]<=yˢᴳ²⁷[t])            # γ_ᵐᵃˣ[2,8]
    @constraint(model, ηₘ_1[9,t]>=yˢᴳ³[t]+yˢᴳ³⁰[t]-1)  # γ_ᵐⁱⁿ[9] 
    @constraint(model, ηₘ_1[9,t]<=yˢᴳ³[t])             # γ_ᵐᵃˣ[1,9]
    @constraint(model, ηₘ_1[9,t]<=yˢᴳ³⁰[t])            # γ_ᵐᵃˣ[2,9]

    @constraint(model, ηₘ_1[10,t]>=yˢᴳ⁴[t]+yˢᴳ⁵[t]-1)   # γ_ᵐⁱⁿ[10]
    @constraint(model, ηₘ_1[10,t]<=yˢᴳ⁴[t])             # γ_ᵐᵃˣ[1,10]
    @constraint(model, ηₘ_1[10,t]<=yˢᴳ⁵[t])             # γ_ᵐᵃˣ[2,10]
    @constraint(model, ηₘ_1[11,t]>=yˢᴳ⁴[t]+yˢᴳ²⁷[t]-1)  # γ_ᵐⁱⁿ[11] 
    @constraint(model, ηₘ_1[11,t]<=yˢᴳ⁴[t])             # γ_ᵐᵃˣ[1,11]
    @constraint(model, ηₘ_1[11,t]<=yˢᴳ²⁷[t])            # γ_ᵐᵃˣ[2,11]
    @constraint(model, ηₘ_1[12,t]>=yˢᴳ⁴[t]+yˢᴳ³⁰[t]-1)  # γ_ᵐⁱⁿ[12] 
    @constraint(model, ηₘ_1[12,t]<=yˢᴳ⁴[t])             # γ_ᵐᵃˣ[1,12]
    @constraint(model, ηₘ_1[12,t]<=yˢᴳ³⁰[t])            # γ_ᵐᵃˣ[2,12]

    @constraint(model, ηₘ_1[13,t]>=yˢᴳ⁵[t]+yˢᴳ²⁷[t]-1)  # γ_ᵐⁱⁿ[13] 
    @constraint(model, ηₘ_1[13,t]<=yˢᴳ⁵[t])             # γ_ᵐᵃˣ[1,13]
    @constraint(model, ηₘ_1[13,t]<=yˢᴳ²⁷[t])            # γ_ᵐᵃˣ[2,13]       
    @constraint(model, ηₘ_1[14,t]>=yˢᴳ⁵[t]+yˢᴳ³⁰[t]-1)  # γ_ᵐⁱⁿ[14]  
    @constraint(model, ηₘ_1[14,t]<=yˢᴳ⁵[t])             # γ_ᵐᵃˣ[1,14]
    @constraint(model, ηₘ_1[14,t]<=yˢᴳ³⁰[t])            # γ_ᵐᵃˣ[2,14]

    @constraint(model, ηₘ_1[15,t]>=yˢᴳ²⁷[t]+yˢᴳ³⁰[t]-1)  # γ_ᵐⁱⁿ[15] 
    @constraint(model, ηₘ_1[15,t]<=yˢᴳ²⁷[t])             # γ_ᵐᵃˣ[1,15]
    @constraint(model, ηₘ_1[15,t]<=yˢᴳ³⁰[t])             # γ_ᵐᵃˣ[2,15]

    @constraint(model, ηₘ_2[1,t]==yˢᴳ²[t]*α_VSG[t])        # γ[1]
    @constraint(model, ηₘ_2[2,t]==yˢᴳ³[t]*α_VSG[t])        # γ[2]
    @constraint(model, ηₘ_2[3,t]==yˢᴳ⁴[t]*α_VSG[t])        # γ[3]
    @constraint(model, ηₘ_2[4,t]==yˢᴳ⁵[t]*α_VSG[t])        # γ[4]
    @constraint(model, ηₘ_2[5,t]==yˢᴳ²⁷[t]*α_VSG[t])       # γ[5]
    @constraint(model, ηₘ_2[6,t]==yˢᴳ³⁰[t]*α_VSG[t])       # γ[6]

    @constraint(model, z_23[t] ==  ( ( K_c_gc[1, 1] *yˢᴳ²[t] +K_c_gc[1, 2] *yˢᴳ³[t] +K_c_gc[1, 3] *yˢᴳ⁴[t] 
    +K_c_gc[1, 4] *yˢᴳ⁵[t] +K_c_gc[1, 5] *yˢᴳ²⁷[t] +K_c_gc[1, 6] *yˢᴳ³⁰[t] )   +α_VSG[t]*K_c_gv[1, 1]
    + sum(K_c_m[1, 1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[1, 16:end] .*ηₘ_2[:,t]   ) ) )     #   ξ_gf_23

    @constraint(model, z_24[t] == ( ( K_c_gc[2, 1] *yˢᴳ²[t] +K_c_gc[2, 2] *yˢᴳ³[t] +K_c_gc[2, 3] *yˢᴳ⁴[t] 
    +K_c_gc[2, 4] *yˢᴳ⁵[t] +K_c_gc[2, 5] *yˢᴳ²⁷[t] +K_c_gc[2, 6] *yˢᴳ³⁰[t] )   +α_VSG[t]*K_c_gv[2, 1]
    + sum(K_c_m[2, 1:15] .*ηₘ_1[:,t] ) +  sum(K_c_m[2, 16:end] .*ηₘ_2[:,t]   ) ) )     #   ξ_gf_24

    @constraint(model, [Qᴳ[8,t] + 1/2*z_23[t]*S_Base,  Pᴳ[8,t],  Qᴳ[8,t]] in SecondOrderCone()  )     # λ_1_23,   λ_2_23,   μ_23
    @constraint(model, [Qᴳ[9,t] + 1/2*z_24[t]*S_Base,  Pᴳ[9,t],  Qᴳ[9,t]] in SecondOrderCone()  )     # λ_1_24,   λ_2_24,   μ_24

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

Pᴳ = value.(Pᴳ)
Qᴳ = value.(Qᴳ)
yˢᴳ² = value.(yˢᴳ²)
yˢᴳ³ = value.(yˢᴳ³)
yˢᴳ⁴ = value.(yˢᴳ⁴)
yˢᴳ⁵ = value.(yˢᴳ⁵)
yˢᴳ²⁷ = value.(yˢᴳ²⁷)
yˢᴳ³⁰ = value.(yˢᴳ³⁰)
Cᵁ² = value.(Cᵁ²)
Cᵁ³ = value.(Cᵁ³)
Cᵁ⁴ = value.(Cᵁ⁴)
Cᵁ⁵ = value.(Cᵁ⁵)
Cᵁ²⁷ = value.(Cᵁ²⁷)
Cᵁ³⁰ = value.(Cᵁ³⁰)
ηₘ_1 = value.(ηₘ_1)
ηₘ_2 = value.(ηₘ_2)


#-------------------------------------------------------
#-------------------Define dual model-------------------
#-------------------------------------------------------

model= Model()  

@variable(model, ψˢᴳ²_ᵐᵃˣ[1:T] >= 0)    
@variable(model, ψˢᴳ³_ᵐᵃˣ[1:T] >= 0)
@variable(model, ψˢᴳ⁴_ᵐᵃˣ[1:T] >= 0)
@variable(model, ψˢᴳ⁵_ᵐᵃˣ[1:T] >= 0)
@variable(model, ψˢᴳ²⁷_ᵐᵃˣ[1:T] >= 0)
@variable(model, ψˢᴳ³⁰_ᵐᵃˣ[1:T] >= 0)   
@variable(model, σˢᴳ²_ˢᵗ[1:T] >= 0)
@variable(model, σˢᴳ³_ˢᵗ[1:T] >= 0)
@variable(model, σˢᴳ⁴_ˢᵗ[1:T] >= 0)
@variable(model, σˢᴳ⁵_ˢᵗ[1:T] >= 0)
@variable(model, σˢᴳ²⁷_ˢᵗ[1:T] >= 0)
@variable(model, σˢᴳ³⁰_ˢᵗ[1:T] >= 0)
@variable(model, τᵖ_ᵘᵇ_ˢᴳ²[1:T] >= 0)
@variable(model, τᵖ_ˡᵇ_ˢᴳ²[1:T] >= 0)
@variable(model, τᵖ_ᵘᵇ_ˢᴳ³[1:T] >= 0)
@variable(model, τᵖ_ˡᵇ_ˢᴳ³[1:T] >= 0)
@variable(model, τᵖ_ᵘᵇ_ˢᴳ⁴[1:T] >= 0)
@variable(model, τᵖ_ˡᵇ_ˢᴳ⁴[1:T] >= 0)
@variable(model, τᵖ_ᵘᵇ_ˢᴳ⁵[1:T] >= 0)
@variable(model, τᵖ_ˡᵇ_ˢᴳ⁵[1:T] >= 0)
@variable(model, τᵖ_ᵘᵇ_ˢᴳ²⁷[1:T] >= 0)
@variable(model, τᵖ_ˡᵇ_ˢᴳ²⁷[1:T] >= 0)
@variable(model, τᵖ_ᵘᵇ_ˢᴳ³⁰[1:T] >= 0)
@variable(model, τᵖ_ˡᵇ_ˢᴳ³⁰[1:T] >= 0)
@variable(model, τᵠ_ᵘᵇ_ˢᴳ²[1:T] >= 0)
@variable(model, τᵠ_ˡᵇ_ˢᴳ²[1:T] >= 0)
@variable(model, τᵠ_ᵘᵇ_ˢᴳ³[1:T] >= 0)
@variable(model, τᵠ_ˡᵇ_ˢᴳ³[1:T] >= 0)
@variable(model, τᵠ_ᵘᵇ_ˢᴳ⁴[1:T] >= 0)
@variable(model, τᵠ_ˡᵇ_ˢᴳ⁴[1:T] >= 0)
@variable(model, τᵠ_ᵘᵇ_ˢᴳ⁵[1:T] >= 0)
@variable(model, τᵠ_ˡᵇ_ˢᴳ⁵[1:T] >= 0)
@variable(model, τᵠ_ᵘᵇ_ˢᴳ²⁷[1:T] >= 0)
@variable(model, τᵠ_ˡᵇ_ˢᴳ²⁷[1:T] >= 0)
@variable(model, τᵠ_ᵘᵇ_ˢᴳ³⁰[1:T] >= 0)
@variable(model, τᵠ_ˡᵇ_ˢᴳ³⁰[1:T] >= 0)
@variable(model, τᵖ_ᵘᵇ_VSG_1[1:T] >= 0)
@variable(model, τᵖ_ᵘᵇ_IBG_23[1:T] >= 0)
@variable(model, τᵖ_ᵘᵇ_IBG_24[1:T] >= 0)
@variable(model, τᵠ_ᵘᵇ_VSG_1[1:T] >= 0)
@variable(model, τᵠ_ˡᵇ_VSG_1[1:T] >= 0)
@variable(model, τᵠ_ᵘᵇ_IBG_23[1:T] >= 0)
@variable(model, τᵠ_ˡᵇ_IBG_23[1:T] >= 0)
@variable(model, τᵠ_ᵘᵇ_IBG_24[1:T] >= 0)
@variable(model, τᵠ_ˡᵇ_IBG_24[1:T] >= 0)

@variable(model, γ_ᵐⁱⁿ[1:15,1:T] >= 0)
@variable(model, γ_ᵐᵃˣ[1:2,1:15,1:T] >= 0)
@variable(model, γ[1:6,1:T] )

@variable(model, υ_1_ˢᴳ²[1:T])    
@variable(model, υ_2_ˢᴳ²[1:T])
@variable(model, υ_3_ˢᴳ²[1:T] >= 0)
@variable(model, υ_1_ˢᴳ³[1:T])
@variable(model, υ_2_ˢᴳ³[1:T])
@variable(model, υ_3_ˢᴳ³[1:T] >= 0)
@variable(model, υ_1_ˢᴳ⁴[1:T])
@variable(model, υ_2_ˢᴳ⁴[1:T])
@variable(model, υ_3_ˢᴳ⁴[1:T] >= 0)
@variable(model, υ_1_ˢᴳ⁵[1:T])
@variable(model, υ_2_ˢᴳ⁵[1:T])
@variable(model, υ_3_ˢᴳ⁵[1:T] >= 0)
@variable(model, υ_1_ˢᴳ²⁷[1:T])
@variable(model, υ_2_ˢᴳ²⁷[1:T])
@variable(model, υ_3_ˢᴳ²⁷[1:T] >= 0)
@variable(model, υ_1_ˢᴳ³⁰[1:T])
@variable(model, υ_2_ˢᴳ³⁰[1:T])
@variable(model, υ_3_ˢᴳ³⁰[1:T] >= 0)
@variable(model, υ_1_VSG_1[1:T])
@variable(model, υ_2_VSG_1[1:T])
@variable(model, υ_3_VSG_1[1:T] >= 0)
@variable(model, υ_1_IBG_23[1:T])
@variable(model, υ_2_IBG_23[1:T])
@variable(model, υ_3_IBG_23[1:T] >= 0)
@variable(model, υ_1_IBG_24[1:T])
@variable(model, υ_2_IBG_24[1:T])
@variable(model, υ_3_IBG_24[1:T] >= 0)

@variable(model, λᴱ[1:T] )
@variable(model, ϕ[1:T] )
    
@variable(model, ξ_gf_23[1:T] )
@variable(model, ξ_gf_24[1:T] )
@variable(model, λ_1_23[1:T] )
@variable(model, λ_2_23[1:T] )
@variable(model, μ_23[1:T] >= 0)
@variable(model, λ_1_24[1:T] )
@variable(model, λ_2_24[1:T] )
@variable(model, μ_24[1:T] >= 0)


#--------------------------------
#  Dual constraints
#--------------------------------

for t in 1:T-1
# ============================================================
# yˢᴳ²...yˢᴳ³⁰

    @constraint(model, cⁿˡ[1] + ψˢᴳ²_ᵐᵃˣ[t] - Pˢᴳₘₐₓ[1]*τᵖ_ᵘᵇ_ˢᴳ²[t] + Pˢᴳₘᵢₙ[1]*τᵖ_ˡᵇ_ˢᴳ²[t] 
                    - Qˢᴳₘₐₓ[1]*τᵠ_ᵘᵇ_ˢᴳ²[t] + Qˢᴳₘᵢₙ[1]*τᵠ_ˡᵇ_ˢᴳ²[t] 
                    + K_c_gc[1, 1]*ξ_gf_23[t] + K_c_gc[2, 1]*ξ_gf_24[t] + α_VSG[t]*γ[1,t]
                    + cˢᵗ[1]*σˢᴳ²_ˢᵗ[t] - cˢᵗ[1]*σˢᴳ²_ˢᵗ[t+1]
                    + sum(γ_ᵐⁱⁿ[1:5,t]) - sum(γ_ᵐᵃˣ[1,1:5,t]) >= 0 )
    @constraint(model, cⁿˡ[2] + ψˢᴳ³_ᵐᵃˣ[t] - Pˢᴳₘₐₓ[2]*τᵖ_ᵘᵇ_ˢᴳ³[t] + Pˢᴳₘᵢₙ[2]*τᵖ_ˡᵇ_ˢᴳ³[t] 
                    - Qˢᴳₘₐₓ[2]*τᵠ_ᵘᵇ_ˢᴳ³[t] + Qˢᴳₘᵢₙ[2]*τᵠ_ˡᵇ_ˢᴳ³[t] 
                    + K_c_gc[1, 2]*ξ_gf_23[t] + K_c_gc[2, 2]*ξ_gf_24[t] + α_VSG[t]*γ[2,t]
                    + cˢᵗ[2]*σˢᴳ³_ˢᵗ[t] - cˢᵗ[2]*σˢᴳ³_ˢᵗ[t+1]
                    + sum(γ_ᵐⁱⁿ[6:9,t]) - sum(γ_ᵐᵃˣ[1,6:9,t]) + γ_ᵐⁱⁿ[1,t] - γ_ᵐᵃˣ[2,1,t] >= 0 )
    @constraint(model, cⁿˡ[3] + ψˢᴳ⁴_ᵐᵃˣ[t] - Pˢᴳₘₐₓ[3]*τᵖ_ᵘᵇ_ˢᴳ⁴[t] + Pˢᴳₘᵢₙ[3]*τᵖ_ˡᵇ_ˢᴳ⁴[t] 
                    - Qˢᴳₘₐₓ[3]*τᵠ_ᵘᵇ_ˢᴳ⁴[t] + Qˢᴳₘᵢₙ[3]*τᵠ_ˡᵇ_ˢᴳ⁴[t] 
                    + K_c_gc[1, 3]*ξ_gf_23[t] + K_c_gc[2, 3]*ξ_gf_24[t] + α_VSG[t]*γ[3,t]
                    + cˢᵗ[3]*σˢᴳ⁴_ˢᵗ[t] - cˢᵗ[3]*σˢᴳ⁴_ˢᵗ[t+1]
                    + sum(γ_ᵐⁱⁿ[10:12,t]) - sum(γ_ᵐᵃˣ[1,10:12,t]) + γ_ᵐⁱⁿ[2,t] + γ_ᵐⁱⁿ[6,t] - γ_ᵐᵃˣ[2,2,t] - γ_ᵐᵃˣ[2,6,t] >= 0 )
    @constraint(model, cⁿˡ[4] + ψˢᴳ⁵_ᵐᵃˣ[t] - Pˢᴳₘₐₓ[4]*τᵖ_ᵘᵇ_ˢᴳ⁵[t] + Pˢᴳₘᵢₙ[4]*τᵖ_ˡᵇ_ˢᴳ⁵[t] 
                    - Qˢᴳₘₐₓ[4]*τᵠ_ᵘᵇ_ˢᴳ⁵[t] + Qˢᴳₘᵢₙ[4]*τᵠ_ˡᵇ_ˢᴳ⁵[t] 
                    + K_c_gc[1, 4]*ξ_gf_23[t] + K_c_gc[2, 4]*ξ_gf_24[t] + α_VSG[t]*γ[4,t]
                    + cˢᵗ[4]*σˢᴳ⁵_ˢᵗ[t] - cˢᵗ[4]*σˢᴳ⁵_ˢᵗ[t+1]
                    + sum(γ_ᵐⁱⁿ[13:14,t]) - sum(γ_ᵐᵃˣ[1,13:14,t]) + γ_ᵐⁱⁿ[3,t] + γ_ᵐⁱⁿ[7,t] + γ_ᵐⁱⁿ[10,t] - γ_ᵐᵃˣ[2,3,t] - γ_ᵐᵃˣ[2,7,t] - γ_ᵐᵃˣ[2,10,t] >= 0 )
    @constraint(model, cⁿˡ[5] + ψˢᴳ²⁷_ᵐᵃˣ[t] - Pˢᴳₘₐₓ[5]*τᵖ_ᵘᵇ_ˢᴳ²⁷[t] + Pˢᴳₘᵢₙ[5]*τᵖ_ˡᵇ_ˢᴳ²⁷[t] 
                    - Qˢᴳₘₐₓ[5]*τᵠ_ᵘᵇ_ˢᴳ²⁷[t] + Qˢᴳₘᵢₙ[5]*τᵠ_ˡᵇ_ˢᴳ²⁷[t] 
                    + K_c_gc[1, 5]*ξ_gf_23[t] + K_c_gc[2, 5]*ξ_gf_24[t] + α_VSG[t]*γ[5,t]
                    + cˢᵗ[5]*σˢᴳ²⁷_ˢᵗ[t] - cˢᵗ[5]*σˢᴳ²⁷_ˢᵗ[t+1]
                    + γ_ᵐⁱⁿ[15,t] + γ_ᵐⁱⁿ[4,t] + γ_ᵐⁱⁿ[8,t] + γ_ᵐⁱⁿ[11,t] + γ_ᵐⁱⁿ[13,t]   
                    - γ_ᵐᵃˣ[1,15,t] - γ_ᵐᵃˣ[2,4,t] - γ_ᵐᵃˣ[2,8,t] - γ_ᵐᵃˣ[2,11,t] - γ_ᵐᵃˣ[2,13,t] >= 0 )
    @constraint(model, cⁿˡ[6] + ψˢᴳ³⁰_ᵐᵃˣ[t] - Pˢᴳₘₐₓ[6]*τᵖ_ᵘᵇ_ˢᴳ³⁰[t] + Pˢᴳₘᵢₙ[6]*τᵖ_ˡᵇ_ˢᴳ³⁰[t] 
                    - Qˢᴳₘₐₓ[6]*τᵠ_ᵘᵇ_ˢᴳ³⁰[t] + Qˢᴳₘᵢₙ[6]*τᵠ_ˡᵇ_ˢᴳ³⁰[t] 
                    + K_c_gc[1, 6]*ξ_gf_23[t] + K_c_gc[2, 6]*ξ_gf_24[t] + α_VSG[t]*γ[6,t]
                    + cˢᵗ[6]*σˢᴳ³⁰_ˢᵗ[t] - cˢᵗ[6]*σˢᴳ³⁰_ˢᵗ[t+1]
                    + γ_ᵐⁱⁿ[15,t] + γ_ᵐⁱⁿ[5,t] + γ_ᵐⁱⁿ[9,t] + γ_ᵐⁱⁿ[12,t] + γ_ᵐⁱⁿ[14,t]   
                    - γ_ᵐᵃˣ[2,15,t] - γ_ᵐᵃˣ[2,5,t] - γ_ᵐᵃˣ[2,9,t] - γ_ᵐᵃˣ[2,12,t] - γ_ᵐᵃˣ[2,14,t] >= 0 )

end

    @constraint(model, cⁿˡ[1] + ψˢᴳ²_ᵐᵃˣ[T] - Pˢᴳₘₐₓ[1]*τᵖ_ᵘᵇ_ˢᴳ²[T] + Pˢᴳₘᵢₙ[1]*τᵖ_ˡᵇ_ˢᴳ²[T] 
                    - Qˢᴳₘₐₓ[1]*τᵠ_ᵘᵇ_ˢᴳ²[T] + Qˢᴳₘᵢₙ[1]*τᵠ_ˡᵇ_ˢᴳ²[T] 
                    + K_c_gc[1, 1]*ξ_gf_23[T] + K_c_gc[2, 1]*ξ_gf_24[T] + α_VSG[T]*γ[1,T]
                    + cˢᵗ[1]*σˢᴳ²_ˢᵗ[T] 
                    + sum(γ_ᵐⁱⁿ[1:5,T]) - sum(γ_ᵐᵃˣ[1,1:5,T]) >= 0 )
    @constraint(model, cⁿˡ[2] + ψˢᴳ³_ᵐᵃˣ[T] - Pˢᴳₘₐₓ[2]*τᵖ_ᵘᵇ_ˢᴳ³[T] + Pˢᴳₘᵢₙ[2]*τᵖ_ˡᵇ_ˢᴳ³[T] 
                    - Qˢᴳₘₐₓ[2]*τᵠ_ᵘᵇ_ˢᴳ³[T] + Qˢᴳₘᵢₙ[2]*τᵠ_ˡᵇ_ˢᴳ³[T] 
                    + K_c_gc[1, 2]*ξ_gf_23[T] + K_c_gc[2, 2]*ξ_gf_24[T] + α_VSG[T]*γ[2,T]
                    + cˢᵗ[2]*σˢᴳ³_ˢᵗ[T] 
                    + sum(γ_ᵐⁱⁿ[6:9,T]) - sum(γ_ᵐᵃˣ[1,6:9,T]) + γ_ᵐⁱⁿ[1,T] - γ_ᵐᵃˣ[2,1,T] >= 0 )
    @constraint(model, cⁿˡ[3] + ψˢᴳ⁴_ᵐᵃˣ[T] - Pˢᴳₘₐₓ[3]*τᵖ_ᵘᵇ_ˢᴳ⁴[T] + Pˢᴳₘᵢₙ[3]*τᵖ_ˡᵇ_ˢᴳ⁴[T] 
                    - Qˢᴳₘₐₓ[3]*τᵠ_ᵘᵇ_ˢᴳ⁴[T] + Qˢᴳₘᵢₙ[3]*τᵠ_ˡᵇ_ˢᴳ⁴[T] 
                    + K_c_gc[1, 3]*ξ_gf_23[T] + K_c_gc[2, 3]*ξ_gf_24[T] + α_VSG[T]*γ[3,T]
                    + cˢᵗ[3]*σˢᴳ⁴_ˢᵗ[T] 
                    + sum(γ_ᵐⁱⁿ[10:12,T]) - sum(γ_ᵐᵃˣ[1,10:12,T]) + γ_ᵐⁱⁿ[2,T] + γ_ᵐⁱⁿ[6,T] - γ_ᵐᵃˣ[2,2,T] - γ_ᵐᵃˣ[2,6,T] >= 0 )
    @constraint(model, cⁿˡ[4] + ψˢᴳ⁵_ᵐᵃˣ[T] - Pˢᴳₘₐₓ[4]*τᵖ_ᵘᵇ_ˢᴳ⁵[T] + Pˢᴳₘᵢₙ[4]*τᵖ_ˡᵇ_ˢᴳ⁵[T] 
                    - Qˢᴳₘₐₓ[4]*τᵠ_ᵘᵇ_ˢᴳ⁵[T] + Qˢᴳₘᵢₙ[4]*τᵠ_ˡᵇ_ˢᴳ⁵[T] 
                    + K_c_gc[1, 4]*ξ_gf_23[T] + K_c_gc[2, 4]*ξ_gf_24[T] + α_VSG[T]*γ[4,T]
                    + cˢᵗ[4]*σˢᴳ⁵_ˢᵗ[T] 
                    + sum(γ_ᵐⁱⁿ[13:14,T]) - sum(γ_ᵐᵃˣ[1,13:14,T]) + γ_ᵐⁱⁿ[3,T] + γ_ᵐⁱⁿ[7,T] + γ_ᵐⁱⁿ[10,T] - γ_ᵐᵃˣ[2,3,T] - γ_ᵐᵃˣ[2,7,T] - γ_ᵐᵃˣ[2,10,T] >= 0 )
    @constraint(model, cⁿˡ[5] + ψˢᴳ²⁷_ᵐᵃˣ[T] - Pˢᴳₘₐₓ[5]*τᵖ_ᵘᵇ_ˢᴳ²⁷[T] + Pˢᴳₘᵢₙ[5]*τᵖ_ˡᵇ_ˢᴳ²⁷[T] 
                    - Qˢᴳₘₐₓ[5]*τᵠ_ᵘᵇ_ˢᴳ²⁷[T] + Qˢᴳₘᵢₙ[5]*τᵠ_ˡᵇ_ˢᴳ²⁷[T] 
                    + K_c_gc[1, 5]*ξ_gf_23[T] + K_c_gc[2, 5]*ξ_gf_24[T] + α_VSG[T]*γ[5,T]
                    + cˢᵗ[5]*σˢᴳ²⁷_ˢᵗ[T] 
                    + γ_ᵐⁱⁿ[15,T] + γ_ᵐⁱⁿ[4,T] + γ_ᵐⁱⁿ[8,T] + γ_ᵐⁱⁿ[11,T] + γ_ᵐⁱⁿ[13,T]   
                    - γ_ᵐᵃˣ[1,15,T] - γ_ᵐᵃˣ[2,4,T] - γ_ᵐᵃˣ[2,8,T] - γ_ᵐᵃˣ[2,11,T] - γ_ᵐᵃˣ[2,13,T] >= 0 )
    @constraint(model, cⁿˡ[6] + ψˢᴳ³⁰_ᵐᵃˣ[T] - Pˢᴳₘₐₓ[6]*τᵖ_ᵘᵇ_ˢᴳ³⁰[T] + Pˢᴳₘᵢₙ[6]*τᵖ_ˡᵇ_ˢᴳ³⁰[T] 
                    - Qˢᴳₘₐₓ[6]*τᵠ_ᵘᵇ_ˢᴳ³⁰[T] + Qˢᴳₘᵢₙ[6]*τᵠ_ˡᵇ_ˢᴳ³⁰[T] 
                    + K_c_gc[1, 6]*ξ_gf_23[T] + K_c_gc[2, 6]*ξ_gf_24[T] + α_VSG[T]*γ[6,T]
                    + cˢᵗ[6]*σˢᴳ³⁰_ˢᵗ[T] 
                    + γ_ᵐⁱⁿ[15,T] + γ_ᵐⁱⁿ[5,T] + γ_ᵐⁱⁿ[9,T] + γ_ᵐⁱⁿ[12,T] + γ_ᵐⁱⁿ[14,T]   
                    - γ_ᵐᵃˣ[2,15,T] - γ_ᵐᵃˣ[2,5,T] - γ_ᵐᵃˣ[2,9,T] - γ_ᵐᵃˣ[2,12,T] - γ_ᵐᵃˣ[2,14,T] >= 0 )
                    
for t in 1:T
# ============================================================
# Cᵁ²...Cᵁ³⁰
    @constraint(model, 1 - σˢᴳ²_ˢᵗ[t] >= 0)
    @constraint(model, 1 - σˢᴳ³_ˢᵗ[t] >= 0)
    @constraint(model, 1 - σˢᴳ⁴_ˢᵗ[t] >= 0)
    @constraint(model, 1 - σˢᴳ⁵_ˢᵗ[t] >= 0)
    @constraint(model, 1 - σˢᴳ²⁷_ˢᵗ[t] >= 0)
    @constraint(model, 1 - σˢᴳ³⁰_ˢᵗ[t] >= 0)

# ============================================================
# ηₘ_1[1:15]
    for k in 1:15
        @constraint(model,
            K_c_m[1,k]*ξ_gf_23[t] + K_c_m[2,k]*ξ_gf_24[t] + γ_ᵐᵃˣ[1,k,t] + γ_ᵐᵃˣ[2,k,t] - γ_ᵐⁱⁿ[k,t] >= 0)
    end

# ============================================================
# ηₘ_2[1:6]
    for j in 1:6
        @constraint(model,
            K_c_m[1, 15+j]*ξ_gf_23[t] + K_c_m[2, 15+j]*ξ_gf_24[t] - γ[j,t] == 0)
    end

# ============================================================
# z_23, z_24
    @constraint(model,  -1/2*S_Base*μ_23[t] - ξ_gf_23[t] == 0)    # SOC sign
    @constraint(model,  -1/2*S_Base*μ_24[t] - ξ_gf_24[t] == 0)

# ============================================================
# Pᴳ[1]...Pᴳ[6]
    @constraint(model, cᵐ[1] + τᵖ_ᵘᵇ_ˢᴳ²[t] - τᵖ_ˡᵇ_ˢᴳ²[t] - υ_1_ˢᴳ²[t] - λᴱ[t] == 0)  # SOC sign
    @constraint(model, cᵐ[2] + τᵖ_ᵘᵇ_ˢᴳ³[t] - τᵖ_ˡᵇ_ˢᴳ³[t] - υ_1_ˢᴳ³[t] - λᴱ[t] == 0)
    @constraint(model, cᵐ[3] + τᵖ_ᵘᵇ_ˢᴳ⁴[t] - τᵖ_ˡᵇ_ˢᴳ⁴[t] - υ_1_ˢᴳ⁴[t] - λᴱ[t] == 0)
    @constraint(model, cᵐ[4] + τᵖ_ᵘᵇ_ˢᴳ⁵[t] - τᵖ_ˡᵇ_ˢᴳ⁵[t] - υ_1_ˢᴳ⁵[t] - λᴱ[t] == 0)
    @constraint(model, cᵐ[5] + τᵖ_ᵘᵇ_ˢᴳ²⁷[t] - τᵖ_ˡᵇ_ˢᴳ²⁷[t] - υ_1_ˢᴳ²⁷[t] - λᴱ[t] == 0)
    @constraint(model, cᵐ[6] + τᵖ_ᵘᵇ_ˢᴳ³⁰[t] - τᵖ_ˡᵇ_ˢᴳ³⁰[t] - υ_1_ˢᴳ³⁰[t] - λᴱ[t] == 0)

# ============================================================
# Pᴳ[7] VSG_1，Pᴳ[8] IBG_23，Pᴳ[9] IBG_24
    @constraint(model, τᵖ_ᵘᵇ_VSG_1[t] - υ_1_VSG_1[t] - λᴱ[t] == 0)   # SOC sign

    @constraint(model, τᵖ_ᵘᵇ_IBG_23[t] - υ_1_IBG_23[t] - λ_1_23[t] - λᴱ[t] == 0)
    @constraint(model, τᵖ_ᵘᵇ_IBG_24[t] - υ_1_IBG_24[t] - λ_1_24[t] - λᴱ[t] == 0)

# ============================================================
# Qᴳ[1]...Qᴳ[6] 
    @constraint(model, τᵠ_ᵘᵇ_ˢᴳ²[t] - τᵠ_ˡᵇ_ˢᴳ²[t] - υ_2_ˢᴳ²[t] - ϕ[t] == 0)  # SOC sign
    @constraint(model, τᵠ_ᵘᵇ_ˢᴳ³[t] - τᵠ_ˡᵇ_ˢᴳ³[t] - υ_2_ˢᴳ³[t] - ϕ[t] == 0)
    @constraint(model, τᵠ_ᵘᵇ_ˢᴳ⁴[t] - τᵠ_ˡᵇ_ˢᴳ⁴[t] - υ_2_ˢᴳ⁴[t] - ϕ[t] == 0)
    @constraint(model, τᵠ_ᵘᵇ_ˢᴳ⁵[t] - τᵠ_ˡᵇ_ˢᴳ⁵[t] - υ_2_ˢᴳ⁵[t] - ϕ[t] == 0)
    @constraint(model, τᵠ_ᵘᵇ_ˢᴳ²⁷[t] - τᵠ_ˡᵇ_ˢᴳ²⁷[t] - υ_2_ˢᴳ²⁷[t] - ϕ[t] == 0)
    @constraint(model, τᵠ_ᵘᵇ_ˢᴳ³⁰[t] - τᵠ_ˡᵇ_ˢᴳ³⁰[t] - υ_2_ˢᴳ³⁰[t] - ϕ[t] == 0)

# ============================================================
# Qᴳ[7] VSG_1，Qᴳ[8] IBG_23，Qᴳ[9] IBG_24
    @constraint(model, τᵠ_ᵘᵇ_VSG_1[t] - τᵠ_ˡᵇ_VSG_1[t] - υ_2_VSG_1[t] - ϕ[t] == 0)  # SOC sign

    @constraint(model, τᵠ_ᵘᵇ_IBG_23[t] - τᵠ_ˡᵇ_IBG_23[t] - υ_2_IBG_23[t]
                   - λ_2_23[t] - μ_23[t] - ϕ[t] == 0)
    @constraint(model, τᵠ_ᵘᵇ_IBG_24[t] - τᵠ_ˡᵇ_IBG_24[t] - υ_2_IBG_24[t]
                   - λ_2_24[t] - μ_24[t] - ϕ[t] == 0)

# ============================================================
# SOC
# ============================================================
    @constraint(model, [υ_3_ˢᴳ²[t], υ_1_ˢᴳ²[t], υ_2_ˢᴳ²[t]] in SecondOrderCone())
    @constraint(model, [υ_3_ˢᴳ³[t], υ_1_ˢᴳ³[t], υ_2_ˢᴳ³[t]] in SecondOrderCone())
    @constraint(model, [υ_3_ˢᴳ⁴[t], υ_1_ˢᴳ⁴[t], υ_2_ˢᴳ⁴[t]] in SecondOrderCone())
    @constraint(model, [υ_3_ˢᴳ⁵[t], υ_1_ˢᴳ⁵[t], υ_2_ˢᴳ⁵[t]] in SecondOrderCone())
    @constraint(model, [υ_3_ˢᴳ²⁷[t], υ_1_ˢᴳ²⁷[t], υ_2_ˢᴳ²⁷[t]] in SecondOrderCone())
    @constraint(model, [υ_3_ˢᴳ³⁰[t], υ_1_ˢᴳ³⁰[t], υ_2_ˢᴳ³⁰[t]] in SecondOrderCone())
    @constraint(model, [υ_3_VSG_1[t], υ_1_VSG_1[t], υ_2_VSG_1[t]] in SecondOrderCone())
    @constraint(model, [υ_3_IBG_23[t], υ_1_IBG_23[t], υ_2_IBG_23[t]] in SecondOrderCone())
    @constraint(model, [υ_3_IBG_24[t], υ_1_IBG_24[t], υ_2_IBG_24[t]] in SecondOrderCone())

    @constraint(model, [μ_23[t], λ_1_23[t], λ_2_23[t]] in SecondOrderCone())
    @constraint(model, [μ_24[t], λ_1_24[t], λ_2_24[t]] in SecondOrderCone())

end

# ============================================================
# Max objective
@variable(model, obj_dual[1:T])
@variable(model, obj_dual_final)

for t in 1:T
    @constraint(model, obj_dual[t] == λᴱ[t]*Pᴰ[t] + ϕ[t]*Qᴰ[t]
        - sum(γ_ᵐⁱⁿ[:,t]) 
        - ( ψˢᴳ²_ᵐᵃˣ[t] + ψˢᴳ³_ᵐᵃˣ[t] + ψˢᴳ⁴_ᵐᵃˣ[t] + ψˢᴳ⁵_ᵐᵃˣ[t] + ψˢᴳ²⁷_ᵐᵃˣ[t] + ψˢᴳ³⁰_ᵐᵃˣ[t] )
  
        - ( υ_3_ˢᴳ²[t]*Sᵐᵃˣ_gc[1]  + υ_3_ˢᴳ³[t]*Sᵐᵃˣ_gc[2]  + υ_3_ˢᴳ⁴[t]*Sᵐᵃˣ_gc[3]
        + υ_3_ˢᴳ⁵[t]*Sᵐᵃˣ_gc[4]  + υ_3_ˢᴳ²⁷[t]*Sᵐᵃˣ_gc[5] + υ_3_ˢᴳ³⁰[t]*Sᵐᵃˣ_gc[6]
        + υ_3_VSG_1[t]*Sᵐᵃˣ_gv[1] + υ_3_IBG_23[t]*Sᵐᵃˣ_c[1] + υ_3_IBG_24[t]*Sᵐᵃˣ_c[2] ) 

        + ( ξ_gf_23[t] * α_VSG[t]*K_c_gv[1,1] + ξ_gf_24[t] * α_VSG[t]*K_c_gv[2,1] )

        - ( α_VSG[t]*Pⱽˢᴳₘₐₓ[1]*τᵖ_ᵘᵇ_VSG_1[t] + α_IBG[1,t]*Pᴵᴮᴳₘₐₓ[1]*τᵖ_ᵘᵇ_IBG_23[t] + α_IBG[2,t]*Pᴵᴮᴳₘₐₓ[2]*τᵖ_ᵘᵇ_IBG_24[t]
        + Qⱽˢᴳₘₐₓ[1]*τᵠ_ᵘᵇ_VSG_1[t] - Qⱽˢᴳₘᵢₙ[1]*τᵠ_ˡᵇ_VSG_1[t] + ratio*Qᴵᴮᴳₘₐₓ[1]*τᵠ_ᵘᵇ_IBG_23[t] - ratio*Qᴵᴮᴳₘᵢₙ[1]*τᵠ_ˡᵇ_IBG_23[t]
        + ratio*Qᴵᴮᴳₘₐₓ[2]*τᵠ_ᵘᵇ_IBG_24[t] - ratio*Qᴵᴮᴳₘᵢₙ[2]*τᵠ_ˡᵇ_IBG_24[t] ) 
                )
end

@constraint(model, obj_dual_final == sum(obj_dual[t] for t in 1:T)
                                    - σˢᴳ²_ˢᵗ[1]  * cˢᵗ[1]*yˢᴳ²_0
                                    - σˢᴳ³_ˢᵗ[1]  * cˢᵗ[2]*yˢᴳ³_0
                                    - σˢᴳ⁴_ˢᵗ[1]  * cˢᵗ[3]*yˢᴳ⁴_0
                                    - σˢᴳ⁵_ˢᵗ[1]  * cˢᵗ[4]*yˢᴳ⁵_0
                                    - σˢᴳ²⁷_ˢᵗ[1] * cˢᵗ[5]*yˢᴳ²⁷_0
                                    - σˢᴳ³⁰_ˢᵗ[1] * cˢᵗ[6]*yˢᴳ³⁰_0  )


# ============================================================
# Nonnegative profit constraint

@variable(model, profit_SG[1:6,1:T] )           # Lower bound for energy price 

@constraint(model, λᴱ.>= 0  )           # Lower bound for energy price 
@constraint(model, λᴱ.<= 2 * maximum(cᵐ)  )           # Upper bound for energy price 

@constraint(model, μ_23.>= 0  )           
@constraint(model, μ_23.<= 1 * (maximum(cⁿˡ) + maximum(cˢᵗ) )  )           

@constraint(model, μ_24.>= 0  )           
@constraint(model, μ_24.<= 1 * (maximum(cⁿˡ) + maximum(cˢᵗ) )  )  

@variable(model, energy_profit_SG[1:6,1:T] )  
@variable(model, vs_revenue_SG[1:6,1:T] )  

@variable(model, profit_wts[1:3,1:T] )  
@variable(model, profit_wts_energy[1:3,1:T] ) 
@variable(model, profit_wts_vs[1:3,1:T] ) 

for t in 1:T

@constraint(model, profit_wts_energy[1,t] == λᴱ[t] * Pᴳ[7,t] )

@constraint(model, profit_wts_energy[2,t] == λᴱ[t] * Pᴳ[8,t] )

@constraint(model, profit_wts_energy[3,t] == λᴱ[t] * Pᴳ[9,t] )

@constraint(model, profit_wts_vs[1,t] == μ_23[t]* ( α_VSG[t]*K_c_gv[1, 1] + sum(K_c_m[1, 16:end] .*ηₘ_2[:,t]   )) *S_Base
                                        + μ_24[t]* ( α_VSG[t]*K_c_gv[2, 1] + sum(K_c_m[2, 16:end] .*ηₘ_2[:,t]   )) *S_Base )

@constraint(model, profit_wts_vs[2,t] == (-λ_2_23[t] + μ_23[t]  )* Qᴳ[8,t] )

@constraint(model, profit_wts_vs[3,t] == (-λ_2_24[t] + μ_24[t]  )* Qᴳ[9,t] )

@constraint(model, profit_wts[:,t] == profit_wts_energy[:,t] +  profit_wts_vs[:,t] )

@constraint(model,  profit_SG[1,t] == energy_profit_SG[1,t] + vs_revenue_SG[1,t] )
@constraint(model,  profit_SG[2,t] == energy_profit_SG[2,t] + vs_revenue_SG[2,t] )
@constraint(model,  profit_SG[3,t] == energy_profit_SG[3,t] + vs_revenue_SG[3,t] )
@constraint(model,  profit_SG[4,t] == energy_profit_SG[4,t] + vs_revenue_SG[4,t] )
@constraint(model,  profit_SG[5,t] == energy_profit_SG[5,t] + vs_revenue_SG[5,t] )
@constraint(model,  profit_SG[6,t] == energy_profit_SG[6,t] + vs_revenue_SG[6,t] )


@constraint(model,  energy_profit_SG[1,t] == λᴱ[t] * Pᴳ[1,t] - cⁿˡ[1]*yˢᴳ²[t] - cᵐ[1]*Pᴳ[1,t] - Cᵁ²[t] )
@constraint(model,  vs_revenue_SG[1,t] == μ_23[t]*K_c_gc[1, 1]*S_Base *yˢᴳ²[t] + μ_24[t]*K_c_gc[2, 1]*S_Base *yˢᴳ²[t]
                                        + μ_23[t]*sum(K_c_m[1, 1:5].*ηₘ_1[1:5,t])*S_Base + μ_24[t]*sum(K_c_m[2, 1:5].*ηₘ_1[1:5,t])*S_Base 
                                        + μ_23[t]*ηₘ_2[1,t]*K_c_m[1, 16]*S_Base + μ_24[t]*ηₘ_2[1,t]*K_c_m[2, 16]*S_Base  )
    
@constraint(model,  energy_profit_SG[2,t] == λᴱ[t] * Pᴳ[2,t] - cⁿˡ[2]*yˢᴳ³[t] - cᵐ[2]*Pᴳ[2,t] - Cᵁ³[t] )
@constraint(model,  vs_revenue_SG[2,t] == μ_23[t]*K_c_gc[1, 2]*S_Base *yˢᴳ³[t] + μ_24[t]*K_c_gc[2, 2]*S_Base *yˢᴳ³[t]
                                        + μ_23[t]*sum(K_c_m[1, 6:9].*ηₘ_1[6:9,t])*S_Base + μ_24[t]*sum(K_c_m[2, 6:9].*ηₘ_1[6:9,t])*S_Base
                                        + μ_23[t]*ηₘ_1[1,t]*K_c_m[1, 1]*S_Base + μ_24[t]*ηₘ_1[1,t]*K_c_m[2, 1]*S_Base 
                                        + μ_23[t]*ηₘ_2[2,t]*K_c_m[1, 17]*S_Base + μ_24[t]*ηₘ_2[2,t]*K_c_m[2, 17]*S_Base  )

@constraint(model,  energy_profit_SG[3,t] == λᴱ[t] * Pᴳ[3,t] - cⁿˡ[3]*yˢᴳ⁴[t] - cᵐ[3]*Pᴳ[3,t] - Cᵁ⁴[t] )
@constraint(model,  vs_revenue_SG[3,t] == μ_23[t]*K_c_gc[1, 3]*S_Base *yˢᴳ⁴[t] + μ_24[t]*K_c_gc[2, 3]*S_Base *yˢᴳ⁴[t]
                                        + μ_23[t]*sum(K_c_m[1, 10:12].*ηₘ_1[10:12,t])*S_Base + μ_24[t]*sum(K_c_m[2, 10:12].*ηₘ_1[10:12,t])*S_Base
                                        + μ_23[t]*ηₘ_1[2,t]*K_c_m[1, 2]*S_Base + μ_23[t]*ηₘ_1[6,t]*K_c_m[1, 6]*S_Base 
                                        + μ_24[t]*ηₘ_1[2,t]*K_c_m[2, 2]*S_Base + μ_24[t]*ηₘ_1[6,t]*K_c_m[2, 6]*S_Base 
                                        + μ_23[t]*ηₘ_2[3,t]*K_c_m[1, 18]*S_Base + μ_24[t]*ηₘ_2[3,t]*K_c_m[2, 18]*S_Base  )

@constraint(model,  energy_profit_SG[4,t] == λᴱ[t] * Pᴳ[4,t] - cⁿˡ[4]*yˢᴳ⁵[t] - cᵐ[4]*Pᴳ[4,t] - Cᵁ⁵[t] )   
@constraint(model,  vs_revenue_SG[4,t] == μ_23[t]*K_c_gc[1, 4]*S_Base *yˢᴳ⁵[t] + μ_24[t]*K_c_gc[2, 4]*S_Base *yˢᴳ⁵[t]
                                        + μ_23[t]*sum(K_c_m[1, 13:14].*ηₘ_1[13:14,t])*S_Base + μ_24[t]*sum(K_c_m[2, 13:14].*ηₘ_1[13:14,t])*S_Base
                                        + μ_23[t]*ηₘ_1[3,t]*K_c_m[1, 3]*S_Base + μ_23[t]*ηₘ_1[7,t]*K_c_m[1, 7]*S_Base + μ_23[t]*ηₘ_1[10,t]*K_c_m[1, 10]*S_Base 
                                        + μ_24[t]*ηₘ_1[3,t]*K_c_m[2, 3]*S_Base + μ_24[t]*ηₘ_1[7,t]*K_c_m[2, 7]*S_Base + μ_24[t]*ηₘ_1[10,t]*K_c_m[2, 10]*S_Base 
                                        + μ_23[t]*ηₘ_2[4,t]*K_c_m[1, 19]*S_Base + μ_24[t]*ηₘ_2[4,t]*K_c_m[2, 19]*S_Base )

@constraint(model,  energy_profit_SG[5,t] == λᴱ[t] * Pᴳ[5,t] - cⁿˡ[5]*yˢᴳ²⁷[t] - cᵐ[5]*Pᴳ[5,t] - Cᵁ²⁷[t] )
@constraint(model,  vs_revenue_SG[5,t] == μ_23[t]*K_c_gc[1, 5]*S_Base *yˢᴳ²⁷[t] + μ_24[t]*K_c_gc[2, 5]*S_Base *yˢᴳ²⁷[t]
                                        + μ_23[t]*K_c_m[1, 15].*ηₘ_1[15,t]*S_Base + μ_24[t]*K_c_m[2, 15].*ηₘ_1[15,t]*S_Base
                                        + μ_23[t]*ηₘ_1[4,t]*K_c_m[1, 4]*S_Base + μ_23[t]*ηₘ_1[8,t]*K_c_m[1, 8]*S_Base + μ_23[t]*ηₘ_1[11,t]*K_c_m[1, 11]*S_Base + μ_23[t]*ηₘ_1[13,t]*K_c_m[1, 13]*S_Base 
                                        + μ_24[t]*ηₘ_1[4,t]*K_c_m[2, 4]*S_Base + μ_24[t]*ηₘ_1[8,t]*K_c_m[2, 8]*S_Base + μ_24[t]*ηₘ_1[11,t]*K_c_m[2, 11]*S_Base + μ_24[t]*ηₘ_1[13,t]*K_c_m[2, 13]*S_Base 
                                        + μ_23[t]*ηₘ_2[5,t]*K_c_m[1, 20]*S_Base + μ_24[t]*ηₘ_2[5,t]*K_c_m[2, 20]*S_Base )

@constraint(model,  energy_profit_SG[6,t] == λᴱ[t] * Pᴳ[6,t] - cⁿˡ[6]*yˢᴳ³⁰[t] - cᵐ[6]*Pᴳ[6,t] - Cᵁ³⁰[t] )                                      
@constraint(model,  vs_revenue_SG[6,t] == μ_23[t]*K_c_gc[1, 6]*S_Base *yˢᴳ³⁰[t] + μ_24[t]*K_c_gc[2, 6]*S_Base *yˢᴳ³⁰[t]
                                        + μ_23[t]*ηₘ_1[5,t]*K_c_m[1, 5]*S_Base + μ_24[t]*ηₘ_1[5,t]*K_c_m[2, 5]*S_Base
                                        + μ_23[t]*ηₘ_1[9,t]*K_c_m[1, 9]*S_Base + μ_23[t]*ηₘ_1[12,t]*K_c_m[1, 12]*S_Base + μ_23[t]*ηₘ_1[14,t]*K_c_m[1, 14]*S_Base 
                                        + μ_24[t]*ηₘ_1[9,t]*K_c_m[2, 9]*S_Base + μ_24[t]*ηₘ_1[12,t]*K_c_m[2, 12]*S_Base + μ_24[t]*ηₘ_1[14,t]*K_c_m[2, 14]*S_Base 
                                        + μ_23[t]*ηₘ_1[15,t]*K_c_m[1, 15]*S_Base + μ_24[t]*ηₘ_1[15,t]*K_c_m[2, 15]*S_Base
                                        + μ_23[t]*ηₘ_2[6,t]*K_c_m[1, 21]*S_Base + μ_24[t]*ηₘ_2[6,t]*K_c_m[2, 21]*S_Base )

#@constraint(model, profit_SG[:,t] .>= 0)

end

@constraint(model, sum(profit_wts[1,:]) >= LCOE_gv)
@constraint(model, sum(profit_wts[2,:]) >= LCOE_gf[1])
@constraint(model, sum(profit_wts[3,:]) >= LCOE_gf[2])

@constraint(model, sum(profit_SG[1,:]) >= 0)
@constraint(model, sum(profit_SG[2,:]) >= 0)
@constraint(model, sum(profit_SG[3,:]) >= 0)
@constraint(model, sum(profit_SG[4,:]) >= 0)
@constraint(model, sum(profit_SG[5,:]) >= 0)
@constraint(model, sum(profit_SG[6,:]) >= 0)


@objective(model, Max, obj_dual_final)  # single-level objective function
#-------Solve and Output Results
set_optimizer(model,  Gurobi.Optimizer)   
optimize!(model)


#-------Analyze results


λᴱ_value = value.(λᴱ)
μ_23_value = value.(μ_23)
μ_24_value = value.(μ_24)
λ_2_23_value = value.(λ_2_23)
λ_2_24_value = value.(λ_2_24)
Pᴳ_value = value.(Pᴳ)

energy_profit_SG = value.(energy_profit_SG)
vs_revenue_SG = value.(vs_revenue_SG)
profit_SG = value.(profit_SG)

profit_wts_energy = value.(profit_wts_energy)
profit_wts_vs = value.(profit_wts_vs)

return λᴱ_value, μ_23_value, μ_24_value, μ_23_value -λ_2_23_value, μ_24_value -λ_2_24_value, energy_profit_SG, vs_revenue_SG, profit_SG, profit_wts_energy, profit_wts_vs

end