

#-----------------
# Define model
#-----------------

model= Model()                     

@variable(model, yˢᴳ²[1:T])                         # status of SGs, buses:2,3,4,5,27,30.    
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



#-----------------------------AC power flow constraints -------------------------------------------------------------------

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

Power_balance_P = [Dict{Int, ConstraintRef}() for _ in 1:bus_num]
Power_balance_Q = [Dict{Int, ConstraintRef}() for _ in 1:bus_num]

for i in 1:bus_num
    for t in 1:T

        Power_balance_P[i][t] = @constraint(model, Pᴳ[i,t] - Pᴰ[i,t]*1.1 ==  sum( Pˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # p_i^g-p_i^d == G_{i,i}*c_{i,i}+...   eq. 50(a)
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




#--------------------------------Voltage stability constraints-----------

ηₘ_1=zeros(15,T)        # product of each pair of SG 
ηₘ_2=zeros(6,T)         # product of each pair of SG and VSG 

for t in 1:T
    # 15+6
    ηₘ_1[1,t] = yˢᴳ²_opt[t] *yˢᴳ³_opt[t]  
    ηₘ_1[2,t] = yˢᴳ²_opt[t] *yˢᴳ⁴_opt[t]
    ηₘ_1[3,t] = yˢᴳ²_opt[t] *yˢᴳ⁵_opt[t]
    ηₘ_1[4,t] = yˢᴳ²_opt[t] *yˢᴳ²⁷_opt[t]
    ηₘ_1[5,t] = yˢᴳ²_opt[t] *yˢᴳ³⁰_opt[t]
    ηₘ_1[6,t] = yˢᴳ³_opt[t] *yˢᴳ⁴_opt[t]
    ηₘ_1[7,t] = yˢᴳ³_opt[t] *yˢᴳ⁵_opt[t]
    ηₘ_1[8,t] = yˢᴳ³_opt[t] *yˢᴳ²⁷_opt[t]
    ηₘ_1[9,t] = yˢᴳ³_opt[t] *yˢᴳ³⁰_opt[t]
    ηₘ_1[10,t] = yˢᴳ⁴_opt[t] *yˢᴳ⁵_opt[t]
    ηₘ_1[11,t] = yˢᴳ⁴_opt[t] *yˢᴳ²⁷_opt[t]
    ηₘ_1[12,t] = yˢᴳ⁴_opt[t] *yˢᴳ³⁰_opt[t]
    ηₘ_1[13,t] = yˢᴳ⁵_opt[t] *yˢᴳ²⁷_opt[t]
    ηₘ_1[14,t] = yˢᴳ⁵_opt[t] *yˢᴳ³⁰_opt[t]
    ηₘ_1[15,t] = yˢᴳ²⁷_opt[t] *yˢᴳ³⁰_opt[t]

    ηₘ_2[1,t] =yˢᴳ²_opt[t]*α_VSG
    ηₘ_2[2,t] =yˢᴳ³_opt[t]*α_VSG
    ηₘ_2[3,t] =yˢᴳ⁴_opt[t]*α_VSG
    ηₘ_2[4,t] =yˢᴳ⁵_opt[t]*α_VSG
    ηₘ_2[5,t] =yˢᴳ²⁷_opt[t]*α_VSG
    ηₘ_2[6,t] =yˢᴳ³⁰_opt[t]*α_VSG

end


@variable(model, P_23_eq[1:T])              # equivalent P and Q of IBGs, buses:23,24    and Γ_c  
@variable(model, P_24_eq[1:T]) 
@variable(model, Q_23_eq[1:T])                 
@variable(model, Q_24_eq[1:T])
@variable(model, Γ_23[1:T])
@variable(model, Γ_24[1:T])

Voltage_constraint_23=Dict()
Voltage_constraint_24=Dict()

for t in 1:T

    @constraint(model, P_23_eq[t]== Pᴳ[23,t]+  ( sh_1- ( K_c_gc[1,1] *yˢᴳ²_opt[t] +K_c_gc[1,2] *yˢᴳ³_opt[t] +K_c_gc[1,3] *yˢᴳ⁴_opt[t] 
    +K_c_gc[1,4] *yˢᴳ⁵_opt[t] +K_c_gc[1,5] *yˢᴳ²⁷_opt[t] +K_c_gc[1,6] *yˢᴳ³⁰_opt[t]    +α_VSG*K_c_gv[1,1] 
    + sum(K_c_m[1,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[1,16:end] .*ηₘ_2[:,t] ) ) )    *Pᴳ[24,t])

    @constraint(model, P_24_eq[t]== Pᴳ[24,t]+  ( sh_2- ( K_c_gc[2,1] *yˢᴳ²_opt[t] +K_c_gc[2,2] *yˢᴳ³_opt[t] +K_c_gc[2,3] *yˢᴳ⁴_opt[t] 
    +K_c_gc[2,4] *yˢᴳ⁵_opt[t] +K_c_gc[2,5] *yˢᴳ²⁷_opt[t] +K_c_gc[2,6] *yˢᴳ³⁰_opt[t]    +α_VSG*K_c_gv[2,1]
    + sum(K_c_m[2,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[2,16:end] .*ηₘ_2[:,t] ) ) )    *Pᴳ[23,t])

    @constraint(model, Q_23_eq[t]== Qᴳ[23,t]+  ( sh_1- ( K_c_gc[1,1] *yˢᴳ²_opt[t] +K_c_gc[1,2] *yˢᴳ³_opt[t] +K_c_gc[1,3] *yˢᴳ⁴_opt[t] 
    +K_c_gc[1,4] *yˢᴳ⁵_opt[t] +K_c_gc[1,5] *yˢᴳ²⁷_opt[t] +K_c_gc[1,6] *yˢᴳ³⁰_opt[t]    +α_VSG*K_c_gv[1,1] 
    + sum(K_c_m[1,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[1,16:end] .*ηₘ_2[:,t] ) ) )    *Qᴳ[24,t])

    @constraint(model, Q_24_eq[t]== Qᴳ[24,t]+  ( sh_2- ( K_c_gc[2,1] *yˢᴳ²_opt[t] +K_c_gc[2,2] *yˢᴳ³_opt[t] +K_c_gc[2,3] *yˢᴳ⁴_opt[t] 
    +K_c_gc[2,4] *yˢᴳ⁵_opt[t] +K_c_gc[2,5] *yˢᴳ²⁷_opt[t] +K_c_gc[2,6] *yˢᴳ³⁰_opt[t]    +α_VSG*K_c_gv[2,1]
    + sum(K_c_m[2,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[2,16:end] .*ηₘ_2[:,t] ) ) )    *Qᴳ[23,t])

    @constraint(model, Γ_23[t]== 1/2 * ( K_c_gc[3,1] *yˢᴳ²_opt[t] +K_c_gc[3,2] *yˢᴳ³_opt[t] +K_c_gc[3,3] *yˢᴳ⁴_opt[t] 
    +K_c_gc[3,4] *yˢᴳ⁵_opt[t] +K_c_gc[3,5] *yˢᴳ²⁷_opt[t] +K_c_gc[3,6] *yˢᴳ³⁰_opt[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[3,16:end] .*ηₘ_2[:,t]   ) )

    @constraint(model, Γ_24[t]== 1/2 * ( K_c_gc[4,1] *yˢᴳ²_opt[t] +K_c_gc[4,2] *yˢᴳ³_opt[t] +K_c_gc[4,3] *yˢᴳ⁴_opt[t] 
    +K_c_gc[4,4] *yˢᴳ⁵_opt[t] +K_c_gc[4,5] *yˢᴳ²⁷_opt[t] +K_c_gc[4,6] *yˢᴳ³⁰_opt[t] )   +α_VSG*K_c_gv[4,1]
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

@objective(model, Min, obj)  # single-level objective function
#-------Solve and Output Results
set_optimizer(model , Gurobi.Optimizer)
set_optimizer_attribute(model, "QCPDual", 1)
optimize!(model)



#------------------------------
# Analyze results
#------------------------------
energy_prices_P=zeros(bus_num,T)
energy_prices_P = [dual(Power_balance_P[i][t]) for i in 1:bus_num, t in 1:T]

mean_energy_prices= sum(energy_prices_P[1,:]) / T 


bar(energy_prices_P')

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
# Calculate commitment revenues of units
#------------------------------
UC_price_yˢᴳ²= [dual(price_yˢᴳ²[t]) for t in 1:T]
UC_price_yˢᴳ³= [dual(price_yˢᴳ³[t]) for t in 1:T]
UC_price_yˢᴳ⁴= [dual(price_yˢᴳ⁴[t]) for t in 1:T]
UC_price_yˢᴳ⁵= [dual(price_yˢᴳ⁵[t]) for t in 1:T]
UC_price_yˢᴳ²⁷= [dual(price_yˢᴳ²⁷[t]) for t in 1:T]
UC_price_yˢᴳ³⁰= [dual(price_yˢᴳ³⁰[t]) for t in 1:T]

yˢᴳ²=JuMP.value.(yˢᴳ²)
yˢᴳ³=JuMP.value.(yˢᴳ³)
yˢᴳ⁴=JuMP.value.(yˢᴳ⁴)
yˢᴳ⁵=JuMP.value.(yˢᴳ⁵)
yˢᴳ²⁷=JuMP.value.(yˢᴳ²⁷)
yˢᴳ³⁰=JuMP.value.(yˢᴳ³⁰)

UC_price_each_SG[1,a]=sum(UC_price_yˢᴳ² .*yˢᴳ²)
UC_price_each_SG[2,a]=sum(UC_price_yˢᴳ³ .*yˢᴳ³)
UC_price_each_SG[3,a]=sum(UC_price_yˢᴳ⁴ .*yˢᴳ⁴)
UC_price_each_SG[4,a]=sum(UC_price_yˢᴳ⁵ .*yˢᴳ⁵)
UC_price_each_SG[5,a]=sum(UC_price_yˢᴳ²⁷ .*yˢᴳ²⁷)
UC_price_each_SG[6,a]=sum(UC_price_yˢᴳ³⁰ .*yˢᴳ³⁰)

profit_each_SG[1,a]=profitˢᴳ² +UC_price_each_SG[1,a]
profit_each_SG[2,a]=profitˢᴳ³ +UC_price_each_SG[2,a]
profit_each_SG[3,a]=profitˢᴳ⁴ +UC_price_each_SG[3,a]
profit_each_SG[4,a]=profitˢᴳ⁵ +UC_price_each_SG[4,a]
profit_each_SG[5,a]=profitˢᴳ²⁷ +UC_price_each_SG[5,a]
profit_each_SG[6,a]=profitˢᴳ³⁰ +UC_price_each_SG[6,a]

println(UC_price_yˢᴳ³⁰ .*yˢᴳ³⁰)

plot(UC_price_yˢᴳ² .*yˢᴳ²)
plot!(UC_price_yˢᴳ³ .*yˢᴳ³)
plot!(UC_price_yˢᴳ⁴ .*yˢᴳ⁴)
plot!(UC_price_yˢᴳ⁵ .*yˢᴳ⁵)
plot!(UC_price_yˢᴳ²⁷ .*yˢᴳ²⁷)
plot!(UC_price_yˢᴳ³⁰ .*yˢᴳ³⁰)




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

bar(-lambda_1[:,1])
bar!(-lambda_1[:,2])

plot!(-lambda_2[:,1]+mu[:,1])
plot!(-lambda_2[:,2]+mu[:,2])


sum(mu[:,1])/T
sum(mu[:,2])/T
sum(-lambda_1[:,1])/T
sum(-lambda_1[:,2])/T
sum(-lambda_2[:,1]+mu[:,1])/T
sum(-lambda_2[:,2]+mu[:,2])/T

println(-lambda_2[:,1]+mu[:,1])

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

#P_23_eqivalent_list=zeros(T)
#P_24_eqivalent_list=zeros(T)
#Q_23_eqivalent_list=zeros(T)
#Q_24_eqivalent_list=zeros(T)

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

sum(revenue_IBG_23_give_active)
sum(revenue_IBG_23_give_reactive)

sum(revenue_IBG_24_give_active)
sum(revenue_IBG_24_give_reactive)




#----revenue for (V)SGs
VS_revenue_VSG_1=zeros(2,T)
VS_revenue_SGs_23=zeros(6,T)
VS_revenue_SGs_24=zeros(6,T)

for i in [1,2,3,4,5,6]
    for t in 1:T
        α_VSG=0.9
        UC_schedule[:,t]=[yˢᴳ²[t];yˢᴳ³[t];yˢᴳ⁴[t];yˢᴳ⁵[t];yˢᴳ²⁷[t];yˢᴳ³⁰[t]]
        zᵟ_23_1, zᵟ_24_1 =  calculate_Tao(α_VSG,Gᵥ,IBG,UC_schedule[:,t]) 
        #UC_schedule[i,t] = 0
        α_VSG=0
        zᵟ_23_1_delta, zᵟ_24_1_delta = calculate_delta_gamma(α_VSG,Gᵥ,IBG,UC_schedule[:,t])
        #VS_revenue_SGs_23[i,t] = mu[t,1]* (zᵟ_23_1 - zᵟ_23_1_delta) *100/2
        #VS_revenue_SGs_24[i,t] = mu[t,2]* (zᵟ_24_1 - zᵟ_24_1_delta) *100/2
        VS_revenue_VSG_1[1,t] =  mu[t,1]* (zᵟ_23_1 - zᵟ_23_1_delta) *100/2
        VS_revenue_VSG_1[2,t] =  mu[t,2]* (zᵟ_24_1 - zᵟ_24_1_delta) *100/2
    end
end

sum(VS_revenue_SGs_23[6,:]) + sum(VS_revenue_SGs_24[6,:])
sum(VS_revenue_VSG_1[1,:]) + sum(VS_revenue_VSG_1[2,:])

plot( VS_revenue_VSG_1[1,:]+ VS_revenue_VSG_1[2,:] )
plot!( VS_revenue_SGs_23[1,:]+ VS_revenue_SGs_24[1,:] )
plot!( VS_revenue_SGs_23[2,:]+ VS_revenue_SGs_24[2,:] )
plot!( VS_revenue_SGs_23[3,:]+ VS_revenue_SGs_24[3,:] )
plot!( VS_revenue_SGs_23[4,:]+ VS_revenue_SGs_24[4,:] )
plot!( VS_revenue_SGs_23[5,:]+ VS_revenue_SGs_24[5,:] )
plot!( VS_revenue_SGs_23[6,:]+ VS_revenue_SGs_24[6,:] )


println(  VS_revenue_VSG_1[1,:]+ VS_revenue_VSG_1[2,:])