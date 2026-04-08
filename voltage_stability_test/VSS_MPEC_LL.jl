# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Pricing voltage stability
# 4.July.2025

function LL(K_c_gc, K_c_gv, K_c_m, sh_1, sh_2, bus_num, T, Pᴰ, Qᴰ, S_lines, Pˢᴳₘₐₓ, Pˢᴳₘᵢₙ   )
    
end


#-------------------------------------------------------
#-------------------Define model-------------------
#-------------------------------------------------------



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


for i in 1:bus_num
    for t in 1:T

        @constraint(model, Pᴳ[i,t] - Pᴰ[i,t] ==  sum( Pˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )   
        @constraint(model, Qᴳ[i,t] - Qᴰ[i,t] ==  sum( Qˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )   

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



#--------------------------Voltage stability constraints--------------------------

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


#------------------------------
# Analyze results
#------------------------------

yˢᴳ²_opt=JuMP.value.(yˢᴳ²)       # store optimal values of decision variables
yˢᴳ³_opt=JuMP.value.(yˢᴳ³)
yˢᴳ⁴_opt=JuMP.value.(yˢᴳ⁴)
yˢᴳ⁵_opt=JuMP.value.(yˢᴳ⁵)
yˢᴳ²⁷_opt=JuMP.value.(yˢᴳ²⁷)
yˢᴳ³⁰_opt=JuMP.value.(yˢᴳ³⁰)

Pᴳ=JuMP.value.(Pᴳ)
sum(Pᴳ[:,1])
