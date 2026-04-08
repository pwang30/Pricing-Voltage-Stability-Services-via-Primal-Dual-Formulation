function MPEC_VS_LL_PD(K_c_gc, K_c_gv, K_c_m, sh_1, sh_2, bus_num, T, Pᴰ, Qᴰ, bus_connection, S_lines, 
                       Pˢᴳₘₐₓ, Pˢᴳₘᵢₙ, Qˢᴳₘₐₓ, Qˢᴳₘᵢₙ, 
                       Pⱽˢᴳₘₐₓ, Pⱽˢᴳₘᵢₙ, Qⱽˢᴳₘₐₓ, Qⱽˢᴳₘᵢₙ,
                       Pᴵᴮᴳₘₐₓ, Pᴵᴮᴳₘᵢₙ, Qᴵᴮᴳₘₐₓ, Qᴵᴮᴳₘᵢₙ,
                       Sᵐᵃˣ_gc, Sᵐᵃˣ_gv, Sᵐᵃˣ_c,
                       α_VSG, α_IBG, ratio,
                       Kˢᵗ, Kˢʰ, Oᵐ, Oⁿˡ,
                       bid_energy_SG_2, bid_energy_SG_3, bid_energy_SG_4, bid_energy_SG_5, bid_energy_SG_27, bid_energy_SG_30,
                       bid_SS_SG_2_bus_23, bid_SS_SG_3_bus_23, bid_SS_SG_4_bus_23, bid_SS_SG_5_bus_23, bid_SS_SG_27_bus_23, bid_SS_SG_30_bus_23,
                       bid_SS_SG_2_bus_24, bid_SS_SG_3_bus_24, bid_SS_SG_4_bus_24, bid_SS_SG_5_bus_24, bid_SS_SG_27_bus_24, bid_SS_SG_30_bus_24,
                       bid_SS_VSG_bus_23, bid_SS_VSG_bus_24,
                       bid_AP_IBG_23_bus_23, bid_AP_IBG_23_bus_24, bid_AP_IBG_24_bus_24, bid_AP_IBG_24_bus_23,
                       c_SS_bus_23, c_SS_bus_24, c_AP_bus_23, c_AP_bus_24, c_RP_bus_23, c_RP_bus_24 )


#------------------------------------------------------------------------------------------------------
# Define Primal in LL model_dual, where the bids of energy and ancillary service could be strategic
#------------------------------------------------------------------------------------------------------

model_primal= Model()                     

@variable(model_primal, yˢᴳ²[1:T],Bin)  # ,Bin                       # status of SGs, buses:2,3,4,5,27,30.    
@variable(model_primal, yˢᴳ³[1:T],Bin)            
@variable(model_primal, yˢᴳ⁴[1:T],Bin)           
@variable(model_primal, yˢᴳ⁵[1:T],Bin)           
@variable(model_primal, yˢᴳ²⁷[1:T],Bin)            
@variable(model_primal, yˢᴳ³⁰[1:T],Bin)  

@variable(model_primal, Cᵁ²[1:T]>=0)                           # startup costs and shutdown costs for SGs          
@variable(model_primal, Cᴰ²[1:T]>=0)                                  
@variable(model_primal, Cᵁ³[1:T]>=0)               
@variable(model_primal, Cᴰ³[1:T]>=0)                
@variable(model_primal, Cᵁ⁴[1:T]>=0)                 
@variable(model_primal, Cᴰ⁴[1:T]>=0)                 
@variable(model_primal, Cᵁ⁵[1:T]>=0)                 
@variable(model_primal, Cᴰ⁵[1:T]>=0)                 
@variable(model_primal, Cᵁ²⁷[1:T]>=0)                 
@variable(model_primal, Cᴰ²⁷[1:T]>=0)                
@variable(model_primal, Cᵁ³⁰[1:T]>=0)                 
@variable(model_primal, Cᴰ³⁰[1:T]>=0) 

@variable(model_primal, Pᴳ[1:30,1:T]) 
@variable(model_primal, Qᴳ[1:30,1:T])  

for t in 2:T                                                # startup costs and shutdown costs occured at hour t
    @constraint(model_primal, Cᵁ²[t]>=(yˢᴳ²[t]-yˢᴳ²[t-1])*Kˢᵗ[1])        
    @constraint(model_primal, Cᴰ²[t]>=(yˢᴳ²[t-1]-yˢᴳ²[t])*Kˢʰ[1])  
    @constraint(model_primal, Cᵁ³[t]>=(yˢᴳ³[t]-yˢᴳ³[t-1])*Kˢᵗ[2]) 
    @constraint(model_primal, Cᴰ³[t]>=(yˢᴳ³[t-1]-yˢᴳ³[t])*Kˢʰ[2])
    @constraint(model_primal, Cᵁ⁴[t]>=(yˢᴳ⁴[t]-yˢᴳ⁴[t-1])*Kˢᵗ[3])
    @constraint(model_primal, Cᴰ⁴[t]>=(yˢᴳ⁴[t-1]-yˢᴳ⁴[t])*Kˢʰ[3])
    @constraint(model_primal, Cᵁ⁵[t]>=(yˢᴳ⁵[t]-yˢᴳ⁵[t-1])*Kˢᵗ[4])
    @constraint(model_primal, Cᴰ⁵[t]>=(yˢᴳ⁵[t-1]-yˢᴳ⁵[t])*Kˢʰ[4])
    @constraint(model_primal, Cᵁ²⁷[t]>=(yˢᴳ²⁷[t]-yˢᴳ²⁷[t-1])*Kˢᵗ[5])
    @constraint(model_primal, Cᴰ²⁷[t]>=(yˢᴳ²⁷[t-1]-yˢᴳ²⁷[t])*Kˢʰ[5])
    @constraint(model_primal, Cᵁ³⁰[t]>=(yˢᴳ³⁰[t]-yˢᴳ³⁰[t-1])*Kˢᵗ[6])
    @constraint(model_primal, Cᴰ³⁰[t]>=(yˢᴳ³⁰[t-1]-yˢᴳ³⁰[t])*Kˢʰ[6])
end   



#---------------------------------------------------------------------------------------------------------------------------
#  AC power flow constraints
#---------------------------------------------------------------------------------------------------------------------------

for t in 1:T                                            # bounds for the active and reactive output of each bus
    @constraint(model_primal, Pᴳ[2,t] <= yˢᴳ²[t]*Pˢᴳₘₐₓ[1])    # bounds for active power of SGs: buses 2 3 4 5 27 30           
    @constraint(model_primal, yˢᴳ²[t]*Pˢᴳₘᵢₙ[1] <= Pᴳ[2,t])       
    @constraint(model_primal, Pᴳ[3,t] <= yˢᴳ³[t]*Pˢᴳₘₐₓ[2])       
    @constraint(model_primal, yˢᴳ³[t]*Pˢᴳₘᵢₙ[2] <= Pᴳ[3,t])         
    @constraint(model_primal, Pᴳ[4,t] <= yˢᴳ⁴[t]*Pˢᴳₘₐₓ[3])       
    @constraint(model_primal, yˢᴳ⁴[t]*Pˢᴳₘᵢₙ[3] <= Pᴳ[4,t])    
    @constraint(model_primal, Pᴳ[5,t] <= yˢᴳ⁵[t]*Pˢᴳₘₐₓ[4])       
    @constraint(model_primal, yˢᴳ⁵[t]*Pˢᴳₘᵢₙ[4] <= Pᴳ[5,t])
    @constraint(model_primal, Pᴳ[27,t] <= yˢᴳ²⁷[t]*Pˢᴳₘₐₓ[5])       
    @constraint(model_primal, yˢᴳ²⁷[t]*Pˢᴳₘᵢₙ[5] <=Pᴳ[27,t])
    @constraint(model_primal, Pᴳ[30,t] <= yˢᴳ³⁰[t]*Pˢᴳₘₐₓ[6])       
    @constraint(model_primal, yˢᴳ³⁰[t]*Pˢᴳₘᵢₙ[6] <= Pᴳ[30,t])

    @constraint(model_primal, Pᴳ[1,t] <= α_VSG*Pⱽˢᴳₘₐₓ[1])            # bounds for active power of VSGs and IBGs: buses 1 23 24           
    @constraint(model_primal, Pⱽˢᴳₘᵢₙ[1] <= Pᴳ[1,t])       
    @constraint(model_primal, Pᴳ[23,t] <= α_IBG*Pᴵᴮᴳₘₐₓ[1])       
    @constraint(model_primal, Pᴵᴮᴳₘᵢₙ[1] <= Pᴳ[23,t])         
    @constraint(model_primal, Pᴳ[24,t] <= α_IBG*Pᴵᴮᴳₘₐₓ[2])       
    @constraint(model_primal, Pᴵᴮᴳₘᵢₙ[2] <= Pᴳ[24,t])   

    @constraint(model_primal, Qᴳ[2,t] <= yˢᴳ²[t]*Qˢᴳₘₐₓ[1])    # bounds for REactive power of SGs: buses 2 3 4 5 27 30           
    @constraint(model_primal, yˢᴳ²[t]*Qˢᴳₘᵢₙ[1] <= Qᴳ[2,t])       
    @constraint(model_primal, Qᴳ[3,t] <= yˢᴳ³[t]*Qˢᴳₘₐₓ[2])       
    @constraint(model_primal, yˢᴳ³[t]*Qˢᴳₘᵢₙ[2] <= Qᴳ[3,t])         
    @constraint(model_primal, Qᴳ[4,t] <= yˢᴳ⁴[t]*Qˢᴳₘₐₓ[3])       
    @constraint(model_primal, yˢᴳ⁴[t]*Qˢᴳₘᵢₙ[3] <= Qᴳ[4,t])    
    @constraint(model_primal, Qᴳ[5,t] <= yˢᴳ⁵[t]*Qˢᴳₘₐₓ[4])       
    @constraint(model_primal, yˢᴳ⁵[t]*Qˢᴳₘᵢₙ[4] <= Qᴳ[5,t])
    @constraint(model_primal, Qᴳ[27,t] <= yˢᴳ²⁷[t]*Qˢᴳₘₐₓ[5])       
    @constraint(model_primal, yˢᴳ²⁷[t]*Qˢᴳₘᵢₙ[5] <=Qᴳ[27,t])
    @constraint(model_primal, Qᴳ[30,t] <= yˢᴳ³⁰[t]*Qˢᴳₘₐₓ[6])       
    @constraint(model_primal, yˢᴳ³⁰[t]*Qˢᴳₘᵢₙ[6] <= Qᴳ[30,t])           

    @constraint(model_primal, Qᴳ[1,t] <= Qⱽˢᴳₘₐₓ[1])       
    @constraint(model_primal, Qⱽˢᴳₘᵢₙ[1] <= Qᴳ[1,t])

    @constraint(model_primal, Qᴳ[23,t] <= ratio*Qᴵᴮᴳₘₐₓ[1])       
    @constraint(model_primal, ratio*Qᴵᴮᴳₘᵢₙ[1] <= Qᴳ[23,t])
    @constraint(model_primal, Qᴳ[24,t] <= ratio*Qᴵᴮᴳₘₐₓ[2])       
    @constraint(model_primal, ratio*Qᴵᴮᴳₘᵢₙ[2] <= Qᴳ[24,t])

    @constraint(model_primal, [ 1, Pᴳ[2,t] /(Sᵐᵃˣ_gc[1]), Qᴳ[2,t] /(Sᵐᵃˣ_gc[1]) ] in SecondOrderCone()  )
    @constraint(model_primal, [ 1, Pᴳ[3,t] /(Sᵐᵃˣ_gc[2]), Qᴳ[3,t] /(Sᵐᵃˣ_gc[2]) ] in SecondOrderCone()  )
    @constraint(model_primal, [ 1, Pᴳ[4,t] /(Sᵐᵃˣ_gc[3]), Qᴳ[4,t] /(Sᵐᵃˣ_gc[3]) ] in SecondOrderCone()  )
    @constraint(model_primal, [ 1, Pᴳ[5,t] /(Sᵐᵃˣ_gc[4]), Qᴳ[5,t] /(Sᵐᵃˣ_gc[4]) ] in SecondOrderCone()  )
    @constraint(model_primal, [ 1, Pᴳ[27,t] /(Sᵐᵃˣ_gc[5]), Qᴳ[27,t] /(Sᵐᵃˣ_gc[5]) ] in SecondOrderCone()  )
    @constraint(model_primal, [ 1, Pᴳ[30,t] /(Sᵐᵃˣ_gc[6]), Qᴳ[30,t] /(Sᵐᵃˣ_gc[6]) ] in SecondOrderCone()  )

    @constraint(model_primal, [ 1, Pᴳ[1,t] /(Sᵐᵃˣ_gv[1]), Qᴳ[1,t] /(Sᵐᵃˣ_gv[1]) ] in SecondOrderCone()  )

    @constraint(model_primal, [ 1, Pᴳ[23,t] /(Sᵐᵃˣ_c[1]), Qᴳ[23,t] /(Sᵐᵃˣ_c[1]) ] in SecondOrderCone()  )
    @constraint(model_primal, [ 1, Pᴳ[24,t] /(Sᵐᵃˣ_c[2]), Qᴳ[24,t] /(Sᵐᵃˣ_c[2]) ] in SecondOrderCone()  )

    gens_SGs = [2, 3, 4, 5, 27, 30]    # buses with SGs
    gens_VSGs = [1]                    # buses with VSGs
    gens_IBGs = [23, 24]               # buses with IBGs gens_all = union(gens_SGs, gens_VSGs, gens_IBGs)

    for i in 1:bus_num                 # no power is injected into the non-generation bus 
        if !(i in gens_SGs) && !(i in gens_VSGs) && !(i in gens_IBGs)
            @constraint(model_primal, Pᴳ[i, t] == 0)
            @constraint(model_primal, Qᴳ[i, t] == 0)
        end
    end
end

@variable(model_primal, Pˡⁱⁿᵉ[1:30,1:30,1:T])
@variable(model_primal, Qˡⁱⁿᵉ[1:30,1:30,1:T])

@variable(model_primal, Pˡⁱⁿᵉ_abs[1:30,1:30,1:T] >=0)
@variable(model_primal, Qˡⁱⁿᵉ_abs[1:30,1:30,1:T] >=0)

for i in 1:bus_num
    for t in 1:T

        @constraint(model_primal, Pᴳ[i,t] - Pᴰ[i,t] ==  sum( Pˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # p_i^g-p_i^d == G_{i,i}*c_{i,i}+...   eq. 50(a)
        @constraint(model_primal, Qᴳ[i,t] - Qᴰ[i,t] ==  sum( Qˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # q_i^g-q_i^d == -B_{i,i}*c_{i,i}+...  eq. 50(b)

        for j in bus_connection[i]
            @constraint(model_primal, Pˡⁱⁿᵉ[i,j,t] + Pˡⁱⁿᵉ[j,i,t] == 0)
            @constraint(model_primal, Qˡⁱⁿᵉ[i,j,t] + Qˡⁱⁿᵉ[j,i,t] == 0)

            @constraint(model_primal, Pˡⁱⁿᵉ[i,j,t] <= Pˡⁱⁿᵉ_abs[i,j,t])
            @constraint(model_primal, -Pˡⁱⁿᵉ[i,j,t] <= Pˡⁱⁿᵉ_abs[i,j,t])
            @constraint(model_primal, Qˡⁱⁿᵉ[i,j,t] <= Qˡⁱⁿᵉ_abs[i,j,t])
            @constraint(model_primal, -Qˡⁱⁿᵉ[i,j,t] <= Qˡⁱⁿᵉ_abs[i,j,t])
            @constraint(model_primal, Pˡⁱⁿᵉ_abs[i,j,t] <= S_lines[i,j])
            @constraint(model_primal, Qˡⁱⁿᵉ_abs[i,j,t] <= S_lines[i,j])
            @constraint(model_primal, Pˡⁱⁿᵉ_abs[i,j,t] + Qˡⁱⁿᵉ_abs[i,j,t] <= S_lines[i,j]*sqrt(2))

        end    
    end
end 




#--------------------------------
#  Voltage stability constraints
#--------------------------------

@variable(model_primal, P_23_eq[1:T])              # equivalent P and Q of IBGs, buses:23,24    and Γ_c  
@variable(model_primal, P_24_eq[1:T]) 
@variable(model_primal, Q_23_eq[1:T])                 
@variable(model_primal, Q_24_eq[1:T])
@variable(model_primal, Γ_23[1:T])
@variable(model_primal, Γ_24[1:T])

@variable(model_primal, ηₘ_1[1:15,1:T],Bin)        # product of each pair of SG 
@variable(model_primal, ηₘ_2[1:6,1:T]>=0)             # product of each pair of SG and VSG 

@variable(model_primal, s_g2_P23[1:T])
@variable(model_primal, s_g3_P23[1:T])
@variable(model_primal, s_g4_P23[1:T])
@variable(model_primal, s_g5_P23[1:T])
@variable(model_primal, s_g27_P23[1:T])
@variable(model_primal, s_g30_P23[1:T])
@variable(model_primal, s_g2_P24[1:T])
@variable(model_primal, s_g3_P24[1:T])
@variable(model_primal, s_g4_P24[1:T])
@variable(model_primal, s_g5_P24[1:T])
@variable(model_primal, s_g27_P24[1:T])
@variable(model_primal, s_g30_P24[1:T])

@variable(model_primal, s_g2_Q23[1:T])
@variable(model_primal, s_g3_Q23[1:T])
@variable(model_primal, s_g4_Q23[1:T])
@variable(model_primal, s_g5_Q23[1:T])
@variable(model_primal, s_g27_Q23[1:T])
@variable(model_primal, s_g30_Q23[1:T])
@variable(model_primal, s_g2_Q24[1:T])
@variable(model_primal, s_g3_Q24[1:T])
@variable(model_primal, s_g4_Q24[1:T])
@variable(model_primal, s_g5_Q24[1:T])
@variable(model_primal, s_g27_Q24[1:T])
@variable(model_primal, s_g30_Q24[1:T])

@variable(model_primal, x_ηₘ_1_P23[1:15,1:T])
@variable(model_primal, x_ηₘ_1_P24[1:15,1:T])
@variable(model_primal, x_ηₘ_1_Q23[1:15,1:T])
@variable(model_primal, x_ηₘ_1_Q24[1:15,1:T])
@variable(model_primal, x_ηₘ_2_P23[1:6,1:T])
@variable(model_primal, x_ηₘ_2_P24[1:6,1:T])
@variable(model_primal, x_ηₘ_2_Q23[1:6,1:T])
@variable(model_primal, x_ηₘ_2_Q24[1:6,1:T])
for t in 1:T
    # 15+6
    @constraint(model_primal, ηₘ_1[1,t]>=yˢᴳ²[t]+yˢᴳ³[t]-1)  
    @constraint(model_primal, ηₘ_1[1,t]<=yˢᴳ²[t])
    @constraint(model_primal, ηₘ_1[1,t]<=yˢᴳ³[t])
    @constraint(model_primal, ηₘ_1[2,t]>=yˢᴳ²[t]+yˢᴳ⁴[t]-1)  
    @constraint(model_primal, ηₘ_1[2,t]<=yˢᴳ²[t])
    @constraint(model_primal, ηₘ_1[2,t]<=yˢᴳ⁴[t])
    @constraint(model_primal, ηₘ_1[3,t]>=yˢᴳ²[t]+yˢᴳ⁵[t]-1)  
    @constraint(model_primal, ηₘ_1[3,t]<=yˢᴳ²[t])
    @constraint(model_primal, ηₘ_1[3,t]<=yˢᴳ⁵[t])
    @constraint(model_primal, ηₘ_1[4,t]>=yˢᴳ²[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model_primal, ηₘ_1[4,t]<=yˢᴳ²[t])
    @constraint(model_primal, ηₘ_1[4,t]<=yˢᴳ²⁷[t])
    @constraint(model_primal, ηₘ_1[5,t]>=yˢᴳ²[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model_primal, ηₘ_1[5,t]<=yˢᴳ²[t])
    @constraint(model_primal, ηₘ_1[5,t]<=yˢᴳ³⁰[t])
    @constraint(model_primal, ηₘ_1[6,t]>=yˢᴳ³[t]+yˢᴳ⁴[t]-1)  
    @constraint(model_primal, ηₘ_1[6,t]<=yˢᴳ³[t])
    @constraint(model_primal, ηₘ_1[6,t]<=yˢᴳ⁴[t])
    @constraint(model_primal, ηₘ_1[7,t]>=yˢᴳ³[t]+yˢᴳ⁵[t]-1)  
    @constraint(model_primal, ηₘ_1[7,t]<=yˢᴳ³[t])
    @constraint(model_primal, ηₘ_1[7,t]<=yˢᴳ⁵[t])
    @constraint(model_primal, ηₘ_1[8,t]>=yˢᴳ³[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model_primal, ηₘ_1[8,t]<=yˢᴳ³[t])
    @constraint(model_primal, ηₘ_1[8,t]<=yˢᴳ²⁷[t])
    @constraint(model_primal, ηₘ_1[9,t]>=yˢᴳ³[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model_primal, ηₘ_1[9,t]<=yˢᴳ³[t])
    @constraint(model_primal, ηₘ_1[9,t]<=yˢᴳ³⁰[t])

    @constraint(model_primal, ηₘ_1[10,t]>=yˢᴳ⁴[t]+yˢᴳ⁵[t]-1)  
    @constraint(model_primal, ηₘ_1[10,t]<=yˢᴳ⁴[t])
    @constraint(model_primal, ηₘ_1[10,t]<=yˢᴳ⁵[t])
    @constraint(model_primal, ηₘ_1[11,t]>=yˢᴳ⁴[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model_primal, ηₘ_1[11,t]<=yˢᴳ⁴[t])
    @constraint(model_primal, ηₘ_1[11,t]<=yˢᴳ²⁷[t])
    @constraint(model_primal, ηₘ_1[12,t]>=yˢᴳ⁴[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model_primal, ηₘ_1[12,t]<=yˢᴳ⁴[t])
    @constraint(model_primal, ηₘ_1[12,t]<=yˢᴳ³⁰[t])

    @constraint(model_primal, ηₘ_1[13,t]>=yˢᴳ⁵[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model_primal, ηₘ_1[13,t]<=yˢᴳ⁵[t])
    @constraint(model_primal, ηₘ_1[13,t]<=yˢᴳ²⁷[t])
    @constraint(model_primal, ηₘ_1[14,t]>=yˢᴳ⁵[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model_primal, ηₘ_1[14,t]<=yˢᴳ⁵[t])
    @constraint(model_primal, ηₘ_1[14,t]<=yˢᴳ³⁰[t])

    @constraint(model_primal, ηₘ_1[15,t]>=yˢᴳ²⁷[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model_primal, ηₘ_1[15,t]<=yˢᴳ²⁷[t])
    @constraint(model_primal, ηₘ_1[15,t]<=yˢᴳ³⁰[t])
    @constraint(model_primal, ηₘ_2[1,t]==yˢᴳ²[t]*α_VSG)
    @constraint(model_primal, ηₘ_2[2,t]==yˢᴳ³[t]*α_VSG)
    @constraint(model_primal, ηₘ_2[3,t]==yˢᴳ⁴[t]*α_VSG)
    @constraint(model_primal, ηₘ_2[4,t]==yˢᴳ⁵[t]*α_VSG)
    @constraint(model_primal, ηₘ_2[5,t]==yˢᴳ²⁷[t]*α_VSG)
    @constraint(model_primal, ηₘ_2[6,t]==yˢᴳ³⁰[t]*α_VSG)

    @constraint(model_primal, s_g2_P23[t] <=yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g2_P23[t] <= Pᴳ[23,t] +yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g2_P23[t] >=Pᴳ[23,t]+yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g2_P23[t] >= yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g3_P23[t] <=yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g3_P23[t] <= Pᴳ[23,t] +yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g3_P23[t] >=Pᴳ[23,t]+yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g3_P23[t] >= yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g4_P23[t] <=yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g4_P23[t] <= Pᴳ[23,t] +yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g4_P23[t] >=Pᴳ[23,t]+yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g4_P23[t] >= yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g5_P23[t] <=yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g5_P23[t] <= Pᴳ[23,t] +yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g5_P23[t] >=Pᴳ[23,t]+yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])   
    @constraint(model_primal, s_g5_P23[t] >= yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g27_P23[t] <=yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g27_P23[t] <= Pᴳ[23,t] +yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g27_P23[t] >=Pᴳ[23,t]+yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])   
    @constraint(model_primal, s_g27_P23[t] >= yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g30_P23[t] <=yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g30_P23[t] <= Pᴳ[23,t] +yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g30_P23[t] >=Pᴳ[23,t]+yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])   
    @constraint(model_primal, s_g30_P23[t] >= yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g2_P24[t] <=yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g2_P24[t] <= Pᴳ[24,t] +yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g2_P24[t] >=Pᴳ[24,t]+yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g2_P24[t] >= yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, s_g3_P24[t] <=yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g3_P24[t] <= Pᴳ[24,t] +yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g3_P24[t] >=Pᴳ[24,t]+yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g3_P24[t] >= yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, s_g4_P24[t] <=yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g4_P24[t] <= Pᴳ[24,t] +yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g4_P24[t] >=Pᴳ[24,t]+yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g4_P24[t] >= yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, s_g5_P24[t] <=yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g5_P24[t] <= Pᴳ[24,t] +yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g5_P24[t] >=Pᴳ[24,t]+yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])   
    @constraint(model_primal, s_g5_P24[t] >= yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, s_g27_P24[t] <=yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g27_P24[t] <= Pᴳ[24,t] +yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g27_P24[t] >=Pᴳ[24,t]+yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])   
    @constraint(model_primal, s_g27_P24[t] >= yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, s_g30_P24[t] <=yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g30_P24[t] <= Pᴳ[24,t] +yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g30_P24[t] >=Pᴳ[24,t]+yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])   
    @constraint(model_primal, s_g30_P24[t] >= yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[2] )

    
    
    @constraint(model_primal, s_g2_Q23[t] <= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g2_Q23[t] <= Qᴳ[23,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g2_Q23[t] >= Qᴳ[23,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g2_Q23[t] >= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g3_Q23[t] <= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g3_Q23[t] <= Qᴳ[23,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g3_Q23[t] >= Qᴳ[23,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g3_Q23[t] >= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g4_Q23[t] <= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g4_Q23[t] <= Qᴳ[23,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g4_Q23[t] >= Qᴳ[23,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g4_Q23[t] >= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g5_Q23[t] <= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g5_Q23[t] <= Qᴳ[23,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g5_Q23[t] >= Qᴳ[23,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])   
    @constraint(model_primal, s_g5_Q23[t] >= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g27_Q23[t] <= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g27_Q23[t] <= Qᴳ[23,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g27_Q23[t] >= Qᴳ[23,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])   
    @constraint(model_primal, s_g27_Q23[t] >= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g30_Q23[t] <= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_primal, s_g30_Q23[t] <= Qᴳ[23,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, s_g30_Q23[t] >= Qᴳ[23,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])   
    @constraint(model_primal, s_g30_Q23[t] >= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, s_g2_Q24[t] <= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g2_Q24[t] <= Qᴳ[24,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g2_Q24[t] >= Qᴳ[24,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g2_Q24[t] >= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )   

    @constraint(model_primal, s_g3_Q24[t] <= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g3_Q24[t] <= Qᴳ[24,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g3_Q24[t] >= Qᴳ[24,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g3_Q24[t] >= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, s_g4_Q24[t] <= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g4_Q24[t] <= Qᴳ[24,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g4_Q24[t] >= Qᴳ[24,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g4_Q24[t] >= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, s_g5_Q24[t] <= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g5_Q24[t] <= Qᴳ[24,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g5_Q24[t] >= Qᴳ[24,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])   
    @constraint(model_primal, s_g5_Q24[t] >= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, s_g27_Q24[t] <= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g27_Q24[t] <= Qᴳ[24,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] ) 
    @constraint(model_primal, s_g27_Q24[t] >= Qᴳ[24,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])   
    @constraint(model_primal, s_g27_Q24[t] >= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, s_g30_Q24[t] <= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_primal, s_g30_Q24[t] <= Qᴳ[24,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, s_g30_Q24[t] >= Qᴳ[24,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])   
    @constraint(model_primal, s_g30_Q24[t] >= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, x_ηₘ_1_P23[:,t] <= ηₘ_1[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model_primal, x_ηₘ_1_P23[:,t] <= Pᴳ[23,t] .+ ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[1] .- Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, x_ηₘ_1_P23[:,t] >= Pᴳ[23,t] .+ ηₘ_1[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[1] .- α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model_primal, x_ηₘ_1_P23[:,t] >= ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, x_ηₘ_2_P23[:,t] <= ηₘ_2[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model_primal, x_ηₘ_2_P23[:,t] <= Pᴳ[23,t] .+ ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[1] .- Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, x_ηₘ_2_P23[:,t] >= Pᴳ[23,t] .+ ηₘ_2[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[1] .- α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model_primal, x_ηₘ_2_P23[:,t] >= ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, x_ηₘ_1_P24[:,t] <= ηₘ_1[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model_primal, x_ηₘ_1_P24[:,t] <= Pᴳ[24,t] .+ ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[2] .- Pᴵᴮᴳₘᵢₙ[2] )    
    @constraint(model_primal, x_ηₘ_1_P24[:,t] >= Pᴳ[24,t] .+ ηₘ_1[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[2] .- α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model_primal, x_ηₘ_1_P24[:,t] >= ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, x_ηₘ_2_P24[:,t] <= ηₘ_2[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model_primal, x_ηₘ_2_P24[:,t] <= Pᴳ[24,t] .+ ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[2] .- Pᴵᴮᴳₘᵢₙ[2] )    
    @constraint(model_primal, x_ηₘ_2_P24[:,t] >= Pᴳ[24,t] .+ ηₘ_2[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[2] .- α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model_primal, x_ηₘ_2_P24[:,t] >= ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, x_ηₘ_1_Q23[:,t] <= ηₘ_1[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[1] )  
    @constraint(model_primal, x_ηₘ_1_Q23[:,t] <= Qᴳ[23,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] .- ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, x_ηₘ_1_Q23[:,t] >= Qᴳ[23,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘₐₓ[1] .- ratio*Qᴵᴮᴳₘₐₓ[1] )
    @constraint(model_primal, x_ηₘ_1_Q23[:,t] >= ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] )
    
    @constraint(model_primal, x_ηₘ_2_Q23[:,t] <= ηₘ_2[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[1] )  
    @constraint(model_primal, x_ηₘ_2_Q23[:,t] <= Qᴳ[23,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] .- ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_primal, x_ηₘ_2_Q23[:,t] >= Qᴳ[23,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘₐₓ[1] .- ratio*Qᴵᴮᴳₘₐₓ[1] )
    @constraint(model_primal, x_ηₘ_2_Q23[:,t] >= ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_primal, x_ηₘ_1_Q24[:,t] <= ηₘ_1[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[2] )  
    @constraint(model_primal, x_ηₘ_1_Q24[:,t] <= Qᴳ[24,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] .- ratio*Qᴵᴮᴳₘᵢₙ[2] )    
    @constraint(model_primal, x_ηₘ_1_Q24[:,t] >= Qᴳ[24,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘₐₓ[2] .- ratio*Qᴵᴮᴳₘₐₓ[2] )
    @constraint(model_primal, x_ηₘ_1_Q24[:,t] >= ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, x_ηₘ_2_Q24[:,t] <= ηₘ_2[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[2] )  
    @constraint(model_primal, x_ηₘ_2_Q24[:,t] <= Qᴳ[24,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] .- ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_primal, x_ηₘ_2_Q24[:,t] >= Qᴳ[24,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘₐₓ[2] .- ratio*Qᴵᴮᴳₘₐₓ[2] )
    @constraint(model_primal, x_ηₘ_2_Q24[:,t] >= ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_primal, P_23_eq[t]== Pᴳ[23,t]+  ( sh_1- ( K_c_gc[1,1] *s_g2_P24[t] +K_c_gc[1,2] *s_g3_P24[t] +K_c_gc[1,3] *s_g4_P24[t] 
    +K_c_gc[1,4] *s_g5_P24[t] +K_c_gc[1,5] *s_g27_P24[t] +K_c_gc[1,6] *s_g30_P24[t]    +α_VSG*K_c_gv[1,1]*Pᴳ[24,t] 
    + sum(K_c_m[1,1:15] .*x_ηₘ_1_P24[:,t]) +  sum(K_c_m[1,16:end] .*x_ηₘ_2_P24[:,t] ) ) )    )

    @constraint(model_primal, P_24_eq[t]== Pᴳ[24,t]+  ( sh_2- ( K_c_gc[2,1] *s_g2_P23[t] +K_c_gc[2,2] *s_g3_P23[t] +K_c_gc[2,3] *s_g4_P23[t] 
    +K_c_gc[2,4] *s_g5_P23[t] +K_c_gc[2,5] *s_g27_P23[t] +K_c_gc[2,6] *s_g30_P23[t]    +α_VSG*K_c_gv[2,1]*Pᴳ[23,t]
    + sum(K_c_m[2,1:15] .*x_ηₘ_1_P23[:,t]) +  sum(K_c_m[2,16:end] .*x_ηₘ_2_P23[:,t] ) ) )  )

    @constraint(model_primal, Q_23_eq[t]== Qᴳ[23,t]+  ( sh_1- ( K_c_gc[1,1] *s_g2_Q24[t] +K_c_gc[1,2] *s_g3_Q24[t] +K_c_gc[1,3] *s_g4_Q24[t] 
    +K_c_gc[1,4] *s_g5_Q24[t] +K_c_gc[1,5] *s_g27_Q24[t] +K_c_gc[1,6] *s_g30_Q24[t]    +α_VSG*K_c_gv[1,1]*Qᴳ[24,t] 
    + sum(K_c_m[1,1:15] .*x_ηₘ_1_Q24[:,t]) +  sum(K_c_m[1,16:end] .*x_ηₘ_2_Q24[:,t] ) ) )    )

    @constraint(model_primal, Q_24_eq[t]== Qᴳ[24,t]+  ( sh_2- ( K_c_gc[2,1] *s_g2_Q23[t] +K_c_gc[2,2] *s_g3_Q23[t] +K_c_gc[2,3] *s_g4_Q23[t] 
    +K_c_gc[2,4] *s_g5_Q23[t] +K_c_gc[2,5] *s_g27_Q23[t] +K_c_gc[2,6] *s_g30_Q23[t]    +α_VSG*K_c_gv[2,1]*Pᴳ[23,t]
    + sum(K_c_m[2,1:15] .*x_ηₘ_1_Q23[:,t]) +  sum(K_c_m[2,16:end] .*x_ηₘ_2_Q23[:,t] ) ) )    )

    @constraint(model_primal, Γ_23[t]== 1/2 * ( K_c_gc[3,1] *yˢᴳ²[t] +K_c_gc[3,2] *yˢᴳ³[t] +K_c_gc[3,3] *yˢᴳ⁴[t] 
    +K_c_gc[3,4] *yˢᴳ⁵[t] +K_c_gc[3,5] *yˢᴳ²⁷[t] +K_c_gc[3,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[3,16:end] .*ηₘ_2[:,t]   ) )

    @constraint(model_primal, Γ_24[t]== 1/2 * ( K_c_gc[4,1] *yˢᴳ²[t] +K_c_gc[4,2] *yˢᴳ³[t] +K_c_gc[4,3] *yˢᴳ⁴[t] 
    +K_c_gc[4,4] *yˢᴳ⁵[t] +K_c_gc[4,5] *yˢᴳ²⁷[t] +K_c_gc[4,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[4,1]
    + sum(K_c_m[4,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[4,16:end] .*ηₘ_2[:,t]   ) )

    @constraint(model_primal, [Q_23_eq[t] + Γ_23[t]*100, P_23_eq[t], Q_23_eq[t]] in SecondOrderCone()  )      
    @constraint(model_primal, [Q_24_eq[t] + Γ_24[t]*100, P_24_eq[t], Q_24_eq[t]] in SecondOrderCone()  )

end


@variable(model_primal, ΔΓ_SG_2_bus_23[1:T] )
@variable(model_primal, ΔΓ_SG_2_bus_24[1:T] )
@variable(model_primal, ΔΓ_SG_3_bus_23[1:T] )
@variable(model_primal, ΔΓ_SG_3_bus_24[1:T] )
@variable(model_primal, ΔΓ_SG_4_bus_23[1:T] )
@variable(model_primal, ΔΓ_SG_4_bus_24[1:T] )
@variable(model_primal, ΔΓ_SG_5_bus_23[1:T] )
@variable(model_primal, ΔΓ_SG_5_bus_24[1:T] )
@variable(model_primal, ΔΓ_SG_27_bus_23[1:T] )
@variable(model_primal, ΔΓ_SG_27_bus_24[1:T] )
@variable(model_primal, ΔΓ_SG_30_bus_23[1:T] )
@variable(model_primal, ΔΓ_SG_30_bus_24[1:T] )
@variable(model_primal, ΔΓ_VSG_bus_23[1:T] )
@variable(model_primal, ΔΓ_VSG_bus_24[1:T] )

for t in 1:T

    @constraint(model_primal, ΔΓ_SG_2_bus_23[t] == Γ_23[t] - ( 1/2 * ( K_c_gc[3,2] *yˢᴳ³[t] +K_c_gc[3,3] *yˢᴳ⁴[t] 
    +K_c_gc[3,4] *yˢᴳ⁵[t] +K_c_gc[3,5] *yˢᴳ²⁷[t] +K_c_gc[3,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,6:15] .*ηₘ_1[6:15 , t]) +  sum(K_c_m[3,17:end] .*ηₘ_2[2:end , t]   ) ) )

    @constraint(model_primal, ΔΓ_SG_2_bus_24[t] == Γ_24[t] - ( 1/2 * ( K_c_gc[4,2] *yˢᴳ³[t] +K_c_gc[4,3] *yˢᴳ⁴[t] 
    +K_c_gc[4,4] *yˢᴳ⁵[t] +K_c_gc[4,5] *yˢᴳ²⁷[t] +K_c_gc[4,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[4,1]
    + sum(K_c_m[4,6:15] .*ηₘ_1[6:15 , t]) +  sum(K_c_m[4,17:end] .*ηₘ_2[2:end , t]   ) ) )

    @constraint(model_primal, ΔΓ_SG_3_bus_23[t] == Γ_23[t] - ( 1/2 * ( K_c_gc[3,1] *yˢᴳ²[t] +K_c_gc[3,3] *yˢᴳ⁴[t] 
    +K_c_gc[3,4] *yˢᴳ⁵[t] +K_c_gc[3,5] *yˢᴳ²⁷[t] +K_c_gc[3,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,2:5] .*ηₘ_1[2:5 , t]) + sum(K_c_m[3,10:15] .*ηₘ_1[10:15 , t]) +  sum(K_c_m[3,16] .*ηₘ_2[1 , t]) + sum(K_c_m[3,18:end] .*ηₘ_2[3:end , t] ) ) )

    @constraint(model_primal, ΔΓ_SG_3_bus_24[t] == Γ_24[t] - ( 1/2 * ( K_c_gc[4,1] *yˢᴳ²[t] +K_c_gc[4,3] *yˢᴳ⁴[t] 
    +K_c_gc[4,4] *yˢᴳ⁵[t] +K_c_gc[4,5] *yˢᴳ²⁷[t] +K_c_gc[4,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[4,1]
    + sum(K_c_m[4,2:5] .*ηₘ_1[2:5 , t]) + sum(K_c_m[4,10:15] .*ηₘ_1[10:15 , t]) +  sum(K_c_m[4,16] .*ηₘ_2[1 , t]) + sum(K_c_m[4,18:end] .*ηₘ_2[3:end , t] ) ) )

    @constraint(model_primal, ΔΓ_SG_4_bus_23[t] == Γ_23[t] - ( 1/2 * ( K_c_gc[3,1] *yˢᴳ²[t] +K_c_gc[3,2] *yˢᴳ³[t] 
    +K_c_gc[3,4] *yˢᴳ⁵[t] +K_c_gc[3,5] *yˢᴳ²⁷[t] +K_c_gc[3,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,1] .*ηₘ_1[1 , t]) + sum(K_c_m[3,3:5] .*ηₘ_1[3:5 , t]) +  sum(K_c_m[3,7:9] .*ηₘ_1[7:9 , t]) + sum(K_c_m[3,13:15] .*ηₘ_1[13:15 , t]) + sum(K_c_m[3,16:17] .*ηₘ_2[1:2 , t]) + sum(K_c_m[3,19:end] .*ηₘ_2[4:end , t] ) ) )

    @constraint(model_primal, ΔΓ_SG_4_bus_24[t] == Γ_24[t] - ( 1/2 * ( K_c_gc[4,1] *yˢᴳ²[t] +K_c_gc[4,2] *yˢᴳ³[t] 
    +K_c_gc[4,4] *yˢᴳ⁵[t] +K_c_gc[4,5] *yˢᴳ²⁷[t] +K_c_gc[4,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[4,1]
    + sum(K_c_m[4,1] .*ηₘ_1[1 , t]) + sum(K_c_m[4,3:5] .*ηₘ_1[3:5 , t]) +  sum(K_c_m[4,7:9] .*ηₘ_1[7:9 , t]) + sum(K_c_m[4,13:15] .*ηₘ_1[13:15 , t]) + sum(K_c_m[4,16:17] .*ηₘ_2[1:2 , t]) + sum(K_c_m[4,19:end] .*ηₘ_2[4:end , t] ) ) )

    @constraint(model_primal, ΔΓ_SG_5_bus_23[t] == Γ_23[t] - ( 1/2 * ( K_c_gc[3,1] *yˢᴳ²[t] +K_c_gc[3,2] *yˢᴳ³[t] 
    +K_c_gc[3,3] *yˢᴳ⁴[t] +K_c_gc[3,5] *yˢᴳ²⁷[t] +K_c_gc[3,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,1:2] .*ηₘ_1[1:2 , t]) + sum(K_c_m[3,4:6] .*ηₘ_1[4:6 , t]) + sum(K_c_m[3,8:9] .*ηₘ_1[8:9 , t]) + sum(K_c_m[3,11:12] .*ηₘ_1[11:12 , t]) 
    + sum(K_c_m[3,15] .*ηₘ_1[15 , t]) +  sum(K_c_m[3,16:18] .*ηₘ_2[1:3 , t]   ) +  sum(K_c_m[3,20:end] .*ηₘ_2[5:end , t]   ) ) )

    @constraint(model_primal, ΔΓ_SG_5_bus_24[t] == Γ_24[t] - ( 1/2 * ( K_c_gc[4,1] *yˢᴳ²[t] +K_c_gc[4,2] *yˢᴳ³[t] 
    +K_c_gc[4,3] *yˢᴳ⁴[t] +K_c_gc[4,5] *yˢᴳ²⁷[t] +K_c_gc[4,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[4,1]
    + sum(K_c_m[4,1:2] .*ηₘ_1[1:2 , t]) + sum(K_c_m[4,4:6] .*ηₘ_1[4:6 , t]) + sum(K_c_m[4,8:9] .*ηₘ_1[8:9 , t]) + sum(K_c_m[4,11:12] .*ηₘ_1[11:12 , t]) 
    + sum(K_c_m[4,15] .*ηₘ_1[15 , t]) +  sum(K_c_m[4,16:18] .*ηₘ_2[1:3 , t]   ) +  sum(K_c_m[4,20:end] .*ηₘ_2[5:end , t]   ) ) )

    @constraint(model_primal, ΔΓ_SG_27_bus_23[t] == Γ_23[t] - ( 1/2 * ( K_c_gc[3,1] *yˢᴳ²[t] +K_c_gc[3,2] *yˢᴳ³[t] 
    +K_c_gc[3,3] *yˢᴳ⁴[t] +K_c_gc[3,4] *yˢᴳ⁵[t] +K_c_gc[3,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,1:3] .*ηₘ_1[1:3 , t]) + sum(K_c_m[3,5:7] .*ηₘ_1[5:7 , t]) + sum(K_c_m[3,9:10] .*ηₘ_1[9:10 , t]) + sum(K_c_m[3,12] .*ηₘ_1[12 , t]) 
    + sum(K_c_m[3,14] .*ηₘ_1[14 , t]) + sum(K_c_m[3,16:19] .*ηₘ_2[1:4 , t])  + sum(K_c_m[3,21] .*ηₘ_2[6 , t] )  ) )

    @constraint(model_primal, ΔΓ_SG_27_bus_24[t] == Γ_24[t] - ( 1/2 * ( K_c_gc[4,1] *yˢᴳ²[t] +K_c_gc[4,2] *yˢᴳ³[t] 
    +K_c_gc[4,3] *yˢᴳ⁴[t] +K_c_gc[4,4] *yˢᴳ⁵[t] +K_c_gc[4,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[4,1]
    + sum(K_c_m[4,1:3] .*ηₘ_1[1:3 , t]) + sum(K_c_m[4,5:7] .*ηₘ_1[5:7 , t]) + sum(K_c_m[4,9:10] .*ηₘ_1[9:10 , t]) + sum(K_c_m[4,12] .*ηₘ_1[12 , t]) 
    + sum(K_c_m[4,14] .*ηₘ_1[14 , t]) + sum(K_c_m[4,16:19] .*ηₘ_2[1:4 , t])  + sum(K_c_m[4,21] .*ηₘ_2[6 , t] )  ) )

    @constraint(model_primal, ΔΓ_SG_30_bus_23[t] == Γ_23[t] - ( 1/2 * ( K_c_gc[3,1] *yˢᴳ²[t] +K_c_gc[3,2] *yˢᴳ³[t] 
    +K_c_gc[3,3] *yˢᴳ⁴[t] +K_c_gc[3,4] *yˢᴳ⁵[t] +K_c_gc[3,5] *yˢᴳ²⁷[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,1:4] .*ηₘ_1[1:4 , t]) + sum(K_c_m[3,6:8] .*ηₘ_1[6:8 , t]) + sum(K_c_m[3,10:11] .*ηₘ_1[10:11 , t]) + sum(K_c_m[3,13] .*ηₘ_1[13 , t]) 
    + sum(K_c_m[3,15] .*ηₘ_1[15 , t]) +  sum(K_c_m[3,16:20] .*ηₘ_2[1:5 , t]   )   ) )

    @constraint(model_primal, ΔΓ_SG_30_bus_24[t] == Γ_24[t] - ( 1/2 * ( K_c_gc[4,1] *yˢᴳ²[t] +K_c_gc[4,2] *yˢᴳ³[t] 
    +K_c_gc[4,3] *yˢᴳ⁴[t] +K_c_gc[4,4] *yˢᴳ⁵[t] +K_c_gc[4,5] *yˢᴳ²⁷[t] )   +α_VSG*K_c_gv[4,1]
    + sum(K_c_m[4,1:4] .*ηₘ_1[1:4 , t]) + sum(K_c_m[4,6:8] .*ηₘ_1[6:8 , t]) + sum(K_c_m[4,10:11] .*ηₘ_1[10:11 , t]) + sum(K_c_m[4,13] .*ηₘ_1[13 , t]) 
    + sum(K_c_m[4,15] .*ηₘ_1[15 , t]) +  sum(K_c_m[4,16:20] .*ηₘ_2[1:5 , t]   )   ) )

    @constraint(model_primal, ΔΓ_VSG_bus_23[t] == Γ_23[t] - ( 1/2 * ( K_c_gc[3,1] *yˢᴳ²[t] +K_c_gc[3,2] *yˢᴳ³[t] +K_c_gc[3,3] *yˢᴳ⁴[t] 
    +K_c_gc[3,4] *yˢᴳ⁵[t] +K_c_gc[3,5] *yˢᴳ²⁷[t] +K_c_gc[3,6] *yˢᴳ³⁰[t] )    + sum(K_c_m[3,1:15] .*ηₘ_1[:,t])  ) )

    @constraint(model_primal, ΔΓ_VSG_bus_24[t] == Γ_24[t] - ( 1/2 * ( K_c_gc[4,1] *yˢᴳ²[t] +K_c_gc[4,2] *yˢᴳ³[t] +K_c_gc[4,3] *yˢᴳ⁴[t] 
    +K_c_gc[4,4] *yˢᴳ⁵[t] +K_c_gc[4,5] *yˢᴳ²⁷[t] +K_c_gc[4,6] *yˢᴳ³⁰[t] )    + sum(K_c_m[4,1:15] .*ηₘ_1[:,t])  ) )



end

#--------------------------------
#  Model solving
#--------------------------------
# obj of SGs for operation cost and ancillary service cost
cost_SG_2 = sum(Cᵁ²)+sum(Cᴰ²) + sum(Oⁿˡ[1].*(yˢᴳ²)) + bid_energy_SG_2 *sum(Oᵐ[1].*Pᴳ[2,:]) + 
            bid_SS_SG_2_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_2_bus_23 ) + bid_SS_SG_2_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_2_bus_24 ) 
cost_SG_3 = sum(Cᵁ³)+sum(Cᴰ³) + sum(Oⁿˡ[2].*(yˢᴳ³)) + bid_energy_SG_3 *sum(Oᵐ[2].*Pᴳ[3,:]) + 
            bid_SS_SG_3_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_3_bus_23 ) + bid_SS_SG_3_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_3_bus_24 )
cost_SG_4 = sum(Cᵁ⁴)+sum(Cᴰ⁴) + sum(Oⁿˡ[3].*(yˢᴳ⁴)) + bid_energy_SG_4 *sum(Oᵐ[3].*Pᴳ[4,:]) + 
            bid_SS_SG_4_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_4_bus_23 ) + bid_SS_SG_4_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_4_bus_24 )
cost_SG_5 = sum(Cᵁ⁵)+sum(Cᴰ⁵) + sum(Oⁿˡ[4].*(yˢᴳ⁵)) + bid_energy_SG_5 *sum(Oᵐ[4].*Pᴳ[5,:]) + 
            bid_SS_SG_5_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_5_bus_23 ) + bid_SS_SG_5_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_5_bus_24 )
cost_SG_27 = sum(Cᵁ²⁷)+sum(Cᴰ²⁷) + sum(Oⁿˡ[5].*(yˢᴳ²⁷)) + bid_energy_SG_27 *sum(Oᵐ[5].*Pᴳ[27,:]) + 
            bid_SS_SG_27_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_27_bus_23 ) + bid_SS_SG_27_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_27_bus_24 )
cost_SG_30 = sum(Cᵁ³⁰)+sum(Cᴰ³⁰) + sum(Oⁿˡ[6].*(yˢᴳ³⁰)) + bid_energy_SG_30 *sum(Oᵐ[6].*Pᴳ[30,:]) + 
            bid_SS_SG_30_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_30_bus_23 ) + bid_SS_SG_30_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_30_bus_24 )

# obj of VSGs for ancillary service cost
cost_VSG = bid_SS_VSG_bus_23 *sum( c_SS_bus_23 .* ΔΓ_VSG_bus_23 ) + bid_SS_VSG_bus_24 *sum( c_SS_bus_24 .* ΔΓ_VSG_bus_24 )

# obj of IBGs for ancillary service cost
@variable(model_primal, cost_IBG_23[1:T] )
@variable(model_primal, cost_IBG_24[1:T] )

            
for t in 1:T
    @constraint(model_primal, cost_IBG_23[t] == bid_AP_IBG_23_bus_23 * c_AP_bus_23[t] * Pᴳ[23,t]  + bid_AP_IBG_23_bus_24 * c_AP_bus_24[t] *  ( sh_2- ( K_c_gc[2,1] *s_g2_P23[t] +K_c_gc[2,2] *s_g3_P23[t] +K_c_gc[2,3] *s_g4_P23[t] 
            +K_c_gc[2,4] *s_g5_P23[t] +K_c_gc[2,5] *s_g27_P23[t] +K_c_gc[2,6] *s_g30_P23[t]    +α_VSG*K_c_gv[2,1]*Pᴳ[23,t]
            + sum(K_c_m[2,1:15] .*x_ηₘ_1_P23[:,t]) +  sum(K_c_m[2,16:end] .*x_ηₘ_2_P23[:,t] ) ) )    ) 

                                        +bid_RP_IBG_23_bus_23 * c_RP_bus_23[t] * Qᴳ[23,t] + bid_RP_IBG_23_bus_24 * c_RP_bus_24[t] *  ( sh_2- ( K_c_gc[2,1] *s_g2_Q23[t] +K_c_gc[2,2] *s_g3_Q23[t] +K_c_gc[2,3] *s_g4_Q23[t] 
            +K_c_gc[2,4] *s_g5_Q23[t] +K_c_gc[2,5] *s_g27_Q23[t] +K_c_gc[2,6] *s_g30_Q23[t]    +α_VSG*K_c_gv[2,1]*Qᴳ[23,t]
            + sum(K_c_m[2,1:15].* x_ηₘ_1_Q23[:,t]) +  sum(K_c_m[2,16:end].* x_ηₘ_2_Q23[:,t]) ) ) 

    @constraint(model_primal, cost_IBG_24[t] == bid_AP_IBG_24_bus_24 * c_AP_bus_24[t] * Pᴳ[24,t]  + bid_AP_IBG_24_bus_23 * c_AP_bus_23[t] *  ( sh_1- ( K_c_gc[1,1] *s_g2_P24[t] +K_c_gc[1,2] *s_g3_P24[t] +K_c_gc[1,3] *s_g4_P24[t] 
            +K_c_gc[1,4] *s_g5_P24[t] +K_c_gc[1,5] *s_g27_P24[t] +K_c_gc[1,6] *s_g30_P24[t]    +α_VSG*K_c_gv[1,1]*Pᴳ[24,t] 
            + sum(K_c_m[1,1:15] .*x_ηₘ_1_P24[:,t]) +  sum(K_c_m[1,16:end] .*x_ηₘ_2_P24[:,t] ) ) )    ) 

                                        +bid_RP_IBG_24_bus_24 * c_RP_bus_24[t] * Qᴳ[24,t] + bid_RP_IBG_24_bus_23 * c_RP_bus_23[t] *  ( sh_1- ( K_c_gc[1,1] *s_g2_Q24[t] +K_c_gc[1,2] *s_g3_Q24[t] +K_c_gc[1,3] *s_g4_Q24[t] 
            +K_c_gc[1,4] *s_g5_Q24[t] +K_c_gc[1,5] *s_g27_Q24[t] +K_c_gc[1,6] *s_g30_Q24[t]    +α_VSG*K_c_gv[1,1]*Qᴳ[24,t] 
            + sum(K_c_m[1,1:15] .*x_ηₘ_1_Q24[:,t]) +  sum(K_c_m[1,16:end] .*x_ηₘ_2_Q24[:,t] ) ) )
end
      

           
            
obj=cost_SG_2 + cost_SG_3 + cost_SG_4 + cost_SG_5 + cost_SG_27 + cost_SG_30 + cost_VSG + sum(cost_IBG_23 + cost_IBG_24)



@objective(model_primal, Min, obj)  

set_optimizer(model_primal , Gurobi.Optimizer)
optimize!(model_primal)


Value_Binary_obj= objective_value(model_primal)















#------------------------------------------------------------------------------------------------------
# Define Dual of LL model
#------------------------------------------------------------------------------------------------------

model_dual= Model()                     

@variable(model_dual, yˢᴳ²[1:T]>=0)          # status of SGs, buses:2,3,4,5,27,30.    
@variable(model_dual, yˢᴳ³[1:T]>=0)            
@variable(model_dual, yˢᴳ⁴[1:T]>=0)           
@variable(model_dual, yˢᴳ⁵[1:T]>=0)           
@variable(model_dual, yˢᴳ²⁷[1:T]>=0)            
@variable(model_dual, yˢᴳ³⁰[1:T]>=0)  

for t in 1:T
    @constraint(model_dual, yˢᴳ²[t]<=1)  
    @constraint(model_dual, yˢᴳ³[t]<=1)  
    @constraint(model_dual, yˢᴳ⁴[t]<=1)  
    @constraint(model_dual, yˢᴳ⁵[t]<=1)  
    @constraint(model_dual, yˢᴳ²⁷[t]<=1)  
    @constraint(model_dual, yˢᴳ³⁰[t]<=1)  
end

@variable(model_dual, Cᵁ²[1:T]>=0)                           # startup costs and shutdown costs for SGs                                            
@variable(model_dual, Cᵁ³[1:T]>=0)                               
@variable(model_dual, Cᵁ⁴[1:T]>=0)                                
@variable(model_dual, Cᵁ⁵[1:T]>=0)                                 
@variable(model_dual, Cᵁ²⁷[1:T]>=0)                                 
@variable(model_dual, Cᵁ³⁰[1:T]>=0)                 

#@variable(model_dual, Cᴰ²[1:T]>=0)  
#@variable(model_dual, Cᴰ³[1:T]>=0)
#@variable(model_dual, Cᴰ⁴[1:T]>=0)
#@variable(model_dual, Cᴰ⁵[1:T]>=0)
#@variable(model_dual, Cᴰ²⁷[1:T]>=0)
#@variable(model_dual, Cᴰ³⁰[1:T]>=0)

@variable(model_dual, Pᴳ[1:30,1:T]) 
@variable(model_dual, Qᴳ[1:30,1:T])  


for t in 2:T                                                # startup costs and shutdown costs occured at hour t
    @constraint(model_dual, Cᵁ²[t]>=(yˢᴳ²[t]-yˢᴳ²[t-1])*Kˢᵗ[1])        
    @constraint(model_dual, Cᵁ³[t]>=(yˢᴳ³[t]-yˢᴳ³[t-1])*Kˢᵗ[2]) 
    @constraint(model_dual, Cᵁ⁴[t]>=(yˢᴳ⁴[t]-yˢᴳ⁴[t-1])*Kˢᵗ[3])
    @constraint(model_dual, Cᵁ⁵[t]>=(yˢᴳ⁵[t]-yˢᴳ⁵[t-1])*Kˢᵗ[4])
    @constraint(model_dual, Cᵁ²⁷[t]>=(yˢᴳ²⁷[t]-yˢᴳ²⁷[t-1])*Kˢᵗ[5])
    @constraint(model_dual, Cᵁ³⁰[t]>=(yˢᴳ³⁰[t]-yˢᴳ³⁰[t-1])*Kˢᵗ[6])

    #@constraint(model_dual, Cᴰ²[t]>=(yˢᴳ²[t-1]-yˢᴳ²[t])*Kˢʰ[1]) 
    #@constraint(model_dual, Cᴰ³[t]>=(yˢᴳ³[t-1]-yˢᴳ³[t])*Kˢʰ[2])
    #@constraint(model_dual, Cᴰ⁴[t]>=(yˢᴳ⁴[t-1]-yˢᴳ⁴[t])*Kˢʰ[3])
    #@constraint(model_dual, Cᴰ⁵[t]>=(yˢᴳ⁵[t-1]-yˢᴳ⁵[t])*Kˢʰ[4])
    #@constraint(model_dual, Cᴰ²⁷[t]>=(yˢᴳ²⁷[t-1]-yˢᴳ²⁷[t])*Kˢʰ[5])
    #@constraint(model_dual, Cᴰ³⁰[t]>=(yˢᴳ³⁰[t-1]-yˢᴳ³⁰[t])*Kˢʰ[6])

end   



#---------------------------------------------------------------------------------------------------------------------------
#  AC power flow constraints
#---------------------------------------------------------------------------------------------------------------------------

for t in 1:T                                            # bounds for the active and reactive output of each bus
    @constraint(model_dual, Pᴳ[2,t] <= yˢᴳ²[t]*Pˢᴳₘₐₓ[1])    # bounds for active power of SGs: buses 2 3 4 5 27 30           
    @constraint(model_dual, yˢᴳ²[t]*Pˢᴳₘᵢₙ[1] <= Pᴳ[2,t])       
    @constraint(model_dual, Pᴳ[3,t] <= yˢᴳ³[t]*Pˢᴳₘₐₓ[2])       
    @constraint(model_dual, yˢᴳ³[t]*Pˢᴳₘᵢₙ[2] <= Pᴳ[3,t])         
    @constraint(model_dual, Pᴳ[4,t] <= yˢᴳ⁴[t]*Pˢᴳₘₐₓ[3])       
    @constraint(model_dual, yˢᴳ⁴[t]*Pˢᴳₘᵢₙ[3] <= Pᴳ[4,t])    
    @constraint(model_dual, Pᴳ[5,t] <= yˢᴳ⁵[t]*Pˢᴳₘₐₓ[4])       
    @constraint(model_dual, yˢᴳ⁵[t]*Pˢᴳₘᵢₙ[4] <= Pᴳ[5,t])
    @constraint(model_dual, Pᴳ[27,t] <= yˢᴳ²⁷[t]*Pˢᴳₘₐₓ[5])       
    @constraint(model_dual, yˢᴳ²⁷[t]*Pˢᴳₘᵢₙ[5] <=Pᴳ[27,t])
    @constraint(model_dual, Pᴳ[30,t] <= yˢᴳ³⁰[t]*Pˢᴳₘₐₓ[6])       
    @constraint(model_dual, yˢᴳ³⁰[t]*Pˢᴳₘᵢₙ[6] <= Pᴳ[30,t])

    @constraint(model_dual, Pᴳ[1,t] <= α_VSG*Pⱽˢᴳₘₐₓ[1])            # bounds for active power of VSGs and IBGs: buses 1 23 24           
    @constraint(model_dual, Pⱽˢᴳₘᵢₙ[1] <= Pᴳ[1,t])       
    @constraint(model_dual, Pᴳ[23,t] <= α_IBG*Pᴵᴮᴳₘₐₓ[1])       
    @constraint(model_dual, Pᴵᴮᴳₘᵢₙ[1] <= Pᴳ[23,t])         
    @constraint(model_dual, Pᴳ[24,t] <= α_IBG*Pᴵᴮᴳₘₐₓ[2])       
    @constraint(model_dual, Pᴵᴮᴳₘᵢₙ[2] <= Pᴳ[24,t])   

    @constraint(model_dual, Qᴳ[2,t] <= yˢᴳ²[t]*Qˢᴳₘₐₓ[1])    # bounds for REactive power of SGs: buses 2 3 4 5 27 30           
    @constraint(model_dual, yˢᴳ²[t]*Qˢᴳₘᵢₙ[1] <= Qᴳ[2,t])       
    @constraint(model_dual, Qᴳ[3,t] <= yˢᴳ³[t]*Qˢᴳₘₐₓ[2])       
    @constraint(model_dual, yˢᴳ³[t]*Qˢᴳₘᵢₙ[2] <= Qᴳ[3,t])         
    @constraint(model_dual, Qᴳ[4,t] <= yˢᴳ⁴[t]*Qˢᴳₘₐₓ[3])       
    @constraint(model_dual, yˢᴳ⁴[t]*Qˢᴳₘᵢₙ[3] <= Qᴳ[4,t])    
    @constraint(model_dual, Qᴳ[5,t] <= yˢᴳ⁵[t]*Qˢᴳₘₐₓ[4])       
    @constraint(model_dual, yˢᴳ⁵[t]*Qˢᴳₘᵢₙ[4] <= Qᴳ[5,t])
    @constraint(model_dual, Qᴳ[27,t] <= yˢᴳ²⁷[t]*Qˢᴳₘₐₓ[5])       
    @constraint(model_dual, yˢᴳ²⁷[t]*Qˢᴳₘᵢₙ[5] <=Qᴳ[27,t])
    @constraint(model_dual, Qᴳ[30,t] <= yˢᴳ³⁰[t]*Qˢᴳₘₐₓ[6])       
    @constraint(model_dual, yˢᴳ³⁰[t]*Qˢᴳₘᵢₙ[6] <= Qᴳ[30,t])           

    @constraint(model_dual, Qᴳ[1,t] <= Qⱽˢᴳₘₐₓ[1])       
    @constraint(model_dual, Qⱽˢᴳₘᵢₙ[1] <= Qᴳ[1,t])

    @constraint(model_dual, Qᴳ[23,t] <= ratio*Qᴵᴮᴳₘₐₓ[1])       
    @constraint(model_dual, Qᴵᴮᴳₘᵢₙ[1] <= Qᴳ[23,t])
    @constraint(model_dual, Qᴳ[24,t] <= ratio*Qᴵᴮᴳₘₐₓ[2])       
    @constraint(model_dual, Qᴵᴮᴳₘᵢₙ[2] <= Qᴳ[24,t])

    @constraint(model_dual, [ 1, Pᴳ[2,t] /(Sᵐᵃˣ_gc[1]), Qᴳ[2,t] /(Sᵐᵃˣ_gc[1]) ] in SecondOrderCone()  )
    @constraint(model_dual, [ 1, Pᴳ[3,t] /(Sᵐᵃˣ_gc[2]), Qᴳ[3,t] /(Sᵐᵃˣ_gc[2]) ] in SecondOrderCone()  )
    @constraint(model_dual, [ 1, Pᴳ[4,t] /(Sᵐᵃˣ_gc[3]), Qᴳ[4,t] /(Sᵐᵃˣ_gc[3]) ] in SecondOrderCone()  )
    @constraint(model_dual, [ 1, Pᴳ[5,t] /(Sᵐᵃˣ_gc[4]), Qᴳ[5,t] /(Sᵐᵃˣ_gc[4]) ] in SecondOrderCone()  )
    @constraint(model_dual, [ 1, Pᴳ[27,t] /(Sᵐᵃˣ_gc[5]), Qᴳ[27,t] /(Sᵐᵃˣ_gc[5]) ] in SecondOrderCone()  )
    @constraint(model_dual, [ 1, Pᴳ[30,t] /(Sᵐᵃˣ_gc[6]), Qᴳ[30,t] /(Sᵐᵃˣ_gc[6]) ] in SecondOrderCone()  )

    @constraint(model_dual, [ 1, Pᴳ[1,t] /(Sᵐᵃˣ_gv[1]), Qᴳ[1,t] /(Sᵐᵃˣ_gv[1]) ] in SecondOrderCone()  )

    @constraint(model_dual, [ 1, Pᴳ[23,t] /(Sᵐᵃˣ_c[1]), Qᴳ[23,t] /(Sᵐᵃˣ_c[1]) ] in SecondOrderCone()  )
    @constraint(model_dual, [ 1, Pᴳ[24,t] /(Sᵐᵃˣ_c[2]), Qᴳ[24,t] /(Sᵐᵃˣ_c[2]) ] in SecondOrderCone()  )

    gens_SGs = [2, 3, 4, 5, 27, 30]    # buses with SGs
    gens_VSGs = [1]                    # buses with VSGs
    gens_IBGs = [23, 24]               # buses with IBGs gens_all = union(gens_SGs, gens_VSGs, gens_IBGs)

    for i in 1:bus_num                 # no power is injected into the non-generation bus 
        if !(i in gens_SGs) && !(i in gens_VSGs) && !(i in gens_IBGs)
            @constraint(model_dual, Pᴳ[i, t] == 0)
            @constraint(model_dual, Qᴳ[i, t] == 0)
        end
    end
end

@variable(model_dual, Pˡⁱⁿᵉ[1:30,1:30,1:T])
@variable(model_dual, Qˡⁱⁿᵉ[1:30,1:30,1:T])

for i in 1:bus_num
    for t in 1:T

        @constraint(model_dual, Pᴳ[i,t] - Pᴰ[i,t] ==  sum( Pˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # p_i^g-p_i^d == G_{i,i}*c_{i,i}+...   eq. 50(a)
        @constraint(model_dual, Qᴳ[i,t] - Qᴰ[i,t] ==  sum( Qˡⁱⁿᵉ[i,j,t] for j in bus_connection[i] ) )    # q_i^g-q_i^d == -B_{i,i}*c_{i,i}+...  eq. 50(b)

        for j in bus_connection[i]
            @constraint(model_dual, Pˡⁱⁿᵉ[i,j,t] + Pˡⁱⁿᵉ[j,i,t] == 0)
            @constraint(model_dual, Qˡⁱⁿᵉ[i,j,t] + Qˡⁱⁿᵉ[j,i,t] == 0)

            @constraint(model_dual, Pˡⁱⁿᵉ[i,j,t] <= S_lines[i,j])
            @constraint(model_dual, -Pˡⁱⁿᵉ[i,j,t] <= S_lines[i,j])
            @constraint(model_dual, Qˡⁱⁿᵉ[i,j,t] <= S_lines[i,j])
            @constraint(model_dual, -Qˡⁱⁿᵉ[i,j,t] <= S_lines[i,j])
            @constraint(model_dual, Pˡⁱⁿᵉ[i,j,t] + Qˡⁱⁿᵉ[i,j,t] <= S_lines[i,j]*sqrt(2))
            @constraint(model_dual, -Pˡⁱⁿᵉ[i,j,t] + Qˡⁱⁿᵉ[i,j,t] <= S_lines[i,j]*sqrt(2))
            @constraint(model_dual, Pˡⁱⁿᵉ[i,j,t] - Qˡⁱⁿᵉ[i,j,t] <= S_lines[i,j]*sqrt(2))
            @constraint(model_dual, -Pˡⁱⁿᵉ[i,j,t] - Qˡⁱⁿᵉ[i,j,t] <= S_lines[i,j]*sqrt(2))

        end    
    end
end 




#--------------------------------
#  Voltage stability constraints
#--------------------------------

@variable(model_dual, P_23_eq[1:T])              # equivalent P and Q of IBGs, buses:23,24    and Γ_c  
@variable(model_dual, P_24_eq[1:T]) 
@variable(model_dual, Q_23_eq[1:T])                 
@variable(model_dual, Q_24_eq[1:T])
@variable(model_dual, Γ_23[1:T])
@variable(model_dual, Γ_24[1:T])

@variable(model_dual, ηₘ_1[1:15,1:T]>=0)        # product of each pair of SG 
@variable(model_dual, ηₘ_2[1:6,1:T]>=0)         # product of each pair of SG and VSG 

@variable(model_dual, s_g2_P23[1:T])
@variable(model_dual, s_g3_P23[1:T])
@variable(model_dual, s_g4_P23[1:T])
@variable(model_dual, s_g5_P23[1:T])
@variable(model_dual, s_g27_P23[1:T])
@variable(model_dual, s_g30_P23[1:T])
@variable(model_dual, s_g2_P24[1:T])
@variable(model_dual, s_g3_P24[1:T])
@variable(model_dual, s_g4_P24[1:T])
@variable(model_dual, s_g5_P24[1:T])
@variable(model_dual, s_g27_P24[1:T])
@variable(model_dual, s_g30_P24[1:T])

@variable(model_dual, s_g2_Q23[1:T])
@variable(model_dual, s_g3_Q23[1:T])
@variable(model_dual, s_g4_Q23[1:T])
@variable(model_dual, s_g5_Q23[1:T])
@variable(model_dual, s_g27_Q23[1:T])
@variable(model_dual, s_g30_Q23[1:T])
@variable(model_dual, s_g2_Q24[1:T])
@variable(model_dual, s_g3_Q24[1:T])
@variable(model_dual, s_g4_Q24[1:T])
@variable(model_dual, s_g5_Q24[1:T])
@variable(model_dual, s_g27_Q24[1:T])
@variable(model_dual, s_g30_Q24[1:T])

@variable(model_dual, x_ηₘ_1_P23[1:15,1:T])
@variable(model_dual, x_ηₘ_1_P24[1:15,1:T])
@variable(model_dual, x_ηₘ_1_Q23[1:15,1:T])
@variable(model_dual, x_ηₘ_1_Q24[1:15,1:T])
@variable(model_dual, x_ηₘ_2_P23[1:6,1:T])
@variable(model_dual, x_ηₘ_2_P24[1:6,1:T])
@variable(model_dual, x_ηₘ_2_Q23[1:6,1:T])
@variable(model_dual, x_ηₘ_2_Q24[1:6,1:T])

for t in 1:T
    # 15+6
    @constraint(model_dual, ηₘ_1[1,t]>=yˢᴳ²[t]+yˢᴳ³[t]-1)  
    @constraint(model_dual, ηₘ_1[1,t]<=yˢᴳ²[t])
    @constraint(model_dual, ηₘ_1[1,t]<=yˢᴳ³[t])
    @constraint(model_dual, ηₘ_1[2,t]>=yˢᴳ²[t]+yˢᴳ⁴[t]-1)  
    @constraint(model_dual, ηₘ_1[2,t]<=yˢᴳ²[t])
    @constraint(model_dual, ηₘ_1[2,t]<=yˢᴳ⁴[t])
    @constraint(model_dual, ηₘ_1[3,t]>=yˢᴳ²[t]+yˢᴳ⁵[t]-1)  
    @constraint(model_dual, ηₘ_1[3,t]<=yˢᴳ²[t])
    @constraint(model_dual, ηₘ_1[3,t]<=yˢᴳ⁵[t])
    @constraint(model_dual, ηₘ_1[4,t]>=yˢᴳ²[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model_dual, ηₘ_1[4,t]<=yˢᴳ²[t])
    @constraint(model_dual, ηₘ_1[4,t]<=yˢᴳ²⁷[t])
    @constraint(model_dual, ηₘ_1[5,t]>=yˢᴳ²[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model_dual, ηₘ_1[5,t]<=yˢᴳ²[t])
    @constraint(model_dual, ηₘ_1[5,t]<=yˢᴳ³⁰[t])

    @constraint(model_dual, ηₘ_1[6,t]>=yˢᴳ³[t]+yˢᴳ⁴[t]-1)  
    @constraint(model_dual, ηₘ_1[6,t]<=yˢᴳ³[t])
    @constraint(model_dual, ηₘ_1[6,t]<=yˢᴳ⁴[t])
    @constraint(model_dual, ηₘ_1[7,t]>=yˢᴳ³[t]+yˢᴳ⁵[t]-1)  
    @constraint(model_dual, ηₘ_1[7,t]<=yˢᴳ³[t])
    @constraint(model_dual, ηₘ_1[7,t]<=yˢᴳ⁵[t])
    @constraint(model_dual, ηₘ_1[8,t]>=yˢᴳ³[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model_dual, ηₘ_1[8,t]<=yˢᴳ³[t])
    @constraint(model_dual, ηₘ_1[8,t]<=yˢᴳ²⁷[t])
    @constraint(model_dual, ηₘ_1[9,t]>=yˢᴳ³[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model_dual, ηₘ_1[9,t]<=yˢᴳ³[t])
    @constraint(model_dual, ηₘ_1[9,t]<=yˢᴳ³⁰[t])

    @constraint(model_dual, ηₘ_1[10,t]>=yˢᴳ⁴[t]+yˢᴳ⁵[t]-1)  
    @constraint(model_dual, ηₘ_1[10,t]<=yˢᴳ⁴[t])
    @constraint(model_dual, ηₘ_1[10,t]<=yˢᴳ⁵[t])
    @constraint(model_dual, ηₘ_1[11,t]>=yˢᴳ⁴[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model_dual, ηₘ_1[11,t]<=yˢᴳ⁴[t])
    @constraint(model_dual, ηₘ_1[11,t]<=yˢᴳ²⁷[t])
    @constraint(model_dual, ηₘ_1[12,t]>=yˢᴳ⁴[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model_dual, ηₘ_1[12,t]<=yˢᴳ⁴[t])
    @constraint(model_dual, ηₘ_1[12,t]<=yˢᴳ³⁰[t])

    @constraint(model_dual, ηₘ_1[13,t]>=yˢᴳ⁵[t]+yˢᴳ²⁷[t]-1)  
    @constraint(model_dual, ηₘ_1[13,t]<=yˢᴳ⁵[t])
    @constraint(model_dual, ηₘ_1[13,t]<=yˢᴳ²⁷[t])
    @constraint(model_dual, ηₘ_1[14,t]>=yˢᴳ⁵[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model_dual, ηₘ_1[14,t]<=yˢᴳ⁵[t])
    @constraint(model_dual, ηₘ_1[14,t]<=yˢᴳ³⁰[t])

    @constraint(model_dual, ηₘ_1[15,t]>=yˢᴳ²⁷[t]+yˢᴳ³⁰[t]-1)  
    @constraint(model_dual, ηₘ_1[15,t]<=yˢᴳ²⁷[t])
    @constraint(model_dual, ηₘ_1[15,t]<=yˢᴳ³⁰[t])

    @constraint(model_dual, ηₘ_2[1,t]==yˢᴳ²[t]*α_VSG)
    @constraint(model_dual, ηₘ_2[2,t]==yˢᴳ³[t]*α_VSG)
    @constraint(model_dual, ηₘ_2[3,t]==yˢᴳ⁴[t]*α_VSG)
    @constraint(model_dual, ηₘ_2[4,t]==yˢᴳ⁵[t]*α_VSG)
    @constraint(model_dual, ηₘ_2[5,t]==yˢᴳ²⁷[t]*α_VSG)
    @constraint(model_dual, ηₘ_2[6,t]==yˢᴳ³⁰[t]*α_VSG)

    @constraint(model_dual, s_g2_P23[t] <=yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g2_P23[t] <= Pᴳ[23,t] +yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g2_P23[t] >=Pᴳ[23,t]+yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g2_P23[t] >= yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g3_P23[t] <=yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g3_P23[t] <= Pᴳ[23,t] +yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g3_P23[t] >=Pᴳ[23,t]+yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g3_P23[t] >= yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g4_P23[t] <=yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g4_P23[t] <= Pᴳ[23,t] +yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g4_P23[t] >=Pᴳ[23,t]+yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g4_P23[t] >= yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g5_P23[t] <=yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g5_P23[t] <= Pᴳ[23,t] +yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g5_P23[t] >=Pᴳ[23,t]+yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])   
    @constraint(model_dual, s_g5_P23[t] >= yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g27_P23[t] <=yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g27_P23[t] <= Pᴳ[23,t] +yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g27_P23[t] >=Pᴳ[23,t]+yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])   
    @constraint(model_dual, s_g27_P23[t] >= yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g30_P23[t] <=yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g30_P23[t] <= Pᴳ[23,t] +yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[1] -Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g30_P23[t] >=Pᴳ[23,t]+yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[1]-α_IBG*Pᴵᴮᴳₘₐₓ[1])   
    @constraint(model_dual, s_g30_P23[t] >= yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g2_P24[t] <=yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g2_P24[t] <= Pᴳ[24,t] +yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g2_P24[t] >=Pᴳ[24,t]+yˢᴳ²[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g2_P24[t] >= yˢᴳ²[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, s_g3_P24[t] <=yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g3_P24[t] <= Pᴳ[24,t] +yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g3_P24[t] >=Pᴳ[24,t]+yˢᴳ³[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g3_P24[t] >= yˢᴳ³[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, s_g4_P24[t] <=yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g4_P24[t] <= Pᴳ[24,t] +yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g4_P24[t] >=Pᴳ[24,t]+yˢᴳ⁴[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g4_P24[t] >= yˢᴳ⁴[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, s_g5_P24[t] <=yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g5_P24[t] <= Pᴳ[24,t] +yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g5_P24[t] >=Pᴳ[24,t]+yˢᴳ⁵[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])   
    @constraint(model_dual, s_g5_P24[t] >= yˢᴳ⁵[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, s_g27_P24[t] <=yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g27_P24[t] <= Pᴳ[24,t] +yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g27_P24[t] >=Pᴳ[24,t]+yˢᴳ²⁷[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])   
    @constraint(model_dual, s_g27_P24[t] >= yˢᴳ²⁷[t]*Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, s_g30_P24[t] <=yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g30_P24[t] <= Pᴳ[24,t] +yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[2] -Pᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g30_P24[t] >=Pᴳ[24,t]+yˢᴳ³⁰[t]*α_IBG*Pᴵᴮᴳₘₐₓ[2]-α_IBG*Pᴵᴮᴳₘₐₓ[2])   
    @constraint(model_dual, s_g30_P24[t] >= yˢᴳ³⁰[t]*Pᴵᴮᴳₘᵢₙ[2] )

    
    
    @constraint(model_dual, s_g2_Q23[t] <= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g2_Q23[t] <= Qᴳ[23,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g2_Q23[t] >= Qᴳ[23,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g2_Q23[t] >= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g3_Q23[t] <= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g3_Q23[t] <= Qᴳ[23,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g3_Q23[t] >= Qᴳ[23,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g3_Q23[t] >= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g4_Q23[t] <= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g4_Q23[t] <= Qᴳ[23,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g4_Q23[t] >= Qᴳ[23,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g4_Q23[t] >= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g5_Q23[t] <= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g5_Q23[t] <= Qᴳ[23,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g5_Q23[t] >= Qᴳ[23,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])   
    @constraint(model_dual, s_g5_Q23[t] >= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g27_Q23[t] <= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g27_Q23[t] <= Qᴳ[23,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g27_Q23[t] >= Qᴳ[23,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])   
    @constraint(model_dual, s_g27_Q23[t] >= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g30_Q23[t] <= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[1])
    @constraint(model_dual, s_g30_Q23[t] <= Qᴳ[23,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] -ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, s_g30_Q23[t] >= Qᴳ[23,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[1]-ratio*Qᴵᴮᴳₘₐₓ[1])   
    @constraint(model_dual, s_g30_Q23[t] >= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, s_g2_Q24[t] <= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g2_Q24[t] <= Qᴳ[24,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g2_Q24[t] >= Qᴳ[24,t] +yˢᴳ²[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g2_Q24[t] >= yˢᴳ²[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )   

    @constraint(model_dual, s_g3_Q24[t] <= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g3_Q24[t] <= Qᴳ[24,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g3_Q24[t] >= Qᴳ[24,t] +yˢᴳ³[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g3_Q24[t] >= yˢᴳ³[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, s_g4_Q24[t] <= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g4_Q24[t] <= Qᴳ[24,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g4_Q24[t] >= Qᴳ[24,t] +yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g4_Q24[t] >= yˢᴳ⁴[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, s_g5_Q24[t] <= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g5_Q24[t] <= Qᴳ[24,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g5_Q24[t] >= Qᴳ[24,t] +yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])   
    @constraint(model_dual, s_g5_Q24[t] >= yˢᴳ⁵[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, s_g27_Q24[t] <= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g27_Q24[t] <= Qᴳ[24,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] ) 
    @constraint(model_dual, s_g27_Q24[t] >= Qᴳ[24,t] +yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])   
    @constraint(model_dual, s_g27_Q24[t] >= yˢᴳ²⁷[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, s_g30_Q24[t] <= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[2])
    @constraint(model_dual, s_g30_Q24[t] <= Qᴳ[24,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] -ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, s_g30_Q24[t] >= Qᴳ[24,t] +yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘₐₓ[2]-ratio*Qᴵᴮᴳₘₐₓ[2])   
    @constraint(model_dual, s_g30_Q24[t] >= yˢᴳ³⁰[t]*ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, x_ηₘ_1_P23[:,t] <= ηₘ_1[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model_dual, x_ηₘ_1_P23[:,t] <= Pᴳ[23,t] .+ ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[1] .- Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, x_ηₘ_1_P23[:,t] >= Pᴳ[23,t] .+ ηₘ_1[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[1] .- α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model_dual, x_ηₘ_1_P23[:,t] >= ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, x_ηₘ_2_P23[:,t] <= ηₘ_2[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model_dual, x_ηₘ_2_P23[:,t] <= Pᴳ[23,t] .+ ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[1] .- Pᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, x_ηₘ_2_P23[:,t] >= Pᴳ[23,t] .+ ηₘ_2[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[1] .- α_IBG*Pᴵᴮᴳₘₐₓ[1] )
    @constraint(model_dual, x_ηₘ_2_P23[:,t] >= ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, x_ηₘ_1_P24[:,t] <= ηₘ_1[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model_dual, x_ηₘ_1_P24[:,t] <= Pᴳ[24,t] .+ ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[2] .- Pᴵᴮᴳₘᵢₙ[2] )    
    @constraint(model_dual, x_ηₘ_1_P24[:,t] >= Pᴳ[24,t] .+ ηₘ_1[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[2] .- α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model_dual, x_ηₘ_1_P24[:,t] >= ηₘ_1[:,t] *Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, x_ηₘ_2_P24[:,t] <= ηₘ_2[:,t] .*α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model_dual, x_ηₘ_2_P24[:,t] <= Pᴳ[24,t] .+ ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[2] .- Pᴵᴮᴳₘᵢₙ[2] )    
    @constraint(model_dual, x_ηₘ_2_P24[:,t] >= Pᴳ[24,t] .+ ηₘ_2[:,t] *α_IBG*Pᴵᴮᴳₘₐₓ[2] .- α_IBG*Pᴵᴮᴳₘₐₓ[2] )
    @constraint(model_dual, x_ηₘ_2_P24[:,t] >= ηₘ_2[:,t] *Pᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, x_ηₘ_1_Q23[:,t] <= ηₘ_1[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[1] )  
    @constraint(model_dual, x_ηₘ_1_Q23[:,t] <= Qᴳ[23,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] .- ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, x_ηₘ_1_Q23[:,t] >= Qᴳ[23,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘₐₓ[1] .- ratio*Qᴵᴮᴳₘₐₓ[1] )
    @constraint(model_dual, x_ηₘ_1_Q23[:,t] >= ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] )
    
    @constraint(model_dual, x_ηₘ_2_Q23[:,t] <= ηₘ_2[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[1] )  
    @constraint(model_dual, x_ηₘ_2_Q23[:,t] <= Qᴳ[23,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] .- ratio*Qᴵᴮᴳₘᵢₙ[1] )
    @constraint(model_dual, x_ηₘ_2_Q23[:,t] >= Qᴳ[23,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘₐₓ[1] .- ratio*Qᴵᴮᴳₘₐₓ[1] )
    @constraint(model_dual, x_ηₘ_2_Q23[:,t] >= ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[1] )

    @constraint(model_dual, x_ηₘ_1_Q24[:,t] <= ηₘ_1[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[2] )  
    @constraint(model_dual, x_ηₘ_1_Q24[:,t] <= Qᴳ[24,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] .- ratio*Qᴵᴮᴳₘᵢₙ[2] )    
    @constraint(model_dual, x_ηₘ_1_Q24[:,t] >= Qᴳ[24,t] .+ ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘₐₓ[2] .- ratio*Qᴵᴮᴳₘₐₓ[2] )
    @constraint(model_dual, x_ηₘ_1_Q24[:,t] >= ηₘ_1[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, x_ηₘ_2_Q24[:,t] <= ηₘ_2[:,t] .*ratio*Qᴵᴮᴳₘₐₓ[2] )  
    @constraint(model_dual, x_ηₘ_2_Q24[:,t] <= Qᴳ[24,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] .- ratio*Qᴵᴮᴳₘᵢₙ[2] )
    @constraint(model_dual, x_ηₘ_2_Q24[:,t] >= Qᴳ[24,t] .+ ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘₐₓ[2] .- ratio*Qᴵᴮᴳₘₐₓ[2] )
    @constraint(model_dual, x_ηₘ_2_Q24[:,t] >= ηₘ_2[:,t] *ratio*Qᴵᴮᴳₘᵢₙ[2] )

    @constraint(model_dual, P_23_eq[t]== Pᴳ[23,t]+  ( sh_1*Pᴳ[24,t] - ( K_c_gc[1,1] *s_g2_P24[t] +K_c_gc[1,2] *s_g3_P24[t] +K_c_gc[1,3] *s_g4_P24[t] 
    +K_c_gc[1,4] *s_g5_P24[t] +K_c_gc[1,5] *s_g27_P24[t] +K_c_gc[1,6] *s_g30_P24[t]    +α_VSG*K_c_gv[1,1]*Pᴳ[24,t] 
    + sum(K_c_m[1,1:15] .*x_ηₘ_1_P24[:,t]) +  sum(K_c_m[1,16:end] .*x_ηₘ_2_P24[:,t] ) ) )    )

    @constraint(model_dual, P_24_eq[t]== Pᴳ[24,t]+  ( sh_2*Pᴳ[23,t] - ( K_c_gc[2,1] *s_g2_P23[t] +K_c_gc[2,2] *s_g3_P23[t] +K_c_gc[2,3] *s_g4_P23[t] 
    +K_c_gc[2,4] *s_g5_P23[t] +K_c_gc[2,5] *s_g27_P23[t] +K_c_gc[2,6] *s_g30_P23[t]    +α_VSG*K_c_gv[2,1]*Pᴳ[23,t]
    + sum(K_c_m[2,1:15] .*x_ηₘ_1_P23[:,t]) +  sum(K_c_m[2,16:end] .*x_ηₘ_2_P23[:,t] ) ) )  )

    @constraint(model_dual, Q_23_eq[t]== Qᴳ[23,t]+  ( sh_1*Qᴳ[24,t] - ( K_c_gc[1,1] *s_g2_Q24[t] +K_c_gc[1,2] *s_g3_Q24[t] +K_c_gc[1,3] *s_g4_Q24[t] 
    +K_c_gc[1,4] *s_g5_Q24[t] +K_c_gc[1,5] *s_g27_Q24[t] +K_c_gc[1,6] *s_g30_Q24[t]    +α_VSG*K_c_gv[1,1]*Qᴳ[24,t] 
    + sum(K_c_m[1,1:15] .*x_ηₘ_1_Q24[:,t]) +  sum(K_c_m[1,16:end] .*x_ηₘ_2_Q24[:,t] ) ) )    )

    @constraint(model_dual, Q_24_eq[t]== Qᴳ[24,t]+  ( sh_2*Qᴳ[23,t] - ( K_c_gc[2,1] *s_g2_Q23[t] +K_c_gc[2,2] *s_g3_Q23[t] +K_c_gc[2,3] *s_g4_Q23[t] 
    +K_c_gc[2,4] *s_g5_Q23[t] +K_c_gc[2,5] *s_g27_Q23[t] +K_c_gc[2,6] *s_g30_Q23[t]    +α_VSG*K_c_gv[2,1]*Qᴳ[23,t]
    + sum(K_c_m[2,1:15] .*x_ηₘ_1_Q23[:,t]) +  sum(K_c_m[2,16:end] .*x_ηₘ_2_Q23[:,t] ) ) )    )

    @constraint(model_dual, Γ_23[t]== 1/2 * ( K_c_gc[3,1] *yˢᴳ²[t] +K_c_gc[3,2] *yˢᴳ³[t] +K_c_gc[3,3] *yˢᴳ⁴[t] 
    +K_c_gc[3,4] *yˢᴳ⁵[t] +K_c_gc[3,5] *yˢᴳ²⁷[t] +K_c_gc[3,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[3,1]
    + sum(K_c_m[3,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[3,16:end] .*ηₘ_2[:,t]   ) )

    @constraint(model_dual, Γ_24[t]== 1/2 * ( K_c_gc[4,1] *yˢᴳ²[t] +K_c_gc[4,2] *yˢᴳ³[t] +K_c_gc[4,3] *yˢᴳ⁴[t] 
    +K_c_gc[4,4] *yˢᴳ⁵[t] +K_c_gc[4,5] *yˢᴳ²⁷[t] +K_c_gc[4,6] *yˢᴳ³⁰[t] )   +α_VSG*K_c_gv[4,1]
    + sum(K_c_m[4,1:15] .*ηₘ_1[:,t]) +  sum(K_c_m[4,16:end] .*ηₘ_2[:,t]   ) )

    @constraint(model_dual, [Q_23_eq[t] + Γ_23[t]*S_base, P_23_eq[t], Q_23_eq[t]] in SecondOrderCone()  )      
    @constraint(model_dual, [Q_24_eq[t] + Γ_24[t]*S_base, P_24_eq[t], Q_24_eq[t]] in SecondOrderCone()  )

end


@variable(model_dual, ΔΓ_SG_2_bus_23[1:T] )
@variable(model_dual, ΔΓ_SG_2_bus_24[1:T] )
@variable(model_dual, ΔΓ_SG_3_bus_23[1:T] )
@variable(model_dual, ΔΓ_SG_3_bus_24[1:T] )
@variable(model_dual, ΔΓ_SG_4_bus_23[1:T] )
@variable(model_dual, ΔΓ_SG_4_bus_24[1:T] )
@variable(model_dual, ΔΓ_SG_5_bus_23[1:T] )
@variable(model_dual, ΔΓ_SG_5_bus_24[1:T] )
@variable(model_dual, ΔΓ_SG_27_bus_23[1:T] )
@variable(model_dual, ΔΓ_SG_27_bus_24[1:T] )
@variable(model_dual, ΔΓ_SG_30_bus_23[1:T] )
@variable(model_dual, ΔΓ_SG_30_bus_24[1:T] )
@variable(model_dual, ΔΓ_VSG_bus_23[1:T] )
@variable(model_dual, ΔΓ_VSG_bus_24[1:T] )

for t in 1:T

    @constraint(model_dual, ΔΓ_SG_2_bus_23[t] == 1/2 * (K_c_gc[3,1] *yˢᴳ²[t] + sum(K_c_m[3,1:5] .*ηₘ_1[1:5 , t]) + sum(K_c_m[3,16] .*ηₘ_2[1 , t] ) ) )
    
    @constraint(model_dual, ΔΓ_SG_2_bus_24[t] == 1/2 * (K_c_gc[4,1] *yˢᴳ²[t] + sum(K_c_m[4,1:5] .*ηₘ_1[1:5 , t]) + sum(K_c_m[4,16] .*ηₘ_2[1 , t] ) ) )

    @constraint(model_dual, ΔΓ_SG_3_bus_23[t] == 1/2 * (K_c_gc[3,2] *yˢᴳ³[t] + sum(K_c_m[3,1] .*ηₘ_1[1 , t]) + sum(K_c_m[3,6:9] .*ηₘ_1[6:9 , t]) +  sum(K_c_m[3,17] .*ηₘ_2[2 , t])  ) )
    
    @constraint(model_dual, ΔΓ_SG_3_bus_24[t] == 1/2 * (K_c_gc[4,2] *yˢᴳ³[t] + sum(K_c_m[4,1] .*ηₘ_1[1 , t]) + sum(K_c_m[4,6:9] .*ηₘ_1[6:9 , t]) +  sum(K_c_m[4,17] .*ηₘ_2[2 , t])  ) )

    @constraint(model_dual, ΔΓ_SG_4_bus_23[t] == 1/2 * (K_c_gc[3,3] *yˢᴳ⁴[t] + sum(K_c_m[3,2] .*ηₘ_1[2 , t]) + sum(K_c_m[3,6] .*ηₘ_1[6 , t]) + sum(K_c_m[3,10:12] .*ηₘ_1[10:12 , t]) +  sum(K_c_m[3,18] .*ηₘ_2[3 , t])  ) )
   
    @constraint(model_dual, ΔΓ_SG_4_bus_24[t] == 1/2 * (K_c_gc[4,3] *yˢᴳ⁴[t] + sum(K_c_m[4,2] .*ηₘ_1[2 , t]) + sum(K_c_m[4,6] .*ηₘ_1[6 , t]) + sum(K_c_m[4,10:12] .*ηₘ_1[10:12 , t]) +  sum(K_c_m[4,18] .*ηₘ_2[3 , t])  ) )
   
    @constraint(model_dual, ΔΓ_SG_5_bus_23[t] == 1/2 * ( K_c_gc[3,4] *yˢᴳ⁵[t] + sum(K_c_m[3,3] .*ηₘ_1[3 , t]) + sum(K_c_m[3,7] .*ηₘ_1[7 , t]) + sum(K_c_m[3,10] .*ηₘ_1[10 , t]) + sum(K_c_m[3,13:14] .*ηₘ_1[13:14 , t]) +  sum(K_c_m[3,19] .*ηₘ_2[4 , t] ) ) )

    @constraint(model_dual, ΔΓ_SG_5_bus_24[t] == 1/2 * ( K_c_gc[4,4] *yˢᴳ⁵[t] + sum(K_c_m[4,3] .*ηₘ_1[3 , t]) + sum(K_c_m[4,7] .*ηₘ_1[7 , t]) + sum(K_c_m[4,10] .*ηₘ_1[10 , t]) + sum(K_c_m[4,13:14] .*ηₘ_1[13:14 , t]) +  sum(K_c_m[4,19] .*ηₘ_2[4 , t] ) ) )

    @constraint(model_dual, ΔΓ_SG_27_bus_23[t] == 1/2 * ( K_c_gc[3,5] *yˢᴳ²⁷[t] + sum(K_c_m[3,4] .*ηₘ_1[4 , t]) + sum(K_c_m[3,8] .*ηₘ_1[8 , t]) + sum(K_c_m[3,11] .*ηₘ_1[11 , t]) + sum(K_c_m[3,13] .*ηₘ_1[13 , t]) + sum(K_c_m[3,15] .*ηₘ_1[15 , t]) +  sum(K_c_m[3,20] .*ηₘ_2[5 , t] ) ) )

    @constraint(model_dual, ΔΓ_SG_27_bus_24[t] == 1/2 * ( K_c_gc[4,5] *yˢᴳ²⁷[t] + sum(K_c_m[4,4] .*ηₘ_1[4 , t]) + sum(K_c_m[4,8] .*ηₘ_1[8 , t]) + sum(K_c_m[4,11] .*ηₘ_1[11 , t]) + sum(K_c_m[4,13] .*ηₘ_1[13 , t]) + sum(K_c_m[4,15] .*ηₘ_1[15 , t]) +  sum(K_c_m[4,20] .*ηₘ_2[5 , t] ) ) )

    @constraint(model_dual, ΔΓ_SG_30_bus_23[t] == 1/2 * ( K_c_gc[3,6] *yˢᴳ³⁰[t] + sum(K_c_m[3,5] .*ηₘ_1[5 , t]) + sum(K_c_m[3,9] .*ηₘ_1[9 , t]) + sum(K_c_m[3,12] .*ηₘ_1[12 , t]) + sum(K_c_m[3,14] .*ηₘ_1[14 , t]) + sum(K_c_m[3,15] .*ηₘ_1[15 , t]) +  sum(K_c_m[3,21] .*ηₘ_2[6 , t] ) ) )

    @constraint(model_dual, ΔΓ_SG_30_bus_24[t] == 1/2 * ( K_c_gc[4,6] *yˢᴳ³⁰[t] + sum(K_c_m[4,5] .*ηₘ_1[5 , t]) + sum(K_c_m[4,9] .*ηₘ_1[9 , t]) + sum(K_c_m[4,12] .*ηₘ_1[12 , t]) + sum(K_c_m[4,14] .*ηₘ_1[14 , t]) + sum(K_c_m[4,15] .*ηₘ_1[15 , t]) +  sum(K_c_m[4,21] .*ηₘ_2[6 , t] ) ) )

    @constraint(model_dual, ΔΓ_VSG_bus_23[t] == 1/2 * ( α_VSG*K_c_gv[3,1] +  sum(K_c_m[3,16:end] .*ηₘ_2[: , t] ) ) )

    @constraint(model_dual, ΔΓ_VSG_bus_24[t] == 1/2 * ( α_VSG*K_c_gv[4,1] +  sum(K_c_m[4,16:end] .*ηₘ_2[: , t] ) ) )

end

#--------------------------------
#  Model solving
#--------------------------------
# obj of SGs for operation cost and ancillary service cost +sum(Cᴰ²)
cost_SG_2 = sum(Cᵁ²) + sum(Oⁿˡ[1].*(yˢᴳ²)) + bid_energy_SG_2 *sum(Oᵐ[1].*Pᴳ[2,:]) + 
            bid_SS_SG_2_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_2_bus_23 ) + bid_SS_SG_2_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_2_bus_24 ) 
cost_SG_3 = sum(Cᵁ³) + sum(Oⁿˡ[2].*(yˢᴳ³)) + bid_energy_SG_3 *sum(Oᵐ[2].*Pᴳ[3,:]) + 
            bid_SS_SG_3_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_3_bus_23 ) + bid_SS_SG_3_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_3_bus_24 )
cost_SG_4 = sum(Cᵁ⁴) + sum(Oⁿˡ[3].*(yˢᴳ⁴)) + bid_energy_SG_4 *sum(Oᵐ[3].*Pᴳ[4,:]) + 
            bid_SS_SG_4_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_4_bus_23 ) + bid_SS_SG_4_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_4_bus_24 )
cost_SG_5 = sum(Cᵁ⁵) + sum(Oⁿˡ[4].*(yˢᴳ⁵)) + bid_energy_SG_5 *sum(Oᵐ[4].*Pᴳ[5,:]) + 
            bid_SS_SG_5_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_5_bus_23 ) + bid_SS_SG_5_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_5_bus_24 )
cost_SG_27 = sum(Cᵁ²⁷) + sum(Oⁿˡ[5].*(yˢᴳ²⁷)) + bid_energy_SG_27 *sum(Oᵐ[5].*Pᴳ[27,:]) + 
            bid_SS_SG_27_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_27_bus_23 ) + bid_SS_SG_27_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_27_bus_24 )
cost_SG_30 = sum(Cᵁ³⁰) + sum(Oⁿˡ[6].*(yˢᴳ³⁰)) + bid_energy_SG_30 *sum(Oᵐ[6].*Pᴳ[30,:]) + 
            bid_SS_SG_30_bus_23 *sum(c_SS_bus_23.* ΔΓ_SG_30_bus_23 ) + bid_SS_SG_30_bus_24 *sum(c_SS_bus_24.* ΔΓ_SG_30_bus_24 )

# obj of VSGs for ancillary service cost
cost_VSG = bid_SS_VSG_bus_23 *sum( c_SS_bus_23 .* ΔΓ_VSG_bus_23 ) + bid_SS_VSG_bus_24 *sum( c_SS_bus_24 .* ΔΓ_VSG_bus_24 )

# obj of IBGs for ancillary service cost
@variable(model_dual, cost_IBG_23[1:T] )
@variable(model_dual, cost_IBG_24[1:T] )

for t in 1:T
    @constraint(model_dual, cost_IBG_23[t] == bid_AP_IBG_23_bus_23 * c_AP_bus_23[t] * Pᴳ[23,t]  + bid_AP_IBG_23_bus_24 * c_AP_bus_24[t] *  ( sh_2*Pᴳ[23,t] - ( K_c_gc[2,1] *s_g2_P23[t] +K_c_gc[2,2] *s_g3_P23[t] +K_c_gc[2,3] *s_g4_P23[t] 
            +K_c_gc[2,4] *s_g5_P23[t] +K_c_gc[2,5] *s_g27_P23[t] +K_c_gc[2,6] *s_g30_P23[t]    +α_VSG*K_c_gv[2,1]*Pᴳ[23,t]
            + sum(K_c_m[2,1:15] .*x_ηₘ_1_P23[:,t]) +  sum(K_c_m[2,16:end] .*x_ηₘ_2_P23[:,t] ) ) )    ) 

                                        +bid_RP_IBG_23_bus_23 * c_RP_bus_23[t] * Qᴳ[23,t] + bid_RP_IBG_23_bus_24 * c_RP_bus_24[t] *  ( sh_2*Qᴳ[23,t] - ( K_c_gc[2,1] *s_g2_Q23[t] +K_c_gc[2,2] *s_g3_Q23[t] +K_c_gc[2,3] *s_g4_Q23[t] 
            +K_c_gc[2,4] *s_g5_Q23[t] +K_c_gc[2,5] *s_g27_Q23[t] +K_c_gc[2,6] *s_g30_Q23[t]    +α_VSG*K_c_gv[2,1]*Qᴳ[23,t]
            + sum(K_c_m[2,1:15].* x_ηₘ_1_Q23[:,t]) +  sum(K_c_m[2,16:end].* x_ηₘ_2_Q23[:,t]) ) ) 

    @constraint(model_dual, cost_IBG_24[t] == bid_AP_IBG_24_bus_24 * c_AP_bus_24[t] * Pᴳ[24,t]  + bid_AP_IBG_24_bus_23 * c_AP_bus_23[t] *  ( sh_1*Pᴳ[24,t] - ( K_c_gc[1,1] *s_g2_P24[t] +K_c_gc[1,2] *s_g3_P24[t] +K_c_gc[1,3] *s_g4_P24[t] 
            +K_c_gc[1,4] *s_g5_P24[t] +K_c_gc[1,5] *s_g27_P24[t] +K_c_gc[1,6] *s_g30_P24[t]    +α_VSG*K_c_gv[1,1]*Pᴳ[24,t] 
            + sum(K_c_m[1,1:15] .*x_ηₘ_1_P24[:,t]) +  sum(K_c_m[1,16:end] .*x_ηₘ_2_P24[:,t] ) ) )    ) 

                                        +bid_RP_IBG_24_bus_24 * c_RP_bus_24[t] * Qᴳ[24,t] + bid_RP_IBG_24_bus_23 * c_RP_bus_23[t] *  ( sh_1*Qᴳ[24,t] - ( K_c_gc[1,1] *s_g2_Q24[t] +K_c_gc[1,2] *s_g3_Q24[t] +K_c_gc[1,3] *s_g4_Q24[t] 
            +K_c_gc[1,4] *s_g5_Q24[t] +K_c_gc[1,5] *s_g27_Q24[t] +K_c_gc[1,6] *s_g30_Q24[t]    +α_VSG*K_c_gv[1,1]*Qᴳ[24,t] 
            + sum(K_c_m[1,1:15] .*x_ηₘ_1_Q24[:,t]) +  sum(K_c_m[1,16:end] .*x_ηₘ_2_Q24[:,t] ) ) )
end
      

           
            
obj=cost_SG_2 + cost_SG_3 + cost_SG_4 + cost_SG_5 + cost_SG_27 + cost_SG_30 + cost_VSG + sum(cost_IBG_23 + cost_IBG_24)



@objective(model_dual, Min, obj)  


dual=dualize(model_dual)

set_optimizer(model_dual, Gurobi.Optimizer)
set_optimizer_attribute(model_dual, "QCPDual", 1)
optimize!(model_dual)


Value_Relaxed_obj= objective_value(model_dual)

DG = (Value_Binary_obj - Value_Relaxed_obj)/Value_Binary_obj

energy_prices_P=zeros(bus_num,T)
energy_prices_P = [dual(Power_balance_P[i][t]) for i in 1:bus_num, t in 1:T]

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

return DG, energy_prices_P, -lambda_1[:,1], -lambda_1[:,2], mu[:,1], mu[:,2], -lambda_2[:,1]+mu[:,1], -lambda_2[:,2]+mu[:,2]

end