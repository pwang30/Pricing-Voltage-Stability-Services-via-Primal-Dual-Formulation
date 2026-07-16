function settlement_dispatchable( Pᴳ,Qᴳ,mu,lambda_2,energy_prices_P, 
                    yˢᴳ²,yˢᴳ³,yˢᴳ⁴,yˢᴳ⁵,yˢᴳ²⁷,yˢᴳ³⁰,
                    Cᵁ²,Cᵁ³,Cᵁ⁴,Cᵁ⁵,Cᵁ²⁷,Cᵁ³⁰,
                    Oᵐ,Oⁿˡ,α_VSG,Gᵥ,IBG,S_Base,
                    K_c_gc, K_c_m, ηₘ_1, ηₘ_2)

μ_23 = mu[:, 1]  
μ_24 = mu[:, 2]  

#------------------------------
# Calculate ENERGY profits of units
#------------------------------

energy_profitˢᴳ² = Pᴳ[1,:] .*energy_prices_P .-Pᴳ[1,:] *Oᵐ[1] .-Cᵁ²  .-yˢᴳ²*Oⁿˡ[1]
energy_profitˢᴳ³ = Pᴳ[2,:] .*energy_prices_P .-Pᴳ[2,:] *Oᵐ[2] .-Cᵁ³  .-yˢᴳ³*Oⁿˡ[2]
energy_profitˢᴳ⁴ = Pᴳ[3,:] .*energy_prices_P .-Pᴳ[3,:] *Oᵐ[3] .-Cᵁ⁴  .-yˢᴳ⁴*Oⁿˡ[3]
energy_profitˢᴳ⁵ = Pᴳ[4,:] .*energy_prices_P .-Pᴳ[4,:] *Oᵐ[4] .-Cᵁ⁵  .-yˢᴳ⁵*Oⁿˡ[4]
energy_profitˢᴳ²⁷ = Pᴳ[5,:] .*energy_prices_P .-Pᴳ[5,:] *Oᵐ[5] .-Cᵁ²⁷  .-yˢᴳ²⁷*Oⁿˡ[5]
energy_profitˢᴳ³⁰ = Pᴳ[6,:] .*energy_prices_P .-Pᴳ[6,:] *Oᵐ[6] .-Cᵁ³⁰  .-yˢᴳ³⁰*Oⁿˡ[6]
energy_profit_VSG_1=Pᴳ[7,:] .*energy_prices_P
energy_profit_IBG_23=Pᴳ[8,:] .*energy_prices_P
energy_profit_IBG_24=Pᴳ[9,:] .*energy_prices_P


#------------------------------
# Calculate VS service revenues of units
#------------------------------

#----revenue for IBGs for providing reactive power support

revenue_IBG_23_give_reactive=zeros(T)
revenue_IBG_24_give_reactive=zeros(T)

for t in 1:T
    revenue_IBG_23_give_reactive[t]=  (-lambda_2[t,1] + mu[t,1]  )* Qᴳ[8,t] 
    revenue_IBG_24_give_reactive[t]=  (-lambda_2[t,2] + mu[t,2]  )* Qᴳ[9,t] 
end

#----revenue for (V)SGs for providing SCR support

VS_revenue_SGs=zeros(6,T)

for t in 1:T
VS_revenue_SGs[1,t] =   μ_23[t]*K_c_gc[1, 1]*S_Base *yˢᴳ²[t] + μ_24[t]*K_c_gc[2, 1]*S_Base *yˢᴳ²[t]
                                        
                                            + μ_23[t]*sum(K_c_m[1, 1:5].*ηₘ_1[1:5,t])*S_Base + μ_24[t]*sum(K_c_m[2, 1:5].*ηₘ_1[1:5,t])*S_Base 
                                        
                                            + μ_23[t]*ηₘ_2[1,t]*K_c_m[1, 16]*S_Base + μ_24[t]*ηₘ_2[1,t]*K_c_m[2, 16]*S_Base  
    
VS_revenue_SGs[2,t] =  μ_23[t]*K_c_gc[1, 2]*S_Base *yˢᴳ³[t] + μ_24[t]*K_c_gc[2, 2]*S_Base *yˢᴳ³[t]

                                            + μ_23[t]*sum(K_c_m[1, 6:9].*ηₘ_1[6:9,t])*S_Base + μ_24[t]*sum(K_c_m[2, 6:9].*ηₘ_1[6:9,t])*S_Base
                                            + μ_23[t]*ηₘ_1[1,t]*K_c_m[1, 1]*S_Base + μ_24[t]*ηₘ_1[1,t]*K_c_m[2, 1]*S_Base 

                                            + μ_23[t]*ηₘ_2[2,t]*K_c_m[1, 17]*S_Base + μ_24[t]*ηₘ_2[2,t]*K_c_m[2, 17]*S_Base   
    
VS_revenue_SGs[3,t] = μ_23[t]*K_c_gc[1, 3]*S_Base *yˢᴳ⁴[t] + μ_24[t]*K_c_gc[2, 3]*S_Base *yˢᴳ⁴[t]

                                            + μ_23[t]*sum(K_c_m[1, 10:12].*ηₘ_1[10:12,t])*S_Base + μ_24[t]*sum(K_c_m[2, 10:12].*ηₘ_1[10:12,t])*S_Base
                                            + μ_23[t]*ηₘ_1[2,t]*K_c_m[1, 2]*S_Base + μ_23[t]*ηₘ_1[6,t]*K_c_m[1, 6]*S_Base 
                                            + μ_24[t]*ηₘ_1[2,t]*K_c_m[2, 2]*S_Base + μ_24[t]*ηₘ_1[6,t]*K_c_m[2, 6]*S_Base 

                                            + μ_23[t]*ηₘ_2[3,t]*K_c_m[1, 18]*S_Base + μ_24[t]*ηₘ_2[3,t]*K_c_m[2, 18]*S_Base 
    
VS_revenue_SGs[4,t] = μ_23[t]*K_c_gc[1, 4]*S_Base *yˢᴳ⁵[t] + μ_24[t]*K_c_gc[2, 4]*S_Base *yˢᴳ⁵[t]

                                            + μ_23[t]*sum(K_c_m[1, 13:14].*ηₘ_1[13:14,t])*S_Base + μ_24[t]*sum(K_c_m[2, 13:14].*ηₘ_1[13:14,t])*S_Base
                                            + μ_23[t]*ηₘ_1[3,t]*K_c_m[1, 3]*S_Base + μ_23[t]*ηₘ_1[7,t]*K_c_m[1, 7]*S_Base + μ_23[t]*ηₘ_1[10,t]*K_c_m[1, 10]*S_Base 
                                            + μ_24[t]*ηₘ_1[3,t]*K_c_m[2, 3]*S_Base + μ_24[t]*ηₘ_1[7,t]*K_c_m[2, 7]*S_Base + μ_24[t]*ηₘ_1[10,t]*K_c_m[2, 10]*S_Base 

                                            + μ_23[t]*ηₘ_2[4,t]*K_c_m[1, 19]*S_Base + μ_24[t]*ηₘ_2[4,t]*K_c_m[2, 19]*S_Base 
    
VS_revenue_SGs[5,t] = + μ_23[t]*K_c_gc[1, 5]*S_Base *yˢᴳ²⁷[t] + μ_24[t]*K_c_gc[2, 5]*S_Base *yˢᴳ²⁷[t]

                                            + μ_23[t]*K_c_m[1, 15].*ηₘ_1[15,t]*S_Base + μ_24[t]*K_c_m[2, 15].*ηₘ_1[15,t]*S_Base
                                            + μ_23[t]*ηₘ_1[4,t]*K_c_m[1, 4]*S_Base + μ_23[t]*ηₘ_1[8,t]*K_c_m[1, 8]*S_Base + μ_23[t]*ηₘ_1[11,t]*K_c_m[1, 11]*S_Base + μ_23[t]*ηₘ_1[13,t]*K_c_m[1, 13]*S_Base 
                                            + μ_24[t]*ηₘ_1[4,t]*K_c_m[2, 4]*S_Base + μ_24[t]*ηₘ_1[8,t]*K_c_m[2, 8]*S_Base + μ_24[t]*ηₘ_1[11,t]*K_c_m[2, 11]*S_Base + μ_24[t]*ηₘ_1[13,t]*K_c_m[2, 13]*S_Base 

                                            + μ_23[t]*ηₘ_2[5,t]*K_c_m[1, 20]*S_Base + μ_24[t]*ηₘ_2[5,t]*K_c_m[2, 20]*S_Base  

VS_revenue_SGs[6,t] = μ_23[t]*K_c_gc[1, 6]*S_Base *yˢᴳ³⁰[t] + μ_24[t]*K_c_gc[2, 6]*S_Base *yˢᴳ³⁰[t]

                                            + μ_23[t]*ηₘ_1[5,t]*K_c_m[1, 5]*S_Base + μ_24[t]*ηₘ_1[5,t]*K_c_m[2, 5]*S_Base
                                            + μ_23[t]*ηₘ_1[9,t]*K_c_m[1, 9]*S_Base + μ_23[t]*ηₘ_1[12,t]*K_c_m[1, 12]*S_Base + μ_23[t]*ηₘ_1[14,t]*K_c_m[1, 14]*S_Base 
                                            + μ_24[t]*ηₘ_1[9,t]*K_c_m[2, 9]*S_Base + μ_24[t]*ηₘ_1[12,t]*K_c_m[2, 12]*S_Base + μ_24[t]*ηₘ_1[14,t]*K_c_m[2, 14]*S_Base 
                                            + μ_23[t]*ηₘ_1[15,t]*K_c_m[1, 15]*S_Base + μ_24[t]*ηₘ_1[15,t]*K_c_m[2, 15]*S_Base

                                            + μ_23[t]*ηₘ_2[6,t]*K_c_m[1, 21]*S_Base + μ_24[t]*ηₘ_2[6,t]*K_c_m[2, 21]*S_Base  

end


# revenue for VSGs for providing SCR support

VS_revenue_VSG_1=zeros(T)

for t in 1:T
    VS_revenue_VSG_1[t] = μ_23[t]* ( α_VSG[t]*K_c_gv[1, 1] + sum(K_c_m[1, 16:end] .*ηₘ_2[:,t]   )) *S_Base
                        + μ_24[t]* ( α_VSG[t]*K_c_gv[2, 1] + sum(K_c_m[2, 16:end] .*ηₘ_2[:,t]   )) *S_Base
end



return energy_profit_VSG_1, energy_profitˢᴳ², energy_profitˢᴳ³, energy_profitˢᴳ⁴, energy_profitˢᴳ⁵, energy_profitˢᴳ²⁷, energy_profitˢᴳ³⁰, energy_profit_IBG_23, energy_profit_IBG_24,
       revenue_IBG_23_give_reactive, revenue_IBG_24_give_reactive,
       VS_revenue_VSG_1, VS_revenue_SGs


end