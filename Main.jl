

import Pkg
using JuMP,Gurobi, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles, Plots, Dualization,Printf,Dates
include("dataset_gene.jl")
include("admittance_matrix_calculation.jl") 
include("offline_trainning.jl")
include("reve_pro_calculation_restricted_approx.jl")
include("reve_pro_calculation_dispatchable.jl")
include("pricing_VS_PD_fixprimal_multi_period.jl")
include("pricing_VS_restricted.jl")
include("pricing_VS_dispatchable.jl")

# SGs_ buses:2_3_4_5_27_30    VSG_ buses:1_     IBGs_ buses:23_24

data_bus = CSV.read("IEEE30_Bus_Data.csv", DataFrame)
Bus_Pd = CSV.read("Bus_Pd_24h.csv", DataFrame)
Bus_Qd = CSV.read("Bus_Qd_24h.csv", DataFrame)



#-----------------------------------Representation and Approximation of VOLTAGE STABILITY Constraints  &  Visualization-----------------------------------

t = @elapsed begin
nᵥ= []              # The interval is evenly divided into nᵥ parts
for i in 0:0.001:1
    push!(nᵥ,i) 
end 
IBG=[23,24]         # GFL.   The location (Bus) of IBGs.
Gᵥ=[1]              # FM.    The location (Bus) of Virtual Synchronous Generators (VSGs).
zᵟ_c, matrix_ω =dataset_gene(nᵥ,Gᵥ,IBG)                                                             # data set generation      +sum(K_c_m[2_:] .*matrix_ω[i_14:91])                
K_c_gc, K_c_gv, K_c_m, MAPE_z_23_1, MAPE_z_24_1 = offline_trainning(zᵟ_c, matrix_ω)  # offline_trainning
end
println("Training time: ", t, " s")


#-----------------------------------Define Parameters for Optimization-----------------------------------  
T=24                                    # number of periods

Pᴰ=  sum( Matrix(Bus_Pd[:,(2:25)]  )  , dims=1)  *1.1         # system demand of active power at each bus
Qᴰ=  sum( Matrix(Bus_Qd[:,(2:25)]  )  , dims=1)   *0.7          # system demand of reactive power at each bus

Pˢᴳₘₐₓ= [ 1317 1152 756 667 650 576] /25
Pˢᴳₘᵢₙ= [  658 576 302 133 130 58 ] /25  

Qˢᴳₘₐₓ= Pˢᴳₘₐₓ*0.6                          # REACtive max generation of SGs,  buses:2,3,4,5,27,30
Qˢᴳₘᵢₙ= -Pˢᴳₘₐₓ*0.6                      # REACtive min generation of SGs,  buses:2,3,4,5,27,30

Pⱽˢᴳₘₐₓ=[100]                           # ACtive max generation of VSGs, buses:1
Pⱽˢᴳₘᵢₙ=[0]                             # ACtive min generation of VSGs,  buses:1
Qⱽˢᴳₘₐₓ= Pⱽˢᴳₘₐₓ *0.6                         # REACtive max generation of SGs,  buses:2,3,4,5,27,30
Qⱽˢᴳₘᵢₙ= -Pⱽˢᴳₘₐₓ*0.6  

Pᴵᴮᴳₘₐₓ= [160, 160]                      # ACtive max generation of IBGs, buses:23,24
Pᴵᴮᴳₘᵢₙ=  [0, 0]                        # ACtive min generation of IBGs, buses:23,24

Qᴵᴮᴳₘₐₓ= Pᴵᴮᴳₘₐₓ*0.6
Qᴵᴮᴳₘᵢₙ= -Pᴵᴮᴳₘₐₓ*0.6

Sᵐᵃˣ_gc=Pˢᴳₘₐₓ*1.2
Sᵐᵃˣ_gv=Pⱽˢᴳₘₐₓ*1.2
Sᵐᵃˣ_c=Pᴵᴮᴳₘₐₓ*1.2
S_Base=90


α_VSG=[0.68,0.72,0.75,0.73,0.69,0.62,0.51,0.38,0.26,0.22,0.24,0.28,0.33,0.39,0.46,0.54,0.61,0.67,0.70,0.69,0.66,0.63,0.60,0.57]
α_IBG= [
    0.62 0.65 0.68 0.66 0.63 0.59 0.55 0.52 0.50 0.49 0.51 0.53 0.55 0.58 0.60 0.63 0.65 0.67 0.66 0.64 0.61 0.58 0.56 0.54
    0.78 0.82 0.75 0.61 0.42 0.35 0.28 0.40 0.55 0.68 0.72 0.65 0.52 0.41 0.33 0.45 0.60 0.74 0.81 0.76 0.64 0.51 0.43 0.38
]
ratio=1

cˢᵗ=[  2000 1250 925 720 550 310 ]  /5            #  Startup cost of SGs                     SGs, buses:2,3,4,5,27,30
cᵐ=[  6.20 6.91 10.47 12.28 13.53 15.36 ]                   #  Marginal generation cost of SGs   in    SGs, buses:2,3,4,5,27,30
cⁿˡ= [ 1743  1501  1376 1093   990    857   ] /25

yˢᴳ²_0 = 1
yˢᴳ³_0 = 1
yˢᴳ⁴_0 = 0
yˢᴳ⁵_0 = 0
yˢᴳ²⁷_0 = 0
yˢᴳ³⁰_0 = 0

LCOE_gv = 12000
LCOE_gf = [13000, 13000]

#-----------------------------------Tests-----------------------------------  

λᴱ_PD_list_against_reactive=zeros(12)
λᴱ_Res_list_against_reactive=zeros(12)
λᴱ_Dis_list_against_reactive=zeros(12)

μ_23_PD_list_against_reactive=zeros(12)
μ_23_Res_list_against_reactive=zeros(12)
μ_23_Dis_list_against_reactive=zeros(12)
μ_24_PD_list_against_reactive=zeros(12)
μ_24_Res_list_against_reactive=zeros(12)
μ_24_Dis_list_against_reactive=zeros(12)

Q_23_PD_list_against_reactive=zeros(12)
Q_23_Res_list_against_reactive=zeros(12)
Q_23_Dis_list_against_reactive=zeros(12)
Q_24_PD_list_against_reactive=zeros(12)
Q_24_Res_list_against_reactive=zeros(12)
Q_24_Dis_list_against_reactive=zeros(12)

profit_PD_list_against_reactive=zeros(9,12)
profit_Res_list_against_reactive=zeros(9,12)
profit_Dis_list_against_reactive=zeros(9,12)


λᴱ_PD_list_against_load=zeros(12)
λᴱ_Res_list_against_load=zeros(12)
λᴱ_Dis_list_against_load=zeros(12)

μ_23_PD_list_against_load=zeros(12)
μ_23_Res_list_against_load=zeros(12)
μ_23_Dis_list_against_load=zeros(12)
μ_24_PD_list_against_load=zeros(12)
μ_24_Res_list_against_load=zeros(12)
μ_24_Dis_list_against_load=zeros(12)
Q_23_PD_list_against_load=zeros(12)
Q_23_Res_list_against_load=zeros(12)
Q_23_Dis_list_against_load=zeros(12)
Q_24_PD_list_against_load=zeros(12)
Q_24_Res_list_against_load=zeros(12)
Q_24_Dis_list_against_load=zeros(12)

uplifts_total_against_reactive=zeros(12)
uplifts_total_against_load=zeros(12)


wind_curtailment_list = zeros(2,12)
SCR_list = zeros(2,12)

ratio_list = [0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]
load_list = [0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1]

profit_PD_list_against_load=zeros(9,12)
profit_Res_list_against_load=zeros(9,12)
profit_Dis_list_against_load=zeros(9,12)

profit_PD_list_against_reactive=zeros(9,12)
profit_Res_list_against_reactive=zeros(9,12)
profit_Dis_list_against_reactive=zeros(9,12)


index=0
for ratio in ratio_list
#for load in load_list
index = index + 1

#Pᴰ= sum( Matrix(Bus_Pd[:,(2:25)]  )  , dims=1)  *load

#------------------------------- Pricing via Primal-Dual -----------------------------------


λᴱ_PD, μ_23_PD, μ_24_PD, Q_23_PD, Q_24_PD, energy_profit_SG_PD, vs_revenue_SG_PD, profit_SG_PD, profit_wts_energy, profit_wts_vs = pricing_VS_PD(K_c_gc, K_c_gv, K_c_m,
                        T,Pᴰ,Qᴰ,Pˢᴳₘₐₓ,Pˢᴳₘᵢₙ,Qˢᴳₘₐₓ,Qˢᴳₘᵢₙ,Pⱽˢᴳₘₐₓ,Pⱽˢᴳₘᵢₙ,Qⱽˢᴳₘₐₓ,Qⱽˢᴳₘᵢₙ,Pᴵᴮᴳₘₐₓ,Pᴵᴮᴳₘᵢₙ,Qᴵᴮᴳₘₐₓ,Qᴵᴮᴳₘᵢₙ,
                        Sᵐᵃˣ_gc,Sᵐᵃˣ_gv,Sᵐᵃˣ_c,S_Base,α_VSG,α_IBG,ratio,
                        cˢᵗ,cᵐ,cⁿˡ,
                        yˢᴳ²_0,yˢᴳ³_0,yˢᴳ⁴_0,yˢᴳ⁵_0,yˢᴳ²⁷_0,yˢᴳ³⁰_0,
                        LCOE_gv, LCOE_gf)

                        

#------------------------------- Pricing via Restricted -----------------------------------

λᴱ_Res, μ_23_Res, μ_24_Res, Q_23_Res, Q_24_Res, pro_SG_no_vs, VS_revenue_SGs_Res, uplifts,
energy_profit_VSG_1_Res, energy_profit_IBG_23_Res, energy_profit_IBG_24_Res, 
VS_revenue_VSG_1_Res, revenue_IBG_23_give_reactive_Res, revenue_IBG_24_give_reactive_Res, energy_profit_SGs_Res,
wind_curtailment, SCR  = pricing_VS_restricted(K_c_gc, K_c_gv, K_c_m,
                        T,Pᴰ,Qᴰ,Pˢᴳₘₐₓ,Pˢᴳₘᵢₙ,Qˢᴳₘₐₓ,Qˢᴳₘᵢₙ,Pⱽˢᴳₘₐₓ,Pⱽˢᴳₘᵢₙ,Qⱽˢᴳₘₐₓ,Qⱽˢᴳₘᵢₙ,Pᴵᴮᴳₘₐₓ,Pᴵᴮᴳₘᵢₙ,Qᴵᴮᴳₘₐₓ,Qᴵᴮᴳₘᵢₙ,
                        Sᵐᵃˣ_gc,Sᵐᵃˣ_gv,Sᵐᵃˣ_c,S_Base,α_VSG,α_IBG,ratio,
                        cˢᵗ,cᵐ,cⁿˡ,
                        yˢᴳ²_0,yˢᴳ³_0,yˢᴳ⁴_0,yˢᴳ⁵_0,yˢᴳ²⁷_0,yˢᴳ³⁰_0)

# wind_curtailment_list[1,index] =   wind_curtailment[1]        
# wind_curtailment_list[2,index] =   wind_curtailment[2]         
# SCR_list[1,index] =   SCR[1]        
# SCR_list[2,index] =   SCR[2]  

# plot(wind_curtailment_list')
# plot(SCR_list')
# println(SCR_list[2,:])


#------------------------------- Pricing via Dispatchable -----------------------------------


λᴱ_Dis, μ_23_Dis, μ_24_Dis, Q_23_Dis, Q_24_Dis, energy_profit_SG_Dis, vs_revenue_SG_Dis, 
energy_profit_VSG_1_Dis, energy_profit_IBG_23_Dis, energy_profit_IBG_24_Dis, 
VS_revenue_VSG_1_Dis, revenue_IBG_23_give_reactive_Dis, revenue_IBG_24_give_reactive_Dis = pricing_VS_dispatchable(K_c_gc, K_c_gv, K_c_m,
                        T,Pᴰ,Qᴰ,Pˢᴳₘₐₓ,Pˢᴳₘᵢₙ,Qˢᴳₘₐₓ,Qˢᴳₘᵢₙ,Pⱽˢᴳₘₐₓ,Pⱽˢᴳₘᵢₙ,Qⱽˢᴳₘₐₓ,Qⱽˢᴳₘᵢₙ,Pᴵᴮᴳₘₐₓ,Pᴵᴮᴳₘᵢₙ,Qᴵᴮᴳₘₐₓ,Qᴵᴮᴳₘᵢₙ,
                        Sᵐᵃˣ_gc,Sᵐᵃˣ_gv,Sᵐᵃˣ_c,S_Base,α_VSG,α_IBG,ratio,
                        cˢᵗ,cᵐ,cⁿˡ,
                        yˢᴳ²_0,yˢᴳ³_0,yˢᴳ⁴_0,yˢᴳ⁵_0,yˢᴳ²⁷_0,yˢᴳ³⁰_0)

#=
λᴱ_PD_list_against_reactive[index]=sum(λᴱ_PD)/T
λᴱ_Res_list_against_reactive[index]=sum(λᴱ_Res)/T
λᴱ_Dis_list_against_reactive[index]=sum(λᴱ_Dis)/T
μ_23_PD_list_against_reactive[index]=sum(μ_23_PD)/T
μ_23_Res_list_against_reactive[index]=sum(μ_23_Res)/T
μ_23_Dis_list_against_reactive[index]=sum(μ_23_Dis)/T
μ_24_PD_list_against_reactive[index]=sum(μ_24_PD)/T
μ_24_Res_list_against_reactive[index]=sum(μ_24_Res)/T
μ_24_Dis_list_against_reactive[index]=sum(μ_24_Dis)/T
Q_23_PD_list_against_reactive[index]=sum(Q_23_PD)/T
Q_23_Res_list_against_reactive[index]=sum(Q_23_Res)/T
Q_23_Dis_list_against_reactive[index]=sum(Q_23_Dis)/T
Q_24_PD_list_against_reactive[index]=sum(Q_24_PD)/T
Q_24_Res_list_against_reactive[index]=sum(Q_24_Res)/T
Q_24_Dis_list_against_reactive[index]=sum(Q_24_Dis)/T

uplifts_total_against_reactive[index]=sum(uplifts)



λᴱ_PD_list_against_load[index]=sum(λᴱ_PD)/T
λᴱ_Res_list_against_load[index]=sum(λᴱ_Res)/T
λᴱ_Dis_list_against_load[index]=sum(λᴱ_Dis)/T
μ_23_PD_list_against_load[index]=sum(μ_23_PD)/T
μ_23_Res_list_against_load[index]=sum(μ_23_Res)/T
μ_23_Dis_list_against_load[index]=sum(μ_23_Dis)/T
μ_24_PD_list_against_load[index]=sum(μ_24_PD)/T
μ_24_Res_list_against_load[index]=sum(μ_24_Res)/T
μ_24_Dis_list_against_load[index]=sum(μ_24_Dis)/T
Q_23_PD_list_against_load[index]=sum(Q_23_PD)/T
Q_23_Res_list_against_load[index]=sum(Q_23_Res)/T
Q_23_Dis_list_against_load[index]=sum(Q_23_Dis)/T
Q_24_PD_list_against_load[index]=sum(Q_24_PD)/T
Q_24_Res_list_against_load[index]=sum(Q_24_Res)/T
Q_24_Dis_list_against_load[index]=sum(Q_24_Dis)/T

uplifts_total_against_load[index]=sum(uplifts)
=#

for i in 1:6
#=
    profit_PD_list_against_load[i,index] = sum(profit_SG_PD[i,:])
    profit_PD_list_against_load[7,index] = sum(profit_wts_energy[1,:]) + sum(profit_wts_vs[1,:]) - LCOE_gv
    profit_PD_list_against_load[8,index] = sum(profit_wts_energy[2,:]) + sum(profit_wts_vs[2,:]) - LCOE_gf[1]
    profit_PD_list_against_load[9,index] = sum(profit_wts_energy[3,:]) + sum(profit_wts_vs[3,:]) - LCOE_gf[2]

    profit_Dis_list_against_load[i,index] = sum(energy_profit_SG_Dis[i,:]) + sum(vs_revenue_SG_Dis[i,:])
    profit_Dis_list_against_load[7,index] = sum(energy_profit_VSG_1_Dis) + sum(VS_revenue_VSG_1_Dis) - LCOE_gv
    profit_Dis_list_against_load[8,index] = sum(energy_profit_IBG_23_Dis) + sum(revenue_IBG_23_give_reactive_Dis) - LCOE_gf[1]
    profit_Dis_list_against_load[9,index] = sum(energy_profit_IBG_24_Dis) + sum(revenue_IBG_24_give_reactive_Dis) - LCOE_gf[2]


    profit_Res_list_against_load[i,index] = pro_SG_no_vs[i] + sum(VS_revenue_SGs_Res[i,:]) 
    profit_Res_list_against_load[7,index] = sum(energy_profit_VSG_1_Res) + sum(VS_revenue_VSG_1_Res) - LCOE_gv
    profit_Res_list_against_load[8,index] = sum(energy_profit_IBG_23_Res) + sum(revenue_IBG_23_give_reactive_Res) - LCOE_gf[1]
    profit_Res_list_against_load[9,index] = sum(energy_profit_IBG_24_Res) + sum(revenue_IBG_24_give_reactive_Res) - LCOE_gf[2]
=#


    profit_PD_list_against_reactive[i,index] = sum(profit_SG_PD[i,:])
    profit_PD_list_against_reactive[7,index] = sum(profit_wts_energy[1,:]) + sum(profit_wts_vs[1,:]) - LCOE_gv
    profit_PD_list_against_reactive[8,index] = sum(profit_wts_energy[2,:]) + sum(profit_wts_vs[2,:]) - LCOE_gf[1]
    profit_PD_list_against_reactive[9,index] = sum(profit_wts_energy[3,:]) + sum(profit_wts_vs[3,:]) - LCOE_gf[2]

    profit_Dis_list_against_reactive[i,index] = sum(energy_profit_SG_Dis[i,:]) + sum(vs_revenue_SG_Dis[i,:])
    profit_Dis_list_against_reactive[7,index] = sum(energy_profit_VSG_1_Dis) + sum(VS_revenue_VSG_1_Dis) - LCOE_gv
    profit_Dis_list_against_reactive[8,index] = sum(energy_profit_IBG_23_Dis) + sum(revenue_IBG_23_give_reactive_Dis) - LCOE_gf[1]
    profit_Dis_list_against_reactive[9,index] = sum(energy_profit_IBG_24_Dis) + sum(revenue_IBG_24_give_reactive_Dis) - LCOE_gf[2]


    profit_Res_list_against_reactive[i,index] = pro_SG_no_vs[i] + sum(VS_revenue_SGs_Res[i,:]) 
    profit_Res_list_against_reactive[7,index] = sum(energy_profit_VSG_1_Res) + sum(VS_revenue_VSG_1_Res) - LCOE_gv
    profit_Res_list_against_reactive[8,index] = sum(energy_profit_IBG_23_Res) + sum(revenue_IBG_23_give_reactive_Res) - LCOE_gf[1]
    profit_Res_list_against_reactive[9,index] = sum(energy_profit_IBG_24_Res) + sum(revenue_IBG_24_give_reactive_Res) - LCOE_gf[2]


end

end

plot(profit_PD_list_against_load[4,:])
plot!(profit_Dis_list_against_load[4,:])
plot!(profit_Res_list_against_load[4,:])
plot!(profit_PD_list_against_load[8,:])
plot!(profit_Dis_list_against_load[8,:])
plot!(profit_Res_list_against_load[8,:])


plot(profit_PD_list_against_reactive[4,:])
plot!(profit_Dis_list_against_reactive[4,:])
plot!(profit_Res_list_against_reactive[4,:])
plot(profit_PD_list_against_reactive[7,:])
plot!(profit_Dis_list_against_reactive[7,:])
plot!(profit_Res_list_against_reactive[7,:])

println(profit_Dis_list_against_reactive[7,:])
println(profit_Res_list_against_reactive[7,:])
println(profit_PD_list_against_reactive[7,:])

#------------------------------- Analyze results -----------------------------------

sum(energy_profit_VSG_1_Res)
sum(energy_profit_IBG_23_Res)
sum(energy_profit_IBG_24_Res)
sum(VS_revenue_VSG_1_Res)
sum(revenue_IBG_23_give_reactive_Res)
sum(revenue_IBG_24_give_reactive_Res)

sum(energy_profit_VSG_1_Dis)
sum(energy_profit_IBG_23_Dis)
sum(energy_profit_IBG_24_Dis)
sum(VS_revenue_VSG_1_Dis)
sum(revenue_IBG_23_give_reactive_Dis)
sum(revenue_IBG_24_give_reactive_Dis)

sum(energy_profit_IBG_24_Dis) - 13000 + sum(revenue_IBG_24_give_reactive_Dis) - 13000





plot(Q_23_PD_list_against_reactive)
plot!(Q_23_Res_list_against_reactive)
plot!(Q_23_Dis_list_against_reactive)

plot(Q_24_PD_list_against_reactive)
plot!(Q_24_Res_list_against_reactive)
plot!(Q_24_Dis_list_against_reactive)

plot(λᴱ_PD_list_against_reactive)
plot!(λᴱ_Res_list_against_reactive)
plot!(λᴱ_Dis_list_against_reactive)
plot!(uplifts_total_against_reactive/400)





plot(Q_23_Dis_list_against_load)
plot!(Q_23_Res_list_against_load)
plot!(Q_23_PD_list_against_load)

plot(Q_24_Dis_list_against_load)
plot!(Q_24_Res_list_against_load)
plot!(Q_24_PD_list_against_load)

plot(λᴱ_Dis_list_against_load)
plot!(λᴱ_Res_list_against_load)
plot!(λᴱ_PD_list_against_load)
plot!(uplifts_total_against_load/400)

println(uplifts_total_against_load/400)



plot(λᴱ_PD)
plot!(λᴱ_Dis)
plot!(λᴱ_Res)

plot(μ_23_Dis)
plot!(μ_23_Res)
plot!(μ_23_PD)
plot(μ_24_Dis)
plot!(μ_24_Res)
plot!(μ_24_PD)

plot(Q_23_Dis)
plot!(Q_23_Res)
plot!(Q_23_PD)

plot(Q_24_Dis)
plot!(Q_24_Res)
plot!(Q_24_PD)

sum(μ_23_Dis)/T
sum(μ_23_Res)/T
sum(μ_23_PD)/T
sum(μ_24_Dis)/T
sum(μ_24_Res)/T
sum(μ_24_PD)/T

sum(Q_23_Dis)/T
sum(Q_23_Res)/T
sum(Q_23_PD)/T
sum(Q_24_Dis)/T
sum(Q_24_Res)/T
sum(Q_24_PD)/T

sum(λᴱ_PD)/T
sum(λᴱ_Res)/T
sum(λᴱ_Dis)/T

final_profit_SG_Dis = energy_profit_SG_Dis + vs_revenue_SG_Dis

sum(energy_profit_SG_Dis[1,:])
sum(energy_profit_SG_Dis[2,:])
sum(energy_profit_SG_Dis[3,:])
sum(energy_profit_SG_Dis[4,:])
sum(energy_profit_SG_Dis[5,:])
sum(energy_profit_SG_Dis[6,:])

sum(vs_revenue_SG_Dis[1,:])
sum(vs_revenue_SG_Dis[2,:])
sum(vs_revenue_SG_Dis[3,:])
sum(vs_revenue_SG_Dis[4,:])
sum(vs_revenue_SG_Dis[5,:])
sum(vs_revenue_SG_Dis[6,:])

sum(final_profit_SG_Dis[1,:])
sum(final_profit_SG_Dis[2,:])
sum(final_profit_SG_Dis[3,:])
sum(final_profit_SG_Dis[4,:])
sum(final_profit_SG_Dis[5,:])
sum(final_profit_SG_Dis[6,:])

final_profit_SG_PD = energy_profit_SG_PD + vs_revenue_SG_PD

sum(energy_profit_SG_PD[1,:])
sum(energy_profit_SG_PD[2,:])
sum(energy_profit_SG_PD[3,:])
sum(energy_profit_SG_PD[4,:])
sum(energy_profit_SG_PD[5,:])
sum(energy_profit_SG_PD[6,:])

sum(vs_revenue_SG_PD[1,:])
sum(vs_revenue_SG_PD[2,:])
sum(vs_revenue_SG_PD[3,:])
sum(vs_revenue_SG_PD[4,:])
sum(vs_revenue_SG_PD[5,:])
sum(vs_revenue_SG_PD[6,:])

sum(final_profit_SG_PD[1,:])
sum(final_profit_SG_PD[2,:])
sum(final_profit_SG_PD[3,:])
sum(final_profit_SG_PD[4,:])
sum(final_profit_SG_PD[5,:])
sum(final_profit_SG_PD[6,:])

sum(profit_wts_energy[1,:]) - 12000 + sum(profit_wts_vs[1,:]) 
sum(profit_wts_energy[2,:]) - 13000 + sum(profit_wts_vs[2,:]) 
sum(profit_wts_energy[3,:]) - 13000 + sum(profit_wts_vs[3,:]) 

sum(profit_wts_vs[1,:]) 
sum(profit_wts_vs[2,:]) 
sum(profit_wts_vs[3,:]) 



final_profit_SG_Res = pro_SG_no_vs + sum(VS_revenue_SGs_Res[:,i] for i in 1:T)

sum(uplifts[1,:])
sum(uplifts[2,:])
sum(uplifts[3,:])
sum(uplifts[4,:])
sum(uplifts[5,:])
sum(uplifts[6,:])

sum(energy_profit_VSG_1_Res) - 12000
sum(energy_profit_IBG_23_Res) - 13000
sum(energy_profit_IBG_24_Res) - 13000

energy_profit_SGs_Res