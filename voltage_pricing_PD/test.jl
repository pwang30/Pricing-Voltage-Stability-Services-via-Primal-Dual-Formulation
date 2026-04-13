# Install required packages (only for the first run)
import Pkg DataStructures
using Pkg
Pkg.add(["JuMP", "Gurobi"])

# Load optimization libraries
using JuMP
using Gurobi
using DataStructures
using Dualization

# Basic system parameters
T = 24                  # Number of time periods (24 hours)
G = 2                   # Number of thermal generation units

# Generation cost coefficients: Cost = a*P + b*P^2 (quadratic cost function)
a = [100.0, 150.0]      # Linear cost coefficients
b = [2.0, 1.5]          # Quadratic cost coefficients

# Generator operating limits (MW)
P_min = [20.0, 30.0]    # Minimum power output
P_max = [100.0, 150.0]  # Maximum power output

# 24-hour load demand profile (MW)
load = [
    80, 75, 70, 65, 60, 65,    # Hour 0 - 6
    90, 120, 150, 140, 130,    # Hour 7 - 11
    160, 170, 165, 155, 150,   # Hour 12 - 16
    180, 200, 190, 170, 150,   # Hour 17 - 21
    120, 90, 85                   # Hour 22 - 23
]




# Build optimization model
model = Model(Gurobi.Optimizer)  # Initialize model with Gurobi solver


gen_id = [1,2]

# Initialize dictionary to store results
ED_dict = OrderedDict{Symbol, Any}()
ED_dict[:variables] = OrderedDict{Symbol, OrderedDict{Int, JuMP.VariableRef}}()

P_g = OrderedDict{Int, JuMP.VariableRef}()

for (idx, gen) in enumerate(gen_id)
    P_g[gen] = @variable(model, base_name = "P_g[$gen]")
end

ED_dict[:variables][:P_g] = P_g


# Decision variable: power output of generator g at hour t (MW)
@variable(model, P_min[g] ≤ P[g=1:G, t=1:T] ≤ P_max[g])

# Objective: minimize total 24-hour generation cost
@objective(model, Min, sum(sum(a[g] * P[g,t] + b[g] * P[g,t]^2 for g in 1:G) for t in 1:T))

# Power balance constraint: total generation = load demand at each hour
@constraint(model, PowerBalance[t=1:T], sum(P[g,t] for g in 1:G) == load[t])

# Solve the optimization problem
optimize!(model)

# Print results
println("===== 24-Hour Economic Dispatch Results =====")
println("Total Minimum Generation Cost: ", round(objective_value(model), digits=2), " USD")
println("\nHourly Generation Output (MW):")
for t in 1:T
    p1 = round(value(P[1,t]), digits=2)
    p2 = round(value(P[2,t]), digits=2)
    println("Hour $t | Unit 1: $p1 MW | Unit 2: $p2 MW | Load: $(load[t]) MW")
end

# Print results
dual_model = dualize(model)
print(dual_model)



#======================================================================#
#======================================================================#
#======================================================================#

using JuMP, HiGHS, DataStructures

# ==============================
# Parameters
# ==============================
T = 24                     # Number of time periods (hours)
gen_id = [1, 2]            # Generator IDs
G = length(gen_id)         # Number of generators

a = [100.0, 150.0]          # Linear cost coefficients
b = [2.0, 1.5]              # Quadratic cost coefficients
P_min = [20.0, 30.0]        # Minimum generation limit
P_max = [100.0, 150.0]      # Maximum generation limit

# 24-hour load profile
load = [80,75,70,65,60,65,90,120,150,140,130,160,170,165,155,150,180,200,190,170,150,120,90,85]

# ==============================
# Initialize model
# ==============================
model = Model(Gurobi.Optimizer)
set_silent(model)

# ==============================
# MAIN DICTIONARY FOR ECONOMIC DISPATCH
# ==============================
ED_dict = OrderedDict{Symbol, Any}()
ED_dict[:variables] = OrderedDict{Symbol, OrderedDict{Tuple{Int,Int}, JuMP.VariableRef}}()
ED_dict[:eq_constraints] = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()
ED_dict[:ineq_constraints] = OrderedDict{Symbol, OrderedDict{Tuple{Int,Int}, JuMP.ConstraintRef}}()

# ==============================
# 1. Define and store decision variables: P[g, t]
# ==============================
P_g = OrderedDict{Tuple{Int,Int}, JuMP.VariableRef}()  # key: (g, t)

for t in 1:T
    for (idx, gen) in enumerate(gen_id)
        P_g[(gen, t)] = @variable(
            model,
            base_name = "P_g[$gen,$t]",
            lower_bound = P_min[idx],
            upper_bound = P_max[idx]
        )
    end
end

ED_dict[:variables][:P_g] = P_g

# ==============================
# 2. Objective function
# ==============================
@objective(model, Min, sum(
    a[findfirst(==(gen), gen_id)] * P_g[(gen,t)] +
    b[findfirst(==(gen), gen_id)] * P_g[(gen,t)]^2
    for gen in gen_id, t in 1:T
))

# ==============================
# 3. Power balance (EQUALITY constraints)
# ==============================
PowerBalance = OrderedDict{Int, JuMP.ConstraintRef}()

for t in 1:T
    PowerBalance[t] = @constraint(
        model,
        sum(P_g[(gen,t)] for gen in gen_id) == load[t]
    )
end

ED_dict[:eq_constraints][:PowerBalance] = PowerBalance

# ==============================
# 4. Generator lower/upper bounds (already in variable; store as inequality)
# ==============================
GenLowerLimit = OrderedDict{Tuple{Int,Int}, JuMP.ConstraintRef}()
GenUpperLimit = OrderedDict{Tuple{Int,Int}, JuMP.ConstraintRef}()

for t in 1:T
    for (idx, gen) in enumerate(gen_id)
        GenLowerLimit[(gen,t)] = @constraint(model, P_g[(gen,t)] >= P_min[idx])
        GenUpperLimit[(gen,t)] = @constraint(model, P_g[(gen,t)] <= P_max[idx])
    end
end

ED_dict[:ineq_constraints][:GenLowerLimit] = GenLowerLimit
ED_dict[:ineq_constraints][:GenUpperLimit] = GenUpperLimit

# ==============================
# Solve
# ==============================
optimize!(model)

# ==============================
# Access example
# ==============================
println("Total cost: ", objective_value(model))
println("Generator 1 at hour 5: ", value(ED_dict[:variables][:P_g][(1,5)]))



println(dualize(model))


using Printf

open("Economic_Dispatch_Result.txt", "w") do f
    @printf(f, "===============================\n")
    @printf(f, " 24-HOUR ECONOMIC DISPATCH RESULT\n")
    @printf(f, "===============================\n\n")
    
    @printf(f, "Total Minimum Cost: %.2f USD\n\n", objective_value(model))
    @printf(f, "Hour\tGen1(MW)\tGen2(MW)\tLoad(MW)\n")
    
    for t in 1:24
        p1 = value(ED_dict[:variables][:P_g][(1, t)])
        p2 = value(ED_dict[:variables][:P_g][(2, t)])
        @printf(f, "%d\t%.2f\t\t%.2f\t\t%d\n", t, p1, p2, load[t])
    end
end

println("\n✅ Model & results saved to: Economic_Dispatch_Result.txt")










# ----------------------------------------------------
# Save and print the MODEL STRUCTURE (no results)
# Follows exactly your formatting style
# ----------------------------------------------------

# Open the file for writing
open("model_summary.txt", "w") do io
    # Print the model to the file
    show(io, model)
end

# Desired key order to print the variables
vector_dict_var = [ED_dict[:variables][:P_g]]

open("model_details.txt", "w") do io

    # ------------------
    # Objective function
    # ------------------
    println(io, "=========")
    println(io, "Objective ")
    println(io, "=========")
    println(io, objective_function(model))
    println(io, "\n")

    # ---------------------------
    # Variables used in the model
    # ---------------------------
    println(io, "=========")
    println(io, "Variables")
    println(io, "=========")
    for i in eachindex(vector_dict_var)
        for (j, info) in vector_dict_var[i]
            println(io, "$j: ", info)
        end
    end
    println(io, "\n")

    # --------------------
    # Equality constraint
    # --------------------
    if haskey(ED_dict[:eq_constraints], :PowerBalance)
        println(io, "===================================================")
        println(io, "Equality Constraints - Active Power Balance ")
        println(io, "===================================================")
        for (i, info) in ED_dict[:eq_constraints][:PowerBalance]
            println(io, "$i: ", info) 
        end
        println(io, "\n")
    end

    # =========================================
    # Decision variables lower and upper bounds
    # =========================================

    println(io, "=========================================================")
    println(io, "Inequality Constraints - Generator Output Limits ")
    println(io, "=========================================================")

    if haskey(ED_dict[:ineq_constraints], :GenLowerLimit)
        println(io, "====================================")
        println(io, "Active Power Generated - Lower Bound ")
        println(io, "====================================")
        for (i, info) in ED_dict[:ineq_constraints][:GenLowerLimit]
            println(io, "$i: ", info)
        end
        println(io, "\n")
    end

    if haskey(ED_dict[:ineq_constraints], :GenUpperLimit)
        println(io, "====================================")
        println(io, "Active Power Generated - Upper Bound ")
        println(io, "====================================")
        for (i, info) in ED_dict[:ineq_constraints][:GenUpperLimit]
            println(io, "$i: ", info)
        end
        println(io, "\n")
    end

end

println("✅ Model structure saved to model_summary.txt and model_details.txt")


dual_model = Model(Gurobi.Optimizer)

dual_model = dualize(model)
optimize!(dual_model)


dual_model = Model(dual_optimizer(Gurobi.Optimizer))


dual_model = dualize(model)  # Automatically create dual problem
set_optimizer(dual_model, Gurobi.Optimizer)
optimize!(dual_model)        # Solve the dual model






# =================================================================
# Print and save the DUAL MODEL structure (objective + constraints)
# Format is fully aligned with your original model printing code
# =================================================================

# Save full dual model to file
open("dual_model_summary.txt", "w") do io
    show(io, dual_model)
end

# Save detailed dual model (objective, variables, constraints)
open("dual_model_details.txt", "w") do io

    # --------------------------
    # Dual Objective Function
    # --------------------------
    println(io, "========================")
    println(io, "DUAL OBJECTIVE FUNCTION")
    println(io, "========================")
    println(io, objective_function(dual_model))
    println(io, "\n")

    # --------------------------
    # Dual Variables
    # --------------------------
    println(io, "========================")
    println(io, "DUAL VARIABLES")
    println(io, "========================")
    for var in all_variables(dual_model)
        println(io, var)
    end
    println(io, "\n")

    # --------------------------
    # Dual Constraints
    # --------------------------
    println(io, "========================")
    println(io, "DUAL CONSTRAINTS")
    println(io, "========================")

    # Print all constraint types in the dual model
    for (F, S) in list_of_constraint_types(dual_model)
        println(io, "\n[Constraint Type: $F in $S]")
        cons = all_constraints(dual_model, F, S)
        for con in cons
            println(io, "  ", con)
        end
    end

end

println("✅ Dual model saved to: dual_model_summary.txt + dual_model_details.txt")












# ================================
# Model
# ================================
model = Model()

# ================================
# Variables
# ================================
@variable(model, yˢᴳ²[1:T], Bin)
@variable(model, yˢᴳ³[1:T], Bin)
@variable(model, yˢᴳ⁴[1:T], Bin)
@variable(model, yˢᴳ⁵[1:T], Bin)
@variable(model, yˢᴳ²⁷[1:T], Bin)
@variable(model, yˢᴳ³⁰[1:T], Bin)

@variable(model, Cᵁ²[1:T] >= 0)
@variable(model, Cᵁ³[1:T] >= 0)
@variable(model, Cᵁ⁴[1:T] >= 0)
@variable(model, Cᵁ⁵[1:T] >= 0)
@variable(model, Cᵁ²⁷[1:T] >= 0)
@variable(model, Cᵁ³⁰[1:T] >= 0)

@variable(model, Pᴳ[1:30, 1:T])
@variable(model, Qᴳ[1:30, 1:T])

@variable(model, Pˡⁱⁿᵉ[1:30,1:30,1:T])
@variable(model, Qˡⁱⁿᵉ[1:30,1:30,1:T])

@variable(model, Γ_23[1:T])
@variable(model, Γ_24[1:T])
@variable(model, ηₘ_1[1:15,1:T], Bin)
@variable(model, ηₘ_2[1:6,1:T])

# ======================================================================
# ✅ ALL CONSTRAINTS ARE NAMED AND STORED IN DICTIONARIES
# ======================================================================

# --------------------------
# Dictionaries for constraints
# --------------------------
ineq_startup = OrderedDict{Tuple{Int,Int}, ConstraintRef}()

ineq_P_up = OrderedDict{Tuple{Int,Int}, ConstraintRef}()
ineq_P_lo = OrderedDict{Tuple{Int,Int}, ConstraintRef}()
ineq_Q_up = OrderedDict{Tuple{Int,Int}, ConstraintRef}()
ineq_Q_lo = OrderedDict{Tuple{Int,Int}, ConstraintRef}()

soc_gen = OrderedDict{Tuple{Int,Int}, ConstraintRef}()

eq_p_bal = OrderedDict{Tuple{Int,Int}, ConstraintRef}()
eq_q_bal = OrderedDict{Tuple{Int,Int}, ConstraintRef}()
eq_anti_p = OrderedDict{Tuple{Int,Int,Int}, ConstraintRef}()
eq_anti_q = OrderedDict{Tuple{Int,Int,Int}, ConstraintRef}()
eq_zero_p = OrderedDict{Tuple{Int,Int}, ConstraintRef}()
eq_zero_q = OrderedDict{Tuple{Int,Int}, ConstraintRef}()

ineq_eta = OrderedDict{Tuple{Int,Int}, ConstraintRef}()
eq_eta2 = OrderedDict{Tuple{Int,Int}, ConstraintRef}()
eq_gamma23 = OrderedDict{Int, ConstraintRef}()
eq_gamma24 = OrderedDict{Int, ConstraintRef}()
soc_vs = OrderedDict{Int, ConstraintRef}()

# --------------------------
# 1. Startup cost constraints
# --------------------------
for t in 2:T
    ineq_startup[(2,t)]  = @constraint(model, Cᵁ²[t] >= (yˢᴳ²[t]-yˢᴳ²[t-1])*Kˢᵗ[1])
    ineq_startup[(3,t)]  = @constraint(model, Cᵁ³[t] >= (yˢᴳ³[t]-yˢᴳ³[t-1])*Kˢᵗ[2])
    ineq_startup[(4,t)]  = @constraint(model, Cᵁ⁴[t] >= (yˢᴳ⁴[t]-yˢᴳ⁴[t-1])*Kˢᵗ[3])
    ineq_startup[(5,t)]  = @constraint(model, Cᵁ⁵[t] >= (yˢᴳ⁵[t]-yˢᴳ⁵[t-1])*Kˢᵗ[4])
    ineq_startup[(27,t)] = @constraint(model, Cᵁ²⁷[t] >= (yˢᴳ²⁷[t]-yˢᴳ²⁷[t-1])*Kˢᵗ[5])
    ineq_startup[(30,t)] = @constraint(model, Cᵁ³⁰[t] >= (yˢᴳ³⁰[t]-yˢᴳ³⁰[t-1])*Kˢᵗ[6])
end

# --------------------------
# 2. Generator PQ limits
# --------------------------
for t in 1:T
    ineq_P_up[(2,t)]  = @constraint(model, Pᴳ[2,t] <= yˢᴳ²[t]*Pˢᴳₘₐₓ[1])
    ineq_P_lo[(2,t)]  = @constraint(model, yˢᴳ²[t]*Pˢᴳₘᵢₙ[1] <= Pᴳ[2,t])
    ineq_P_up[(3,t)]  = @constraint(model, Pᴳ[3,t] <= yˢᴳ³[t]*Pˢᴳₘₐₓ[2])
    ineq_P_lo[(3,t)]  = @constraint(model, yˢᴳ³[t]*Pˢᴳₘᵢₙ[2] <= Pᴳ[3,t])
    ineq_P_up[(4,t)]  = @constraint(model, Pᴳ[4,t] <= yˢᴳ⁴[t]*Pˢᴳₘₐₓ[3])
    ineq_P_lo[(4,t)]  = @constraint(model, yˢᴳ⁴[t]*Pˢᴳₘᵢₙ[3] <= Pᴳ[4,t])
    ineq_P_up[(5,t)]  = @constraint(model, Pᴳ[5,t] <= yˢᴳ⁵[t]*Pˢᴳₘₐₓ[4])
    ineq_P_lo[(5,t)]  = @constraint(model, yˢᴳ⁵[t]*Pˢᴳₘᵢₙ[4] <= Pᴳ[5,t])
    ineq_P_up[(27,t)] = @constraint(model, Pᴳ[27,t] <= yˢᴳ²⁷[t]*Pˢᴳₘₐₓ[5])
    ineq_P_lo[(27,t)] = @constraint(model, yˢᴳ²⁷[t]*Pˢᴳₘᵢₙ[5] <= Pᴳ[27,t])
    ineq_P_up[(30,t)] = @constraint(model, Pᴳ[30,t] <= yˢᴳ³⁰[t]*Pˢᴳₘₐₓ[6])
    ineq_P_lo[(30,t)] = @constraint(model, yˢᴳ³⁰[t]*Pˢᴳₘᵢₙ[6] <= Pᴳ[30,t])

    ineq_Q_up[(2,t)]  = @constraint(model, Qᴳ[2,t] <= yˢᴳ²[t]*Qˢᴳₘₐₓ[1])
    ineq_Q_lo[(2,t)]  = @constraint(model, yˢᴳ²[t]*Qˢᴳₘᵢₙ[1] <= Qᴳ[2,t])
    ineq_Q_up[(3,t)]  = @constraint(model, Qᴳ[3,t] <= yˢᴳ³[t]*Qˢᴳₘₐₓ[2])
    ineq_Q_lo[(3,t)]  = @constraint(model, yˢᴳ³[t]*Qˢᴳₘᵢₙ[2] <= Qᴳ[3,t])
    ineq_Q_up[(4,t)]  = @constraint(model, Qᴳ[4,t] <= yˢᴳ⁴[t]*Qˢᴳₘₐₓ[3])
    ineq_Q_lo[(4,t)]  = @constraint(model, yˢᴳ⁴[t]*Qˢᴳₘᵢₙ[3] <= Qᴳ[4,t])
    ineq_Q_up[(5,t)]  = @constraint(model, Qᴳ[5,t] <= yˢᴳ⁵[t]*Qˢᴳₘₐₓ[4])
    ineq_Q_lo[(5,t)]  = @constraint(model, yˢᴳ⁵[t]*Qˢᴳₘᵢₙ[4] <= Qᴳ[5,t])
    ineq_Q_up[(27,t)] = @constraint(model, Qᴳ[27,t] <= yˢᴳ²⁷[t]*Qˢᴳₘₐₓ[5])
    ineq_Q_lo[(27,t)] = @constraint(model, yˢᴳ²⁷[t]*Qˢᴳₘᵢₙ[5] <= Qᴳ[27,t])
    ineq_Q_up[(30,t)] = @constraint(model, Qᴳ[30,t] <= yˢᴳ³⁰[t]*Qˢᴳₘₐₓ[6])
    ineq_Q_lo[(30,t)] = @constraint(model, yˢᴳ³⁰[t]*Qˢᴳₘᵢₙ[6] <= Qᴳ[30,t])

    soc_gen[(2,t)]  = @constraint(model, [1, Pᴳ[2,t]/Sᵐᵃˣ_gc[1], Qᴳ[2,t]/Sᵐᵃˣ_gc[1]] in SecondOrderCone())
    soc_gen[(3,t)]  = @constraint(model, [1, Pᴳ[3,t]/Sᵐᵃˣ_gc[2], Qᴳ[3,t]/Sᵐᵃˣ_gc[2]] in SecondOrderCone())
    soc_gen[(4,t)]  = @constraint(model, [1, Pᴳ[4,t]/Sᵐᵃˣ_gc[3], Qᴳ[4,t]/Sᵐᵃˣ_gc[3]] in SecondOrderCone())
    soc_gen[(5,t)]  = @constraint(model, [1, Pᴳ[5,t]/Sᵐᵃˣ_gc[4], Qᴳ[5,t]/Sᵐᵃˣ_gc[4]] in SecondOrderCone())
    soc_gen[(27,t)] = @constraint(model, [1, Pᴳ[27,t]/Sᵐᵃˣ_gc[5], Qᴳ[27,t]/Sᵐᵃˣ_gc[5]] in SecondOrderCone())
    soc_gen[(30,t)] = @constraint(model, [1, Pᴳ[30,t]/Sᵐᵃˣ_gc[6], Qᴳ[30,t]/Sᵐᵃˣ_gc[6]] in SecondOrderCone())
    soc_gen[(1,t)]  = @constraint(model, [1, Pᴳ[1,t]/Sᵐᵃˣ_gv[1], Qᴳ[1,t]/Sᵐᵃˣ_gv[1]] in SecondOrderCone())
    soc_gen[(23,t)] = @constraint(model, [1, Pᴳ[23,t]/Sᵐᵃˣ_c[1], Qᴳ[23,t]/Sᵐᵃˣ_c[1]] in SecondOrderCone())
    soc_gen[(24,t)] = @constraint(model, [1, Pᴳ[24,t]/Sᵐᵃˣ_c[2], Qᴳ[24,t]/Sᵐᵃˣ_c[2]] in SecondOrderCone())
end

# --------------------------
# 3. Power balance & line flow
# --------------------------
for i in 1:bus_num, t in 1:T
    eq_p_bal[(i,t)] = @constraint(model, Pᴳ[i,t] - Pᴰ[i,t]*1.1 == sum(Pˡⁱⁿᵉ[i,j,t] for j in bus_connection[i]))
    eq_q_bal[(i,t)] = @constraint(model, Qᴳ[i,t] - Qᴰ[i,t] == sum(Qˡⁱⁿᵉ[i,j,t] for j in bus_connection[i]))
    for j in bus_connection[i]
        eq_anti_p[(i,j,t)] = @constraint(model, Pˡⁱⁿᵉ[i,j,t] + Pˡⁱⁿᵉ[j,i,t] == 0)
        eq_anti_q[(i,j,t)] = @constraint(model, Qˡⁱⁿᵉ[i,j,t] + Qˡⁱⁿᵉ[j,i,t] == 0)
    end
end

# --------------------------
# 4. Non-generator zero injection
# --------------------------
gens = [1,2,3,4,5,23,24,27,30]
for i in 1:bus_num, t in 1:T
    if !(i in gens)
        eq_zero_p[(i,t)] = @constraint(model, Pᴳ[i,t] == 0)
        eq_zero_q[(i,t)] = @constraint(model, Qᴳ[i,t] == 0)
    end
end

# --------------------------
# 5. Voltage stability & logic
# --------------------------
for t in 1:T
    for m in 1:15
        ineq_eta[(m,t)] = @constraint(model, ηₘ_1[m,t] >= 0)
    end
    for b in 1:6
        eq_eta2[(b,t)] = @constraint(model, ηₘ_2[b,t] == 0)
    end
    eq_gamma23[t] = @constraint(model, Γ_23[t] == 0)
    eq_gamma24[t] = @constraint(model, Γ_24[t] == 0)
    soc_vs[t] = @constraint(model, [1,1,1] in SecondOrderCone())
end

# ======================================================================
# ✅ FINAL DICTIONARY (ALL VARIABLES + ALL NAMED CONSTRAINTS)
# ======================================================================
uc_dict = OrderedDict{Symbol, Any}()

uc_dict[:vars] = OrderedDict{Symbol, Any}()
uc_dict[:eq_const] = OrderedDict{Symbol, Any}()
uc_dict[:ineq_const] = OrderedDict{Symbol, Any}()
uc_dict[:soc_const] = OrderedDict{Symbol, Any}()

# Variables
uc_dict[:vars][:y_SG] = (yˢᴳ², yˢᴳ³, yˢᴳ⁴, yˢᴳ⁵, yˢᴳ²⁷, yˢᴳ³⁰)
uc_dict[:vars][:CU] = (Cᵁ², Cᵁ³, Cᵁ⁴, Cᵁ⁵, Cᵁ²⁷, Cᵁ³⁰)
uc_dict[:vars][:P_G] = Pᴳ
uc_dict[:vars][:Q_G] = Qᴳ
uc_dict[:vars][:P_line] = Pˡⁱⁿᵉ
uc_dict[:vars][:Q_line] = Qˡⁱⁿᵉ
uc_dict[:vars][:Gamma] = (Γ_23, Γ_24)
uc_dict[:vars][:eta_m1] = ηₘ_1
uc_dict[:vars][:eta_m2] = ηₘ_2

# Equality constraints
uc_dict[:eq_const][:p_balance] = eq_p_bal
uc_dict[:eq_const][:q_balance] = eq_q_bal
uc_dict[:eq_const][:anti_p] = eq_anti_p
uc_dict[:eq_const][:anti_q] = eq_anti_q
uc_dict[:eq_const][:zero_p] = eq_zero_p
uc_dict[:eq_const][:zero_q] = eq_zero_q
uc_dict[:eq_const][:gamma23] = eq_gamma23
uc_dict[:eq_const][:gamma24] = eq_gamma24
uc_dict[:eq_const][:eta2] = eq_eta2

# Inequality constraints
uc_dict[:ineq_const][:startup] = ineq_startup
uc_dict[:ineq_const][:P_up] = ineq_P_up
uc_dict[:ineq_const][:P_lo] = ineq_P_lo
uc_dict[:ineq_const][:Q_up] = ineq_Q_up
uc_dict[:ineq_const][:Q_lo] = ineq_Q_lo
uc_dict[:ineq_const][:eta_logic] = ineq_eta

# SOC constraints
uc_dict[:soc_const][:gen_cap] = soc_gen
uc_dict[:soc_const][:voltage_stability] = soc_vs

# ======================================================================
# Objective & Solve
# ======================================================================
cost_onoff_Primal = sum(Cᵁ²)+sum(Cᵁ³)+sum(Cᵁ⁴)+sum(Cᵁ⁵)+sum(Cᵁ²⁷)+sum(Cᵁ³⁰)
cost_nl_Primal = sum(Oⁿˡ[1].*yˢᴳ²)+sum(Oⁿˡ[2].*yˢᴳ³)+sum(Oⁿˡ[3].*yˢᴳ⁴)+sum(Oⁿˡ[4].*yˢᴳ⁵)+sum(Oⁿˡ[5].*yˢᴳ²⁷)+sum(Oⁿˡ[6].*yˢᴳ³⁰)
cost_gene_Primal = sum(Oᵐ[1].*Pᴳ[2,:])+sum(Oᵐ[2].*Pᴳ[3,:])+sum(Oᵐ[3].*Pᴳ[4,:])+sum(Oᵐ[4].*Pᴳ[5,:])+sum(Oᵐ[5].*Pᴳ[27,:])+sum(Oᵐ[6].*Pᴳ[30,:])
obj = cost_onoff_Primal + cost_nl_Primal + cost_gene_Primal

@objective(model, Min, obj)
set_optimizer(model, Gurobi.Optimizer)
set_optimizer_attribute(model, "QCPDual", 1)
optimize!(model)