# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Pricing voltage stability
# 4.July.2025


#------------------------------------------------------------------------------
#--------Price voltage stability with restricted method------------------------
#------------------------------------------------------------------------------



import Pkg
Pkg.add("SpecialFunctions")
include("levy.jl")
include("HO.jl")
include("fun_info.jl")

using JuMP,Gurobi, CSV,DataFrames,LinearAlgebra, XLSX, IterTools, DelimitedFiles,Plots,MAT,Distributions,Random,SpecialFunctions 
#Optimization, OptimizationMetaheuristics


#-------------------------------------------------------
#-------------------Define model-------------------
#-------------------------------------------------------


#-------------------Loop for Heuristic-------------------

#.  demo
rosenbrock(x, p) = -(p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2 +5^x[1]
x0 = zeros(2)
p = [1.0, 100.0]
f = OptimizationFunction(rosenbrock)
prob = SciMLBase.OptimizationProblem(f, x0, p, lb = [-1.0, -1.0], ub = [1.0, 1.0])
sol = solve(prob, WOA(), maxiters = 100000, maxtime = 1000.0)



SearchAgents=16;                     # number of Hippopotamus (population members)
Max_iterations=100;                     # maximum number of iteration
lowerbound,upperbound,dimension,fitness=fun_info("F1");                     # Object function
Best_score,Best_pos,HO_curve =HO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);



obj(x)= x[1] + x[2]
x0 = zeros(2)
p = [1.0, 100.0]
f = OptimizationFunction(obj)
prob = SciMLBase.OptimizationProblem(f, x0, p, lb = [-1.0, -1.0], ub = [1.0, 1.0])
sol = solve(prob, ECA(), maxiters = 100000, maxtime = 1000.0)

