function fun_info(name)
    if name == "F1"
        fitness(x) = sum(x.^3)
        return -100, 100, 30, fitness
    else
        error("Unknown function name")
    end
end