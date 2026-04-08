function HO(SearchAgents, Max_iterations, lowerbound, upperbound, dimension, fitness)

    lowerbound_vec = fill(lowerbound, dimension)
    upperbound_vec = fill(upperbound, dimension)

    # 初始化
    X = [lowerbound_vec[j] + rand() * (upperbound_vec[j] - lowerbound_vec[j])
         for i in 1:SearchAgents, j in 1:dimension]

    fit = [fitness(X[i, :]) for i in 1:SearchAgents]

    best_so_far = zeros(Float64, Max_iterations)
    Xbest = zeros(dimension)
    fbest = Inf

    for t in 1:Max_iterations

        best, location = findmin(fit)
        if t == 1 || best < fbest
            fbest = best
            Xbest = copy(X[location, :])
        end

        # -------------------------
        # Phase 1 – Exploration
        # -------------------------
        for i in 1:Int(SearchAgents/2)

            Dominant = Xbest
            I1 = rand(1:2)
            I2 = rand(1:2)
            Ip1 = rand(0:1, 2)

            RandGroupNumber = rand(1:SearchAgents)
            RandGroup = randperm(SearchAgents)[1:RandGroupNumber]

            if length(RandGroup) == 1
                MeanGroup = X[RandGroup[1], :]
            else
                MeanGroup = sum(X[RandGroup, :], dims=1)[:] / length(RandGroup)
            end

            # Alfa 参数集
            Alfa = Dict{Int,Any}()
            Alfa[1] = I2 .* rand(1, dimension) .+ (1 - Ip1[1])
            Alfa[2] = 2 .* rand(1, dimension) .- 1
            Alfa[3] = rand(1, dimension)
            Alfa[4] = I1 .* rand(1, dimension) .+ (1 - Ip1[2])
            Alfa[5] = rand()

            A = Alfa[rand(1:5)]
            B = Alfa[rand(1:5)]

            # X_P1
            X_P1 = X[i, :] .+ rand() .* (Dominant .- I1 .* X[i, :])

            # X_P2
            T = exp(-t / Max_iterations)
            if T > 0.6
                X_P2 = X[i, :] .+ A .* (Dominant .- I2 .* MeanGroup)
            else
                if rand() > 0.5
                    X_P2 = X[i, :] .+ B .* (MeanGroup .- Dominant)
                else
                    X_P2 = lowerbound_vec .+ rand() .* (upperbound_vec .- lowerbound_vec)
                end
            end

            X_P2 = clamp.(X_P2, lowerbound_vec, upperbound_vec)

            # 更新
            F_P1 = fitness(X_P1)
            if F_P1 < fit[i]
                X[i, :] = X_P1
                fit[i] = F_P1
            end

            F_P2 = fitness(X_P2)
            if F_P2 < fit[i]
                X[i, :] = X_P2
                fit[i] = F_P2
            end
        end

        # -------------------------
        # Phase 2 – Predator Defense
        # -------------------------
        RL = 0.05 .* levy(SearchAgents, dimension, 1.5)

        for i in Int(SearchAgents/2)+1:SearchAgents
            predator = lowerbound_vec .+ rand(dimension) .* (upperbound_vec .- lowerbound_vec)
            F_HL = fitness(predator)

            dist = abs.(predator .- X[i, :])

            b = rand(Uniform(2, 4))
            c = rand(Uniform(1, 1.5))
            d = rand(Uniform(2, 3))
            l = rand(Uniform(-2*pi, 2*pi))

            if fit[i] > F_HL
                X_P3 = RL[i, :] .* predator + (b / (c - d * cos(l))) .* (1 ./ dist)
            else
                X_P3 = RL[i, :] .* predator + (b / (c - d * cos(l))) .* (1 ./ (2 .* dist .+ rand(dimension)))
            end

            X_P3 = clamp.(X_P3, lowerbound_vec, upperbound_vec)

            F_P3 = fitness(X_P3)
            if F_P3 < fit[i]
                X[i, :] = X_P3
                fit[i] = F_P3
            end
        end

        # -------------------------
        # Phase 3 – Escape (Exploitation)
        # -------------------------
        for i in 1:SearchAgents
            LO_LOCAL = lowerbound_vec ./ t
            HI_LOCAL = upperbound_vec ./ t

            Alfa = Dict{Int,Any}()
            Alfa[1] = 2 .* rand(1,dimension) .- 1
            Alfa[2] = rand()
            Alfa[3] = randn()

            D = Alfa[rand(1:3)]
            X_P4 = X[i, :] .+ rand() .* (LO_LOCAL .+ D .* (HI_LOCAL .- LO_LOCAL))
            X_P4 = clamp.(X_P4, lowerbound_vec, upperbound_vec)

            F_P4 = fitness(X_P4)
            if F_P4 < fit[i]
                X[i, :] = X_P4
                fit[i] = F_P4
            end
        end

        best_so_far[t] = fbest
        println("Iteration $t: Best Cost = $(best_so_far[t])")
    end

    return fbest, Xbest, best_so_far
end