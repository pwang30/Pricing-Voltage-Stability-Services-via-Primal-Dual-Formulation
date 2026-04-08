function build_GB_matrices(file_path::String)



    data = CSV.read(file_path, DataFrame)

    nbus = maximum(union(data.fbus, data.tbus))
    Ybus = zeros(ComplexF64, nbus, nbus)

    for row in eachrow(data)
        f = row.fbus
        t = row.tbus
        r = row.r
        x = row.x

        # 跳过关闭的线路或非法阻抗
        if row.status != 1 || (r == 0 && x == 0)
            continue
        end

        z = r + im * x
        y = 1 / z

        # 构建 Ybus（无 shunt）
        Ybus[f, f] += y
        Ybus[t, t] += y
        Ybus[f, t] -= y
        Ybus[t, f] -= y
    end

    # 提取实部与虚部
    G = real.(Ybus)
    B = imag.(Ybus)

    return G, B
end