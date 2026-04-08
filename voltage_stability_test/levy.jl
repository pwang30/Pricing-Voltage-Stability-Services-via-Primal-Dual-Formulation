

function levy(n, m, beta)
    num = gamma(1 + beta) * sin(pi * beta / 2)
    den = gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)
    sigma_u = (num / den)^(1 / beta)

    u = rand.(Normal(0, sigma_u), n, m)
    v = rand.(Normal(0, 1), n, m)

    z = u ./ (abs.(v) .^ (1 / beta))
    return z
end