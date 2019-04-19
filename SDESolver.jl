using Plots, Random

#solve SDE using schemeEuler–Maruyama dX = a(t,x)*dt + b(t,x)*dw
function SDEsolver(fa, fb, t0, tn, N, y0)
    t_init = t0
    t_end = tn
    dt = (t_end - t_init) / N
    y_init = y0

    function dW(delta_t)
        return (rand()*2 - 1)* sqrt(delta_t)
    end

    ts = t_init : dt : t_end - dt
    ys = zeros(N)
    ys[1] = y_init

    for i in 2 : N
        t = ts[i - 1]
        y = ys[i - 1]
        ys[i] = y + fa(y, t) * dt + fb(y, t) * dW(dt)
    end
    return (ts, ys)
end


c_theta = 0.7
c_mu = 0.2
c_sigma = 0.06

function mu(y, t)
    return c_theta * (c_mu - y)
end

function sigma(y, t)
    return c_sigma
end

sol = SDEsolver(mu, sigma, -5, 5, 10000, 0)
plot(sol)
