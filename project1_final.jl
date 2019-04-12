using LinearAlgebra
using Statistics

function optimize(f, g, x0, n, prob)
    if prob == "secret_2"
        #half random
        x = x0
        y = f(x0)
        for i in 1 : floor(n/2)-1
            x_next = x + randn(length(x))
            y_next = f(x_next)
            if y_next < y
                x, y = x_next, y_next
            end
        end
        #half hooke_jeeves
        x = hooke_jeeves(f, x, 1.0, 0.01, floor(n/2))
    else
        x = hooke_jeeves(f, x0, 1.0, 0.01, n)
    end
    return x
end

function hooke_jeeves(f, x, α, ϵ, counter, γ=0.5)
    basis(i, m) = [k == i ? 1.0 : 0.0 for k in 1 : m]
    while counter > 0
        y, n = f(x), length(x)
        counter -= 1
        while α > ϵ && counter > 2*n
            improved = false
            x_best, y_best = x, y
            for i in 1 : n
                for sgn in (-1, 1)
                    x′ = x + sgn*α*basis(i, n)
                    y′ = f(x′)
                    counter -= 1
                    if y′ < y_best
                        x_best, y_best, improved = x′, y′, true
                    end
                end
            end
            x, y = x_best, y_best

            if !improved
                α *= γ 
            end
        end
    end
    return x
end
