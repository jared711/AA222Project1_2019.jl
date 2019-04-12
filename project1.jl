using Pkg
Pkg.add("Plots")
Pkg.add("PyPlot")
using Plots
pyplot(size = (900,900), legend = false)

using LinearAlgebra
using Statistics

X = -15:0.1:15
Y = -15:0.1:15
Z = zeros(length(X), length(Y))
for i in 1: length(X)
    for j in 1: length(Y)
        x_state = [X[i], Y[j]]
        Z[i,j] = rosenbrock(x_state)
    end
end

contour(X,Y,Z)
scatter!([1,1],[1,1],marker = (:circle, 5, 0.6, :green, stroke(3, 0.2, :green, :dot)))

x0 = [0.0,0.0]
optimize(rosenbrock,myfun,x0,100,"hello")

function myfun(x)
    return x
end

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
        scatter!([x[1],x[1]],[x[2],x[2]],marker = (:circle, 5, 0.6, :magenta, stroke(3, 0.2, :magenta, :dot)))
        #x_list[list_counter,:], f_list[list_counter] = x, y
        #list_counter += 1

        counter -= 1
        while α > ϵ && counter > 2*n
            improved = false
            x_best, y_best = x, y
            for i in 1 : n
                for sgn in (-1, 1)
                    x′ = x + sgn*α*basis(i, n)
                    y′ = f(x′)
                    scatter!([x'[1],x'[1]],[x'[2],x'[2]],marker = (:circle, 5, 0.6, :blue, stroke(3, 0.2, :blue, :dot)))
                    #x_list[list_counter,:], f_list[list_counter] = x, y
                    #list_counter += 1

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

    scatter!([x[1],x[1]],[x[2],x[2]],marker = (:circle, 5, 0.6, :red, stroke(3, 0.2, :red, :dot)))
    return x
end

function rosenbrock(x)
    x_star = ([1,1])
    f_star = 0
    n = 20
    x1 = x[1]
    x2 = x[2]
    f = 100*(x2-x1^2)^2 + (1-x1)^2;
    return f
end


