using LinearAlgebra
using Statistics

function optimize(f, g, x0, n, prob)
    if prob == "secret_2" || prob == "secret_1"
        basis(i, m) = [k == i ? 1.0 : 0.0 for k in 1 : m]
    
        S = [basis(i, length(x0)) * x0[i] for i in 1:length(x0)]
        push!(S, randn(length(x0)))
    
        x = nelder_mead(f, S, 1e-8, n)
        return x
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
                                  

function nelder_mead(f, S, ϵ, n; α=1.0, β=2.0, γ=0.5)
    counter = n
    Δ, y_arr = Inf, f.(S)
    counter -= length(S)
    while Δ > ϵ && counter > length(S) #+ 1 #maybe not necessary
        p = sortperm(y_arr) #sort lowest to highest
        S, y_arr = S[p], y_arr[p]
        xl, yl = S[1], y_arr[1] #lowest
        xh, yh = S[end], y_arr[end] #highest
        xs, ys = S[end-1], y_arr[end-1] #second-highest
        xm = mean(S[1:end-1]) #centroid
        xr = xm + α*(xm - xh) #reflection point
        yr = f(xr)
        counter -= 1 #MINUS 1
        
        if yr < yl
            xe = xm + β*(xr-xm) #expansion point
            ye = f(xe)
            counter -= 1 #MINUS 1
            S[end], y_arr[end] = ye < yr ? (xe, ye) : (xr, yr)
        elseif yr > ys
            if yr ≤ yh
                xh, yh, S[end], y_arr[end] = xr, yr, xr, yr
            end
            xc = xm + γ*(xh-xm) #contraction point
            yc = f(xc)
            counter -= 1 #MINUS 1
            if yc > yh
                for i in 2 : length(y_arr)
                    S[i] = (S[i] + xl)/2
                    y_arr[i] = f(S[i])
                    counter -= length(y_arr)-1 #MINUS (length(y_arr)-1)
                end
            else
                S[end], y_arr[end] = xc, yc
            end
        else
            S[end], y_arr[end] = xr, yr
        end

        Δ = std(y_arr, corrected=false)
    end
    return S[argmin(y_arr)]
end