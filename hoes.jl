# File contains an implementation of polynomial annihilation based High Order Edge Sensors (HOES)

function hoes(grid, AR, i)
	if i < 2*HOESorder 
                localgrid = grid[1: 2*HOESorder+1]
                ls = AR[1:2*HOESorder+1]
        elseif i > length(grid) - 2*HOESorder
                localgrid = grid[end-2*HOESorder:end]
                ls = AR[end-2*HOESorder:end]
        else
                localgrid = grid[i-HOESorder:i+HOESorder]
                ls = AR[i-HOESorder:i+HOESorder]
        end

        w = ones(HOESorder*2+1)
        c = zeros(HOESorder*2+1)
        for j in 1:(2*HOESorder+1)
                for i in 1:j-1
                        w[j] = w[j]*(localgrid[j]-localgrid[i])
                end
                for i in (j+1):(2*HOESorder+1)
                        w[j] = w[j]*(localgrid[j]-localgrid[i])
                end
        end
        c = factorial(2*HOESorder)./w
        q = [sum(c[j:end]) for j in 1:(2*HOESorder+1)]
        L = sum(c.*ls) / q[HOESorder+1]
	return L
end

function hoes(u; grid = collect(1:length(u)))
	umax = maximum(u)
	umin = minimum(u)
	utest = vcat(ones(10)*umax, ones(10)*umin)
	bn = maximum(abs.(map(i->hoes(grid[1:20], utest, i), 1:20)))
	b = abs.(map(i->hoes(grid, u,i), 1:length(u)))	
        return b, bn
end

