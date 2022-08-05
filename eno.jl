# File contains ENO style reconstructions of order 2


# Performs classic ENO reconstruction of order 2, u is the vector [u_{i-1}, u_i, u_{i+1}]
# returns mu, du/dx
@inline function reconENO2(u)
        a1 = u[2] - u[1]
        a2 = u[3] - u[2]
        if a1^2 < a2^2
                return u[2], a1
        else
                return u[2], a2
        end
end

@inline function classicENO2(u)
        mean, a = reconENO2(u)
        return mean - a/2, mean + a/2
end

@inline function ENO2Mat(u)
        N, K = size(u)
        ul = zeros(N, K)
        ur = zeros(N, K)
        for i=2:N-1
                for k = 1:K
                        ul[i, k], ur[i, k] = classicENO2(u[i-1:i+1, k])
                end
        end
	# Alternative: Reconstruction with only the inner stencil
	ul[1, :] = u[1, :]
	ur[1, :] = u[1, :]
	ul[N, :] = u[N, :]
	ur[N, :] = u[N, :]
        return ul, ur
end

@inline minmod(a::Float64, b::Float64) = min(abs(a), abs(b))*(sign(a) + sign(b))/2.0 :: Float64


@inline function MUSCLMat(u::Matrix{Float64}) :: Tuple{Matrix{Float64}, Matrix{Float64}}
	N, K = size(u)
        ul = zeros(N, K)
        ur = zeros(N, K)
#Threads.@threads  
	for i=1:N
		for k=1:K
			a1 = u[mod1(i, N), k] - u[mod1(i-1, N), k]
			a2 = u[mod1(i+1, N), k] - u[mod1(i, N), k]
			a = minmod(a1, a2)
               		ul[i, k], ur[i, k] = u[i, k] - a/2.0, u[i, k] + a/2.0
		end
        end
        # Alternative: Reconstruction with only the inner stencil
     #   ul[1, :] = u[1, :]
      #  ur[1, :] = u[1, :]
      #  ul[N, :] = u[N, :]
      #  ur[N, :] = u[N, :]
        return ul, ur
end

function MUSCLMat!(u, ul, ur)
        N, K = size(u)
        for i=1:N
                for k = 1:K
                        a1 = u[mod1(i, N), k] - u[mod1(i-1, N), k]
                        a2 = u[mod1(i+1, N), k] - u[mod1(i, N), k]
                        a = minmod(a1, a2)
                        ul[i, k], ur[i, k] = u[i, k] - a/2, u[i, k] + a/2
                end
        end
end

MUSCLr(u) = u[2, :] + 0.5*minmod.(u[2, :] - u[1, :], u[3, :] - u[2, :])
MUSCLl(u) = u[2, :] - 0.5*minmod.(u[2, :] - u[1, :], u[3, :] - u[2, :])


