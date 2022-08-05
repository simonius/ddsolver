using Polynomials

# File contains initial conditions that can be used to start the testcases and produces training initial data

# Function takes physical values like density, speed and pressure and returns the conserved variables.
function cons_var(rho, v, p)
    E = p / (gamma-1) + 0.5*v^2*rho
    rhov = rho*v
    return [rho, rhov, E]
end


function toro123(grid)
        N = size(grid)[1]
        u0 = zeros(N, 3)
        for i=1:N
                if grid[i] < 5.0
                        u0[i, :] = cons_var(1.0, -2.0, 0.4)
                else
                        u0[i, :] = cons_var(1.0, 2.0, 0.4)
                end
        end
        return u0
end


function shuosh6a(grid)
	N = size(grid)[1]
	u0 = zeros(N, 3)
	for i=1:N
		if grid[i] < 5.0
			u0[i, :] = cons_var(1.0, 0.0, 1.0)
		else
			u0[i, :] = cons_var(0.125, 0.0, 0.10)
		end
	end
	return u0
end

function shuosh6b(grid)
	 N = size(grid)[1]
        u0 = zeros(N, 3)
        for i=1:N
                if grid[i] < 5.0
                        u0[i, :] = cons_var(0.445, 0.698, 3.528)
                else
                        u0[i, :] = cons_var(0.5, 0.0, 0.571)
                end
        end
	return u0
end

function shuosh7(grid)
	N = size(grid)[1]
	u0 = zeros(N, 3)
	ul = cons_var(1.0, 0.0, 1.0E3)
	um = cons_var(1.0, 0.0, 1.0E-2)
	ur = cons_var(1.0, 0.0, 1.0E2)
	for i=1:N
		if grid[i] < 1.0
			u0[i, :] = ul
		elseif grid[i] < 9.0
			u0[i, :] = um
		else
			u0[i, :] = ur
		end
	end
	return u0
end

function shuosh8(grid)
        eps = 0.2
	N = size(grid)[1]
        u0 = zeros(N, 3)
        for i=1:N
                if grid[i] < 1.0
                        u0[i, :] = cons_var(3.857143, 2.629369, 10.33333333333)
                else
                        u0[i, :] = cons_var(1.0 + eps*sin(5*grid[i]), 0.0, 1.0)
                end
        end
        return u0
end


function rrinit(grid)
	rhol = 0.01+5.0*rand()
	vl = 5.0*(rand()-0.5)
	pl = 0.01+5.0*rand()

	rhoml = 0.01+5.0*rand()
        vml = 5.0*(rand()-0.5)
        pml = 0.01+5.0*rand()

        rhomr = 0.01+5.0*rand()
        vmr = 5.0*(rand()-0.5)
        pmr = 0.01+5.0*rand()


        rhor = 0.01+5.0*rand()
        vr = 5.0*(rand()-0.5)
        pr = 0.01+5.0*rand()

	N = size(grid)[1]
        u0 = zeros(N, 3)

        for i=1:N
		if grid[i] < 2.5
                	u0[i, :] = cons_var(rhol, vl, pl)
        	elseif grid[i] < 5.0 
			u0[i, :] = cons_var(rhoml, vml, pml)
		elseif grid[i] < 7.5
                        u0[i, :] = cons_var(rhomr, vmr, pmr)
                else
                        u0[i, :] = cons_var(rhor, vr, pr)
                end
	end
        return u0
end

uex(x, t) = cons_var(3.85 + exp(-(x-2*t-3.0)^2)*sin(2*(x-2*t)), 2.0, 10.3)
	
function smootheuler(grid)
        eps = 0.2
    	N = length(grid)
        u0 = zeros(N, 3)
        for i=1:N
                eps = exp(-(grid[i]-3.0)^2)
                u0[i, :] = cons_var(3.85 + eps*sin(2*grid[i]), 2.0, 10.3)
        end
        return u0
end
		
