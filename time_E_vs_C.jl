using EfficientLorentz
using ClassicalLorentz

#x = [0.532282, 0.531928, 0.267589]; v = [0.665341, 0.668661, 0.331986] - bad velocity (too symmetric). Below is another.
#x = [0.532282, 0.531928, 0.267589]; v = [cos(1)*sin(0.8), sin(1)*sin(0.8), cos(0.8)]

# Measuring first collision execution time as a function of radius (general function)
# x, v - init.cond, N - count down to r = 10^-N, alg - 'E' or 'C', dim - 2 or 3
function exectime_vs_r(x, v, N0, N, alg, dim, prec)
	if alg == 'E'
		if dim == 2
			f(r, prec) = collisions(x[1], x[2], v[1], v[2], r, 1, prec)
		elseif dim == 3
			f(r, prec) = collisions3d_time(x, v, r, 1, prec)
		else
			error("Dimensions must be 2 or 3")
		end
	elseif alg == 'C'
		if dim == 2
			f(r, prec) = first_collision_classical(x, v, r, prec)
		elseif dim == 3
			f(r, prec) = first_collision3d_classical(x, v, r, prec)
		else
			error("Dimensions must be 2 or 3")
		end
	else 
		error("Algorithm must be 'E' or 'C'")
	end
	time_vs_r = Array{Real, 1}[]
	@elapsed f(0.11, prec) #For some reason, first calc. takes much longer (probably initialization)
	push!(time_vs_r, [0.1, @elapsed f(0.1, prec)])
	for n = N0:N-1
	    r0 = 1/10^n
	    for i = 1:9
		    r = r0 - i/10^(n+1)
		    time_E = (@elapsed f(r, prec))
		    #println(r, " ", time_E)
		    push!(time_vs_r, [r, time_E])
	    end
    end
    time_vs_r
end


# Average runner to run many times the exec. time functions with many random init. cond. and take average results
# N - count down to r = 10^-N, alg - 'E' or 'C', dim - 2 or 3, times - how many init. cond.
function averagerunner(N0, N, alg, dim, times, prec)
	time_vs_r_many = Array{Array{Real, 1}, 1}[]
	for i = 1:times
		if dim == 2
			x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
			phi = 2pi*rand()
	        v = [cos(phi), sin(phi)]
		elseif dim == 3
			x = [0.8*rand() + 0.1, 0.8*rand() + 0.1, 0.8*rand() + 0.1]
			#phi = 2pi*rand(); theta = pi*rand()
			#v = [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
			phi = (1 + sqrt(5))/2
			v = [1/(phi + 2), phi/(phi + 2), phi/sqrt(phi + 2)]
		else 
			error("Only works for 2 or 3 dimensions")
		end
		
		push!(time_vs_r_many, exectime_vs_r(x, v, N0, N, alg, dim, prec))
	end
	averaged_points = Array{Real, 1}[]
	l = length(time_vs_r_many) #How many runs (random speeds/coords) we made
	l1 = length(time_vs_r_many[1]) #How many radii we ran through for each random speed/coord
	for j = 1:l1
		s = 0
		for i = 1:l
		    s += time_vs_r_many[i][j][2]    
		end
		point = [time_vs_r_many[1][j][1], s/l]
		push!(averaged_points, point)
	end

	averaged_points
end








###########################################################################################

### Below is old code that has been improved and cleaned up above


## Old separate functions that have been joined to exectime_vs_r()

#= Measuring first collision execution time as a function of radius (2D)
function exec_time_vs_r_E2D(x, v)
    time_vs_r = Array{Real, 1}[]
    @elapsed collisions(x[1], x[2], v[1], v[2], 0.11, 1, 64)
    push!(time_vs_r, [0.1, @elapsed collisions(x[1], x[2], v[1], v[2], 0.1, 1, 64)])
    for n = 1:7
	    r0 = 1/10^n
	    for i = 1:9
		    r = r0 - i/10^(n+1)
		    time_E = (@elapsed collisions(x[1], x[2], v[1], v[2], r, 1, 64))
		    #println(r, " ", time_E)
		    push!(time_vs_r, [r, time_E])
	    end
    end
    time_vs_r
end


# Measuring first collision execution time (classical) as a function of radius (2D)
function exec_time_vs_r_C2D(x, v)
    time_vs_r = Array{Real, 1}[]
    @elapsed first_collision_classical(x, v, 0.11, 64)
    push!(time_vs_r, [0.1, @elapsed first_collision_classical(x, v, 0.11, 64)])
    for n = 1:6
	    r0 = 1/10^n
	    for i = 1:9
		    r = r0 - i/10^(n+1)
		    time_E = (@elapsed first_collision_classical(x, v, r, 64))
		    #println(r, " ", time_E)
		    push!(time_vs_r, [r, time_E])
	    end
	end
    time_vs_r
end

# Measuring first collision execution time as a function of radius
function exec_time_vs_r_E3D(x, v)
	time_vs_r = Array{Real, 1}[]
	@elapsed collisions3d_time(x, v, 0.11, 1, 64) #For some reason, first calc. takes much longer (probably initialization)
	push!(time_vs_r, [0.1, @elapsed collisions3d_time(x, v, 0.1, 1, 64)]) 
	for n = 1:4
		r0 = 1/10^n
		for i = 1:9
			r = r0 - i/10^(n+1)
			time_E = (@elapsed collisions3d_time(x, v, r, 1, 64))
			println(r, " ", time_E)
			push!(time_vs_r, [r, time_E])
		end
	end
	time_vs_r
end


# Measuring first collision execution time (classical) as a function of radius
function exec_time_vs_r_C3D(x, v)
	time_vs_r = Array{Real, 1}[]
	@elapsed first_collision3d_classical(x, v, 0.11, 64) #For some reason, first calc. takes much longer (probably initialization)
	push!(time_vs_r, [0.1, @elapsed first_collision3d_classical(x, v, 0.1, 64)]) 
	for n = 1:4
		r0 = 1/10^n
		for i = 1:9
			r = r0 - i/10^(n+1)
			time_E = (@elapsed first_collision3d_classical(x, v, r, 64))
			println(r, " ", time_E)
			push!(time_vs_r, [r, time_E])
		end
	end
	time_vs_r
end
=#

## Code to average over fixed values

#= Code to average all points over angles 1..44, 46..89 with fixed x0
x = [0.01, 0.445]
time_vs_r_deg = Array{Array{Real, 1}, 1}[]
for deg = 1:89
    if deg != 45
    	#x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
        phi = (pi/180)*deg
        v = [cos(phi), sin(phi)]
        push!(time_vs_r_deg, exec_time_vs_r_E2D(x, v))
    end        
end
=#
#= Code to average all points over random x0 and v0
time_vs_r_deg = Array{Array{Real, 1}, 1}[]
for i = 1:100
    	x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
        phi = 2pi*rand()
        v = [cos(phi), sin(phi)]
        push!(time_vs_r_deg, exec_time_vs_r_E2D(x, v)) 
end
averaged_points = Array{Real, 1}[]
l = length(time_vs_r_deg) #How many init.cond.
l1 = length(time_vs_r_deg[1]) #How many radii for each init.cond., e.g., first
for j = 1:l1
    s = 0
    for i = 1:l
        s += time_vs_r_deg[i][j][2]    
    end
    point = [time_vs_r_deg[1][j][1], s/l]
    push!(averaged_points, point)
end
println(averaged_points)
=#
#= Code to average all points over angles 1..44, 46..89
x = [0.01, 0.445]
time_vs_r_deg = Array{Array{Real, 1}, 1}[]
for deg = 1:89
    if deg != 45
        phi = (pi/180)*deg
        v = [cos(phi), sin(phi)]
        push!(time_vs_r_deg, exec_time_vs_r_C2D(x, v))
    end        
end
averaged_points = Array{Real, 1}[]
for j = 1:60
    s = 0
    l = length(time_vs_r_deg)
    for i = 1:l
        s += time_vs_r_deg[i][j][2]    
    end
    point = [time_vs_r_deg[1][j][1], s/l]
    push!(averaged_points, point)
end
println(averaged_points)
=#

