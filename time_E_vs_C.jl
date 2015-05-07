#using EfficientLorentz
using ClassicalLorentz

#x = [0.532282, 0.531928, 0.267589]; v = [0.665341, 0.668661, 0.331986] - bad velocity (too symmetric). Below is another.
#x = [0.532282, 0.531928, 0.267589]; v = [cos(1)*sin(0.8), sin(1)*sin(0.8), cos(0.8)]


################ E vs C 3D ################

#= Measuring first collision execution time as a function of radius
time_vs_r = Array{Real, 1}[]

@elapsed collisions3d_time(x, v, 0.11, 1, 64)
for n = 1:7
	r0 = 1/10^n
	for i = 0:9
		r = r0 - i/10^(n+1)
		time_E = (@elapsed collisions3d_time(x, v, r, 1, 64))
		println(r, " ", time_E)
		push!(time_vs_r, [r, time_E])
	end
	time_vs_r
end
=#

#= Measuring first collision execution time (classical) as a function of radius
@elapsed first_collision3d_classical(x, v, 0.11, 64)
for n = 1:7
	r0 = 1/10^n
	for i = 0:9
		r = r0 - i/10^(n+1)
		time_E = (@elapsed first_collision3d_classical(x, v, r, 64))
		println(r, " ", time_E)
		push!(time_vs_r, [r, time_E])
	end
	time_vs_r
end
=#


################ E vs C 2D ################

# Measuring first collision execution time as a function of radius (2D)
function exec_time_vs_r_E(x, v)
    time_vs_r = Array{Real, 1}[]
    @elapsed collisions(x[1], x[2], v[1], v[2], 0.11, 1, 64)
    for n = 1:7
	    r0 = 1/10^n
	    for i = 0:9
		    r = r0 - i/10^(n+1)
		    time_E = (@elapsed collisions(x[1], x[2], v[1], v[2], r, 1, 64))
		    #println(r, " ", time_E)
		    push!(time_vs_r, [r, time_E])
	    end
    end
    time_vs_r
end

#= Code to average all points over angles 1..44, 46..89
x = [0.01, 0.445]
time_vs_r_deg = Array{Array{Real, 1}, 1}[]
for deg = 1:89
    if deg != 45
        phi = (pi/180)*deg
        v = [cos(phi), sin(phi)]
        push!(time_vs_r_deg, exec_time_vs_r_E(x, v))
    end        
end

averaged_points = Array{Real, 1}[]
for j = 1:70
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

# Measuring first collision execution time (classical) as a function of radius (2D)
#x = [0.01, 0.445]; v = [cos(1), sin(1)]

function exec_time_vs_r_C(x, v)
    time_vs_r = Array{Real, 1}[]
    #@elapsed first_collision_classical(x, v, 0.11, 64)
    for n = 1:6
	    r0 = 1/10^n
	    for i = 0:9
		    r = r0 - i/10^(n+1)
		    time_E = (@elapsed first_collision_classical(x, v, r, 64))
		    #println(r, " ", time_E)
		    push!(time_vs_r, [r, time_E])
	    end
	end
    time_vs_r
end

# Code to average all points over angles 1..44, 46..89
x = [0.01, 0.445]
time_vs_r_deg = Array{Array{Real, 1}, 1}[]
for deg = 1:89
    if deg != 45
        phi = (pi/180)*deg
        v = [cos(phi), sin(phi)]
        push!(time_vs_r_deg, exec_time_vs_r_C(x, v))
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


