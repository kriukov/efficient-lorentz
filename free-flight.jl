using EfficientLorentz

################ Free flight time ################

## 2D

# Fix the velocity and vary r, checking the longest free flight as a function of r

# Measuring free flight time until the first collision in 2D
#x = [0.01, 0.445]
#phi = (1 + sqrt(5))/2

function freeflight2d(x, v)
    time_to_1st = Array{Real, 1}[]
    for n = 1:8
	    r0 = 1/10^n
	    for i = 1:9
	        #for j = 1:9 # This was for more data points
		    r = r0 - i/10^(n+1) # - j/10^(n+2)
		    first_place = collisions(x[1], x[2], v[1], v[2], r, 1, 256)[1][1]
		    dist_to_1st = norm(first_place - x)
		    #println(r, " ", dist_to_1st)
		    push!(time_to_1st, [r, dist_to_1st])
		    #end
	    end
	end
    time_to_1st
end

#= Code to average points over velocities with angles 1..44, 46..89
x = [0.01, 0.445]
freeflight2d_deg = Array{Array{Real, 1}, 1}[]
for deg = 1:89
    if deg != 45
        phi = (pi/180)*deg
        v = [cos(phi), sin(phi)]
        push!(freeflight2d_deg, freeflight2d(x, v))
    end        
end
averaged_points = Array{Real, 1}[]
for j = 1:60
    s = 0
    l = length(freeflight2d_deg)
    for i = 1:l
        s += freeflight2d_deg[i][j][2]    
    end
    point = [freeflight2d_deg[1][j][1], s/l]
    push!(averaged_points, point)
end
println(averaged_points)
=#

# Fix the velocity and randomize initial position, then choose the maximum free flight of all
#v = [sqrt(1/3), sqrt(2/3)] # Slope of velocity = sqrt(2)
#k = 1 - 0.2/sqrt(2); b = 0.1/sqrt(2)
#=
function freeflight2d_randomxmax(v)
    time_to_1st = Array{Real, 1}[]
    for n = 6:8
	    r0 = 1/10^n
	    for i = 1:9
		    r = r0 - i/10^(n+1)
	        # For this given radius, find free flight for random x (in reduced 0-1 square) and take the maximum
	        freeflights = Real[]
	        for i = 1:1000
	            #x = [k*rand() + b, k*rand() + b]
	            x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
		        first_place = collisions(x[1], x[2], v[1], v[2], r, 1, 256)[1][1]
		        dist_to_1st = norm(first_place - x)
		        push!(freeflights, dist_to_1st)
		    end
		    max_dist_to_1st = findmax(freeflights)[1]
		    push!(time_to_1st, [r, max_dist_to_1st])
	    end
	end
    time_to_1st
end

phi = (1 + sqrt(5))/2
v = [cos(phi), sin(phi)]
data = freeflight2d_randomxmax(v)
println(data)
=#


## 3D

#= Measuring free flight time until the first collision (basically distance between x0 and first place of collision)

time_to_1st = Array{Real, 1}[]
for n = 1:5
	r0 = 1/10^n
	for i = 0:9
		r = r0 - i/10^(n+1)
		first_place = collisions3d_time(x, v, r, 1, 64)[2][1]
		dist_to_1st = norm(first_place - x)
		println(r, " ", dist_to_1st)
		push!(time_to_1st, [r, dist_to_1st])
	end
	time_to_1st
end
=#

# Trying out different velocities, at least in the first octant
time_to_1st = Array{Real, 1}[]

#for phi1 = 1:9
#for theta1 = 1:9

#phi = phi1*(pi/2)/10
#theta = theta1*(pi/2)/10

#phi = 2pi*rand()
#theta = pi*rand()
#=
phi = sqrt(2)
theta = (1 + sqrt(5))/2
v = [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
=#

phi = (1 + sqrt(5))/2
v = [1/(phi + 2), phi/(phi + 2), phi/sqrt(phi + 2)]

for n = 3:4
	r0 = 1/10^n
	for i = 1:9
	    for k = 1:2
		    r = r0 - i/10^(n+1) + (k-1)/(2*10^(n+1))
		    freeflights = Real[]
		    exec_t = 0
		    N = 100
		    for j = 1:N
        		x = [0.8*rand() + 0.1, 0.8*rand() + 0.1, 0.8*rand() + 0.1]
        		tic()
	        	first_place = collisions3d_time(x, v, r, 1, 64)[2][1]
	        	exec_t += toq()
		        dist_to_1st = norm(first_place - x)
		        push!(freeflights, dist_to_1st)
		    end
		    exec_t /= N
		    max_dist_to_1st = findmax(freeflights)[1]
		    sigma = sqrt(var(freeflights))
		    println(r, " ", max_dist_to_1st, " ", sigma, "; avg. exectime = ", exec_t)
		    push!(time_to_1st, [r, max_dist_to_1st, sigma, exec_t])
        end
	end
	time_to_1st
end

#end
#end


#= Collecting all points instead of taking the maximum
phi = (1 + sqrt(5))/2
v = [1/(phi + 2), phi/(phi + 2), phi/sqrt(phi + 2)]

freeflightsdist = Array{BigFloat, 1}[]
for n = 3:4
	r0 = 1/10^n
	for i = 1:9
		r = r0 - i/10^(n+1)
		for i = 1:100
    		x = [0.8*rand() + 0.1, 0.8*rand() + 0.1, 0.8*rand() + 0.1]
	    	first_place = collisions3d_time(x, v, r, 1, 64)[2][1]
		    dist_to_1st = norm(first_place - x)
		    println(r, " ", dist_to_1st)
		    push!(freeflightsdist, [r, dist_to_1st])
		end
	end
	println(freeflightsdist)
end
=#

