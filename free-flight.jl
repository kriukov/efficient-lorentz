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

# Code to average points over velocities with angles 1..44, 46..89
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


## 3D

#= Measuring free flight time until the first collision (basically distance between x0 and first place of collision)
time_to_1st = Array{Real, 1}[]
for n = 1:7
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

#= Trying out different velocities, at least in the first octant
time_to_1st = Array{Real, 1}[]

for phi1 = 1:9
for theta1 = 1:9

phi = phi1*(pi/2)/10
theta = theta1*(pi/2)/10

v = [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
for n = 1:4
	r0 = 1/10^n
	for i = 0:9
		r = r0 - i/10^(n+1)
		first_place = collisions3d_time(x, v, r, 1, 64)[2][1]
		dist_to_1st = norm(first_place - x)
		println(v, " ", r, " ", dist_to_1st)
		push!(time_to_1st, [v, r, dist_to_1st])
	end
	time_to_1st
end

end
end
=#
