using EfficientLorentz
#using ClassicalLorentz

#x = [0.532282, 0.531928, 0.267589]; v = [0.665341, 0.668661, 0.331986] - bad velocity (too symmetric). Below is another.
x = [0.532282, 0.531928, 0.267589]; v = [cos(1)*sin(0.8), sin(1)*sin(0.8), cos(0.8)]

#= Measuring first collision time as a function of radius
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

#= Measuring first collision time (classical) as a function of radius
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

# Trying out different velocities, at least in the first octant
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
