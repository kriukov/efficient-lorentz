using EfficientLorentz
#using ClassicalLorentz

#x = [0.532282, 0.531928, 0.267589]; v = [0.665341, 0.668661, 0.331986]
x = [0.532282, 0.531928, 0.267589]; v = [cos(1)*sin(0.8), sin(1)*sin(0.8), cos(0.8)]

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


#=
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
