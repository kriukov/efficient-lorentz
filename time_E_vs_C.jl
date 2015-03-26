using EfficientLorentz

x = [0.532282, 0.531928, 0.267589]; v = [0.665341, 0.668661, 0.331986]


time_vs_r = Array{Real, 1}[]


for n = 1:7
	r0 = 1/10^n
	for i = 0:9
		r = r0 - i/10^(n+1)
		time_E = (tic(); collisions3d(x, v, r, 1, 64); toc())
		push!(time_vs_r, [r, time_E])
	end
	time_vs_r
end



