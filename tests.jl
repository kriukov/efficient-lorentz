

# using EfficientLorentz; x = 0; y = 0.445; vx = cos(1); vy = sin(1); r = 0.1
# using EfficientLorentz; x = 0; y = 0.35; vx = cos(1); vy = sin(1); r = 0.1


#c = sqrt(1.001) # Irrationalizing coefficient
# Split the 1st quadrant into N = 30 equal angular parts
# For each n < N, slope is tg(pi n/(2N))

#N = 30

#=
for n = 1:N-1
	loops = 0
	for m = 1:N
		for b = 1:N
			x = eff(tan(pi*n/(2N))*c, b/N*c, n/N)
			if x == 3000
				loops += 1
				break
			end
		end
	end
	println(n/N, " ", loops/N^2)
end
=#

#=
for n = 1:N-1
	loops = 0
	@show n
	for m = 1:N
		for b = 1:N
			m = tan(pi*rand()/2)
			x = eff(m, rand(), n/N*sqrt(1 + m^2))
			if x == 300
				loops += 1
				break
			end
		end
	end
	println(n/N, " ", loops/N^2)
end
=#

# Testing

#=
r = 0.001
x = 0
y = 0.445
for deg = 1:89
	if deg != 45
		vx = cos(deg*pi/180)
		vy = sin(deg*pi/180)
		println(deg, " ", first_collision(x, y, vx, vy, r/vx))
	end
end
=#

# Measuring the times

#= This one still gives error due to ContinuedFractions limits. Trying to work around
function measuring()
for i = 1:100
	r = i/1000
	for j = 1:100
		for k = 1:100
			vx = j*sqrt(0.999)
			vy = k*sqrt(1.002)
			#println(i, " ", j, " ", k, " ", vy/vx)
			first_collision(0, 0.445, vx, vy, r/vx)
		end
	end
end
end

@time measuring()
=#

#=
x = 0
y = 0.445


function measuring()
for i = 1:100
	r = i/1000
	for deg = 1:89
		if deg != 45
			vx = cos(deg*pi/180)
			vy = sin(deg*pi/180)
			#println(r, " ", deg, " ", first_collision(x, y, vx, vy, r/vx))
			first_collision(x, y, vx, vy, r/vx)
		end
	end
end
end

@time measuring()
=#

# with println: 112.495258954 s
# without println: 109.723781459 s

# After changing the frac() function into Atahualpa's one

# with println: 62.02397 s
# without println: 20.11622 s
# deg 1-44: 11.938 s
# deg 46-89: 7.825 s

#=
function measuring2()
	r = 0.000001
	for deg = 46:89
		if deg != 45
			vx = cos(deg*pi/180)
			vy = sin(deg*pi/180)
			#println(r, " ", deg, " ", first_collision(x, y, vx, vy, r/vx))
			first_collision(x, y, vx, vy, r/vx)
		end
	end
end

@time measuring2()
=#
