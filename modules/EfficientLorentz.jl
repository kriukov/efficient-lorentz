module EfficientLorentz
export frac, efficient_algorithm, first_collision, collide, post_collision, dist_point_line, dist_point_line_sign, collisions, collisions3d

function frac(alpha, epsilon)
    p = 0; q = 0
    h1 = 0
    h2 = 1
    k1 = 1
    k2 = 0
    mm1 = alpha
    a1 = ifloor(mm1)
    while abs(k2*alpha - h2) > epsilon
        mm2 = 1/(mm1 - a1)
        a2 = ifloor(mm2)
        h3 = a1*h2 + h1
        k3 = a1*k2 + k1
        p = h3
        q = k3
        h1 = h2
        h2 = h3
        a1 = a2
        k1 = k2
        k2 = k3
        mm1 = mm2
    end
    q, p
end 

function efficient_algorithm(m, b, epsilon)
	kn = 0
	#i = 0

	while b > epsilon && 1 - b > epsilon
		if b < 0.5
			(q, p) = frac(m, 2b)
		else
			(q, p) = frac(m, 2*(1 - b))
		end
		b = mod(m*q + b, 1)
		kn += q
		#i += 1
		#if i >= 3000
			#println("cycle break: too many iterations")
			#break
		#end
	end
	#println("cycle ended")
	q = kn
	p = ifloor(m*q) + 1
	return (q, p)	
	#return i
end

#-> Finds the coordinates of the center of the first-collision obstacle
function first_collision(x, y, vx, vy, delta)

	# Normalize velocity if it wasn't normalized
	v = sqrt(vx^2 + vy^2)
	vx1 = vx/v
	vy1 = vy/v
	vx = vx1
	vy = vy1
		
	m = vy/vx
	b = y - m*x
	
	if b > delta && 1 - b > delta
	
		if vx > 0 && vy > 0
			(q, p) = efficient_algorithm(m, b, delta)
			p = ifloor(m*q) + 1
		elseif vx < 0 && vy > 0
			m = -m
			(q, p) = efficient_algorithm(m, b, delta)
			p = ifloor(m*q) + 1
			m = -m
			q = -q
		elseif vx < 0 && vy < 0
			b = 1 - b
			(q, p) = efficient_algorithm(m, b, delta)
			b = 1 - b
			p = -ifloor(m*q) #p = -(ifloor(m*q) + 1)
			q = -q
		elseif vx > 0 && vy < 0
			b = 1 - b
			m = -m
			(q, p) = efficient_algorithm(m, b, delta)
			b = 1 - b
			#m = -m # this line should not be here
			p = -ifloor(m*q) #p = -(ifloor(m*q) + 1)
		end
		
	# Added code for cases when b or 1-b < delta
	# The problem with efficient_algorithm() is that whenever b<delta or 1-b<delta, it always outputs (0, 1) - see plot b_or_1-b_lt_delta.png: even though the speed is negative in both directions (starting point is (0, b)) and the first collision is (-2, -2), the algorithm outputs (0, 1), which is in the other direction.
	# We have to make special cases to fix it
	
	elseif b <= delta
		# Four situations possible: leaves through top, leaves through right/left, hits (1, 1)/(-1, 1), hits (0, 1)
		r = abs(delta*vx)
		# Three critical values of m, ascending
		m1 = (b - 1 + r*sqrt(2 + b*(b - 2) - r^2))/(r^2 - 1)
		m2 = (r - (b - 1)*sqrt(2 + b*(b - 2) - r^2))/((b - 1)*r + sqrt(2 + b*(b - 2) - r^2))
		m3 = sqrt(((b - 1)/r)^2 - 1)
		
		# Leaves through top
		if abs(m) < m3 && abs(m) > m2
			if vx > 0 && vy > 0
				(q, p) = efficient_algorithm(1/m, (1-b)/m, delta/m)
				q, p = p, q
				p += 1
			elseif vx < 0 && vy > 0
				m = -m
				(q, p) = efficient_algorithm(1/m, (1-b)/m, delta/m)
				q, p = p, q
				p += 1
				q = -q
			elseif vy < 0
				q = 0
				p = 0
			end
		elseif abs(m) > m3
			q = 0
			p = 1
		elseif m < m2 && m > m1
			q, p = 1, 1
		elseif m > -m2 && m < -m1
			q, p = -1, 1
		elseif abs(m) < m1
			if vx > 0 && vy > 0
				(q, p) = efficient_algorithm(m, m + b, delta)
				q += 1
			elseif vx < 0 && vy > 0
				m = -m
				(q, p) = efficient_algorithm(m, m + b, delta)
				q += 1
				q = -q
			elseif vy < 0
				q = 0
				p = 0
			end		
		end
		
	elseif 1 - b <= delta
		
		# It is just transformation y -> 1 - y, vy -> -vy and the previous situation applies
		
		m = -m
		b = 1 - b
		vy = -vy
		
		r = abs(delta*vx)
		# Three critical values of m, ascending
		m1 = (b - 1 + r*sqrt(2 + b*(b - 2) - r^2))/(r^2 - 1)
		m2 = (r - (b - 1)*sqrt(2 + b*(b - 2) - r^2))/((b - 1)*r + sqrt(2 + b*(b - 2) - r^2))
		m3 = sqrt(((b - 1)/r)^2 - 1)
		
		# Leaves through top
		if abs(m) < m3 && abs(m) > m2
			if vx > 0 && vy > 0
				(q, p) = efficient_algorithm(1/m, (1-b)/m, delta/m)
				q, p = p, q
				p += 1
			elseif vx < 0 && vy > 0
				m = -m
				(q, p) = efficient_algorithm(1/m, (1-b)/m, delta/m)
				q, p = p, q
				p += 1
				q = -q
			elseif vy < 0
				q = 0
				p = 0
			end
		elseif abs(m) > m3
			q = 0
			p = 1
		elseif m < m2 && m > m1
			q, p = 1, 1
		elseif m > -m2 && m < -m1
			q, p = -1, 1
		elseif abs(m) < m1
			if vx > 0 && vy > 0
				(q, p) = efficient_algorithm(m, m + b, delta)
				q += 1
			elseif vx < 0 && vy > 0
				m = -m
				(q, p) = efficient_algorithm(m, m + b, delta)
				q += 1
				q = -q
			elseif vy < 0
				q = 0
				p = 0
			end		
		end		
		
		p = 1 - p
		
	end	
	return q, p
end

#-> Gives the new velocity and the new coordinates (which are a point of collision) having the initial position/velocity, the obstacle coordinates and radius (this one is based on the classic algorithm)
# Tested: OK
function collide(q, p, x, y, vx, vy, r)
	r0 = [x, y]
	v0 = [vx, vy]
	R = [q, p]
	crossz(x, y) = x[1]*y[2] - x[2]*y[1]
	discr = norm(v0)^2*r^2 - (crossz(v0, r0-R))^2
	t1 = (-dot(v0, r0-R) - sqrt(discr))/norm(v0)^2
	N0 = r0 + v0*t1 - R
	N = N0/norm(N0)			
	v1 = v0 - 2*dot(v0, N)*N
	r1 = r0 + v0*t1
	return r1[1], r1[2], v1[1], v1[2]
end

# Distance from point (x, y) to line y = kx + b
dist_point_line(x, y, k, b) = abs(y - k*x - b)/sqrt(k^2 + 1)
dist_point_line_sign(x, y, k, b) = (y - k*x - b)/sqrt(k^2 + 1)

#-> Finds the trajectory
function collisions(x, y, vx, vy, r, maxsteps, precision::Integer=64)
	# Normalize velocity if it wasn't normalized
	v = sqrt(vx^2 + vy^2)
	vx1 = vx/v
	vy1 = vy/v
	vx = vx1
	vy = vy1
	
	set_bigfloat_precision(precision)
	x = big(x); y = big(y); vx = big(vx); vy = big(vy)
	steps = 1
	places = Vector{BigInt}[]
	coords = Vector{BigFloat}[]
	speeds = Vector{BigFloat}[]
	# Push a dummy place to "places" - it cannot be empty for array_corners; also, it can't be near [x, y] and it can't be NaN or Inf
	push!(places, [-10^9, -10^9])
	
	push!(coords, [x, y])
	push!(speeds, [vx, vy])
	
	while steps <= maxsteps
	
		# Will the particle exit the square?
		n = ifloor(x)
		m = ifloor(y)
		k = vy/vx
		b = - k*x + y
		
		# If exits, not counting the obstacle that has just experienced collision (of course distance to it < r)
		# Make an array of corners and mark the corner where the collision just happened (which is the last item in "places" array)
		array_corners = Array{Int, 1}[]
		push!(array_corners, [n, m], [n, m+1], [n+1, m], [n+1, m+1])
		
		d(i) = dist_point_line(array_corners[i][1], array_corners[i][2], k, b)
		
		j = 0
		for i = 1:length(array_corners)
			if array_corners[i] == places[length(places)]
				j = i
			end
		end
		
		
		# Array of the rest of the corners (3 other which may or may not experience collision); does not work at 1st step
		# The first two numbers are the corner coords, and the third is true/false (1/0) whether it leaves square
		array_rest_corners = Array{Int, 1}[]
		for i = 1:4
			if i != j
				push!(array_rest_corners, [array_corners[i], d(i) >= r])
			end
		end
		
		# 1st step: delete the obstacle that is the closest backwards
		# Determine through which wall it would exit if moved backwards
   		# Times to each wall (vertical, horizontal)
   		if steps == 1	
	   		
	   		tv1 = (n - x)/vx 
			tv2 = (n - x + 1)/vx
			th1 = (m - y)/vy 
			th2 = (m - y + 1)/vy

			array_times_back = BigFloat[]
			push!(array_times_back, tv1, tv2, th1, th2)

			# Extract the maximum negative time value (closest wall backwards)
			maxneg = -Inf
			number = 0
			for i = 1:length(array_times_back)
				if array_times_back[i] > maxneg && array_times_back[i] < 0
					maxneg = array_times_back[i]
					number = i
				end				
			end
			# number: 1 = left, 2 = right, 3 = bottom, 4 = top
			if number == 1
				if k*n + b < m + 0.5
					deleteat!(array_rest_corners, 1)
				else
					deleteat!(array_rest_corners, 2)
				end
			elseif number == 2
				if k*(n + 1) + b < m + 0.5
					deleteat!(array_rest_corners, 3)
				else
					deleteat!(array_rest_corners, 4)
				end
			elseif number == 3
				if (m - b)/k < n + 0.5
					deleteat!(array_rest_corners, 1)
				else
					deleteat!(array_rest_corners, 3)
				end
			elseif number == 4
				if (m + 1 - b)/k < n + 0.5
					deleteat!(array_rest_corners, 2)
				else
					deleteat!(array_rest_corners, 4)
				end
			end
		end
		
		# Whether the particle leaves square or collides in the same square
		leaves_square = array_rest_corners[1][3] == 1 && array_rest_corners[2][3] == 1 && array_rest_corners[3][3] == 1



		# There is a difficulty: at the first step, we have a dummy place in places, which is obviously not the previous collision, and the corner which is behind the initial position of the particle doesn't get deleted and the algorithm may mistakenly count it as the first collision. We have to detect and delete it.
		# The algorithm below didn't help. There is another one above, at steps == 1
		#=
		if steps == 1 && !leaves_square
			hit_corners = Array{Int, 1}[]
			if length(array_rest_corners) == 4
				for i = 1:4
					if array_rest_corners[i][3] == 0
						push!(hit_corners, [array_rest_corners[i][1], array_rest_corners[i][2], i])
					end
				end
				# Time to each hit corner: if it is negative, it is behind and needs to be deleted
				for i = 1:length(hit_corners)
					t1 = dot([hit_corners[i][1], hit_corners[i][2]] - [x, y], [vx, vy])
					if t1 < 0
						deleteat!(array_rest_corners, hit_corners[i][3])
					end
				end
			
			end
					
		end
		=#


		# If the particle exits the unit square without another collision
		if leaves_square
		
	   		# determine through which wall it will exit
	   		# times to each wall (vertical, horizontal)
	   		tv1 = (n - x)/vx 
			tv2 = (n - x + 1)/vx
			th1 = (m - y)/vy 
			th2 = (m - y + 1)/vy

			array_times = BigFloat[]
			push!(array_times, tv1, tv2, th1, th2)

			# Extract the minimum positive time value
			minpos = Inf
			number = 0
			for i = 1:length(array_times)
				if array_times[i] < minpos && array_times[i] > 0
					minpos = array_times[i]
					number = i
				end				
			end
			# number: 1 = left, 2 = right, 3 = bottom, 4 = top
			
		   	# if the walls are vertical, we don't need to rotate coords for first_collision(). We just run first_collision(), get the next obstacle and collide()
		   	if number == 1 || number == 2
		   	
		   		if number == 1
					posx = n
					posy = m
					x1 = 0
					y1 = y + k*(n - x) - m
				elseif number == 2
					posx = n + 1
					posy = m
					x1 = 0
					y1 = y + k*(n + 1 - x) - m
				end
				
				q, p = first_collision(x1, y1, vx, vy, abs(r/vx))
				push!(places, [q + posx, p + posy])
	
				x, y, vx, vy = collide(q + posx, p + posy, x1 + posx, y1 + posy, vx, vy, r)
				push!(coords, [x, y])
				push!(speeds, [vx, vy])
				# And we obtain coords of collision and new velocity, and then cycle continues
				
		   	elseif number == 3 || number == 4
		   		# Here we have to rotate
		   	
				if number == 3
					posx = m
					posy = n
					y1 = x + (m - y)/k - n
					x1 = 0
				elseif number == 4
					posx = m + 1
					posy = n
					y1 = x + (m + 1 - y)/k - n
					x1 = 0
				end
				
				# Now the coords are rotated. In this rotated system, find the coords of the next obstacle
				q, p = first_collision(x1, y1, vy, vx, abs(r/vy))
				# Record the coords rotated back
				push!(places, [p + posy, q + posx])
				# Find the next point of collision
				x, y, vx, vy = collide(p + posy, q + posx, y1 + posy, x1 + posx, vx, vy, r)
				push!(coords, [x, y])
				push!(speeds, [vx, vy])
		   	
		   	end
		   	
		# Now what if the particle doesn't exit the square? Obtain where it doesn't (coords of obstacle), collide there and continue cycle   	
		
		elseif array_rest_corners[1][3] == 0
			place = [array_rest_corners[1][1], array_rest_corners[1][2]]
			push!(places, place)
			x, y, vx, vy = collide(place[1], place[2], x, y, vx, vy, r)
			push!(coords, [x, y])
			push!(speeds, [vx, vy])
		elseif array_rest_corners[2][3] == 0
			place = [array_rest_corners[2][1], array_rest_corners[2][2]]
			push!(places, place)		
			x, y, vx, vy = collide(place[1], place[2], x, y, vx, vy, r)
			push!(coords, [x, y])
			push!(speeds, [vx, vy])
		elseif array_rest_corners[3][3] == 0
			place = [array_rest_corners[3][1], array_rest_corners[3][2]]
			push!(places, place)
			x, y, vx, vy = collide(place[1], place[2], x, y, vx, vy, r)
			push!(coords, [x, y])
			push!(speeds, [vx, vy])
		end
		steps += 1
	end
	# Delete the dummy place at position 1
	deleteat!(places, 1)
	return places, coords, speeds
end


## 3D version



# to calculate the place where a particle will collide, if the obstacle have center x2, and radius r, and the particle have velocity v and initial position x1
function collide3d(x1, x2, v, r)
    b = BigFloat(dot((x1 - x2), v))/norm(v)^2
    c = BigFloat(norm(x1 - x2)^2 - r^2)
    if BigFloat(b^2 - c) < BigFloat(0)   # if there is no collision, return false
    return false
    end
    t = -b - sqrt(BigFloat(b^2 - c))
    x = v*t + x1
    return x
end

# And, to calculate the velocity after the collision at the point x1, and x2 is the center of the sphere
function v_new(x1, x2, v)
	n = x1 - x2
	n = n/norm(n)
	vn = dot(n, v)*n
	v = v - 2vn
	v = v/norm(v)
	return v
end

function collisions3d(x, v, r, maxsteps)

	places = Vector{BigInt}[]
	coords = Vector{BigFloat}[]
	speeds = Vector{BigFloat}[]

	steps = 1
	
	while steps <= maxsteps

		# Projection functions that convert 3D vectors into 2D without dot and cross
		Pxy(x) = [x[1], x[2]]
		Pyz(x) = [x[2], x[3]]
		Pxz(x) = [x[1], x[3]]

		# Calculating first collision in all 3 planes
		# Normalize velocities first
		v_xy = [Pxy(v)[1], Pxy(v)[2]]
		v_yz = [Pyz(v)[1], Pyz(v)[2]]
		v_xz = [Pxz(v)[1], Pxz(v)[2]]
	
		v_xy = v_xy/norm(v_xy)
		v_yz = v_yz/norm(v_yz)
		v_xz = v_xz/norm(v_xz)
		
		# First collision in 2D (function first_collision() only works for a special set of initial conditions, that's why collisions() has so much code generalizing it)
		first(x, v, r, precision::Integer=64) = collisions(x[1], x[2], v[1], v[2], r, 1, precision)[1][1]

		@show q1, p1 = first(Pxy(x), v_xy, r) #x, y
		@show q2, p2 = first(Pyz(x), v_yz, r) #y, z
		@show q3, p3 = first(Pxz(x), v_xz, r) #x, z
	
		# Condition that all of them correspond to the same point (x, y, z) of collision
		hit = q1 == q3 && p1 == q2 && p2 == p3
		
		x_new = collide3d(x, [q1, p1, p2], v, r)
	
		# If the collision happens
		if @show hit && x_new != false
			@show v = v_new(x_new, x, v)
			@show x = x_new
			@show push!(places, [q1, p1, p2])
			@show push!(coords, x)
			@show push!(speeds, v)
		# If it doesn't happen
		elseif @show hit && x_new == false || !hit && x_new == false
			@show q1, x[1], v[1]
			@show t1 = (q1 - x[1])/v[1]
			@show t2 = (p1 - x[2])/v[2]
			@show t3 = (p2 - x[3])/v[3]
			t = max(t1, t2, t3)
			# When the coordinate is integer, sometimes the 2d version can't decide in which square it is. We need to extend a little further from the possible integer
			@show x += v*(t + 0.01)
			
		end
	
	@show steps += 1
	end
	return places, coords, speeds
end


# End of module
end

#= Testing: using EfficientLorentz; x = [0.2, 0.3, 0.2]; v = [cos(1), sin(1), sin(0.4)]; r = 0.1


Line that will cross [3, 4, 5] with r = 0.1:

using EfficientLorentz; x = [0, 0.445, 0.342]; v = [0.451207, 0.541449, 0.709398]; r = 0.1
collisions3d(x, v, r, 5)


x = [-2.125652728087777,-2.1057861003273635,-3.0]
v = [-0.4512069325775071,-0.5414489190929208,-0.7093978939968069]

x = [-37.06038678750248, -44.027513426667056, -57.92519059385232]; v = [0.4512069325775123, 0.5414489190929171, 0.7093978939968061]

2d version gives strange results too... Test:

using EfficientLorentz; x = [0.445, 0.342]; v = [0.60672, 0.794915]; r = 0.2
collisions(x[1], x[2], v[1], v[2], r, 3)

using EfficientLorentz; x = [0.0,0.445]; v = [0.640184,0.768222]; r = 0.2

Next test 3D

using EfficientLorentz; x = [-1.96460353861745473485,12.0000000000739345148,-19.089208334145857042]; v = [0.492197878294012166407,0.739345146123457909362,-0.459467086423560331943]; r = 0.2

This line really hits [0, 15, -21]

!!!!!!!!!!!!!!
Similar test with the start not from the inside of an obstacle:
using EfficientLorentz; x = [-1.81694, 12.2218, -19.227]; v = [0.492197878294012166407,0.739345146123457909362,-0.459467086423560331943]; r = 0.2
The algorithm gets stuck at xy, yz, xz = [0, 15], [15, -21], [-1, 20]

using EfficientLorentz; x = [-31.5646794874135777267, 10.0000000000071421636, -8.86136400925839916656]; v = [-0.995170604930726595226, 0.0714216332842911117977, 0.0673380826933461717346]; r = 0.2
collisions3d(x, v, r, 5)



=#






