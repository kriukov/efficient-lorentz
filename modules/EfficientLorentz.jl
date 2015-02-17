module EfficientLorentz
export frac, efficient_algorithm, first_collision, collide, post_collision, dist_point_line, dist_point_line_sign, collisions, collisions2

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
dist_point_line(x, y, k, b) = abs(y - k*x - b)/sqrt(k^2 + b^2)
dist_point_line_sign(x, y, k, b) = (y - k*x - b)/sqrt(k^2 + b^2)

function collisions2(x, y, vx, vy, r, maxsteps)
	steps = 1
	places = Vector[]
	# Push a dummy vector to "places" - it cannot be empty for array_corners
	push!(places, [-Inf, -Inf])
	
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
		
		# Array of the rest of the corners (3 other which may or may not experience collision)
		# The first two  numbers are the corner coords, and the third is true/false (1/0) whether it leaves square
		array_rest_corners = Array{Int, 1}[]
		for i = 1:4
			if i != j
				push!(array_rest_corners, [array_corners[i], d(i) >= r])
			end
		end
		
		leaves_square = array_rest_corners[1][3] == 1 && array_rest_corners[2][3] == 1 && array_rest_corners[3][3] == 1

		# If the particle exits the unit square without another collision
		if leaves_square
	   
	   		# determine through which wall it will exit
	   		# times to each wall (vertical, horizontal)
	   		tv1 = (n - x)/vx 
			tv2 = (n - x + 1)/vx
			th1 = (m - y)/vy 
			th2 = (m - y + 1)/vy

			array_times = Real[]
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
				@show push!(places, [q + posx, p + posy])
				@show x, y, vx, vy = collide(q + posx, p + posy, x1 + posx, y1 + posy, vx, vy, r) # And we obtain coords of collision and new velocity, and then cycle continues
				
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
				@show push!(places, [p + posy, q + posx])
				# Find the next point of collision
				@show x, y, vx, vy = collide(p + posy, q + posx, y1 + posy, x1 + posx, vx, vy, r)
		   	
		   	end
		   	
		# Now what if the particle doesn't exit the square? Obtain where it doesn't (coords of obstacle), collide there and continue cycle   	
		elseif array_rest_corners[1][3] == 0
			place = [array_rest_corners[1][1], array_rest_corners[1][2]]
			@show push!(places, place)
			@show x, y, vx, vy = collide(place[1], place[2], x, y, vx, vy, r)
		elseif array_rest_corners[2][3] == 0
			place = [array_rest_corners[2][1], array_rest_corners[2][2]]
			@show push!(places, place)		
			@show x, y, vx, vy = collide(place[1], place[2], x, y, vx, vy, r)
		elseif array_rest_corners[3][3] == 0
			place = [array_rest_corners[3][1], array_rest_corners[3][2]]
			@show push!(places, place)
			@show x, y, vx, vy = collide(place[1], place[2], x, y, vx, vy, r)
		end
		steps += 1
	end
	deleteat!(places, 1)
	return places
end






















#-> Takes the point inside a square formed by the four closest obstacles and determines which wall will be crossed, thus converting the conditions back to x0 = (0, b) to use with first_collision() (with possible turning the coordinates)
function post_collision(x, y, vx, vy, r)
	n = ifloor(x)
	m = ifloor(y)
	
	k = vy/vx
	b = k*x + y
	
	if dist_point_line(n, m, k, b) >= r &&
	   dist_point_line(n, m+1, k, b) >= r &&
	   dist_point_line(n+1, m, k, b) >= r &&
	   dist_point_line(n+1, m+1, k, b) >= r	


		tv1 = (n - x)/vx 
		tv2 = (n - x + 1)/vx
		th1 = (m - y)/vy 
		th2 = (m - y + 1)/vy

		array_times = Real[]
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
		# number = 1 or 2 - vertical crossing, 3 or 4 - horizontal

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
		elseif number == 3
			posx = n
			posy = m
			x1 = x + (m - y)/k - n
			y1 = 0
		elseif number == 4
			posx = n
			posy = m + 1
			x1 = x + (m + 1 - y)/k - n
			y1 = 0
		end

		return x1, y1, posx, posy, number
	# It returns the coordinates (posx, posy) of the bottom-left corner of the unit square and coordinates (x1, y1) within this unit square
	# The exact coordinates of the particle at crossing are posx+x1, posy+y1

	# What if one of the distances to the four obstacles is less than r? Collide and reflect, then recursively apply the function until it goes out of the square
	#=
	elseif dist_point_line(n, m, k, b) < r
		println("1 more collision, case $n $m")
		x1, y1, vx1, vy1 = collide(n, m, x, y, vx, vy, r)
		post_collision(x1, y1, vx1, vy1, r)
	elseif dist_point_line(n, m+1, k, b) < r
		println("1 more collision, case $n, $(m+1)")		
		x1, y1, vx1, vy1 = collide(n, m+1, x, y, vx, vy, r)
		post_collision(x1, y1, vx1, vy1, r)
	elseif dist_point_line(n+1, m, k, b) < r
		println("1 more collision, case $(n+1), $m")		
		x1, y1, vx1, vy1 = collide(n+1, m, x, y, vx, vy, r)
		post_collision(x1, y1, vx1, vy1, r)
	elseif dist_point_line(n+1, m+1, k, b) < r
		println("1 more collision, case $(n+1), $(m+1)")
		x1, y1, vx1, vy1 = collide(n+1, m+1, x, y, vx, vy, r)
		post_collision(x1, y1, vx1, vy1, r)
		=#
	end

end

#-> Calculates the coordinates of the obstacles for a given number of collisions (trajectory)
function collisions(x, y, vx, vy, r, maxsteps)
	steps = 1
	places = Vector[]
	q, p = first_collision(x, y, vx, vy, abs(r/vx))
	push!(places, [q, p])
	#b = y
	posx = 0; posy = 0
	while steps <= maxsteps
		println("step ", steps)
	
		x1, y1, vx1, vy1 = collide(q, p, x, y, vx, vy, r) # And we obtain coords of collision and new velocity
		println("Coords of collision and new speed: x1 = $x1, y1 = $y1; posx = $posx, posy = $posy; vx1 = $vx1, vy1 = $vy1")
		
		# What if no exit?
		k = vy1/vx1
		b = k*x1 + y1
		
		if dist_point_line(q + posx, p + posy, k, b) >= r &&
	   dist_point_line(q + posx, p+posy+1, k, b) >= r &&
	   dist_point_line(q+posx+1,p+posy, k, b) >= r &&
	   dist_point_line(q+posx+1, p+posy+1, k, b) >= r	
			x2, y2, posx, posy, number = post_collision(x1 + posx, y1 + posy, vx1, vy1, r)
			println("Exiting square: x2 = $x2, y2 = $y2, posx = $posx, posy = $posy, wall = $number")
		end
		
		
		
		if number == 1 || number == 2 # Left/right wall crossed
			q, p = first_collision(x2, y2, vx1, vy1, abs(r/vx1))
			@show push!(places, [q + posx, p + posy])
		elseif number == 3 || number == 4 # Bottom/top wall crossed
			# Rotate the axes, x->y, y->x and then interchange q and p
			q, p = first_collision(y2, x2, vy1, vx1, abs(r/vy1))
			#q = q1
			#p = p1
			@show push!(places, [p + posx, q + posy])
			y = x2
			x = y2
			vy = vx1
			vx = vy1
		end
			
		@show steps += 1
	end
	places
end

# End of module
end
