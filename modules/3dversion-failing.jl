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

		@show x1, y1 = first(Pxy(x), v_xy, r) #x, y
		@show y2, z2 = first(Pyz(x), v_yz, r) #y, z
		@show x3, z3 = first(Pxz(x), v_xz, r) #x, z
	
		# Condition that all of them correspond to the same point (x, y, z) of collision
		#hit = q1 == q3 && p1 == q2 && p2 == p3
		
		condition1 = x1 == x3
		condition2 = y1 == y2
		condition3 = z2 == z3
		
		# Possible hit
		possible_hit = condition1 || condition2 || condition3
		
		if condition1
			x0 = [x1, y1, z3]
		elseif condition2
			x0 = [x1, y1, z2]
		elseif condition3
			x0 = [x3, y2, z2]
		#else
			
			
			#=
			t11 = (x1 - x[1])/v[1]
			t12 = (x3 - x[1])/v[1]
			t21 = (y1 - x[2])/v[2]
			t22 = (y2 - x[2])/v[2]
			t31 = (z2 - x[3])/v[3]
			t32 = (z3 - x[3])/v[3]
			@show tmax = max(t11, t12, t21, t22, t31, t32)
			@show x_new = x + v*tmax
			=#
		end
		
		
		if possible_hit
			x_new = collide3d(x, x0, v, r)
			if x_new != false
				@show v = v_new(x_new, x0, v)
				@show x = x_new
				@show push!(places, x0)
				@show push!(coords, x)
				@show push!(speeds, v)
			end
		
		else
			newplacexy(x, v, r, precision::Integer=64) = collisions(x[1], x[2], v[1], v[2], r, 1, precision)[2][1]
			x2, y2 = newplacexy(Pxy(x), v_xy, r)
			z2 = x[3] + v[3]/v[2]*(y2 - x[2])
			x = [x2, y2, z2] + 0.001v
		end
		
		# If the collision happens
		#=
		if @show hit && x_new != false

		elseif @show !hit && x_new == false || !hit && x_new == false
			@show x1, x[1], v[1]
			@show t1 = (x0[1] - x[1])/v[1]
			@show t2 = (x0[2] - x[2])/v[2]
			@show t3 = (x0[3] - x[3])/v[3]
			t = max(t1, t2, t3)
			@show x += v*(t + 0.1)
			
		end
		=#
		
	@show steps += 1
	end
	return places, coords, speeds
end
