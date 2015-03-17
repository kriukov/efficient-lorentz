			
			
			
			
			
			
			
				# Find the coordinates where the straight line ends on one 2d circle, let's say xy, and continue from a little further from there - wrong idea: misses balls
				# We need to find nearest place and continue from there
				place_xy = collisions(x[1], x[2], v_xy[1], v_xy[2], r, 1, precision)[2][2]
				place_yz = collisions(x[2], x[3], v_yz[1], v_yz[2], r, 1, precision)[2][2]
				place_xz = collisions(x[1], x[3], v_xz[1], v_xz[2], r, 1, precision)[2][2]
				#place_z = x[3] + v[3]/v[2]*(place_xy[2] - x[2])
				
				l1 = norm(place_xy - [x[1], x[2]])
				l2 = norm(place_yz - [x[2], x[3]])
				l3 = norm(place_xz - [x[1], x[3]])
				
				l = max(l1, l2, l3)
				
				# Need to advance from the false collision point at least by a distance larger than a ball diameter, otherwise the algorithm may get stuck at the same point
				
				if l == l1
					place_z = x[3] + v[3]/v[2]*(place_xy[2] - x[2])
					x = [place_xy[1], place_xy[2], place_z]	+ v*(2r + 0.1) 
				elseif l == l2
					place_x = x[1] + v[1]/v[2]*(place_xy[2] - x[2])
					x = [place_x, place_yz[1], place_yz[2]]	+ v*(2r + 0.1) 
				elseif l == l3
					place_y = x[2] + v[2]/v[1]*(place_xy[1] - x[1])
					x = [place_xz[1], place_y, place_xz[2]]	+ v*(2r + 0.1) 
				end





		
		if possible_hit 
		
			if suspected_times(times1)
				ball = [x1, y1, z3]
				x_new = collide3d(x, ball, v, r)
				if x_new == false
					continue
				else
					# Definitely a collision
					v = v_new(x_new, ball, v)
					x = x_new
					push!(places, ball)
					push!(coords, x)
					push!(speeds, v)
					steps += 1
				end
			elseif suspected_times(times2)
				ball = [x1, y1, z2]
				x_new = collide3d(x, ball, v, r)
				if x_new == false
					continue
				else
					# Definitely a collision
					v = v_new(x_new, ball, v)
					x = x_new
					push!(places, ball)
					push!(coords, x)
					push!(speeds, v)
					steps += 1
				end
			elseif suspected_times(times3)
				ball = [x3, y2, z2]
				x_new = collide3d(x, ball, v, r)
				if x_new == false
					continue
				else
					# Definitely a collision
					v = v_new(x_new, ball, v)
					x = x_new
					push!(places, ball)
					push!(coords, x)
					push!(speeds, v)
					steps += 1
				end
			else
				advance(x, v, x1, y1, y2, z2, x3, z3)
			end

			
		# And if possible_hit returns false, i.e., none of the 3 conditions is satisfied, we take the farthest (!!! nearest - see timeline) point it went and continue from there
		else
			advance(x, v, x1, y1, y2, z2, x3, z3)
		end	
