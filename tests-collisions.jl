using ClassicalLorentz
using EfficientLorentz

#=
x = [0.,0.445]; v = [cos(1), sin(1)]; r = 0.1
t = 30
precision = 64
=#

function compare_lorentz2d(x, v, r, t, precision)
	places_classical = collisions_classical(x, v, r, t, precision)[2]
	l = length(places_classical)
	places_efficient = collisions(x[1], x[2], v[1], v[2], r, l, precision)[1]
	for i = 1:l
		if places_classical[i] == places_efficient[i]
			println("Element $i passed")
		else
			println("Diverged at element $i - increase precision")
			break
		end
		
	end
end
