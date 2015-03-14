using ClassicalLorentz
using EfficientLorentz

#=
x = [0, 0.445]; v = [cos(1), sin(1)]; r = 0.1; t = 30; precision = 64
=#

function compare_lorentz2d(x, v, r, t, precision)
	elements_classical = collisions_classical(x, v, r, t, precision)[2]
	l = length(elements_classical)
	elements_efficient = collisions(x[1], x[2], v[1], v[2], r, l, precision)[1]
	for i = 1:l
		if elements_classical[i] == elements_efficient[i]
			println("Element $i passed: $(elements_efficient[i])")
		else
			println("Diverged at element $i - increase precision")
			println("Classical: $(elements_classical[i]); Efficient: $(elements_efficient[i])")
			break
		end
		
	end
end


for deg = 0:359
	phi = deg*pi/180
	v = [cos(phi), sin(phi)]
	println(deg)
	compare_lorentz2d([0.01, 0.445], v, 0.1, 200, 2048)
end

