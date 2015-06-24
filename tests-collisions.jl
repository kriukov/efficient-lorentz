using ClassicalLorentz
using EfficientLorentz

#=
x = [0.01, 0.445]; v = [cos(1), sin(1)]; r = 0.1; t = 30; prec = 64
=#

function compare_lorentz2d(x, v, r, t, prec)
	classical = collisions_classical(x, v, r, t, prec)
	elements_classical = classical[2]
	places_classical = classical[1]
	speeds_classical = classical[3]
	l = length(elements_classical)
	
	efficient = collisions(x[1], x[2], v[1], v[2], r, l, prec)
	elements_efficient = efficient[1]
	places_efficient = efficient[2]
	speeds_efficient = efficient[3]
	div_place = 0
	div_element = 0
	for i = 1:l
		if elements_classical[i] == elements_efficient[i]
			#println("Element $i passed: $(elements_efficient[i])")
			#println("Place: \n efficient  $(places_efficient[i+1]), \n classical  $(places_classical[i+1]), \n difference $(places_efficient[i+1] - places_classical[i+1])")
			#println("Speed: \n efficient  $(speeds_efficient[i+1]), \n classical  $(speeds_classical[i+1]), \n difference $(speeds_efficient[i+1] - speeds_classical[i+1])")
			#println("---------------------------")
		else
			println("Diverged at element $i - increase precision")
			println("Classical: $(elements_classical[i]); Efficient: $(elements_efficient[i])")
			div_place = places_efficient[i+1] - places_classical[i+1]
			div_element = i
			break
		end
		
	end
	div_place, div_element
end


set_bigfloat_precision(512)
diverging_points = Array{BigFloat, 1}[]
for times = 0:1000
	deg = times*360/1000
	phi = deg*pi/180
	v = [cos(phi), sin(phi)]
	compare = compare_lorentz2d([0.01, 0.445], v, 0.1, 700, 512)
	println(deg)
	push!(diverging_points, [deg, norm(compare[1]), compare[2]])
end
println(diverging_points)

