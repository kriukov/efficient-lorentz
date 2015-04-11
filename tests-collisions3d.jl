using ClassicalLorentz
using EfficientLorentz

diverding_points = Real[]

function compare_lorentz3d(x, v, r, t, precision)
	classical = collisions3d_classical(x, v, r, t, precision)
	elements_classical = classical[2]
	places_classical = classical[1]
	speeds_classical = classical[3]
	l = length(places_classical)
	
	efficient = collisions3d_time(x, v, r, l, precision)
	elements_efficient = efficient[1]
	places_efficient = efficient[2]
	speeds_efficient = efficient[3]
	
	div_place = 0
	
	for i = 1:l
		if elements_classical[i] == elements_efficient[i]
			println("Element $i passed: $(elements_efficient[i])")
			println("Place: \n efficient  $(places_efficient[i]), \n classical  $(places_classical[i+1]), \n difference $(places_efficient[i] - places_classical[i+1])")
			println("Speed: \n efficient  $(speeds_efficient[i]), \n classical  $(speeds_classical[i+1]), \n difference $(speeds_efficient[i] - speeds_classical[i+1])")
			println("---------------------------")
		else
			println("Diverged at element $i")
			println("Classical: $(elements_classical[i]); Efficient: $(elements_efficient[i])")
			div_place = places_efficient[i-1] - places_classical[i]
			break
		end
		
	end
	div_place
end

#=
x = [0.532282, 0.531928, 0.267589]; v = [0.665341, 0.668661, 0.331986]; r = 0.1
compare_lorentz3d(x, v, r, 1000, 1024)
=#

#=
diverging_points = Array{BigFloat, 1}[]
x = [0.532282, 0.531928, 0.267589]; v = [0.665341, 0.668661, 0.331986]; r = 0.05
for deg = 1:10
	for deg1 = 1:180
		phi = deg*pi/180
		theta = deg1*pi/180
		v = [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
		println("$deg, $deg1")
		@show push!(diverging_points, [deg, deg1, norm(compare_lorentz3d(x, v, r, 20000, 256))])
	end
end

println(diverding_points)
=#


