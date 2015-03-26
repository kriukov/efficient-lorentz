include("ClassicalLorentz.jl")
include("EfficientLorentz.jl")

x = [0.,0.445]
v = [cos(1), sin(1)]
r = 0.1
t = 30
precision = 64

function compare_lorentz2d(x, v, r, t, prec)
	elements_classical = collisions_classical(x, v, r, t, prec)[2]
	places_classical = collisions_classical(x, v, r, t, prec)[1]
	speeds_classical = collisions_classical(x, v, r, t, prec)[3]
	l = length(elements_classical)
	elements_efficient = collisions(x[1], x[2], v[1], v[2], r, l, prec)[1]
	places_efficient = collisions(x[1], x[2], v[1], v[2], r, l, prec)[2]
	speeds_efficient = collisions(x[1], x[2], v[1], v[2], r, l, prec)[3]
	div_place = 0
	div_element = 0
	for i = 1:l
		if elements_classical[i] == elements_efficient[i]
			println("Element $i passed: $(elements_efficient[i])")
			println("Place: \n efficient  $(places_efficient[i+1]), \n classical  $(places_classical[i+1]), \n difference $(places_efficient[i+1] - places_classical[i+1])")
			println("Speed: \n efficient  $(speeds_efficient[i+1]), \n classical  $(speeds_classical[i+1]), \n difference $(speeds_efficient[i+1] - speeds_classical[i+1])")
			println("---------------------------")
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
