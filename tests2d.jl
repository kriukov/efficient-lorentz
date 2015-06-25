include("efficient-ata.jl")
using ClassicalLorentz
using EfficientLorentz

for j = 10:50
    for i = 1:100
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
	    phi = 2pi*rand()
	    v = [cos(phi), sin(phi)]
        r = 0.8^j
        
        effic = Lorentz(x, v, r)
        class = int(ClassicLorentz(x, v, r))
        
        my_effic = collisions(x[1], x[2], v[1], v[2], r, 1, 64)[1][1]
        
        
        
        if effic == my_effic
            #println("Passed with x = $x; v = $v; r = $r")
        else
            t = norm(my_effic - x)/norm(v)
            my_class = collisions_classical(x, v, r, t + 0.1, 64)[2]
            println("Failed with x = $x; v = $v; r = $r : E = $effic, MyE = $my_effic, C = $class, MyC = $my_class")
        end
    end
end
