include("eff-ata-latest.jl")
using EfficientLorentz
using ClassicalLorentz

# Ata's E vs Ata's C
#=
for j = 10:50
    for i = 1:10
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1, 0.8*rand() + 0.1]
	    phi = 2pi*rand(); theta = pi*rand()
		v = [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
        r = 0.8^j
        
        effic = int(Lorentz3D(x, v, r))
        class = int(ClassicLorentz(x, v, r))
        
        if effic == class
            #println("Passed with x = $x; v = $v; r = $r")
        else
            #my_effic = collisions3d_time(x, v, r, 1, 64)[1]
            println("Failed with x = $x; v = $v; r = $r : E = $effic, C = $class")
        end
    end
end
=#

# Ata's E vs my E
for j = 10:30
    for i = 1:10
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1, 0.8*rand() + 0.1]
	    phi = 2pi*rand(); theta = pi*rand()
		v = [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
        r = 0.8^j
        
        effic = int(Lorentz3D2(x, v, r))
        my_effic = collisions3d_time(x, v, r, 1, 64)[1][1]
        
        if effic == my_effic
            #println("Passed with x = $x; v = $v; r = $r")
        else
            t = norm(my_effic - x)/norm(v)
            my_class = collisions3d_classical(x, v, r, t + 0.1, 64)[2][1]
            println("Failed with x = $x; v = $v; r = $r : E = $effic, myE = $my_effic, myC = $my_class")
        end
    end
end
