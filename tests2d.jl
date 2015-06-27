include("eff-ata-latest.jl")
using ClassicalLorentz
using EfficientLorentz

# Ata's E vs my E

set_bigfloat_precision(512)

for j = 10:70
    for i = 1:1000
        x = BigFloat[0.8*rand() + 0.1, 0.8*rand() + 0.1]
	    phi = 2pi*rand()
	    v = BigFloat[cos(phi), sin(phi)]
        r = 0.8^j
        
        effic = Lorentz2(x, v, r)
        class = int(ClassicLorentz(x, v, r))
        
        my_effic = collisions(x[1], x[2], v[1], v[2], r, 1, 512)[1][1]
                
        if effic == my_effic
            #println("Passed with x = $x; v = $v; r = $r")
        else
            t = norm(my_effic - x)/norm(v)
            my_class = collisions_classical(x, v, r, t + 0.1, 512)[2]
            println("Failed with x = $x; v = $v; r = $r : E = $effic, MyE = $my_effic, C = $class, MyC = $my_class")
        end
    end
end


# Ata's E vs Ata's C
#=
for j = 10:50
    for i = 1:100
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
	    phi = 2pi*rand()
	    v = [cos(phi), sin(phi)]
        r = 0.8^j
        
        effic = Lorentz2(x, v, r)
        class = int(ClassicLorentz(x, v, r))
        
        if effic == class
            #println("Passed with x = $x; v = $v; r = $r")
        else
            my_effic = collisions(x[1], x[2], v[1], v[2], r, 1, 64)[1][1]
            println("Failed with x = $x; v = $v; r = $r : E = $effic, C = $class, myE = $my_effic")
        end
    end
end
=#