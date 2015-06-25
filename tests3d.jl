include("efficient-ata.jl")

for j = 10:50
    for i = 1:10
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
	    phi = 2pi*rand()
	    v = [cos(phi), sin(phi)]
        r = 0.8^j
        if Lorentz(x, v, r) == int(ClassicLorentz(x, v, r))
            println("Passed with x = $x; v = $v; r = $r")
        else
            println("Failed with x = $x; v = $v; r = $r")
        end
    end
end
