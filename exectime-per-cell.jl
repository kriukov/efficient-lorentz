include("eff-ata-latest.jl")
exectimepercell_vs_r = Array{Real, 1}[]

# E 2D
#=
for i = 10:75
    r = 0.8^i
    
    # Average exectimepercell for a given r
    sum = 0
    N = 1000
    for j = 1:N    
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
        phi = 2pi*rand()
        v = [cos(phi), sin(phi)]
        tic()
        x1 = Lorentz2(x, v, r)
        exectime = toq()
	    Ncells = abs(x1[1] - x[1]) + abs(x1[2] - x[2])
	    sum += exectime/Ncells
	end
	exectimepercell = sum/N
	push!(exectimepercell_vs_r, [r, exectimepercell])
end

# The first exectime value is always unusually large, so discard it
deleteat!(exectimepercell_vs_r, 1)
println(exectimepercell_vs_r)
=#

# C 2D

using ClassicalLorentz
for i = 40:60
    r = 0.8^i

    sum = 0
    N = 50
    for j = 1:N    
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
        phi = 2pi*rand()
        v = [cos(phi), sin(phi)]
        tic()
        x1 = first_collision_classical(x, v, r, 64)
        exectime = toq()
	    Ncells = abs(x1[1] - x[1]) + abs(x1[2] - x[2])
	    sum += exectime/Ncells
	end
	exectimepercell = sum/N
	push!(exectimepercell_vs_r, [r, exectimepercell])
end
deleteat!(exectimepercell_vs_r, 1)
println(exectimepercell_vs_r)


# 3D E
#=
for i = 10:50
    r = 0.8^i

    sum = 0
    N = 100
    for j = 1:N    
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1, 0.8*rand() + 0.1]
		phi = 2pi*rand(); theta = pi*rand()
		v = [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
        tic()
        x1 = Lorentz3D2(x, v, r)
        exectime = toq()
	    Ncells = abs(x1[1] - x[1]) + abs(x1[2] - x[2]) + abs(x1[3] - x[3])
	    sum += exectime/Ncells
	end
	exectimepercell = sum/N
	push!(exectimepercell_vs_r, [r, exectimepercell])
end

deleteat!(exectimepercell_vs_r, 1)

println(exectimepercell_vs_r)
=#

# 3D C
#=
using ClassicalLorentz
for i = 10:40
    r = 0.8^i
    sum = 0
    N = 10
    for j = 1:N    
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1, 0.8*rand() + 0.1]
		phi = 2pi*rand(); theta = pi*rand()
		v = [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
        tic()
        x1 = first_collision3d_classical(x, v, r, 64)
        exectime = toq()
	    Ncells = abs(x1[1] - x[1]) + abs(x1[2] - x[2]) + abs(x1[3] - x[3])
	    sum += exectime/Ncells
	end
	exectimepercell = sum/N
	@show push!(exectimepercell_vs_r, [r, exectimepercell])
end
deleteat!(exectimepercell_vs_r, 1)
println(exectimepercell_vs_r)
=#

