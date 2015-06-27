#include("eff-ata-latest.jl")
exectimepercell_vs_r = Array{Real, 1}[]

#=
for i = 10:75
    r = 0.8^i
    
    # Average exectimepercell for a given r
    sum = 0
    N = 100
    for j = 1:N    
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
        phi = 2pi*rand()
        v = [cos(phi), sin(phi)]
        exectime = @elapsed Lorentz2(x, v, r)
	    Ncells = abs(x[1]) + abs(x[2])
	    sum += exectime/Ncells
	end
	exectimepercell = sum/N
	push!(exectimepercell_vs_r, [r, exectimepercell])
end

# The first exectime value is always unusually large, so discard it
deleteat!(exectimepercell_vs_r, 1)

println(exectimepercell_vs_r)
=#

# Classical algorithm

using ClassicalLorentz
for i = 10:40
    r = 0.8^i
    
    # Average exectimepercell for a given r
    sum = 0
    N = 10
    for j = 1:N    
        x = [0.8*rand() + 0.1, 0.8*rand() + 0.1]
        phi = 2pi*rand()
        v = [cos(phi), sin(phi)]
        exectime = @elapsed first_collision_classical(x, v, r, 64)
	    Ncells = abs(x[1]) + abs(x[2])
	    sum += exectime/Ncells
	end
	exectimepercell = sum/N
	push!(exectimepercell_vs_r, [r, exectimepercell])
end

# The first exectime value is always unusually large, so discard it
deleteat!(exectimepercell_vs_r, 1)

println(exectimepercell_vs_r)
