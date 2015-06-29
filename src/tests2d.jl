include("efficient-Lorentz.jl")
#include("eff-ata-latest-old.jl")
using ClassicalLorentz
#using EfficientLorentz

# Ata's E vs my E

#set_bigfloat_precision(512)
precision = get_bigfloat_precision()

srand(10)

#using Base.Test

#using ImmutableArrays
#typealias Vec2d Vector2{Float64}

#using Vector2d

#using Docile

typealias Vector2D Vector{BigFloat}

function compare()
    for j = 10:30
        @show j
        for i = 1:100
            x = big([0.8*rand() + 0.1, 0.8*rand() + 0.1])

            phi = big(2pi*rand())

            v = [cos(phi), sin(phi)]
            r = big(0.8)^j

            efficient = Lorentz2(x, v, r)
            #class = int(ClassicLorentz(x, v, r))

            #my_effic = collisions(x[1], x[2], v[1], v[2], r, 1, precision)[1][1]

            #if effic == my_effic
                #println("Passed with x = $x; v = $v; r = $r")
            #else
            t = norm(efficient - x) / norm(v)
            classical = collisions_classical(x, v, r, t+1e-10, precision)[2][1]
              #  println("Failed with x = $x; v = $v; r = $r : E = $effic, MyE = $my_effic, C = $class, MyC = $my_class")
            #@show (efficient, classical)
            @assert efficient == classical
            end
        end
    end


@time compare()






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
