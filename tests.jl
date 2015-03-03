

# using EfficientLorentz; x = 0; y = 0.445; vx = cos(1); vy = sin(1); r = 0.1
# using EfficientLorentz; x = 0; y = 0.35; vx = cos(1); vy = sin(1); r = 0.1


#c = sqrt(1.001) # Irrationalizing coefficient
# Split the 1st quadrant into N = 30 equal angular parts
# For each n < N, slope is tg(pi n/(2N))

#N = 30

#=
for n = 1:N-1
	loops = 0
	for m = 1:N
		for b = 1:N
			x = eff(tan(pi*n/(2N))*c, b/N*c, n/N)
			if x == 3000
				loops += 1
				break
			end
		end
	end
	println(n/N, " ", loops/N^2)
end
=#

#=
for n = 1:N-1
	loops = 0
	@show n
	for m = 1:N
		for b = 1:N
			m = tan(pi*rand()/2)
			x = eff(m, rand(), n/N*sqrt(1 + m^2))
			if x == 300
				loops += 1
				break
			end
		end
	end
	println(n/N, " ", loops/N^2)
end
=#

# Testing

#=
r = 0.001
x = 0
y = 0.445
for deg = 1:89
	if deg != 45
		vx = cos(deg*pi/180)
		vy = sin(deg*pi/180)
		println(deg, " ", first_collision(x, y, vx, vy, r/vx))
	end
end
=#

# Measuring the times

#= This one still gives error due to ContinuedFractions limits. Trying to work around
function measuring()
for i = 1:100
	r = i/1000
	for j = 1:100
		for k = 1:100
			vx = j*sqrt(0.999)
			vy = k*sqrt(1.002)
			#println(i, " ", j, " ", k, " ", vy/vx)
			first_collision(0, 0.445, vx, vy, r/vx)
		end
	end
end
end

@time measuring()
=#

#=
x = 0
y = 0.445


function measuring()
for i = 1:100
	r = i/1000
	for deg = 1:89
		if deg != 45
			vx = cos(deg*pi/180)
			vy = sin(deg*pi/180)
			#println(r, " ", deg, " ", first_collision(x, y, vx, vy, r/vx))
			first_collision(x, y, vx, vy, r/vx)
		end
	end
end
end

@time measuring()
=#

# with println: 112.495258954 s
# without println: 109.723781459 s

# After changing the frac() function into Atahualpa's one

# with println: 62.02397 s
# without println: 20.11622 s
# deg 1-44: 11.938 s
# deg 46-89: 7.825 s

#=
function measuring2()
	r = 0.000001
	for deg = 46:89
		if deg != 45
			vx = cos(deg*pi/180)
			vy = sin(deg*pi/180)
			#println(r, " ", deg, " ", first_collision(x, y, vx, vy, r/vx))
			first_collision(x, y, vx, vy, r/vx)
		end
	end
end

@time measuring2()
=#

############ Tests for 3D version


# Testing: using EfficientLorentz; x = [0.2, 0.3, 0.2]; v = [cos(1), sin(1), sin(0.4)]; r = 0.1


Line that will cross [3, 4, 5] with r = 0.1:

using EfficientLorentz; x = [0, 0.445, 0.342]; v = [0.451207, 0.541449, 0.709398]; r = 0.1
collisions3d(x, v, r, 5)


x = [-2.125652728087777,-2.1057861003273635,-3.0]
v = [-0.4512069325775071,-0.5414489190929208,-0.7093978939968069]

x = [-37.06038678750248, -44.027513426667056, -57.92519059385232]; v = [0.4512069325775123, 0.5414489190929171, 0.7093978939968061]

2d version gives strange results too... Test:

using EfficientLorentz; x = [0.445, 0.342]; v = [0.60672, 0.794915]; r = 0.2
collisions(x[1], x[2], v[1], v[2], r, 3)

using EfficientLorentz; x = [0.0,0.445]; v = [0.640184,0.768222]; r = 0.2

Next test 3D

using EfficientLorentz; x = [-1.96460353861745473485,12.0000000000739345148,-19.089208334145857042]; v = [0.492197878294012166407,0.739345146123457909362,-0.459467086423560331943]; r = 0.2

This line really hits [0, 15, -21]

!!!!!!!!!!!!!!
Similar test with the start not from the inside of an obstacle:
using EfficientLorentz; x = [-1.81694, 12.2218, -19.227]; v = [0.492197878294012166407,0.739345146123457909362,-0.459467086423560331943]; r = 0.2
The algorithm gets stuck at xy, yz, xz = [0, 15], [15, -21], [-1, 20]

using EfficientLorentz; x = [-31.5646794874135777267, 10.0000000000071421636, -8.86136400925839916656]; v = [-0.995170604930726595226, 0.0714216332842911117977, 0.0673380826933461717346]; r = 0.2
collisions3d(x, v, r, 5)

###############

using EfficientLorentz; x = [0, 0.445, 0.342]; v = [0.451207, 0.541449, 0.709398]; r = 0.2
collisions3d(x, v, r, 500)

Additional step at the end: 0.1
 BigInt[3,4,5]            
 BigInt[-5377,-1447,-1929]
 BigInt[-5386,-1453,-1927]
 BigInt[-5383,-1455,-1932]
 BigInt[-5432,-1490,-1909]
 BigInt[-5422,-1500,-1916]

Additional step at the end: 0.2
 BigInt[3,4,5]          
 BigInt[-709,-188,-251] 
 BigInt[-682,-320,-160] 
 BigInt[-685,-322,-158] 
 BigInt[-670,-310,-162] 
 BigInt[-759,-311,-187] 
 BigInt[-748,-287,-189] 
 BigInt[-739,-285,-189] 
 BigInt[-704,-286,-187] 
 BigInt[-634,-274,-230] 
 BigInt[-563,-246,-193] 
 BigInt[-580,-251,-214] 
 BigInt[-569,-233,-223] 
 BigInt[-468,-335,-167] 
 BigInt[-459,-326,-162] 


Additional step at the end: 0.5
 BigInt[3,4,5]          
 BigInt[-798,-212,-283] 
 BigInt[-798,-214,-282] 
 BigInt[-813,-182,-321] 
 BigInt[-768,-290,-320] 
 BigInt[-674,-254,-473] 
 BigInt[-307,-284,-529] 
 BigInt[-279,-323,-531] 
 BigInt[-280,-320,-528] 
 BigInt[-241,-341,-505] 
 BigInt[-139,-289,-486] 
 BigInt[-98,-311,-450]  
 BigInt[-57,-329,-428]  
 BigInt[-46,-300,-391]  
 BigInt[-95,-875,-623]  
 BigInt[-96,-877,-622]  
 BigInt[-135,-843,-640] 
 BigInt[-140,-857,-601] 
 BigInt[-287,-562,-163] 
 
No additional step
 BigInt[3,4,5]
 (none new even after 1000 iterations)
 
 
Probably this step is wrong


Testing again with r = 0.2, because 0.1 seems to fall in channels often

No step:
 BigInt[3,4,5]     
 BigInt[-14,-6,-8] 
 BigInt[-4,9,-17]  
 BigInt[-3,10,-20] 
 BigInt[-3,9,-19]  
 BigInt[-3,12,-18] 
 BigInt[-1,9,-19]  
 BigInt[-8,-11,47] 
 BigInt[-13,-15,44]
 BigInt[-20,-5,48] 
 BigInt[-26,-8,47] 
 BigInt[-20,-1,52] 
 BigInt[-19,-1,52] 
 BigInt[-19,-3,53] 
 BigInt[-20,-5,52] 

 
Step 0.1
 BigInt[3,4,5]      
 BigInt[-14,-6,-8]  
 BigInt[-4,9,-17]   
 BigInt[-3,10,-20]  
 BigInt[-2,11,-20]  
 BigInt[-3,10,-22]  
 BigInt[1,0,-11]    
 BigInt[0,-1,-17]   
 BigInt[-21,-17,-20]
 BigInt[-20,-14,-20]
 BigInt[-26,-12,-21]
 BigInt[6,-20,-101] 
 BigInt[5,-19,-102] 
 BigInt[-2,-20,-90] 
 BigInt[1,-20,-91]  
 BigInt[14,-2,-108] 

Step 0.3
 BigInt[3,4,5]       
 BigInt[-14,-6,-8]   
 BigInt[-4,9,-17]    
 BigInt[-3,10,-19]   
 BigInt[-32,7,-7]    
 BigInt[-27,16,-8]   
 BigInt[-28,16,-8]   
 BigInt[-69,18,-50]  
 BigInt[-74,12,-52]  
 BigInt[-74,12,-51]  
 BigInt[-73,21,-43]  
 BigInt[-78,15,-46]  
 BigInt[-77,18,-49]  
 BigInt[-79,15,-51]  
 BigInt[-70,16,-85]  
 BigInt[-71,12,-83]  
 BigInt[-77,6,-84]   
 BigInt[-81,5,-85]   
 BigInt[-85,4,-84]   
 BigInt[-84,3,-85]   
 BigInt[-79,11,-79]  

Maybe the variations are caused by chaoticity and not by this step?
Step is unnecessary

## ND tests
Example: N = 5

using EfficientLorentz; x = [0.213, 0.445, 0.331, 0.254, 0.188]; v = [0.347678, 0.541476, 0.233173, 0.414546, 0.599755]; N = 5; r = 0.1
collisionsnd(x, v, r, 5)



