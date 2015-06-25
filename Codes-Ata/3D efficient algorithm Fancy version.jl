
function collision(x1,x2,r,v)
    v2=dot(v,v)
    b=dot((x1-x2),v)/v2
    c=(dot((x1-x2),(x1-x2))-r^2)/v2
    if(BigFloat(b^2-c)<BigFloat("0"))
        return false, 0
    end
    t=-b-sqrt(b^2-c)
    x=v*t+x1
    return x, t

end

function velo_col(x1,x2,v)
    n=(x1-x2)
    n=n/norm(n)
    vn=dot(n,v)*n
    v=v-2*vn
    v=v/norm(v)
    return v
end

function frac(x, epsilon)
   h1, h2 = 1, 0
   k1, k2 = 0, 1
   b = x
   while abs(k1*x - h1) > epsilon
       a = ifloor(b)
       h1, h2 = a*h1 + h2, h1
       k1, k2 = a*k1 + k2, k1
       b = 1/(b - a)
   end
   return k1, h1
end



function nextt1(x,v,r) 
    e1=[1,0]
    e2=[0,1]
    nn=int(x)
    x1=x-nn      
    t=(-x1[1])/v[1]
    test=1
    if(t<0)
        t=(1-x1[1])/v[1] 
        test=2
    end
    xx1=x1+v*t
    b=xx1[2]    
    epsilon=r/v[1]
    testt=0
    if(test==1)
        if((b>-epsilon && b<epsilon))
            return nn, testt
        elseif((1-b>-epsilon) && (1-b<epsilon))
            return nn+e2, testt
        end       
    elseif(test==2)
        if((b>-epsilon && b<epsilon))
            return nn+e1, testt
        elseif(1-b>-epsilon && 1-b<epsilon) 
            return nn+e1+e2,testt
        end
    end
    testt=1
    return nn+xx1, testt
end

function eff(m, b, r)
	kn = 0
    b1=b
    epsilon=r*sqrt(m*m+1) 
    if(b < epsilon || 1 - b < epsilon)
        if b < BigFloat("0.5")
			(q, p) = frac(m, 2.*b)
		else
			(q, p) = frac(m, 2.*(1. - b))
		end
		b = mod(m*q + b, 1)
		kn += q
    end  
	while b > epsilon && 1 - b > epsilon
		if b < BigFloat("0.5")
			(q, p) = frac(m, 2.*b)
		else
			(q, p) = frac(m, 2.*(1. - b))
		end
		b = mod(m*q + b, 1)
		kn += q
	end
	q = kn
    p = int(m*q+b1)
    return [q, p]
end

function Lor(x,v,r)   #This function obtain the coordenates of the obstacle with which 
                    #a particle at a initial position x, and intial velocity v collids 
                    #if both components of the velocity are positive, and the slope is less than 1.
    x1,test=nextt1(x,v,r)
    if(test==0)
        return x1
    end
    v1=v[1]
    v2=v[2]
    m=v2/v1
    b=x1[2]
    b=b-floor(b)
    de=[0, int(b)]
    centro=eff(m,b,r)
    x2=int(x1)-de+centro
    return x2
end

function Lorentz1(x,v,r)  
    ROT=[[0 1],
        [-1 0]]  #Rotational matrix π/2 radians
    REF=[[0 1],
        [1 0]]   #Reflect matrix, change y->x and x->y
    ROT2=ROT^2
    ROT3=ROT^3
    v1=v[1]
    v2=v[2]
    vv=copy(v)
    xx=copy(x)
    m1=v2/v1
    if(norm(int(x)-x)<r)   #if a particle begin inside an obstacle, then the first collision 
                          #is considered with the same obstacle.
        return int(x)
    end
    if(m1>=0 && v2>=0)  # if the velocity is in the quadrant I 
        if(m1<=1)
            xx=Lor(xx,vv,r)
        elseif(m1>1)  
            xx=REF*xx
            vv=REF*vv
            xx=Lor(xx,vv,r)
            xx=REF*xx
            vv=REF*vv            
        end
        return xx
    elseif(m1>=0 && v2<0) #if the velocity is in the quadrant III
        xx=ROT2*xx
        vv=ROT2*vv
        if(m1<=1)
            xx=Lor(xx,vv,r)  
        elseif(m1>1)
            xx=REF*xx
            vv=REF*vv
            xx=Lor(xx,vv,r)
            xx=REF*xx
            vv=REF*vv   
        end
        xx=ROT2*xx
        vv=ROT2*vv
        return xx
    elseif(m1<0 && v2>=0) #if the velocity is in the quadrant II
        xx=ROT*xx
        vv=ROT*vv
        if(m1<-1)
            xx=Lor(xx,vv,r)            
        elseif(m1>=-1)
            xx=REF*xx
            vv=REF*vv
            xx=Lor(xx,vv,r)
            xx=REF*xx   
            vv=REF*vv 
        end     
        xx=ROT3*xx
        vv=ROT3*vv
        return xx
    elseif(m1<0 && v2<0) #if the velocity is in the quadrant IV
        xx=ROT3*xx
        vv=ROT3*vv
        if(m1<-1)
            xx=Lor(xx,vv,r)        
        elseif(m1>=-1)
            xx=REF*xx
            vv=REF*vv
            xx=Lor(xx,vv,r)
            xx=REF*xx
            vv=REF*vv 
        end   
        xx=ROT*xx
        vv=ROT3*vv
        return xx
    end
end

function LorentzGas1(x,v,r,steps)
    for i=1:steps
        center=Lorentz1(x,v,r) 
        x, t=collision(x,center,r,v)
        v=velo_col(x,center,v)
    end
    return x, v
end 

r=BigFloat("0.00001")
x=rand(2).*BigFloat[1,1]
angle=BigFloat(rand(1)[1])*pi*2.
v=[cos(angle), sin(angle)]
v=v/norm(v)
steps=1000
x,v=LorentzGas1(x,v,r,steps)


function LorentzGas2(x,v,r,time)
    t=0
    v1=copy(v)
    while(t<time)
        x1=copy(x)
        v1=copy(v)
        center=Lorentz1(x,v,r) 
        x, tt=collision(x,center,r,v)
        v=velo_col(x,center,v)
        t=t+tt #norm(x-x1)
    end
    t=t-time
    if(t==0)
        v1=v
    end
    x=x-t*v1
    return x, v1
end 

function Lorentz3D(x,v,r)
    if(norm(int(x)-x)<r)   #if a particle begin inside an obstacle, then the first collision 
                          #is considered with the same obstacle.
    xx, t=collision(x,int(x),r,v)
        return xx
    end
    x1=zeros(2)
    x2=zeros(2)
    x3=zeros(2)
    xr=false
    while (xr==false)
        v_xy=[v[1], v[2]]
        v_xz=[v[1], v[3]]
        v_yz=[v[2], v[3]]                #Velocities in the planes xy, xz, and yz
        vn_xy=v_xy./BigFloat(sqrt(dot(v_xy,v_xy)))
        vn_xz=v_xz./BigFloat(sqrt(dot(v_xz,v_xz)))
        vn_yz=v_yz./BigFloat(sqrt(dot(v_yz,v_yz))) #normalized velocities in the planes xy, xz, and yz
        tmax=BigFloat("0")
        tmin=BigFloat("-2")
        x1[1]=1
        x1[2]=2
        x2[1]=3
        x2[2]=4
        x3[1]=5 
        x3[2]=6
        #        while (int(tmax)>int(tmin+BigFloat("0.5")) )   #it can be considerin time conincidence
        while(x1[1]!=x2[1] || x1[2]!=x3[1] || x2[2]!=x3[2])    #or considering the obstacle coincidence.
            x_xy=[x[1], x[2]]
            x_xz=[x[1], x[3]]
            x_yz=[x[2], x[3]]                #positions in the planes xy, xz, and yz

            x1=Lorentz1(x_xy,vn_xy,r)#calculate the position of obstacle that collide with the particle 
            x2=Lorentz1(x_xz,vn_xz,r)# in the 2D Lorentz gas for the 3 planes
            x3=Lorentz1(x_yz,vn_yz,r) 
            t1=(sqrt(dot(x1-x_xy,x1-x_xy))-r)/sqrt(dot(v_xy,v_xy))
            t2=(sqrt(dot(x2-x_xz,x2-x_xz))-r)/sqrt(dot(v_xz,v_xz))
            t3=(sqrt(dot(x3-x_yz,x3-x_yz))-r)/sqrt(dot(v_yz,v_yz))  #calculate the time used to reach the obstacle. 

            tmax=maximum([t1 t2 t3])    # if the minumum and maximum valueds of time are almost the same, then 
            tmin=minimum([t1 t2 t3])    # there is very high probable a collision  
            x=x+v*(tmax-0.1*r)   

        end
        xr, tt=collision(x,int(x),r,v) #calculates the collision. In case that there is not, then, it advance the particle to the 
                                   #next cell
        if(tt<0)
            xr=false
        end
        x=x+v*(1-2*r)
    end
    return xr
end

function ClassicLorentz(x,v,r)    #This is a simpler and inefficient version of the Classic algorithm for ND. Just to 
                                  #check if the 3D version works correctly.
   if(norm(int(x)-x)<r)   #if a particle begin inside an obstacle, then the first collision 
                          #is considered with the same obstacle.
        xx,tt=collision(x,int(x),r,v)
        return xx
    end
    xr=false
    i=0.0
    while (xr==false)
        while (xr==false)
            xr, tt=collision(x,int(x+v*(i)),r,v)
            i=i+1-2*r
        end
        t=sqrt(dot(xr-x,xr-x)/dot(v,v))
        if(dot(x+t*v-xr,x+t*v-xr)>dot(xr-x,xr-x))
            xr=false
        end
    end
    return xr
end


