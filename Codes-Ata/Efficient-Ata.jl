function frac(alpha, epsilon)
    p = 0; q = 0
	h1 = 0
	h2 = 1
	k1 = 1
	k2 = 0
	mm1 = BigFloat(alpha)
	a1 = ifloor(mm1)
    while abs(k2*alpha-h2) > epsilon
        mm2 = BigFloat("1")/(mm1 - a1)
        a2 = ifloor(mm2)
		h3 = a1*h2 + h1
		k3 = a1*k2 + k1
		p = h3
		q = k3
        h1=h2
        h2=h3
        a1=a2
        k1=k2
        k2=k3
        mm1=mm2
    end
    return q,p
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

function collision(x1,x2,r,v)
    v2=dot(v,v)
    b=dot((x1-x2),v)/v2
    c=(dot((x1-x2),(x1-x2))-r^2)/v2
    if(BigFloat(b^2-c)<BigFloat("0"))
        return false
    end
    t=-b-sqrt(b^2-c)
    x=v*t+x1
    return x

end

function velo_col(x1,x2,v)
    n=(x1-x2)
    n=n/sqrt(dot(n,n))
    vn=dot(n,v)*n
    v=v-2*vn
    v=v/sqrt(dot(v,v))
    return v
end

function nextt(x,v,r) 
    e1=[1,0]
    e2=[0,1]
    
    nn=int(x)
    x1=x-nn      
    t=(-dot(x1,e1))/dot(v,e1)
    test=1
    if(t<0)
        t=(1-dot(x1,e1))/dot(v,e1) #o =(1,b)
        test=2
    end
    xx1=x1+v*t
    b=dot(xx1,e2)    
    epsilon=r/dot(v,e1)
    if(test==1)
        if((b>-epsilon && b<epsilon))
            return nn
        elseif((1-b>-epsilon) && (1-b<epsilon))
            return nn+e2
        end       
    elseif(test==2)
        if((b>-epsilon && b<epsilon))
            return nn+e1
        elseif(1-b>-epsilon && 1-b<epsilon) 
            return nn+e1+e2
        end
    end
    return nn+xx1
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

function Lorentz(x,v,r)   # I know this is not beautifull :-P
    e1={1,0}
    e2={0,1}
    x1=dot(x,e1)
    x2=dot(x,e2)
    v1=dot(v,e1)
    v2=dot(v,e2)
    m=v2/v1
    if(dot(int(x)-x,int(x)-x)<r)
      return int(x)
    end
        if(m>=1 && v2>0)
            xx=x2 
            yy=x1
            vx=v2
            vy=v1
            x=[xx,yy]
            v=[vx,vy]
            x=nextt(x,v,r)
            if(typeof(x[2])==Int64)
                xx2=dot(x,e2)
                xx1=dot(x,e1)
                xx=xx2 
                yy=xx1
                x=[xx,yy]
                return x
            end
            b=dot(x,e2)
            b=b-int(b)
            de=[0, 0]
            if b<0 
                b=b+1
                de=[0, 1] 
            end
            m=vy/vx
            centro=eff(m,b,r)
            x=int(x)-de+centro
            xx2=dot(x,e2)
            xx1=dot(x,e1)
            xx=xx2
            yy=xx1
            x=[xx,yy]
            return x 
        elseif (m>=1 && v2<0)
            xx=-x2
            yy=-x1
            vx=-v2
            vy=-v1
            x=[xx,yy]
            v=[vx,vy]
            x=nextt(x,v,r)
            if(typeof(x[2])==Int64)
                xx2=dot(x,e2) 
                xx1=dot(x,e1)
                xx=xx2
                yy=xx1
                x=[xx,yy]
                return x
            end
            b=dot(x,e2)
            b=b-int(b)
            de=[0, 0]
            if b<0 
                b=b+1
                de=[0, 1]
            end
            m=vy/vx
            centro=eff(m,b,r)
            x=int(x)-de+centro  
            xx2=dot(x,e2)
            xx1=dot(x,e1)
            xx=-xx2
            yy=-xx1 
            x=[xx,yy]
            return x
        elseif (m>=0 && m<1 && v2>0)
        #ideal
            xx=x1
            yy=x2
            vx=v1
            vy=v2
            x=[xx,yy]
            v=[vx,vy] 
            x=nextt(x,v,r)

            if(typeof(x[2])==Int64)
                xx2=dot(x,e2)
                xx1=dot(x,e1)
                xx=xx1
                yy=xx2
                x=[xx,yy]
                return x
            end
            b=dot(x,e2) 
            b=b-int(b)
            de=[0, 0]
            if b<0 
                b=b+1
                de=[0, 1]
            end
            m=vy/vx
            centro=eff(m,b,r)
 
            x=int(x)-de+centro 
            xx2=dot(x,e2) 
            xx1=dot(x,e1)
            xx=xx1
            yy=xx2
            x=[xx,yy]
            return x
        elseif (m>=0 && m<1 && v2<0)
            xx=-x1
            yy=-x2
            vx=-v1
            vy=-v2
            x=[xx,yy]
            v=[vx,vy]
            x=nextt(x,v,r)
            if(typeof(x[2])==Int64)
                xx2=dot(x,e2)
                xx1=dot(x,e1)
                xx=xx2
                yy=xx1
                x=[xx,yy]
                return x 
            end
            b=dot(x,e2)
            b=b-int(b)
            de=[0, 0]
            if b<0 
                b=b+1
                de=[0, 1]
            end
            m=vy/vx
            centro=eff(m,b,r)
            x=int(x)-de+centro 
            xx2=dot(x,e2)
            xx1=dot(x,e1)
            xx=-xx1
            yy=-xx2
            x=[xx,yy]
            return x
        elseif (m<=-1 && v2>0)
            xx=x2
            yy=-x1
            vx=v2
            vy=-v1   
            x=[xx,yy]
            v=[vx,vy]
            x=nextt(x,v,r)
            if(typeof(x[2])==Int64)
                xx2=dot(x,e2)
                xx1=dot(x,e1)
                xx=xx2
                yy=xx1 
                x=[xx,yy]
                return x
            end
            b=dot(x,e2)
            b=b-int(b)
            de=[0, 0]
            if b<0 
                b=b+1
                de=[0, 1]
            end #160
            m=vy/vx
            centro=eff(m,b,r)
            x=int(x)-de+centro 
            xx2=dot(x,e2)
            xx1=dot(x,e1)
            xx=-xx2
            yy=xx1
            x=[xx,yy]
            return x
        elseif (m<=-1 && v2<0) 
            xx=-x2
            yy=x1
            vx=-v2
            vy=v1   
            x=[xx,yy]
            v=[vx,vy]
            x=nextt(x,v,r)
            if(typeof(x[2])==Int64)
                xx2=dot(x,e2)
                xx1=dot(x,e1)
                xx=xx2
                yy=xx1
                x=[xx,yy]
                return x
            end
            b=dot(x,e2)
            b=b-int(b)
            de=[0, 0]
            if b<0 
                b=b+1 
                de=[0, 1]
            end
            m=vy/vx
            centro=eff(m,b,r)
            x=int(x)-de+centro 
            xx2=dot(x,e2)
            xx1=dot(x,e1)
            xx=xx2
            yy=-xx1
            x=[xx,yy] 
            return x
        elseif (m<0 && m>-1 && v2>0)
            xx=-x1
            yy=x2
            vx=-v1
            vy=v2  
            x=[xx,yy]
            v=[vx,vy]
            x=nextt(x,v,r)
            if(typeof(x[2])==Int64) 
                xx2=dot(x,e2)
                xx1=dot(x,e1)
                xx=xx2
                yy=xx1
                x=[xx,yy]
                return x
            end
            b=dot(x,e2)
            b=b-int(b)
            de=[0, 0] 
            if b<0 
                b=b+1
                de=[0, 1]
            end
            m=vy/vx
            centro=eff(m,b,r)
            x=int(x)-de+centro 
            xx2=dot(x,e2)
            xx1=dot(x,e1)
            xx=-xx1 
            yy=xx2
            x=[xx,yy]
            return x
        elseif (m<0 && m>-1 && v2<0)
            xx=x1
            yy=-x2
            vx=v1
            vy=-v2
            x=[xx,yy]
            v=[vx,vy]
            x=nextt(x,v,r)
            if(typeof(x[2])==Int64)
                xx2=dot(x,e2)
                xx1=dot(x,e1)
                xx=xx2
                yy=xx1
                x=[xx,yy]
                return x
            end
            b=dot(x,e2) 
            b=b-int(b)
            de=[0, 0]
            if b<0 
                b=b+1
                de=[0, 1]
            end
            m=vy/vx
            centro=eff(m,b,r)
            x=int(x)-de+centro 
            xx2=dot(x,e2) 
            xx1=dot(x,e1)
            xx=xx1
            yy=-xx2
            x=[xx,yy]
            return x
        end
end

function Lorentz3D(x,v,r)
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
#        while (int(tmax)>int(tmin+BigFloat("0.5")) )
        while(x1[1]!=x2[1] || x1[2]!=x3[1] || x2[2]!=x3[2])
            x_xy=[x[1], x[2]]
            x_xz=[x[1], x[3]]
            x_yz=[x[2], x[3]]                #positions in the planes xy, xz, and yz

            x1=Lorentz(x_xy,vn_xy,r)#calculate the position of obstacle that collide with the particle 
            x2=Lorentz(x_xz,vn_xz,r)# in the 2D Lorentz gas for the 3 planes
            x3=Lorentz(x_yz,vn_yz,r) 

            t1=(sqrt(dot(x1-x_xy,x1-x_xy))-r)/sqrt(dot(v_xy,v_xy))
            t2=(sqrt(dot(x2-x_xz,x2-x_xz))-r)/sqrt(dot(v_xz,v_xz))
            t3=(sqrt(dot(x3-x_yz,x3-x_yz))-r)/sqrt(dot(v_yz,v_yz))  #calculate the time used to reach the obstacle. 

            tmax=maximum([t1 t2 t3])    # if the minumum and maximum valueds of time are almost the same, then 
            tmin=minimum([t1 t2 t3])    # there is very probable a collision  
            x=x+v*(tmax-0.1*r)   
        end
        xr=collision(x,int(x),r,v)
        x=x+v*(1-2*r)
    end
    return xr
end

function ClassicLorentz(x,v,r)
    xr=false
    i=0.0
    while (xr==false)
        while (xr==false)
            xr=collision(x,int(x+v*(i)),r,v)
            i=i+1-2*r
        end
        t=sqrt(dot(xr-x,xr-x)/dot(v,v))
        if(dot(x+t*v-xr,x+t*v-xr)>dot(xr-x,xr-x))
            xr=false
        end
    end
    return xr
end