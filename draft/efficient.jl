normsq(v) = dot(v,v)

function collision(x1, x2, r, v)
    v2 = normsq(v)
    b = dot((x1 - x2),v) / v2
    c = normsq(x1 - x2) - r^2
    c /= v2

    if b^2 - c < 0
        return false
    end

    t = -b - sqrt(b^2-c)
    x = v*t + x1

    x
end

function velo_col(x1, x2, v)
    n = x1 - x2
    n /= norm(n)
    vn = dot(n,v)*n
    v = v - 2*vn
    v /= norm(v)

    v
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

# function frac(alpha, epsilon)
#     p = 0; q = 0
# 	h1 = 0
# 	h2 = 1
# 	k1 = 1
# 	k2 = 0
# 	mm1 = alpha
# 	a1 = ifloor(mm1)
#     while abs(k2*alpha-h2) > epsilon
#         mm2 = 1/(mm1 - a1)
#         @show mm2
#         a2 = ifloor(mm2)
# 		h3 = a1*h2 + h1
# 		k3 = a1*k2 + k1
# 		p = h3
# 		q = k3
#         h1=h2
#         h2=h3
#         a1=a2
#         k1=k2
#         k2=k3
#         mm1=mm2
#     end
#     return q,p
# end

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
        if b < 0.5
			(q, p) = frac(m, 2.*b)
		else
			(q, p) = frac(m, 2.*(1. - b))
		end
		b = mod(m*q + b, 1)
		kn += q
    end

	while b > epsilon && 1 - b > epsilon
		if b < 0.5
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

normalise(x) = x / norm(x)

function Lorentz3D(x,v,r)
    xr = false  # whether have found a collision

    while xr == false

        # velocities in the planes xy, xz, and yz:
        v_xy = [v[1], v[2]]
        v_xz = [v[1], v[3]]
        v_yz = [v[2], v[3]]

        # normalise:
        vn_xy = normalise(v_xy)
        vn_xz = normalise(v_xz)
        vn_yz = normalise(v_yz)

        tmax=0
        tmin=-2

        while int(tmax) > int(tmin+0.5)

            x_xy = [x[1], x[2]]
            x_xz = [x[1], x[3]]
            x_yz = [x[2], x[3]]                #positions in the planes xy, xz, and yz

            x1 = Lorentz(x_xy, vn_xy, r) # calculate the position of obstacle that collide with the particle
            x2 = Lorentz(x_xz, vn_xz, r) # in the 2D Lorentz gas for the 3 planes
            x3 = Lorentz(x_yz, vn_yz, r)

            t1_squared = normsq(x1 - x_xy) / normsq(v_xy)
            t2_squared = normsq(x2 - x_xz) / normsq(v_xz)
            t3_squared = normsq(x3 - x_yz) / normsq(v_yz)  #calculate the time used to reach the obstacle.

            tmax = sqrt(max(t1_squared, t2_squared, t3_squared))    # if the minumum and maximum values of time are almost the same, then
            tmin = sqrt(min(t1_squared, t2_squared, t3_squared))    # there is very probably a collision

            x += v*(tmax-0.5)
   #         println(x1,x2,x3," ", int(x), " ", (tmax), " ", tmin)

        end

        xr = collision(x, int(x), r, v)
#        println(xr)
        if(sqrt(dot(x,x))>10^20)
#            xr=[0,0,0]
            println(sqrt(dot(x,x)))
        end
        x += v
    end
    return xr
end
