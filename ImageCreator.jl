const c=299792458.

function IC(l,p,h,X,Y,Z,tilt,azimut,phase,amplitude,order)
  const Rx=vcat([-1 0 0],[0 1 0],[0 0 1])
  const Ry=vcat([1 0 0],[0 -1 0],[0 0 1])
  const Rz=vcat([1 0 0],[0 1 0],[0 0 -1])
  const i0=amplitude*[sin(tilt).*cos(azimut), sin(tilt).*sin(azimut), cos(tilt)]
  POS=[0 0 0 0 0 0 0]
  for i=0:order
    perc=round(i/order*1000)/10
    println("$perc %")
    for j=0:order
      for k=0:order
        n = int(i+j+k);
        ip = (-1)^n*Rx^abs(i)*Ry^abs(j)*Rz^abs(k)*i0;
        x = 2*i*l+(-1)^abs(i)*X;
        y =	2*j*p+(-1)^abs(j)*Y;
        z = 2*k*h+(-1)^abs(k)*Z;
        POS=[POS; ip' x  y  z n ]
        POS=[POS; (-Rz*ip)' x  y  -z n+1 ]
        POS=[POS; (-Ry*ip)' x  -y  z n+1 ]
        POS=[POS; (-Ry*-Rz*ip)' x  -y  -z n+2 ]
        POS=[POS; (-Rx*ip)' -x  y  z n+1 ]
        POS=[POS; (-Rx*-Rz*ip)' -x  y  -z n+2 ]
        POS=[POS; (-Rx*-Ry*ip)' -x  -y  z n+2 ]
        POS=[POS; (-Rx*-Ry*-Rz*ip)' -x  -y  -z n+3 ]
      end
    end
  end
  POS = POS[[2:end], :]
  perm=sortperm(POS[:,7])
  return POS[perm,:]
end







