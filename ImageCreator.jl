const c=299792458.

function IC(l,p,h,X,Y,Z,tilt,azimut,phase,amplitude,order)
  const Rx=vcat([-1 0 0],[0 1 0],[0 0 1])
  const Ry=vcat([1 0 0],[0 -1 0],[0 0 1])
  const Rz=vcat([1 0 0],[0 1 0],[0 0 -1])
  const i0=amplitude*[sin(tilt).*cos(azimut), sin(tilt).*sin(azimut), cos(tilt)]
  POS=zeros((order+1)^3*8,7)#[0 0 0 0 0 0 0]
  count=1
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
        POS[count,:]=[ip' x  y  z n ]
        POS[count+1,:]=[(-Rz*ip)' x  y  -z n+1 ]
        POS[count+2,:]=[(-Ry*ip)' x  -y  z n+1 ]
        POS[count+3,:]=[(-Ry*-Rz*ip)' x  -y  -z n+2 ]
        POS[count+4,:]=[(-Rx*ip)' -x  y  z n+1 ]
        POS[count+5,:]=[(-Rx*-Rz*ip)' -x  y  -z n+2 ]
        POS[count+6,:]=[(-Rx*-Ry*ip)' -x  -y  z n+2 ]
        POS[count+7,:]=[(-Rx*-Ry*-Rz*ip)' -x  -y  -z n+3 ]
	count+=8
      end
    end
  end
  perm=sortperm(POS[:,7])
  return POS[perm,:]
end
