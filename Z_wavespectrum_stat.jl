include("dipole.jl")

const c = 299792458.
const mu0 = 4*pi*1e-7
const eps0 = 1/(mu0*c^2)

#dipoles:
const N=1000;
const distance=1000+1000*rand(N)
const theta=2*acos(rand(N)-1)
const phi=2*pi*rand(N)
xyz = hcat(distance.*sin(theta).*cos(phi), distance.*sin(theta).*sin(phi), distance.*cos(theta))
#moments:
const th_dip=2*acos(rand(N)-1)
const phi_dip=2*pi*rand(N)

#dipole moment
#total time averaged radiated power P= 1 W dipole moment => |p|=sqrt(12πcP/µOω⁴)
const Pow=1
const amp_max=sqrt(12*pi*c*Pow/(mu0*(2*pi*freq)^4))/N
const amplitude=amp_max*rand(N)
const phase=2*pi*rand(N)
I=hcat(amplitude.*sin(th_dip).*cos(phi_dip), amplitude.*sin(th_dip).*sin(phi_dip), amplitude.*cos(th_dip))

const freq=1e9



#observation points
M=500
#observation points
x=rand(M)*(10)-5
y=rand(M)*(10)-5
z=rand(M)*(10)-5

println("Computing the radiation...")
Z=zeros(M)
for i=1:M
  perc=round(i/M*100)
  println("$perc %")
  r=[x[i],y[i],z[i]]
  E=zeros(Complex128,3)
  B=zeros(Complex128,3)
  for m=1:N
    p=vec(I[m,:])
    R=vec(xyz[m,:])
    Ed,Bd=Hertz_dipole (r, p, R, phase[m], freq)
    E+=Ed
    B+=Bd
  end
  Z[i]=sqrt(sum(abs(E).^2))/sqrt(sum(abs(B).^2))*mu0 #real(E).^2#0.5*numpy.cross(E.T,conjugate(B.T))
end

using PyPlot
plt.hist(Z)
