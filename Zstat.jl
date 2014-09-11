const c = 299792458.
const mu0 = 4*pi*1e-7
const eps0 = 1/(mu0*c^2)

include("ImageCreator.jl")
include("dipole.jl")

#cavity
L=5
l=4
h=3
order=100

#dipole
X=1.2
Y=1.3
Z=1.0
tilt=0.3
azimut=2*pi*rand()
phase=0

const freq=1e9
#dipole moment
#total time averaged radiated power P= 1 W dipole moment => |p|=sqrt(12πcP/µOω⁴)
Pow=1
amplitude=sqrt(12*pi*c*Pow/(mu0*(2*pi*freq)^4))

#POS=IC(L,l,h,X,Y,Z,tilt,azimut,phase,amplitude,order)
using HDF5,JLD

@load "POS100.jld" POS

numberofimages=1+2*order+(2*order*(order+1)*(2*order+1))/3
POS=POS[1:numberofimages,:]

M=500
#observation points
x=rand(M)*(L-0.5)+0.5
y=rand(M)*(l-0.5)+0.5
z=rand(M)*(h-0.5)+0.5


println("Computing the radiation...")
Z=zeros(M)
for i=1:M
  perc=round(i/M*1000)/10
  println("$perc %")
  r=[x[i],y[i],z[i]]
  E=zeros(Complex128,3)
  B=zeros(Complex128,3)
  for m=1:numberofimages
    p=vec(POS[m,1:3])
    R=vec(POS[m,4:6])
    Ed,Bd=Hertz_dipole (r, p, R, phase, freq)
    E+=Ed
    B+=Bd
  end
  Z[i]=sqrt(sum(abs(E).^2))/sqrt(sum(abs(B).^2))*mu0 #real(E).^2#0.5*numpy.cross(E.T,conjugate(B.T))
end

@save "Z.jld" Z=Z

using PyPlot
plt.hist(Z,30)
xlabel("\$ Z\$ in \$\\Omega\$")
ylabel("Probability")
