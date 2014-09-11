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
dx=0.01
dy=0.01
dz=0.01
nx=int(2/dx)
ny=int(2/dy)
const  x=linspace(0.,2.,nx)

const  y=linspace(0.,2.,ny)

z=0
println("Computing the radiation...")
Z=zeros(nx,ny)
Ez=zeros(nx,ny)
Bz=zeros(nx,ny)
for i=1:nx
  perc=round(i/nx*1000)/10
  println("$perc %")
  for j=1:ny
    r=[x[i],y[j],z]
    E=zeros(Complex128,3)
    B=zeros(Complex128,3)
    for m=1:N
      p=vec(I[m,:])
      R=vec(xyz[m,:])
      Ed,Bd=Hertz_dipole (r, p, R, phase[m], freq)
      E+=Ed
      B+=Bd
    end
    Z[i,j]=sqrt(sum(abs(E).^2))/sqrt(sum(abs(B).^2))*mu0 #real(E).^2#0.5*numpy.cross(E.T,conjugate(B.T))
    Ez[i,j]=abs(E[3])
    Bz[i,j]=abs(B[3])
  end
end

using PyPlot
figure(1)
title("\$E_z\$ (V/m)")
pcolor(x,y,Ez,cmap="jet")
clim(0,20)
colorbar()
axis("scaled")
xlim(0,L)
ylim(0,l)
xlabel("\$x\$/m")
ylabel("\$y\$/m")
savefig("Ez.png",bbox="tight")
clf()


figure(2)
title("\$H_z\$ (A/m)")
pcolor(x,y,Bz/mu0,cmap="jet")
#clim(0,.2)
colorbar()
axis("scaled")
xlim(0,L)
ylim(0,l)
xlabel("x/m")
ylabel("\$y\$/m")
savefig("Hz.png",bbox="tight")
clf()


figure(3)
title("\$Z\$ (\$ \\Omega\$)")
pcolor(x,y,Z,cmap="jet")
clim(0,1000)
colorbar(ticks=[0,100,200,300,120*pi,500,60*pi^2,700,800,900,1000])

axis("scaled")
xlim(0,L)
ylim(0,l)
#xlabel("x/m")
ylabel("\$y\$/m")
savefig("Z.png",bbox="tight")
clf()
