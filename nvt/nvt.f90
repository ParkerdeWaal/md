program nvt
implicit none
real*8 :: dt, H, Hs, px, py, vx, vy, fx, fy, fxold, fyold, potential
real*8 :: gmma, kbt, sigma1, sigma2, crv, s0, s1, s2, dx1, dx2, dy1, dy2
integer :: i,n,nstep,natoms

!user settings
nstep = 1000000000
dt = 0.04d0
gmma = 1.5d0
kbt = 1d0
H = 0d0
Hs = 0d0
call R250Init(2018) !2018 is the thermostat seed

!integration scheme constants
sigma1 = (dt*kbt)/gmma * (2d0-(3d0-4d0*exp(-gmma*dt)+exp(-2d0*gmma*dt))/(gmma*dt))
sigma2 = kbt*(1d0-exp(-2d0*gmma*dt))
crv = (kbt*(1d0-exp(-gmma*dt))**2d0)/(gmma*sqrt(sigma1*sigma2))
s0 = exp(-gmma*dt)
s1 = (1d0-s0)/(gmma*dt)
s2 = (1d0-s1)/(gmma*dt)

!initial coordinates, velocity and forces
px = 0.5d0
py = 0.5d0
vx = 0.2d0
vy = 0.2d0
fx = 0d0
fy = 0d0

dx1 = 0d0
dx2 = 0d0
dy1 = 0d0
dy2 = 0d0

!calculate initial forces at t0
call calculate_Forces(px,py,fx,fy)

do n=1,nstep
	call boxmul(dx1,dy1)
	call boxmul(dx2,dy2)
	px = px + s1*vx*dt + s2*dt**2d0*fx + sigma1*dx1
	py = py + s1*vy*dt + s2*dt**2d0*fy + sigma1*dy1
	fxold = fx
	fyold = fy
	call calculate_Forces(px,py,fx,fy)
	vx = s0*vx + (s1-s2)*dt*fxold + s2*dt*fx + sigma2*(crv*dx1+sqrt(1d0-crv)*dx2)
	vy = s0*vy + (s1-s2)*dt*fyold + s2*dt*fy + sigma2*(crv*dy1+sqrt(1d0-crv)*dy2)
	
	!call calculate_Forces(px,py,fx,fy)
    if(MOD(n, 500).eq.0)then
    	H = H + (vx**2d0+vy**2d0)/(2*nstep)
    	Hs = sqrt(H)
    	write(*,*) n, Hs, px, py
	endif
enddo

end program nvt

subroutine calculate_Forces(px,py,fx,fy)
implicit none
real*8 :: px, py, fx, fy
real*8 :: ff1x, ff2x, ff3x, ff4x, ff5x, ff6x, ff7x
real*8 :: ff1y, ff2y, ff3y, ff4y, ff5y, ff6y, ff7y

fx = 0d0
fy = 0d0

ff1x = 16d0*px*py**2d0*(px**2d0+py**2d0-1d0)
ff1y = 8d0*py*(px**2d0+py**2d0-1d0)**2d0+16d0*py**3d0*(px**2d0+py**2d0-1d0)

ff2x = -8d0*(px-1d0)*exp(-4d0*((px-1d0)**2d0+py**2d0))
ff2y = -8d0*py*exp(-4d0*((px-1d0)**2d0+py**2d0))

ff3x = -8d0*(px+1d0)*exp(-4d0*((px+1d0)**2d0+py**2d0))
ff3y = -8d0*py*exp(-4d0*((px+1d0)**2d0+py**2d0))

ff4x = 8d0*exp(8d0*(px-1.5d0))
ff4y = 0d0

ff5x = -8d0*exp(-8d0*(px+1.5d0))
ff5y = 0d0

ff6x = 0d0
ff6y = -4d0*exp(-4d0*(py+0.25d0))

ff7x = (-16d0*px*exp(-8d0*px**2d0))/5d0
ff7y = 0d0

fx = ff1x - ff2x - ff3x + ff4x + ff5x + ff6x + ff7x
fx = -fx
fy = ff1y - ff2y - ff3y + ff4y + ff5y + ff6y + ff7y
fy = -fy
return
end subroutine calculate_Forces