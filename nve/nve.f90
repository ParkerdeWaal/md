program nve
implicit none
real*8 :: dt, H, px, py, vx, vy, fx, fy, potential
integer :: i,n,nstep,natoms

!user settings
nstep = 5000000
dt = 0.02d0

!initial coordinates, velocity and forces
px = 0.5d0
py = 0.5d0
vx = 0.2d0
vy = 0.2d0
fx = 0d0
fy = 0d0

!calculate initial forces at t0
call calculate_Forces(px,py,fx,fy)

do n=1,nstep
	!velocity verlet
	vx = vx - (dt/2d0) * fx
	vy = vy - (dt/2d0) * fy

	px = px + dt * vx
	py = py + dt * vy
	
	call calculate_Forces(px,py,fx,fy)

	vx = vx - (dt/2d0) * fx
	vy = vy - (dt/2d0) * fy


	if(MOD(n, 500).eq.0)then
		potential = 0d0

		potential = 4d0*(px**2d0+py**2d0-1d0)**2d0*py**2d0 - &
				exp(-4d0*((px-1d0)**2d0+py**2d0)) - &
				exp(-4d0*((px+1d0)**2d0+py**2d0)) + &
				exp(8d0*(px-1.5d0)) + &
				exp(-8d0*(px+1.5d0)) + &
				exp(-4d0*(py+0.25d0)) + &
				0.2d0 * exp(-8d0*px**2d0)

		H = 0.5d0 * (vx**2+vy**2) + potential

		write(*,*) n, H, px, py
	endif
enddo

end program nve

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
fy = ff1y - ff2y - ff3y + ff4y + ff5y + ff6y + ff7y

return
end subroutine calculate_Forces