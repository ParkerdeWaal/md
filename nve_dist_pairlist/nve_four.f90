program nve
implicit none
real*8 :: xx(4*2), ff(4*2), vv(4*2)
real*8 :: xij,yij,rij,dif,potential,H
real*8 :: dt,d0
integer :: i,j,k,n,nstep,natoms,pairlist(8)

!settings
nstep = 500000000
dt = 0.0025d0
d0 = 1d0

!atom count and pairlists
natoms = 4
pairlist = (/1,2,2,3,3,4,4,1/)

!initial coordinates and velocities
xx(1) = 0.0d0
xx(2) = 0.0d0
xx(3) = 1.0d0
xx(4) = 0.0d0
xx(5) = 1.0d0
xx(6) = 1.0d0
xx(7) = 0.0d0
xx(8) = 1.0d0

vv = 0d0
ff= 0d0

vv(1) = 0.0050d0
vv(2) = -0.0050d0
vv(3) = -0.0050d0
vv(4) = 0.0050d0
vv(5) = 0.0050d0
vv(6) = 0.0050d0
vv(7) = -0.0050d0
vv(8) = 0.0050d0

! initial forces
call calculate_Forces(xx,ff,pairlist,d0)

do n=1,nstep
	! velocity verlet
	do i=1,natoms
		do k=1,2
			vv(i*2-(2-k)) = vv(i*2-(2-k)) + (dt/2d0) * ff(i*2-(2-k))
			xx(i*2-(2-k)) = xx(i*2-(2-k)) + dt * vv(i*2-(2-k))
		enddo
	enddo
	call calculate_Forces(xx,ff,pairlist,d0)
	do i=1,natoms
		do k=1,2
			vv(i*2-(2-k)) = vv(i*2-(2-k)) + (dt/2d0) * ff(i*2-(2-k))
		enddo
	enddo
	!calculate energy
	if(MOD(n, 500).eq.0)then
		call calculate_Energy(n,xx,vv,pairlist,d0)
		write(81,*) natoms
		write(81,*)
	 	do i=1,natoms
			write(81,*) 'C ', xx(i*2-(2-1))*10d0, xx(i*2-(2-2))*10d0, 1d0
		enddo
		call flush(81)
	endif
enddo
end program nve

subroutine calculate_Energy(n,xx,vv,pairlist,d0)
implicit none
real*8 :: xx(4*2),vv(4*2),potential,v,H,d0,dif,xij,yij,rij
integer :: i,j,k,n,pairlist(8)

potential = 0d0
v = 0d0 
H = 0d0

do i=1,size(pairlist),2
	xij = xx(pairlist(i)*2-(2-1)) - xx(pairlist(i+1)*2-(2-1))
	yij = xx(pairlist(i)*2-(2-2)) - xx(pairlist(i+1)*2-(2-2))

	rij = sqrt(xij*xij+yij*yij)
	dif = rij - d0

	potential = potential + 0.5d0 * 0.62d0 * dif ** 2d0
	do k=1,2
		v = v + vv(pairlist(i)*2-(2-k))**2d0
	enddo
enddo

H = 0.5d0 * v + potential
write(*,*) n,H    

end subroutine calculate_Energy

subroutine calculate_Forces(xx,ff,pairlist,d0)
implicit none
real*8 :: xx(4*2), ff(4*2)
real*8 :: xij, yij, rij, d0, dif, fx, fy
integer :: i,j,k,n,pairlist(8)

ff = 0 

do i=1,size(pairlist),2
	xij = xx(pairlist(i)*2-(2-1)) - xx(pairlist(i+1)*2-(2-1))
	yij = xx(pairlist(i)*2-(2-2)) - xx(pairlist(i+1)*2-(2-2))

	rij = sqrt(xij*xij+yij*yij)
	dif = rij - d0

	fx = 0.62d0*(dif)*xij/rij
	fy = 0.62d0*(dif)*yij/rij
	ff(pairlist(i)*2-(2-1)) = ff(pairlist(i)*2-(2-1)) - fx
	ff(pairlist(i)*2-(2-2)) = ff(pairlist(i)*2-(2-2)) - fy

	fx = -0.62d0*(dif)*xij/rij
	fy = -0.62d0*(dif)*yij/rij
	ff(pairlist(i+1)*2-(2-1)) = ff(pairlist(i+1)*2-(2-1)) - fx
	ff(pairlist(i+1)*2-(2-2)) = ff(pairlist(i+1)*2-(2-2)) - fy
enddo
return
end subroutine calculate_Forces