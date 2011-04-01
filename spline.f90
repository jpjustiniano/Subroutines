Program testSpline
	integer, parameter:: n=15,no=30
	integer:: i, j
	real,dimension(1:n)::x,y,p2
	real,dimension(1:no)::xo,yo
	integer :: horaf, horaini, horasdia, diasmes
	
	Print *, 'Datos de entrada:'

	horaini=10
	horasdia= no
	 
	 
	x= (/10.1,10.8,11.0,11.6,12.0,12.6,13.5,14.0,15.5,15.9,16.4,17.1,18.0,19.0,19.2/)

	y= (/120.0, , , , , , , /)
	 
	xo(1)= 0.5+int(x(1))
	Do i=2,horasdia
		xo(i)= xo(i-1)+0.2
	end do

	call spline(n,x,y,1.,1. ,no,xo,yo)
	Do i = 1,n
		Print *, 'x',i,':',x(i), '  f(x)',i,':',y(i)
	end do
	Print *, 'Datos de salida:'
	Do i = 1,no
		Print *, 'x',i,':',xo(i), '  f(x)',i,':',yo(i)
	end do
	
	open (12,file='spline.txt') 
	write (12, *) 'Datos de entrada:'
	write (12, 100) (i,x(i),y(i), i=1,n)
	write (12, *) 'Datos de salida:'
	write (12, 100) (i,xo(i),yo(i), i=1,no)
	close (12)
100 Format ('X(',I2')',F6.1,' ;',F7.1)

	call CUBIC_SPLINE(n, x, y, P2)





contains

subroutine spline(n,x,y,yp1,ypn,no,xo,yo)
!***************************************************************
! This subroutine will calculate the cubic-spline intepolated value of
! given any function for single variable of any dimension 
!
! inputs/output: 
! n : (input) [scalar] dimension of input array
! x : (input) [vector of length 'n'], independent variable of dimension 'n'
! y : (input)[vector of length 'n']' the function of x (an array of dimension 'n'
! yp1 : (input) [scalar] =1.e+30 :: the routine is signaled to set
! the corresponding boundary condition for a natural spline, with zero
! second derivative on that boundary
! ypn : (input) same to yp1 = 1.e+30
! no : (input) [scalar] dimension output array
! xo : (input) [vector of length 'no'], 'x' values at which you want to spline
! yo : (output) [vector of length 'no'], splined value
!
! Note: This works properly in any inter fortran compiler and f90-compiler
!***************************************************************
! Tanmoy Das.
! Feb 13, 2008.
!***************************************************************
	integer::n,no,io
	real,dimension(1:n)::x,y,y2,u
	real,dimension(1:no)::xo,yo

	if (yp1.gt.0.99e30) then
		y2(1)=0.
		u(1)=0.
	else
		y2(1)=-0.5
		u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	endif
	do  i=2,n-1
		sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
		p=sig*y2(i-1)+2.
		y2(i)=(sig-1.)/p
		u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))& 
		/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
	end do
	if (ypn.gt.0.99e30) then
		qn=0.
		un=0.
	else
		qn=0.5
		un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
	endif
	y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
	do k=n-1,1,-1
		y2(k)=y2(k)*y2(k+1)+u(k)
	end do
!***************************************************************
! now for any arbitrary x
!***************************************************************
	do io = 1,no
		klo=1
		khi=n
		1 if (khi-klo.gt.1) then
			k=(khi+klo)/2.
				if (x(k).gt.xo(io)) then
					khi=k
				else
					klo=k
				endif
			goto 1
		endif
		h=x(khi)-x(klo)
		if (h.eq.0.) then 
			write (*,*) 'bad XA input'
			read (*,*)
		end if
		a=(x(khi)-xo(io))/h
		b=(xo(io)-x(klo))/h
		yo(io)=a*y(klo)+b*y(khi)+(a*(a*a-1.)*y2(klo)+b*(b*b-1.)*y2(khi))*h*h/6.
	enddo

	return
end subroutine spline

SUBROUTINE CUBIC_SPLINE (N, XI, FI, P2)
!
! Function to carry out the cubic-spline approximation
! with the second-order derivatives returned.
!
  INTEGER :: I
  INTEGER, INTENT (IN) :: N
  REAL, INTENT (IN), DIMENSION (N):: XI, FI
  REAL, INTENT (OUT), DIMENSION (N):: P2
  REAL, DIMENSION (N):: G, H
  REAL, DIMENSION (N-1):: D, B, C
!
! Assign the intervals and function differences
!
  DO I = 1, N
    H(I) = XI(I+1) - XI(I)
    G(I) = FI(I+1) - FI(I)
  END DO
!
! Evaluate the coefficient matrix elements
  DO I = 1, N-1
    D(I) = 2*(H(I+1)+H(I))
    B(I) = 6*(G(I+1)/H(I+1)-G(I)/H(I))
    C(I) = H(I+1)
  END DO
!
! Obtain the second-order derivatives
!
  CALL TRIDIAGONAL_LINEAR_EQ (N-1, D, C, C, B, G)
  P2(1) = 0
  P2(N+1) = 0
  DO I = 2, N 
    P2(I) = G(I-1)
  END DO
END SUBROUTINE CUBIC_SPLINE
!
SUBROUTINE TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)
!
! Functione to solve the tridiagonal linear equation set.
!
  INTEGER, INTENT (IN) :: L
  INTEGER :: I
  REAL, INTENT (IN), DIMENSION (L):: D, E, C, B
  REAL, INTENT (OUT), DIMENSION (L):: Z
  REAL, DIMENSION (L):: Y, W
  REAL, DIMENSION (L-1):: V, T
!
! Evaluate the elements in the LU decomposition
!
  W(1) = D(1)
  V(1)  = C(1)
  T(1)  = E(1)/W(1)
  DO I = 2, L - 1
    W(I) = D(I)-V(I-1)*T(I-1)
    V(I) = C(I)
    T(I) = E(I)/W(I)
  END DO
  W(L) = D(L)-V(L-1)*T(L-1)
!
! Forward substitution to obtain y
!
  Y(1) = B(1)/W(1)
  DO I = 2, L
    Y(I) = (B(I)-V(I-1)*Y(I-1))/W(I)
  END DO
!
! Backward substitution to obtain z
  Z(L) = Y(L)
  DO I = L-1, 1, -1
    Z(I) = Y(I) - T(I)*Z(I+1)
  END DO
END SUBROUTINE TRIDIAGONAL_LINEAR_EQ

End program
