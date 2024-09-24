module m
	implicit none
	
	contains
	!https://stackoverflow.com/questions/43728087/computing-inverse-of-a-matrix-using-lapack-blas-in-fortran90
function inv(A) result(Ainv)
    implicit none
    real,intent(in) :: A(:,:)
    real            :: Ainv(size(A,1),size(A,2))
    real            :: work(size(A,1))            ! work array for LAPACK
    integer         :: n,info,ipiv(size(A,1))     ! pivot indices

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call SGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call SGETRI(n,Ainv,n,ipiv,work,n,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
end function inv

end module

program matt
    use m
    implicit none
    integer, parameter :: n=20
    real, dimension(n,n) :: I, A, Ainv, B, C, R, Rinv
    real, dimension(1,n) :: f
    real, dimension(1,n) :: C1
    real :: alpha=1, h=0.1, gamma=1, pi=acos(-1.), x
    integer :: k, j

    print*, "start"

    !la matrice identité
    do k=1,n
        do j=1,n
            if (k.eq.j) then
            	I(k,j)=1
            else
            	I(k,j)=0
            end if
        end do
    end do

    ! la matrice A et son inverse
    do k=1,n
        do j=1,n
            if (k.eq.j) then
                A(k,j)=-2/(h**2)
            !else if ((k.eq.j-1) .or. (k.eq.j+1)) then !ajoute les CL !!!!
            else if ((k.eq.j-1) .or. (k.eq.j+1) .or. ((k.eq.n) .and. (j.eq.1)) .or. ((j.eq.n) .and. (k.eq.1))) then
                A(k,j)=1/(h**2)
            else
            A(k,j)=0
            end if
        end do
    end do

    Ainv=inv(A)

    !la matrice B
    do k=1,n
        do j=1,n
            if (k.eq.j) then
                B(k,j)=5./6.
            !else if ((k.eq.j-1) .or. (k.eq.j+1)) then !ajoute les CL !!!!
            else if ((k.eq.j-1) .or. (k.eq.j+1) .or. ((k.eq.n) .and. (j.eq.1)) .or. ((j.eq.n) .and. (k.eq.1))) then
                B(k,j)=1./12.
            else
            B(k,j)=0
            end if
        end do
    end do
    
    C = matmul(Ainv, B)
    
    R=alpha*I-C
    
    !matrice globale resultante
    !do k=1,n
    !	print*, R(k,:)
    !end do 
    
    x=2./n
    
    do j=1,n
    	f(1,j)=(alpha-(gamma*pi)**2)*cos(gamma*pi*real(j)*x)
    !	print*, f(1,j)
    end do
    
    Rinv=inv(R)
    
    open (unit=20, file="sol.txt")
    open(unit=30, file="sol_exacte.txt")
    
    C1=matmul(f, Rinv)
    do j=1,n
    	write(20,*) C1(:,j)
    end do
    
    do j=1,n
    	write(30,*) cos(gamma*pi*real(j)*x)
    end do

	close(20)
	close(30)
    print*, "end"

end program
