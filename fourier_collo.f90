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

program fourier_colloc
    use m
    implicit none
    
    integer, parameter :: n=16
    real, dimension(n,n) :: B, D, D2, R, beta
    real, dimension(n) :: U, f
    integer :: i, j
    real :: be=1.5, gamma=1., pi=acos(-1.), x, hij=0.5, alpha=2

    print*, "start"

	x=2*pi/n

    do i=1,n
            f(i)=cos(alpha*real(i)*x)
    end do

    do i=1,n
        do j=1,n
            if (i.eq.j) then
                D(i,j)=0
            else
            	if (modulo(n,2).ne.0) then
            		D(i,j)=(-1)**(i+j)/(2*sin(hij))
            	else 
                	D(i,j)=(-1)**(i+j)/(2*tan(hij))
                end if
            end if
        end do
    end do
    
    D2=matmul(D,D)

    do i=1,n
    	do j=1,n
    		if (i==j) then
    			beta(i,j) = be
    		else 
    			beta(i,j)=0
    		end if
    	end do
    end do

    do i=1,n
        do j=1,n
            R(i,j)=D2(i,j)+beta(i,j)
    !        R(i,j)=(D(i,j)**2)-beta
        end do
    end do

    U=matmul(inv(R), f)

    open (unit=50, file="result.txt")
    
    do i=1,n
        write(50,*) U(i)
    end do

    close (50)
    print*, "end"

end program


