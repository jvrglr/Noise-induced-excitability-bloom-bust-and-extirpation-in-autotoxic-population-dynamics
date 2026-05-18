module functions
  !Module with public functions.

  implicit none

  double precision, public, parameter :: pi = &
   3.1415926535897932384626433832795028841972

  double precision, private :: f
contains

  pure function rev(string) result(reverse_string)
      character(len=*), intent(in) :: string
      character(len=len(string)) :: reverse_string
      integer :: i, n

      n = len(string)
      do i = 1, n
          reverse_string(n-i+1:n-i+1) = string(i:i)
      end do

  end function rev

  character(len=30) function str(k) !fUNCTION
    implicit none
!   "Convert an integer*4 to string."
    integer, intent(in) :: k
    write (str, *) k !write to a string
    str = adjustl(str)
  end function str

  character(len=30) function str8(k) !fUNCTION
    implicit none
!   "Convert an integer*4 to string."
    integer*8, intent(in) :: k
    write (str8, *) k !write to a string
    str8 = adjustl(str8) !Adding trim() here doesn't avoid //trim(str8(...))//
  end function str8

  function round(val, n) result(a)
    implicit none
    double precision :: val, a
    integer :: n
    a = anint(val*10.0**n)/10.0**n
  end function round

  function dble_mean(x) result(a)
    !compute average of double precision vector
    implicit none
    double precision :: a
    double precision, dimension (:), intent(in) :: x
    a = sum(x)/size(x)
  end function dble_mean

  function dble_var(x) result(a)
    !compute variance of double precision vector
    implicit none
    double precision :: a
    double precision, dimension (:), intent(in) :: x
    double precision, dimension (:), allocatable :: dummy
    integer*4 :: i, len

    len=size(x)
    allocate(dummy(len))

    dummy=x-dble_mean(x)
    a=sum(dummy*dummy)/len
    deallocate(dummy)
  end function dble_var

  function dble_err(x) result(a)
    !compute standard error of double precision vector
    implicit none
    double precision :: a
    double precision, dimension (:), intent(in) :: x
    a=sqrt(dble_var(x)/size(x))
  end function dble_err

  function G(x,mu,sig) result(a)
    !Gaussian function
    implicit none
    double precision :: a
    double precision, intent (in):: x,mu,sig
    a=(1.0/(sqrt(2.0*pi)*sig))*exp(-(x - mu)*(x - mu) / (2.0d0 *sig*sig))
  end function G

end module functions