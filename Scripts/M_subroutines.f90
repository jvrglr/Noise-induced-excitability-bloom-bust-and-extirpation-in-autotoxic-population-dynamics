module Subroutines

    !Module with public subroutines.
    !Depends on declarations_module and functions.
  
    use declarations_module
    use functions
  
  contains

    subroutine Draw_trajectory(until_t,init_t,dt,init_x,init_y,points,name_traj,until_extinction,print_theta_w)
      implicit none
      double precision, intent(in) :: until_t,init_t,dt,init_x,init_y
      double precision :: x0,y0,theta,theta_old,dtheta,theta_w
      integer*8, intent(in) :: points
      character (len=*),intent(in) :: name_traj
      logical,intent(in) :: until_extinction
      double precision :: t0,tf
      double precision,dimension(:),allocatable :: xs,ts,ys,thetas_w
      integer*8 :: ii,points_before_abs
      logical, optional :: print_theta_w
      logical :: print_theta_w_

      x0 = init_x
      y0 = init_y
      t0 = init_t
      theta_w = 0.0d0
      print_theta_w_ = .false.
      if ( present(print_theta_w) ) then
        print_theta_w_ = print_theta_w
      end if
      allocate(xs(0:points),ys(0:points),ts(0:points),thetas_w(0:points))
      thetas_w(0)=0.0d0
      xs(0) = x0
      ys(0) = y0
      ts(0) = t0
      theta_old = ATAN2(y0-1.0,x0-1.0)
      print*, "theta_old=",theta_old/pi
      points_before_abs = points
      do ii = 1, points, 1
        tf = init_t+(until_t-init_t)*ii/points
        ! print *, "tf =",tf,"t0=",t0
        call trajectory_milstein(t0,x0,y0,tf,dt)
        
        if ( x<=0 ) then
          points_before_abs = ii
          x = 0.0d0
          if ( until_extinction ) then
            exit
          end if 
        end if
        t0 = tf
        x0 = x
        y0 = y
        theta = ATAN2(y-1.0d0,x-1.0d0)
        dtheta = theta - theta_old
        dtheta = modulo(dtheta + pi, 2.0d0*pi) - pi
        theta_w = theta_w + dtheta
        theta_old = theta
        thetas_w(ii) = theta_w
        xs(ii) = x
        ys(ii) = y
        ts(ii) = t
      end do
      open(unit=1001, file=trim(name_traj),iostat=ios, status="unknown", action="write")
      print*, "Trajectory saved in: ",trim(name_traj)
      if ( ios /= 0 ) stop "Error opening scratch file on unit 1001"
        write(1001, *) ts(0:points_before_abs)
        write(1001, *) xs(0:points_before_abs)
        write(1001, *) ys(0:points_before_abs)
        if ( print_theta_w_ ) then
          write(1001, *) thetas_w(0:points_before_abs)
        end if
      close(1001)      
      print*, "final t=",t,"xf=",x,"yf=",y

    
    deallocate(ts,xs,ys)
    end subroutine Draw_trajectory

    subroutine trajectory_milstein(t0,x0,y0,tf,dt)
    !Evolve trajectory from (x0,t0) to (x,tf) using milstein method
    !if dt bigger than tf-t0, or dt negative and tf<t0, or dt>0 and t0>td then evolve only one step (t->t+dt).
      implicit none
      double precision,intent (in) :: t0,x0,y0,tf,dt
      double precision :: D2dt,Dsqdt,u
      double precision :: drift_x,drift_y,diff_x,milstein_x
      integer*8 :: steps,ii
      double precision :: dran_g

      t = t0
      x = x0
      y = y0
      D2dt = D*D*dt
      Dsqdt = D*sqrt(dt)
      steps = int8((tf-t)/dt)
      if ( steps<1 ) then
        steps = 1
      end if
      do ii = 1, steps, 1
        t = t+dt
        if ( x<=0.0d0 ) then
          x = 0.0d0
        else
          drift_x = (x - x*y)*dt
          u = dran_g()
          diff_x = Dsqdt*sqrt((x + x*y))*u
          milstein_x = (1.0d0-u**2.0d0)*(1+y)*D2dt/4.0d0 !1/2 * g(x)*g'(x)*dt*(1.0d0-u**2.0d0)
          x = x + drift_x + diff_x + milstein_x 
        end if
        drift_y = r*(x-y)*dt
        y = y + drift_y
        
      end do
    end subroutine trajectory_milstein

    subroutine read_xF_tf()
      !Example of subroutine, read data file
      implicit none
      integer*4 :: ios,i,dum,datapoints
      double precision :: dummy
  
      open(unit=1001, file="data/F_Target_test3.dat", iostat=ios, status="old", action="read")
      if ( ios /= 0 ) stop "Error opening file "
      do i = 1, datapoints, 1
        read(1001,*) dummy,dummy,dummy,dum
      end do
  
  
    end subroutine read_xF_tf

  subroutine logspace(start, end, num, x, endpoint, base)
  implicit none
  real(kind=8), intent(in)  :: start, end
  integer*8,      intent(in)  :: num
  real(kind=8), intent(out) :: x(num)
  logical,      intent(in), optional :: endpoint
  real(kind=8), intent(in), optional :: base

  logical :: endpoint_
  real(kind=8) :: base_, step
  integer :: i

  if (num <= 0) return   ! nothing to do (just exits)

  endpoint_ = .true.
  if (present(endpoint)) endpoint_ = endpoint

  base_ = 10.0d0
  if (present(base)) base_ = base

  if (num == 1) then
    x(1) = base_**start
    return
  end if

  if (endpoint_) then
    step = (end - start) / real(num-1,kind=8)
  else
    step = (end - start) / real(num,kind=8)
  end if

  do i = 1, num
    x(i) = base_**( start + step*real(i-1,kind=8) )
  end do

end subroutine logspace

    
    subroutine search_list_binary_algorithm(list,position,p)
      !Rafle event:
      !Given a list of probabilities called "list" such that sum(list)=1 and a probability p.
      !Look for "position" such that C(position)>=p and C(j)<p for all j in [1,position[.
      !Where C is the cumulative of list: C(i)=list(1)+list(2)+...+list(i)-
      !REFERENCE:Brainerd, W. S. (2015). Guide to Fortran 2008 programming (p. 141). Berlin: Springer.
  
      implicit none
      double precision, dimension(:), intent (in) :: list
      double precision, intent (in) :: p
      integer*4, intent(out) :: position
      double precision, dimension(size(list)) :: C
      integer*4 ii,N,first,last,half
  
      N=size(list) !It would be cool to define this as a parameter (constant), I don't know how...
      C(1)=list(1)
      do ii = 2, N, 1 !Compute cumulative of list
        c(ii)=C(ii-1)+list(ii)
      end do
  
      first=1;last=N
      do while ( first.ne.last )
        half=(first+last)/2
        if ( p>C(half) ) then
          first=half+1
        else
          last=half
        end if
      end do
  
      position=first
  
    end subroutine search_list_binary_algorithm
  
    subroutine float_to_string(a, n, result)
      !Convert float number a to string with format int(a)//"d"//dec(a).
      !E.g. 0.00543--> "0d00543"
      !For dec(a) select first n decimals
      !If number of digits in dec(a)>n then add zeros to the left
      double precision, intent(in) :: a
      integer, intent(in) :: n
      character(len=*), intent(out) :: result
      character(len=n) :: dumc
      integer :: digit,zeros
  
    
      integer :: int_part
      double precision :: dec_part
    
      ! Get the integer and decimal parts
      int_part = int(a)
      dec_part = a - int_part
  
      ! Convert the decimal part to a string with n digits
      zeros=0
      digit=int(dec_part*10**(zeros+1))
      do while ((digit==0).and.(zeros<n))
        zeros=zeros+1
        digit=int(dec_part*10**(zeros+1))
      enddo
      dumc=repeat("0",zeros)//trim(str(int(dec_part*10**(n)))) 
  
  
      result=trim(str(int_part))//"d"//trim(dumc)
  
    
    end subroutine float_to_string
  
  end module subroutines