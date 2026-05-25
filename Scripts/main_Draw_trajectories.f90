program main_program
    !Main program, depends on the rest of modules (declarations_module, functions and subroutines))
    use declarations_module
    use functions
    use subroutines
    implicit none
    double precision :: dt

    call assingments()
    call dran_ini(time())
 
    dt = 0.001d0
    ! r = 0.08d0
    ! D = 0.4d0
    r =   1.000000000000000      
    D = 0.14384498882876628    
    call Draw_trajectory(until_t=100.0d0,init_t=0.0d0,dt=dt,init_x=0.01d0,init_y=0.0d0,points=int8(10000),name_traj = "filename.dat",until_extinction=.false.,print_theta_w=.true.)
    
  101 end program main_program
  
