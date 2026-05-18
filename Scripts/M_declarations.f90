module declarations_module
    !Module for public variables to be used in the rest of the modules
    implicit none
    public
    
    integer*4 :: ios
    double precision :: x,y,t,D,r
  
  contains
  
    subroutine assingments()
      implicit none
      ! Here you read and assign values
    end subroutine assingments
  
  end module declarations_module
