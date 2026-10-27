module Subs
  integer                           :: Id
  integer    (kind=4), dimension(8) :: vt
  real       (kind=8), dimension(8) :: dv
  character (len=100), dimension(8) :: sv 
  contains
  
  subroutine ExtractWrite(cell)
    use IPhreeqc
    IMPLICIT NONE
    integer :: j
    integer :: Ierr
    integer    (kind=4), intent(in) :: cell
    do j = 1, 8
      ! Headings are on row 0
      Ierr = GetSelectedOutputValue(Id,1,j,vt(j),dv(j),sv(j))
      if(Ierr .ne. IPQ_OK) call EHandler()
    enddo
    write(*,"(a,i2,a,i4,a/2(5x,a,f7.2))") "Cell",cell," ",int(dv(8))," ", &
        "pH:",dv(6),"SR(calcite):",dv(7) 
  end subroutine ExtractWrite
  
  subroutine EHandler()
    use IPhreeqc
    IMPLICIT NONE
    call OutputErrorString(Id)
    stop
  end subroutine EHandler    
end module Subs
    
module MyData
    double precision year
end module MyData    

module Callback    
    contains
    REAL(kind=C_DOUBLE)  function MyCallback(x1, x2, str, l) BIND(C, name='MyCallback')
    USE ISO_C_BINDING
    use MyData, only : year
!       
!   Use of a callback is optional.
!   
!   The callback provides a way to obtain data from a Basic program
!   through the variables x1, x2, and str, and send data to a 
!   Basic program through the return value of the callback.
!   
!   The callback function is called whenever CALLBACK(x1, x2, str$)  
!   is used in a Basic program (usually USER_PUNCH). See file "ic".
!   
        REAL(kind=C_DOUBLE),    INTENT(in)        :: x1, x2
        CHARACTER(kind=C_CHAR), INTENT(in)        :: str(*)
        INTEGER(kind=C_INT),    INTENT(in), value :: l
        character(len=l) fstr
        MyCallback = -1.0
	    ! An example of a datum available in main program    
        do i = 1, l
            fstr(i:i) = str(i)
        enddo
        if (trim(fstr) .eq. "Year") then
            write (*,"(/a,i2,a,f8.2)") "Callback for cell ", int(x1), ": pH ", x2
            MyCallback = year
        endif
        return
    end function MyCallback
end module Callback
    
program Advect
  use Subs
  use Callback
  use MyData
  use IPhreeqc
  IMPLICIT NONE
  integer :: Ierr
  character(len=1024) Istring 
  
!Create module, load database, define initial conditions and selected output
  year = 2012.
  Id = CreateIPhreeqc()
  if (LoadDatabase(Id, "phreeqc.dat") .ne. 0) call EHandler()
  if (SetBasicFortranCallback(id, MyCallback) .ne. 0) call EHandler()
  If (RunFile(Id, "ic") .ne. 0) call EHandler()

!Run cell 1, extract/write result
  if (RunString(Id, "RUN_CELLS; -cells; 1; END") .ne. 0) call EHandler()
  call ExtractWrite(1)

!Advect cell 1 solution to cell 2, run cell 2, extract/write results
  Ierr = AccumulateLine(Id, "SOLUTION_MODIFY 2")    
  Ierr = AccumulateLine(Id, "   -cb      " // sv(1))
  Ierr = AccumulateLine(Id, "   -total_h " // sv(2))
  Ierr = AccumulateLine(Id, "   -total_o " // sv(3))
  Ierr = AccumulateLine(Id, "   -totals  ")
  Ierr = AccumulateLine(Id, "      C     " // sv(4))
  Ierr = AccumulateLine(Id, "      Ca    " // sv(5))
  Ierr = AccumulateLine(Id, "RUN_CELLS; -cells; 2; END")
  year = year + 1
  if (RunAccumulated(Id) .ne. 0) call EHandler()
  call ExtractWrite(2)
  
 !Destroy module 
  if (DestroyIPhreeqc(Id) .ne. 0) call EHandler()
end program Advect