program writedenshell
implicit none

integer, parameter :: Ngmin = 19
integer, parameter :: Ngmax = 21
integer, parameter :: Ntmin = 1
integer, parameter :: Ntmax = 24
character(300), parameter :: view = "topview"

  real, parameter :: ozonesizes(23) = (/  5.27, 14.4, 41.0, & 
                                          2.67, 7.76, 21.5, & 
                                          13.0, 25.7, 72.5, &
                                          4.48, 13.2, 36.7, &
                                          5.60, 15.6, 0.00, &
                                          5.90, 15.6, 0.00, &
                                          5.51, 5.34, 5.55, &
                                          14.6, 44.1 /)  

  real, parameter ::  zonesizesX(23) = (/ 8.00, 14.4, 30.0, & 
                                          10.0, 20.0, 30.0, & 
                                          8.00, 14.0, 30.0, &
                                          10.0, 14.0, 30.0, &
                                          12.0, 14.0, 0.00, &
                                          10.0, 14.0, 0.00, &
                                          8.00, 8.00, 8.00, &
                                          14.0, 30.0 /)  

  real, parameter ::  zonesizesY(23) = (/ 8.00, 14.4, 30.0, & 
                                          10.0, 20.0, 30.0, & 
                                          8.00, 14.0, 30.0, &
                                          10.0, 14.0, 30.0, &
                                          12.0, 14.0, 0.00, &
                                          10.0, 14.0, 0.00, &
                                          8.00, 8.00, 8.00, &
                                          14.0, 30.0 /)  

  real, parameter ::  zonesizesZ(23) = (/ 8.00, 14.4, 30.0, & 
                                          10 .0, 20.0, 30.0, & 
                                          8.00, 14.0, 30.0, &
                                          10.0, 14.0, 30.0, &
                                          12.0, 14.0, 0.00, &
                                          10.0, 14.0, 0.00, &
                                          8.00, 8.00, 8.00, &
                                          14.0, 30.0 /)  


character(300) :: command
character(300) :: prefix
character(300) :: sout
character(300) :: rhofile
character(300) :: ionfile
character(300) :: tempfile
character(300) :: line


integer*8 :: lun

real :: zonesize(3)
integer :: i,j

  open(unit=10,file="density_test.sh")

  command = "./density_test"
  prefix = "~/DwarfGalAnalysis" 
 
  101 format(A,"/DwarfGals/set2-gal",I1,"_",I3.3,".1")
  102 format(A,"/DwarfGals/set2-gal",I2,"_",I3.3,".1")

  201 format(A,"/",A,"/set2-gal",I1,"_",I3.3,".rho")
  202 format(A,"/",A,"/set2-gal",I2,"_",I3.3,".rho")

  301 format(A,"/",A,"/set2-gal",I1,"_",I3.3,".ion")
  302 format(A,"/",A,"/set2-gal",I2,"_",I3.3,".ion")

  401 format(A,"/",A,"/set2-gal",I1,"_",I3.3,".temp")
  402 format(A,"/",A,"/set2-gal",I2,"_",I3.3,".temp")

  500 format(A," ",A," ",A," ",A," ",A)

  do i = Ngmin,Ngmax

        if (i.ge.0 .and. i.le.9) then
           write(sout,101) trim(prefix), i, 1
        else if (i.ge.10 .and. i.le.99) then
           write(sout,102) trim(prefix), i, 1
        end if
 

        if (i==3 .or. i==6 .or. i==15 .or. i==18) cycle


     do j = Ntmin,Ntmax

        if (i.ge.0 .and. i.le.9) then
           write(sout,101) trim(prefix), i, j
           write(rhofile,201) trim(prefix), trim(view), i, j
           write(ionfile,301) trim(prefix), trim(view), i, j
           write(tempfile,401) trim(prefix), trim(view), i, j
        else if (i.ge.10 .and. i.le.99) then
           write(sout,102) trim(prefix), i, j
           write(rhofile,202) trim(prefix), trim(view), i, j
           write(ionfile,302) trim(prefix), trim(view), i, j
           write(tempfile,402) trim(prefix), trim(view), i, j
        end if

        zonesize = (/ zonesizesX(i), zonesizesY(i), zonesizesZ(i) /)

        write(line,500) trim(command), trim(sout), &
                        trim(rhofile), trim(ionfile), trim(tempfile)

        write(10,'(A,3F15.5)') trim(line), zonesize

     end do
  end do




end program writedenshell
