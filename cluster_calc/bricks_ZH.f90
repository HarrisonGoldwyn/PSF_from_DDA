program shapemaker

    ! Here is where we declared variables that this shapemaker
    ! will use. 
    integer :: radius
    integer :: x,y,z,diel,idx
    character(len=200)::row

    ! These are the values that determine the character of your 
    ! system. The separation is a the face to face separation 
    ! distance in your dimer. If you want the two spheres to 
    ! have the same composition, change both diel and diel2 to
    ! hold a value of 1. Otherwise, set one at 1 and the other 
    ! at 2 to allow for different materials. 

    diel = 1
    npts = 1
    ! The centers for your spheres are calculated based on the
    ! input parameters. 

    idx = 0

    ! Here we open a temporary file, that will be deleted later. 
    ! This is my rigged way of getting around allocating the size
    ! of an array for our shape
    open(12, file='temp1',status='replace')
    !open(14, file='ddfield.in',status='new')
    !write(14,*) 'w000r000k000.pol1  =name of file with stored polarization'
    !write(14,*) '5.00e-3 = gamma (interaction cutoff parameter'

    ! Right brick
    do x = 1, 9
       do y = -15, 15
          do z = -8, 8
             idx = idx + 1
              write(12,*) idx ,x ,y,z ,diel, diel, diel
             !write(14,*) x, y, z, x, y, z, npts
           end do
       end do
    end do

    ! This closes the temporary file
    close(12)
    !close(14)
    ! We now reopen the file
    open(12, file='temp1',status='old')
    ! We now create our shape file
    open(13, file='shape_temp.dat',status='replace')
    
    ! And write the header of the file
    write(13,*) 'Sphere shape'
    write(13,*) idx, '= number of dipoles in target'
    write(13,*) '1.000000 0.000000 0.000000 = A_1 vector'
    write(13,*) '0.000000 1.000000 0.000000 = A_2 vector'
    write(13,*) '1.000000 1.000000 1.000000 = (d_x,d_y,d_z)/d'
    write(13,*) '0.000000 0.000000 0.000000 = (x,y,z)/d'
    write(13,*) 'JA  IX  IY  IZ ICOMP(x,y,z)'
    
! And here we read the information of the temporary file, and 
! rewrite it after the header.
do x = 1,idx
   read(12,'(a)') row 
   write(13,'(a)') trim(row)
end do

!Here we delete the temporary file and close the shape file.
close(12,status='delete')
close(13)


end program shapemaker
