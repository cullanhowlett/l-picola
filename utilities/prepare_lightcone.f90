program prepare_lightcone
  
  ! ============================================================================================================ !
  ! This routine reads in the unformatted binary output from L-PICOLA's lightcone mode and converts it to ASCII  !
  ! It has similar behaviour to 'prepare_gadget.f90', and assumes a contiguously numbered set of unformatted     !
  ! input files, in L-PICOLA's unformatted lightcone output format, with the same base filename, .i.e.,          !
  ! inputfile.0, inputfile.1 inputfile.2, etc. The output ASCII files are named as inputfile_reform.*, with      !
  ! the same numbers as the corresponding input file.                                                            !
  ! Compile with 'ifort prepare_lightcone.f90 -o prepare_lightcone                                               !
  ! ============================================================================================================ !

  ! Input parameters
  ! ================
  character(len=5) :: ninput_str
  character(len=1000) :: inputfile

  ! Other variables
  ! ===============
  character(len=5) :: file_str
  character(len=1000) :: foutname, finname
  logical*1 :: errorflag
  integer*4 :: i, j, EOF, nchunk, ninput_int
  real*4, allocatable, dimension(:,:) :: r, v

  ! Reads the command line arguments
  ! ================================
  if (command_argument_count().eq.2) then
    call get_command_argument(1,inputfile)
    call get_command_argument(2,ninput_str)
  else 
    write(*,*) 'Usage: prepare_lightcone, [Input_Filename, Number_of_Input_Files]'
    STOP
  end if

  read(ninput_str,'(i4.4)') ninput_int

  ! Loop over the input files
  do i = 0, ninput_int-1

    write(file_str,'(i5)') i

    foutname=trim(adjustl(inputfile))//'_reform.'//trim(adjustl(file_str))
    open(10,file=foutname,form='formatted')

    finname=trim(adjustl(inputfile))//'.'//trim(adjustl(file_str))
    open(unit=11,file=finname,form='unformatted')
    print*, 'Opened: '//trim(adjustl(finname))

    ! Loop over each 'chunk' of particles in the file
    errorflag = .false.
    do 
      read(11, IOSTAT=EOF) nchunk
      print*, nchunk
      if (EOF .gt. 0) then
        print*, 'Read error in file: ', finname
        print*, 'Exiting program' 
        errorflag = .true.
      else if (EOF .lt. 0) then
        ! If EOF < 0 then we have read all the chunks
        exit
      else
        ! Reformat the individual particles in the chunk
        allocate(r(3,nchunk))
        allocate(v(3,nchunk))
        read(11)(r(1,j),r(2,j),r(3,j),v(1,j),v(2,j),v(3,j),j=1,nchunk)
        do j = 1, nchunk
          write(10,'(6F13.6)') r(1,j), r(2,j), r(3,j), v(1,j), v(2,j), v(3,j)
        enddo
        deallocate(r, v)
      endif
    enddo
 
    close(10)
    if (errorflag) then 
      exit
    endif

  enddo

end program prepare_lightcone
