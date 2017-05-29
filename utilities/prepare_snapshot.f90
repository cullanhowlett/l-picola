program prepare_snapshot
  
  ! =========================================================================================================================== !
  ! This example Fortran90 routine reads in data output in GADGET-2's standard 'unformatted' format, which is also shared by    !
  ! the 2LPT code (Manera et al. 2012, Scoccimarro et al. 2012) and L-PICOLA snapshot outputs, and then rewrites them in ASCII  !
  ! format. It assumes a set of files named 'inputfile.*', where the files are contiguously numbered from 0 to 'ninput_str',    !
  ! i.e, inputfile.0, inputfile.1, inputfile.2, etc. The ASCII files are output as inputfile_z?p? reform.*, with the same       !
  ! numbers as the input unformatted files, and where the 'z?p?' is the redshift of the snapshot as defined in the header       !   
  ! Compile with 'ifort prepare_gadget.f90 -o prepare_gadget                                                                    !
  ! =========================================================================================================================== !

  ! The required input parameters
  ! =============================
  character(len=5) :: ninput_str
  character(len=1000) :: inputfile

  ! The header used by GADGET-2, 2LPT and L-PICOLA, for unformatted snapshots
  ! =========================================================================
  type io_gheader
    integer :: npart(0:5)  !should be unsigned integer but not supported by fortran90
    real*8  :: mass(0:5)
    real*8  :: time
    real*8  :: redshift
    integer :: flag_sfr
    integer :: flag_feedback
    integer :: npartTotal(0:5) !should be unsigned integer but not supported by fortran90
    integer :: flag_cooling
    integer :: num_files
    real*8  :: BoxSize   !in Gadget code units
    real*8  :: Omega0
    real*8  :: OmegaL0
    real*8  :: HubbleParam
    integer :: flag_stellarage
    integer :: flag_metals
    integer :: npartTotalHighWord(0:5) !should be unsigned integer but not supported by fortran90
    integer :: flag_entropy_instead_u
    character :: fill*60   !to fill the header
  end type io_gheader

  ! Other variables
  ! ===============
  character(len=2) :: red_str1, red_str2
  character(len=5) :: file_str
  character(len=1000) :: filename
  integer*4 :: i, j, k, ninput_int
  integer*8, allocatable, dimension(:) :: id
  real*4, allocatable, dimension(:,:) :: r, v
  type(io_gheader) :: gheader

  ! Reads the command line arguments
  ! ================================
  if (command_argument_count().eq.2) then
    call get_command_argument(1,inputfile)
    call get_command_argument(2,ninput_str)
  else 
    write(*,*) 'Usage: prepare_gadget, [Input_Filename, Number_of_Input_Files]'
    STOP
  end if

  read(ninput_str,'(i4.4)') ninput_int

  ! Loops over the input files
  ! ==========================
  do i = 0, ninput_int-1
    write(file_str,'(i5)') i
    filename=trim(adjustl(inputfile))//'.'//trim(adjustl(file_str))
    open(unit=11,file=filename,form='unformatted')
    print*, 'Opened: '//trim(adjustl(filename))

    ! Read in the header
    read(11) (gheader%npart(j),j=0,5),(gheader%mass(j),j=0,5),gheader%time,gheader%redshift, &
       & gheader%flag_sfr,gheader%flag_feedback,(gheader%npartTotal(j),j=0,5),gheader%flag_cooling,& 
       & gheader%num_files,gheader%BoxSize,gheader%Omega0,gheader%OmegaL0,gheader%HubbleParam, &
       & gheader%flag_stellarage,gheader%flag_metals,(gheader%npartTotalHighWord(j),j=0,5), &
       & gheader%flag_entropy_instead_u,gheader%fill
    print*, gheader%redshift

    ! Read in the particle data
    allocate(r(3,gheader%npart(1)))
    allocate(v(3,gheader%npart(1)))
    allocate(id(gheader%npart(1)))
    read(11)((r(j,k),j=1,3),k=1,gheader%npart(1))
    read(11)((v(j,k),j=1,3),k=1,gheader%npart(1))
    !read(11)(id(j),j=1,gheader%npart(1)) 
    write(red_str1,'(i1)') int(gheader%redshift)
    write(red_str2,'(i1)') nint(10.0*(gheader%redshift-int(gheader%redshift)))
    close(11)
    
    ! Output the data in ASCII
    filename=trim(adjustl(inputfile))//'_z'//trim(adjustl(red_str1))//'p' & 
            &   //trim(adjustl(red_str2))//'_reform.'//trim(adjustl(file_str))
    open(10,file=filename,form='formatted')
    do j = 1, gheader%npart(1)
      write(10,'(6F15.6)') r(1,j), r(2,j), r(3,j), v(1,j), v(2,j), v(3,j)
      !write(10,'(I15, 6F15.6)') id(j), r(1,j), r(2,j), r(3,j), v(1,j), v(2,j), v(3,j)
    end do

    deallocate(id, r, v)
    close(10)
  end do

end program prepare_snapshot
