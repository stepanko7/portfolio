program main
implicit none
integer*8 fftw_plan_r2c_1d, fftw_plan_c2r_1d
integer*8 fftw_plan_r2c_2d, fftw_plan_c2r_2d
integer nx, ny
double precision, allocatable, dimension(:) :: rho_1d, conv_1d, w3_1d
double complex, allocatable, dimension(:) :: w3k_1d, rhok_1d, tmpk ! Yes! Fortran supports complex numbers out of the box!
double precision, allocatable, dimension(:,:) :: rho_2d, conv_2d, w3_2d
double complex, allocatable, dimension(:,:) :: w3k_2d, rhok_2d, tmpk_2d 
double precision :: sgm = 1d0, pi = 3.141592654d0, Lx=5d0, Ly=5d0
double precision sn,cs,kx,ky,k,t
integer i, j, xi, yi


  !Lets start with Fourier 1D-transform!
  
  nx = 100

  allocate( rho_1d(nx), w3_1d(nx), conv_1d(nx), w3k_1d(nx/2+1), rhok_1d(nx/2+1), tmpk(nx/2+1) )
  call initialize_fftw_plans_1d(fftw_plan_r2c_1d, fftw_plan_c2r_1d, nx)    


  !It is also possible to get k-images of "square" functions analytically, however it is confusing,
  !so I have commented the code below.
  
  ! weight function is simple square, i.e. constant in some region and 0 outside!
!  do i=1, nx/2+1 ! prepare the k-space image of the weight function! It is analytical!
!    if( i==1 )then
!      w3k_1d(i) = pi*sgm**3/6
!    else
!      xi = i-1
!      if( xi > nx/2) xi = xi - nx
!      kx = 2*pi*xi/Lx
!      k  = dabs( kx )
!      t = k*sgm/2
!      sn = dsin(t)
!      cs = dcos(t)
!      w3k_1d(i) = dcmplx( (pi*sgm**3/2) * (sn-t*cs)/t**3, 0d0)
!    endif
!  enddo

! Instead, we will get the images by the forward fftw transform
  w3_1d = 0d0        ! Fortran supports vectorized operations. Just like NumPy!
  w3_1d(1:9) = 1d0   ! Slicing is also supported!
  w3_1d(90:100) = 1d0 
  ! forward r->k transform
  call dfftw_execute_dft_r2c(fftw_plan_r2c_1d, w3_1d, w3k_1d) !get k-image of the rho_1d


  !similarly as the weight function, this function is constant inside some region and 0 everywhere else
  rho_1d = 0d0        
  rho_1d(40:59) = 1d0 

  ! forward r->k transform
  call dfftw_execute_dft_r2c(fftw_plan_r2c_1d, rho_1d, rhok_1d) !get k-image of the rho_1d

  tmpk = w3k_1d*rhok_1d ! elementwise multiplication of the images to get the image of the convolution

  ! backward k->r transform
  call dfftw_execute_dft_c2r(fftw_plan_c2r_1d, tmpk, conv_1d)
  conv_1d = conv_1d/Nx ! we have to divide by array size!

  conv_1d = conv_1d/(2*pi**2)
  
  ! as you might see by plotting this file with "python plot.py",
  ! by convoluting two "squares" you will get you a triangle! :)
  print *, "Writing convolution_1d.dat ..."
  open(unit=1, file='convolution_1d.dat', status='replace')
  do i=1, nx
    write(1,fmt='(1i4,3e16.6)') i, rho_1d(i), w3_1d(i), conv_1d(i)
  enddo
  close(1)
!  Please, plot the above file with plot.py script supplied with this source code.
!  simply run:
!  python plot.py

  call destroy_fftw_plans_1d(fftw_plan_r2c_1d, fftw_plan_c2r_1d)
  deallocate( rho_1d, w3_1d, conv_1d, w3k_1d, rhok_1d, tmpk )




  !Now the 2D-convolutions ---------------------------------------------------------------------------
  
  ! We will convolute two functions that are nonzero and constant on a square support region in two dimensions.
  ! Outside the region the functions are 0. 

  ny = nx

  allocate( rho_2d(nx,ny), w3_2d(nx,ny), conv_2d(nx,ny), w3k_2d(nx/2+1,ny), rhok_2d(nx/2+1,ny), tmpk_2d(nx/2+1,ny) )
  call initialize_fftw_plans_2d(fftw_plan_r2c_2d, fftw_plan_c2r_2d, nx, ny)    

  w3_2d = 0d0       
  w3_2d(1:9,1:9) = 1d0         ! functions (box) are assumed to be periodic, so we define the 
  w3_2d(90:100,1:9) = 1d0      ! the square support region of the function to be at the border of the box.
  w3_2d(1:9,90:100) = 1d0      
  w3_2d(90:100,90:100) = 1d0 
  ! forward r->k transform
  call dfftw_execute_dft_r2c(fftw_plan_r2c_2d, w3_2d, w3k_2d) !get k-image of the rho_2d


  !This function is constant (eqaul to 1) inside the square region locted in the center of the box, and 0 everywhere else
  rho_2d = 0d0        
  rho_2d(40:59,40:59) = 1d0 

  ! forward r->k transform
  call dfftw_execute_dft_r2c(fftw_plan_r2c_2d, rho_2d, rhok_2d) !get k-image of the rho_2d

  tmpk_2d = w3k_2d*rhok_2d ! elementwise multiplication of the images to get the image of the convolution

  ! backward k->r transform
  call dfftw_execute_dft_c2r(fftw_plan_c2r_2d, tmpk_2d, conv_2d)
  conv_2d = conv_2d/(Nx*Ny) ! we have to divide by 2D array size, and we are done!

  conv_2d = conv_2d/(4*pi**4)
  ! In the 2D case instead of triangle, you will get something that looks like pyramid,
  ! but decays nonlinearly ... so, ot is not pyramid, but some kid of peak.
  print *, "Writing convolution_2d.dat ..."
  open(unit=1, file='convolution_2d.dat', status='replace')
  write(1,fmt='(2i4)') nx, ny
  do j=1, ny
  do i=1, nx
    write(1,fmt='(3e16.6)') rho_2d(i,j), w3_2d(i,j), conv_2d(i,j)
  enddo
  enddo
  close(1)
!  Please, plot the above file with plot.py script supplied with this source code

  call destroy_fftw_plans_2d(fftw_plan_r2c_2d, fftw_plan_c2r_2d)
  deallocate( rho_2d, w3_2d, conv_2d, w3k_2d, rhok_2d, tmpk_2d )
  
  print *, "Now run:  python plot.py"
  print *, "To plot the output files!"
  
end program



subroutine initialize_fftw_plans_1d(fftw_plan_r2c_1d, fftw_plan_c2r_1d, nx)
  use fftw_const
  implicit none
  integer*8 fftw_plan_r2c_1d, fftw_plan_c2r_1d
  integer error, nx
  double precision, dimension(nx) :: in
  double complex, dimension(nx/2+1) :: out
  integer isuccess, numompthreads, resinitthr, omp_get_max_threads

  write(*,*) 'init 1d plans'
  
  call dfftw_plan_dft_r2c_1d(fftw_plan_r2c_1d,nx,in,out,fftw_estimate)
  call dfftw_plan_dft_c2r_1d(fftw_plan_c2r_1d,nx,out,in,fftw_estimate)
end subroutine

subroutine destroy_fftw_plans_1d(fftw_plan_r2c_1d, fftw_plan_c2r_1d)
  use fftw_const
  implicit none
  integer*8 fftw_plan_r2c_1d, fftw_plan_c2r_1d

  write(*,*) 'destroy 1d plans'
  call dfftw_destroy_plan(fftw_plan_r2c_1d)
  call dfftw_destroy_plan(fftw_plan_c2r_1d)
end subroutine




SUBROUTINE initialize_fftw_plans_2d(fftw_plan_r2c_2d, fftw_plan_c2r_2d, nx, ny)
  use fftw_const
  implicit none
  integer*8 fftw_plan_r2c_2d, fftw_plan_c2r_2d
  integer error, Nx, Ny
  double precision, dimension(nx,ny) :: in
  double complex, dimension(nx/2+1,ny) :: out

  write(*,*) 'init 2d plans'
  call dfftw_plan_dft_r2c_2d(fftw_plan_r2c_2d,Nx,Ny,in,out,FFTW_MEASURE)
  call dfftw_plan_dft_c2r_2d(fftw_plan_c2r_2d,Nx,Ny,out,in,FFTW_MEASURE)
end subroutine


subroutine destroy_fftw_plans_2d(fftw_plan_r2c_2d, fftw_plan_c2r_2d)
  USE fftw_const
  implicit none
  integer*8 fftw_plan_r2c_2d, fftw_plan_c2r_2d

  write(*,*) 'destroy 2d plans'
  call dfftw_destroy_plan(fftw_plan_r2c_2d)
  call dfftw_destroy_plan(fftw_plan_c2r_2d)
end subroutine


