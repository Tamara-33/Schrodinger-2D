!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       																       !
!     RESOLUCIÓN DE LA ECUACIÓN DE SCHRODINGER DEPENDIENTE DEL TIEMPO          !
!      EN DOS DIMENSIONES, UTILIZANDO EL MÉTODO DEL OPERADOR DIVIDIDO          !
!              Y LA TRANSFORMADA RÁPIDA DE FOURIER (FFT)                       !
!                             GPU-CUFFT                                        !
!                                                                              !
!                                                                              !
!  																		       !
! PROFESOR: ALEJANDRO BENEDYKT KOLTON                                          !
! ALUMNO: TAMARA GUARDA                                                        !
! MATERIA: INTRODUCCIÓN AL CÁLCULO NUMÉRICO EN PROCESADORES GRÁFICOS           !
! FECHA: 27/06/2021                                                              !
! LUGAR: INSTITUTO BALSEIRO-CAB                                                !
! 																			   !
!                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         MÓDULO DE PRECISIÓN (Single)                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module precision
integer, parameter, public :: Single = kind(0.0) ! Precisión simple 
integer, parameter, public :: fp_kind = Single
end module precision

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                            MÓDULO CUFFT                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cufft                                                                                                   
                                                                                                               
integer, public :: CUFFT_FORWARD = -1                                                                          
integer, public :: CUFFT_INVERSE = 1                                                                           
integer, public :: CUFFT_R2C = Z'2a' ! Real to Complex (interleaved)                                           
integer, public :: CUFFT_C2R = Z'2c' ! Complex (interleaved) to Real                                           
integer, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved                                         
integer, public :: CUFFT_D2Z = Z'6a' ! Double to Double-Complex                                                
integer, public :: CUFFT_Z2D = Z'6c' ! Double-Complex to Double                                                
integer, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                          
!                                                                              !                               
! cufftPlan1d(cufftHandle *plan, int nx,cufftType type,int batch)              !                                
!                                                                              !                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                          
                                                                                                                                                                                                                     
interface cufftPlan1d                                                                                          
subroutine cufftPlan1d(plan, nx, type, batch) bind(C,name='cufftPlan1d')                                       
use iso_c_binding                                                                                              
integer(c_int):: plan                                                                                          
integer(c_int),value:: nx, batch,type                
end subroutine cufftPlan1d                                                                                     
end interface cufftPlan1d                                                                                      
                                                                                                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                          
!                                                                              !                               
! cufftPlan2d(cufftHandle *plan, int nx,int ny, cufftType type,int batch)      !                               
!                                                                              !                                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
                                                                                                               
interface cufftPlan2d                                                                                          
subroutine cufftPlan2d(plan, nx,ny, type) bind(C,name='cufftPlan2d')                                         
use iso_c_binding                                                                                            
integer(c_int):: plan                                                                                        
integer(c_int),value:: nx,ny,type                                                                            
end subroutine cufftPlan2d                                                                                   
end interface cufftPlan2d         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                     
!                                                                              !                               
! cufftDestroy(cufftHandle plan)                                               !                                
!                                                                              !                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                     
                                                                                                               
      interface cufftDestroy                                                                                   
subroutine cufftDestroy(plan) bind(C,name='cufftDestroy')                                                      
use iso_c_binding                                                                                              
integer(c_int),value:: plan                                                                                    
end subroutine cufftDestroy                                                                                    
end interface cufftDestroy                                                                                     
                                                                                                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                     
!                                                                              !                               
! cufftExecC2C(cufftHandle plan,                                               !                                
! cufftComplex *idata,                                                         !                               
! cufftComplex *odata,                                                         !                               
! int direction)                                                               !                               
!                                                                              !                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                     
                                                                                                               
interface cufftExecC2C                                                                                         
subroutine cufftExecC2C(plan, idata, odata, direction) &                                                       
& bind(C,name='cufftExecC2C')                                                                                  
use iso_c_binding 
use precision                                                                                             
integer(c_int),value:: direction                                                                           
integer(c_int),value:: plan                                                                                    
complex(fp_kind),device:: idata(*),odata(*)                                                                    
end subroutine cufftExecC2C                                                                                    
end interface cufftExecC2C                                                                                     
                                                                                                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                     
!                                                                              !                                
! cufftExecZ2Z(cufftHandle plan,                                               !                                
! cufftDoubleComplex *idata,                                                   !                                
! cufftDoubleComplex *odata,                                                   !                                
! int direction);                                                              ! 
                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          

interface cufftExecZ2Z                                                                                         
subroutine cufftExecZ2Z(plan, idata, odata, direction) &                                                         
& bind(C,name='cufftExecZ2Z')                                                                                    
use iso_c_binding                                                                                                
use precision                                                                                                    
integer(c_int),value:: direction                                                                                 
integer(c_int),value:: plan                                                                                      
complex(fp_kind),device:: idata(*),odata(*)                                                                      
end subroutine cufftExecZ2Z                                                                                      
end interface cufftExecZ2Z                                                                                                                                                                               
end module cufft             

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                 
!                  KERNEL DE MULTIPLICACIÓN MATRICIAL                          !
!                        ELEMENTO A ELEMENTO                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mathOps                                                                                                   
use precision                                                                                                    
contains                                                                                                         
    attributes(global) subroutine mmul(expv,u0v,fin,npx,npy)                                                       
        implicit none                                                                                                
        complex(fp_kind), dimension(:,:) :: expv, u0v, fin                                                           
        integer, value :: npx,npy                                                                                    
        integer :: i,j                                                                                               
        i = blockDim%x * (blockIdx%x - 1) + threadIdx%x                                                              
        j= blockDim%y * (blockIdx%y - 1) + threadIdx%y                                                               
                                                                                                                 
        if (i <= npx.and.j<=npy) fin(i,j) = expv(i,j)*u0v(i,j)                                                       
                                                                                                                 
    end subroutine mmul                                                                                            
end module mathOps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                 
!       KERNEL DE MULTIPLICACIÓN MATRICIAL ELEMENTO A ELEMENTO                 !
!                  MULTIPLICADO POR UN ESCALAR                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mathOps_3                                                                             
use precision                                                                                                    
contains                                                                                                         
    attributes(global) subroutine mmul_3(expv,fout,unew,npx,npy)                                                   
        implicit none                                                                                                
        complex(fp_kind), dimension(:,:) :: expv, fout, unew                                                         
        integer, value :: npx,npy                                                                                    
        integer :: i,j                                                                                               
        i = blockDim%x * (blockIdx%x - 1) + threadIdx%x                                                              
        j= blockDim%y * (blockIdx%y - 1) + threadIdx%y                                                               
                                                                                                                 
        if (i <= npx.and.j<=npy) unew(i,j) = expv(i,j)*fout(i,j)/(npx*npy)                                           
                                                                                                                 
    end subroutine mmul_3                                                                                          
end module mathOps_3                

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     DEFINIMOS MATRICES EN CPU                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                
module dyn_array                                                                                           
    use precision                                                                                         
    integer :: n
    complex(fp_kind), dimension(:,:), allocatable :: unew,u0v,expv,expt
    complex(fp_kind), dimension(:,:), allocatable :: fin,fout
    real(fp_kind), dimension(:), allocatable :: xgrid,ygrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      DEFINIMOS MATRICES EN GPU                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    complex(fp_kind), dimension(:,:), allocatable, device :: d_fin,d_fout
    complex(fp_kind), dimension(:,:), allocatable, device :: d_unew,d_u0v
    complex(fp_kind), dimension(:,:), allocatable, device :: d_expv,d_expt                                            
end module dyn_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                       MÓDULO PAQUETE DE ONDA                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                          
module wavep                                                                                               
    use precision                                                                                          
    real(fp_kind) sx,sy,x0,y0,kx,ky                                                                            
end module 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                        MÓDULO DE CONSTANTES                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                            
module constants                                                                                           
    use precision                                                                                          
    real(fp_kind), parameter :: a_0=5.291772108d-11,t_0=24.188843265d-18                                       
    real(fp_kind), parameter :: i_0=3.51d16,v_c=137
    real(fp_kind), parameter :: pi=3.141592653589793238,eps0=1d-3
    real(fp_kind), parameter :: ry=27.21                                                          
    integer, parameter :: ix=2  ! ix-times the electron mass                                                   
    real(fp_kind), parameter :: mass=1.0  ! mass=1 au for electrons                                          
    complex(fp_kind), parameter :: ci=(0,1),c1=(1,0),c0=(0,0)                                      
end module                                              

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        INICIALIZAMOS EL PROGRAMA PARA RESOLVER LA TSDE  2D                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                               
 program split2d  
		use mathOps                                                                                        
    use mathOps_3                                                                                          
    use precision                                                                                          
    use cufft                                                                                              
    use cudafor                                                                                            
    use wavep                                                                                              
    use constants                                                                                          
    use dyn_array                                                                                          
    type(dim3)::grid, tBlock                                                                                                             
    real(fp_kind) t,tmin,deltat                                                                                
    real(fp_kind) x,xmin,xmax,deltax,deltakx,deltaky                                                           
    real(fp_kind) y,ymin,ymax,deltay,vpot,v0,ke,kix,kiy                                                        
    complex(fp_kind) u                                                                                         
    integer :: npx,npy                                                                                          
    integer :: nt,i,j,jj,jx,jy,tstep,frame,totframes                                                         
    integer :: planB, planF,allstatus                                                                          
    character*16 filename        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               LLAMAMOS AL POTENCIAL Y AL PAQUETE DE ONDA                     !
!                    DEFINIDO FUERA DEL PROGRAMA                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                                                            
      external vpot,u                                                                                            

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          INICIALIZAMOS LOS PARÁMETROS DEL PAQUETE DE ONDA                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      
                                                                               
      x0=35       !Centro paquete de onda, (parametro de impacto)                                            
      y0=5                                                                                                  
      sx=2        ! Ancho paquete de onda                                                                                 
      sy=2                                                                                                     
      kx=ix*2.0   ! momento                                                                                   
      ky=0        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     DEFINIMOS LA GRILLA ESPACIAL                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                                 
      npx=ix*128                                                                                 
      npy=128 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           ALOCAMOS MEMORIA PARA LAS MATRICES EN CPU                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      allocate(unew(npx,npy),u0v(npx,npy),stat=allstatus)
      allocate(expv(npx,npy),expt(npx,npy),stat=allstatus)
      allocate(fin(npx,npy),fout(npx,npy),stat=allstatus)
      allocate(xgrid(npx),ygrid(npy),stat=allstatus)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           ALOCAMOS MEMORIA PARA LAS MATRICES EN GPU                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
    
      allocate(d_unew(npx,npy),d_u0v(npx,npy),stat=allstatus)
      allocate(d_expv(npx,npy),d_expt(npx,npy),stat=allstatus)
      allocate(d_fin(npx,npy),d_fout(npx,npy),stat=allstatus)

      xmin=-60                  ! size of box                                                                                  
      xmax=60                                                                                                  
      deltax=(xmax-xmin)/npx    !spatial step                                                                    
      deltakx=2*pi/npx/deltax                                                                                    
      ymin=-60                                                                                                 
      ymax=60                                                                                                  
      deltay=(ymax-ymin)/npy                                                                                     
      deltaky=2*pi/npy/deltay                                                                                  
      tBlock = dim3(32,32,1)                                                                                 
      grid = dim3( ceiling(real(npx)/tBlock%x), ceiling(real(npy)/tBlock%y),1)                                   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Potential step height                                                                                          
      v0=5                      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    INSTANTÁNEAS EN EL TIEMPO                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      tstep=2400     !número de pasos de tiempo después de los cuales 
      frame=0        !exportamos los datos
      totframes=100  !total de archivos a imprimir 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   DEFINIMOS PARÁMETROS DE TIEMPO                             !     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      tmin=0                   
      nt=tstep*totframes    !número total de pasos de tiempo
      deltat=1.25d-4        !delta de tiempo entre cada paso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  INICIALIZACIÓN DE LA FUNCIÓN DE ONDA                        !
!                  Y OPERADOR DE ESPACIO DE POSICIÓN "V"                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                   
      do i=1,npx                                                                                                 
      x=xmin+(i-1)*deltax                                                                                        
      xgrid(i)=x                                                                                                 
      do j=1,npy                                                                                                 
      y=ymin+(j-1)*deltay                                                                                        
      ygrid(j)=y 
	    u0v(i,j)=u(x,y)         !llamamos a la función: Función de onda Gaussiana                                                      
      expv(i,j)=exp(-ci*vpot(v0,x,y)*deltat/2d0)                                                          
      enddo                                                                                               
      enddo                                                                                               

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  OPERADOR DE ESPACIO DE IMPULSO "T"                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
 
      do i=1, npx/2
      do j=1, npy/2                                                                                                                                 
      kix=(i-1)*deltakx                                                                                   
      kiy=(j-1)*deltaky                                                                                   
      ke=0.5*(kix**2+kiy**2)/(ix*mass)                                                                  
      expt(i,j)=exp(-ci*ke*deltat)                                                                        
      kix=-i*deltakx                                                                                      
      ke=0.5*(kix**2+kiy**2)/(ix*mass)                                                                  
      expt(npx+1-i,j)=exp(-ci*ke*deltat)                                                                  
      kiy=-j*deltaky                                                                                      
      kix=(i-1)*deltakx                                                                                   
      expt(i,npy+1-j)=exp(-ci*ke*deltat)                                                                  
      kix=-i*deltakx                                                                                      
      expt(npx+1-i,npy+1-j)=exp(-ci*ke*deltat)                                                            
      enddo    
      enddo                                                                                    
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              COPIAMOS LAS MATRICES DE HOST A DEVICE                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     
      d_expt = expt                                                                                       
      d_expv = expv                                                                                   
      d_u0v  = u0v  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      CREAMOS LOS PLANES PARA LA TRANSFORMADA Y LA TRANSFORMADA INVERSA       !
!                                CUFFT                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                           
                                                                                                                                                  
      call cufftPlan2d(planB,npy,npx,CUFFT_C2C)                                                           
      call cufftPlan2d(planF,npy,npx,CUFFT_C2C)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!               LLAMAMOS A LOS KERNELS DE MULTIPLICACIÓN Y                     !
!      CUFFT PARA EJECUTAR LA TRANSFORMADA DIRECTA E INVERSA DE FOURIER        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                  
    do jj=1,nt                                                                                        
      call mmul<<<grid, tBlock>>>(d_expv, d_u0v, d_fin, npx, npy)    !Llamado al kernel de multiplicación                                       
      call cufftExecC2C(planB,d_fin,d_fout,CUFFT_FORWARD)            !Llamado a ejecución planB FFT                                                    
      call mmul<<<grid, tBlock>>>(d_expt, d_fout, d_fin, npx, npy)   !Llamado al kernel de multiplicación                                  
      call cufftExecC2C(planF,d_fin,d_fout,CUFFT_INVERSE)            !Llamado a ejecución planB FFT inversa                                      
      call mmul_3<<<grid, tBlock>>>(d_expv, d_fout, d_unew, npx, npy)!Llamado al kernel de multiplicación     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!           ESCRIBIMOS LOS FRAMES CON LA POSICIÓN EN X,Y,                      !
!          DISTRUBUCIÓN DE PROBABILIDAD UNEW(JX,JY)                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                            
       if(mod(jj-1,tstep).eq.0)then                                                                        
         unew = d_unew                                                                                   
         frame=frame+1         
	       write(*,*) frame                                                                   
         write(filename,'("frame",I3.3,".txt")') frame                                                              
         open(file=filename,unit=9)                                                                                 
         write(9,200) ((xgrid(jx),ygrid(jy),abs(unew(jx,jy))**2,jx=1,npx,ix),jy=1,npy)                              
         close(9)                                                                                                   
       endif                                                                                                      
      d_u0v=d_unew                                                                                                        
    enddo                                                                                                      
200   format(3es23.14e3)                                                                                         
          call cufftDestroy(planB)                                                                               
          call cufftDestroy(planF)                                                                               
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               FUNCIÓN PARA DEFINIR EL PAQUETE DE ONDA INICIAL                ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

complex(fp_kind) function u(x,y)                                                                             
      use wavep                                                                                                  
      use precision                                                                                          
      use constants                                                                                              
      real(fp_kind) x,y,norm                                                                                     
      norm=1/sqrt(2*pi*sx*sy)                                                                               
      u=norm*exp(-((x+x0)**2/(4*sx**2)+(y+y0)**2/(4*sy**2)))*exp(ci*(kx*x+ky*y))                            
      return        
end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            FUNCIÓN PARA DEFINIR EL POTENCIAL CATIÓN DIHIDRÓGENO "H2+"        ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                          
      real(fp_kind) function vpot(v0,x,y)                                                                  
      use precision
      implicit none
      real(fp_kind), x,v0,t,a,vx,vy,y,aux,vxy1,vxy2                                                                   
      real(fp_kind) z,r,a0,a1,a2,r1,b,c                                                                              
      z=-1.0                                                                                                
      r=2.0                                                                                                 
      t=0.01                                                                                                
      a0=1.56                                                                                               
      a1=0.91                                                                                               
      a2=0.9                                                                                                
      b=5.0                                                                                                 
      r1=5.5   
	    a=2                                                                                    
      c = a1+(a0 - a1)*((R1-R)/R1)**3                                                                      
      aux=1.0*x*x +(c/b)**2                                                                              
      vxy1 = z/(1/c-c/b +sqrt((y-R/2)**2 +aux ))                                                    
      vxy2 = z/(1/c-c/b+sqrt((y+R/2)**2+ aux))                                                      
      vpot=vxy1+vxy2                                                                                          
      return                                                                                                  
      end function                
