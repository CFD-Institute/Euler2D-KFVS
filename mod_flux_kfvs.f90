module mod_flux_kfvs
  
contains

  function fluxp(rho,ux,uy,e,p,t,a,b,nx,ny) 
    
    implicit none

    real(8), intent(in) :: rho,ux,uy,e,p,t,a,b,nx,ny
    real(8), dimension(1:4) :: fluxp

    real(8) :: un,ut,f1,f2,f3,f4

    un=nx*ux+ny*uy
    ut=-ny*ux+nx*uy

    if (un>=b) then
       
       f1=rho*un
       f2=rho*un**2+p
       f3=ut*f1
       f4=(e+p)*un

    elseif (un<=-b) then

       f1=0.0d0
       f2=0.0d0
       f3=0.0d0
       f4=0.0d0

    else

       f1=2.0d0*a*b**2*(un+b)**2
       f2=4.0d0/3.0d0*a*b**2*(un+b)**3
       f3=ut*f1
       f4=a*b**2/2.0d0*(un+b)**2*(un**2+2.0d0*un*b+2*ut**2+7.0d0/3.0d0*b**2)

    end if

    fluxp(1)=f1
    fluxp(2)=nx*f2-ny*f3
    fluxp(3)=ny*f2+nx*f3
    fluxp(4)=f4

  end function fluxp

  function fluxm(rho,ux,uy,e,p,t,a,b,nx,ny) 
    
    implicit none

    real(8), intent(in) :: rho,ux,uy,e,p,t,a,b,nx,ny
    real(8), dimension(1:4) :: fluxm

    fluxm=-fluxp(rho,ux,uy,e,p,t,a,b,-nx,-ny) 

  end function fluxm
  

end module mod_flux_kfvs


