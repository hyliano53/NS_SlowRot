program TOV

implicit none 
      
! Constant variables
double precision pi,mevtokmm2,kmpersunsmass,k,NuR,cspeed,xi0,p2,xi2,rp,re,ellip,Jtot,deltam

! Variables for the Runge-Kutta iteration:
double precision ka1,ka2,ka3,ka4,kb1,kb2,kb3,kb4,kc1,kc2,kc3,kc4,h,kd1,kd2,kd3,kd4,ke1,ke2,ke3,ke4
double precision kf1,kf2,kf3,kf4,kg1,kg2,kg3,kg4,kh1,kh2,kh3,kh4,ki1,ki2,ki3,ki4
double precision Maux,Edensaux2,Paux,drhodpaux,dmdraux,jinterpol,derjinterpol,aaux,waux,dnudraux

! Arrays to contain the dynamical functions (pressure, mass, density, r)
integer i,j,N,Nit,len
parameter(N=10000)
parameter(Nit=2000)
double precision M(N),dmdr(N),Edens(N),P(N),r(N),Nu(N),dnudr(N),w(N),a(N),jjj(N),derj(N),drhodp(N),m0(N),p0(N),h2(N),v2(N)
      
! Function for the equation of state
double precision Edens_interpol,Edensaux,badcs

external Edens_interpol
common/unitconversion/mevtokmm2
common/pinumber/pi
pi=3.14159265d0

! A file_to_save the mass of radius relation
open(unit=7,file='Mass_of_Radius.dat')
open(unit=17,file='Mass_Profile.dat')
open(unit=19,file='Pressure_Profile.dat')
open(unit=21,file='Nu_Profile.dat')
open(unit=23,file='Omega_Profile.dat')
open(unit=24,file='testing.dat')
open(unit=30,file='m0_profile.dat')
open(unit=31,file='p0_profile.dat')
open(unit=40,file='w_profile.dat')
open(unit=50,file='h2_profile.dat')
open(unit=51,file='v2_profile.dat')


! Unit conversion from MeV/fm3 to km**-2
mevtokmm2=1.3152d-6

! Unit coversion from s to km
cspeed=2.99792d5

! Unit conversion from km to solar masses
kmpersunsmass=1.47d0

! The energy density for which the speed of sound exceeds the light speed
! in the effective theory. For the free neutron gas this does not occur.
badcs=1d20*mevtokmm2

! Prepare a table rho(P) from a file
call init_eqofstate()

! Initialize the radial variable from 0 to 50 kilometers
h=5d1/dble(N)
do i=1,N
	r(i)=dble(i)*h
enddo

!Initialize the mass to zero (boundary condition)
M(1)=0d0
! Initialize the nu fuction to zero (boundary condition)
Nu(1)=0d0


! The typical cetral pressure of a NS is about hundreds of MeV/fm^3. We need it in km**-2
      P(1)=140*mevtokmm2
	Edens(1)=Edens_interpol(P(1))
! Auxiliary iteration where_we integrate out in r
	do i=2,N
            len=i
		call auxders(r(i-1),Edens(i-1),M(i-1),P(i-1),ka1,kb1,kc1)
        dmdr(i-1)=ka1
        dnudr(i-1)=kc1
! Construct the Runge-Kutta derivatives
		Edensaux=Edens_interpol(P(i-1)+h*kb1/2d0)
        call auxders(r(i-1)+h/2d0, Edensaux, M(i-1)+h/2d0*ka1,P(i-1)+h/2d0*kb1, ka2, kb2, kc2)
        Edensaux=Edens_interpol(P(i-1)+h*kb2/2d0)
        
        call auxders(r(i-1)+h/2d0, Edensaux, M(i-1)+h/2d0*ka2,P(i-1)+h/2d0*kb2, ka3, kb3, kc3)
        Edensaux=Edens_interpol(P(i-1)+h*kb3)
        
        call auxders(r(i-1)+h, Edensaux, M(i-1)+h*ka3,P(i-1)+h*kb3, ka4, kb4, kc4)
! advance_ the pressure with Runge-Kutta
        P(i)=P(i-1)+(h/6d0)*(kb1+2d0*kb2+2d0*kb3+kb4)
        write(19,*) r(i), P(i)/mevtokmm2
! advance_the mass with Runge-Kutta
        M(i)=M(i-1)+(h/6d0)*(ka1+2d0*ka2+2d0*ka3+ka4)
		write(17,*) r(i), M(i)/kmpersunsmass
! advance_the Nu fuction with Runge-Kutta
        Nu(i)=Nu(i-1)+(h/6d0)*(kc1+2d0*kc2+2d0*kc3+kc4)
        
        Edens(i)=Edens_interpol(P(i))
        drhodp(i-1)=(Edens(i)-Edens(i-1))/(P(i)-P(i-1))

! Conditions to break from the loop: reached the_end of the star
! or reached an energy density for which causality breaks down
        if(P(i).le.5d-10) then 
        	goto 130
        endif

        if(Edens(i).gt.badcs) then
        	print*, 'This point discarded: causality breach, P high'
        	print*, 'Energy density =',Edens(i),'Pressure=',P(i)
        	goto 140 
        endif
	enddo
	
! Keep in the disk the mass/suns mass and radius in km**-2
 130   write(7,*) i, r(i),M(i)/kmpersunsmass
       print*, 'end results ',i-1,r(i),M(i)/kmpersunsmass


 140 continue

! len is the length of r and M vectors
 drhodp(len)=(Edens(len)-Edens(len-1))/(P(len)-P(len-1))


! Initialize the M derivative at the edge of the star
 dmdr(len)=4d0*pi*r(len)**2*Edens(len)
 dnudr(len)=-(2d0*kb4)/(Edens(len)+P(len))


 NuR=Nu(len)
 k=M(len)/r(len) ! Compactness at the star's edge
do i=1,len
      ! Modify the Nu function to have the proper boundary condition, matching with the Swarzschild sol in R
      Nu(i)=Nu(i)-NuR+dlog(1d0-2d0*k)
      write(21,*) r(i), Nu(i)
      ! Save the j function and its derivative (necessary to solve the rotating star)
      jjj(i)=dexp(-Nu(i)/2d0)*dsqrt(1d0-2d0*M(i)/r(i))
      derj(i)=-4d0*pi*r(i)*(Edens(i)+P(i))*dexp(-Nu(i))/jjj(i)
enddo


! Keplerian frequency
print*, 'The max ang velocity is (ms^-1)', 240d0*dsqrt((M(len)/kmpersunsmass)/(r(len)**3))


! Now we solve the rotation
 a(1)=0d0
 w(1)=0.003
 m0(1)=0d0
 p0(1)=0d0
 h2(1)=0d0
 v2(1)=0d0

 do i=2,len
    ! First of all we calculate the frequency
    call auxders2(r(i-1),jjj(i-1),derj(i-1),a(i-1),w(i-1),kd1,ke1)
    jinterpol=(h/2d0)*(jjj(i)-jjj(i-1))/(r(i)-r(i-1))+jjj(i-1)
    derjinterpol=(h/2d0)*(derj(i)-derj(i-1))/(r(i)-r(i-1))+derj(i-1)
    call auxders2(r(i-1)+h/2d0,jinterpol,derjinterpol,a(i-1)+h/2d0*kd1,w(i-1)+h/2d0*ke1,kd2,ke2)
    call auxders2(r(i-1)+h/2d0,jinterpol,derjinterpol,a(i-1)+h/2d0*kd2,w(i-1)+h/2d0*ke2,kd3,ke3)
    call auxders2(r(i-1)+h,jjj(i),derj(i),a(i-1)+h*kd3,w(i-1)+h*ke3,kd4,ke4)
    a(i)=a(i-1)+(h/6d0)*(kd1+2d0*kd2+2d0*kd3+kd4)
    write(24,*) r(i), a(i)*cspeed**2
    w(i)=w(i-1)+(h/6d0)*(ke1+2d0*ke2+2d0*ke3+ke4)
    write(23,*) r(i), w(i)*cspeed



    ! Now we solve the spherical deformations l=0
    call auxders3(r(i-1),M(i-1),Edens(i-1),P(i-1),drhodp(i-1),dmdr(i-1),jjj(i-1),derj(i-1),a(i-1),w(i-1),m0(i-1),p0(i-1),kf1,kg1)
    Maux=M(i-1)+h/2d0*dmdr(i-1)
    Edensaux2=(h/2d0)*(Edens(i)-Edens(i-1))/(r(i)-r(i-1))+Edens(i-1)
    Paux=(h/2d0)*(P(i)-P(i-1))/(r(i)-r(i-1))+P(i-1)
    drhodpaux=(h/2d0)*(drhodp(i)-drhodp(i-1))/(r(i)-r(i-1))+drhodp(i-1)
    dmdraux=(h/2d0)*(dmdr(i)-dmdr(i-1))/(r(i)-r(i-1))+dmdr(i-1)
    aaux=(h/2d0)*(a(i)-a(i-1))/(r(i)-r(i-1))+a(i-1)
    waux=w(i-1)+h/2d0*a(i-1)
    call auxders3(r(i-1)+h/2d0,Maux,Edensaux2,Paux,drhodpaux,dmdraux,jinterpol,derjinterpol,aaux,waux, &
    m0(i-1)+h/2d0*kf1,p0(i-1)+h/2d0*kg1,kf2,kg2)
    call auxders3(r(i-1)+h/2d0,Maux,Edensaux2,Paux,drhodpaux,dmdraux,jinterpol,derjinterpol,aaux,waux, &
    m0(i-1)+h/2d0*kf2,p0(i-1)+h/2d0*kg2,kf3,kg3)
    call auxders3(r(i-1)+h,M(i),Edens(i),P(i),drhodp(i),dmdr(i),jjj(i),derj(i),a(i),w(i), &
    m0(i-1)+h*kf3,p0(i-1)+h*kg3,kf4,kg4)

    m0(i)=m0(i-1)+(h/6d0)*(kf1+2d0*kf2+2d0*kf3+kf4)
    write(30,*) r(i), m0(i)/kmpersunsmass
    p0(i)=p0(i-1)+(h/6d0)*(kg1+2d0*kg2+2d0*kg3+kg4)
    write(31,*) r(i), p0(i)



    ! Now we solve the l=2 equations
    call auxders4(r(i-1),M(i-1),Edens(i-1),P(i-1),drhodp(i-1),dnudr(i-1),jjj(i-1),derj(i-1),a(i-1),w(i-1),h2(i-1),v2(i-1),kh1,ki1)
    dnudraux=(h/2d0)*(dnudr(i)-dnudr(i-1))/(r(i)-r(i-1))+dnudr(i-1)
    call auxders4(r(i-1)+h/2d0,Maux,Edensaux2,Paux,drhodpaux,dnudraux,jinterpol,derjinterpol,aaux,waux, &
    h2(i-1)+h/2d0*kh1,v2(i-1)+h/2d0*ki1,kh2,ki2)
    call auxders4(r(i-1)+h/2d0,Maux,Edensaux2,Paux,drhodpaux,dnudraux,jinterpol,derjinterpol,aaux,waux, &
    h2(i-1)+h/2d0*kh2,v2(i-1)+h/2d0*ki2,kh3,ki3)
    call auxders4(r(i-1)+h,M(i),Edens(i),P(i),drhodp(i),dnudr(i),jjj(i),derj(i),a(i),w(i), &
    h2(i-1)+h*kh3,v2(i-1)+h*ki3,kh4,ki4)

    h2(i)=h2(i-1)+(h/6d0)*(kh1+2d0*kh2+2d0*kh3+kh4)
    write(50,*) r(i), h2(i)
    v2(i)=v2(i-1)+(h/6d0)*(ki1+2d0*ki2+2d0*ki3+ki4)
    write(51,*) r(i), v2(i)
 enddo

  ! We calculate xi0, p2 and xi2
 xi0=-p0(len)*(Edens(len)+P(len))/kb4 ! kb1 is the last value that dpdr takes in the TOV loop (dpdr in R)
 p2=-h2(len)-r(len)**2*w(len)**2*dexp(-nu(len))/3d0
 xi2=-p2*(Edens(len)+P(len))/kb4

! We calculate the ellipticity
 rp=r(len)+xi0+xi2
 re=r(len)+xi0-xi2/2d0
 ellip=dsqrt(1d0-(rp/re)**2)

 print*, 'La excentricidad es:', ellip


print*, 'Omega (s^-1):', (w(len)+r(len)*a(len)/3)*cspeed
print*, 'Omega (ms^-1):', (w(len)+r(len)*a(len)/3)*cspeed/1000d0
print*, r(len), w(len), (w(len)+r(len)*a(len)/3)*r(len)


!This loop calculates the angular velocity \omega(r)
do i=2,len
      write(40,*) r(i)/r(len), ((w(len)+r(len)*a(len)/3)-w(i))/(w(len)+r(len)*a(len)/3)
      !write(40,*) r(i), (r(i)*a(i)/3)
enddo


! We compute the total angular momentum and the total mass
Jtot=r(len)**4*a(len)/6
deltam=m0(len)+Jtot**2/r(len)**3

print*, 'La masa total es ', (M(len)+deltam)/kmpersunsmass



close(7)
close(17)
close(19)
close(21)
close(23)
close(24)
close(30)
close(31)
close(40)
close(50)
close(51)
end




! SUBROUTINES
subroutine auxders(r,Edens,M,P,derM,derP,derNu)
!Calculates the derivatives of the pressure, the mass ad the Nu function
      implicit none
      double precision pi,r,Edens,M,P,derM,derP,derNu
      common/pinumber/pi
      derM=4d0*pi*r**2*Edens
      derP=-(1d0/r**2)*(Edens+P)*(M+4d0*pi*P*r**3)
      derP=derP/(1d0-2d0*M/r)
      derNu=-(2d0*derP)/(Edens+P)
      return
      end

subroutine auxders2(r,jjj,derj,a,w,dera,derw)
!Calculates the derivatives of the frequency and its derivative
      implicit none
      double precision pi,r,jjj,derj,a,w,dera,derw
      common/pinumber/pi
      dera=-derj/jjj*(a+4d0*w/r)-4d0*a/r
      derw=a
      return
      end


subroutine auxders3(r,M,Edens,P,drhodp,dmdr,jjj,derj,a,w,m0,p0,derm0,derp0)
!Calculates the derivatives of the mass and the pressure perturbation factors
      implicit none
      double precision pi,r,M,Edens,P,drhodp,dmdr,jjj,derj,a,w,m0,p0,derm0,derp0,K
      common/pinumber/pi
      derm0=4d0*pi*r**2*(Edens+P)*drhodp*p0+(1d0/12d0)*r**4*jjj**2*a**2-(2d0/3d0)*jjj*r**3*w**2*derj
      K=(r**3*jjj**2*w**2/(r-2d0*M)**2)*((r-2d0*M)*(3d0/r+2d0*derj/jjj+2d0*a/w)+2d0*dmdr-1d0)
      derp0=-m0*(1d0+8d0*pi*r**2*P)/(r-2d0*M)**2-(4d0*pi*p0*r**2*(Edens+P)/(r-2d0*M))+(r**4*jjj**2*a**2)/(12d0*(r-2d0*M))+K/3d0
      return
      end


subroutine auxders4(r,M,Edens,P,drhodp,dnudr,jjj,derj,a,w,h2,v2,derh2,derv2)
!Calculates the derivatives of the h2 and the v2 functions
      implicit none
      double precision pi,r,M,Edens,P,drhodp,dnudr,jjj,derj,a,w,h2,v2,derh2,derv2,aa,bb,cc,dd
      common/pinumber/pi
      aa=h2*(-dnudr+r*(8d0*pi*(Edens+P)-4d0*M/r**3)/((r-2d0*M)*dnudr))
      bb=4d0*v2/(r*(r-2d0*M)*dnudr)
      cc=(r*dnudr/2d0-1d0/((r-2d0*M)*dnudr))*(r**3*jjj**2*a**2)/6d0
      dd=(r*dnudr/2d0+1d0/((r-2d0*M)*dnudr))*(2*r**2*jjj*derj*w**2)/3d0
      derh2=aa-bb+cc-dd
      derv2=-dnudr*h2+(1d0/r+dnudr/2d0)*(r**4*jjj**2*a**2/6d0-2d0*jjj*r**3*w**2*derj/3d0)
      return
      end


subroutine init_eqofstate()
      implicit none
      integer i,N
      parameter(N=121)
      double precision rhoinit(N),Pinit(N),mevtokmm2
      common/eqstate/rhoinit,Pinit
      common/unitconversion/mevtokmm2
      open(unit=5,file='P_of_rho.dat')

      do i=1,N
        read(5,*) rhoinit(i),Pinit(i)  
        !Convert to the correct units, from MeV/fm3 to km**-2
        rhoinit(i)=rhoinit(i)*mevtokmm2
        Pinit(i)=Pinit(i)*mevtokmm2
      enddo
      close(5)
      return
      end



double precision function Edens_interpol(P)
! This is to obtain the energy density from the pressure for a loaded
! equation of state in a table
      implicit none
      double precision rho,P,mevtokmm2
      integer N,j,jstart,jup,jm
      parameter(N=121)
      
      double precision rhoinit(N),Pinit(N)
      common/eqstate/rhoinit,Pinit
      common/unitconversion/mevtokmm2
      if(P.LT.Pinit(1)) then
! Linear interpolation between the last points
       rho=rhoinit(1)/Pinit(1)*P
       goto 200
      endif
      if(P.GT.Pinit(N)) then
! Relativistic equation of state for extrapolating
       rho=3d0*P
       goto 200
      endif
      jstart=1
      jup=N
 270  if((jup-jstart).gt.1) then
        jm=(jup+jstart)/2
        if(P.GT.Pinit(jm)) then
          jstart=jm
        else
          jup=jm
        endif   
        goto 270
      endif   
      j=jstart
      rho=rhoinit(j)+(P-Pinit(j))*(rhoinit(j+1)-rhoinit(j))/(Pinit(j+1)-Pinit(j))
 200  Edens_interpol=rho
      return 
      end


double precision function Const_density(P)
      implicit none
      double precision P
! This equation of state is just a constant density for neutron matter
! We employ the Gaussian system of units, so that energy density is erg/cm**3
      Const_density=1E15
      return
      end
