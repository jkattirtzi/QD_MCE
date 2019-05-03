
!%      Author: John Kattirtzi, Date 23/04/09
!%                                                      Module Integration
!%      
!%      This module deals with the Integration part of the program. 
!%      It will consist of two subroutines:
!%      1. Derivs routine- a routine that will give the derivatives of
!%              the basis set      
!%      2. RK4ccs- a 4th order Runge-kutta subroutine that will use the
!%      type made in the bset module
!%      The rk4 routine in numerical recipes uses assert_eq function to
!%      check that the arrays are the same size- can't copy this easily
!%      so might use an IF statement instead
!%      
!%      Links with other modules:
!%      Need to link with BSET to understand type and overload the
!%      operators
!%      Need to link with program/module that calls the integration
!%      Need to link with GenHam to get the Hord and DHord that is needed 
!%        
!%      Global parameters for this module:
!%      Img=0.0d0,1.0d0
!%********************************************************************************************************************!!
  Module Integration
    !use BSET, only: csbasisfn,ndim,nbf
    use BSET2!, only: operator(+) , operator(*)  
    use GenHam
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Global parameters:
    complex(kind=8),PARAMETER::Img=cmplx(0.0d0,1.0d0)
    contains

!$********************************************************************************************************************!!
!$                                                      Subroutine ReadTime
!$      This subroutine will read in the time parameters:
!$      t0- initial time
!$      tmax- maximum time you want the trajectory to run for
!$      dt- time increment
!$********************************************************************************************************************!!
  
      Subroutine ReadTime(Line1,t0,tmax,dt)
        real(kind=8)::t0
        real(kind=8)::tmax
        real(kind=8)::dt
        character(LEN=100)::Line1
        integer::ierr=0
        OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierr)
        read(30,*,iostat=ierr)LINE1 ! unit 30 was chosen randomly
        If (ierr .ne. 0) then
           Print *, 'WARNING: Not reading Time in Int mod'
           Print *, 'ierr', ierr
        End If
        do while (ierr==0)
           if(Line1 == "t0") then
              backspace(30)
              read(30,*,iostat=ierr)Line1,t0
              If (ierr .ne. 0) then
                 Print *, 'Error reading in t0 '
                 Print *, 'ierr', ierr
              End If
           else if(Line1 == "tmax") then
              backspace(30)
              read(30,*,iostat=ierr)Line1,tmax
              Print *, 'tmax'
              If (ierr .ne. 0) then
                 Print *, 'Error reading in tmax '
                 Print *, 'ierr', ierr
              End If
           else if(Line1 == "dt") then
              backspace(30)
              read(30,*,iostat=ierr)LINE1,dt
              If (ierr .ne. 0) then
                 Print *, 'Error reading in dt '
                 Print *, 'ierr', ierr
              End If
           end if

           read(30,*, iostat=ierr)LINE1
        end do
        close(30)
      End Subroutine ReadTime


!$********************************************************************************************************************!!
!$                                                      Subroutine Derivs
!$      This subroutine will take in the basis set and give out the derivative
!$      of the z,s and ds. 
!$      The arguments are:
!$      t- time
!$,     bs- basis set(input)
!$      dbsdt- derivative of basis set(output)
!$       
!$      Equation for zdot:
!$      zdot=-i/hbar*dhord(zic,zi)/dzic
!$      Equation for Sdot=L
!$      Sdot=L=ihbar/2.0sum_ndim((zc*zdot-zdotc*z)-Hord(zc,z))
!$
!$********************************************************************************************************************!!
 Subroutine Derivs(t,bs,dbsdt) 
   Implicit none
   real(kind=8)::t
   type(csbasisfn),dimension(:)::bs,dbsdt
   type(csbasisfn)::bf1,bf1cz,bf2,dbf1,dbf1cz
   complex(kind=8)::z,zc,zdt,zcdt,d
   complex(kind=8),dimension(nconf,nconf)::HMCE
   integer::i,j,k,m
   call INBASISFN(bf1)
   call INCONF(bf1)
   call INBASISFN(bf1cz)
   call INCONF(bf1cz)
   call INBASISFN(bf2)
   call INCONF(bf2)
   call INBASISFN(dbf1)
   call INCONF(dbf1)
   call INBASISFN(dbf1cz)
   call INCONF(dbf1cz)
   !z derivs
   do i=1,nbf
      !print *, 'nbf', nbf
      print *, 'bsi(1)%z', bs(i)%z
      bf1=bs(i)
      print *, 'bf1%z', bf1%z
      print *, 'bf1cz%z', bf1cz%z
      call bfconz(bf1,bf1cz)
      print *, 'bf1%z', bf1%z
      print *, 'bf1cz%z', bf1cz%z
      !STOP
      m=0
      do m=1,ndim
         dbsdt(i)%z(m)=-1.0d0*Img*DHehrNew(bf1cz,bf1,m)
      end do
   end do  
      print *, dbsdt(nbf)%z
   !print *, 'dbsdtZ', dbsdt(nbf)%z
   !STOP
   !s derivs
   i=0
   bf1%z=0.0d0
   bf1cz%z=0.0d0
   HMCE(:,:)=(0.d0,0.0d0)
   do i=1,nbf 
      zdt=product(dbsdt(i)%z(1:ndim))
      zcdt=dconjg(zdt)
      z=product(bs(i)%z(1:ndim))
      zc=dconjg(z)
      bf1=bs(i)
      call bfconz(bf1,bf1cz)
      call H_MCESpin(bf1cz,bf1,HMCE)
      do j=1,nconf
         dbsdt(i)%allconf%s_config(j)=(Img/2.0d0)*(zdt*zc-zcdt*z)-HMCE(j,j)
      end do
   end do 
   i=0
   j=0
   !d derivs
   bf1%allconf%d_config(:)=0.0d0
   bf1%allconf%s_config(:)=0.0d0
   bf1cz%allconf%s_config(:)=0.0d0
   bf1cz%allconf%d_config(:)=0.0d0
   HMCE(:,:)=(0.d0,0.0d0)
   do i=1,nbf
      bf1=bs(i)
      call bfconz(bf1,bf1cz)
      call H_MCESpin(bf1cz,bf1,HMCE)
      do j=1,nconf
            IF (j .eq. 1) then
              k=2  
            ELSE IF (j .eq. 2) then 
              k=1 
            ELSE 
              Print *, 'fatal error in derivs, stop'
              STOP
            END IF 
            d=-1.0d0*Img*HMCE(j,k)*bf1%allconf%d_config(k)*&
               &cdexp(Img*(bf1%allconf%s_config(k)-bf1%allconf%s_config(j)))
            dbsdt(i)%allconf%d_config(j)=d
      end do
   end do
   !a derivs
   !D derivs
   do i=1,nbf
      dbsdt(i)%d=0.0d0
      dbsdt(i)%allconf%a(1:nconf)=0.0d0
   end do
   Print *, 'finished derivs, deallocating'
   call OUTBASISFN(bf1)
   call OUTBASISFN(bf1cz)
   call OUTBASISFN(bf2)
   call OUTBASISFN(dbf1)
   call OUTBASISFN(dbf1cz)

 End Subroutine Derivs

!$********************************************************************************************************************!!
!$                                                      Subroutine Ampdmce
!$      This subroutine will do what Amp2d for the ccs code but for the mce code
!$      The main problem with this for the ccs was its length so it will be split up into smaller routines
!$      
!$      Arguments:
!$      bs- input-  bs(i)%d will be used
!$      ddot-output- which will go into dbsdt(i)%d
!$      external:
!$      need to call overlap matrix and  Hord (Delta2Hord) functions
!$      need to call LAPACK
!$      Equation solved in element form:
!$      Sum_i D_i^dot *e^i/hbar*S_i*Omega_ji=i/hbar*Sum_i d_i e^(i/hbar*S_i)*Omega_ji*D2Hordji -old eqn
!$      Sum_i <phi_j|phi_i> D_i^dot= -i*Sum_i D2Hordji* D_i    
!$      In LAPACK I will solve
!$      ==> Omega*x=b
!$      where Omega is the overlap matrix for the phis
!$      x is the vector D_idot
!$      b is the vector:
!$      b=D2Hordji*Di - matrix*vector=vector
!
!$********************************************************************************************************************!!

    Subroutine Ampdmce (bs, ddot,dbsdt)
      Implicit none
      type(csbasisfn),dimension(:),intent(in)::bs 
      type(csbasisfn),dimension(:),intent(in)::dbsdt
      type(csbasisfn)::bf1,bf2
      complex(kind=8),dimension(:),intent(out)::ddot
      complex(kind=8),dimension(nbf)::D
      complex(kind=8)::el
      integer::i, j
      complex(kind=8),dimension(:,:),allocatable::omphi, dhord2
      complex(kind=8),dimension(:),allocatable::rhs
      ! variables for LAPACK
      integer, dimension(nbf) ::IPIV
      integer:: LDA,LDB, N, NRHS, INFO
      complex(kind=8),dimension(:,:),allocatable::omphiin
!      complex(kind=8),dimension(nbf)::bout
      !INTERFACE
       include "./LAPACK/ZGESV.f90"
      !END INTERFACE
      allocate(omphi(nbf,nbf))
      allocate(omphiin(nbf,nbf))
      allocate(dhord2(nbf,nbf))
      allocate(rhs(nbf))
      call ovlphimat2(omphi,bs)
      D(1:nbf)=bs(1:nbf)%d
      !LHS = omphi*ddot so have that sorted now
      !move on to RHS=-i*D2Hordji * Di - call D2hord then do mat mul times D
      ! then put into lapack and have got it solved
      !complex(kind=8),dimension(nbf:nconf)::ai,aj
      !LHS= <phi_j|phi_i> D_i
      call D2HordMCE(dhord2,bs,dbsdt)
      rhs(:)=-1.0d0*Img*MATMUL(dhord2(:,:),D(:))
     ! do i=1,nbf
     !    do j=1,nbf
     !       bf1=bs(i)
     !       bf2=bs(j)
     !       !call ovlphiel(bf1,bf2,el)
     !    end do
     ! end do
!      !Print *, 'b before LAPACK', b
!      
!      ! matrix equation now is : omphi * ddot=rhs ,
!      ! variables for LAPACK****************************************************************
      omphiin(:,:)=omphi(:,:)
!      omegain(:,:)=OMEGA(:,:)
      LDA=nbf
      LDB=nbf
      NRHS=1
      N=nbf
      INFO=0
      IPIV=0
!      !Print *, 'LDA', LDA, 'LDB', LDB,'NRHS', NRHS, 'N', N, 'INFO', INFO, 'IPIV', IPIV
!      !***********************************************************************************
      CALL ZGESV(N, NRHS, omphiin , LDA, IPIV, D , LDB, INFO )
!      bout(:)=b(:)
      If (INFO .ne. 0) then
         Print *, 'WARNING INFO =', INFO
      End If
      ddot(:)=D(:)
!      do j=1,nbf
!        ddot(j)=bout(j)/E(j,j)
!      end do
!      b(:)=MATMUL(Omega,bout)! omegain gets changed by LAPACK so need to use OMEGA
!      !Print *, 'b after LAPACK', b
!      !STOP
      deallocate(omphiin)
      deallocate(omphi)
      deallocate(dhord2)
      deallocate(rhs)
    End Subroutine Ampdmce

    Subroutine D2hordMCE(omega,bs, dbsdt)
      Implicit none
      type(csbasisfn),dimension(:),intent(in)::bs
      type(csbasisfn),dimension(:),intent(in)::dbsdt
      complex(kind=8),dimension(:),allocatable::zdoti
      complex(kind=8),dimension(:,:),intent(out)::omega
      complex(kind=8)::dh2elji
      integer ::i,j,k,n, nbf2
      nbf2=nbf**2 ! nbf is global variable
      If (size(omega) .ne. nbf2 ) then
         Print *, 'error in size of omega, stop!'
         Print *, 'size of omega', size(omega)
         STOP
      End if
      n=size(bs)
      if (n .ne. nbf) then
         Print *, "error in OVRLMAT routine"
      end if
      allocate(zdoti(ndim))
      do i=1,n-1
         zdoti(1:ndim)=dbsdt(i)%z(1:ndim)
         omega(i,i)=cmplx(0.0d0,0.0d0)! ensures diagonal elements are 0           
         do j=i+1,n
            call d2hel(dh2elji,bs(j),bs(i),zdoti)
            omega(i,j)=dh2elji
            !omega(j,i)=dconjg(omega(i,j))
         end do
      end do
      omega(n,n)=cmplx(0.0d0,0.0d0)      
      deallocate(zdoti)
    End Subroutine D2hordMCE

    Subroutine d2hel(dh2elji,bf1,bf2,zdoti)
      Implicit none
      type(csbasisfn),intent(in)::bf1,bf2
      type(csbasisfn)::bf1cz,bf2cz
      complex(kind=8),intent(out)::dh2elji
      complex(kind=8),dimension(:),intent(in)::zdoti
      complex(kind=8)::ovphi,ovz
      complex(kind=8),dimension(nconf,nconf)::HMCE1cz2,HMCE2cz2
      integer::i,j,k
      complex(kind=8),dimension(nconf)::a1,a2
      complex(kind=8)::A,B,C ! split up the D2Hord into 3 terms called A,B,C
      complex(kind=8)::zdst ! term for zdot_i*(zj*-zi*)
      ovphi=ovlphiel(bf1,bf2)
      ovz=OVRLPELMT(bf1,bf2)
      call bfconz(bf1,bf1cz)
      call bfconz(bf2,bf2cz)
      call H_MCEspin(bf1cz,bf2,HMCE1cz2)
      call H_MCEspin(bf2cz,bf2,HMCE2cz2)
      dh2elji=cmplx(0.0d0,0.0d0)
      call Mka(bf1,a1)
      call Mka(bf2,a2)
      zdst=sum((bf1%z(1:ndim)-bf2%z(1:ndim))*zdoti(1:ndim))
      A=cmplx(0.0d0,0.0d0)
      B=cmplx(0.0d0,0.0d0)
      do i=1,nconf
         do j=1,nconf
            A=A+ovz*HMCE1cz2(i,j)*dconjg(a1(i))*a2(j)
            B=B+ovz*HMCE2cz2(i,j)*dconjg(a1(i))*a2(j)
!            do k=1,ndim
!            dh2elji=dh2elji+(ovz*HMCE1cz2(j,i)*dconjg(aj(j))*ai(i))-&
!                   &(ovz*HMCE2cz2(j,i)*dconjg(aj(j))*ai(i))-Img*ovphi!*(bf1cz%z(k)-bf2cz%z(k))*z1dot(k)
!         end do
         end do
      end do
        C=Img*ovphi*zdst
        dh2elji=A-B-C
    End Subroutine d2hel

! Subroutine bfconz takes an allocated basis function bf1 and returns the complex conjugate bf1cz
! bf1cz must be unallocated! - now in bset
!$********************************************************************************************************************!!
!$                                                      Subroutine Amp2d
!$      This subroutine will do what Ampsd was supposed to do but failed
!$      I have decided to keep this routine out of the derivs routine to keep the derivs routine more simple
!$      Arguments:
!$      bs- input-  bs(i)%d will be used
!$      ddot-output- which will go into dbsdt(i)%d
!$      external:
!$      need to call overlap matrix and  Hord (Delta2Hord) functions
!$      need to call LAPACK
!$      Equation solved in element form:
!$      Sum_i D_i^dot *e^i/hbar*S_i*Omega_ji=i/hbar*Sum_i d_i e^(i/hbar*S_i)*Omega_ji*D2Hordji
!$      In LAPACK I will solve
!$      ==> Omega*x=b
!$      where Omega is the overlap matrix
!$      x is the vector ddot*exp(i/hbar*Si)                     - ddot *exp is a scalar multiplication
!$      b is the vector:
!$      b=OmegaD2*diexp(i/hbar*Si)                              -Matrix vector multiplication
!$      where OmegaD2 is a matrix where the elements are constructed from a scalar multiplication of the ij elements 
!$      diexp(i/hbar*si) is a vector where the element di is multiplied by exp(i/hbar*si) - also called vector a!
!$      After solved with LAPACK need to divide x by exp(i/hbar*Si) to get ddot
!$********************************************************************************************************************!!


!    Subroutine Amp2d(bs,ddot)
!      Implicit none
!      type(csbasisfn),dimension(:),intent(in)::bs
!      type(csbasisfn)::bf1,bf2
!      complex(kind=8),dimension(:),intent(out)::ddot
!      complex(kind=8),dimension(:),allocatable::a
!      complex(kind=8),dimension(:),allocatable::b
!      complex(kind=8),dimension(:,:),allocatable::OMEGA
!      complex(kind=8),dimension(:,:),allocatable::D2
!      complex(kind=8),dimension(:,:),allocatable::OMD2
!      complex(kind=8),dimension(nbf,nbf)::E ! new phase matrix
!      integer::i,j
!      ! variables for LAPACK
!      integer, dimension(nbf) ::IPIV
!      integer:: LDA,LDB, N, NRHS, INFO
!      complex(kind=8),dimension(nbf,nbf)::Omegain
!      complex(kind=8),dimension(nbf)::bout
!      !INTERFACE
!       include "./LAPACK/ZGESV.f90"
!      !END INTERFACE
!      allocate(OMEGA(nbf,nbf))
!      allocate(OMD2(nbf,nbf))
!      allocate(D2(nbf,nbf))
!      allocate(a(nbf))
!      allocate(b(nbf))
!      call OVRLMAT(OMEGA,bs)
!      CALL INBASISFN(bf1)
!      CALL INBASISFN(bf2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Right Hand Side!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! sort out right hand side matrix
!      do i=1,nbf
!         bf2=bs(i)
!         do j=1,nbf
!            bf1=bs(j)
!            IF ( i .eq. j) then
!               D2(j,i)=cmplx(0.0d0, 0.0d0)
!            ELSE IF (i .ne. j) then
!                D2(j,i)=Delta2Hord(bf1,bf2)
!            END IF
!         end do
!      end do
!      !Print *, 'D2', D2(1:nbf,1:nbf)
!      do i=1,nbf 
!         do j=1,nbf
!                IF (j .ne. i) then
!                        E(j,i)=0.0d0
!                Else if (j .eq. i) then
!                        E(i,j)=exp(Img/hbar*(bs(i)%s))!-bs(j)%s))  - find it confusing as it matches what i have in notes but differs to paper
!                        !Print *, 'e', E(i,j)
!                End IF
!         end do
!      end do
!      do i=1,nbf
!         do j=1,nbf
!            OMD2(j,i)=OMEGA(j,i)*D2(j,i)
!         end do
!      end do
!      ! sort out right hand side vector a
!      do j=1,nbf
!        a(j)=bs(j)%d*E(j,j)                              
!      end do
!      b(:)=MATMUL(OMD2,a)
!      b(:)=-1.0d0*(Img/hbar)*b(:)
!      !Print *, 'b before LAPACK', b
!      ! matrix equation now is : Omega * x=b , where x is ddot*exp(i/hbar*si) and b is vector above- solve with LAPACK
!      ! variables for LAPACK****************************************************************
!      omegain(:,:)=OMEGA(:,:)
!      LDA=nbf
!      LDB=nbf
!      NRHS=1
!      N=nbf
!      INFO=0
!      IPIV=0
!      !Print *, 'LDA', LDA, 'LDB', LDB,'NRHS', NRHS, 'N', N, 'INFO', INFO, 'IPIV', IPIV
!      !***********************************************************************************
!      CALL ZGESV(N, NRHS, omegain , LDA, IPIV, b , LDB, INFO )
!      bout(:)=b(:)
!      If (INFO .ne. 0) then
!         Print *, 'WARNING INFO =', INFO
!      End If
!      ! ddot=b/exp(i/hbar*si)
!      do j=1,nbf
!        ddot(j)=bout(j)/E(j,j)
!      end do
!      b(:)=MATMUL(Omega,bout)! omegain gets changed by LAPACK so need to use OMEGA
!      !Print *, 'b after LAPACK', b
!      !STOP
!      deallocate(D2)
!      deallocate(OMD2)
!      deallocate(OMEGA)
!      deallocate(b)
!      deallocate(a)
!      deallocate(bf1%z)
!      deallocate(bf2%z)
!    End Subroutine Amp2d


!$********************************************************************************************************************!!
!$                                                      Subroutine Ampdconfig
!$      Arguments:
!$      bs- input-  bs(i)%allconfig%d_config will be used
!$      ddot-output- which will go into dbsdt(i)%allconfig%d_config
!$      external:
!$      #need to call overlap matrix and  Hord (Delta2Hord) functions
!$      #need to call LAPACK
!$      Equation solved in element form:
!$      Sum_i D_i^dot *e^i/hbar*S_i*Omega_ji=i/hbar*Sum_i d_i e^(i/hbar*S_i)*Omega_ji*D2Hordji
!$      D_i^dot =i/hbar*Omega_ji*Hordji*D_j*exp(i*sj-si)
!$********************************************************************************************************************!!
    Subroutine Ampdconfig(bf1,dconfdot)
      Implicit none
      type(csbasisfn),intent(in)::bf1
      type(csbasisfn)::bf1cz
      complex(kind=8),dimension(:),intent(out)::dconfdot
      complex(kind=8),dimension(nconf)::a
      complex(kind=8),dimension(nconf)::ac
      complex(kind=8),dimension(nconf,nconf)::Hmce
      integer::i,j,k
      call bfconz(bf1,bf1cz)
      call Mka(bf1,a)
      call Mka(bf1cz,ac)
      call H_MCESpin(bf1cz,bf1,Hmce)
      do i=1,nconf
         IF (i .eq. 1) then
         j=2
         else if (i .eq. 2) then
         j=1
         else 
         print *, "do not understand j in Ampdconfig"
         STOP
         end if
          !dconfdot(i)=-1.0d0*Img*HMce(i,j)*bf1%allconf%d_config(j)&
          !           &*cdexp(Img*(bf1%allconf%s_config(j)))
          dconfdot(i)=-1.0d0*Img*HMce(i,j)*ac(i)*a(j)*bf1%allconf%d_config(j)&
                     &*cdexp(Img*(bf1%allconf%s_config(j)-bf1%allconf%s_config(i)))
      end do
    End Subroutine Ampdconfig

!$********************************************************************************************************************!!
!$                                                      Subroutine RK4ccs
!$      This subroutine will be heavily based on the one in the numberical
!$      recipes book, p1297 (second volume) but it will work for the type
!$      defined in the basis set module.
!$      Double precision will be used instead of SP and DP to be consistent
!$      No interface will be needed as the derivs routine will also be in this
!$      module. 
!$      The arguments are:
!$      y - this is bs
!$      dydx- this is bsdot (the deriv of bs)
!$      x- time t- I will actually change this to t so not to be confusing
!$      h- time step, I will call this dt 
!$      yout- the integrated output of the bs
!$      The RK4 routine has derivs as an argument but I dont think it is needed
!$      if it is in the same module.
!$********************************************************************************************************************!!

    Subroutine RK4ccs(y,dydx,x,h,yout)
      Implicit none
      type(csbasisfn),dimension(:),intent(in)::y,dydx
      type(csbasisfn),dimension(:),intent(out)::yout
      real(kind=8),intent(in)::x,h
      real(kind=8)::h6,hh,xh
      type(csbasisfn),dimension(:),allocatable::dym,dyt,yt
      !if(.not. (size(y) == size (dydx) .OR. size(y)== size(yout) .or. size(yout)==(size(dydx)) ))then
      !   Print *, 'error in size of bs arrays in RK4ccs'
      !end if
      !Print *, ' starting derivs'
      CALL ALCBASIS(dym)
      CALL ALCBASIS(dyt)
      CALL ALCBASIS(yt)
      !Print *, 'allocated memory in derivs'
      hh=h*0.5d0
      h6=h/6.0d0
      xh=x+hh
      yt(:)=hh*dydx(:)
      yt(:)=yt(:)+y(:)
      !yt=y+tt*dydx ! first step
      !Print *, 'in rk4 calling derivs'
      call derivs(xh,yt,dyt)!second
      ! equation was: yt=y+tt*dyt but had to split it up
      yt(:)=hh*dyt(:)!+y
      yt(:)=yt(:)+y(:)
      !Print *, 'in rk4 calling derivs'
      call derivs(xh,yt,dym)! third step
      !equation was: yt= y+tt*dym
      yt(:)=h*dym(:) ! changed this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      yt(:)=yt(:)+y(:)
      dym(:)=dyt(:)+dym(:)
      !Print *, 'in rk4 calling derivs'
      call derivs(x+h,yt,dyt)! fourth step
      ! equation was: yout=y+t6*(dydx+dyt+2.0d0*dym)! accumulate the proper weights
      dym(:)=2.0d0*dym(:)
      yout(:)=dydx(:)+dyt(:)
      yout=yout(:)+dym(:)
      yout=h6*yout(:)
      yout=yout(:)+y(:)
      DEALLOCATE(dym)
      DEALLOCATE(dyt)
      DEALLOCATE(yt)
      !Print *, 'finished derivs'
    End Subroutine RK4ccs
      

  End Module Integration

