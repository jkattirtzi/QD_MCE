
!% Author: John Kattirtzi, Date: 21/04/09      
!%                                    Module General Hamiltonian
!%
!%      This module is written to read in the type of Hamiltonian and to
!%      assign the correct one to Hord (i.e HordFree, HordHarm, HordMorse)
!%      Subroutines/functions that are needed
!%      read in name of hamiltonian needed
!%      take in name and assign the correct Hord
!%      take in name and assign the correct DHorDz1c
!%      calculate Delta2
!%      calculate total energy that needs to be conserved (classical)
!%      calculate total quantum energy
!%      calculate expectation energy (quantum)
!%
!%      Links with other modules:
!%      Need to link with HamiltonianSpec as that has Hord for each system
!%      Need to link with the integration module as derivs routine calls Hord
!%      Might need to link with BSET to get parameters like ndig
!$
!%********************************************************************************************************************!!

Module GenHam
use BSET2, only:csbasisfn!, ndim, nbf
use BSET2, only:ndim,nbf, DTOC, nconf! does it get ndim from hamspec? it might be worth investigating as it seems to get it anyway
use BSET2, only:INBASISFN, INCONF, OUTBASISFN
use BSET2, only:Mka, bfconz
!use BSET, only: gam, hbar ! needs them for the Echeck mod
use HamiltonianSpec
!$********************************************************************************************************************!!
!$      Global Variables
!$ Name is a global variable that will read in the name of the system (Free_Particle, Harmonic, Morse)
!$********************************************************************************************************************!!

Character(LEN=20)::NameSys='Notread'
!real(kind=8)::Ecut=0.0d0
!integer ntries
contains


!$********************************************************************************************************************!!
!$                                                      Subroutine ReadSysName
!$      This subroutine will read in the name of the system
!$      It will work by doing a do loop where it searches for System: Name
!$      where System: is a character called Line1 and Name is the actual name of the system
!$********************************************************************************************************************!!


  Subroutine ReadSysName(LINE1)
    Implicit None
    Character(Len=10)::LINE1
    !System is the character that should read System: whilst Name is the name of the system
    integer::ierr
    ierr=0
    OPEN(UNIT=30, FILE='inham.dat',STATUS='OLD', iostat=ierr)
    read(30,*,iostat=ierr)LINE1 ! unit 30 was chosen randomly
    do while (ierr==0)
       if(LINE1 == "System:") then
          backspace(30)
          read(30,*,iostat=ierr)LINE1,NameSys
       end if
       read(30,*, iostat=ierr)LINE1
    end do
    close(30)
  End Subroutine ReadSysName

!$********************************************************************************************************************!!
!$                                                      Function H_ord
!$      This is a function that will take in the name of the system read in the routine above and select the
!$      appropriate H_ord. The advantage of doing it this way is that you can use H_ord in an equation without having
!$      to call a subroutine. 
!$      Obviously have to call ReadSysName first otherwise you won't have a system
!$      I'm going to use IF statements here rather than CASE statements
!$      I've decided this becausse it will mean I don't have to write case(1) when I want to use it 
!$      and I can have an if statement if nothing is read to tell me that
!$      Arguments:
!$      bf1cz and bf2- these are the same as the arguments in Hamiltonian Specific module
!$********************************************************************************************************************!!

  Function H_ord(bf1cz,bf2)
    Implicit none
    type(csbasisfn),intent(in)::bf1cz,bf2
    !integer,intent(in)::i,j
    complex(kind=8)::H_ord
    if (NameSys == 'Notread')then ! Name was initialised to notread 
       Print *, 'havent read in System Name in, have you called ReadSysName? '
    else if (NameSys=='SpinBose') then
    !   call H_MCE(bf1cz,bf2,HmceMat)
    H_ord=Hehr(bf1cz,bf2) ! this is going to be used for derivs of z- no i,j
    else if (.not. (NameSys=='SpinBose' )) then
       Print *, 'WARNING, I havent understood name of system!!!!!!!!!!!,'
       Print *, 'I only understand SpinBose'
    end if
  End Function H_ord
  

!  Subroutine H_MCE(bf1cz,bf2,HmceMat)
!    Implicit None
!    type(csbasisfn),intent(in)::bf1cz,bf2
!    complex(kind=8),dimension(nconf:nconf)::HmceMat
!    integer::i
!    do i=1,nconf
!       do j=1,nconf
!          HmceMat(i,j)=H_ord(bf1cz,bf2,i,j)
!       end do
!    end do
!  End Subroutine H_MCE

!H_MCESpin is a subroutine that will get the hamiltonian matrix for the 2d spin boson model


     Subroutine H_MCESpin(bf1cz,bf2,HMceMat)
        Implicit None
        type(csbasisfn),intent(in)::bf1cz,bf2
        complex(kind=8),dimension(nconf,nconf)::HmceMat
        integer::i,j
        IF (nconf .ne. 2) then
           Print *, 'Error in H_MCESpin, I only know how to do a 2 spin boson model'
           STOP
        END IF
        HmceMat(1,2)=cmplx(spdelta,0.0d0)
        HmceMat(2,1)=cmplx(spdelta,0.0d0)
        HmceMat(1,1)=HordHarmB(bf1cz,bf2)+HordHarmC(bf1cz,bf2)+speps
        HmceMat(2,2)=HordHarmB(bf1cz,bf2)-HordHarmC(bf1cz,bf2)-speps
        ! do j=1,nconf
        ! IF (AIMAG(HMCEMat(j,j)) .lt. 0.0000000001) then
        !     print *, 'checking this hmce spin'
        !     print *, 'HMceMat(j,j)', HMCEMAT(j,j)
        !     HmceMat(j,j)=cmplx(DBLE(HMCEMat(j,j)),0.0d0)
        !     print *, 'HMceMat(j,j)', HMCEMAT(j,j)
        !     STOP
        ! END IF
        ! end do
     !   Print *, 'HMCE Spin, 11, 12, 21, 22,', HmceMat(1,1), HmceMat(1,2), HmceMat(2,1), HmceMat(2,2)
      END SUBROUTINE H_MCESpin

  Subroutine HehrA(bs,ener)
  Implicit None
  type(csbasisfn),dimension(:),intent(in)::bs
  complex(kind=8),intent(out)::ener
  type(csbasisfn)::bf1,bf1c
  complex(kind=8),dimension(nconf)::a,ac
  complex(kind=8),dimension(nconf,nconf)::HMCE
  complex(kind=8)::H11, H12, H21, H22
  integer::i,j,k
  ener=(0.0d0,0.0d0)
  CALL INBASISFN(bf1) 
  CALL INBASISFN(bf1c) 
  CALL INCONF(bf1) 
  CALL INCONF(bf1c) 
  do i=1,nbf
     bf1=bs(i)
     call bfconz(bf1,bf1c)
     call Mka(bf1,a)
     ac(1:nconf)=dconjg(a(1:nconf))
     !Print *, 'a', a
     !Print *, 'dconjga', dconjg(a)
     !Print *, 'ac', ac
     call H_MCESpin(bf1c,bf1,HMCE)
     H11=HMCE(1,1)*a(1)*ac(1)
     H22=HMCE(2,2)*a(2)*ac(2)
     H21=HMCE(2,1)*a(2)*ac(1)
     H12=HMCE(1,2)*a(1)*ac(2)
     ener=ener+H11 + H22 + H12 + H21 
  end do
     !Print *, 'H11', H11, 'H22', H22, 'H21', H21, 'H12', H12
     !STOP
  CALL OUTBASISFN(bf1) 
  CALL OUTBASISFN(bf1c) 
  End subroutine


 Subroutine Hehr_tot(bs, ehr)
    ! this function is like the classical one (one do loop) but with amplitudes
    Implicit none
    type(csbasisfn),dimension(:),allocatable,intent(in)::bs
    type(csbasisfn)::bf1cz,bf2
    complex(kind=8), dimension(nconf)::a,ac, act
    integer::i,j
  !  complex(kind=8),dimension(:),allocatable::c
    complex(kind=8)::ehr, HMCE, H1, H2, HMCE2
   ! allocate(c(nbf))
    call INBASISFN(bf1cz)
    call INCONF(bf1cz)
    call INBASISFN(bf2)
    call INCONF(bf2)
    bf2=bs(1)
    call bfconz(bf2,bf1cz)
    call Mka(bf2,a)
    call Mka(bf1cz,ac)
    !call DTOC(c,bs1)
    ehr=cmplx(0.0d0,0.0d0)
    !do i=1,nbf
    !   bf2=bs(i)
    !   call bfconz(bf2,bf1cz)
    !   ehr= ehr + Hehr(bf1cz,bf2)
    !end do
    H1=HordHarmB(bf1cz,bf2)
    H2=HordHarmC(bf1cz,bf2)
    HMCE=HordHarmB(bf1cz,bf2)+HordHarmC(bf1cz,bf2)+speps
    HMCE2=H1+H2+speps
    act=ac(1)*a(1)
    print *, 'H1', H1
    print *, 'H2', H1
    print *, 'HMCE', HMCE
    print *, 'HMCE2', HMCE2
    !print *, 'speps', speps
    print *, 'act', act 
    ehr=HMCE
    print *, 'in hehr-tot', ehr
   ! STOP
    deallocate(bf1cz%z)
    deallocate(bf2%z)
    deallocate(bf1cz%allconf%d_config)
    deallocate(bf2%allconf%d_config)
  End Subroutine Hehr_tot

! Function thats calculates the Ehr Ham <H>=<phi|H|phi>
      Function Hehr(bf1cz,bf)
        Implicit None
        type(csbasisfn),intent(in)::bf1cz,bf
        !integer, intent(in)::i,j
        complex(kind=8)::Hehr, Img
        complex(kind=8),dimension(nconf)::a,ac
        complex(kind=8),dimension(nconf,nconf)::HMceMat
        integer::i,j,k
        Img=(0.0d0,1.0d0)
        !do k=1,nconf
        !   a(k)=bf1cz%allconf%d_config(k)*exp(Img*bf1cz%allconf%s_config(k))
        !   print *, '1s(k)',bf1cz%allconf%s_config(k)
        !end do
      !  print *, 'need to make a, here are bits, d=', bf%d
        call Mka(bf,a)
        call Mka(bf1cz,ac)
        print *, 'a', a
        call H_MCESpin(bf1cz,bf,HMceMat)
        Hehr=cmplx(0.0d0,0.0d0)
        do i=1,nconf
           do j=1,nconf 
             ! print *, 'i=', i, 'j=', j
             ! print *, 'made a- here it is, a(i)', a(i),'aj', a(j)
             ! print *, 'dconjai*aj', dconjg(a(i))*a(j)
              Hehr=Hehr+(HMceMat(i,j)*ac(i)*a(j))!+HMceMat(2,2)+HMceMat(1,2)+HMceMat(2,1)
             ! Print *, 'HMCEMAT_IJ', HMceMat(i,j), 'i=', i, 'j=', j
           end do
        end do
      End Function Hehr

     Subroutine DH_MCESpinDz1c(bf1cz,bf2,DHMceDz1cMat,m)
        Implicit None
        type(csbasisfn),intent(in)::bf1cz,bf2
        integer,intent(in)::m
        complex(kind=8),dimension(nconf,nconf)::DHmceDz1cMat
        integer::i,j
        IF (nconf .ne. 2) then
           Print *, 'Error, I only know how to do a 2 spin boson model'
           STOP
        END IF
        do i=1,nconf
           do j=1,nconf
              IF (i .ne. j) then
                 DHMceDz1cMat(i,j)=cmplx(0.0d0,0.0d0)
          !       Print *, 'i=', i, 'j=', j, 'DHMceDZ1cMAt', DHMceDz1cMat(i,j)
              else if ((i .eq. 1) .AND. (j .eq. 1)) then
                 DHMceDz1cMat(i,j)=DHordDz1cHarmB(bf1cz,bf2,m)+DHarmCDz1c(bf1cz,bf2,m)
         !        Print *, 'i=', i, 'j=', j, 'DHMceDZ1cMAt', DHMceDz1cMat(i,j)
              else if ((i .eq. 2) .AND. (j .eq. 2)) then
                 DHMceDz1cMat(i,j)=DHordDz1cHarmB(bf1cz,bf2,m)-DHarmCDz1c(bf1cz,bf2,m)
        !         Print *, 'i=', i, 'j=', j, 'DHMceDZ1cMAt', DHMceDz1cMat(i,j)
              else  
               Print *, 'Error in DH_MCESpinDz1c, Somehow I didnt get placed'
               Print *, 'i=',i, 'j=', j, 'nconf=',nconf
               Print *, 'Is nconf > 2? It shouldnt be!'
               STOP
              end if
           end do
        end do
        !Print *, 'DHMCEZz1cMat', DHMceDz1cMat
!        !STOP
      END SUBROUTINE DH_MCESpinDz1c

!      Function DHehr2(bf1cz,bf)
!        Implicit None
!        type(csbasisfn),intent(in)::bf1cz,bf
!!        integer,intent(in)::m
!        complex(kind=8)::DHehr2, Img
!        real(kind=8)::test1,test2
!        complex(kind=8),dimension(nconf,nconf)::DHmceMat
!        complex(kind=8),dimension(:),allocatable::a1cz,a2
!        integer::i,j
!        allocate (a1cz(nconf))
!        allocate (a2(nconf))
!        call Mka(bf,a2)
!        call Mka(bf1cz,a1cz)
!        call DH_MCESpinDz1c(bf1cz,bf,DHMceMat)
!      !print *, 'a2', a2(1)
!      !print *, 'DHMat', DHMceMat
!      !test1= AIMAG(a2(1))*AIMAG(DHMCeMat(1,1))
!      !test2= REAL(a2(1))*REAL(DHMCeMat(1,1))
!      !print *, 'AIMAG(A2)*DH', test1
!      !print *, 'REAL(A2)*DH',  test2
!      !print *, 'AIMAGs', AIMAG(a1cz(2)), AIMAG(a2(2))
!      !print *, 'AIMAGs1', AIMAG(a1cz(1)), AIMAG(a2(1))
!      !  DHehr2=cmplx(0.0d0,0.0d0)
!        DHehr2=(DHMceMat(1,1))+(DHMceMat(2,2))
!        !DHehr2=(DHMceMat(1,1)*a1cz(1)*a2(1))+(DHMceMat(2,2)*a1cz(2)*a2(2))
!        !print *, 'bf%s', bf%allconf%s_config
!        !print *, 'bf%d', bf%allconf%d_config
!        !print *, 'a1cz', a1cz
!        !print *, 'a2', a2
!        !STOP
!      !     do j=1,nconf
!      !        do i=1,nconf
!      !        IF (i .eq. 1 .AND. j .eq. 1) then
!      !           Dhehr2=DHMceMat(1,1)*a1cz(1)*a2(1)
!      !        DHehr2=DHehr2+DHMceMat(j,i)*a1cz(i)*a2(j)
!      !     end do
!      !  end do
!        print *, 'DHMCEMat(1,1)', DHMceMat(1,1)
!        print *, 'a1cz', a1cz(1)
!        print *, 'a2', a2(1)
!        print *, 'DHehr2', DHehr2 
!        !STOP
!        deallocate (a1cz)
!        deallocate (a2)
!      End Function DHehr2
      
      Function DHehrNew(bf1cz,bf,j)
        Implicit None
        type(csbasisfn),intent(in)::bf1cz,bf
        integer,intent(in)::j
        complex(kind=8)::DHehrNew
        complex(kind=8),dimension(nconf,nconf)::DHmceMat
        complex(kind=8),dimension(:),allocatable::acz,a
        print *, 'in dhehr new', j
        allocate (acz(nconf))
        allocate (a(nconf))
        call Mka(bf,a)
        acz(1:nconf)=dconjg(a(1:nconf))
        !print *, 'In new',  'a', a, 'ac', acz
        !print *, 'In Dhehrnew, made ac', acz
        !print *, 'bf%co%d', bf1cz%allconf%d_config
        !print *, 'bf%co%d', bf1cz%allconf%s_config
        print *, 'bf1cz', bf1cz%z
        print *, 'bf1z', bf%z
        !STOP
        DHMceMat(1,1)=DHordDz1cHarmB(bf1cz,bf,j)+DHarmCDz1c(bf1cz,bf,j)
        DHMceMat(1,2)=cmplx(0.0d0,0.0d0)
        DHMceMat(2,1)=cmplx(0.0d0,0.0d0)
        DHMceMat(2,2)=DHordDz1cHarmB(bf1cz,bf,j)-DHarmCDz1c(bf1cz,bf,j)
        DHehrNew=DHMceMat(1,1)*a(1)*acz(1)+DHMceMat(2,2)*a(2)*acz(2) + DHMceMat(2,1)*acz(2)*a(1)+DHMceMat(1,2)*acz(1)*a(2)
        !DHehrNew=DHMceMat(1,1)+DHMceMat(2,2)
        IF (j .eq. ndim) then
           print *, 'In DHeHrNew'
           print *, 'DHehrNew', DHehrNew, j
           print *, 'DHMCEMat', DHMceMat
        END IF
        deallocate (acz)
        deallocate (a)
      End Function DHehrNew
 



!  Function DH_ehr_Dz1c(bf1cz,bf2,i,j)
!    Implicit None
!    type(csbasisfn),intent(in)::bf1cz,bf2
!    integer::i,j,k,l
!    complex(kind=8)::Img
!    complex(kind=8),dimension(nconf)::a
!    complex(kind=8):: DH_ehr_Dz1c
!    Img=cmplx(0.0d0,1.0d0)
!    DH_ehr_Dz1c=DHord
!  END Function DH_ehr_Dz1c

 ! Function Hehr(bf1,bf2)
 ! Implicit NONE
 ! type(csbasisfn),intent(in)::bf1,bf2
 ! complex(kind=8)::Hehr, Img
 ! complex(kind=8),dimension(:), allocatable:: a1,a2
 ! integer::i,j
 ! Img=cmplx(0.0d0,1.0d0)
 ! allocate(a1(nconf))
 ! allocate(a2(nconf))
 ! a1=bf1%allconf%d_config*exp(Img*bf1%allconf%d_config)
 ! a2=bf2%allconf%d_config*exp(Img*bf2%allconf%d_config)
 ! Hehr= H_ord(bf1,bf1)*dconjg(a1)*a1 + H_ord(bf2,bf2)*dconjg(a2)*a2 + H_ord(bf1,bf2)*dconjg(a1)*a2+ H_ord(bf2,bf1)*dconjg(a2)*a1
 ! deallocate(a1)
 ! deallocate(a2)
  !END FUNCTION 

!$********************************************************************************************************************!!
!$                                                      Function DH_ord_Dz1c
!$      This does the same as the H_ord function but for DHord/Dz1*
!$      Name is global so isn't an argument
!$      I like how each of the functions takes the same two arguments
!$********************************************************************************************************************!!

!  Function DHehr_Dz1c(bf1cz,bf2,i)
!    Implicit none
!    type(csbasisfn),intent(in)::bf1cz,bf2
!    integer,intent(in)::i
!    complex(kind=8)::DH_ord_Dz1c
!    if (NameSys == 'Notread')then ! Name was initialised to notread 
!       Print *, 'havent read in System Name in, have you called ReadSysName? '
!    else if (NameSys =='Free') then 
!       DH_ord_Dz1c=DHordDz1cFree(bf1cz,bf2,i)
!    else if (NameSys=='Harmonic') then
!       DH_ord_Dz1c=DHordDz1cHarm(bf1cz,bf2,i)
!    else if (NameSys=='Morse') then
!       DH_ord_Dz1c= DHordDz1cMorse2(bf1cz,bf2,i) 
!    else if (NameSys=='Hubbard') then
!       DH_ord_Dz1c=DHordDz1cHubb(bf1cz,bf2,i)     
!    else if (.not. (NameSys=='Morse' .OR. NameSys=='Harmonic' .OR. NameSys=='Free')) then
!       Print *, 'WARNING, I havent understood name of system!!!!!!!!!!!'
!    end if
!  End Function DH_ord_Dz1c

!$********************************************************************************************************************!!
!$                                                      Function DHordDz1
!$      This function gives the derivative of z by dz whilst the function above gives it by z*
!$********************************************************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!I don't actually need the function below anymore
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Function DH_ord_Dz1(bf1cz,bf1,i)
!    Implicit none
!    type(csbasisfn),intent(in)::bf1cz,bf1
!    integer,intent(in)::i
!    complex(kind=8)::DH_ord_Dz1
!    if (NameSys == 'Notread')then ! Name was initialised to notread 
!       Print *, 'havent read in System Name in, have you called ReadSysName? '
!    else if (NameSys =='Free') then
!       DH_ord_Dz1=DHordDz1Free(bf1cz,bf1,i)
!    else if (NameSys=='Harmonic') then
!       !Print *, 'havent written this yet'
!       !STOP
!       DH_ord_Dz1=DHordDz1Harm(bf1cz,bf1,i)
!    else if (NameSys=='Morse') then
!       !Print *, 'havent written this yet'
!       !STOP
!    !   DH_ord_Dz1= DHordDz1Morse(bf1cz,bf1,i)
!    else if (.not. (NameSys=='Morse' .OR. NameSys=='Harmonic' .OR. NameSys=='Free')) then
!       Print *, 'WARNING, I havent understood name of system!!!!!!!!!!!'
!    end if
!  End Function DH_ord_Dz1

!$********************************************************************************************************************!!
!$                                                      Function Delta2Hord
!$      This value uses the Hord and DhordDz values calculated above
!$      It will be used when calculating the amplitudes
!$      This function looks a bit more complicated then others because it is
!$      split into two parts. 
!$      The DHordBit does the sum of the DH_ord function over all the degrees of
!$      freedom. This could not be done in the DH_ord function itself as the
!$      DH_ord(zic,zj) must be multiplied by zi-zj then summed over all the
!$      degrees of freedom. For the H_ord function the sum is written in the
!$      function as it is simpler and this doesn't need to be multiplied.  
!$********************************************************************************************************************!!


!  Function Delta2Hord(bfi,bfj)
!    Implicit none
!    type(csbasisfn),intent(in)::bfi,bfj
!    type(csbasisfn)::bfic, bfjc
!    integer::i,j
!    complex(kind=8)::Delta2Hord
!    complex(kind=8)::DHordBit
!    CALL INBASISFN(bfic)
!    CALL INBASISFN(bfjc)
!    bfic%z(1:ndim)=conjg(bfi%z(1:ndim))
!    bfjc%z(1:ndim)=conjg(bfj%z(1:ndim))
!    DHordBit=cmplx(0.0d0,0.0d0)
!    !do i=1,ndim  !                      - might need to change this for 2d prob
!       do j=1,ndim
!          DHordBit=DHordBit+(DH_ord_Dz1c(bfjc,bfj,j)*(bfic%z(j)-bfjc%z(j)))
!          !DHordBit=DHordBit+(DH_ord_Dz1c(bfjc,bfj,j)*(bfjc%z(j)-bfic%z(i)))
!       end do   
!    !end do
!    !Delta2Hord=H_ord(bfic,bfi)-H_ord(bfic,bfj)+DHordBit
!    Delta2Hord=(H_ord(bfic,bfj)-H_ord(bfjc,bfj)-DHordBit)
!    !Print *, 'Delta2Hord', Delta2Hord
!    deallocate(bfic%z)
!    deallocate(bfjc%z)
!    ! The DH_ord_Dz1c should actually be by dz2c so might have to check if this
!    ! makes a difference
!  End Function Delta2Hord


!       Write Energy functions:

!$********************************************************************************************************************!!
!$                                                      Function Energyclass
!$      This function is to calculate the classical energy for a basis set(bs)
!$      It will take the whole basis set as the argument
!$      and sum over all the degrees of freedom and basis functions
!$      I think the expression for the classical energy is:
!$      classical energy=Sum_i Hord(zic,zi)
!$      I think it doesn't involve zj and there won't be any overlap
!$********************************************************************************************************************!!

  Function Energyclass(bs)
    Implicit none
    type(csbasisfn),dimension(:), intent(in)::bs
    type(csbasisfn)::bf1cz,bf2
    integer::i
    complex(kind=8)::Energyclass ! makes sense for this to be a real number so check if imag is small
    !real(kind=8)::Energyclass
    Energyclass=cmplx(0.0d0,0.0d0)
    call INBASISFN(bf1cz)
    call INBASISFN(bf2)
    do i=1,nbf
       bf2=bs(i)
       bf1cz%z(1:ndim)=conjg(bf2%z(1:ndim))
       Energyclass=Energyclass + (H_ord(bf1cz,bf2))
    end do
    deallocate(bf1cz%z)
    deallocate(bf2%z)
  End Function Energyclass

!  Function EnergyclassBF(bf)
!    Implicit none
!    type(csbasisfn)::bf, bf1c
!    complex(kind=8)::EnergyclassBF
!    call INBASISFN(bf1c)
!    bf1c%s=bf%s
!    bf1c%d=bf%d
!    bf1c%z(1:ndim)=dconjg(bf%z(1:ndim))
!    EnergyclassBF=H_ord(bf1c,bf)
!    deallocate(bf1c%z) 
!  End Function EnergyclassBF

!$********************************************************************************************************************!!
!$                                                      Function Energyquantum
!$      This function calculates the total quantum energy of the basis set (bs)
!$      It will take the whole basis set as the argument, bs also has the
!$      amplitudes in it. 
!$      The expression I will programme is:
!$      Energy=Sum_i Sum_j( D_iconjg * D_j*Hord(zic,zj)
!$********************************************************************************************************************!!

  Function Energyquantum(bs)
    Implicit none
    type(csbasisfn),dimension(:),allocatable,intent(in)::bs
    type(csbasisfn)::bf1cz,bf2
    integer::i,j
    complex(kind=8)::Energyquantum
    call INBASISFN(bf1cz)
    call INBASISFN(bf2)
    Energyquantum=cmplx(0.0d0,0.0d0)
    do i=1,nbf
       bf1cz%z(1:ndim)=conjg(bs(i)%z(1:ndim))
       do j=1,nbf
          bf2=bs(j)
          Energyquantum=Energyquantum + conjg(bs(i)%d)*bs(j)%d*H_ord(bf1cz,bf2)
       end do
    end do
    deallocate(bf1cz%z)
    deallocate(bf2%z)
  End Function Energyquantum

!$********************************************************************************************************************!!
!$                                                      Function Energyexpct
!$      This function calculates the quantum expectation energy of the basis set
!$      It will take the whole bais set (bs) as argument, like the other energy
!$      functions. 
!$      The expression I will programme is:
!$      Total Energy Expectation: Sum_i C_ic * Di * Hord(zic,zi)
!$      This argument needs C so need to call DTOC routine in basis set
!$      I have changed the DTOC routine so that it calls overlap rather than
!$      needing to call overlap in this function too
!$      This keeps the overlap up to date. 
!$********************************************************************************************************************!!
  
  Function Energyexpct(bs)
    ! this function is like the classical one (one do loop) but with amplitudes
    Implicit none
    type(csbasisfn),dimension(:),allocatable,intent(in)::bs
    type(csbasisfn)::bf1cz,bf2
    integer::i
    complex(kind=8),dimension(:),allocatable::c
    complex(kind=8)::Energyexpct
    allocate(c(nbf))
    call INBASISFN(bf1cz)
    call INBASISFN(bf2)
    call DTOC(c,bs)
    Energyexpct=cmplx(0.0d0,0.0d0)
    do i=1,nbf
       bf2=bs(i)
       bf1cz%z(1:ndim)=conjg(bf2%z(1:ndim))
       Energyexpct= Energyexpct + c(i)*bs(i)%d*H_ord(bf1cz,bf2)
    end do
    deallocate(bf1cz%z)
    deallocate(bf2%z)
    deallocate(c)
  End Function Energyexpct


End Module GenHam
