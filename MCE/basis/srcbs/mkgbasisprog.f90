!%**************************************************************
!%                      Program MkGBasis
!%      This program will make a basis set chosen from 
!%      a random Gaussian distribution     
!%**************************************************************

Program MkBasis
use BSETMK
use GenHam
use Check
Implicit none
type(csbasisfn),dimension(:),allocatable::bs                 ! basis set
type(csbasisfn)::bf,bfconjz
character(LEN=100)::LINE1                                    ! character
character(LEN=100)::LINE2
! parameters for reading in z (except for LINE1 above):
integer(kind=8):: myseed
integer::outp
real(kind=8):: sigp,sigq
real(kind=8),dimension(:), allocatable::mup, muq
real(kind=8),external::ZBQLNOR
real(kind=8)RLnbf,  Modnbf
integer::i,j, nbfhalf
integer::nout
integer::nbfcount, ntimes, dim1, dim2
real(kind=8)::Totpop, Pop1, Pop2
complex(kind=8)::EHord

! extra variables for initialising d
complex(kind=8), dimension(:),allocatable::c
complex(kind=8), dimension(:),allocatable::zp
complex(kind=8), dimension(:,:),allocatable::omega
complex(kind=8)::psi
! variable for norm
real(kind=8)::norm
! variable for time read in
real(kind=8)::t0,tmax,dt
real(kind=8) :: starttime,stoptime
OPEN(unit=130, FILE='General.out')
OPEN(unit=135, FILE='basisset.dat')
Write(130,*), 'Program is starting'
call CPU_TIME(starttime)
!$**************************************************************
!$      Read in parameters nbf and ndim
!$**************************************************************
call READBSPAR(LINE1)
!$**************************************************************
!$      Allocate memory for the basis set
!$**************************************************************
call alcbasis(bs)
allocate(mup(ndim))
allocate(muq(ndim))
allocate(c(nbf))
call ALCPAR ! allocates memory for gamma, m and w

!$**************************************************************
!$      Read in the parameters that are arrays(ndim)
!$**************************************************************
call READARR(LINE1)
call Readz(myseed,outp, mup, muq,sigp, sigq, LINE1)

!$**************************************************************
!$      Do some quick checks:
!$      Check that gamma given is equal to m and w calc
!$      Check that the input.dat file was read in correctly
!$**************************************************************
call BSGAMCHK
call ReadChk(Line1)

!$**************************************************************
!$      Need to read in the Hamiltonian and its parameters
!$      Otherwise can't initialise the bs with E/pop cut off 
!$      Also for MCE- need to read PES and nconf
!$**************************************************************
call ALCHamPar ! allocates memory to the hamiltonian par
call ReadSysName(LINE1)
Write(130, *), 'Namesys', NameSys
call ReadHamSpec(LINE1, NameSys)
If (NameSys=='Harmonic') then
   Write(130,*), 'checking harmonic parameters'
   call HarmChk
End If

!$**************************************************************
!$      Read in time variables so you know if you are 
!$      initialising at t0
!$**************************************************************

call ReadTime(Line1,t0,tmax,dt)
!$**************************************************************
!$      Need to read in values for the energy and pop cut off
!$**************************************************************
call ReadEcut(Line1)
Write(130,*), 'Emax', Ecut, 'Emin',Emin,'Popmax',  Popmax, 'Popmin',Popmin
!$**************************************************************
!$      Start to initialise
!$      Initialise s, z and then d
!$**************************************************************
IF (t0 .ne. 0.0d0) then
   Write(130,*),  't0=', t0, 'this is not 0.0d0'
   Write(130,*), 'I havent programmed how to initialise at not t0=0'
   Write(130,*), 'Stopping'
   STOP
ELSE
   Write(130,*), 'Initialising S'        
   call IntS(bs)
END IF
!!! IntZ routine!!
Write(130,*), 'Initialising z'
RLnbf=REAL(nbf)
IF (.not. allocated(mup)) then
   Write(130, *), 'WARNING not allocated mup in intz'
END IF
IF (.not. allocated(muq)) then
   Write(130,*), 'WARNING not allocated muq in intz'
END IF
CALL INBASISFN(bfconjz)
CALL INBASISFN(bf)
CALL ZBQLINI(myseed,outp)! this will automatically write the SEED to an outp
nbfhalf=0.5d0*nbf
dim1=1
dim2=2
Write(130,*), 'Pcheck ', PCheck
! symm check
IF (symop =="YES") then
   IF (nbf ==1) then
      Write(130,*), 'This isnt a warning, it is more of a reminder'
      Write(130,*), 'You have asked for a symmetrical one basis fnc basis set'
      Write(130,*), 'Either use more basis function or put sym NO'
   END IF
   Modnbf=mod(RLnbf,2.0d0)
   IF (Modnbf > 0.0d0) then
      Print *, 'WARNING How can I ensure a symmetrical odd numbered basis set?'
      STOP
   END IF
   do i=1,nbfhalf
      bf=bs(i)
      do j=1,ndim
         bf%z(j)=cmplx(SQRT(gam(j)/2.0d0)*(ZBQLNOR(muq(j),sigq))&            
              ,(1.0d0/hbar*SQRT(1.0d0/(2.0d0*gam(j)))*ZBQLNOR(mup(j),sigp)))
      end do
      bs(i)=bf
      bfconjz%s=bf%s
      bfconjz%d=bf%d
      bfconjz%z(1:ndim)=dconjg(bf%z(1:ndim))
      bs(i+nbfhalf)=bfconjz
   end do

ELSE IF (symop =="NO") then

   nbfcount=1
   i=1
   do while (nbfcount .le. nbf)
      nout=0
      bf=bs(i)
      ntimes=0
      do while(nout .eq. 0)
         ntimes=ntimes+1
         IF (ntimes .gt. ntries) then
            Write(130,*), 'ntimes greater than ntries'
            Write(130,*), 'ntimes', ntimes, 'ntries', ntries
            Write(130,*), 'nbfcoutnt', nbfcount, 'nbf', nbf
            Write(130,*), 'Giving up here'
            STOP
            nout=1
         END IF
         do j=1,ndim
            bf%z(j)=cmplx(SQRT(gam(j)/2.0d0)*(ZBQLNOR(muq(j),sigq))&
                 ,(1.0d0/hbar*SQRT(1.0d0/(2.0d0*gam(j)))*ZBQLNOR(mup(j),sigp)))
            bfconjz%z(j)=dconjg(bf%z(j))
         end do
         EHord=H_ord(bfconjz,bf)
         IF (.not. (Ecut .eq. 0.0d0 .AND. Emin .eq. 0.0d0)) then 
            IF ((Ecut .le. Real(EHord)) .OR. (Emin .gt. Real(EHord))) then
               Write(130,*), 'Rejecting bf because of energy, ntime=', ntimes, 'for nbf=', nbfcount, 'EHord', EHord
               Write(130,*), 'Emax=', Ecut, 'Emin=', Emin, 'ntries=', ntries
               nout=0
            ELSE IF (PCheck =="NO") then
               Write(130,*), 'Accepted bf'
               Write(130,*), 'Using bf', nbfcount, 'Ehord', Real(Ehord)
               nout=1
            ELSE IF (.not.(PCheck=="NO")) then
               IF (.not.(NameSys =="Hubbard")) then
                  Write(130,*), 'You want a pop test for a non Hubbard system?'
                  Write(130,*), 'I think this is wrong...stopping'
                  STOP
               END IF
               Pop1=Npopbf(bf,dim1)
               Pop2=Npopbf(bf,dim2)
               Totpop=Pop1+Pop2
               IF ((Totpop .gt. Popmax) .OR. (Totpop .lt. Popmin)) then
                  Write(130,*), 'Rejecting Bf because of pop, ntime=', ntimes, 'for nbf=', nbfcount, 'Totpop', Totpop
                  Write(130,*), 'Popmax=', Popmax, 'Popmin=', Popmin, 'ntries=', ntries
                  nout=0
               ELSE
                  IF ((Ecut .le. Real(EHord)) .OR. (Emin .gt. Real(EHord))) then
                     Write(130,*), 'Rejecting bf because of energy,ntime=', ntimes, 'for nbf=', nbfcount,'EHord', Real(EHord)
                     Write(130,*), 'Emax=', Ecut,'Emin=', Emin, 'ntries=',ntries
                     nout=0
                  ELSE                    
                     Write(130,*), 'Accepted bf for Pop and energy'
                     Write(130,*),  'Using bf', nbfcount, 'Total Pop', Totpop, 'pop1', Pop1, 'pop2', Pop2, 'EHord', Real(EHord)
                     nout=1
                  END IF
               END IF
            END IF
            IF (.not.((nout .eq. 0 ).OR. (nout .eq. 1))) then
               Write(130,*), 'Error with nout, nout=', nout
               STOP
            END IF
         ELSE
            nout=1
         End IF
         
         bs(i)=bf
      END do
      nbfcount=nbfcount+1
      i=i+1
   end do
ELSE 
Print *, 'Havent understood symop, should be YES or NO ', symop
END IF
!ELSE IF (symop =="YES") then
!   Print *, 'Symmetry in the basis set has been specified'
!   Modnbf=mod(RLnbf,2.0d0)
!   IF (Modnbf > 0.0d0) then
!      Print *, 'WARNING How can I ensure a symmetrical odd numbered basis set?'
!      STOP
!   END IF
!   do i=1,nbfhalf
!      bf=bs(i)
!      do j=1,ndim
!         bf%z(j)=cmplx(SQRT(gam(j)/2.0d0)*(ZBQLNOR(muq(j),sigq))&
!              ,(1.0d0/hbar*SQRT(1.0d0/(2.0d0*gam(j)))*ZBQLNOR(mup(j),sigp)))
!      end do
!      bs(i)=bf
!      bfconjz%s=bf%s
!      bfconjz%d=bf%d
!      bfconjz%z(1:ndim)=dconjg(bf%z(1:ndim))
!      bs(i+nbfhalf)=bfconjz
!   end do
!ELSE
!   Print *, 'Havent understood symop'
!   Print *, 'Should be YES or NO in input.dat'
!   Print *, 'Symop=', Symop
!END IF
deallocate (bf%z)
deallocate (bfconjz%z)
!!!! Now initialise D - need to change this to c
Write(130, *), 'Initialising c'
call IntC(c, mup, muq, bs)
! write out basisset
Write(135,*), 'ndof  ', ndim
Write(135,*), 'nconf ', nconf
write(135,*), 'nbasisfns  ', nbf
Write (135,*), 'PES   ', pes
do i=1,nbf
   write(135,*),'basis ', i
   write(135,*), 'C', Real(c(i)), AIMAG(c(i))
   do j=1,ndim
      write(135,*) REAL(bs(i)%z(j)), AIMAG(bs(i)%z(j))
   end do  
end do
close(135)

!!!finish
call CPU_TIME(stoptime)
Write (130,*), 'Finished basis set program'
IF ((stoptime-starttime)/3600.0d0 .gt. 1.0d0)then 
Write (130,*), 'Time taken', (stoptime-starttime)/3600.0d0, 'hours' 
ELSE IF ((stoptime-starttime)/60.0d0 .gt. 1.0d0)then
Write (130,*), 'Time taken', (stoptime-starttime)/60 , 'mins'
ELSE
Write(130,*), 'Time taken', stoptime-starttime, 's'        
END IF
close(130)
END Program MkBasis
