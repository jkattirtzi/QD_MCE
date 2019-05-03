!%********************************************************************************************************************!!
!%                                                      Program TestEverything
!%      This is a program this is written to test the modules written
!%      It will test:
!%      Initialisation of the basis set
!%      Running the trajectory (Integration)
!%      If the classical energy is being conserved
!%      If the norm is being conserved
!%********************************************************************************************************************!!             


Program TestEverything
use BSET2
use HamiltonianSpec                                          
use GenHam                                                   ! namesys is declared here   
use Integration
use OutMod
!use ACFMOD
!use Check
!use Hubbard
Implicit none
type(csbasisfn),dimension(:),allocatable::bs                 ! basis set
character(LEN=100)::LINE1                                    ! character variable that is read in
character(LEN=100)::LINE2
! parameters for reading in z (except for LINE1 above):
integer(kind=8):: myseed
integer::outp
!real(kind=8):: sigp,sigq
real(kind=8),dimension(:), allocatable::mup, muq
!real(kind=8),external::ZBQLNOR
! extra variables for initialising d
complex(kind=8), dimension(:),allocatable::c
complex(kind=8), dimension(:),allocatable::zp
complex(kind=8), dimension(:,:),allocatable::omega
complex(kind=8)::psi
! variable for norm
real(kind=8)::norm
! variables for hamiltonian tests
complex(kind=8):: classical
complex(kind=8):: quantum
complex(kind=8):: total
! variables for integration test
real(kind=8)::t0
real(kind=8)::dt
real(kind=8)::tmax
real(kind=8)::t
type(csbasisfn),dimension(:),allocatable::bsout
type(csbasisfn),dimension(:),allocatable::dbsdt
integer::i
! variable for output
!real(kind=8)::t2
integer::n
integer::j,ierr
integer::dn
! variable for ACF
complex(kind=8)::ACFval
integer::jcount
real(kind=8) :: starttime,stoptime
Print *, 'Hello, programme is starting'
call CPU_TIME(starttime)
!$********************************************************************************************************************!!
!$                                                      Initialise basis set
!$      Read in the input file
!$      Allocate the memory
!$      Initialise z,s, d -call initialisation routines
!$********************************************************************************************************************!!
Print *, 'Read In Basis Set'
call ReadBSIN(bs,c)
call CTODin(c,bs)
call Norm2MCE(bs,norm)
Print *, 'initial norm2', norm
CALL ReadSysName(LINE1)
CALL  ReadHamSpec(LINE1, NameSys)
Print *, 'Namesys', Namesys
call HehrA(bs, quantum)
Print *, 'initial ehr', quantum

!$********************************************************************************************************************!!
!$                                                      Integration part testing
!$      Read in the time parameters
!$      write a do loop to calculate the trajectory
!$      In the do loop print off:
!$      the classical energy
!$      the norm
!$      maybe also look at the quantum and expectation energy
!$********************************************************************************************************************!!
call ReadTime(Line1,t0,tmax,dt)
Print *, 'number of basis functions', nbf, 'dimensions', ndim
Print *, 'initial time', t0, 'tmax', tmax, 'dt', dt
Print *,'*************************************************************************'
call alcbasis(bsout)
call alcbasis(dbsdt)
t=t0
!            Print *, 'z=', bs(1)%z 
!t2=t
!bsout(1:nbf)%z(1:ndim)=(0.0d0,0.0d0)                                       ! might need it later but now giving me problems
!bsout(1:nbf)%d=0.0d0                                       ! might need it later but now giving me problems
!bsout(1:nbf)%allconf%s_config(1:nconf)=0.0d0                                       ! might need it later but now giving me problems
!bsout(1:nbf)%allconf%d_config(1:nconf)=(0.0d0,0.0d0)                                       ! might need it later but now giving me problems
n=0
call ParName(Line1,bs, t0, tmax, dt,norm,classical,quantum, total ) ! output mod
Dn=0
call ReadDn(Dn)
            Print *, 'z=', bs(1)%z 
!            !STOP
!jcount=1
Print *, 'Ready to start main do loop'
Print *, 'starting main do loop'
     do while (ABS(t)<=ABS(tmax))
            print *, 'before derivs'
        call derivs(t,bs,dbsdt)
        ! do i=1,nbf
          Print *, 'dbsdt=', dbsdt(nbf)%z
        !end do 
       !STOP
      ! Print *, 'bsout=', bsout%z(1)
       ! Print *, 'Calling rk4, t=', t, 'dt=', dt
        call RK4ccs(bs,dbsdt,t,dt,bsout)
       call Norm2MCE(bsout,norm)
       !Print *, 'norm2 after derivs before rk4, dbsdt', norm 
       Print *, 'norm after rk4 calculated with bsout', norm
       call HehrA(bsout, quantum) 
      Print *, 'Hehr_tot after rk4 calculated with bsout', quantum
       bs=bsout
       t=t+dt
!       STOP
     end do
     deallocate(bsout)
     deallocate(dbsdt)
          !IF (size(muq) > 0 .AND. size(mup) > 0) then
     IF (allocated(muq)) then
       deallocate(muq)
     end if
     If (allocated(mup)) then
      deallocate(mup)
     end if
     deallocate(c)
    call CPU_TIME(stoptime)
    Print *,  '************************************************************************'
    Print *,  '************************************************************************'
    Print *, 'Time taken:', (stoptime-starttime)/60.0d0, 'mins'
    Print *, 'Outputs:',Trajname, Consvname, Outpname ,  bfoutname
    Print *, 'program terminated ok'
    Print *,  '************************************************************************'
    Print *,  '************************************************************************'

Print *, 'wohooooooooooooooooooo'
END Program TestEverything
