! program getPartialInbreeding.f90
! ifort -O3 -heap-arrays -fopenmp getPartialInbreeding.f90 -o getPartialInbreeding
! ./getPartialInbreeding pedigree_name first_animal last_animal (default pedigree)
! it reads a pedigree file and writes a file with partial inbreeding coefficients per ancestor
! A Legarra, andres.legarra@uscdcb.com, 6 June 2023


module kinds
  integer, parameter :: single = SELECTED_REAL_KIND( 6, 37 )
  integer, parameter :: double = SELECTED_REAL_KIND( 15, 307 )
  integer, parameter :: extended = SELECTED_REAL_KIND( 18, 4931 )
  integer, parameter :: r4 = SELECTED_REAL_KIND( 6, 37 )
  integer, parameter :: r8 = SELECTED_REAL_KIND( 15, 307 )
  integer, parameter :: r16 = SELECTED_REAL_KIND( 18, 4931 )

  ! current precison for hash storage
  integer, parameter :: rh=r8
end module kinds


module pedigree
use kinds

double precision, allocatable::F(:)
logical:: lmeuw=.true.
contains



subroutine meuw(sss,ddd,fout) 
   !  Meuwissen & Luo
   use iso_fortran_env
   implicit none
   integer :: sss(:),ddd(:)
   integer :: ss(size(sss)),dd(size(sss))
   double precision:: f(0:size(sss)),fout(:)
   integer, allocatable:: point(:)
   double precision, allocatable:: l(:),d(:)
   integer :: is,id,i,j,k,n,ks,kd
   double precision :: r,fi
   real :: t1

   !print*,'Calculating Inbreeding by M&L function'
   n=size(sss)
   ss=sss
   dd=ddd
  t1=seconds()
   allocate(point(n),l(n),d(n))
   point=0;l=0;d=0
   f(0)=-1
   do i=1,n
      if(n>100 .and. mod(i,int(n/10.0))==0) then
          ! if I don't use flush() this print is not synchronized with code !!!
        write(error_unit,*) ' at',i,'from ',n,' animals'
        call flush(error_unit)
      endif
      is=ss(i); id=dd(i)
      ss(i)=max(is,id)
      dd(i)=min(is,id)
      d(i)=.5-.25*(f(is)+f(id))
      if (is==0 .or. id==0) then
         f(i)=0
         cycle
      endif
      fi=-1
      l(i)=1.0
      j=i
      do while (j/=0)
         k=j
         r=0.5*l(k)
         ks=ss(k)
         kd=dd(k)
         if (ks>0) then
            do while(point(k)>ks)
               k=point(k)
            enddo
            l(ks)=l(ks)+r
            if (ks/=point(k))then
               point(ks)=point(k)
               point(k)=ks
            endif
            if (kd>0)then
               do while (point(k)>kd)
                  k=point(k)
               enddo
               l(kd)=l(kd)+r
               if (kd/=point(k)) then
                  point(kd)=point(k)
                  point(k)=kd
               endif
            endif
         endif

         fi=fi+l(j)*l(j)*d(j)
         l(j)=0
         k=j
         j=point(j)
         point(k)=0
      enddo
      f(i)=fi
   enddo
   f(0)=-1
   fout(:)=f(0:n)
   write(error_unit,*)'   Calculating Inbreeding by M&L, elapsed time',seconds()-t1

end subroutine
 real function seconds()
     call cpu_time(seconds)
 end function


function A_times_v(s,d,v,anc,bound) result(x)
!computes x=A*v using Colleau 2002 indirect method 
! modified to compute partial A using ancestor anc (if present)
! and to compute only up to bound (if present) i.e. A(1:bound,1:bound)
implicit none
integer:: i,nanim
integer,optional:: anc,bound
integer:: s(:),d(:)
double precision::v(:),bigD(size(s))
double precision::a(size(s)),b(size(s)),dii
double precision,allocatable::x(:)
integer:: last

allocate(x(size(s)))
if(lmeuw) then
  allocate(F(0:size(s)))
  F=0d0
  call meuw(s,d,F)
  lmeuw=.false.
  ! write inbreeding in file 100 for checking
  do i=0,size(s)
    write(100,*)i, F(i)
  enddo
endif


nanim=size(s)
a=0
if(present(bound))nanim=bound
do i=nanim,1,-1
      a(i)=a(i)+v(i)
      if (s(i)>0) a(s(i))=a(s(i))+0.5*a(i)
      if (d(i)>0) a(d(i))=a(d(i))+0.5*a(i)
enddo


b=0

!if D is precomputed
!b=D*a

do i=1,nanim
      !dii=0.5-0.25*(F(s(i))+F(d(i)))
      ! we can use the simple rule because F(0)=-1
      bigD(i)=0.5-0.25*(F(s(i))+F(d(i)))
      !b(i)=a(i)*dii
enddo

if(present(anc)) then
    if(anc /=0) then
        ! delete all other mendelian sampling variance
        dii=bigD(anc)
        bigD=0d0
        bigD(anc)=dii
    endif
endif
b=bigD*a


x=0
do i=1,nanim
      x(i)=b(i)
      if(s(i)/=0) then
          x(i)=x(i)+0.5*x(s(i))
      endif
      if(d(i)/=0) then
          x(i)=x(i)+0.5*x(d(i))
      endif
enddo
end function

end module


program test
!$ use omp_lib

use iso_fortran_env
use pedigree
implicit none
character(len=200):: pedfile
character(len=100):: dummy(2)
integer:: i,j,k,n,io,cnt,nanc
integer,allocatable:: ped(:,:)
double precision, allocatable:: v(:),x(:),w(:)
double precision:: Fp,val,t1
logical,allocatable:: is_anc(:)
integer,allocatable:: list_anc(:)
integer::first,last

call get_command_argument(1,value=pedfile,status=io)
if(io/=0) then
        write(error_unit,*) ('what pedigree file?')
        read(input_unit,*)pedfile
endif

open(1,file=pedfile,status='old')
n=0
do 
  read(1,*,iostat=io)
  if(io/=0) exit 
  n=n+1
enddo
write(error_unit,*)'pedfile with',n,'lines'
rewind(1)
allocate(ped(n,3))
do i=1,n
  read(1,*) ped(i,:)
enddo

first=1
last=n
call get_command_argument(2,value=dummy(1),status=io)
if(io==0) then
        call get_command_argument(3,value=dummy(2),status=io)
        if(io==0) then
            read(dummy(1),*)first
            read(dummy(2),*)last
        endif
endif       
write(dummy(1),'(i0)')first
write(dummy(2),'(i0)')last
write(*,*)'first= ',first,'last= ',last


!allocate(v(n),x(n),w(n))
allocate(F(0:n))
! compute whole pedigree
F=0d0
call meuw(ped(:,2),ped(:,3),F)
lmeuw=.false.

! two uses of Colleau
! this below is all for illustration !!
!print *,'first task'
! compute the first column of A
!v=0
!v(1)=1
!x=A_times_v(ped(:,2),ped(:,3),v)
!print *,'compute the first column of A'
!do i=1,n
!  write (10,*)i,x(i)
!enddo

!print *,'2nd task'
! sum of all even elements of A
!v=0
!do i=2,n,2
!  v(i)=1
!enddo
!x=A_times_v(ped(:,2),ped(:,3),v)
! here we get v'x=v'Av
!print *,'sum of even elements of A'
!print *,sum(v*x)

open(unit=2,file='PartialInbreeding.txt'//trim(adjustl(dummy(1)))//'to'//trim(adjustl(dummy(2))),status='replace')
!open(unit=2,file='PartialInbreeding.txt',status='replace')
write (2,*)'ancestor individual Fpartial numberOfAncestors'
cnt=0


! Method flagging ancestors

allocate(is_anc(0:n))
t1=seconds()
cnt=0
do i=first,last
    is_anc=.false.
    if (ped(i,2)/=0 .and. ped(i,3)/=0) then
        call flag_ancestors(i)
        nanc=count(is_anc)-1-1
        allocate(list_anc(nanc))
        call list_ancestors(i)
        !$OMP parallel &
        !$OMP          default(none) &
        !$OMP          shared(i,ped,f,n,is_anc,list_anc,nanc) &
        !$OMP          private(k,j,w,v,x,Fp) &
        !$OMP          reduction(+:cnt)
        allocate(v(n),x(n),w(n))
        !$OMP do
        !do j=1,i
        !    if (is_anc(j)) then
        do k=1,nanc
            j=list_anc(k)
                    ! we compute how much inbreeding in i was generated by j using
                    ! Colleau (2002) indirect method
                    v=0; v(ped(i,2))=1
                    w=0; w(ped(i,3))=1
                    ! extract A(ped(i,2),:) for A generated by mendelian sampling variance of j
                    !x=A_times_v(ped(:,2),ped(:,3),v,anc=j)
                    x=A_times_v(ped(:,2),ped(:,3),v,anc=j,bound=maxval([ped(i,2),ped(i,3)]))
                    Fp=x(ped(i,3))/2
                    !Fp=sum(w*x)/2 ![the same bcs w=1]
                    if(Fp>0) then
                        cnt=cnt+1
                        write (2,'(2i9,f12.8,i19)')j,i,Fp,nanc
                    endif
            !endif
        enddo
        !$OMP end do
        deallocate(v,x,w)
        !$OMP end parallel
        deallocate(list_anc)
    endif
    if(mod(i,10000)==0 .and. n>2000) print *,'in : ',i,' number of ancestors: ',count(is_anc),'cumulated time ', &
           seconds()-t1 ,'non zero partialF so far ',cnt
enddo
write(error_unit,*)'found ',cnt,' non zero Partial F'
write(error_unit,*)'time for 2nd approach ',seconds()-t1



contains
    recursive subroutine flag_ancestors(i)
        integer:: i
        if(i==0) then
            is_anc(0)=.true.
            return
        else
            is_anc(i)=.true.
            if (.not. is_anc(ped(i,2))) call flag_ancestors(ped(i,2))
            if (.not. is_anc(ped(i,3))) call flag_ancestors(ped(i,3))
        endif
    end subroutine

    subroutine list_ancestors(i)
        integer:: j,pos,i
        pos=0
        ! skip 0
        do j=1,i
            if (is_anc(j) .and. j/=i)then
                pos=pos+1
                list_anc(pos)=j
            endif
        enddo
    end subroutine

end





