program monte_carlo_dataset
implicit none
include 'mpif.h'
integer i,j,k,l,step,loop,loop2,loop0,init
integer n1,dummyi,info
integer ransu1,ransu2
integer,parameter :: natom = 32
integer,parameter :: dataset = 2000
integer,parameter :: totalset = 2000*7
integer,parameter :: numset = 1000
integer,parameter :: num_loop = 1000
integer ransu(numset),ransu_tmp(numset)
real lattice,a(3),b(3),c(3),x(4)!x(numset)
real posi(totalset,natom,3),posi_set(numset,natom,3)
real pse_ene,pse_old,pse_t,delta_ene,r_ij,prob,pse_ene_tmp
real temperature,energy,free_energy,t1,t2,t_req
character*10 system,dummyc
character*3  e1,cT,cE,cF
integer petot,my_rank,ierr,nodename_len
character*8 nodename
!integer myid

call MPI_INIT      (ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD,petot,ierr)
call MPI_COMM_RANK (MPI_COMM_WORLD,my_rank,ierr)
call MPI_GET_PROCESSOR_NAME(nodename,nodename_len,ierr)
print*,"Hello world",my_rank,petot,nodename!,nodename_len

if( my_rank .eq. 0 )then

  open(unit=999,file="./output.dat",action="write")
  open(unit=1,file="../XDATCAR_bulk_au_1000k_2000steps",action="read")
  open(unit=2,file="../XDATCAR_bulk_au_2000k_2000steps",action="read")
  open(unit=3,file="../XDATCAR_bulk_au_3000k_2000steps",action="read")
  open(unit=4,file="../XDATCAR_bulk_au_4000k_2000steps",action="read")
  open(unit=5,file="../XDATCAR_bulk_au_5000k_2000steps",action="read")
  open(unit=6,file="../XDATCAR_bulk_au_300k_2000steps",action="read")
  open(unit=7,file="../XDATCAR_bulk_au_500k_2000steps",action="read")

  read(1,*)system
  read(1,*)lattice
  read(1,*)a(1),a(2),a(3)
  read(1,*)b(1),b(2),b(3)
  read(1,*)c(1),c(2),c(3)
  read(1,*)e1
  read(1,*)n1
  do i=1,dataset
      read(1,*)dummyc,dummyc,step
    do j=1,natom
      read(1,*)posi(i,j,1),posi(i,j,2),posi(i,j,3)
    enddo
  enddo

  read(2,*)system
  read(2,*)lattice
  read(2,*)a(1),a(2),a(3)
  read(2,*)b(1),b(2),b(3)
  read(2,*)c(1),c(2),c(3)
  read(2,*)e1
  read(2,*)n1
  do i=dataset+1,2*dataset
      read(2,*)dummyc,dummyc,step
    do j=1,natom
      read(2,*)posi(i,j,1),posi(i,j,2),posi(i,j,3)
    enddo
  enddo

  read(3,*)system
  read(3,*)lattice
  read(3,*)a(1),a(2),a(3)
  read(3,*)b(1),b(2),b(3)
  read(3,*)c(1),c(2),c(3)
  read(3,*)e1
  read(3,*)n1
  do i=2*dataset+1,3*dataset
      read(3,*)dummyc,dummyc,step
    do j=1,natom
      read(3,*)posi(i,j,1),posi(i,j,2),posi(i,j,3)
    enddo
  enddo

  read(4,*)system
  read(4,*)lattice
  read(4,*)a(1),a(2),a(3)
  read(4,*)b(1),b(2),b(3)
  read(4,*)c(1),c(2),c(3)
  read(4,*)e1
  read(4,*)n1
  do i=3*dataset+1,4*dataset
      read(4,*)dummyc,dummyc,step
    do j=1,natom
      read(4,*)posi(i,j,1),posi(i,j,2),posi(i,j,3)
    enddo
  enddo

  read(5,*)system
  read(5,*)lattice
  read(5,*)a(1),a(2),a(3)
  read(5,*)b(1),b(2),b(3)
  read(5,*)c(1),c(2),c(3)
  read(5,*)e1
  read(5,*)n1
  do i=4*dataset+1,5*dataset
      read(5,*)dummyc,dummyc,step
    do j=1,natom
      read(5,*)posi(i,j,1),posi(i,j,2),posi(i,j,3)
    enddo
  enddo

  read(6,*)system
  read(6,*)lattice
  read(6,*)a(1),a(2),a(3)
  read(6,*)b(1),b(2),b(3)
  read(6,*)c(1),c(2),c(3)
  read(6,*)e1
  read(6,*)n1
  do i=5*dataset+1,6*dataset
      read(6,*)dummyc,dummyc,step
    do j=1,natom
      read(6,*)posi(i,j,1),posi(i,j,2),posi(i,j,3)
    enddo
  enddo

  read(7,*)system
  read(7,*)lattice
  read(7,*)a(1),a(2),a(3)
  read(7,*)b(1),b(2),b(3)
  read(7,*)c(1),c(2),c(3)
  read(7,*)e1
  read(7,*)n1
  do i=6*dataset+1,7*dataset
      read(7,*)dummyc,dummyc,step
    do j=1,natom
      read(7,*)posi(i,j,1),posi(i,j,2),posi(i,j,3)
    enddo
  enddo


!Initial choice of "numset"
  do i=1,numset
      ransu(i) = (totalset/dataset)*i
!      ransu(i) = dint( x(i)*totalset ) + 1
    do j=1,natom
      posi_set(i,j,1) = posi(ransu(i),j,1)
      posi_set(i,j,2) = posi(ransu(i),j,2)
      posi_set(i,j,3) = posi(ransu(i),j,3)
    enddo
  enddo

endif! my_rank .eq. 0


call MPI_BCAST(lattice,1,mpi_real,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(a,3,mpi_real,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(b,3,mpi_real,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(c,3,mpi_real,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(posi_set,numset*natom*3,mpi_real,0,MPI_COMM_WORLD,ierr)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)


!Pseudo energy evaluation
pse_ene = 0.0
do i=(numset/petot)*(my_rank)+1,(numset/petot)*(my_rank+1)
!do i=1,numset
  do l=1,numset
    do j=1,natom
      do k=1,natom
        if( i .ne. l )then
          r_ij =   ( lattice*a(1)*( posi_set(i,j,1) - posi_set(l,k,1) ) )**2 &
                 + ( lattice*b(2)*( posi_set(i,j,2) - posi_set(l,k,2) ) )**2 &
                 + ( lattice*c(3)*( posi_set(i,j,3) - posi_set(l,k,3) ) )**2 
          pse_ene = pse_ene + (1.0/r_ij)
        elseif( i .eq. l .and. j .ne. k )then
          r_ij =   ( lattice*a(1)*( posi_set(i,j,1) - posi_set(l,k,1) ) )**2 &
                 + ( lattice*b(2)*( posi_set(i,j,2) - posi_set(l,k,2) ) )**2 &
                 + ( lattice*c(3)*( posi_set(i,j,3) - posi_set(l,k,3) ) )**2 
          pse_ene = pse_ene + (1.0/r_ij)
        endif
      enddo
    enddo
  enddo
enddo
          pse_ene = pse_ene/(natom*numset)

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!if( my_rank .eq. 0 )print*,pse_ene,"my_rank",my_rank
!if( my_rank .eq. 1 )print*,pse_ene,"my_rank",my_rank
!if( my_rank .eq. 2 )print*,pse_ene,"my_rank",my_rank
!if( my_rank .eq. 3 )print*,pse_ene,"my_rank",my_rank
call MPI_BARRIER(MPI_COMM_WORLD,ierr)
call MPI_REDUCE(pse_ene,pse_ene_tmp,1,mpi_real,mpi_sum,0,MPI_COMM_WORLD,ierr)
!if( my_rank .eq. 0 )print*,pse_ene_tmp,"my_rank",my_rank,"sum"

if( my_rank .eq. 0 )then
  pse_old = pse_ene_tmp
endif


!**********************
!Monte-Carlo loop start
!**********************

if( my_rank .eq. 0 )then
   call pre_random
   pse_t = 1000.0
endif

do loop=1,num_loop

  if( my_rank .eq. 0 )then
    pse_t = pse_t*0.990
    print*,"loop:",loop
  endif

  if( my_rank .eq. 0 )call cpu_time(t1)

  do loop2=1,numset

    if( my_rank .eq. 0 )then
!      print*,loop,loop2

      call random_number(x)
      ransu1 = int( x(2)*numset ) + 1
      ransu2 = int( x(3)*totalset ) + 1

      info = 0
      do i=1,numset
        if( ransu(i) .eq. ransu2 )info = 1
      enddo
      if( info .eq. 1 )then
!        print*,"cycle",loop,loop2
        cycle
      endif

    endif! my_rank .eq. 0

!Exchange data
    if( my_rank .eq. 0 )then
      do j=1,natom
        posi_set(ransu1,j,1) = posi(ransu2,j,1)
        posi_set(ransu1,j,2) = posi(ransu2,j,2)
        posi_set(ransu1,j,3) = posi(ransu2,j,3)
      enddo
    endif

    call MPI_BCAST(posi_set,numset*natom*3,mpi_real,0,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    pse_ene = 0.0
    do i=(numset/petot)*(my_rank)+1,(numset/petot)*(my_rank+1)
!    do i=1,numset
      do l=1,numset
        do j=1,natom
          do k=1,natom
            if( i .ne. l )then
              r_ij =   ( lattice*a(1)*( posi_set(i,j,1) - posi_set(l,k,1) ) )**2 &
                     + ( lattice*b(2)*( posi_set(i,j,2) - posi_set(l,k,2) ) )**2 &
                     + ( lattice*c(3)*( posi_set(i,j,3) - posi_set(l,k,3) ) )**2
              pse_ene = pse_ene + (1.0/r_ij)
            elseif( i .eq. l .and. j .ne. k )then
              r_ij =   ( lattice*a(1)*( posi_set(i,j,1) - posi_set(l,k,1) ) )**2 &
                     + ( lattice*b(2)*( posi_set(i,j,2) - posi_set(l,k,2) ) )**2 &
                     + ( lattice*c(3)*( posi_set(i,j,3) - posi_set(l,k,3) ) )**2
              pse_ene = pse_ene + (1.0/r_ij)
            endif
          enddo
        enddo
      enddo
    enddo
              pse_ene = pse_ene/(natom*numset)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(pse_ene,pse_ene_tmp,1,mpi_real,mpi_sum,0,MPI_COMM_WORLD,ierr)

    if( my_rank .eq. 0 )then

      pse_ene = pse_ene_tmp

      if( pse_ene < pse_old )then

          ransu(ransu1) = ransu2
!          print*,"Exchange!"
          pse_old = pse_ene

      elseif( pse_ene > pse_old )then

        delta_ene = pse_ene - pse_old
        prob = exp( -(delta_ene/pse_t) )

        if( x(4) < prob )then
          ransu(ransu1) = ransu2
!          print*,"Exchange prob!"
          pse_old = pse_ene
        elseif( x(4) > prob )then
          do j=1,natom
            posi_set(ransu1,j,1) = posi(ransu(ransu1),j,1)
            posi_set(ransu1,j,2) = posi(ransu(ransu1),j,2)
            posi_set(ransu1,j,3) = posi(ransu(ransu1),j,3)
          enddo
!          print*,"Not exchange!"
        endif

      endif! pse_ene < pse_old

    endif! my_rank .eq. 0

    call MPI_BCAST(posi_set,numset*natom*3,mpi_real,0,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  enddo!loop2

    if( my_rank .eq. 0 )then
      call cpu_time(t2)
      t_req = t2 - t1
      print*,"Times:",t_req,"[s]"
    endif

  if( my_rank .eq. 0 )then
    write(999,*)loop,pse_ene,pse_t
  endif

enddo!loop

if( my_rank .eq. 0 )then
  print*,"MC loop finished"
endif


!XDATCAR_selected
if( my_rank .eq. 0 )then
  call bubblesort(numset,ransu)

  open(unit=101,file="./XDATCAR_bulk_au_1000k_2000steps_s",action="write")
  open(unit=102,file="./XDATCAR_bulk_au_2000k_2000steps_s",action="write")
  open(unit=103,file="./XDATCAR_bulk_au_3000k_2000steps_s",action="write")
  open(unit=104,file="./XDATCAR_bulk_au_4000k_2000steps_s",action="write")
  open(unit=105,file="./XDATCAR_bulk_au_5000k_2000steps_s",action="write")
  open(unit=106,file="./XDATCAR_bulk_au_300k_2000steps_s",action="write")
  open(unit=107,file="./XDATCAR_bulk_au_500k_2000steps_s",action="write")

  write(101,*)system
  write(101,*)lattice
  write(101,*)a(1),a(2),a(3)
  write(101,*)b(1),b(2),b(3)
  write(101,*)c(1),c(2),c(3)
  write(101,*)e1
  write(101,*)n1

  write(102,*)system
  write(102,*)lattice
  write(102,*)a(1),a(2),a(3)
  write(102,*)b(1),b(2),b(3)
  write(102,*)c(1),c(2),c(3)
  write(102,*)e1
  write(102,*)n1

  write(103,*)system
  write(103,*)lattice
  write(103,*)a(1),a(2),a(3)
  write(103,*)b(1),b(2),b(3)
  write(103,*)c(1),c(2),c(3)
  write(103,*)e1
  write(103,*)n1

  write(104,*)system
  write(104,*)lattice
  write(104,*)a(1),a(2),a(3)
  write(104,*)b(1),b(2),b(3)
  write(104,*)c(1),c(2),c(3)
  write(104,*)e1
  write(104,*)n1

  write(105,*)system
  write(105,*)lattice
  write(105,*)a(1),a(2),a(3)
  write(105,*)b(1),b(2),b(3)
  write(105,*)c(1),c(2),c(3)
  write(105,*)e1
  write(105,*)n1

  write(106,*)system
  write(106,*)lattice
  write(106,*)a(1),a(2),a(3)
  write(106,*)b(1),b(2),b(3)
  write(106,*)c(1),c(2),c(3)
  write(106,*)e1
  write(106,*)n1

  write(107,*)system
  write(107,*)lattice
  write(107,*)a(1),a(2),a(3)
  write(107,*)b(1),b(2),b(3)
  write(107,*)c(1),c(2),c(3)
  write(107,*)e1
  write(107,*)n1


  do i=1,numset

    if( ransu(i) <= dataset )then
      write(101,*)"Direct configuration=",ransu(i)
      do j=1,natom
        write(101,*)posi_set(i,j,1),posi_set(i,j,2),posi_set(i,j,3)
      enddo
    elseif( dataset < ransu(i) .and. ransu(i) <= 2*dataset )then
      write(102,*)"Direct configuration=",ransu(i)
      do j=1,natom
        write(102,*)posi_set(i,j,1),posi_set(i,j,2),posi_set(i,j,3)
      enddo
    elseif( 2*dataset < ransu(i) .and. ransu(i) <= 3*dataset )then
      write(103,*)"Direct configuration=",ransu(i)
      do j=1,natom
        write(103,*)posi_set(i,j,1),posi_set(i,j,2),posi_set(i,j,3)
      enddo
    elseif( 3*dataset < ransu(i) .and. ransu(i) <= 4*dataset )then
      write(104,*)"Direct configuration=",ransu(i)
      do j=1,natom
        write(104,*)posi_set(i,j,1),posi_set(i,j,2),posi_set(i,j,3)
      enddo
    elseif( 4*dataset < ransu(i) .and. ransu(i) <= 5*dataset )then
      write(105,*)"Direct configuration=",ransu2
      do j=1,natom
        write(105,*)posi_set(i,j,1),posi_set(i,j,2),posi_set(i,j,3)
      enddo
    elseif( 5*dataset < ransu(i) .and. ransu(i) <= 6*dataset )then
      write(106,*)"Direct configuration=",ransu(i)
      do j=1,natom
        write(106,*)posi_set(i,j,1),posi_set(i,j,2),posi_set(i,j,3)
      enddo
    elseif( 6*dataset < ransu(i) .and. ransu(i) <= 7*dataset )then
      write(107,*)"Direct configuration=",ransu(i)
      do j=1,natom
        write(107,*)posi_set(i,j,1),posi_set(i,j,2),posi_set(i,j,3)
      enddo
    endif

  enddo!i


!Corresponding OSZICAR_selected
  open(unit=301,file="./OSZICAR_bulk_au_1000k_2000steps_s",action="write")
  open(unit=302,file="./OSZICAR_bulk_au_2000k_2000steps_s",action="write")
  open(unit=303,file="./OSZICAR_bulk_au_3000k_2000steps_s",action="write")
  open(unit=304,file="./OSZICAR_bulk_au_4000k_2000steps_s",action="write")
  open(unit=305,file="./OSZICAR_bulk_au_5000k_2000steps_s",action="write")
  open(unit=306,file="./OSZICAR_bulk_au_300k_2000steps_s",action="write")
  open(unit=307,file="./OSZICAR_bulk_au_500k_2000steps_s",action="write")

  do i=1,numset

    if( ransu(i) <= dataset )then
      open(unit=201,file="../oszicar_bulk_au_1000k_2000steps",action="read")
      do j=1,ransu(i)-1
        read(201,*)dummyi
      enddo
        read(201,*)step,cT,temperature,cE,energy,cF,free_energy
        write(301,*)step,cT,temperature,cE,energy,cF,free_energy
      close(201)
    elseif( dataset < ransu(i) .and. ransu(i) <= 2*dataset )then
      open(unit=202,file="../oszicar_bulk_au_2000k_2000steps",action="read")
      ransu2 = ransu(i) - dataset
      do j=1,ransu2-1
        read(202,*)dummyi
      enddo
        read(202,*)step,cT,temperature,cE,energy,cF,free_energy
        write(302,*)step,cT,temperature,cE,energy,cF,free_energy
      close(202)
    elseif( 2*dataset < ransu(i) .and. ransu(i) <= 3*dataset )then
      open(unit=203,file="../oszicar_bulk_au_3000k_2000steps",action="read")
      ransu2 = ransu(i) - dataset*2
      do j=1,ransu2-1
        read(203,*)dummyi
      enddo
        read(203,*)step,cT,temperature,cE,energy,cF,free_energy
        write(303,*)step,cT,temperature,cE,energy,cF,free_energy
      close(203)
    elseif( 3*dataset < ransu(i) .and. ransu(i) <= 4*dataset )then
      open(unit=204,file="../oszicar_bulk_au_4000k_2000steps",action="read")
      ransu2 = ransu(i) - dataset*3
      do j=1,ransu2-1
        read(204,*)dummyi
      enddo
        read(204,*)step,cT,temperature,cE,energy,cF,free_energy
        write(304,*)step,cT,temperature,cE,energy,cF,free_energy
      close(204)
    elseif( 4*dataset < ransu(i) .and. ransu(i) <= 5*dataset )then
      open(unit=205,file="../oszicar_bulk_au_5000k_2000steps",action="read")
      ransu2 = ransu(i) - dataset*4
      do j=1,ransu2-1
        read(205,*)dummyi
      enddo
        read(205,*)step,cT,temperature,cE,energy,cF,free_energy
        write(305,*)step,cT,temperature,cE,energy,cF,free_energy
      close(205)
    elseif( 5*dataset < ransu(i) .and. ransu(i) <= 6*dataset )then
      open(unit=206,file="../oszicar_bulk_au_300k_2000steps",action="read")
      ransu2 = ransu(i) - dataset*5
      do j=1,ransu2-1
        read(206,*)dummyi
      enddo
        read(206,*)step,cT,temperature,cE,energy,cF,free_energy
        write(306,*)step,cT,temperature,cE,energy,cF,free_energy
      close(206)
    elseif( 6*dataset < ransu(i) .and. ransu(i) <= 7*dataset )then
      open(unit=207,file="../oszicar_bulk_au_500k_2000steps",action="read")
      ransu2 = ransu(i) - dataset*6
      do j=1,ransu2-1
        read(207,*)dummyi
      enddo
        read(207,*)step,cT,temperature,cE,energy,cF,free_energy
        write(307,*)step,cT,temperature,cE,energy,cF,free_energy
      close(207)
    endif

  enddo!i

endif! my_rank .eq. 0




call MPI_FINALIZE  (ierr)

end

!********************
subroutine pre_random
  implicit none
  integer::seedsize,c
  integer,allocatable::seed(:)
  call random_seed(size=seedsize)
  allocate(seed(1:seedsize))
  call system_clock(count=c)
  seed=c
  call random_seed(put=seed)
  return
end subroutine pre_random

subroutine bubblesort(N,array)
  implicit none
  integer,intent(in)::N
  integer,intent(inout)::array(N)
!  double precision,intent(inout)::array(N)
  integer::i,j
  real::t
  do i=1,N-1
     do j=i+1,N
        if(array(i) .gt. array(j))then
           t=array(i)
           array(i)=array(j)
           array(j)=t
        end if
     end do
  end do
  return
end subroutine bubblesort
