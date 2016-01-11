      program parcheck
!  computes the principal axis of rotation for an ensemble of polymer
!  structures and computes the alignment of the transition dipole
!  moments in relation to this principal component axis

      implicit none
      integer :: i,j,k,l,m,n
      real*8 :: x,y,z
      integer,parameter :: nsnap=500,natom=392,nexst=50
      real*8,dimension(nsnap,natom,3) :: q  ! coordinate matrix
      real*8,dimension(nsnap,3) :: irot  ! principal axis of rotation
      real*8,dimension(nsnap,nexst,5) :: exst ! excited state data
      real*8,dimension(nsnap,nexst) :: align ! alignment of trdip & irot
      character,dimension(natom) :: al  ! atom labels
      character*50 :: c1              ! placeholder for data prefix
      
! read in coordinates
      c1='p3et'
      call readxyz(c1,nsnap,natom,q,al)
! compute principal axis of rotation
      call compiaxes(nsnap,natom,q,al,irot)
! read in excited state data
      call readexst(c1,nsnap,nexst,exst)
! compute alignment of transition dipoles with principal axes
      call alignment(nsnap,nexst,exst,irot,align)

      end 
     
      subroutine alignment(nsnap,nexst,exst,irot,align)
! computes the dot product of the principal axis of rotation and each of
! the transition dipole moments of the excited states
      implicit none
      integer :: i,j,k,l,m,n
      real*8 :: x,y,z
      integer,intent(in) :: nsnap,nexst
      real*8,dimension(nsnap,nexst,5),intent(in) :: exst  ! excited state data
      real*8,dimension(nsnap,3),intent(in) :: irot
      real*8,dimension(nsnap,nexst),intent(inout) :: align 
      real*8,dimension(3) :: v

      do i=1,nsnap
        do j=1,10 !up to nexst
! normalize the transition dipole vectors
          x=0.0
          do k=1,3
            x=x+exst(i,j,k)**2
          end do
          x=1/sqrt(x)
          do k=1,3
            v(k)=exst(i,j,k)*x
          end do
          align(i,j)=dot_product(irot(i,:),v)
          print *, align(i,j)
        end do
      end do
      end
          

      subroutine readexst(c1,nsnap,nexst,exst)
! reads in excited state data
      implicit none
      integer :: i,j,k,l,m,n
      real*8 :: x,y,z
      integer,intent(in) :: nsnap,nexst
      character*50,intent(in) :: c1              ! placeholder for data prefix
      real*8,dimension(nsnap,nexst,5),intent(inout) :: exst  ! excited state data
      character*50 :: fn  ! filename for excited state data
! first read in excited state data
      fn=trim(adjustl(c1))//'data/'//trim(adjustl(c1))//'.trdip.dat'
      open(10,file=trim(adjustl(fn)))
      do i=1,nsnap
        do j=1,nexst
          read(10,*) n,(exst(i,j,k),k=1,5)
! data structure:
! State   X   Y   Z   Dip.Str   Osc.Str
        end do
      end do
      close(10)
      end

      subroutine compiaxes(nsnap,natom,q,al,irot)
! computes the principal axis of rotation for a set of polymer
! coordinates
      implicit none
      integer :: i,j,k,l,m,n
      real*8 :: x,y,z,r,ix,iy,iz
      integer,intent(in) :: nsnap,natom
      real*8,dimension(nsnap,natom,3),intent(inout) :: q  ! coordinate matrix
      character,dimension(natom),intent(in) :: al  ! atom labels
      real*8,dimension(3,3) :: itens ! inertia tensor
      real*8,dimension(nsnap,3),intent(inout) :: irot  ! princ axis of rot
      real*8,dimension(natom) :: mass  ! atomic weight of atoms
      real*8,dimension(9) :: work  ! array for lapack
      real*8, dimension(3) :: s,ieigen ! array for maths

! populate the atomic mass array
      do i=1,natom
        if(al(i).eq.'H') then
          mass(i)=1.0079
        elseif(al(i).eq.'C') then
          mass(i)=12.011
        else if(al(i).eq.'S') then
          mass(i)=32.06
        else
          print *, 'mass of atom '//al(i)//' not found'
        endif
      end do

! center the structure on the center of mass
      do i=1,nsnap
        r=1.0/sum(mass)
        do j=1,3  ! compute center of mass
          s(j)=r*dot_product(mass,q(i,:,j))
        end do
        do j=1,3  ! center the structure
          q(i,:,j)=q(i,:,j)-s(j)
        end do
      end do
      
! compute the intertia tensor         
      q=q**2  ! the r**2 term
      do i=1,nsnap
        do j=1,3
          s(j)=dot_product(mass,q(i,:,j)) ! m*r**2
        end do
        do j=1,3
          do k=1,3 ! fill intertia tensor
            itens(j,k)=s(j)*s(k) 
          end do
        end do
! now diag inertia tensor
        call dsyev('V','U',3,itens,3,ieigen,work,9,m)
        irot(i,:)=itens(:,3)
      end do        

      end

      subroutine readxyz(c1,nsnap,natom,q,al)
! reads in atomic coordinates of polymer structures
      implicit none
      integer :: i,j,k,l,m,n
      real*8 :: x,y,z
      integer,intent(in) :: nsnap,natom
      real*8,dimension(nsnap,natom,3),intent(inout) :: q  ! coordinate matrix
      character*50,intent(in) :: c1              ! placeholder for data prefix
      character,dimension(natom),intent(inout) :: al  ! atom labels
      character*50 :: c5,fn   ! character for sequential data files

      do i=1,nsnap
        write(c5,'(I5)') i
        fn=trim(adjustl(c1))//'data/'//trim(adjustl(c5))//'.xyz'
        open(10,file=trim(adjustl(fn)))
        do j=1,natom
          read(10,*) al(j),(q(i,j,k),k=1,3)
        end do
        close(10)
      end do

      end
