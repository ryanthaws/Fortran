      program metrics 
! a series of structure evals for a triblock co-polymer
      implicit none
      integer i,j,k,l,m,n
      real x,y,z,q
      integer ntri,nsnap,natom
      parameter(ntri=1,nsnap=5000,natom=1560)
      character*2 ci,cj
      character*5 prefix
      dimension q(nsnap,ntri*natom,3)

      do i=40,40,2
        write(ci,'(I2)') i
        do j=40,40,2
          write(cj,'(I2)') j
          prefix=ci//'-'//cj
          print *, prefix
          call readcoord(prefix,q,ntri,nsnap,natom)
          call minpair(prefix,q,ntri,nsnap,natom)
          call dihedral(prefix,q,ntri,nsnap,natom)
!          call contact pair angle
        end do
      end do
      
      end
      
      subroutine dihedral(prefix,q,ntri,nsnap,natom) 
! calculate dihedral and skew properties
      implicit none
      integer i,j,k,l,m,n,ii,jj,kk,ll,mm,nn
      real x,y,z,q,cen,r,s,nv,tv
      integer ntri,nsnap,natom
      dimension q(nsnap,ntri*natom,3),nv(nsnap,ntri*60,3),tv(3,3)
      character*5 prefix

      open(10,file='mpd.'//adjustl(trim(prefix))//'.dat')
      open(21,file='dihed.'//adjustl(trim(prefix))//'.dat')
      open(21,file='skew.'//adjustl(trim(prefix))//'.dat')
      do i=1,nsnap
        do j=1,ntri ! calc normal vectors
          do k=1,30
            n=60*(j-1)+k
            l=11*(k-1)+1+natom*(j-1)
            m=11*(k-1)+1231+natom*(j-1)
            do ii=1,3 ! 1st chain
              x=(q(i,l,ii)+q(i,l+1,ii)+q(i,l+2,ii))/3.
              tv(1,ii)=q(i,l+4,ii)-x
              tv(2,ii)=q(i,l+5,ii)-x
            end do
            tv(3,1)=tv(1,2)*tv(2,3)-tv(1,3)*tv(2,2)
            tv(3,2)=tv(1,3)*tv(2,1)-tv(1,1)*tv(2,3)
            tv(3,3)=tv(1,1)*tv(2,2)-tv(1,2)*tv(2,1)
            r=sqrt(tv(3,1)**2+tv(3,2)**2+tv(3,3)**2)
            tv=tv/r
            nv(i,n,1)=tv(3,1)
            nv(i,n,2)=tv(3,2)
            nv(i,n,3)=tv(3,3)
            do ii=1,3 ! 2nd chain
              x=(q(i,m,ii)+q(i,m+1,ii)+q(i,m+2,ii))/3.
              tv(1,ii)=q(i,m+4,ii)-x
              tv(2,ii)=q(i,m+5,ii)-x
            end do
            tv(3,1)=tv(1,2)*tv(2,3)-tv(1,3)*tv(2,2)
            tv(3,2)=tv(1,3)*tv(2,1)-tv(1,1)*tv(2,3)
            tv(3,3)=tv(1,1)*tv(2,2)-tv(1,2)*tv(2,1)
            r=sqrt(tv(3,1)**2+tv(3,2)**2+tv(3,3)**2)
            tv=tv/r
            nv(i,n+30,1)=tv(3,1)

            nv(i,n+30,2)=tv(3,2)
            nv(i,n+30,3)=tv(3,3)
          end do
        end do
        do j=1,ntri*60 ! backbone dihedrals
          
      end do
      close(21)
      close(10)
      end


      subroutine minpair(prefix,q,ntri,nsnap,natom) 
! calculate minimum pair distance for thiophenes on unique chains
      implicit none
      integer i,j,k,l,m,n,id1,id2
      real x,y,z,q,cen,r,s
      integer ntri,nsnap,natom
      dimension q(nsnap,ntri*natom,3),cen(nsnap,ntri*60,3)
      character*5 prefix

      open(21,file='mpd.'//adjustl(trim(prefix))//'.dat')
      do i=1,nsnap
        do j=1,ntri
          do k=1,30
            n=60*(j-1)+k
            l=11*(k-1)+1+natom*(j-1)
            m=11*(k-1)+1231+natom*(j-1)
            cen(i,n,1)=q(i,l,1)+q(i,l+1,1)+q(i,l+2,1)+q(i,l+3,1)+q(i,l+4,1)
            cen(i,n,2)=q(i,l,2)+q(i,l+1,2)+q(i,l+2,2)+q(i,l+3,2)+q(i,l+4,2)
            cen(i,n,3)=q(i,l,3)+q(i,l+1,3)+q(i,l+2,3)+q(i,l+3,3)+q(i,l+4,3)
            cen(i,n+30,1)=q(i,m,1)+q(i,m+1,1)+q(i,m+2,1)+q(i,m+3,1)+q(i,m+4,1)
            cen(i,n+30,2)=q(i,m,2)+q(i,m+1,2)+q(i,m+2,2)+q(i,m+3,2)+q(i,m+4,2)
            cen(i,n+30,3)=q(i,m,3)+q(i,m+1,3)+q(i,m+2,3)+q(i,m+3,3)+q(i,m+4,3)
          end do
        end do
        cen=cen*0.2
        do k=1,2*ntri
          do l=30*(k-1)+1,30*k
            r=100.
            id1=9999
            do m=1,ntri*60
              if(m.lt.30*(k-1)+1.or.m.gt.30*k)then
                x=cen(i,l,1)-cen(i,m,1)
                y=cen(i,l,2)-cen(i,m,2)
                z=cen(i,l,3)-cen(i,m,3)
                s=sqrt(x**2+y**2+z**2)
                if(s.lt.r) then
                  r=s
                  id1=m
                endif
              endif
            end do
            write(21,*) l,id1,r,i
          end do
        end do

      end do
      close(21)
      end

      subroutine readcoord(prefix,q,ntri,nsnap,natom) 
! read coordinates of trajectories
      implicit none
      integer i,j,k,l,m,n,id1,id2
      real x,y,z,q
      integer ntri,nsnap,natom
      dimension q(nsnap,ntri*natom,3)
      character*5 prefix,c1,c2
  87  format(I5,2A5,I5,3F8.3)

      open(10,file='./trajectories/d'//adjustl(trim(prefix))//'.gro')
      do i=1,nsnap
        read(10,*)
        read(10,*)
        do j=1,natom*ntri
          read(10,87) id1,c1,c2,id2,q(i,j,1),q(i,j,2),q(i,j,3)
!          write(6,87) id1,c1,c2,id2,q(i,j,1),q(i,j,2),q(i,j,3)
        end do
        read(10,*)
      end do
      q=q*10.0
      close(10)
      end
