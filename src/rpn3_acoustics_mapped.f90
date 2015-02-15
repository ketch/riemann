! ==================================================================
subroutine rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! ==================================================================

! Riemann solver for the acoustics equations in 3d, with varying
! material properties and mapped grids.
!
! waves: 2
! equations: 4
! aux fields:
!
! Conserved quantities:
!       1 pressure
!       2 x_velocity
!       3 y_velocity
!       4 z_velocity
!
! Auxiliary variables:
!       1 a_x
!       2 a_y
!       3 a_z
!       4 area_ratio_left
!       5 b_x
!       6 b_y
!       7 b_z
!       8 area_ratio_front
!       9  c_x
!       10 c_y
!       11 c_z
!       12 area_ratio_bottom
!       13 volume_ratio
!       14 impedance
!       15 sound_speed

! Note that although there are 4 eigenvectors, two eigenvalues are
! always zero and so we only need to compute 2 waves.

! Solve Riemann problems along one slice of data.
! This data is along a slice in the x-direction if ixyz=1
!                               the y-direction if ixyz=2.
!                               the z-direction if ixyz=3.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr


    implicit none
    double precision, intent(out) :: wave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision, intent(out) ::    s(mwaves,1-mbc:maxm+mbc)
    double precision, intent(in)  ::   ql(meqn,1-mbc:maxm+mbc)
    double precision, intent(in)  ::   qr(meqn,1-mbc:maxm+mbc)
    double precision, intent(out) :: amdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(out) :: apdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(in)  :: auxl(maux,1-mbc:maxm+mbc)
    double precision, intent(in)  :: auxr(maux,1-mbc:maxm+mbc)
    double precision :: delta(3), zi, zim, a1, a2, nx, ny, nz, ci, cim
    double precision :: u_normal_left, u_normal_right
    integer, intent(in) :: ixyz, maxm, meqn, mwaves, mbc, mx, maux
    integer :: i, m, inx, iny, inz, iratio, mw

    if (ixyz == 1) then
        inx = 1
        iny = 2
        inz = 3
        iratio = 4
    else if (ixyz == 2) then
        inx = 5
        iny = 6
        inz = 7
        iratio = 8
    else if (ixyz == 3) then
        inx = 9
        iny = 10
        inz = 11
        iratio = 12
    endif

    do i = 2-mbc, mx+mbc
        nx = auxl(inx,i)
        ny = auxl(iny,i)
        nz = auxl(inz,i)
        u_normal_left  = nx*ql(2  ,i) + ny*ql(3  ,i) + nz*ql(4  ,i)
        u_normal_right = nx*qr(2,i-1) + ny*qr(3,i-1) + nz*qr(4,i-1)

        delta(1) = ql(1,i) - qr(1,i-1)
        delta(2) = u_normal_left - u_normal_right

        zi  = auxl(14,i)
        zim = auxr(14,i-1)
        ci  = auxl(15,i)
        cim = auxl(15,i-1)

        a1 = (-delta(1) + zi*delta(2)) / (zim + zi)
        a2 =  (delta(1) + zim*delta(2)) / (zim + zi)

        wave(1,1,i) = -a1*zim
        wave(2,1,i) = a1 * nx
        wave(3,1,i) = a1 * ny
        wave(4,1,i) = a1 * nz
        s(1,i) = -cim * auxl(iratio,i)
    
        wave(1,2,i) = a2*zi
        wave(2,2,i) = a2 * nx
        wave(3,2,i) = a2 * ny
        wave(4,2,i) = a2 * nz
        s(2,i) = ci * auxl(iratio,i)
    end do


    do m=1,meqn
        do i = 2-mbc, mx+mbc
            amdq(m,i) = s(1,i)*wave(m,1,i)
            apdq(m,i) = s(2,i)*wave(m,2,i)
        end do
    end do


    return
    end subroutine rpn3

