! =====================================================
subroutine rptt3(ixy,icoor,imp,impt,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,bsasdq,cmbsasdq,cpbsasdq)
! =====================================================
    implicit double precision (a-h,o-z)

!     # Dummy transverse Riemann solver, for use in dimensionally-split algorithm.

    write(*,*) 'Error: Dummy transverse Riemann solver called!'
    return
    end subroutine rptt3
