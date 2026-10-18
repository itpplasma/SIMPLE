program test_fo_symplectic
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use field_spectre, only: spectre_field_t, create_spectre_field
    use magfie_sub, only: set_magfie_spectre_field, spectre_field
    use spectre_fo_gauss2, only: gauss2_prepare, gauss2_geometry, &
        gauss2_potential, gauss2_encode, shipped_step => gauss2_step, &
        gauss2_entry_step
    use reference_fo_gauss2, only: reference_prepare => gauss2_prepare, &
        gauss2_hamiltonian, gauss2_step
    use util, only: c, p_mass, e_charge, ev
    use generated_fo_stepper, only: generated_eval, generated_step
    implicit none
    real(dp)::width=0.06_dp,ro0=1.0e5_dp,h=0.02_dp
    real(dp)::omega(6,6),eye(6,6),efixed=0,lfixed=1,gfixed=1,vscale=1,mfixed=0
    type(spectre_field_t)::field
    character(1024)::filename,arg
    integer::case_id,i,nlong,ierr,diagnostic_step=0,only_marker=-1
    real(dp)::step_factor=0.0_dp,energy_tol=1.0e-12_dp
    logical::failure_mode=.false.,energy_failed=.false.
    nlong=1000
    if(command_argument_count()<3.or.command_argument_count()>7) &
        error stop 'two fixtures, step count, optional velocity scale required'
    call get_command_argument(3,arg)
    read(arg,*) nlong
    if(command_argument_count()>=4) then
        call get_command_argument(4,arg)
        if(trim(arg)=='failure') then
            failure_mode=.true.; vscale=100
        else
            read(arg,*) vscale
        end if
    end if
    if(command_argument_count()>=5) then
        call get_command_argument(5,arg)
        read(arg,*) only_marker
    end if
    if(command_argument_count()>=6) then
        call get_command_argument(6,arg)
        read(arg,*) step_factor
    end if
    if(command_argument_count()>=7) then
        call get_command_argument(7,arg)
        read(arg,*) energy_tol
    end if
    print '(a,es12.4)', 'velocity multiplier=',vscale
    omega=0; eye=0
    do i=1,6
        eye(i,i)=1
    end do
    omega(1:3,4:6)=eye(1:3,1:3)
    omega(4:6,1:3)=-eye(1:3,1:3)
    do case_id=1,2
        if((failure_mode.or.only_marker>=0).and.case_id==1) cycle
        call get_command_argument(case_id,filename)
        call create_spectre_field(field,trim(filename),ierr)
        if(ierr/=0) error stop 'fixture load failed'
        call set_magfie_spectre_field(field)
        call gauss2_prepare
        call reference_prepare
        print '(a,i0,2a)', 'CASE ',case_id,' ',trim(filename)
        if(failure_mode) then
            call check_newton_boundary
            cycle
        end if
        if(only_marker>=0) then
            call energy_run(nlong,only_marker)
            cycle
        end if
        call check_derivatives
        call check_monitor_derivatives
        efixed=0; lfixed=1; gfixed=1; mfixed=0
        call check_map(.false.,0.65_dp)
        call check_map(.true.,0.65_dp)
        if(case_id==1) then
            call energy_run(min(nlong,1000),0)
        else
            call energy_run(nlong,0)
            call energy_run(nlong,1)
            call energy_run(nlong,2)
            if(index(trim(filename),'ncsx')>0.and.abs(vscale-1)<0.1_dp) then
                call energy_run(nlong,9)
                call energy_run(nlong,20)
                call energy_run(nlong,21)
                call energy_run(nlong,120)
            end if
        end if
        width=0.06_dp;ro0=1.0e5_dp;h=0.02_dp
        lfixed=0.12_dp; gfixed=0.1_dp
        call check_map(.false.,0.995_dp)
        call check_map(.true.,0.995_dp)
    end do
    if(energy_failed) error stop 'one or more energy runs failed; see diagnostics'
contains
    subroutine initial(z,rho)
        real(dp),intent(out)::z(6)
        real(dp),intent(in)::rho
        real(dp)::x(3),e(3,3),g(3,3),dg(3,3,3),hg(3,3,3,3)
        real(dp)::av(3),da(3,3),ha(3,3,3)
        z(1:3)=[rho,0.7_dp,0.3_dp]
        call gauss2_geometry(z(1:3),width,x,e,g,dg,hg)
        call gauss2_potential(z(1:3),width,ro0,av,da,ha)
        z(4:6)=av+matmul(transpose(e),vscale*[0.001_dp,0.002_dp,0.001_dp])
    end subroutine initial

    subroutine check_derivatives
        real(dp)::z(6),zp(6),zm(6),fp(6),fm(6),f(6),df(6,6),dfr(6,6)
        real(dp)::en,er,grad(6),hh(6,6),errj,errh,eps,dnum(6,6)
        real(dp)::ep,em,ep2,em2,eg(6),gradfd(6),errgrad,errvalue,errrhs
        real(dp)::rhos(5)=[0.65_dp,0.98_dp,1.0_dp,1.02_dp,1.35_dp]
        integer::r,j,refinement
        do r=1,size(rhos)
            call initial(z,rhos(r))
            call generated_eval(z,width,ro0,en,f,df)
            call gauss2_hamiltonian(z,width,ro0,er,grad,hh)
            dfr(1:3,:)=hh(4:6,:); dfr(4:6,:)=-hh(1:3,:)
            eg(1:3)=grad(4:6); eg(4:6)=-grad(1:3)
            errrhs=maxval(abs(f-eg))/max(1.0e-10_dp,maxval(abs(f)))
            errvalue=abs(en-er)/max(1.0e-10_dp,abs(en))
            print '(a,2es12.4)', 'VALUE/RHS generated-reference=',errvalue,errrhs
            if(max(errvalue,errrhs)>1.0e-10_dp) &
                error stop 'generated/reference value or RHS mismatch'
            errh=maxval(abs(df-dfr))/max(1.0_dp,maxval(abs(df)))
            if(errh>1.0e-10_dp) error stop 'generated/reference Hessian mismatch'
            do refinement=1,3
            do j=1,6
                eps=1.0e-4_dp*0.25_dp**(refinement-1)*max(1.0_dp,abs(z(j)))
                zp=z; zm=z; zp(j)=zp(j)+eps; zm(j)=zm(j)-eps
                call generated_eval(zp,width,ro0,ep,fp,dfr)
                call generated_eval(zm,width,ro0,em,fm,dfr)
                dnum(:,j)=8*(fp-fm)
                zp=z; zm=z; zp(j)=zp(j)+2*eps; zm(j)=zm(j)-2*eps
                call generated_eval(zp,width,ro0,ep2,fp,dfr)
                call generated_eval(zm,width,ro0,em2,fm,dfr)
                gradfd(j)=(8*(ep-em)-(ep2-em2))/(12*eps)
                dnum(:,j)=(dnum(:,j)-(fp-fm))/(12*eps)
            end do
            eg(1:3)=-f(4:6); eg(4:6)=f(1:3)
            errgrad=maxval(abs(gradfd-eg))/max(1.0e-10_dp,maxval(abs(eg)))
            print '(a,i0,a,es12.4)', 'FD refinement=',refinement, &
                ' energy/gradient=',errgrad
            if(errgrad>1.0e-6_dp.and.refinement==3) &
                error stop 'energy/RHS inconsistency'
            errj=maxval(abs(dnum-df))/max(1.0_dp,maxval(abs(df)))
            print '(a,f5.2,2(a,es12.4))','DERIV rho=',rhos(r), &
                ' generated/reference=',errh,' FD/RHS=',errj
            if(errj>1.0e-6_dp.and.refinement==3) &
                error stop 'field-jet/RHS derivative mismatch'
            end do
        end do
    end subroutine check_derivatives

    real(dp) function reference_k(z)
        real(dp),intent(in)::z(6)
        real(dp)::en,grad(6),hh(6,6),radial,t
        integer::k
        call gauss2_hamiltonian(z,width,ro0,en,grad,hh,radial_metric=radial)
        t=1
        do k=1,spectre_field%data%Mvol-1
            t=t*(1-(1-gfixed)*exp(-((z(1)-real(k,dp))/lfixed)**2))
        end do
        reference_k=t*sqrt(mfixed/(mfixed+radial))*(en-efixed)
    end function reference_k

    subroutine check_monitor_derivatives
        real(dp)::z(6),zp(6),zm(6),f(6),df(6,6),fp(6),fm(6),scratch(6,6)
        real(dp)::en,ke,ep,em,ep2,em2,gradfd(6),dnum(6,6),grad(6)
        real(dp)::eps,radial,hh(6,6),eg(6),errgrad,errj,errvalue
        real(dp),parameter::rhos(3)=[0.98_dp,0.995_dp,1.02_dp]
        integer::r,j,refinement
        lfixed=0.12_dp;gfixed=0.1_dp
        do r=1,size(rhos)
            call initial(z,rhos(r))
            call gauss2_hamiltonian(z,width,ro0,en,grad,hh,radial_metric=radial)
            ! Exercise all second-monitor terms, including (H-E0) != 0.
            efixed=0.9_dp*en;mfixed=radial
            call generated_eval(z,width,ro0,en,f,df,efixed,lfixed,gfixed,mfixed,ke)
            errvalue=abs(ke-reference_k(z))/max(1.0e-10_dp,abs(ke))
            if(errvalue>1.0e-10_dp) error stop 'monitored K value mismatch'
            do refinement=1,4
                do j=1,6
                    eps=1.0e-4_dp*0.25_dp**(refinement-1)*max(1.0_dp,abs(z(j)))
                    zp=z;zm=z;zp(j)=zp(j)+eps;zm(j)=zm(j)-eps
                    ep=reference_k(zp);em=reference_k(zm)
                    call generated_eval(zp,width,ro0,en,fp,scratch, &
                        efixed,lfixed,gfixed,mfixed)
                    call generated_eval(zm,width,ro0,en,fm,scratch, &
                        efixed,lfixed,gfixed,mfixed)
                    dnum(:,j)=8*(fp-fm)
                    zp=z;zm=z;zp(j)=zp(j)+2*eps;zm(j)=zm(j)-2*eps
                    ep2=reference_k(zp);em2=reference_k(zm)
                    call generated_eval(zp,width,ro0,en,fp,scratch, &
                        efixed,lfixed,gfixed,mfixed)
                    call generated_eval(zm,width,ro0,en,fm,scratch, &
                        efixed,lfixed,gfixed,mfixed)
                    gradfd(j)=(8*(ep-em)-(ep2-em2))/(12*eps)
                    dnum(:,j)=(dnum(:,j)-(fp-fm))/(12*eps)
                end do
                eg(1:3)=-f(4:6);eg(4:6)=f(1:3)
                errgrad=maxval(abs(gradfd-eg))/max(1.0e-10_dp,maxval(abs(eg)))
                errj=maxval(abs(dnum-df))/max(1.0_dp,maxval(abs(df)))
                print '(a,f6.3,a,i0,3(a,es12.4))', 'MONITOR rho=',rhos(r), &
                    ' refinement=',refinement,' value=',errvalue, &
                    ' K/gradient=',errgrad,' RHS/Jacobian=',errj
                if(refinement==4.and.max(errgrad,errj)>1.0e-6_dp) &
                    error stop 'monitored kernel derivative mismatch'
            end do
        end do
    end subroutine check_monitor_derivatives

    subroutine step(z,tol,reference,out,m,res)
        real(dp),intent(in)::z(6),tol
        logical,intent(in)::reference
        real(dp),intent(out)::out(6),m(6,6),res
        real(dp)::dt
        logical::ok
        if(reference) then
            call gauss2_step(z,width,ro0,efixed,lfixed,gfixed,h,out,dt,ok,m,radial_scale=mfixed)
            res=-1 ! Reference API does not expose stage residual.
        else
            call generated_step(z,h,width,ro0,tol,out,m,res,ok,efixed,lfixed,gfixed,mfixed)
        end if
        if(.not.ok) then
            print '(a,i0,a,l1,2(a,es16.8))','NEWTON_FAILURE step=',diagnostic_step, &
                ' reference=',reference,' tol=',tol,' residual=',res
            print '(a,6es25.17)', 'NEWTON_FAILURE z=',z
            print '(a,6(1x,z16.16))', 'NEWTON_FAILURE bits=',transfer(z,[0_int64],6)
            call diagnose_newton(z)
            error stop 'Newton failure'
        end if
    end subroutine step

    subroutine check_newton_boundary
        real(dp)::z(6),start(6),out(6),m(6,6),res,f(6),df(6,6)
        logical::ok
        integer(int64),parameter::bits(6)=[int(z'3FEFA155692841DE',int64), &
            int(z'C004455BB6C61D05',int64),int(z'40208B3C9ECD7D3C',int64), &
            int(z'BFD9B36E018FD46F',int64),int(z'4012DCDA91AFF80F',int64), &
            int(z'40354681F5654DA0',int64)]
        call initial(start,0.65_dp)
        call generated_eval(start,width,ro0,efixed,f,df)
        z=transfer(bits,z)
        lfixed=1; gfixed=1
        mfixed=(lfixed/(16*h*sqrt(2*efixed)))**2
        call diagnose_newton(z)
        call diagnose_chart(z)
        ! A prescribed half-sized step is a different map; this test does not
        ! change the failing long run or use state-dependent adaptation.
        call generated_step(z,h/2,width,ro0,1.0e-12_dp,out,m,res,ok, &
            efixed,lfixed,gfixed,mfixed)
        if(.not.ok) error stop 'prescribed half step failed at saved state'
        if(res>1.0e-12_dp) error stop 'half step accepted excessive residual'
        if(defect(m)>1.0e-7_dp) error stop 'half-step tangent defect'
        call generated_step(z,gauss2_entry_step(h,width,lfixed,gfixed),width,ro0, &
            1.0e-12_dp,out,m,res,ok,efixed,lfixed,gfixed,mfixed)
        print '(a,l1,a,es16.8)', 'ENTRY_POLICY saved-state converged=',ok, &
            ' residual=',res
        if(.not.ok.or.res>1.0e-12_dp) error stop 'entry policy saved-state failure'
        if(defect(m)>1.0e-7_dp) error stop 'entry policy tangent defect'
    end subroutine check_newton_boundary

    subroutine diagnose_newton(z)
        real(dp),intent(in)::z(6)
        real(dp)::out(6),m(6,6),res,dt
        logical::ok
        integer::j
        do j=1,3
            call generated_step(z,h/real(2**(j-1),dp),width,ro0,1.0e-12_dp, &
                out,m,res,ok,efixed,lfixed,gfixed,mfixed)
            print '(a,es12.4,a,l1,a,es12.4)', 'NEWTON_PROBE h=', &
                h/real(2**(j-1),dp),' converged=',ok,' residual=',res
        end do
        call diagnose_kernel(z)
        call gauss2_step(z,width,ro0,efixed,lfixed,gfixed,h,out,dt,ok,m,radial_scale=mfixed)
        print '(a,l1)', 'NEWTON_PROBE frozen reference converged=',ok
    end subroutine diagnose_newton

    subroutine diagnose_chart(z)
        real(dp),intent(in)::z(6)
        real(dp)::q(3),x(3),basis(3,3),metric(3,3),dg(3,3,3),d2g(3,3,3,3)
        real(dp)::det,min_det,max_det,minrho,energy,f(6),df(6,6)
        integer::j
        min_det=huge(1.0_dp);max_det=-huge(1.0_dp);minrho=0
        do j=0,1000
            q=z(1:3);q(1)=1-width/2+width*real(j,dp)/1000
            call gauss2_geometry(q,width,x,basis,metric,dg,d2g)
            det=dot_product(basis(:,1),cross(basis(:,2),basis(:,3)))
            if(det<min_det) minrho=q(1)
            min_det=min(min_det,det);max_det=max(max_det,det)
        end do
        call generated_eval(z,width,ro0,energy,f,df)
        print '(a,4es24.16)', 'CHART min_det max_det rho_at_min H=', &
            min_det,max_det,minrho,energy
    end subroutine diagnose_chart

    subroutine diagnose_kernel(z)
        real(dp),intent(in)::z(6)
        real(dp)::en,f(6),df(6,6),fp(6),fm(6),scratch(6,6),fd(6,6)
        real(dp)::zp(6),zm(6),eps,ep,em,ep2,em2,gradfd(6),eg(6),errg,errj
        integer::j,k
        call generated_eval(z,width,ro0,en,f,df,efixed,lfixed,gfixed,mfixed)
        eg(1:3)=-f(4:6);eg(4:6)=f(1:3)
        do k=1,3
            do j=1,6
                eps=1.0e-5_dp*0.25_dp**(k-1)*max(1.0_dp,abs(z(j)))
                zp=z;zm=z;zp(j)=zp(j)+eps;zm(j)=zm(j)-eps
                ep=reference_k(zp);em=reference_k(zm)
                call generated_eval(zp,width,ro0,en,fp,scratch, &
                    efixed,lfixed,gfixed,mfixed)
                call generated_eval(zm,width,ro0,en,fm,scratch, &
                    efixed,lfixed,gfixed,mfixed)
                fd(:,j)=8*(fp-fm)
                zp=z;zm=z;zp(j)=zp(j)+2*eps;zm(j)=zm(j)-2*eps
                ep2=reference_k(zp);em2=reference_k(zm)
                call generated_eval(zp,width,ro0,en,fp,scratch, &
                    efixed,lfixed,gfixed,mfixed)
                call generated_eval(zm,width,ro0,en,fm,scratch, &
                    efixed,lfixed,gfixed,mfixed)
                gradfd(j)=(8*(ep-em)-(ep2-em2))/(12*eps)
                fd(:,j)=(fd(:,j)-(fp-fm))/(12*eps)
            end do
            errg=maxval(abs(gradfd-eg))/max(1.0e-10_dp,maxval(abs(eg)))
            errj=maxval(abs(fd-df))/max(1.0_dp,maxval(abs(df)))
            print '(a,i0,2(a,es12.4))', 'FAILURE_KERNEL refinement=',k, &
                ' K/gradient=',errg,' RHS/Jacobian=',errj
        end do
        if(max(errg,errj)>1.0e-6_dp) &
            error stop 'kernel derivative mismatch at failed Newton state'
    end subroutine diagnose_kernel

    subroutine map(z,tol,reference,n,out,m,res)
        real(dp),intent(in)::z(6),tol
        logical,intent(in)::reference
        integer,intent(in)::n
        real(dp),intent(out)::out(6),m(6,6),res
        real(dp)::current(6),next(6),mi(6,6),ri
        integer::s
        current=z; m=eye; res=0
        do s=1,n
            call step(current,tol,reference,next,mi,ri)
            m=matmul(mi,m); res=max(res,ri); current=next
        end do
        out=current
        if(reference) res=-1 ! Unavailable from frozen reference API.
    end subroutine map

    real(dp) function defect(m)
        real(dp),intent(in)::m(6,6)
        defect=sqrt(sum((matmul(transpose(m),matmul(omega,m))-omega)**2))
    end function defect

    subroutine finite_difference(z,tol,reference,n,eps0,fd)
        real(dp),intent(in)::z(6),tol,eps0
        logical,intent(in)::reference
        integer,intent(in)::n
        real(dp),intent(out)::fd(6,6)
        real(dp)::zp(6),zm(6),fp(6),fm(6),m(6,6),res,eps
        integer::j
        do j=1,6
            eps=eps0*max(1.0_dp,abs(z(j)))
            zp=z; zm=z; zp(j)=zp(j)+eps; zm(j)=zm(j)-eps
            call map(zp,tol,reference,n,fp,m,res)
            call map(zm,tol,reference,n,fm,m,res)
            fd(:,j)=8*(fp-fm)
            zp=z; zm=z; zp(j)=zp(j)+2*eps; zm(j)=zm(j)-2*eps
            call map(zp,tol,reference,n,fp,m,res)
            call map(zm,tol,reference,n,fm,m,res)
            fd(:,j)=(fd(:,j)-(fp-fm))/(12*eps)
        end do
    end subroutine finite_difference

    subroutine check_map(reference,rho)
        logical,intent(in)::reference
        real(dp),intent(in)::rho
        real(dp)::energy,f0(6),d0(6,6),outref(6),mref(6,6),refres,comparison
        real(dp)::z(6),out(6),m(6,6),fd(6,6),fd2(6,6),res,tol,noise
        real(dp),parameter::tols(3)=[1.0e-7_dp,1.0e-10_dp,1.0e-12_dp]
        integer::k,n,ncase
        h=0.02_dp
        call initial(z,rho)
        call generated_eval(z,width,ro0,energy,f0,d0)
        efixed=energy
        mfixed=(lfixed/(16*h*sqrt(2*efixed)))**2
        if(step_factor>0) then
            h=h*step_factor
        else
            h=gauss2_entry_step(h,width,lfixed,gfixed)
        end if
        print '(a,f7.4,a,f6.3)', 'MAP configuration rho=',rho,' gmin=',gfixed
        do k=1,3
            if(reference.and.k==2) cycle
            tol=tols(k)
            if(reference) tol=-1 ! Reference Newton tolerance is fixed internally.
            do ncase=1,2
            n=1
            if(ncase==2) n=1000
            call map(z,tol,reference,n,out,m,res)
            if(.not.reference) then
                call map(z,tol,.true.,n,outref,mref,refres)
                comparison=maxval(abs(out-outref))/max(1.0_dp,maxval(abs(out)))
                print '(a,i0,a,es12.4)', 'MAP correspondence steps=',n, &
                    ' relative max error=',comparison
                if(n==1.or.k==3) then
                    if(comparison>1.0e-8_dp) error stop 'compiled step correspondence failure'
                end if
            end if
            call finite_difference(z,tol,reference,n,1.0e-4_dp,fd)
            call finite_difference(z,tol,reference,n,5.0e-5_dp,fd2)
            noise=sqrt(sum((fd-fd2)**2))
            if(n==1000.and.k==3.and.(noise>1.0e-5_dp.or.defect(fd2)>1.0e-5_dp)) then
                print '(a,2es16.8)', 'FD_REFINEMENT coarse defect/disagreement=', &
                    defect(fd2),noise
                call finite_difference(z,tol,reference,n,2.5e-5_dp,fd)
                call finite_difference(z,tol,reference,n,1.25e-5_dp,fd2)
                noise=sqrt(sum((fd-fd2)**2))
                if(noise>1.0e-5_dp) error stop 'FD refinement unresolved'
            end if
            print '(a,l1,a,i0,5(a,es12.4))','MAP reference=',reference, &
                ' steps=',n,' tol=',tol,' residual=',res,' tangent defect=', &
                defect(m),' FD defect=',defect(fd2),' FD disagreement=',noise
            if(defect(m)>1.0e-7_dp) error stop 'tangent symplecticity defect'
            ! This finite-difference threshold is fixed, separately reported
            ! from the Newton tolerance and refinement disagreement.
            if(n==1000.and.k<3) then
                ! Extra tolerance sweep: preserve failures of the strict criterion
                ! as observations. The original strict gates apply to 1-step
                ! maps and the 1000-step map at 1e-12, unchanged.
                if(defect(fd2)>1.0e-5_dp) &
                    print '(a)', 'LOOSE_TOLERANCE_MAP exceeds strict FD criterion'
            else
                if(defect(fd2)>1.0e-5_dp) error stop 'compiled map symplecticity defect'
            end if
            end do
        end do
    end subroutine check_map

    subroutine energy_run(n,marker)
        integer,intent(in)::n,marker
        real(dp)::z(6),out(6),res,e0,en,f(6),df(6,6)
        real(dp)::emin,emax,eblock_min,eblock_max,rhomin,rhomax
        real(dp)::xx,yy,mx,my,sxx,sxy,syy,dx,dy,slope,se
        real(dp)::maxres
        integer::s
        logical::ok
        width=0.06_dp;ro0=1.0e5_dp;h=0.02_dp;lfixed=1;gfixed=1
        select case(marker)
        case(0)
            call initial(z,0.65_dp)
        case(1)
            call initial(z,0.95_dp)
        case(2)
            call initial(z,0.995_dp)
            lfixed=0.12_dp;gfixed=0.1_dp
        case(9,20,21)
            call validation_initial(marker,z)
        case(120)
            ro0=sqrt(2*3500.0_dp*ev/(4*p_mass))* &
                4*p_mass*c/(2*e_charge)/sqrt(2.0_dp)
            width=1.8214223148575103e-3_dp
            h=3.7552110327543201e-2_dp
            lfixed=5.8285514075440331e-2_dp;gfixed=1.953125e-3_dp
            z=[1.0025024863675955_dp,3.2758489329512472_dp, &
                -3.9738344980129323_dp,-13.708571643795706_dp, &
                84.369367121529450_dp,-16.358194458166761_dp]
        end select
        call generated_eval(z,width,ro0,e0,f,df)
        efixed=e0
        mfixed=(lfixed/(16*h*sqrt(2*efixed)))**2
        if(step_factor>0) then
            h=h*step_factor
        else
            h=gauss2_entry_step(h,width,lfixed,gfixed)
        end if
        print '(a,i0,a,6es24.16)', 'START marker=',marker,' z=',z
        print '(a,6es24.16)', 'PARAM width ro0 h length gmin mref=', &
            width,ro0,h,lfixed,gfixed,mfixed
        emin=e0; emax=e0; eblock_min=e0; eblock_max=e0
        rhomin=z(1);rhomax=z(1);mx=0;my=0;sxx=0;sxy=0;syy=0;maxres=0
        do s=1,n
            diagnostic_step=s
            call shipped_step(z,width,ro0,efixed,lfixed,gfixed,h,out,xx,ok, &
                radial_scale=mfixed,tolerance=energy_tol,residual=res)
            if(.not.ok) then
                print '(a,i0,a,i0,2(a,es16.8))', 'ENERGY_FAILURE marker=',marker, &
                    ' step=',s,' rho=',z(1),' residual=',res
                print '(a,6es25.17)', 'ENERGY_FAILURE z=',z
                call diagnose_newton(z)
                energy_failed=.true.
                return
            end if
            maxres=max(maxres,res);z=out
            call generated_eval(z,width,ro0,en,f,df)
            emin=min(emin,en); emax=max(emax,en)
            rhomin=min(rhomin,z(1));rhomax=max(rhomax,z(1))
            eblock_min=min(eblock_min,en); eblock_max=max(eblock_max,en)
            xx=real(s,dp)/real(n,dp);yy=en/e0-1
            dx=xx-mx;dy=yy-my
            mx=mx+dx/s;my=my+dy/s
            sxx=sxx+dx*(xx-mx);sxy=sxy+dx*(yy-my);syy=syy+dy*(yy-my)
            if(mod(s,max(1,n/10))==0) then
                print '(a,i0,a,i0,3(a,es12.4),a,f8.4)', &
                    'ENERGY marker=',marker,' steps=',s, &
                    ' min/E0-1=',eblock_min/e0-1,' max/E0-1=',eblock_max/e0-1, &
                    ' end/E0-1=',yy,' rho=',z(1)
                eblock_min=en; eblock_max=en
            end if
        end do
        slope=0;se=0
        if(n>2) then
            slope=sxy/sxx/real(n,dp)
            se=sqrt(max(0.0_dp,syy-sxy*sxy/sxx)/real(n-2,dp)/sxx)/real(n,dp)
        end if
        print '(a,i0,a,i0,7(a,es16.8))','ENVELOPE marker=',marker,' steps=',n, &
            ' min/E0-1=',emin/e0-1,' max/E0-1=',emax/e0-1, &
            ' slope/step=',slope,' OLS_SE/step=',se, &
            ' rho_min=',rhomin,' rho_max=',rhomax,' max_residual=',maxres
        if(max(abs(emin/e0-1),abs(emax/e0-1))>1.0e-5_dp) then
            print '(a)', 'ENERGY_FAILURE envelope exceeds unchanged 1e-5 criterion'
            energy_failed=.true.
        end if
    end subroutine energy_run

    subroutine validation_initial(marker,z)
        ! Recover the validation GC seeds; construct a gyrophase-zero FO
        ! entry using the same displacement and parameter formulas as production.
        ! These are starts at the seed, not later hybrid excursion checkpoints.
        integer,intent(in)::marker
        real(dp),intent(out)::z(6)
        real(dp)::x(3),basis(3,3),g(3,3),dg(3,3,3),d2g(3,3,3,3)
        real(dp)::av(3),da(3,3),d2a(3,3,3),b(3),axis(3),perp(3),v(3)
        real(dp)::pitch,bmod,det,radial,energy,grad(6),hh(6,6),speed
        logical::ok
        z(1)=0.95_dp
        select case(marker)
        case(9)
            z(2:3)=[5.7432240698438353_dp,4.0742529726242589_dp]
            pitch=0.34030985461878194_dp
        case(20)
            z(2:3)=[3.2888548092268115_dp,0.73631077818510704_dp]
            pitch=-0.18109296850536349_dp
        case(21)
            z(2:3)=[5.9395736106931967_dp,4.4669520543229826_dp]
            pitch=0.31788627649272905_dp
        end select
        ro0=sqrt(2*3500.0_dp*ev/(4*p_mass))* &
            4*p_mass*c/(2*e_charge)/sqrt(2.0_dp)
        call gauss2_geometry(z(1:3),0.0_dp,x,basis,g,dg,d2g)
        call gauss2_potential(z(1:3),0.0_dp,ro0,av,da,d2a)
        det=dot_product(basis(:,1),cross(basis(:,2),basis(:,3)))
        b=ro0*matmul(basis,[da(3,2)-da(2,3),da(1,3)-da(3,1), &
            da(2,1)-da(1,2)])/det
        bmod=norm2(b);b=b/bmod
        axis=[0.0_dp,0.0_dp,1.0_dp]
        if(abs(b(3))>=0.9_dp) axis=[1.0_dp,0.0_dp,0.0_dp]
        perp=axis-dot_product(axis,b)*b;perp=perp/norm2(perp)
        speed=sqrt(2.0_dp)
        v=speed*(pitch*b+sqrt(1-pitch**2)*perp)
        x=x+ro0/bmod*cross(b,speed*sqrt(1-pitch**2)*perp)
        z(4:6)=av
        call gauss2_hamiltonian(z,0.0_dp,ro0,energy,grad,hh,radial_metric=radial)
        radial=speed*sqrt(radial)
        width=0.05_dp*ro0/bmod*radial;h=0.1_dp*ro0/bmod
        lfixed=max(16*h*radial,16*width)
        gfixed=min(0.1_dp,width/(256*h*radial))
        call gauss2_encode(x,v,width,ro0,z,ok)
        if(.not.ok) error stop 'validation seed FO encode failed'
    end subroutine validation_initial

    pure function cross(a,b) result(out)
        real(dp),intent(in)::a(3),b(3)
        real(dp)::out(3)
        out=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
    end function cross
end program test_fo_symplectic
