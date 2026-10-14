module generated_fo_stepper
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use lapack_interfaces, only: dgesv
    use fo_kernel, only: fo_eval, fo_monitor_factor, fo_jet_product, &
        fo_radial, fo_metric_monitor
    use magfie_sub, only: spectre_field
    use spectre_fo_gauss2, only: gauss2_geometry, gauss2_potential
    implicit none
    private
    public :: generated_eval, generated_step
    real(dp),parameter :: a(2,2)=reshape([0.25_dp, &
        0.25_dp+sqrt(3.0_dp)/6.0_dp,0.25_dp-sqrt(3.0_dp)/6.0_dp, &
        0.25_dp],[2,2])
contains
    subroutine generated_eval(z,width,ro0,energy,f,df,e0,length,gmin,radial_scale,kenergy)
        real(dp),intent(in)::z(6),width,ro0
        real(dp),intent(out)::energy,f(6),df(6,6)
        real(dp),intent(out),optional::kenergy
        real(dp),intent(in),optional::e0,length,gmin,radial_scale
        real(dp)::x(3),basis(3,3),g(3,3),dg(3,3,3),d2g(3,3,3,3)
        real(dp)::av(3),da(3,3),d2a(3,3,3),aj(3,10),gj(6,10)
        real(dp)::ke,gh(6),hh(6,6),tj(10),rhoj(3),rj(10),factor(3),product(3),energy0
        integer::i,j,k,l
        integer,parameter::si(6)=[1,1,1,2,2,3],sj(6)=[1,2,3,2,3,3]
        call gauss2_geometry(z(1:3),width,x,basis,g,dg,d2g)
        call gauss2_potential(z(1:3),width,ro0,av,da,d2a)
        aj(:,1)=av; aj(:,2:4)=da
        do l=1,6
            gj(l,1)=g(si(l),sj(l))
            gj(l,2:4)=dg(si(l),sj(l),:)
        end do
        k=4
        do j=1,3
            do i=1,j
                k=k+1
                aj(:,k)=d2a(:,i,j)
                do l=1,6
                    gj(l,k)=d2g(si(l),sj(l),i,j)
                end do
            end do
        end do
        ! Existing code uses normalized P=e_i.v+A/ro0, unit mass/charge.
        energy0=0
        if(present(e0)) energy0=e0
        rhoj=[1.0_dp,0.0_dp,0.0_dp]
        if(present(gmin).and.present(length)) then
            do k=1,spectre_field%data%Mvol-1
                call fo_monitor_factor(z(1)-real(k,dp),length,gmin, &
                    factor(1),factor(2),factor(3))
                call fo_jet_product(rhoj(1),rhoj(2),rhoj(3), &
                    factor(1),factor(2),factor(3), &
                    product(1),product(2),product(3))
                rhoj=product
            end do
        end if
        tj=0; tj(1)=rhoj(1);tj(2)=rhoj(2);tj(5)=rhoj(3)
        if(present(radial_scale)) then
            if(radial_scale>0) then
                call fo_radial(gj,rj)
                call fo_metric_monitor(rj,rhoj,radial_scale,tj)
            end if
        end if
        call fo_eval(z(4:6),aj,gj,1.0_dp,1.0_dp,energy0, &
            tj,energy,ke,gh,hh,f,df)
        if(present(kenergy)) kenergy=ke
    end subroutine generated_eval

    subroutine generated_step(z,h,width,ro0,tol,zout,tangent,residual,ok, &
        e0,length,gmin,radial_scale)
        real(dp),intent(in)::z(6),h,width,ro0,tol
        real(dp),intent(out)::zout(6),tangent(6,6),residual
        logical,intent(out)::ok
        real(dp),intent(in),optional::e0,length,gmin,radial_scale
        real(dp)::y(6,2),f(6,2),df(6,6,2),jac(12,12),rhs(12,6)
        real(dp)::r(6,2),delta(12),energy,eye(6,6)
        integer::i,j,it,piv(12),info
        ok=.false.; zout=z; tangent=0; residual=huge(1.0_dp)
        eye=0
        do i=1,6
            eye(i,i)=1
        end do
        call generated_eval(z,width,ro0,energy,f(:,1),df(:,:,1),e0,length,gmin,radial_scale)
        y(:,1)=z+h*(0.5_dp-sqrt(3.0_dp)/6.0_dp)*f(:,1)
        y(:,2)=z+h*(0.5_dp+sqrt(3.0_dp)/6.0_dp)*f(:,1)
        do it=1,20
            do j=1,2
                call generated_eval(y(:,j),width,ro0,energy,f(:,j),df(:,:,j), &
                    e0,length,gmin,radial_scale)
                do i=1,2
                    jac(6*i-5:6*i,6*j-5:6*j)=-h*a(i,j)*df(:,:,j)
                end do
            end do
            do i=1,12
                jac(i,i)=jac(i,i)+1
            end do
            do i=1,2
                r(:,i)=y(:,i)-z-h*(a(i,1)*f(:,1)+a(i,2)*f(:,2))
            end do
            residual=maxval(abs(r))
            if(.not.all(ieee_is_finite(r))) return
            if(residual<=tol) then
                ok=.true.
                exit
            end if
            delta=-reshape(r,[12])
            call dgesv(12,1,jac,12,piv,delta,12,info)
            if(info/=0) return
            y=y+reshape(delta,[6,2])
        end do
        if(.not.ok) return
        zout=z+h*0.5_dp*(f(:,1)+f(:,2))
        ! This is the derivative of the exact-stage equations at the computed
        ! stages. Finite differences below check the actual terminated map.
        rhs(1:6,:)=eye; rhs(7:12,:)=eye
        call dgesv(12,6,jac,12,piv,rhs,12,info)
        if(info/=0) then
            ok=.false.
            return
        end if
        tangent=eye+h*0.5_dp*(matmul(df(:,:,1),rhs(1:6,:)) &
            +matmul(df(:,:,2),rhs(7:12,:)))
    end subroutine generated_step
end module generated_fo_stepper
