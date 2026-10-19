program gen_fo_hamiltonian
    use fortsym_string, only: str, chars, str_t
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**), func, exp, sqrt
    use fortsym_diff, only: diff, partial_derivative
    use fortsym_subs, only: subs_many
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE, &
        cse_result_t, cse_analyse
    use fortsym_fparse, only: parse_fortran_expr
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: eng
    type(expr_t) :: q(3), p(3), z(6), av(3), g(3,3)
    type(expr_t) :: gi(3,3), w(3), det, ham, kham, monitor, e0, mass, charge
    type(expr_t) :: old(100), replacement(100)
    type(expr_t),allocatable :: temp_symbols(:),temp_values(:)
    type(cse_result_t) :: cse
    type(expr_t) :: first_raw(12)
    type(expr_t) :: roots(86), raw(86), tmp, dg(6)
    type(engine_result_t) :: res
    type(kernel_spec_t) :: spec
    type(suite_t) :: suite
    type(str_t) :: args(106), outputs(86)
    character(1024)::outroot
    character(32) :: name
    character(:),allocatable::message,first_code
    logical::ok
    integer :: i,j,k,l,n,u,ii,jj,idx,nrep,order
    integer, parameter :: si(6)=[1,1,1,2,2,3], sj(6)=[1,2,3,2,3,3]
    outroot='.'
    call get_environment_variable('FO_GENERATOR_OUTPUT',outroot,status=i)
    if(i/=0) outroot='.'
    call arena%init()
    eng=make_symengine_engine(arena)
    if(command_argument_count()>0) stop
    n=0
    do i=1,3
        write(name,'(a,i0)') 'q',i
        q(i)=sym(arena,trim(name))
        write(name,'(a,i0)') 'p',i
        p(i)=input(trim(name))
    end do
    call suite_begin(suite,"FO scalar generator exact identities")
    call emit_monitor
    mass=input('mass'); charge=input('charge'); e0=input('e0')
    nrep=0
    do i=1,3
        write(name,'(a,i0)') 'Afield',i
        av(i)=func(trim(name),q)
        call jet_replacements(av(i),'a',i)
        w(i)=p(i)-charge*av(i)
    end do
    do l=1,6
        write(name,'(a,i0)') 'Gfield',l
        tmp=func(trim(name),q)
        call jet_replacements(tmp,'g',l)
        g(si(l),sj(l))=tmp
        g(sj(l),si(l))=tmp
    end do
    ! Adjugate/determinant inversion of the symmetric covariant metric.
    det=g(1,1)*g(2,2)*g(3,3)+2*g(1,2)*g(1,3)*g(2,3) &
        -g(1,1)*g(2,3)**2-g(2,2)*g(1,3)**2-g(3,3)*g(1,2)**2
    do j=1,3
        do i=1,3
            ii=mod(i,3)+1; jj=mod(j,3)+1
            gi(j,i)=(g(ii,jj)*g(mod(ii,3)+1,mod(jj,3)+1) &
                -g(ii,mod(jj,3)+1)*g(mod(ii,3)+1,jj))/det
        end do
    end do
    ham=num(arena,0)
    do j=1,3
        do i=1,3
            ham=ham+w(i)*gi(i,j)*w(j)/(2*mass)
        end do
    end do
    call emit_radial
    monitor=func('monitor',q)
    call jet_replacements(monitor,'t',1)
    if(n/=106) error stop 'input count'
    kham=monitor*(ham-e0)
    z(1:3)=q; z(4:6)=p
    raw(1)=ham; outputs(1)=str('energy')
    raw(2)=kham; outputs(2)=str('kenergy')
    do i=1,6
        print *, "first derivative",i
        raw(2+i)=diff(ham,z(i))
        first_raw(i)=raw(2+i)
        res=eng%simplify(raw(2+i))
        if(.not.res%ok) error stop "first derivative simplify refused"
        raw(2+i)=res%value
        write(name,'(a,i0)') 'gh',i
        outputs(2+i)=str(trim(name))
        dg(i)=diff(kham,z(i))
        first_raw(6+i)=dg(i)
        res=eng%simplify(dg(i))
        if(.not.res%ok) error stop "K gradient simplify refused"
        dg(i)=res%value
    end do
    do i=1,6
        k=mod(i+2,6)+1
        raw(44+i)=dg(k)
        if(i>3) raw(44+i)=-raw(44+i)
        write(name,'(a,i0)') 'f',i
        outputs(44+i)=str(trim(name))
    end do
    do j=1,6
        print *, "second derivative column",j
        do i=1,6
            k=8+(j-1)*6+i
            raw(k)=diff(raw(2+i),z(j))
            write(name,'(a,i0,a,i0)') 'hh',i,'_',j
            outputs(k)=str(trim(name))
            k=50+(j-1)*6+i
            raw(k)=diff(raw(44+i),z(j))
            write(name,'(a,i0,a,i0)') 'df',i,'_',j
            outputs(k)=str(trim(name))
        end do
    end do
    do k=1,86
        raw(k)=subs_many(raw(k),old,replacement)
        res=eng%simplify(raw(k))
        if(.not.res%ok) error stop 'simplification refused'
        roots(k)=res%value
        write(name,'(a,i0)') 'simplification ',k
        call check_identity(suite,eng,trim(name),roots(k)-raw(k))
    end do
    do i=1,6
        call check_identity(suite,eng,'gradient simplification', &
            roots(2+i)-subs_many(first_raw(i),old,replacement))
        k=mod(i+2,6)+1
        tmp=subs_many(first_raw(6+k),old,replacement)
        if(i>3) tmp=-tmp
        call check_identity(suite,eng,'RHS equals J grad K',roots(44+i)-tmp)
    end do
    do j=1,6
        do i=1,j
            call check_identity(suite,eng,'Hessian symmetry', &
                roots(8+(j-1)*6+i)-roots(8+(i-1)*6+j))
        end do
    end do
    spec%name=str('fo_scalar'); spec%mode=KERNEL_SUBROUTINE
    spec%temp_prefix=str('work')
    spec%regenerate_command=str('bash gen/regenerate_fo.sh')
    spec%generator=str('gen/gen_fo_hamiltonian.f90')
    spec%args=args; spec%outputs=outputs
    open(newunit=u,file=trim(outroot)//'/src/generated/fo_scalar.inc',status='replace')
    write(u,'(a)') chars(emit_kernel(roots,spec))
    close(u)
    ! The same Hamiltonian and differentiated expressions, with unused Hessian
    ! outputs omitted. FortSym eliminates their dependencies in this leaf.
    spec%name=str('fo_scalar_first')
    spec%outputs=[outputs(1:8),outputs(45:50)]
    open(newunit=u,file=trim(outroot)//'/src/generated/fo_scalar_first.inc', &
        status='replace')
    first_code=chars(emit_kernel([roots(1:8),roots(45:50)],spec))
    if(first_code(len(first_code):)==new_line('a')) &
        first_code=first_code(:len(first_code)-1)
    write(u,'(a)') first_code
    close(u)
    ! Parse the actual emitted assignments back through FortSym. Expand their
    ! straight-line temporaries in order; no copied expression transcription.
    cse=cse_analyse(roots,'work')
    allocate(temp_symbols(cse%n),temp_values(cse%n))
    do i=1,cse%n
        name=chars(cse%names(i))
        temp_symbols(i)=sym(arena,trim(name))
        tmp=parse_fortran_expr(arena,trim(outroot)//'/src/generated/fo_scalar.inc', &
            trim(name),ok,message)
        if(.not.ok) then
            print *,message
            error stop 'emitted temporary parse refused'
        end if
        if(i>1) tmp=subs_many(tmp,temp_symbols(:i-1),temp_values(:i-1))
        temp_values(i)=tmp
    end do
    do i=1,86
        tmp=parse_fortran_expr(arena,trim(outroot)//'/src/generated/fo_scalar.inc', &
            chars(outputs(i)),ok,message)
        if(.not.ok) then
            print *,message
            error stop 'emitted output parse refused'
        end if
        tmp=subs_many(tmp,temp_symbols,temp_values)
        call check_identity(suite,eng,'emitted '//chars(outputs(i)),tmp-roots(i))
    end do
    call check_first_order_leaf
    ! Array facade is generated from the same ordered interface, never edited.
    open(newunit=u,file=trim(outroot)//'/src/generated/fo_kernel.f90',status='replace')
    write(u,'(a)') 'module fo_kernel', &
        'use, intrinsic :: iso_fortran_env, only: dp => real64', &
        'implicit none','contains'
    do order=1,2
        if(order==1) then
            write(u,'(a)') &
                'subroutine fo_eval(p,aj,gj,mass,charge,e0,tj,energy,kenergy,gh,hh,f,df)'
        else
            write(u,'(a)') &
                'subroutine fo_eval_first(p,aj,gj,mass,charge,e0,tj,energy,kenergy,gh,hh,f,df)'
        end if
        write(u,'(a)') &
            'real(dp),intent(in)::p(3),aj(3,10),gj(6,10),mass,charge,e0,tj(10)', &
            'real(dp),intent(out)::energy,kenergy,gh(6),hh(6,6),f(6),df(6,6)'
        if(order==1) then
            write(u,'(a)') 'call fo_scalar( &'
        else
            write(u,'(a)') 'hh=0; df=0', 'call fo_scalar_first( &'
        end if
        do i=1,3
            write(u,'(a,i0,a)') 'p(',i,'), &'
        end do
        write(u,'(a)') 'mass, charge, e0, &'
        do i=1,3
            do k=1,10
                write(u,'(a,i0,a,i0,a)') 'aj(',i,',',k,'), &'
            end do
        end do
        do i=1,6
            do k=1,10
                write(u,'(a,i0,a,i0,a)') 'gj(',i,',',k,'), &'
            end do
        end do
        do k=1,10
            write(u,'(a,i0,a)') 'tj(',k,'), &'
        end do
        write(u,'(a)') 'energy,kenergy, &'
        do i=1,6
            write(u,'(a,i0,a)') 'gh(',i,'), &'
        end do
        if(order==1) then
        do j=1,6
            do i=1,6
                write(u,'(a,i0,a,i0,a)') 'hh(',i,',',j,'), &'
            end do
        end do
        end if
        do i=1,6
            if(order==2.and.i==6) then
                write(u,'(a)') 'f(6))'
            else
                write(u,'(a,i0,a)') 'f(',i,'), &'
            end if
        end do
        if(order==1) then
        do j=1,6
            do i=1,6
                if(i==6.and.j==6) then
                    write(u,'(a)') 'df(6,6))'
                else
                    write(u,'(a,i0,a,i0,a)') 'df(',i,',',j,'), &'
                end if
            end do
        end do
            write(u,'(a)') 'end subroutine fo_eval'
        else
            write(u,'(a)') 'end subroutine fo_eval_first'
        end if
    end do
    write(u,'(a)') &
        'subroutine fo_radial(gj,rj)', &
        'real(dp),intent(in)::gj(6,10)', 'real(dp),intent(out)::rj(10)', &
        'call fo_radial_scalar( &'
    do i=1,6
        do k=1,10
            write(u,'(a,i0,a,i0,a)') 'gj(',i,',',k,'), &'
        end do
    end do
    do k=1,9
        write(u,'(a,i0,a)') 'rj(',k,'), &'
    end do
    write(u,'(a)') 'rj(10))','end subroutine fo_radial', &
        'subroutine fo_metric_monitor(rj,tj,mref,out)', &
        'real(dp),intent(in)::rj(10),tj(3),mref', &
        'real(dp),intent(out)::out(10)', 'call fo_metric_monitor_scalar( &'
    do k=1,10
        write(u,'(a,i0,a)') 'rj(',k,'), &'
    end do
    write(u,'(a)') 'tj(1),tj(2),tj(3),mref, &'
    do k=1,9
        write(u,'(a,i0,a)') 'out(',k,'), &'
    end do
    write(u,'(a)') 'out(10))','end subroutine fo_metric_monitor', &
        "include 'fo_scalar.inc'", "include 'fo_scalar_first.inc'", &
        "include 'fo_monitor.inc'","include 'fo_product.inc'", &
        "include 'fo_radial.inc'", "include 'fo_metric_monitor.inc'", &
        'end module fo_kernel'
    close(u)
    call suite_end(suite,trim(outroot)//'/concept/fo_symplectic_proof/symbolic.json')
contains
    subroutine check_first_order_leaf
        type(cse_result_t) :: first_cse
        type(expr_t), allocatable :: names(:), values(:)
        type(expr_t) :: leaves(14), parsed
        type(str_t) :: labels(14)
        integer :: j

        leaves = [roots(1:8), roots(45:50)]
        labels = [outputs(1:8), outputs(45:50)]
        first_cse = cse_analyse(leaves, 'work')
        allocate(names(first_cse%n), values(first_cse%n))
        do j = 1, first_cse%n
            names(j) = sym(arena, chars(first_cse%names(j)))
            parsed = parse_fortran_expr(arena, &
                trim(outroot)//'/src/generated/fo_scalar_first.inc', &
                chars(first_cse%names(j)), ok, message)
            if (.not. ok) error stop 'first-order temporary readback refused'
            if (j > 1) parsed = subs_many(parsed, names(:j-1), values(:j-1))
            values(j) = parsed
        end do
        do j = 1, size(leaves)
            parsed = parse_fortran_expr(arena, &
                trim(outroot)//'/src/generated/fo_scalar_first.inc', &
                chars(labels(j)), ok, message)
            if (.not. ok) error stop 'first-order output readback refused'
            parsed = subs_many(parsed, names, values)
            call check_identity(suite, eng, 'first-order '//chars(labels(j)), &
                parsed - leaves(j))
        end do
    end subroutine check_first_order_leaf

    subroutine emit_radial
        type(expr_t)::leaves(10),rfield,tfield,mref,scalar
        type(expr_t)::rold(13),rnew(13)
        type(str_t)::rargs(14)
        integer::k
        leaves=jet(gi(1,1))
        do k=1,10
            leaves(k)=subs_many(leaves(k),old(:90),replacement(:90))
        end do
        call checked_leaf(leaves,args(37:96),'fo_radial_scalar','fo_radial.inc')
        rfield=func('radial',q); tfield=func('rhomonitor',q(1:1))
        rold(1:10)=jet(rfield)
        rold(11)=tfield
        rold(12)=partial_derivative(tfield,1)
        rold(13)=partial_derivative(rold(12),1)
        do k=1,13
            write(name,'(a,i0)') 'rinput',k
            rnew(k)=sym(arena,trim(name)); rargs(k)=str(trim(name))
        end do
        mref=sym(arena,'mref');rargs(14)=str('mref')
        scalar=tfield*sqrt(mref/(mref+rfield))
        leaves=jet(scalar)
        do k=1,10
            leaves(k)=subs_many(leaves(k),rold,rnew)
        end do
        call checked_leaf(leaves,rargs,'fo_metric_monitor_scalar', &
            'fo_metric_monitor.inc')
    end subroutine emit_radial

    function jet(scalar) result(leaves)
        type(expr_t),intent(in)::scalar
        type(expr_t)::leaves(10)
        integer::a1,b1,slot
        leaves(1)=scalar
        do a1=1,3
            leaves(a1+1)=diff(scalar,q(a1))
        end do
        slot=4
        do b1=1,3
            do a1=1,b1
                slot=slot+1
                leaves(slot)=diff(leaves(a1+1),q(b1))
            end do
        end do
    end function jet

    subroutine checked_leaf(leaves,leafargs,leafname,filename)
        type(expr_t),intent(inout)::leaves(:)
        type(str_t),intent(in)::leafargs(:)
        character(*),intent(in)::leafname,filename
        type(kernel_spec_t)::lspec
        type(cse_result_t)::lcse
        type(expr_t),allocatable::lsymbols(:),lvalues(:)
        type(str_t),allocatable::loutputs(:)
        type(expr_t)::parsed
        integer::lu,lk
        allocate(loutputs(size(leaves)))
        do lk=1,size(leaves)
            res=eng%simplify(leaves(lk))
            if(.not.res%ok) error stop 'leaf simplification refused'
            call check_identity(suite,eng,leafname//' simplify',res%value-leaves(lk))
            leaves(lk)=res%value
            write(name,'(a,i0)') 'out',lk
            loutputs(lk)=str(trim(name))
        end do
        lspec%name=str(leafname);lspec%mode=KERNEL_SUBROUTINE
        lspec%args=leafargs;lspec%outputs=loutputs;lspec%temp_prefix=str('work')
        lspec%generator=str('gen/gen_fo_hamiltonian.f90')
        lspec%regenerate_command=str('bash gen/regenerate_fo.sh')
        open(newunit=lu,file=trim(outroot)//'/src/generated/'//filename,status='replace')
        write(lu,'(a)') chars(emit_kernel(leaves,lspec))
        close(lu)
        lcse=cse_analyse(leaves,'work')
        allocate(lsymbols(lcse%n),lvalues(lcse%n))
        do lk=1,lcse%n
            name=chars(lcse%names(lk));lsymbols(lk)=sym(arena,trim(name))
            parsed=parse_fortran_expr(arena,trim(outroot)//'/src/generated/'//filename, &
                trim(name),ok,message)
            if(.not.ok) error stop 'leaf temporary readback refused'
            if(lk>1) parsed=subs_many(parsed,lsymbols(:lk-1),lvalues(:lk-1))
            lvalues(lk)=parsed
        end do
        do lk=1,size(leaves)
            parsed=parse_fortran_expr(arena,trim(outroot)//'/src/generated/'//filename, &
                chars(loutputs(lk)),ok,message)
            if(.not.ok) error stop 'leaf output readback refused'
            parsed=subs_many(parsed,lsymbols,lvalues)
            call check_identity(suite,eng,leafname//' emitted',parsed-leaves(lk))
        end do
    end subroutine checked_leaf

    subroutine emit_monitor
        type(expr_t)::r,le,gm,scalar,mroots(3),rr,aa,bb
        type(expr_t)::jold(6),jnew(6)
        type(kernel_spec_t)::mspec
        type(suite_t)::msuite
        integer::kk
        call suite_begin(msuite,'monitor and jet product identities')
        r=sym(arena,'r'); le=sym(arena,'length'); gm=sym(arena,'gmin')
        scalar=1-(1-gm)*exp(-(r/le)**2)
        mroots(1)=scalar
        mroots(2)=diff(scalar,r)
        mroots(3)=diff(mroots(2),r)
        do kk=1,3
            res=eng%simplify(mroots(kk))
            if(.not.res%ok) error stop 'monitor simplify'
            call check_identity(msuite,eng,'monitor simplify',res%value-mroots(kk))
            mroots(kk)=res%value
        end do
        mspec%name=str('fo_monitor_factor'); mspec%mode=KERNEL_SUBROUTINE
        mspec%temp_prefix=str('work')
        mspec%generator=str('gen/gen_fo_hamiltonian.f90')
        mspec%regenerate_command=str('bash gen/regenerate_fo.sh')
        mspec%args=[str('r'),str('length'),str('gmin')]
        mspec%outputs=[str('t0'),str('t1'),str('t2')]
        call checked_leaf(mroots,mspec%args,'fo_monitor_factor','fo_monitor.inc')
        rr=sym(arena,'r')
        aa=func('leftjet',[rr]); bb=func('rightjet',[rr])
        jold(1)=aa; jold(2)=partial_derivative(aa,1)
        jold(3)=partial_derivative(jold(2),1)
        jold(4)=bb; jold(5)=partial_derivative(bb,1)
        jold(6)=partial_derivative(jold(5),1)
        do kk=1,6
            write(name,'(a,i0)') 'j',kk
            jnew(kk)=sym(arena,trim(name))
        end do
        scalar=aa*bb
        mroots(1)=scalar; mroots(2)=diff(scalar,rr)
        mroots(3)=diff(mroots(2),rr)
        do kk=1,3
            mroots(kk)=subs_many(mroots(kk),jold,jnew)
            res=eng%simplify(mroots(kk))
            if(.not.res%ok) error stop 'product simplify'
            call check_identity(msuite,eng,'product simplify',res%value-mroots(kk))
            mroots(kk)=res%value
        end do
        mspec%name=str('fo_jet_product')
        mspec%args=[str('j1'),str('j2'),str('j3'),str('j4'),str('j5'),str('j6')]
        call checked_leaf(mroots,mspec%args,'fo_jet_product','fo_product.inc')
        call suite_end(msuite,trim(outroot)//'/concept/fo_symplectic_proof/monitor.json')
    end subroutine emit_monitor
    subroutine jet_replacements(field,prefix,component)
        type(expr_t),intent(in)::field
        character(*),intent(in)::prefix
        integer,intent(in)::component
        type(expr_t)::der
        integer::slot,r,t
        character(32)::label
        do slot=1,10
            nrep=nrep+1
            old(nrep)=field
            if(slot>=2.and.slot<=4) &
                old(nrep)=partial_derivative(field,slot-1)
            if(slot>=5) then
                idx=4
                do t=1,3
                    do r=1,t
                        idx=idx+1
                        if(idx==slot) then
                            der=partial_derivative(field,r)
                            old(nrep)=partial_derivative(der,t)
                        end if
                    end do
                end do
            end if
            write(label,'(a,i0,a,i0)') prefix,component,'_',slot
            replacement(nrep)=input(trim(label))
        end do
    end subroutine jet_replacements
    function input(label) result(x)
        character(*),intent(in)::label
        type(expr_t)::x
        n=n+1
        args(n)=str(label)
        x=sym(arena,label)
    end function input
end program gen_fo_hamiltonian
