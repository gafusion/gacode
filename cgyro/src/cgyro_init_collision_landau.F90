module cgyro_init_collision_landau
  implicit none
  logical, parameter :: DIFF_OFF=.false.
  ! ^^ Switch off energy restoration for new Sugama op
  logical, parameter :: KPERP_FIELD_OFF=.false.
  ! ^^ Switch off gyrotrafo of field part, i.e., treat that one drift-kinetic
  integer, parameter :: verbose=1000
  real, parameter :: eps=1e-13
  integer,parameter :: ng2=8 ! number of Gauss points for inner integration
  real,private,parameter :: pi1=atan(1.)*4
  private :: qsort,old_sugama
!!$  interface
!!$     function qsort ( base, nmemb, size,cmp  ) bind ( c,name="qsort")
!!$       use iso_c_binding
!!$       type (*) :: base
!!$       integer :: nmemb,size
!!$       interface
!!$          integer,value :: function cmp(a,b) !<-- this is probably illegal which is why one can't  use this.
!!$            type (*),intent(in) :: a,b
!!$          end function cmp
!!$       end interface
!!$     end function qsort
!!$  end interface
contains
  subroutine cgyro_init_landau()
    ! populate cmat with Galerkin based gyrokinetic Landau operator.
    ! cmat1 is only for comparison purposes
    use cgyro_globals, only : CGYRO_COMM_WORLD,vth,temp,mass,dens,temp_ele,mass_ele,dens_ele,rho,z,&
         n_energy,e_max,n_xi,n_radial,n_theta,n_species,n_toroidal,nt1,nt2,nc_loc_coll,nc_cl1,nc_cl2,nc,nv,&
         nu_ee,&
         xi,w_xi,& !needed for projleg calc
         collision_model,&   ! if this is 7, we switch to calculating Sugama.
         collision_field_max_l,collision_test_max_l,&
         collision_test_mode,&  !if this is 1, we print out comparisons with traditional matrix element
         ic_c,it_c,iv_v,&
         k_perp,bmag,&
         alpha_poly,&
         i_proc,n_proc
    use cgyro_coll_data, only : cmat, copy_into_cmat, copy_from_cmat_all
    use landau, landauvb=>verbose
    use gyrotransformation, gtvb=>verbose
    use half_hermite
    use mpi
    ! need to calculate up to lphys=n_xi-1, n_energy polynomials and for n_species.
    use, intrinsic :: ieee_exceptions
    implicit none
    !character(*),parameter :: sr='init_collision_landau: '
    real, dimension(:), allocatable :: a1,b1,c1,lg,a,bsq,sp,sw,gw,gp,gw2,gp2
    real, dimension(:,:), allocatable :: projsteen,projleg,L2xi,xi2L,poly2v,v2poly,Landau2v,&
         lor,lor1,lor_self,dif,dif1,dif_self,t1t2,t1t21,energymatrix
    real, dimension(:,:,:), allocatable :: field,field2
    real, dimension(:,:,:,:), allocatable :: Landauop,dk_self_field !<- only for KPERP_FIELD_OFF
    real, dimension(:,:,:,:,:,:,:), allocatable :: gyrocolmat

    ! for comparison purposes only for collision_test_mode==1:
    real, allocatable :: cmat1(:,:,:,:)
    real, target, allocatable :: cmat2(:,:,:,:)
    
    ! for Lapack
    real, dimension(:),allocatable :: work
    integer,dimension(:),allocatable :: iwork,isuppz,ifail
    integer lwork,liwork,info

    ! for MPI
    !type(MPI_Status) :: status !mpi_f08
    integer status(MPI_STATUS_SIZE)
    integer ierror

    ! for IEEE
    type(ieee_status_type) ieee_status

    ! for processor allocation
    integer,allocatable:: proc(:,:,:),load(:),gtcost(:,:,:),sortidx(:)
    integer idx,max_load,min_load,cost

    ! for timing
    real t(11),t1,t2

    ! for calculating Sugama op
    real,allocatable :: polyrep(:,:)

    integer ngauss
    integer it,ic,itor,ic_loc,ia,ib,ik,nkmax,is,is1,ns
    integer js
    integer,allocatable :: nk(:,:)
    real, allocatable :: kperp_arr(:,:,:),loss(:),dist(:),id(:,:)
    real,allocatable:: AF(:,:),pv(:,:),em(:,:),Sc(:),ferr(:),berr(:)
    real rcond
    character equed
    real, allocatable, dimension(:,:,:,:) :: m1,m2
    real, dimension(:,:), target, allocatable :: cmat_loc
    real, allocatable :: chebweightarr(:)
    real target_k,target_ik
    integer ix,ie,jx,je,iv,jv
    integer i,j,k,l,m
    real xmax,kperp_bmag_max,rhomax,kperprhomax
    real kperp_bmag,md,d,md1,d1,devi,s,s1
    integer nmaxpoly,lmax,lmax_field,lmax_test
    integer totalcost,lneeded
    real beta,t1t2ratio,normalization,fieldnormalization,fieldnormalization_ba,testnormalization
    logical t1t2flag
    character(1000) :: fn
    real, allocatable, dimension(:,:,:,:,:) :: c
    real, allocatable, dimension(:,:) :: v2polytimesemat,mommat,v2momtimesemat,mom2v
    integer, allocatable, dimension(:) :: nc1_proc,nc2_proc,nt1_proc,nt2_proc
    integer, allocatable, dimension(:,:) :: proc_c
1   format ("init_collision_landau: ",9A)
    if (i_proc==0) print 1,'WARNING: dens_rot not yet implemented!!'
    if (i_proc==0) print 1,'WARNING: nu_global not yet implemented!!'
    if (i_proc==0 .and. collision_test_mode==1) print 1,'collision_test_mode==1, comparing ...'
    if (i_proc==0 .and. maxval(abs(temp(1:n_species)-temp(1)))/=0) then
       print 1,'Warning: Landau not yet working for different species temperatures!!'
    end if


    if (i_proc==0 .and. verbose>0) then ! verboseness settings of modules
      landauvb=1
      gtvb=1
    end if

    !$    call MPI_Barrier(CGYRO_COMM_WORLD,ierror) ! may improve timing
    call cpu_time(t1)
    ns=ispec(n_species,n_species) !number of non-redundant species pairs
    xmax=sqrt(e_max) !cut off at exp(-xmax^2)
    ! find kperp_max in the system
    ! 1st local:
    kperp_bmag_max=0
    if (i_proc==0) then
       do itor=nt1,nt2
          print 9,"itor=",itor
          do i=1,nc
9            format ("init_collision_landau: ",A,I0,5(A,G24.16E3))
             print 9,'kperp(',i,')=',k_perp(i,itor)/bmag(it_c(i))*rho,"bmag(it)=",bmag(it_c(i))
          end do
       enddo
    end if
    do itor=nt1,nt2
       do i=1,n_theta
          do j=1,n_radial
             kperp_bmag=k_perp(ic_c(j,i),itor)/bmag(i)
             !        global array k_perp(1 ... nc,nt1 ... nt2) nc=n_theta*n_radial
             kperp_bmag_max=max(kperp_bmag_max,kperp_bmag)
          end do
       end do
    end do
    ! kperp_bmag_max is not completely global, there is still the n dependence.
    ! we need to maximize over the toroidal mode numbers:
    call MPI_ALLREDUCE(MPI_IN_PLACE,kperp_bmag_max,1,MPI_REAL8,MPI_MAX,CGYRO_COMM_WORLD,ierror)
    rhomax=maxval(abs(rho_spec([(i,i=1,n_species)])))*xmax
    kperprhomax=kperp_bmag_max*rhomax
    if (verbose>0 .and. i_proc==0) print 6,'using kperprhomax=',kperprhomax

    ! first calculate the required polynomials
    nmaxpoly=n_energy+est_extradegree(kperprhomax,eps)
    lmax=n_xi+est_mpullback(kperprhomax,eps=eps)
    lmax=lmax+mod(lmax,2)
    if (verbose>1 .and. i_proc==0) print 7,'using nmaxpoly',nmaxpoly,'lmax=',lmax

    lmax_field=lmax
    if (collision_field_max_l>=-1) lmax_field=min(lmax_field,collision_field_max_l+1)
    lmax_test=lmax
    if (collision_test_max_l>=-1) lmax_test=min(lmax_test,collision_test_max_l+1)
    allocate(a1(nmaxpoly+1),b1(nmaxpoly+1),c1(nmaxpoly+1),a(nmaxpoly+1),bsq(nmaxpoly+1),lg(nmaxpoly+1))
    call half_hermite_norm(nmaxpoly+1,0.,xmax,1.,alpha_poly,a1,b1,c1,a,bsq,lg)
    allocate(sp(nmaxpoly),sw(nmaxpoly),projsteen(nmaxpoly,nmaxpoly))

    !dstevr: liwork>=1,10*nmaxpoly
    liwork=10*nmaxpoly
    !dstevr: lwork >=1,20*nmaxpoly
    lwork=20*nmaxpoly
    !dstevr: isuppz: dim>=2*nmaxpoly

    !dstein: ifail: dim>=nmaxpoly
    !dstein: work: dim>=5*nmaxpoly
    !dstein: iwork: dim>=nmaxpoly
    
    !later need for dposvx below:
    !lowest dimensions: work(3*n_energy), iwork(n)
    
    allocate(iwork(liwork),work(lwork),isuppz(2*nmaxpoly),ifail(nmaxpoly))
    call ieee_set_flag(ieee_all,.false.)
    call ieee_get_status(ieee_status)
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call dstevr('N','A',nmaxpoly,a(1:nmaxpoly),bsq(2:nmaxpoly),0.,0.,0,0,0.,m,sp,&
         projsteen,nmaxpoly,isuppz,work,lwork,iwork,liwork,info)
    if (info/=0 .or. m/=nmaxpoly) then
       if (i_proc==0) print '("init_collision_landau: ",A,9(I3))','dstevr error, here is info,m,nmaxpol',info,m,nmaxpoly
       stop
    end if
    call dstein(nmaxpoly,a(1:nmaxpoly),bsq(2:nmaxpoly),nmaxpoly,sp,(/(1,i=1,nmaxpoly)/),&
         (/(nmaxpoly,i=1,nmaxpoly)/),projsteen,nmaxpoly,work ,iwork,ifail,info)
    call ieee_set_status(ieee_status)
    if (info/=0) then
       if (i_proc==0) print '("init_collision_landau: ",A,2I3," ",A,I3)','dstein error, here is i,info',i,info,'and ifail',ifail
       stop
    end if
    do i=1,nmaxpoly
       if (projsteen(1,i)<0) projsteen(:,i)=-projsteen(:,i)
    end do
    sw=projsteen(1,:)**2/c1(1)**2 ! Steen weights

    if (collision_test_mode==1 .or. collision_model==7) then
       ! Calculate representation of the monomials in terms of the
       ! fine grained half_hermite polynomials.
       ! This is useful for moment calculations, needed for Sugama and diagnostics.
       ! Not more than 3 moments are actually needed at present.
       allocate(polyrep(nmaxpoly,nmaxpoly))
       polyrep=0
       i=1
       polyrep(i,1)=c1(1)
       do i=2,nmaxpoly
          polyrep(i,2:i)=c1(i)*polyrep(i-1,1:i-1)
          polyrep(i,1:i-1)=polyrep(i,1:i-1)-polyrep(i-1,1:i-1)*a1(i)
          if (i>2) polyrep(i,1:i-2)=polyrep(i,1:i-2)-polyrep(i-2,1:i-2)*b1(i)
       enddo
       ! polyrep(i,j) significant for j<=i otherwise zero.
       ! invert Lower triangular matrix, diagonal is Non-unit i.e. significant
       call dtrtri('L','N',nmaxpoly,polyrep,nmaxpoly,info)
       if (info /=0) then
          print 7,'dtrtriinfo',info
          stop
       end if
    end if

    ngauss=max(70,ceiling(15*xmax))+nmaxpoly+ceiling(alpha_poly/2) ! this is a good empirical value for sqrt(2)*xmax
7   format ("init_collision_landau: ",9(A,I0,' '))
    if (i_proc==0 .and. verbose>1) print 7,'using ng2=',ng2,'and ngauss=',ngauss
    allocate(gw(ngauss),gp(ngauss),gw2(ng2),gp2(ng2))
    call gauss_legendre(0.,1.,gp,gw,ngauss)
    call gauss_legendre(0.,1.,gp2,gw2,ng2)

    ! now we go for the real space Landauop collision matrix (Galerkin coefficients)
    allocate(Landauop(nmaxpoly,nmaxpoly,lmax,ns))
    if (KPERP_FIELD_OFF) allocate(dk_self_field(nmaxpoly,nmaxpoly,lmax,n_species))
    allocate(lor_self(nmaxpoly,nmaxpoly),dif_self(nmaxpoly,nmaxpoly),&
         lor(nmaxpoly,nmaxpoly),dif(nmaxpoly,nmaxpoly),&
         lor1(nmaxpoly,nmaxpoly),dif1(nmaxpoly,nmaxpoly),field(nmaxpoly,nmaxpoly,lmax))
    t1t2flag=maxval(abs(temp(1:n_species)-temp(1)))>0 ! need to take different T into account.
    if (t1t2flag) then
       allocate(t1t2(nmaxpoly,nmaxpoly),t1t21(nmaxpoly,nmaxpoly),field2(nmaxpoly,nmaxpoly,lmax))
    end if
    ! first the self-collision matrices
    ! (obviously could save here, since these are symmetric)

    call cpu_time(t2)
    t(6)=t2-t1
    t1=t2
    
    call gentestkernel(nmaxpoly,a1,b1,c1,xmax,1.,gp,gw,ngauss,lor_self,dif_self)
    if (collision_model==7) then
       ! calc. mock up Sugama field op.
       field=0
       ! cancel l=2 (lphys=1) v^1 and l=1 (lphys=0) v^2 polynomial.
       if (.not. DIFF_OFF .and. lmax_field>=1) then
          l=1;k=3
          do i=1,nmaxpoly
             field(:,i,l)=matmul(dif_self,polyrep(k,:))*dot_product(polyrep(k,:),dif_self(:,i))&
                  /dot_product(polyrep(k,:),matmul(dif_self,polyrep(k,:)))
          end do
       end if
       if (lmax_field>=2) then
          l=2;k=2
          do i=1,nmaxpoly
             field(:,i,l)=matmul(dif_self+lor_self*l*(l-1),polyrep(k,:))*&
                  dot_product(polyrep(k,:),dif_self(:,i)+lor_self(:,i)*l*(l-1))&
                  /dot_product(polyrep(k,:),matmul(dif_self+lor_self*l*(l-1),polyrep(k,:)))
          end do
       end if
    else
       if (lmax_field>=1) &
            call genfieldkernel(nmaxpoly,lmax_field,a1,b1,c1,xmax,1.,1.,gp,gw,ngauss,gp2,gw2,ng2,field)
       field(:,:,lmax_field+1:lmax)=0
    end if
    ialoop: do ia=1,n_species
       is=ispec(ia,ia)
       normalization=z(ia)**4*dens(ia)**2/&
            (sqrt(mass(ia))*temp(ia)**1.5) !*normcol <- this we save for later.
       lor=lor_self*normalization
       dif=dif_self*normalization
       if (t1t2flag) t1t2=0
       ibloop: do ib=1,n_species
          if (ib==ia) cycle
          t1t2ratio=temp(ia)/temp(ib)
          beta=sqrt(mass(ib)/mass(ia)*t1t2ratio)
          testnormalization=(z(ia)*z(ib))**2*dens(ia)*dens(ib)*mass(ib)/&
               (mass(ia)**1.5*sqrt(temp(ia))*temp(ib)) !*normcol <- this we save for later.
          ! ^^^ cgyrolandauop.tex Eq. (82) before "Analogously ..."
          if (temp(ia)==temp(ib)) then
             call gentestkernel(nmaxpoly,a1,b1,c1,xmax,beta,gp,gw,ngauss,lor1,dif1)
          else
             call gentestkernel(nmaxpoly,a1,b1,c1,xmax,beta,gp,gw,ngauss,lor1,dif1,t1t2_int=t1t21)
             dif1=dif1+(t1t2ratio-1)*t1t21
          end if
          lor=lor+testnormalization*lor1
          dif=dif+testnormalization*dif1
          if (collision_model==7) then
             ! calc. mock up Sugama field op.
             ! left out for the moment t1t2ratio !!!!
             ! cancel l=2 (lphys=1) v^1 and l=1 (lphys=0) v^2 polynomial.
             ! at this point for every pair a,b with b>a we first have (ia,ib)=(b,a) and then (a,b)
             if (ib>ia) then
                Landauop(:,:,:,ispec(ia,ib))=0
                Landauop(:,:,:,ispec(ib,ia))=0
             end if
             if (.not. DIFF_OFF .and. lmax_field>=1) then
                l=1;k=3
                if (ib>ia) then !first half
                   Landauop(:,1,l,ispec(ia,ib))=-matmul(dif1,polyrep(k,:)-polyrep(1,:)*1.5)/&
                        dot_product(polyrep(k,:)-polyrep(1,:)*1.5,matmul(dif1,polyrep(k,:)-polyrep(1,:)*1.5))&
                        /(dens(ia)*temp(ia))
                   Landauop(1,:,l,ispec(ib,ia))=-matmul(polyrep(k,:)-polyrep(1,:)*1.5,dif1(:,:))&
                        *testnormalization*dens(ia)*temp(ia)
                else
                   do i=nmaxpoly,1,-1
                      Landauop(i,:,l,ispec(ia,ib))=Landauop(1,:,l,ispec(ia&
                           ,ib))*dot_product(dif1(i,:),polyrep(k,:)-polyrep(1&
                           ,:)*1.5)/ dot_product(polyrep(k,:)-polyrep(1,:)&
                           *1.5,matmul(dif1,polyrep(k,:)-polyrep(1,:)*1.5)) &
                           /(dens(ia)*temp(ia))
                      Landauop(:,i,l,ispec(ib,ia))=dot_product(polyrep(k,:)&
                           -polyrep(1,:)*1.5,dif1(:,i))*Landauop(:,1,l&
                           ,ispec(ib,ia)) *testnormalization*dens(ia)*temp(ia)
                   end do
                end if
             end if
             if (.not. DIFF_OFF .and. lmax_field>=2) then
                l=2;k=2
                if (ib>ia) then !first half
                   Landauop(:,1,l,ispec(ia,ib))=-matmul(dif1+lor1*l*(l-1),polyrep(k,:))/&
                        dot_product(polyrep(k,:),matmul(dif1+lor1*l*(l-1),polyrep(k,:)))&
                        /(dens(ia)*sqrt(temp(ia)*mass(ia)))
                   Landauop(1,:,l,ispec(ib,ia))=-matmul(polyrep(k,:),dif1+lor1*l*(l-1))&
                        *testnormalization*dens(ia)*sqrt(temp(ia)*mass(ia))
                else
                   do i=nmaxpoly,1,-1
                      Landauop(i,:,l,ispec(ia,ib))=Landauop(1,:,l,ispec(ia,ib))*&
                           dot_product(dif1(i,:)+lor1(i,:)*l*(l-1),polyrep(k,:))/&
                           dot_product(polyrep(k,:),matmul(dif1+lor1*l*(l-1),polyrep(k,:)))&
                           /(dens(ia)*sqrt(temp(ia)*mass(ia)))
                      Landauop(:,i,l,ispec(ib,ia))=dot_product(polyrep(k,:),dif1(:,i)+lor1(:,i)*l*(l-1))&
                           *Landauop(:,1,l,ispec(ib,ia))&
                           *testnormalization*dens(ia)*sqrt(temp(ia)*mass(ia))
                   end do
                end if
             end if
          end if
       end do ibloop
       Landauop(:,:,lmax_test+1:,is)=0
       do l=1,lmax_test
          Landauop(:,:,l,is)=((l-1)*l)*lor+dif
       end do
       if (KPERP_FIELD_OFF) then
          dk_self_field(:,:,:lmax_field,ia)=-field(:,:,:lmax_field)*normalization
          dk_self_field(:,:,lmax_field+1:,ia)=0
       else
          Landauop(:,:,:lmax_field,is)=Landauop(:,:,:lmax_field,is)-field(:,:,:lmax_field)*normalization
       end if
    end do ialoop
    if (collision_model==7) then
!!$       ! **Now** normalize interspecies Sugama field terms
!!$       do ia=1,n_species
!!$          do ib=1,n_species
!!$             if (ib==ia) cycle
!!$             l=1
!!$             Landauop(:,:,l,ispec(ia,ib))=Landauop(:,:,l,ispec(ia,ib))&
!!$                  /(temp(ia)*dens(ia))*temp(ib)*dens(ib)
!!$             l=2
!!$             Landauop(:,:,l,ispec(ia,ib))=Landauop(:,:,l,ispec(ia,ib))&
!!$                  /(sqrt(temp(ia)*mass(ia))*dens(ia))*sqrt(temp(ib)*mass(ib))*dens(ib)
!!$          end do
!!$       end do
       !Now check momentum and energy conservation (purely for diagnostics)
       if (i_proc==0) then
          allocate(loss(nmaxpoly),dist(nmaxpoly))
          l=1;k=3 !energy
          do ib=1,n_species
             dist=0
             do ia=1,n_species
                if (ib==ia) then
                   loss=matmul(polyrep(k,:),Landauop(:,:,l,ispec(ia,ib)))*temp(ia)*dens(ia)
                else
                   dist=dist+matmul(polyrep(k,:),Landauop(:,:,l,ispec(ia,ib)))*temp(ia)*dens(ia)
                end if
             end do
6            format ("init_collision_landau: ",A,5(G24.16))

             print 9,'spec',ib,' energy:'
             print 6,'loss',loss
             print 6,'dist',dist
          end do
          l=2;k=2 !momentum
          do ib=1,n_species
             dist=0
             do ia=1,n_species
                if (ib==ia) then
                   loss=matmul(polyrep(k,:),Landauop(:,:,l,ispec(ia,ib)))*&
                        sqrt(temp(ia)*mass(ia))*dens(ia)
                else
                   dist=dist+matmul(polyrep(k,:),Landauop(:,:,l,ispec(ia,ib)))*&
                        sqrt(temp(ia)*mass(ia))*dens(ia)
                end if
             end do
             print 9,'spec',ib,' momentum:'
             print 6,'loss',loss
             print 6,'dist',dist
          end do
          deallocate(loss,dist)
       end if
    end if
    ! Now the interspecies collision matrices
    call cpu_time(t2)
    t(5)=t2-t1
    t1=t2
    if (collision_model==6) then
       do ib=1,n_species
          do ia=1,n_species
             if (ia==ib) cycle
             is=ispec(ia,ib)
             is1=ispec(ib,ia)
             t1t2ratio=temp(ia)/temp(ib)
             beta=sqrt(mass(ib)/mass(ia)*t1t2ratio)
             fieldnormalization=(z(ia)*z(ib))**2*dens(ia)*dens(ib)/&
                  (sqrt(mass(ia)*temp(ia))*temp(ib)) !*normcol <- this we save for later.
             ! fieldnormalization for a<->b is this /(beta*t1t2ratio)
             fieldnormalization_ba=fieldnormalization/(beta*t1t2ratio)
             Ta_eq_Tb:          if (temp(ia)==temp(ib)) then
!!$             if (ia>ib .or. (i_proc==0 .and. verbose>4)) then
                if (ia>ib) then
                   if (lmax_field>=1) then
                      call genfieldkernel(nmaxpoly,lmax_field,a1,b1,c1,xmax,beta,1.,&
                           gp,gw,ngauss,gp2,gw2,ng2,Landauop(:,:,:,is))
                      call dscal(nmaxpoly**2*lmax_field,-fieldnormalization,Landauop(:,:,:,is),1)
                   end if
                   Landauop(:,:,lmax_field+1:,is)=0
                end if
!!$             if (ia>ib .and. .not. (i_proc==0 .and. verbose>4)) then
                if (ia>ib) then
                   do l=1,lmax
                      ! in this case Landauop is self-adjoint, i.e. for given l and ia and ib
                      ! Landauop=Landauop^T
#ifdef DIMATCOPY
#ifdef __INTEL_COMPILER
                      !use MKL
                      call mkl_domatcopy('c','t',nmaxpoly,nmaxpoly,1.,Landauop(:,:,l,is),&
                           nmaxpoly,Landauop(:,:,l,is1),nmaxpoly)
#else
                      !use openBLAS
                      call domatcopy('c','t',nmaxpoly,nmaxpoly,1.,Landauop(:,:,l,is),&
                           nmaxpoly,Landauop(:,:,l,is1),nmaxpoly)
#endif
#else
                      !alternative to domatcopy
                      do i=1,nmaxpoly
                         do j=1,nmaxpoly
                            Landauop(j,i,l,is1)=Landauop(i,j,l,is)
                         end do
                      end do
#endif
                   end do
                end if
                !Now check self-adjointness, when possible:
                if (verbose>4 .and. i_proc==0 .and. ia<ib) then
3                  format("init_collision_landau: ",2(A,I3),A,2G24.16)
                   print 3,'Checking self-adjointness at ia',ia,'ib',ib,'Ta=Tb',temp(ia),temp(ib)
                   md=-1
                   do l=1,lmax
                      do i=1,nmaxpoly
                         do j=1,nmaxpoly
                            d=abs(Landauop(i,j,l,is)-Landauop(j,i,l,is1))
                            if (d>1e-15) then
4                              format("init_collision_landau: ",A,5(I0," "),A,2(I0," "),A,3G24.16E3)
                               print 4,'SA-error',i,j,l,ia,ib,'is',is,is1,&
                                    'cm',Landauop(i,j,l,is),Landauop(j,i,l,is1)
                            elseif (d>md) then
                               md=d
                               print 4,'so far max SA deviation',i,j,l,ia,ib,'is',is,is1,&
                                    'cm',Landauop(i,j,l,is),Landauop(j,i,l,is1),d
                            end if
                         end do
                      end do
                   end do
                end if
             else if (ia>ib) then ! Ta /= Tb
                ! in this case genfieldkernel calculates ab and ba field simultaneously
                if (lmax_field>=1) then
                   call genfieldkernel(nmaxpoly,lmax_field,a1,b1,c1,xmax,beta,t1t2ratio,&
                        gp,gw,ngauss,gp2,gw2,ng2,Landauop(:,:,:,is),intkernel2=Landauop(:,:,:,is1))
                end if
                Landauop(:,:,lmax_field+1:,is)=0
                Landauop(:,:,lmax_field+1:,is1)=0
                if (verbose>4 .and. i_proc==0 .and. lmax_field>=1) then
                   print '("init_collision_landau: ",2(A,I0),2(A,G24.16))',&
                        'Checking symmetries of field kernel at ia=',ia,' ib=',ib,' Ta',temp(ia),' Tb',temp(ib)
                   call genfieldkernel(nmaxpoly,lmax_field,a1,b1,c1,xmax,1./beta,1./t1t2ratio,&
                        gp,gw,ngauss,gp2,gw2,ng2,field,intkernel2=field2)
                   md=-1
                   md1=-1
                   do l=1,lmax_field
                      do i=1,nmaxpoly
                         do j=1,nmaxpoly
                            d=abs(Landauop(i,j,l,is)-field2(i,j,l))
                            if (d>1e-15) then
                               print 4,'symm-error',i,j,l,ia,ib,'is',is,is1,&
                                    'cm',Landauop(i,j,l,is),field2(i,j,l),d
                            elseif (d>md) then
                               md=d
                               print 4,'so far max symm error',i,j,l,ia,ib,'is',is,is1,&
                                    'cm',Landauop(i,j,l,is),field2(i,j,l),d
                            end if
                            d1=abs(Landauop(i,j,l,is1)-field(i,j,l))
                            if (d1>1e-15) then
                               print 4,'symm-error (2)',i,j,l,ia,ib,'is',is,is1,&
                                    'cm',Landauop(i,j,l,is1),field(i,j,l),d1
                            elseif (d1>md1) then
                               md1=d1
                               print 4,'so far max symm error (2)',i,j,l,ia,ib,'is',is,is1,&
                                    'cm',Landauop(i,j,l,is1),field(i,j,l),d1
                            end if
                         end do
                      end do
                   end do
                end if

                call dscal(nmaxpoly**2*lmax_field,-fieldnormalization,Landauop(:,:,:,is),1)
                call dscal(nmaxpoly**2*lmax_field,-fieldnormalization_ba,Landauop(:,:,:,is1),1)

             end if Ta_eq_Tb
          end do
       end do
    end if
    ! Final touch to the normalization: (cgyrolandauop.tex eq. 122,124)
    call cpu_time(t2)
    t(4)=t2-t1
    t1=t2
    call dscal(nmaxpoly**2*lmax*ns,normcolcgyro*nu_ee*sqrt(mass_ele)*temp_ele**1.5/dens_ele,Landauop,1)
    if (KPERP_FIELD_OFF) &
         call dscal(nmaxpoly**2*lmax*n_species,normcolcgyro*nu_ee*sqrt(mass_ele)*temp_ele**1.5/dens_ele,dk_self_field,1)
    ! Got the Landauop
    ! The only thing left is its gyro-transformation.
    ! now need to find out how many Chebyshev points for k we need to compute for each combination of species.
    if (i_proc==0) then
       print 6,'landauop(1,1,2,1)',Landauop(1,1,2,1)
    end if
    call cpu_time(t2)
    t(1)=t2-t1
    t1=t2

    allocate(nk(n_species,n_species))
    do ia=1,n_species
       do ib=1,n_species
          nk(ia,ib)=est_k_sampling(kperp_bmag_max*xmax*&
               (abs(rho_spec(ia))+abs(rho_spec(ib))),eps)
          if (i_proc==0 .and. verbose>0) then
             print 7,'number Chebyshev k-points nk(',ia,',',ib,')=',nk(ia,ib)
          end if
       end do
    end do

    ! Now distribute the species pairs, the k, the l, and the n onto the processors
    ! 
    ! Noting that for equal temperatures, we have self-adjointness.
    ! Higher k require more computer time.
    ! In particular: A measure is nmaxpoly**2*mpullback*lmax
    !  which for large krho is propto krho**4.
    !    nmaxpoly=n_energy+est_extradegree(kperprho,eps)
    !        lmax=n_xi+est_mpullback(kperprho,eps=eps)
    !        mpullback=est_mpullback(kperprho,eps=eps)
    ! Then one should not distribute over k but rather over
    ! n and l.
    ! and only as a last resort over k in a sort of round robin fashion.
    ! problem is: for every k we'd need a different set of
    !  - associated Legendre polys
    !  - Steen polys
    !  - would hate this.
    ! So (at first) *only* distribute over k, but using a cost estimate to achieve
    ! balanced distribution. We use the multiplication with the collision matrix
    ! in the inner loop as cost estimate.
    nkmax=maxval(nk)
    allocate(gtcost(nkmax,n_species,n_species),kperp_arr(nkmax,n_species,n_species))
    do ia=1,n_species
       do ib=1,n_species
          do ik=1,nk(ia,ib)
             !chebyshev points to be used for species ia in this pairing (ia,ib):
             kperp_arr(ik,ia,ib)=kperp_bmag_max*rho_spec(ia)*sin(.5*pi1*(ik-.5)/nk(ia,ib))**2
          end do
       end do
    end do
    gtvb=0
    gtcost=0
    do ia=1,n_species
       do ib=1,n_species
          if (ia>=ib .or. temp(ia)/=temp(ib)) then
             do ik=1,nk(ia,ib)
                lneeded=n_xi+est_mpullback(kperp_arr(ik,ia,ib)*xmax,kboundb=kperp_arr(ik,ib,ia)*xmax,eps=eps)
                ! npolyneeded=n_energy+est_extradegree(kbound=max(kperp_arr(ik,ia,ib),kperp_arr(ik,ib,ia))*xmax,eps=eps)
                ! for simplicity, always the maximum available Steen(-like) polys are used.
                gtcost(ik,ia,ib)=(lneeded-n_xi+1)*lneeded !*nmaxpoly**2 !only need relative estimate
                if (gtcost(ik,ia,ib)==0) then
                   print '("init_collision_landau: ",A,I0,A,3I3)','ERROR iproc',i_proc,' cost=0, (ik,ia,ib)=',ik,ia,ib
                   stop
                end if
                if (i_proc==0 .and. verbose>0) then
                   print '("init_collision_landau: ",3(A,I0),A,G24.16E3)',&
                        'cost estimate for ia=',ia,' ib=',ib,' ik=',ik,' is',gtcost(ik,ia,ib)
                end if
                if (i_proc==0 .and. gtcost(ik,ia,ib)==0) then
                   print 7,'WARNING: gtcost for ia=',ia,'ib=',ib,'ik=',ik,'is zero!'
                end if
             end do
          end if
       end do
    end do
    if (i_proc==0 .and. verbose>0) gtvb=1

    ! Now distribute the ia,ib,k *fairly* among the processors, and
    ! when Ta=Tb only do the ia>=ib case.
    ! I.e. minimize max. cost on single processor.
    totalcost=sum(gtcost)
    !sort gtcost:
    allocate(sortidx(n_species**2*nkmax))
    sortidx=[(i,i=1,n_species**2*nkmax)]
    if (i_proc==0) print 7,'qsort:'
    call qsort(sortidx,gtcost,n_species**2*nkmax)
!!$    if (i_proc==0) then
!!$       do i=1,n_species**2*nkmax
!!$          idx=sortidx(i)-1
!!$          ik=mod(idx,nkmax)+1
!!$          idx=idx/nkmax
!!$          ia=mod(idx,n_species)+1
!!$          idx=idx/n_species
!!$          ib=idx+1   
!!$          print 1,'i',i,sortidx(i),'cost',gtcost(ik,ia,ib)
!!$       end do
!!$    end if
    ! Now distribute the work on the processors in a simple heuristic form.
    allocate(load(n_proc),proc(nkmax,n_species,n_species))
    load=0
    proc=0
    allocate(gyrocolmat(n_energy,n_xi,n_energy,n_xi,n_species,n_species,nkmax))
    ! could be optimized
    gyrocolmat=1e300

    min_load=0
    max_load=0
    do i=1,n_species**2*nkmax
       idx=sortidx(i)-1
       ik=mod(idx,nkmax)+1
       idx=idx/nkmax
       ia=mod(idx,n_species)+1
       idx=idx/n_species
       ib=idx+1
       cost=gtcost(ik,ia,ib)
       if (cost==0) exit !done.
       !if (i_proc==0) print *,'init_collision_landau: i=',i,sortidx(i),ik,ia,ib,cost
       do j=1,n_proc
          if (load(j)==min_load .or. load(j)+cost<=max_load) then
             if (i_proc==0) then
                print '("init_collision_landau: ",A,3I3,2(" ",A,I0))','assigning (ik,ia,ib)=',ik,ia,ib,'to proc',j-1,'cost',cost
             end if
             load(j)=load(j)+cost
             proc(ik,ia,ib)=j
             max_load=max(max_load,load(j))
             min_load=minval(load)
             exit
          end if
       end do
    enddo
    call cpu_time(t2)
    t(2)=t2-t1
    call cpu_time(t1)
    do i=1,n_species**2*nkmax
       idx=sortidx(i)-1
       ik=mod(idx,nkmax)+1
       idx=idx/nkmax
       ia=mod(idx,n_species)+1
       idx=idx/n_species
       ib=idx+1
       if (proc(ik,ia,ib)==0) exit
       if (i_proc==proc(ik,ia,ib)-1) then
          ! NOTE: Somehow writing these routines the notion of a and b species flipped.
          ! in gyrotrafo the a species are the source of the perturbation, and the b species
          ! the recipient. But b is on the left of the operator and a on the right.
          ! in the landauoperator here it is oposite.
          if (.not. KPERP_FIELD_OFF .or. ia==ib) then
             call gyrotrafo(gyrocolmat(:,:,:,:,ia,ib,ik),n_energy,n_xi,&
                  Landauop(:,:,:,ispec(ia,ib)),nmaxpoly,lmax,&
                  projsteen,sp,nmaxpoly,xmax,&
                  kperp_arr(ik,ib,ia),kperp_arr(ik,ia,ib),eps=eps)
          end if
          if (KPERP_FIELD_OFF) then
             if (ia==ib) then
                do l=1,n_xi
                   gyrocolmat(:,l,:,l,ia,ia,ik)=gyrocolmat(:,l,:,l,ia,ia,ik)+&
                        dk_self_field(:n_energy,:n_energy,l,ia)
                end do
             else
                gyrocolmat(:,:,:,:,ia,ib,ik)=0
                do l=1,n_xi
                   gyrocolmat(:,l,:,l,ia,ib,ik)=Landauop(:n_energy,:n_energy,l,ispec(ia,ib))
                end do
             end if
          end if
       end if
    enddo
    call cpu_time(t2)
    t(3)=t2-t1
    t1=t2
    ! Now do the fill in of potentially left out symmetric halfs
    do i=1,n_species**2*nkmax
       idx=sortidx(i)-1
       ik=mod(idx,nkmax)+1
       idx=idx/nkmax
       ia=mod(idx,n_species)+1
       idx=idx/n_species
       ib=idx+1
       if (proc(ik,ia,ib)==0) exit
       if (i_proc==proc(ik,ia,ib)-1 .and. ia>ib .and. temp(ia)==temp(ib)) then
#ifdef DIMATCOPY
#ifdef __INTEL_COMPILER
          !use MKL
          call mkl_domatcopy('c','t',n_xi*n_energy,n_xi*n_energy,1.,&
               gyrocolmat(:,:,:,:,ia,ib,ik),n_xi*n_energy,&
               gyrocolmat(:,:,:,:,ib,ia,ik),n_xi*n_energy)
#else
          !use openBLAS
          call domatcopy('c','t',n_xi*n_energy,n_xi*n_energy,1.,&
               gyrocolmat(:,:,:,:,ia,ib,ik),n_xi*n_energy,&
               gyrocolmat(:,:,:,:,ib,ia,ik),n_xi*n_energy)
#endif
#else
          !alternative to domatcopy
          do k=1,n_xi
             do l=1,n_energy
                do j=1,n_xi
                   do m=1,n_energy
                      gyrocolmat(m,j,l,k,ib,ia,ik)=gyrocolmat(l,k,m,j,ia,ib,ik)
                   end do
                end do
             end do
          end do
#endif
       end if
    end do
    call cpu_time(t2)
    t(9)=t2-t1
    t1=t2
    ! next: transform from polynomial into vertex spaced for angular and
    ! radial direction
    allocate(projleg(n_xi,n_xi),L2xi(n_xi,n_xi),xi2L(n_xi,n_xi))
!!$    if (i_proc==0) gtvb=5
    call calc_projleg(projleg,n_xi,n_xi,xi,2*w_xi) !interval=[-1,1] so w_xi*=2
!!$    if (i_proc==0) gtvb=1
    do i=1,n_xi
       L2xi(i,:)=projleg(i,:)/sqrt(w_xi(i))
       xi2L(:,i)=projleg(i,:)*sqrt(w_xi(i))
    end do
    if (verbose>2 .and. i_proc==0) then
       allocate(id(n_xi,n_xi))
       call dgemm('N','N',n_xi,n_xi,n_xi,1.,xi2L,n_xi,L2xi,n_xi,0.&
            ,id,n_xi)
       devi=0
       do i=1,n_xi
          do j=1,n_xi
             if (i==j) then
                d=abs(id(i,j)-1)
             else
                d=abs(id(i,j))
             end if
             if (d>devi) then
                print '("init_collision_landau: ",A,G24.16E3,A,2I4,A,G23.16E3)',&
                     'so far max xi-mat deviation',d,' at i,j',i,j,' mat=',id(i,j)
                devi=d
             end if
          end do
       end do
       deallocate(id)
    end if
    ! do the same thing from above again for the lesser n_energy vertices:
    deallocate(projsteen)

    allocate(projsteen(n_energy,n_energy)) 
    call ieee_get_status(ieee_status)
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call dstevr('N','A',n_energy,a(1:n_energy),bsq(2:n_energy),0.,0.,0,0,0.,m,sp,&
         projsteen,n_energy,isuppz,work,lwork,iwork,liwork,info)
    if (info/=0 .or. m/=n_energy) then
       if (i_proc==0) print 4,'dstevr error, here is info,m,nmaxpol',info,m,n_energy
       stop
    end if
    call dstein(n_energy,a(1:n_energy),bsq(2:n_energy),n_energy,sp,(/(1,i=1,n_energy)/),&
         (/(n_energy,i=1,n_energy)/),projsteen,n_energy,work ,iwork,ifail,info)
    if (info/=0) then
       if (i_proc==0) print '(A,2I0," ",A,I0)','dstein error, here is i,info',i,info,'and ifail',ifail
       stop
    end if
    call ieee_set_status(ieee_status)
    do i=1,n_energy
       if (projsteen(1,i)<0) projsteen(:,i)=-projsteen(:,i)
    end do
    !sw=projsteen(1,:)**2/c1(1)**2 ! Steen weights
    ! now polycoeff=projsteen*vertexcoeff*sqrt(weights)
    !   vertexcoeff*sqrt(weights)=projsteen^T*polycoeff.
    allocate(v2poly(n_energy,n_energy),poly2v(n_energy,n_energy))
    do i=1,n_energy
       v2poly(:,i)=projsteen(:,i)*(projsteen(1,i)/c1(1))
       poly2v(:,i)=projsteen(:,i)/(projsteen(1,i)/c1(1))
       ! ^^ purposely calculating the "transposed" matrix
    end do
    ! now polycoeff=v2poly*vertexcoeff
    !   vertexcoeff=poly2v^T*polycoeff.
    if (verbose>2 .and. i_proc==0) then
       allocate(id(n_energy,n_energy))
       call dgemm('T','N',n_energy,n_energy,n_energy,1.,poly2v,n_energy,v2poly,n_energy,0.&
            ,id,n_energy)
       devi=0
       do i=1,n_energy
          do j=1,n_energy
             if (i==j) then
                d=abs(id(i,j)-1)
             else
                d=abs(id(i,j))
             end if
             if (d>devi) then
                print '("init_collision_landau: ",A,G24.16,A,2I4,A,G23.16E3)',&
                     'so far max "e"-mat deviation',d,' at i,j',i,j,' mat=',id(i,j)
                devi=d
             end if
          end do
       end do
       deallocate(id)
    end if
    call ieee_set_flag(ieee_all,.false.)

    allocate(energymatrix(n_energy,n_energy))
    call genenergymatrix(n_energy,a1,b1,c1,xmax,gp,gw,ngauss,energymatrix)
    !   polycoeff=emat^-1*Landauop. ; vertexcoeff=poly2v^T*emat^-1*Landauop
    !   ==> Landau2v=poly2v^T*emat^-1=(emat^-1*poly2v)^T
    allocate(Landau2v(n_energy,n_energy))
    allocate(AF(n_energy,n_energy),Sc(n_energy),ferr(n_energy),&
         berr(n_energy))
    ! lowest dimensions: work(3*n_energy), iwork(n_energy)
    ! solve for vertex projection matrix
    allocate(pv(n_energy,n_energy),em(n_energy,n_energy))
    em=energymatrix
    pv=poly2v
    call dposvx('E','U', n_energy,n_energy,em,n_energy, AF,n_energy,&
         EQUED, Sc, pv,n_energy,Landau2v,n_energy,&
         RCOND, FERR, BERR, WORK, IWORK, INFO)
    ! note: energymatrix, poly2v are potentially destroyed by this call.
    if (verbose>1 .and. i_proc==0) print 6,'energy matrix dposvx rcond=',rcond
    if (info/=0) then
       print 7,'dposvx info=',info,'i_proc',i_proc
       stop
    end if
    deallocate(AF,Sc,ferr,berr,pv,em)
    ! *** Note that Landau2v needs to be transposed upon application.
    allocate(m1(n_energy,n_xi,n_energy,n_xi),m2(n_energy,n_xi,n_energy,n_xi))
    do i=1,n_species**2*nkmax
       idx=sortidx(i)-1
       ik=mod(idx,nkmax)+1
       idx=idx/nkmax
       ia=mod(idx,n_species)+1
       idx=idx/n_species
       ib=idx+1
       if (proc(ik,ia,ib)==0) exit
       if (i_proc==proc(ik,ia,ib)-1) then
          call dgemm('T','N',n_energy,n_xi**2*n_energy,n_energy,1.,Landau2v&
               ,n_energy,gyrocolmat(:,:,:,:,ia,ib,ik),n_energy,0.,&
               m1,n_energy)
          do k=1,n_xi
             do j=1,n_energy
                call dgemm('N','T',n_energy,n_xi,n_xi,1.,m1(:,:,j,k),n_energy,&
                     L2xi,n_xi,0.,m2(:,:,j,k),n_energy)
             end do
          end do
          do j=1,n_xi
             call dgemm('N','N',n_xi*n_energy,n_energy,n_energy,1.,m2(:,:,:&
                  ,j),n_xi*n_energy,v2poly,n_energy,0.,&
                  m1(:,:,:,j),n_xi*n_energy)
          end do
          call  dgemm('N','N',n_xi*n_energy**2,n_xi,n_xi,1.,m1,n_xi*n_energy&
               **2,xi2L,n_xi,0.,gyrocolmat(:,:,:,:,ia,ib,ik),n_xi*n_energy**2)
          if (ia>ib .and. temp(ia)==temp(ib)) then
             call dgemm('T','N',n_energy,n_xi**2*n_energy,n_energy,1.,Landau2v&
                  ,n_energy,gyrocolmat(:,:,:,:,ib,ia,ik),n_energy,0.,&
                  m1,n_energy)
             do k=1,n_xi
                do j=1,n_energy
                   call dgemm('N','T',n_energy,n_xi,n_xi,1.,m1(:,:,j,k),n_energy,&
                        L2xi,n_xi,0.,m2(:,:,j,k),n_energy)
                end do
             end do
             do j=1,n_xi
                call dgemm('N','N',n_xi*n_energy,n_energy,n_energy,1.,m2(:,:,:&
                     ,j),n_xi*n_energy,v2poly,n_energy,0.,&
                     m1(:,:,:,j),n_xi*n_energy)
             end do
             call  dgemm('N','N',n_xi*n_energy**2,n_xi,n_xi,1.,m1,n_xi*n_energy&
                  **2,xi2L,n_xi,0.,gyrocolmat(:,:,:,:,ib,ia,ik),n_xi*n_energy**2)
          end if
       end if
    enddo
    !deallocate(m1)
    deallocate(m2)
    call cpu_time(t2)
    t(7)=t2-t1
    if (i_proc==0) then
       print 7,'pre_scatter timing:'
       do i=1,n_proc
          if (i>1) then
             call MPI_Recv(t,11,MPI_REAL8,i-1,i-1,CGYRO_COMM_WORLD,status,ierror)
          end if
5         format("init_collision_landau: ",A,I0,A,7G24.16E3,A,I0,A,G24.16E3)
          print 5,'i_proc=',i-1,' took',t(1:7),' load ',load(i),' rel',t(3)/load(i)
       end do
    else
       call MPI_Send(t,11,MPI_REAL8,0,i_proc,CGYRO_COMM_WORLD,ierror)
    end if
    call cpu_time(t1)
    ! Now do the scatter

    do i=1,n_species**2*nkmax
       idx=sortidx(i)-1
       ik=mod(idx,nkmax)+1
       idx=idx/nkmax
       ia=mod(idx,n_species)+1
       idx=idx/n_species
       ib=idx+1
       if (proc(ik,ia,ib)/=0) then
!!$          do j=1,n_proc
!!$             call MPI_BARRIER(CGYRO_COMM_WORLD,ierror)
!!$!          if (i_proc==0 .and. verbose>100) then
!!$             if (i_proc==j-1) print *,'bcasting (ik,ia,ib,proc)=',ik,ia,ib,proc(ik,ia,ib)-1,'ip',i_proc
!!$          !          end if
!!$             call MPI_BARRIER(CGYRO_COMM_WORLD,ierror)
!!$          end do
          call MPI_Bcast(gyrocolmat(:,:,:,:,ia,ib,ik),n_xi**2*n_energy**2,&
               MPI_REAL8,proc(ik,ia,ib)-1,CGYRO_COMM_WORLD,ierror)
          if (ia>ib .and. temp(ia)==temp(ib)) then
             call MPI_Bcast(gyrocolmat(:,:,:,:,ib,ia,ik),n_xi**2*n_energy**2,&
                  MPI_REAL8,proc(ik,ia,ib)-1,CGYRO_COMM_WORLD,ierror)
          end if
       end if
    enddo
    call cpu_time(t2)
    t(8)=t2-t1
    t1=t2

    allocate(cmat_loc(nv,nv))
    ! now to the k-interpolation of cmat
    allocate(chebweightarr(nkmax))
    !allocate(m1(n_energy,n_xi,n_energy,n_xi))
    do itor=nt1,nt2
       do ic_loc=1,nc_loc_coll
          ic=ic_loc-1+nc_cl1
          it=it_c(ic)
          target_k=k_perp(ic,itor)/bmag(it)/kperp_bmag_max
          !target_k=sin(.5*pi1*(target_ik-.5)/nk(ia,ib))**2
          do ia=1,n_species
             do ib=1,n_species
                target_ik=asin(sqrt(target_k))*(nk(ia,ib)/(.5*pi1))+.5
                ! for sinc --> see below
                chebweightarr(1:nk(ia,ib))=[(sinc(target_ik-i,nk(ia,ib))+sinc(target_ik+i-1,nk(ia,ib)),&
                     i=1,nk(ia,ib))]
                if (verbose>4 .and. i_proc==0) then
                   print *,'interpolation for krel=',target_k,'target_ik=',target_ik,&
                        'weightsum',sum(chebweightarr(1:nk(ia,ib))),&
                        'ipoltest',sum(chebweightarr(1:nk(ia,ib))*kperp_arr(1:nk(ia,ib),ia,ib)),&
                        'should be',target_k*rho_spec(ia)*kperp_bmag_max
                   print *,'interpolationtest',sum(chebweightarr(1:nk(ia,ib))*&
                        [(sin(.5*pi1*(i-.5)/nk(ia,ib))**2,i=1,nk(ia,ib))]),&
                        'should be',sin(.5*pi1*(target_ik-.5)/nk(ia,ib))**2,nk(ia,ib)
                end if
                call dgemm('N','N',n_xi**2*n_energy**2,1,nk(ia,ib),1.,&
                     gyrocolmat(:,:,:,:,ia,ib,1),n_xi**2*n_energy**2*n_species**2,chebweightarr,nk(ia,ib),&
                     0.,m1,n_xi**2*n_energy**2)
                do jx=1,n_xi
                   do je=1,n_energy
                      jv=iv_v(je,jx,ib)
                      do ix=1,n_xi
                         do ie=1,n_energy
                            iv=iv_v(ie,ix,ia)
                            cmat_loc(iv,jv)=-m1(ie,ix,je,jx)
!!$                           if (gyrocolmat(ie,ix,je,jx,ia,ib,1)==0 .and. ic_loc==1) &
!!$                                print *,'zero gc i_proc',i_proc,&
!!$                                '(ie,ix,ia,je,jx,ib,ic_loc)',ie,ix,ia,je,jx,ib,ic_loc
                         end do
                      end do
                   end do
                end do
                call copy_into_cmat(cmat_loc,ic_loc,itor)
             end do
          end do
       end do
    end do
    deallocate(chebweightarr)
    deallocate(cmat_loc)
    deallocate(m1)
    call cpu_time(t2)
    t(10)=t2-t1
    t1=t2

    if (i_proc==0) then
       do i=1,n_proc
          if (i>1) then
             call MPI_Recv(t,11,MPI_REAL8,i-1,i-1,CGYRO_COMM_WORLD,status,ierror)
          end if
          print *,'i_proc=',i-1,'took',t(1:10),'load',load(i),'rel',t(3)/load(i)
       end do
    else
       call MPI_Send(t,11,MPI_REAL8,0,i_proc,CGYRO_COMM_WORLD,ierror)
    end if
    
    coltestmode: if(collision_test_mode==1) then
       allocate(cmat1(nv,nv,nc_loc_coll,nt1:nt2))
       allocate(cmat2(nv,nv,nc_loc_coll,nt1:nt2))
       ! compute old collision operator in cmat1
       call old_sugama(cmat1) !<--- this is in this module.
       ! Make a copy of the main cmat into cmat2, to ensure it is CPU-accessible
       call copy_from_cmat_all(cmat2)
!!$    do ia=1,n_species
!!$       do ib=1,n_species
!!$          if (ia>=ib .and. temp(ia)/=temp(ib)) then
!!$             ! need to do it, otherwise use self-adjointness
!!$             is=ispec(ia,ib)
!!$             do ik=1,nk(is)
!!$                do l=1,n_xi
!!$                   do i=1,n_energy !polynomial index

       !compare with supplied cmat1
       allocate(m1(n_energy,n_xi,n_energy,n_xi),m2(n_energy,n_xi,n_energy,n_xi),&
            c(n_energy,n_xi,n_energy,n_xi,2))
       allocate(v2polytimesemat(n_energy,n_energy))
       ! I do not know how everything is distributed. I only know the collision matrix 
       allocate (nc1_proc(n_proc),nc2_proc(n_proc),nt1_proc(n_proc),nt2_proc(n_proc),proc_c(nc,0:n_toroidal))
       nc1_proc(i_proc+1)=nc_cl1
       nc2_proc(i_proc+1)=nc_cl2
       nt1_proc(i_proc+1)=nt1
       nt2_proc(i_proc+1)=nt2
       proc_c=0 ! dummy value if no processor is responsible
       do i=1,n_proc
          call MPI_BCAST(nc1_proc(i),1,MPI_INTEGER,i-1,CGYRO_COMM_WORLD,ierror)
          call MPI_BCAST(nc2_proc(i),1,MPI_INTEGER,i-1,CGYRO_COMM_WORLD,ierror)
          call MPI_BCAST(nt1_proc(i),1,MPI_INTEGER,i-1,CGYRO_COMM_WORLD,ierror)
          call MPI_BCAST(nt2_proc(i),1,MPI_INTEGER,i-1,CGYRO_COMM_WORLD,ierror)
          ! this assigns the processor with the highest number to the respective
          ! nc and nt range
          proc_c(nc1_proc(i):nc2_proc(i),nt1_proc(i):nt2_proc(i))=i
       end do
       if (i_proc==0) then
          print *,'nc1',nc1_proc
          print *,'nc2',nc2_proc
          print *,'nt1',nt1_proc
          print *,'nt2',nt2_proc
          print *,'proc_c',proc_c
       end if
       call dgemm('n','n',n_energy,n_energy,n_energy,1.,energymatrix,n_energy,v2poly,n_energy,&
            0.,v2polytimesemat,n_energy)
       md=-1
       do itor=1,n_toroidal ! don't know what toroidal indices actually exist. We don't try 0 here.
          do ic=1,nc,nc-1
             ! If proc_c==0 then there is nothing to be done.
             ! Otherwise we deal with it, if we are either the responsible processor or proc 0
             if (.not. (proc_c(ic,itor)==i_proc+1 .or. i_proc==0 .and. proc_c(ic,itor)/=0)) continue
             ic_loc=ic+1-nc_cl1
             do ia=1,n_species
                specbloop: do ib=1,n_species
                   if (proc_c(ic,itor)==i_proc+1) then
                      do jx=1,n_xi
                         do je=1,n_energy
                            jv=iv_v(je,jx,ib)
                            do ix=1,n_xi
                               do ie=1,n_energy
                                  iv=iv_v(ie,ix,ia)
                                  c(ie,ix,je,jx,1)=cmat2(iv,jv,ic_loc,itor)
                                  c(ie,ix,je,jx,2)=cmat1(iv,jv,ic_loc,itor)
                               end do
                            end do
                         end do
                      end do
                      do l=1,2
                         call dgemm('N','N',n_energy,n_xi**2*n_energy,n_energy,1.,v2polytimesemat&
                              ,n_energy,c(:,:,:,:,l),n_energy,0.,&
                              m1,n_energy)
                         do k=1,n_xi
                            do j=1,n_energy
                               call dgemm('N','T',n_energy,n_xi,n_xi,1.,m1(:,:,j,k),n_energy,&
                                    xi2L,n_xi,0.,m2(:,:,j,k),n_energy)
                            end do
                         end do
                         do j=1,n_xi
                            call dgemm('N','T',n_xi*n_energy,n_energy,n_energy,1.,m2(:,:,:&
                                 ,j),n_xi*n_energy,poly2v,n_energy,0.,&
                                 m1(:,:,:,j),n_xi*n_energy)
                         end do
                         call  dgemm('N','N',n_xi*n_energy**2,n_xi,n_xi,1.,m1,n_xi*n_energy&
                              **2,L2xi,n_xi,0.,c(:,:,:,:,l),n_xi*n_energy**2)
                      end do
                      if (i_proc/=0) then
                         call MPI_SEND(c,size(c),MPI_REAL8,0,1234,CGYRO_COMM_WORLD,ierror)
                      end if
                   else
                      if (i_proc==0) then
                         call MPI_RECV(c,size(c),MPI_REAL8,proc_c(ic,itor)-1,1234,CGYRO_COMM_WORLD,status,ierror)
                      end if
                   end if
!!$                  block
!!$                    real s(n_energy,n_xi),s1(n_energy,n_xi)
!!$                    if (i_proc==0) then
!!$                       s=0
!!$                       s1=0
!!$                       do ix=1,n_xi
!!$                          do ie=1,n_energy
!!$                             do jx=1,n_xi
!!$                                do je=1,n_energy
!!$                                   s(je,jx)=s(je,jx)+cmat2(iv_v(je,jx,ia),iv_v(ie,ix,ib),ic_loc,itor)
!!$                                   s1(je,jx)=s1(je,jx)+cmat1(iv_v(je,jx,ia),iv_v(ie,ix,ib),ic_loc,itor)
!!$                                end do
!!$                             end do
!!$                          end do
!!$                       end do
!!$
!!$                       print *,'particle conservation',ia,ib,'s',s
!!$                       print *,'particle conservation',ia,ib,'s1',s1
!!$                    end if
!!$                  end block
                   if (i_proc==0) then
                      jxloop: do jx=1,n_xi !n_xi
                         do je=1,n_energy/2
                            do ix=1,n_xi
                               do ie=1,n_energy/2 !n_energy
                                  s=c(ie,ix,je,jx,1)
                                  s1=c(ie,ix,je,jx,2)
                                  d=abs(s-s1)
                                  if (d>md) then
                                     md=d
                                  end if
                                  print '(A,4I3,G14.6,2(A,I3,I4),A,G25.16E3,A,2G25.16E3)','ia,ib,itor,ic_loc,kp',ia,&
                                       ib,itor,ic_loc,k_perp(ic,itor)/bmag(it_c(ic))*rho*sqrt(2.),&
                                       'jx,je',jx,je,' ix,ie',ix,ie,' d',d,' c,c1',s,s1
                               end do
                            end do
                         end do
                      end do jxloop
                   end if
                end do specbloop
             end do
          end do
       end do
       ! print out all kradial for certain moments:
       allocate(mommat(3,n_energy))
       mommat(1,:)=polyrep(1,:n_energy)
       mommat(2,:)=polyrep(2,:n_energy)
       mommat(3,:)=polyrep(3,:n_energy)-1.5*polyrep(1,:n_energy) !no units attached
       allocate(v2momtimesemat(3,n_energy))
       call dgemm('n','n',3,n_energy,n_energy,1.,mommat,3,v2polytimesemat,n_energy,0.,v2momtimesemat,3)
       allocate(mom2v(3,n_energy))
       call dgemm('n','n',3,n_energy,n_energy,1.,mommat,3,poly2v,n_energy,0.,mom2v,3)
       if (i_proc==0) then
          do l=1,2
             do ia=1,n_species
                do ib=1,n_species
                   if (collision_model==6) then
                      ! V3=Landau V2=old sugama V1=new sugama
                      write(fn,"('cgyro.moments.V',I1,'a',I1,'b',I1)") 4-l,ia,ib
                   else
                      write(fn,"('cgyro.moments.V',I1,'a',I1,'b',I1)") l,ia,ib
                   end if
                   open(6+ib+n_species*(ia+n_species*l),file=fn)
                end do
             enddo
          end do
       end if
       do itor=1,n_toroidal ! don't know what toroidal indices actually exist. We don't try 0 here.
          do ic=1,nc
             ! If proc_c==0 then there is nothing to be done.
             ! Otherwise we deal with it, if we are either the responsible processor or proc 0
             if (.not. (proc_c(ic,itor)==i_proc+1 .or. i_proc==0 .and. proc_c(ic,itor)/=0)) continue
             ic_loc=ic+1-nc_cl1
             do ia=1,n_species
                specbloop2: do ib=1,n_species
                   if (ic>=nc_cl1 .and. ic<=nc_cl2) then
                      do jx=1,n_xi
                         do je=1,n_energy
                            jv=iv_v(je,jx,ib)
                            do ix=1,n_xi
                               do ie=1,n_energy
                                  iv=iv_v(ie,ix,ia)
                                  c(ie,ix,je,jx,1)=cmat2(iv,jv,ic_loc,itor)
                                  c(ie,ix,je,jx,2)=cmat1(iv,jv,ic_loc,itor)
                               end do
                            end do
                         end do
                      end do
!!$                     do l=1,2
!!$                        call dgemm('N','N',n_energy,n_xi**2*n_energy,n_energy,1.,v2polytimesemat&
!!$                             ,n_energy,c(:,:,:,:,l),n_energy,0.,&
!!$                             m1,n_energy)
!!$                        do k=1,n_xi
!!$                           do j=1,n_energy
!!$                              call dgemm('N','T',n_energy,n_xi,n_xi,1.,m1(:,:,j,k),n_energy,&
!!$                                   xi2L,n_xi,0.,m2(:,:,j,k),n_energy)
!!$                           end do
!!$                        end do
!!$                        do j=1,n_xi
!!$                           call dgemm('N','T',n_xi*n_energy,n_energy,n_energy,1.,m2(:,:,:&
!!$                                ,j),n_xi*n_energy,poly2v,n_energy,0.,&
!!$                                m1(:,:,:,j),n_xi*n_energy)
!!$                        end do
!!$                        call  dgemm('N','N',n_xi*n_energy**2,n_xi,n_xi,1.,m1,n_xi*n_energy&
!!$                             **2,L2xi,n_xi,0.,c(:,:,:,:,l),n_xi*n_energy**2)
!!$                     end do
                      
                      do l=1,2
                         call dgemm('N','N',3,n_xi**2*n_energy,n_energy,1.,v2momtimesemat&
                              ,3,c(:,:,:,:,l),n_energy,0.,&
                              m1,n_energy)
                         do k=1,n_xi
                            do j=1,n_energy
                               call dgemm('N','T',n_energy,n_xi,n_xi,1.,m1(:,:,j,k),n_energy,&
                                    xi2L,n_xi,0.,m2(:,:,j,k),n_energy)
                            end do
                         end do
                         do j=1,n_xi
                            call dgemm('N','T',n_xi*n_energy,3,n_energy,1.,m2(:,:,:&
                                 ,j),n_xi*n_energy,mom2v,3,0.,&
                                 m1(:,:,:,j),n_xi*n_energy)
                         end do
                         call  dgemm('N','N',n_xi*n_energy*3,n_xi,n_xi,1.,m1,n_xi*n_energy&
                              **2,L2xi,n_xi,0.,c(:,:,:,:,l),n_xi*n_energy**2)
                      end do
                      if (i_proc/=0) then
                         call MPI_SEND(c,size(c),MPI_REAL8,0,1234,CGYRO_COMM_WORLD,ierror)
                      end if
                   else
                      if (i_proc==0) then
                         call MPI_RECV(c,size(c),MPI_REAL8,proc_c(ic,itor)-1,1234,CGYRO_COMM_WORLD,status,ierror)
                      end if
                   end if

                   if (i_proc==0) then
                      do l=1,2
                         ! kperp n-diff,mom-diff heat-diff T->n n->T
                         write(6+ib+n_species*(ia+n_species*l),'(6(G24.16))') &
                              k_perp(ic,itor)/bmag(it_c(ic))*rho,c(1,1,1,1,l),c(2,2,2,2,l),&
                              c(3,1,3,1,l),c(1,1,3,1,l),c(3,1,1,1,l)
                      end do
                   end if
                end do specbloop2
             end do
          end do
       end do
       if (i_proc==0) then
          do i=1,2*n_species**2
             close(6+i)
          end do
       end if
       deallocate(m1,m2,c,v2polytimesemat,nc1_proc,nc2_proc,nt1_proc,nt2_proc,proc_c)

       ! now let's compare on processor 0
       md=-1
       if (i_proc==0) then
          do itor=nt1,nt2
             do ic_loc=1,nc_loc_coll
                do is=1,n_species
                   do ix=1,n_xi
                      do ie=1,n_energy
                         iv=iv_v(ie,ix,is)
                         do js=1,n_species
                            do jx=1,n_xi
                               do je=1,n_energy
                                  jv=iv_v(je,jx,js)
                                  d=abs(cmat2(iv,jv,ic_loc,itor)-cmat1(iv,jv,ic_loc,itor))
                                  if (d>md) then
                                     md=d
77                                   format (A,2I3,2(A,3I3),3(A,G23.16))
                                     !if (i_proc==0) &
                                     print 77,'so far max cmat diff @ (itor,ic_loc)',itor,ic_loc,'(is,ix,ie)=',is,ix,ie,&
                                          '(js,jx,je)=',js,jx,je,'d=',d,'c=',cmat2(iv,jv,ic_loc,itor),&
                                          'c1=',cmat1(iv,jv,ic_loc,itor)
                                  end if
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          enddo
       end if
       call MPI_reduce(md,d,1,MPI_REAL8,MPI_MAX,0,CGYRO_COMM_WORLD,ierror)
       if (i_proc==0) print 11,'Max. deviation over all processors:',d
11      format ('cgyro_in._col.: ',A,G23.16)

       call MPI_Barrier(CGYRO_COMM_WORLD,ierror)
       call MPI_finalize(ierror)
       stop
    end if coltestmode

!!$    call MPI_Barrier(CGYRO_COMM_WORLD,ierror)
!!$    print *,'i_proc',i_proc,'done with init_landau'
!!$    call MPI_Barrier(CGYRO_COMM_WORLD,ierror)

  contains
    elemental real function sinc(target_k,halfperiod)
      ! sinc function for chebyshev interpolation
      real target_k
      integer halfperiod
      intent(in):: target_k,halfperiod
      if (target_k==0.) then
         sinc=1.
      else
         sinc=sin(pi1*target_k)/tan((.5*pi1/halfperiod)*target_k)/(2*halfperiod)
      end if
    end function sinc
    elemental real function rho_spec(ia)
      implicit none
      integer,intent(in) :: ia
      rho_spec=vth(ia)*mass(ia)/z(ia)*sqrt(2.)*rho
    end function rho_spec
    integer function ispec(ia,ib)
      implicit none
      integer ia,ib
      ! unfortunately, for the case of unequal temperatures we do not have a
      ! self-adjoint (symmetric) matrix
      ! maybe we will switch at a later time to a self-adjoint formulation
      ! even for unequal temperatures.
      ispec=ia+(ib-1)*n_species
!!$      !1<=ib>=ia>=1
!!$      ! E.g. (ia,ib)= (1,1) (1,2) (2,2) (1,3) (2,3)  (3,3)
!!$      !                 1    2      3     4     5     6
!!$      !(1,ib) at 1 2 4 for ib=1,2,3
!!$      !        at 1 , 1+1 , 1+1+2 ... for 1,2,3
!!$      !        at 1+ (ib-1)*ib/2
!!$      !(ia,ib) at ia+(ib-1)*ib/2
!!$      ! max at ia=ib=n_species: n_species+(n_species-1)n_species/2
!!$      !                       =(n_species+1)n_species/2
!!$      ispec=ia+((ib-1)*ib)/2
    end function ispec
  end subroutine cgyro_init_landau
  
  subroutine old_sugama(cmat)

    use cgyro_globals, only : vth,temp,mass,dens,temp_ele,dens_ele,rho,z,&
         n_energy,e_max,n_xi,n_species,nt1,nt2,nc_loc_coll,nc_cl1,nc_cl2,nv,&
         xi,w_e,w_exi,e_deriv1_mat,xi_lor_mat,energy,dens2_rot,dens_rot,vel,vel2,nu,&
         collision_ene_diffusion,collision_mom_restore,collision_kperp,&
         collision_ene_restore,collision_ion_model,&
         is_v,ix_v,ie_v,ir_c,it_c,k_perp,bmag,is_ele,pi
    
    use mpi

    implicit none
    
    real,intent(out) :: cmat(:,:,:,:)
    real, dimension(:,:,:), allocatable :: nu_d, nu_par
    real, dimension(:,:), allocatable :: rs
    real, dimension(:,:,:,:), allocatable :: rsvec, rsvect0
    real, dimension(:,:), allocatable :: klor_fac, kdiff_fac
    real :: rsvtmp, rsvtmp0
    real :: my_dens2_rot,my_bj0,my_bj1

    real :: arg,xi_s1s,xi_prop
    real :: xa, xb, tauinv_ab
    real :: rval
    integer :: jv
    integer :: is,ir,it,ix,ie,iv,js,je,jx,ks
    integer :: ic,ic_loc,itor
    real, dimension(:,:), allocatable :: cmat_base1,cmat_base2
    real, dimension(:,:), allocatable :: cmat_loc
    real, dimension(:,:,:,:,:,:), allocatable :: ctest
    real, dimension(:,:,:,:), allocatable :: bessel
    allocate(nu_d(n_energy,n_species,n_species))
    allocate(nu_par(n_energy,n_species,n_species))
    allocate(klor_fac(n_species,n_species))
    allocate(kdiff_fac(n_species,n_species))
    nu_d(:,:,:) = 0.0
    nu_par(:,:,:) = 0.0
    klor_fac(:,:) = 0.0
    kdiff_fac(:,:) = 0.0

    do js=1,n_species
       do is=1,n_species
          do ie=1,n_energy
             xa = vel(ie)
             xb = xa*vth(is)/vth(js)
             tauinv_ab = nu(is)*z(js)**2/z(is)**2*dens(js)/dens(is)
             ! Ad hoc op
             ! (Fix for underflow)
             nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                  * (exp(-xb*xb)/(xb*sqrt(pi)) &
                  + (1.0-1.0/(2.0*xb*xb)) * erf(xb))
             ! No i-e Lorentz
             !if (is /= is_ele .and. js == is_ele) then
             !   nu_d(ie,is,js) = 0.0
             !endif
             if(collision_kperp == 1) then
                klor_fac(is,js) = 1.0
             endif

             ! Diffusion 
             nu_par(ie,is,js) = tauinv_ab * (2.0/xa**3) &
                  * (-exp(-xb*xb)/(xb*sqrt(pi)) &
                  + (1.0/(2.0*xb*xb)) * erf(xb))
             if(collision_kperp == 1) then
                kdiff_fac(is,js) = 1.0
             endif
          enddo
       enddo
    enddo

     if (collision_ion_model == 1) then
        do is=1,n_species
           if(is /= is_ele) then
              do js=1,n_species
                 nu_d(:,is,js) = 0.0
                 nu_par(:,is,js) = 0.0
              enddo
           endif
        enddo
     endif

     if (collision_kperp == 1 .and. &
         (collision_mom_restore == 1 .or. collision_ene_restore == 1)) then
       allocate(bessel(0:1,nv,nc_loc_coll,nt1:nt2))
!$omp parallel do collapse(2) private(ic_loc,it,ie,ix,is,iv,arg,xi_s1s,xi_prop)
       do itor=nt1,nt2
          do ic=nc_cl1,nc_cl2
             ic_loc = ic-nc_cl1+1
             it = it_c(ic)
             do iv=1,nv
                is = is_v(iv)
                ix = ix_v(iv)
                ie = ie_v(iv)
                xi_s1s = sqrt(1.0-xi(ix)**2)
                xi_prop = xi_s1s / xi(ix)
                arg = k_perp(ic,itor)*rho*vth(is)*mass(is)&
                     /(z(is)*bmag(it)) * vel2(ie) * xi_s1s
                bessel(0,iv,ic_loc,itor) = bessel_j0(arg)
                ! always used with the correction, so do it here
                bessel(1,iv,ic_loc,itor) = bessel_j1(arg) * xi_prop
             enddo
          enddo
       enddo
    endif
    allocate(ctest(n_species,n_species,n_xi,n_xi,n_energy,n_energy))
    ! Collision test particle component
    ! ctest

    ! Lorentz
!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
    do je=1,n_energy
       do ie=1,n_energy
          if (je == ie ) then
             do jx=1,n_xi
                do ix=1,n_xi
                   do js=1,n_species
                      do is=1,n_species
                         ctest(is,js,ix,jx,ie,je) = &
                              + xi_lor_mat(ix,jx) *0.5*nu_d(ie,is,js)
                      enddo
                   enddo
                enddo
             enddo
          else
             ctest(:,:,:,:,ie,je) = 0
          endif
       enddo
    enddo

    ! Diffusion
    if (collision_ene_diffusion == 1) then
!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
       do je=1,n_energy
          do ie=1,n_energy
             do ix=1,n_xi
                jx = ix
                do js=1,n_species
                   do is=1,n_species
                      ! From K. Hallatschek
                      ! self-adjoint part of ctest written self-adjointly
                      ctest(is,js,ix,jx,ie,je) = ctest(is,js,ix,jx,ie,je) &
                           -0.5 / w_e(ie) &
                           *sum(w_e(:)*e_deriv1_mat(:,ie)*energy(:) &
                           *nu_par(:,is,js) *e_deriv1_mat(:,je))/(1.0*e_max)
                      ! non-self-adjoint part proportional 1-Ta/Tb written
                      ! in a way that supports inherent particle number 
                      ! conservation for small kperp
                      ctest(is,js,ix,jx,ie,je) = ctest(is,js,ix,jx,ie,je) &
                           + (1-temp(is)/temp(js)) / sqrt(1.0*e_max)/w_e(ie) &
                           * w_e(je)*e_deriv1_mat(je,ie) &
                           * nu_par(je,is,js)*energy(je)**1.5
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    allocate(cmat_base1(nv,nv))
    allocate(cmat_base2(nv,nv))

    allocate(rs(n_species,n_species))
    allocate(rsvec(n_species,n_species,n_xi,n_energy))
    allocate(rsvect0(n_species,n_species,n_xi,n_energy))
    ! Momentum Restoring

    if (collision_mom_restore == 1) then

       ! C_test_ab(v_par f0a,f0b) and w_e v_par C_test_ab(H_a)
       rs(:,:) = 0.0

       do ie=1,n_energy
          do ix=1,n_xi
             do js=1,n_species
                do is=1,n_species
                   rsvtmp = 0
                   rsvtmp0 = 0
                   do je=1,n_energy
                      do jx=1,n_xi
                         rsvtmp = rsvtmp &
                              + ctest(is,js,ix,jx,ie,je) &
                              * vel2(je) * xi(jx)
                         rsvtmp0 = rsvtmp0 &
                              + ctest(is,js,jx,ix,je,ie) &
                              * vel2(je) * xi(jx) &
                              * w_exi(je,jx)
                      enddo
                   enddo
                   rsvtmp = rsvtmp * vth(is)
                   rsvtmp0 = rsvtmp0 * vth(is)
                   rsvec(is,js,ix,ie) = rsvtmp
                   rsvect0(is,js,ix,ie) = rsvtmp0
                   ! int v_par C_test_ab(v_par f0a,f0b) / n_0a
                   rs(is,js) = rs(is,js) + w_exi(ie,ix) * dens(is) &
                        * rsvtmp * vel2(ie) * xi(ix) &
                        * vth(is)
                enddo
             enddo
          enddo
       enddo
       
!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
       do jv=1,nv
          do iv=1,nv
             js = is_v(jv)
             jx = ix_v(jv)
             je = ie_v(jv)

             is = is_v(iv)
             ix = ix_v(iv)
             ie = ie_v(iv)

             if (abs(rs(is,js))>epsilon(0.0)) then
                cmat_base1(iv,jv) = &
                     - mass(js)/mass(is) &
                     * rsvec(is,js,ix,ie) / rs(is,js) &
                     * rsvect0(js,is,jx,je)
             else
                cmat_base1(iv,jv) = 0
             endif
          enddo
       enddo
    endif

    ! Energy Restoring
    
    if (collision_ene_restore == 1) then

       ! C_test_ab(u_a^2 f0a,f0b) and w_e u_a^2 C_test_ab(H_a)
       rs(:,:) = 0.0

       do ie=1,n_energy
          do ix=1,n_xi
             do js=1,n_species
                do is=1,n_species
                   rsvtmp = 0
                   rsvtmp0 = 0
                   do je=1,n_energy
                      do jx=1,n_xi
                         rsvtmp = rsvtmp &
                              + ctest(is,js,ix,jx,ie,je) * energy(je)
                         rsvtmp0 = rsvtmp0 &
                              + ctest(is,js,jx,ix,je,ie) * energy(je) &
                              * w_exi(je,jx)
                      enddo
                   enddo
                   rsvec(is,js,ix,ie) = rsvtmp
                   rsvect0(is,js,ix,ie) = rsvtmp0
                   ! int v^2 C_test_ab(u_a^2 f0a,f0b) 
                   rs(is,js) = rs(is,js) + w_exi(ie,ix) &
                        * dens(is) * rsvec(is,js,ix,ie) * (energy(ie)-1.5)
                   !Using temperature, not energy here. Matters for Ta/=Tb.
                enddo
             enddo
          enddo
       enddo
!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
       do jv=1,nv
          do iv=1,nv
             js = is_v(jv)
             jx = ix_v(jv)
             je = ie_v(jv)

             is = is_v(iv)
             ix = ix_v(iv)
             ie = ie_v(iv)

             if (abs(rs(is,js)) > epsilon(0.0)) then
                cmat_base2(iv,jv) = &
                     - temp(js)/temp(is) &
                     * rsvec(is,js,ix,ie) / rs(is,js) &
                     * rsvect0(js,is,jx,je)
             else
                cmat_base2(iv,jv) = 0
             endif
          enddo
       enddo
    endif
    deallocate(rs)
    deallocate(rsvec)
    deallocate(rsvect0)
    allocate(cmat_loc(nv,nv))
!$omp  parallel do collapse(2) default(none) &
!$omp& shared(nc_cl1,nc_cl2,nt1,nt2,nv,n_species,rho,is_ele,n_energy,n_xi) &
!$omp& shared(collision_kperp) &
!$omp& firstprivate(collision_mom_restore,collision_ene_restore) &
!$omp& shared(dens_ele,temp_ele,dens_rot,dens2_rot) &
!$omp& shared(cmat_base1,cmat_base2,bessel) &
!$omp& shared(it_c,ir_c,is_v,ix_v,ie_v,ctest) &
!$omp& shared(temp,dens,energy,vel,vel2) &
!$omp& shared(xi_lor_mat) &
!$omp& shared(k_perp,vth,mass,z,bmag,nu_d,xi,nu_par,w_e,w_exi) &
!$omp& shared(klor_fac,kdiff_fac) &
!$omp& private(ic,ic_loc,it,ir,rval,my_dens2_rot,my_bj0,my_bj1) &
!$omp& private(iv,is,ix,ie,jv,js,jx,je,ks) &
!$omp& private(cmat_loc) &
!$omp& shared(cmat)
    do itor=nt1,nt2
       do ic=nc_cl1,nc_cl2

          ic_loc = ic-nc_cl1+1

          it = it_c(ic)
          ir = ir_c(ic)
          ! Momentum Restoring

          if (collision_mom_restore == 1) then

             if (collision_kperp == 0) then
                do jv=1,nv
                   js = is_v(jv)
                   my_dens2_rot = dens2_rot(it,js)

                   do iv=1,nv

                      ! we know we are the first one to modify cmat_loc, so no need for +=
                      cmat_loc(iv,jv) = &
                           cmat_base1(iv,jv) &
                           * my_dens2_rot
                   enddo
                enddo

             else
                do jv=1,nv
                   js = is_v(jv)
                   my_dens2_rot = dens2_rot(it,js)
                   my_bj0 = bessel(0,jv,ic_loc,itor)
                   my_bj1 = bessel(1,jv,ic_loc,itor)

                   do iv=1,nv
                      
                      ! we know we are the first one to modify cmat_loc, so no need for +=
                      cmat_loc(iv,jv) = &
                           cmat_base1(iv,jv) &
                           * (  bessel(0,iv,ic_loc,itor) * my_bj0 &
                           + bessel(1,iv,ic_loc,itor) * my_bj1 ) &
                           * my_dens2_rot
                   enddo
                enddo
             endif

          else
             cmat_loc(:,:) = 0.0
          endif

          ! Energy Restoring

          if (collision_ene_restore == 1) then
             
             if (collision_kperp == 0) then
                do jv=1,nv
                   js = is_v(jv)
                   my_dens2_rot = dens2_rot(it,js)

                   do iv=1,nv

                      ! likely not the first, use +=
                      cmat_loc(iv,jv) &
                           = cmat_loc(iv,jv) &
                           + cmat_base2(iv,jv) &
                           * my_dens2_rot
                   enddo
                enddo

             else
                do jv=1,nv
                   js = is_v(jv)
                   my_dens2_rot = dens2_rot(it,js)
                   my_bj0 = bessel(0,jv,ic_loc,itor)

                   do iv=1,nv
                      ! likely not the first, use +=
                      cmat_loc(iv,jv) &
                           = cmat_loc(iv,jv) &
                           + cmat_base2(iv,jv) &
                           * bessel(0,iv,ic_loc,itor) * my_bj0 &
                           * my_dens2_rot
                   enddo
                enddo
             endif

          endif
          ! Change necessary to be able to compare with Landau operator: combine
          ! cmat and ctest right now here into one matrix cmat, including the kperp components.

          do iv=1,nv

             is = is_v(iv)
             ix = ix_v(iv)
             ie = ie_v(iv)

             do jv=1,nv

                js = is_v(jv)
                jx = ix_v(jv)
                je = ie_v(jv)

                ! Collision component: Test particle
                if (is == js) then
                   do ks=1,n_species
                      cmat_loc(iv,jv) = cmat_loc(iv,jv) &
                           + ctest(is,ks,ix,jx,ie,je) &
                           * dens_rot(it,ks)
                   enddo
                endif
                ! Finite-kperp test particle corrections 
                if (collision_kperp == 1) then
                   if (is == js .and. jx == ix .and. je == ie) then
                      do ks=1,n_species
                         cmat_loc(iv,jv) = cmat_loc(iv,jv) &
                              + (-0.25*(k_perp(ic,itor)*rho*vth(is)*mass(is) &
                              / (z(is)*bmag(it)))**2 * 2.0*energy(ie) &
                              * (klor_fac(is,ks)*nu_d(ie,is,ks) * (1+xi(ix)**2) &
                              + kdiff_fac(is,ks)*nu_par(ie,is,ks)* (1-xi(ix)**2)))
                      enddo
                   endif
                endif
             end do
          end do
          cmat(:,:,ic_loc,itor) = cmat_loc(:,:)
       enddo ! ic loop
    enddo ! itor loop
    deallocate(cmat_loc)

  end subroutine old_sugama

  subroutine qsort(sortidx,val,n)
    use cgyro_globals, only : i_proc
    implicit none
    integer,intent(in):: n,val(n)
    integer,intent(inout) :: sortidx(n)
    integer is,ie,i,j,ic,vic,s,isstack(n),iestack(n),sp
!!$    integer v
!!$    v(i)=val(sortidx(i)) !forced by fortran standard (LOL) to use obsolete statement function here.
    is=1
    ie=n
    sp=0
    do
       ic=(is+ie)/2
       vic=v(ic)
       i=is
       j=ie
!!$       if (i_proc==0 .and. verbose>1000) then
!!$          print *,'Sorting',is,ie
!!$       end if
       do
          do while (v(i)>vic)
             i=i+1
          end do
          do while (v(j)<vic)
             j=j-1
          end do
!!$          if (i_proc==0 .and. verbose>1000) then
!!$             print *,'found',i,j,v(i),v(j),vic
!!$          end if
          if (i<j) then
             s=sortidx(i)
             sortidx(i)=sortidx(j)
             sortidx(j)=s
             i=i+1
             j=j-1
          else
!!$             if (i_proc==0 .and. verbose>1000) then
!!$                print *,'stopped at',i,j,v(i),v(j)
!!$                do k=is,j
!!$                   if (v(k)<vic) then
!!$                      print *,'qsort intermediate error error at ',k,'v(k)',v(k),'vic',vic
!!$                      stop
!!$                   end if
!!$                end do
!!$                do k=i,ie
!!$                   if (v(k)<vic) then
!!$                      print *,'qsort intermediate error at ',k,'v(k)',v(k),'vic',vic
!!$                      stop
!!$                   end if
!!$                end do
!!$             end if
             exit
          end if
       end do
       ! now i>=j and all <i are <=vic all >j are >=vic
       ! the element at i is >= vic and the one at j is <=vic

       if (i-1>is) then
          if (j+1<ie) then
             sp=sp+1
             isstack(sp)=j+1
             iestack(sp)=ie
          end if
          ie=i-1
       else if (j+1<ie) then
          is=j+1
       else if (sp>0) then
          is=isstack(sp)
          ie=iestack(sp)
          sp=sp-1
       else
          exit
       end if
    end do
    if (verbose>=1000 .and. i_proc==0) then
       do i=1,n-1
          if (v(i)<v(i+1)) then
             print *,'qsort error at ',i,'v(i)',v(i),'v(i+1)',v(i+1)
             stop
          end if
       end do
       print *,'qsort correctly done',n
       do i=1,n
          print *,'i',i,sortidx(i),'v',v(i)
          if (v(i)==0) exit
       end do
    end if
  contains
    integer function v(i)
      implicit none
      integer, intent(in) ::i
      v=val(sortidx(i))
    end function v
  end subroutine qsort
end module cgyro_init_collision_landau
  
