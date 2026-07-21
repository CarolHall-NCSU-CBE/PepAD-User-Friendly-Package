module constant

    implicit none
! The Data type "groupdetails" is used to define an array to store the information on the residues and atoms.
	type groupdetails
		integer				:: cnum1, cnum2, cnum3
		character*4			:: atype1(20), atype2(60), atype3(20)
		character*4			:: gtype
		real				:: coo1(20,3), coo2(60,3), coo3(20,3)
	end type

	type energyparameters
		integer				:: iac1(20), iac2(60), iac3(20), atomid1(20), atomid2(60), atomid3(20)
		real				:: charge1(20), epsion1(20), r1(20), rborn1(20), fs1(20), dielecons1(20)
		real				:: charge2(60), epsion2(60), r2(60), rborn2(60), fs2(60), dielecons2(60)
		real				:: charge3(20), epsion3(20), r3(20), rborn3(20), fs3(20), dielecons3(20)
	end type
	
	type lib4aminoacid
		integer				:: cnum1, cnum2, cnum3
		integer				:: dihedralangle(34,4)		
		character*4			:: atype1(20), atype2(60), atype3(20)
		character*4			:: gtype
		real				:: coo1(20,3), coo2(60,3), coo3(20,3)
		integer				:: grade, rotanum
	end type
	
	type index4sidechain
		integer				:: class_No, member_No
	end type

	type conformer4sidechain
		real				:: member(15,3)
	end type

	type dihedralparameters
		integer*8 			:: iph(36), jph(36), kph(36), lph(36), multiply(36)
		real*8    			:: pk(36,4), pn(36,4), phase(36,4)
	end type

	type peptideassembly
		integer			    :: num4peptides
		integer, allocatable :: peptideID(:)
    end type

	type sasa_cache_type
		logical :: valid=.false.
		integer :: n_atoms=0
		real(kind=8) :: area_complex=0.0d0
		real(kind=8), allocatable :: area_group(:)
		real(kind=8), allocatable :: atom_area_complex(:)
		real(kind=8), allocatable :: atom_area_group(:)
		real(kind=8), allocatable :: coord(:,:),sasa_radius(:)
		integer, allocatable :: chain_id(:),site_id(:),group_id(:)
		character(len=4), allocatable :: atom_name(:)
	end type sasa_cache_type
	
    type hydration_types      ! new data type, contains hydration properties
        integer				:: hnum, pnum, cnum, onum
        character*4			:: hgtype(10), pgtype(10), cgtype(10), ogtype(10)
    end type
    
    type :: AminoRestriction
		character*4 :: gtype
		integer		:: max
        integer     :: min
        integer		:: count
        integer		:: clas
	end type AminoRestriction

	integer, parameter		:: maximum_aa_restrictions=20
	type(AminoRestriction)  :: aa_restrictions(maximum_aa_restrictions), &
							 aa_morethan(maximum_aa_restrictions)
	integer					:: n_restrictions, n_morethan		! number of restriction of AA
    
	integer					:: gnum, repeated_unit, num4category, num4betasheets
	integer					:: amino_acid_residue_count
	integer					:: atom_num

	integer, parameter		:: num4pdbsave=10
	integer, parameter		:: maximum_nmr_site_num=10
	integer, parameter		:: maximum_void_site_num=30
	integer					:: void_site(maximum_void_site_num)
	character*4				:: void_site_fixed_name(maximum_void_site_num)
	integer					:: void_site_num
	integer					:: nmr_site_num, nmr_site_ID(maximum_nmr_site_num)
    integer					:: restraint_num, exclusion_num, NMR_pool_size
	integer					:: nstep_start, nstep_terminal, idum, recalcu_switch, flag4sheet
	integer					:: fpho,fpol,fchg,foth
	integer					:: fpho_min, fpho_max, fpol_min, fpol_max
	integer					:: fchg_min, fchg_max, foth_min, foth_max
	logical					:: composition_min_set(4), composition_max_set(4)
	integer					:: anneal_stages, steps_before_sheetmove
	integer					:: sitenum4mutation
    integer, allocatable    :: sitenum4mutation_group(:)
	integer					:: sheetmove_interval
    integer					:: seed_switch, generate_traj_switch
	integer					:: sheetmove_flag, ESURF_flag, CONROT_flag
    !integer					:: ASN_limit, ASN_count

	real, parameter			:: scmfswitch1=0.6
!	real, parameter			:: scmfswitch2=0.8
	real, parameter			:: dihedral_weighting_factor=0.20
	!real, parameter        :: propensity_weighting_factor=3.0 ! Fixed-value alternative.
    real					:: propensity_weighting_factor
	real(kind=8), parameter	:: surftens=0.0072d0
	real(kind=8), parameter	:: offsetvalue=0.0d0
	real, parameter			:: vdw14_coeff=2.0  !1-4 VDW scale factor
	real, parameter			:: ele14_coeff=1.2  !1-4 ELE scale factor

	real					:: ekt_seq, ekt_seq_low, ekt_seq_high
    real					:: ekt_sheet, ekt_sheet_low, ekt_sheet_high
	real					:: energy_min(num4pdbsave)
    real         			:: sheet_switch
    real, parameter			:: CONROT_switch=0.6
	real					:: rmsd_max, rmsd_max_x, rmsd_max_y, rmsd_max_z
	real					:: displacement_factor, displacement_factor_x, displacement_factor_y, displacement_factor_z
	real					:: CONROT_angle_limit, CONROT_step, CONROT_closure_tol, CONROT_rmsd_limit
	!type RES4chain
	!	integer				:: num
	!	integer				:: IDs(gnum)
	!end type

	type(lib4aminoacid)		:: aa_lib(23)
	type(hydration_types)   :: hydrationprop
	type(peptideassembly), allocatable 	:: selfassembly(:)
    type(peptideassembly), allocatable 	:: betasheets(:)
	type(sasa_cache_type) :: sasa_current,sasa_trial
	integer :: sasa_validation_count=0
	character*50			:: filename							! The input PDB file
	character*4             :: NMR_restraint_AA(10), exclusion_AA(10), NMR_AA_pool(10)
	character(len=4), allocatable :: initial_residue_name(:)
    character(len=:), allocatable :: mydir						! The exe file located directory
    type(groupdetails),  allocatable  :: original_group(:,:)

	! DSSP hydrogen bonds and beta-sheet information.
	! Residue matrices use the index (chain-1)*gnum + residue.
	logical, allocatable :: hbond(:,:)
	logical, allocatable :: parallel_bridge(:,:), antiparallel_bridge(:,:)
	logical, allocatable :: parallel_ladder(:,:), antiparallel_ladder(:,:)
	logical, allocatable :: sheet_neighbor(:,:)
	logical, allocatable :: parallel_sheet_neighbor(:,:), antiparallel_sheet_neighbor(:,:)
	integer, allocatable :: chain_hbond_count(:,:)

end module constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module randomgenerator
	use constant
    contains
! Generate the random number
    
	subroutine init_random_seed()
        implicit none
        integer :: kkk, i, clock_count
        integer, allocatable :: seed(:)

        real :: ran

        call random_seed(size = kkk)
        allocate(seed(kkk))

		! inputfile normally generates the clock seed before printing its
		! summary.  Generate it here only if this routine is called directly.
        if(seed_switch.ne.1 .and. idum.eq.0) then
			call system_clock(count=clock_count)
			idum=modulo(clock_count,1000000000)
			if(idum.eq.0) idum=123457
		endif

        do i = 1, kkk
            seed(i) = modulo(idum, 1000000000)
            seed(i) = modulo(seed(i)+104729*i+7919*i*i, 1000000000)
            if(seed(i).eq.0) seed(i)=i+37
        end do
        
		call random_seed(put = seed)
		deallocate(seed)
        call random_number(ran) !discard first random value
        write(*,'(A,I0)') " Random seed = ", idum
        write(*,*) "the first random number:", ran
    end subroutine init_random_seed 
 
	subroutine ran_gen(ran2,flag_ran)
		implicit none
		real			:: ran2
		integer			:: flag_ran
    
		call random_number(ran2)
		return
    end subroutine ran_gen

end module randomgenerator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module pdbfile

	use	constant

	contains
	subroutine analyzepdb(user_gnum,user_chains)
	implicit none
	integer, intent(in), optional	:: user_gnum, user_chains
	integer						:: anum, status, ic, previous_ic
	integer						:: nres, candidate, repeat_count, i, j
	real						:: x, y, z
	character(len=6)			:: record_name
	character(len=4)			:: atype, name, aa_name
	character(len=150)			:: line
	character(len=4), allocatable	:: residue_names(:), enlarged_names(:)
	logical						:: first_residue, same_pattern, pattern_found

    if (recalcu_switch == 1) filename = 'pdbfiles/backup4pdb.pdb'
	nres=0
	first_residue=.true.
	allocate(residue_names(0))

	open(10, file=trim(filename), status='old', action='read', iostat=status)
	if(status.ne.0) stop "ERROR: cannot open the PDB file."
	do
		read(10,'(A)',iostat=status) line
		if(status.ne.0) exit
		if(len_trim(line).lt.3) cycle
		if(line(1:3)=='END') exit

		read(line,*,iostat=status) record_name, anum, atype, name, ic, x, y, z
		if(status.ne.0) cycle
		if(trim(record_name)/='ATOM' .and. trim(record_name)/='HETATM') cycle

		! Use the residue number because consecutive residues can have the
		! same name, for example the PHE-PHE pair in comp1.pdb.
		if(first_residue.or.ic.ne.previous_ic) then
			nres=nres+1
			allocate(enlarged_names(nres))
			if(nres.gt.1) enlarged_names(1:nres-1)=residue_names
			enlarged_names(nres)=name
			call move_alloc(enlarged_names,residue_names)
			previous_ic=ic
			first_residue=.false.
		endif
	enddo
	close(10)

	if(nres.eq.0) stop "ERROR: no amino-acid residues were found in the PDB file."

	if(present(user_gnum).neqv.present(user_chains)) then
		stop "ERROR: N_RESIDUE and N_CHAINS must be supplied together."
	endif

	if(present(user_gnum)) then
		! Use the dimensions supplied by the user.
		if(user_gnum.le.0.or.user_chains.le.0) then
			stop "ERROR: N_RESIDUE and N_CHAINS must be positive."
		endif
		if(user_gnum*user_chains.ne.nres) then
			stop "ERROR: N_RESIDUE*N_CHAINS does not match the PDB residue count."
		endif
		gnum=user_gnum
		repeated_unit=user_chains
	else
		! Without user dimensions, every chain must have the same residue pattern.
		pattern_found=.false.
		do candidate=1,nres/2
			if(mod(nres,candidate).ne.0) cycle
			repeat_count=nres/candidate
			same_pattern=.true.
			do i=candidate+1,nres
				j=mod(i-1,candidate)+1
				if(residue_names(i).ne.residue_names(j)) then
					same_pattern=.false.
					exit
				endif
			enddo
			if(same_pattern) then
				gnum=candidate
				repeated_unit=repeat_count
				pattern_found=.true.
				exit
			endif
		enddo
		if(.not.pattern_found) then
			stop "ERROR: PDB chains are different; supply N_RESIDUE and N_CHAINS."
		endif
	endif

	if (allocated(initial_residue_name)) deallocate(initial_residue_name)
	allocate(initial_residue_name(gnum))
	initial_residue_name=residue_names(1:gnum)

	! Record cap positions separately because they are also fixed mutation sites.
	do i=1,gnum
		aa_name=adjustl(residue_names(i))
		if (trim(aa_name)=='ACE'.or.trim(aa_name)=='NME'.or.trim(aa_name)=='NHE') then
			void_site_num=void_site_num+1
			if (void_site_num > maximum_void_site_num) &
				stop 'ERROR: too many cap and single-site constraints.'
			void_site(void_site_num)=i
		endif
	enddo

	! Use the first peptide chain to set omitted composition bounds.
	call count_amino_acid_categories(residue_names(1:gnum),fpho,fpol,fchg,foth)
	amino_acid_residue_count=fpho+fpol+fchg+foth

	if (.not.any(composition_min_set) .and. &
		.not.any(composition_max_set)) then
		! No composition flags: preserve the first-chain composition exactly.
		fpho_min=fpho; fpho_max=fpho
		fpol_min=fpol; fpol_max=fpol
		fchg_min=fchg; fchg_max=fchg
		foth_min=foth; foth_max=foth
	else
		! Custom mode: every category needs at least one supplied bound.
		if (.not.(composition_min_set(1).or.composition_max_set(1))) &
			stop 'ERROR: custom composition requires hydrophobic MIN or MAX.'
		if (.not.(composition_min_set(2).or.composition_max_set(2))) &
			stop 'ERROR: custom composition requires polar MIN or MAX.'
		if (.not.(composition_min_set(3).or.composition_max_set(3))) &
			stop 'ERROR: custom composition requires charged MIN or MAX.'
		if (.not.(composition_min_set(4).or.composition_max_set(4))) &
			stop 'ERROR: custom composition requires other MIN or MAX.'

		if (.not.composition_min_set(1)) fpho_min=0
		if (.not.composition_max_set(1)) fpho_max=gnum
		if (.not.composition_min_set(2)) fpol_min=0
		if (.not.composition_max_set(2)) fpol_max=gnum
		if (.not.composition_min_set(3)) fchg_min=0
		if (.not.composition_max_set(3)) fchg_max=gnum
		if (.not.composition_min_set(4)) foth_min=0
		if (.not.composition_max_set(4)) foth_max=gnum
	endif
    
	deallocate(residue_names)
	return
	end subroutine analyzepdb

	subroutine count_amino_acid_categories(residue_names,n_hydrophobic,n_polar,n_charged,n_other)
	implicit none
	character(len=*), intent(in) :: residue_names(:)
	integer, intent(out) :: n_hydrophobic, n_polar, n_charged, n_other
	integer :: residue
	character(len=4) :: aa_name

	n_hydrophobic=0
	n_polar=0
	n_charged=0
	n_other=0

	do residue=1,size(residue_names)
		aa_name=adjustl(residue_names(residue))

		! Caps are not amino acids and do not contribute to composition.
		if (trim(aa_name) == 'ACE' .or. trim(aa_name) == 'NME' .or. &
			trim(aa_name) == 'NHE') cycle

		! Normalize terminal residue names such as NLYS and CGLU.
		if (len_trim(aa_name) == 4) then
			if (aa_name(1:1) == 'N' .or. aa_name(1:1) == 'C') &
				aa_name=aa_name(2:4)
		endif

		select case(trim(aa_name))
		case('ALA','LEU','VAL','ILE','MET','PHE','TYR','TRP','GLY')
			n_hydrophobic=n_hydrophobic+1
		case('ASN','GLN','SER','THR','HIE','HIS','HID','HIP')
			n_polar=n_polar+1
		case('ARG','LYS','GLU','ASP')
			n_charged=n_charged+1
		case('PRO','CYS')
			n_other=n_other+1
		case default
			write(*,'(3A)') 'ERROR: unknown residue ',trim(aa_name), &
				' while determining sequence composition.'
			stop
		end select
	enddo

	end subroutine count_amino_acid_categories



! Read an initial file.
	subroutine readpdb(group,structure_only)
	implicit none
	integer						:: anum, status, chainID, ic, i, flag
	real						:: x, y, z
	character*4					:: char, atype, name
	type(groupdetails)			:: group(repeated_unit,gnum)
	logical, intent(in), optional :: structure_only
	logical :: load_structure_only
	character*150               :: line

	load_structure_only=.false.
	if (present(structure_only)) load_structure_only=structure_only
    
	group%cnum1=0
	group%cnum2=0
	group%cnum3=0
	chainID=1
    flag=0

    if (recalcu_switch == 1) filename = 'pdbfiles/backup4pdb.pdb'
    
	open(10, file=filename)
	do while(.true.)

        read(10, '(A)', iostat=status) line
		if(status.ne.0) exit
        if(line(1:3) == 'TER') then
            goto 10
        else if (line(1:3) == 'END') then
            exit
        endif
        read(line, *, iostat=status) char, anum, atype, name, ic, x, y, z
        if(status.ne.0) exit
		if(ic.gt.(chainID*gnum)) chainID=chainID+1 ! Identify the peptide chain.
		ic=ic-(chainID-1)*gnum						! Convert to the residue position within the chain.
        
        if(name /= "ACE" .and. name /= "NME" .and. name /= "NHE" .and. len_trim(name)==3) then !judge if current AA is not label as N or C
			if (ic==1) then 
                name = "N"//name 
            elseif (ic==gnum) then
                name = "C"//name 
            endif
        endif
        
		group(chainID,ic)%gtype=name
		if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
		   .or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then  
			group(chainID,ic)%cnum1=group(chainID,ic)%cnum1+1
			group(chainID,ic)%atype1(group(chainID,ic)%cnum1)=atype
			group(chainID,ic)%coo1(group(chainID,ic)%cnum1,1)=x
			group(chainID,ic)%coo1(group(chainID,ic)%cnum1,2)=y
			group(chainID,ic)%coo1(group(chainID,ic)%cnum1,3)=z
		elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
			group(chainID,ic)%cnum3=group(chainID,ic)%cnum3+1
			group(chainID,ic)%atype3(group(chainID,ic)%cnum3)=atype
			group(chainID,ic)%coo3(group(chainID,ic)%cnum3,1)=x
			group(chainID,ic)%coo3(group(chainID,ic)%cnum3,2)=y
			group(chainID,ic)%coo3(group(chainID,ic)%cnum3,3)=z
		else
			group(chainID,ic)%cnum2=group(chainID,ic)%cnum2+1
			group(chainID,ic)%atype2(group(chainID,ic)%cnum2)=atype
			group(chainID,ic)%coo2(group(chainID,ic)%cnum2,1)=x
			group(chainID,ic)%coo2(group(chainID,ic)%cnum2,2)=y
			group(chainID,ic)%coo2(group(chainID,ic)%cnum2,3)=z
		endif

		if (.not.load_structure_only) then
			if(chainID==1) then
				if(ic.ne.flag) then
					flag=ic
					do i=1, void_site_num
						if(ic == void_site(i)) then
							goto 10
						endif
					enddo
					sitenum4mutation=sitenum4mutation+1				! record how many sites for mutation 
					sitenum4mutation_group(sitenum4mutation)=ic     ! Record the mutable residue position in the peptide.
				endif
			endif
		endif
10		continue
		!write(*,*) chainID, ic, group(chainID,ic)%gtype
	enddo
	
	close(10)
	
	if (.not.load_structure_only) then
		if (sitenum4mutation /= gnum-void_site_num) then
			stop 'ERROR: mutation-site count is inconsistent with positional constraints.'
		endif
		original_group=group
	endif
	
	return
	end subroutine readpdb


	subroutine generatepdb(step, attempt, group)
	implicit none
	integer							:: step, attempt
	integer							:: i, j, k, atomnum
	character*5						:: stepchar, attemptchar
	type(groupdetails)				:: group(repeated_unit,gnum)

	write(stepchar,"(i5)") step
	write(attemptchar,"(i4)") attempt
            open(10,file='pdbfiles/'//trim(adjustl(stepchar))//'_'//trim(adjustl(attemptchar))//'.pdb', status="replace", access="sequential")
		atomnum=1
		do i=1, repeated_unit
			do j=1, gnum
				do k=1, group(i,j)%cnum1
					if(len_trim(group(i,j)%atype1(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype1(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo1(k,1), group(i,j)%coo1(k,2), group(i,j)%coo1(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype1(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo1(k,1), group(i,j)%coo1(k,2), group(i,j)%coo1(k,3)
					endif
					atomnum=atomnum+1
				enddo
				do k=1, group(i,j)%cnum2
					if(len_trim(group(i,j)%atype2(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype2(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo2(k,1), group(i,j)%coo2(k,2), group(i,j)%coo2(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype2(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo2(k,1), group(i,j)%coo2(k,2), group(i,j)%coo2(k,3)
					endif
					atomnum=atomnum+1
				enddo
				do k=1, group(i,j)%cnum3
					if(len_trim(group(i,j)%atype3(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype3(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo3(k,1), group(i,j)%coo3(k,2), group(i,j)%coo3(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype3(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo3(k,1), group(i,j)%coo3(k,2), group(i,j)%coo3(k,3)
					endif
					atomnum=atomnum+1
				enddo
			enddo
		enddo
1		format(a4,i7,a6,a4,i5,f12.3,2f8.3,2f6.2,a16)
2		format(a4,i7,a5,a1,a4,i5,f12.3,2f8.3,2f6.2,a16)
	close(10)

	return
    end subroutine generatepdb


	subroutine generatepdb_debug(prefix, step, attempt, group)
	implicit none
	integer							:: step, attempt
	integer							:: i, j, k, atomnum, model_id
	character(len=*)				:: prefix
	character*5						:: stepchar, attemptchar
	type(groupdetails)				:: group(repeated_unit,gnum)

	write(stepchar,"(i5)") step
	write(attemptchar,"(i4)") attempt
	if (trim(adjustl(prefix)) == 'scmf_failure_before') then
		open(10,file='pdbfiles/scmf_failure_before.pdb',status='replace',access='sequential')
		model_id=1
		write(10,'(a5,i9)') 'MODEL',model_id
		write(10,'(a,i0,a,i0)') 'REMARK SCMF FAILURE BEFORE ATTEMPT; CHAIN SITE ', &
			step,'; FAILED CHAIN ',attempt
	elseif (trim(adjustl(prefix)) == 'scmf_failure_candidate') then
		open(10,file='pdbfiles/scmf_failure_candidate.pdb',status='replace',access='sequential')
		model_id=1
		write(10,'(a5,i9)') 'MODEL',model_id
		write(10,'(a,i0,a,i0)') 'REMARK SCMF FINAL REJECTED ROTAMER; CHAIN SITE ', &
			step,'; FAILED CHAIN ',attempt
	else
		open(10,file='pdbfiles/conrot_debug.pdb', access="append")
		model_id=1
		if(trim(adjustl(prefix))=="conrot_repack_failed") model_id=2
		if(trim(adjustl(prefix))=="conrot_repack") model_id=3
		write(10,'(a5,i9)') "MODEL", model_id
		write(10,'(a,1x,a)') "REMARK CONROT_DEBUG", trim(adjustl(prefix))
	endif
	!if(trim(adjustl(prefix))=="conrot_failed".or.trim(adjustl(prefix))=="conrot_raw") then
	!	open(10,file='pdbfiles/conrot_debug_'//trim(adjustl(stepchar))//'_'// &
	!		trim(adjustl(attemptchar))//'.pdb', status="replace", access="sequential")
	!else
	!	open(10,file='pdbfiles/conrot_debug_'//trim(adjustl(stepchar))//'_'// &
	!		trim(adjustl(attemptchar))//'.pdb', status="unknown", access="sequential", position="append")
	!endif
		atomnum=1
		do i=1, repeated_unit
			do j=1, gnum
				do k=1, group(i,j)%cnum1
					if(len_trim(group(i,j)%atype1(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype1(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo1(k,1), group(i,j)%coo1(k,2), group(i,j)%coo1(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype1(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo1(k,1), group(i,j)%coo1(k,2), group(i,j)%coo1(k,3)
					endif
					atomnum=atomnum+1
				enddo
				do k=1, group(i,j)%cnum2
					if(len_trim(group(i,j)%atype2(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype2(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo2(k,1), group(i,j)%coo2(k,2), group(i,j)%coo2(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype2(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo2(k,1), group(i,j)%coo2(k,2), group(i,j)%coo2(k,3)
					endif
					atomnum=atomnum+1
				enddo
				do k=1, group(i,j)%cnum3
					if(len_trim(group(i,j)%atype3(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype3(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo3(k,1), group(i,j)%coo3(k,2), group(i,j)%coo3(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype3(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo3(k,1), group(i,j)%coo3(k,2), group(i,j)%coo3(k,3)
					endif
					atomnum=atomnum+1
				enddo
			enddo
		enddo
	write(10,'(a6)') "ENDMDL"
1		format(a4,i7,a6,a4,i5,f12.3,2f8.3,2f6.2,a16)
2		format(a4,i7,a5,a1,a4,i5,f12.3,2f8.3,2f6.2,a16)
	close(10)

	return
	end subroutine generatepdb_debug
    
    subroutine generatebackup4pdb(group)
	implicit none
	integer							:: i, j, k, atomnum
	character*5						:: stepchar, attemptchar
	type(groupdetails)				:: group(repeated_unit,gnum)

	open(10,file='pdbfiles/backup4pdb.pdb', access="sequential", status="replace")
		atomnum=1
		do i=1, repeated_unit
			do j=1, gnum
				do k=1, group(i,j)%cnum1
					if(len_trim(group(i,j)%atype1(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype1(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo1(k,1), group(i,j)%coo1(k,2), group(i,j)%coo1(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype1(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo1(k,1), group(i,j)%coo1(k,2), group(i,j)%coo1(k,3)
					endif
					atomnum=atomnum+1
				enddo
				do k=1, group(i,j)%cnum2
					if(len_trim(group(i,j)%atype2(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype2(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo2(k,1), group(i,j)%coo2(k,2), group(i,j)%coo2(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype2(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo2(k,1), group(i,j)%coo2(k,2), group(i,j)%coo2(k,3)
					endif
					atomnum=atomnum+1
				enddo
				do k=1, group(i,j)%cnum3
					if(len_trim(group(i,j)%atype3(k))==4) then
						write(10,2) "ATOM", atomnum, group(i,j)%atype3(k), " ", group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo3(k,1), group(i,j)%coo3(k,2), group(i,j)%coo3(k,3)
					else
						write(10,1) "ATOM", atomnum, group(i,j)%atype3(k), group(i,j)%gtype, (i-1)*gnum+j, group(i,j)%coo3(k,1), group(i,j)%coo3(k,2), group(i,j)%coo3(k,3)
					endif
					atomnum=atomnum+1
				enddo
			enddo
		enddo
1		format(a4,i7,a6,a4,i5,f12.3,2f8.3,2f6.2,a16)
2		format(a4,i7,a5,a1,a4,i5,f12.3,2f8.3,2f6.2,a16)
	close(10)

	return
    end subroutine generatebackup4pdb
    
    subroutine Hbond_matrix
	implicit none
	integer :: nres, chain_i, chain_j, site_i, site_j
	integer :: residue_i, residue_j, atom, donor_hydrogen
	integer :: incomplete_backbone_count
	integer :: queue_head, queue_tail, current_chain, neighbor_chain, sheet_id_number
	integer, allocatable :: sheet_id(:), energy_group(:), chain_queue(:)
	real(kind=8) :: rON, rCH, rOH, rCN, energy
	real(kind=8), parameter :: minimum_distance=1.0d-6
	type(groupdetails) :: group(repeated_unit,gnum)
	real(kind=8), allocatable :: n_coo(:,:), h_coo(:,:,:), c_coo(:,:), o_coo(:,:)
	integer, allocatable :: h_count(:)
	logical, allocatable :: has_n(:), has_c(:), has_o(:), is_cap(:), is_n_terminal(:)
	logical :: is_parallel, is_antiparallel
	character(len=4) :: residue_name

	! At this stage PepAD only needs coordinates for DSSP analysis.
	! Mutation and simulation arrays are initialized later.
	call readpdb(group,.true.)
	nres=repeated_unit*gnum

	allocate(n_coo(nres,3),h_coo(nres,3,3),c_coo(nres,3),o_coo(nres,3))
	allocate(h_count(nres))
	allocate(has_n(nres),has_c(nres),has_o(nres),is_cap(nres),is_n_terminal(nres))
	n_coo=0.0d0; h_coo=0.0d0; c_coo=0.0d0; o_coo=0.0d0
	h_count=0
	has_n=.false.; has_c=.false.; has_o=.false.
	is_cap=.false.; is_n_terminal=.false.

	! Put the DSSP donor and acceptor atoms into residue-indexed coordinate arrays.
	do chain_i=1,repeated_unit
		do site_i=1,gnum
			residue_i=(chain_i-1)*gnum+site_i
			residue_name=adjustl(group(chain_i,site_i)%gtype)
			select case(trim(residue_name))
			case('ACE','NME','NHE')
				is_cap(residue_i)=.true.
			end select
			if (len_trim(residue_name) == 4 .and. residue_name(1:1) == 'N') &
				is_n_terminal(residue_i)=.true.
			do atom=1,group(chain_i,site_i)%cnum1
				select case(trim(adjustl(group(chain_i,site_i)%atype1(atom))))
				case('N')
					n_coo(residue_i,:)=group(chain_i,site_i)%coo1(atom,:)
					has_n(residue_i)=.true.
				case('H')
					if (.not.is_n_terminal(residue_i)) then
						h_count(residue_i)=1
						h_coo(residue_i,1,:)=group(chain_i,site_i)%coo1(atom,:)
					endif
				case('H1','H2','H3')
					if (is_n_terminal(residue_i) .and. h_count(residue_i) < 3) then
						h_count(residue_i)=h_count(residue_i)+1
						h_coo(residue_i,h_count(residue_i),:)= &
							group(chain_i,site_i)%coo1(atom,:)
					endif
				end select
			enddo
			do atom=1,group(chain_i,site_i)%cnum3
				select case(trim(adjustl(group(chain_i,site_i)%atype3(atom))))
				case('C')
					c_coo(residue_i,:)=group(chain_i,site_i)%coo3(atom,:)
					has_c(residue_i)=.true.
				case('O')
					o_coo(residue_i,:)=group(chain_i,site_i)%coo3(atom,:)
					has_o(residue_i)=.true.
				end select
			enddo
		enddo
	enddo

	incomplete_backbone_count=count(.not.is_cap .and. &
		.not.(has_n.and.(h_count > 0).and.has_c.and.has_o))
	if (incomplete_backbone_count > 0) then
		write(*,'(A,I0,A)') 'WARNING: ',incomplete_backbone_count, &
			' non-cap residues are missing DSSP atoms (N, C, O, and H or H1/H2/H3).'
		write(*,'(A)') '         Beta-sheet detection cannot be performed reliably for this PDB file.'
	endif

	if (allocated(hbond)) deallocate(hbond)
	if (allocated(parallel_bridge)) deallocate(parallel_bridge)
	if (allocated(antiparallel_bridge)) deallocate(antiparallel_bridge)
	if (allocated(parallel_ladder)) deallocate(parallel_ladder)
	if (allocated(antiparallel_ladder)) deallocate(antiparallel_ladder)
	if (allocated(sheet_neighbor)) deallocate(sheet_neighbor)
	if (allocated(parallel_sheet_neighbor)) deallocate(parallel_sheet_neighbor)
	if (allocated(antiparallel_sheet_neighbor)) deallocate(antiparallel_sheet_neighbor)
	if (allocated(chain_hbond_count)) deallocate(chain_hbond_count)

	allocate(hbond(nres,nres))
	allocate(parallel_bridge(nres,nres),antiparallel_bridge(nres,nres))
	allocate(parallel_ladder(nres,nres),antiparallel_ladder(nres,nres))
	allocate(sheet_neighbor(repeated_unit,repeated_unit))
	allocate(parallel_sheet_neighbor(repeated_unit,repeated_unit))
	allocate(antiparallel_sheet_neighbor(repeated_unit,repeated_unit))
	allocate(chain_hbond_count(repeated_unit,repeated_unit))
	hbond=.false.
	parallel_bridge=.false.; antiparallel_bridge=.false.
	parallel_ladder=.false.; antiparallel_ladder=.false.
	sheet_neighbor=.false.
	parallel_sheet_neighbor=.false.; antiparallel_sheet_neighbor=.false.
	chain_hbond_count=0

	! hbond(donor,acceptor) is directed.  Caps and directly connected
	! same-chain residues are excluded from the DSSP calculation.
	do residue_i=1,nres
		if (is_cap(residue_i)) cycle
		if (.not.has_n(residue_i) .or. h_count(residue_i) == 0) cycle
		chain_i=(residue_i-1)/gnum+1
		site_i=mod(residue_i-1,gnum)+1
		do residue_j=1,nres
			if (is_cap(residue_j)) cycle
			if (.not.(has_c(residue_j).and.has_o(residue_j))) cycle
			chain_j=(residue_j-1)/gnum+1
			site_j=mod(residue_j-1,gnum)+1
			if (chain_i == chain_j .and. abs(site_i-site_j) <= 2) cycle

			rON=sqrt(sum((o_coo(residue_j,:)-n_coo(residue_i,:))**2))
			rCN=sqrt(sum((c_coo(residue_j,:)-n_coo(residue_i,:))**2))
			if (min(rON,rCN) <= minimum_distance) cycle

			! Ordinary residues use H.  An N*** terminal residue uses each of
			! H1, H2, and H3; any one passing the cutoff defines the directed bond.
			do donor_hydrogen=1,h_count(residue_i)
				rCH=sqrt(sum((c_coo(residue_j,:)- &
					h_coo(residue_i,donor_hydrogen,:))**2))
				rOH=sqrt(sum((o_coo(residue_j,:)- &
					h_coo(residue_i,donor_hydrogen,:))**2))
				if (min(rCH,rOH) <= minimum_distance) cycle
				energy=27.888d0*(1.0d0/rON+1.0d0/rCH-1.0d0/rOH-1.0d0/rCN)
				if (energy < -0.5d0) then
					hbond(residue_i,residue_j)=.true.
					exit
				endif
			enddo
		enddo
	enddo

	! Detect DSSP parallel and antiparallel bridges between different chains.
	! The local-site checks prevent every +/-1 offset from crossing a chain end.
	do chain_i=1,repeated_unit-1
		do chain_j=chain_i+1,repeated_unit
			do site_i=1,gnum
				residue_i=(chain_i-1)*gnum+site_i
				do site_j=1,gnum
					residue_j=(chain_j-1)*gnum+site_j
					is_antiparallel=hbond(residue_i,residue_j).and. &
						hbond(residue_j,residue_i)
					if (site_i > 1 .and. site_i < gnum .and. &
						site_j > 1 .and. site_j < gnum) then
						is_antiparallel=is_antiparallel.or. &
							hbond(residue_i-1,residue_j+1).and. &
							hbond(residue_j-1,residue_i+1)
					endif

					is_parallel=.false.
					if (site_i > 1 .and. site_i < gnum) then
						is_parallel=hbond(residue_i-1,residue_j).and. &
							hbond(residue_j,residue_i+1)
					endif
					if (site_j > 1 .and. site_j < gnum) then
						is_parallel=is_parallel.or. &
							hbond(residue_j-1,residue_i).and. &
							hbond(residue_i,residue_j+1)
					endif

					if (is_antiparallel) then
						antiparallel_bridge(residue_i,residue_j)=.true.
						antiparallel_bridge(residue_j,residue_i)=.true.
					endif
					if (is_parallel) then
						parallel_bridge(residue_i,residue_j)=.true.
						parallel_bridge(residue_j,residue_i)=.true.
					endif
				enddo
			enddo
		enddo
	enddo

	! Two adjacent bridges are the minimum ladder.  Longer runs are merged
	! because every neighboring bridge pair marks both of its bridge positions.
	do chain_i=1,repeated_unit-1
		do chain_j=chain_i+1,repeated_unit
			do site_i=1,gnum-1
				residue_i=(chain_i-1)*gnum+site_i
				do site_j=1,gnum-1
					residue_j=(chain_j-1)*gnum+site_j
					if (parallel_bridge(residue_i,residue_j).and. &
						parallel_bridge(residue_i+1,residue_j+1)) then
						parallel_ladder(residue_i,residue_j)=.true.
						parallel_ladder(residue_j,residue_i)=.true.
						parallel_ladder(residue_i+1,residue_j+1)=.true.
						parallel_ladder(residue_j+1,residue_i+1)=.true.
						parallel_sheet_neighbor(chain_i,chain_j)=.true.
						parallel_sheet_neighbor(chain_j,chain_i)=.true.
					endif
				enddo

				do site_j=2,gnum
					residue_j=(chain_j-1)*gnum+site_j
					if (antiparallel_bridge(residue_i,residue_j).and. &
						antiparallel_bridge(residue_i+1,residue_j-1)) then
						antiparallel_ladder(residue_i,residue_j)=.true.
						antiparallel_ladder(residue_j,residue_i)=.true.
						antiparallel_ladder(residue_i+1,residue_j-1)=.true.
						antiparallel_ladder(residue_j-1,residue_i+1)=.true.
						antiparallel_sheet_neighbor(chain_i,chain_j)=.true.
						antiparallel_sheet_neighbor(chain_j,chain_i)=.true.
					endif
				enddo
			enddo
		enddo
	enddo

	! Count all directed backbone H-bonds between each pair of chains.
	! Three or more H-bonds are sufficient to connect the two strands, even
	! when the strict DSSP bridge pattern does not form a two-bridge ladder.
	do chain_i=1,repeated_unit-1
		do chain_j=chain_i+1,repeated_unit
			do site_i=1,gnum
				residue_i=(chain_i-1)*gnum+site_i
				do site_j=1,gnum
					residue_j=(chain_j-1)*gnum+site_j
					if (hbond(residue_i,residue_j)) &
						chain_hbond_count(chain_i,chain_j)= &
						chain_hbond_count(chain_i,chain_j)+1
					if (hbond(residue_j,residue_i)) &
						chain_hbond_count(chain_i,chain_j)= &
						chain_hbond_count(chain_i,chain_j)+1
				enddo
			enddo
			chain_hbond_count(chain_j,chain_i)=chain_hbond_count(chain_i,chain_j)
		enddo
	enddo

	sheet_neighbor=parallel_sheet_neighbor.or.antiparallel_sheet_neighbor
	do chain_i=1,repeated_unit-1
		do chain_j=chain_i+1,repeated_unit
			if (chain_hbond_count(chain_i,chain_j) >= 3) then
				sheet_neighbor(chain_i,chain_j)=.true.
				sheet_neighbor(chain_j,chain_i)=.true.
			endif
		enddo
	enddo

	! Layer construction: every connected component of the
	! strand-neighbor map is one beta sheet.  At the same time, two-color
	! the graph so directly neighboring strands enter opposite energy groups.
	allocate(sheet_id(repeated_unit),energy_group(repeated_unit))
	allocate(chain_queue(repeated_unit))
	sheet_id=0
	energy_group=0
	num4betasheets=0

	do chain_i=1,repeated_unit
		if (sheet_id(chain_i) /= 0) cycle
		num4betasheets=num4betasheets+1
		queue_head=1
		queue_tail=1
		chain_queue(1)=chain_i
		sheet_id(chain_i)=num4betasheets
		energy_group(chain_i)=1

		do while(queue_head <= queue_tail)
			current_chain=chain_queue(queue_head)
			queue_head=queue_head+1
			do neighbor_chain=1,repeated_unit
				if (.not.sheet_neighbor(current_chain,neighbor_chain)) cycle
				if (sheet_id(neighbor_chain) == 0) then
					sheet_id(neighbor_chain)=num4betasheets
					energy_group(neighbor_chain)=3-energy_group(current_chain)
					queue_tail=queue_tail+1
					chain_queue(queue_tail)=neighbor_chain
				elseif (energy_group(neighbor_chain) == energy_group(current_chain)) then
					write(*,'(A,I0,A,I0,A)') 'ERROR: chains ',current_chain,' and ', &
						neighbor_chain,' make a non-alternating beta-sheet connectivity map.'
					stop
				endif
			enddo
		enddo
	enddo

	if (allocated(betasheets)) deallocate(betasheets)
	allocate(betasheets(num4betasheets))
	do sheet_id_number=1,num4betasheets
		betasheets(sheet_id_number)%num4peptides=0
		allocate(betasheets(sheet_id_number)%peptideID(repeated_unit))
		betasheets(sheet_id_number)%peptideID=0
	enddo

	if (allocated(selfassembly)) deallocate(selfassembly)
	allocate(selfassembly(2))
	do chain_i=1,2
		selfassembly(chain_i)%num4peptides=0
		allocate(selfassembly(chain_i)%peptideID(repeated_unit))
		selfassembly(chain_i)%peptideID=0
	enddo

	do chain_i=1,repeated_unit
		sheet_id_number=sheet_id(chain_i)
		betasheets(sheet_id_number)%num4peptides= &
			betasheets(sheet_id_number)%num4peptides+1
		betasheets(sheet_id_number)%peptideID( &
			betasheets(sheet_id_number)%num4peptides)=chain_i

		chain_j=energy_group(chain_i)
		selfassembly(chain_j)%num4peptides=selfassembly(chain_j)%num4peptides+1
		selfassembly(chain_j)%peptideID(selfassembly(chain_j)%num4peptides)=chain_i
	enddo

	write(6,*)
	write(6,'(A)') 'Beta-strand connectivity map'
	do chain_i=1,repeated_unit
		write(6,'(A,I0,A)',advance='no') ' Chain ',chain_i,' ->'
		if (.not.any(sheet_neighbor(chain_i,:))) then
			write(6,'(A)') ' none'
		else
			do chain_j=1,repeated_unit
				if (sheet_neighbor(chain_i,chain_j)) &
					write(6,'(1X,I0)',advance='no') chain_j
			enddo
			write(6,*)
		endif
	enddo

	deallocate(sheet_id,energy_group,chain_queue)
	deallocate(n_coo,h_coo,c_coo,o_coo,h_count)
	deallocate(has_n,has_c,has_o,is_cap,is_n_terminal)

    end subroutine Hbond_matrix
    

end module pdbfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module input

	use constant
	use pdbfile, only: analyzepdb, Hbond_matrix

    contains
    
	subroutine inputfile(input_name)
	implicit none
	integer, parameter :: maximum_user_sheets=100
	integer, parameter :: number_of_counted_flags=37
	character(len=3), parameter :: amino_acid_names(20) = (/ &
		'GLY','ALA','VAL','LEU','ILE','MET','PHE','TYR','TRP','SER', &
		'ASN','GLN','THR','HIE','ARG','LYS','GLU','ASP','CYS','PRO' /)
	integer, parameter :: amino_acid_categories(20) = (/ &
		1,1,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4 /)
	character(len=32), parameter :: counted_flag_names(number_of_counted_flags) = &
		(/ character(len=32) :: &
		'PDBFILE','N_RESIDUES','N_CHAINS','GROUP_1','GROUP_2','RESTART', &
		'RANSEED','N_STEPS','KBT_SEQ','KBT_SEQ_HIGH','KBT_SEQ_LOW', &
		'KBT_SHEETMOVE','KBT_SHEETMOVE_HIGH','KBT_SHEETMOVE_LOW', &
		'RMSD_MAX','PAGG_WEIGHT','SHEETMOVE','ANNEAL_STAGES','ESURF_MODE', &
		'STEPS_BEFORE_SHEETMOVE','STEPS_BETWEEN_SHEETMOVE', &
		'N_HYDROPHOBIC_MIN','N_HYDROPHOBIC_MAX','N_POLAR_MIN','N_POLAR_MAX', &
		'N_CHARGED_MIN','N_CHARGED_MAX','N_OTHER_MIN','N_OTHER_MAX', &
		'SINGLE_SITE_CONSTRAINTS','GROUPED_SITE_AA_POOL', &
		'GROUPED_SITE_CONSTRAINTS','CONROT_FLAG','CONROT_ANGLE_LIMIT', &
		'CONROT_STEP','CONROT_CLOSURE_TOL','CONROT_RMSD_LIMIT' /)
	character(len=*), intent(in), optional :: input_name
	integer :: unit, ios, line_number, equals_at, expected_gnum, expected_chains
	integer :: clock_count
	integer :: i, j, k, min_value, max_value, group_index, sheet_number, aa_index
	integer :: counted_flag_index
	integer :: max_user_sheet, token_count, chain_id
	integer :: specific_minimums(4), specific_maximum_capacity(4)
	integer :: category_minimum_requirements(4), category_maximum_capacity(4)
	integer :: fixed_aa_counts(20), user_aa_minimums(20), user_aa_maximums(20)
	integer :: effective_aa_minimums(20), cap_site_count
	integer :: fixed_category_counts(4), grouped_pool_category_counts(4)
	integer :: grouped_category_minimums(4), grouped_category_maximums(4)
	integer :: positional_category_minimums(4), positional_category_maximums(4)
	integer :: expected_mutation_sites
	integer :: flag_count(number_of_counted_flags)
	integer, allocatable :: chain_counts(:)
	logical :: matching_minimum, in_token, user_groups_used, user_sheets_used
	logical :: amino_acid_flag, aa_min_defined(20), aa_max_defined(20)
	character(len=256) :: input_path, line, key, value, iomsg
	character(len=256) :: single_site_value, group_value(2)
	character(len=256) :: sheet_value(maximum_user_sheets)
	logical :: sheet_defined(maximum_user_sheets)
	character(len=256) :: grouped_pool_value, grouped_site_value
	character(len=20) :: restriction_label
	character(len=40) :: sheet_label
	character(len=4) :: fixed_residue_name

	input_path='input.txt'
	if (present(input_name)) input_path=trim(input_name)

    ! Initialize flag_count
    flag_count=0
	! Default and initial values.
	nstep_start=1
	num4category=2
	seed_switch=0
    sheetmove_flag=1
	sheet_switch=0.6
	ESURF_flag=0
	CONROT_flag=0
	CONROT_angle_limit=2.0
	CONROT_step=0.1
	CONROT_closure_tol=0.10
	CONROT_rmsd_limit=0.5
	displacement_factor_x=0.06
	displacement_factor_y=0.06
	displacement_factor_z=0.06

	filename=''
	expected_gnum=-1
	expected_chains=-1
	recalcu_switch=0
	idum=0
	nstep_terminal=-1
	ekt_seq=1.0
    ekt_seq_high=2.0
    ekt_seq_low=0.5
	ekt_sheet=0.6
    ekt_sheet_high=1.2
    ekt_sheet_low=0.3
	rmsd_max=3.0
	propensity_weighting_factor=2.0
	anneal_stages=0                   ! =0 no anneal, >=1 number of anneal stages
	steps_before_sheetmove=500
	sheetmove_interval=200
	! Sentinels used until analyzepdb selects automatic or custom mode.
	fpho_min=-1; fpho_max=-1
	fpol_min=-1; fpol_max=-1
	fchg_min=-1; fchg_max=-1
	foth_min=-1; foth_max=-1
	composition_min_set=.false.
	composition_max_set=.false.
    
    
	single_site_value=''
	grouped_pool_value=''
	grouped_site_value=''
	void_site_num=0
	void_site=0
	void_site_fixed_name=''
	group_value=''
	sheet_value=''
	sheet_defined=.false.
	max_user_sheet=0
	user_groups_used=.false.
	user_sheets_used=.false.

	! Specific-amino-acid count limits are appended in input-file order.
	n_morethan=0
	n_restrictions=0
	aa_min_defined=.false.
	aa_max_defined=.false.

	open(newunit=unit,file=trim(input_path),status='old',action='read', &
		iostat=ios,iomsg=iomsg)
	if (ios /= 0) then
		write(*,'(2A)') 'ERROR: cannot open input file: ',trim(iomsg)
		stop
	endif

	! Read each FLAG = value line. Flag position does not matter.
	line_number=0
	do
		read(unit,'(A)',iostat=ios) line
		if (ios < 0) exit
		line_number=line_number+1
		if (ios > 0) then
			write(*,'(A,I0)') 'ERROR: cannot read input line ',line_number
			stop
		endif

		call strip_comment(line)
		if (len_trim(line) == 0) cycle
		do i=1,len_trim(line)
			if (line(i:i) == achar(9)) line(i:i)=' '
		enddo

		equals_at=index(line,'=')
		if (equals_at <= 1) then
			write(*,'(A,I0)') 'ERROR: missing = on input line ',line_number
			stop
		endif

		key=adjustl(line(:equals_at-1))
		value=adjustl(line(equals_at+1:))
		call upcase(key)
		if (trim(key) == 'N_RESIDUE') key='N_RESIDUES'
		ios=0

		! All specific-amino-acid restrictions have the form N_<AA>_MIN/MAX.
		! Keeping them here avoids forty nearly identical CASE blocks and makes
		! duplicate MIN or MAX flags straightforward to reject.
		amino_acid_flag=.false.
		if (len_trim(key) == 9 .and. key(1:2) == 'N_' .and. &
			(key(6:9) == '_MIN' .or. key(6:9) == '_MAX')) then
			aa_index=0
			do i=1,size(amino_acid_names)
				if (key(3:5) == amino_acid_names(i)) then
					aa_index=i
					exit
				endif
			enddo

			if (aa_index > 0) then
				amino_acid_flag=.true.
				if (key(6:9) == '_MIN') then
					if (aa_min_defined(aa_index)) &
						call input_error('duplicate flag '//trim(key),line_number)
					if (n_morethan >= maximum_aa_restrictions) &
						call input_error('too many amino-acid minimum restrictions',line_number)
					n_morethan=n_morethan+1
					read(value,*,iostat=ios) aa_morethan(n_morethan)%min
					aa_morethan(n_morethan)%gtype=amino_acid_names(aa_index)
					aa_min_defined(aa_index)=.true.
				else
					if (aa_max_defined(aa_index)) &
						call input_error('duplicate flag '//trim(key),line_number)
					if (n_restrictions >= maximum_aa_restrictions) &
						call input_error('too many amino-acid maximum restrictions',line_number)
					n_restrictions=n_restrictions+1
					read(value,*,iostat=ios) aa_restrictions(n_restrictions)%max
					aa_restrictions(n_restrictions)%gtype=amino_acid_names(aa_index)
					aa_max_defined(aa_index)=.true.
				endif
				if (ios /= 0) call input_error('invalid value for '//trim(key),line_number)
			endif
		endif
		if (amino_acid_flag) cycle

		! Every ordinary input flag may appear only once.  The specific-amino-acid
		! flags are counted separately above, and SHEET_n flags are checked below.
		do counted_flag_index=1,number_of_counted_flags
			if (trim(key) == trim(counted_flag_names(counted_flag_index))) then
				flag_count(counted_flag_index)=flag_count(counted_flag_index)+1
				if (flag_count(counted_flag_index) > 1) &
					call input_error('duplicate flag '//trim(key),line_number)
				exit
			endif
		enddo

		! SHEET_1, SHEET_2, ... are optional user-defined beta sheets.
		if (index(trim(key),'SHEET_') == 1) then
			read(key(7:),*,iostat=ios) sheet_number
			if (ios /= 0 .or. sheet_number < 1 .or. &
				sheet_number > maximum_user_sheets) then
				write(*,'(A,I0,A)') 'ERROR on input line ',line_number, &
					': invalid SHEET number.'
				stop
			endif
			if (sheet_defined(sheet_number)) &
				call input_error('duplicate flag '//trim(key),line_number)
			sheet_value(sheet_number)=value
			sheet_defined(sheet_number)=.true.
			max_user_sheet=max(max_user_sheet,sheet_number)
			cycle
		endif

		select case(trim(key))
		case('PDBFILE')
			read(value,*,iostat=ios) filename
		case('N_RESIDUES')
			read(value,*,iostat=ios) expected_gnum
            
		case('N_CHAINS')
			read(value,*,iostat=ios) expected_chains
		case('GROUP_1')
			group_value(1)=value
		case('GROUP_2')
			group_value(2)=value
            
		case('RESTART')
			read(value,*,iostat=ios) recalcu_switch
		case('RANSEED')
			read(value,*,iostat=ios) idum
			seed_switch=1
		case('N_STEPS')
			read(value,*,iostat=ios) nstep_terminal
		case('KBT_SEQ')
			read(value,*,iostat=ios) ekt_seq
		case('KBT_SEQ_HIGH')
			read(value,*,iostat=ios) ekt_seq_high            
		case('KBT_SEQ_LOW')
			read(value,*,iostat=ios) ekt_seq_low          
		case('KBT_SHEETMOVE')
			read(value,*,iostat=ios) ekt_sheet
		case('KBT_SHEETMOVE_HIGH')
			read(value,*,iostat=ios) ekt_sheet_high
		case('KBT_SHEETMOVE_LOW')
			read(value,*,iostat=ios) ekt_sheet_low         
		case('RMSD_MAX')
			read(value,*,iostat=ios) rmsd_max
		case('PAGG_WEIGHT')
			read(value,*,iostat=ios) propensity_weighting_factor
		case('SHEETMOVE')
			read(value,*,iostat=ios) sheetmove_flag           
		case('ANNEAL_STAGES')
			read(value,*,iostat=ios) anneal_stages
		case('ESURF_MODE')
			read(value,*,iostat=ios) ESURF_flag
		case('STEPS_BEFORE_SHEETMOVE')
			read(value,*,iostat=ios) steps_before_sheetmove
		case('STEPS_BETWEEN_SHEETMOVE')
			read(value,*,iostat=ios) sheetmove_interval

		!case('N_HYDROPHOBIC')
		!	read(value,*,iostat=ios) fpho
		!case('N_POLAR')
		!	read(value,*,iostat=ios) fpol
		!case('N_CHARGED')
		!	read(value,*,iostat=ios) fchg
		!case('N_OTHER')
		!	read(value,*,iostat=ios) foth
		case('N_HYDROPHOBIC_MIN')
			read(value,*,iostat=ios) fpho_min
			composition_min_set(1)=.true.
		case('N_HYDROPHOBIC_MAX')
			read(value,*,iostat=ios) fpho_max
			composition_max_set(1)=.true.
		case('N_POLAR_MIN')
			read(value,*,iostat=ios) fpol_min
			composition_min_set(2)=.true.
		case('N_POLAR_MAX')
			read(value,*,iostat=ios) fpol_max
			composition_max_set(2)=.true.
		case('N_CHARGED_MIN')
			read(value,*,iostat=ios) fchg_min
			composition_min_set(3)=.true.
		case('N_CHARGED_MAX')
			read(value,*,iostat=ios) fchg_max
			composition_max_set(3)=.true.
		case('N_OTHER_MIN')
			read(value,*,iostat=ios) foth_min
			composition_min_set(4)=.true.
		case('N_OTHER_MAX')
			read(value,*,iostat=ios) foth_max
			composition_max_set(4)=.true.

		case('SINGLE_SITE_CONSTRAINTS')
			single_site_value=value
		case('GROUPED_SITE_AA_POOL')
			grouped_pool_value=value
		case('GROUPED_SITE_CONSTRAINTS')
			grouped_site_value=value

		case('CONROT_FLAG')
			read(value,*,iostat=ios) CONROT_flag
		case('CONROT_ANGLE_LIMIT')
			read(value,*,iostat=ios) CONROT_angle_limit
		case('CONROT_STEP')
			read(value,*,iostat=ios) CONROT_step
		case('CONROT_CLOSURE_TOL')
			read(value,*,iostat=ios) CONROT_closure_tol
		case('CONROT_RMSD_LIMIT')
			read(value,*,iostat=ios) CONROT_rmsd_limit

		case default
			write(*,'(A,I0,2A)') 'ERROR on input line ',line_number, &
				': unknown flag ',trim(key)
			stop
		end select

		if (ios /= 0) then
			write(*,'(A,I0,2A)') 'ERROR on input line ',line_number, &
				': invalid value for ',trim(key)
			stop
		endif
	enddo
	close(unit)

	! Generate the clock seed now so the displayed value is the seed that
	! init_random_seed will actually use later.
	if (seed_switch == 0) then
		call system_clock(count=clock_count)
		idum=modulo(clock_count,1000000000)
		if (idum == 0) idum=123457
	endif

	! Check the required values.
	if (len_trim(filename) == 0) then
		write(*,*) 'ERROR: PDBFILE is required.'
		stop
	endif
	if (nstep_terminal < 1) then
		write(*,*) 'ERROR: N_STEPS is required and must be positive.'
		stop
	endif
	if (anneal_stages < 0 .or. rmsd_max < 0.0 .or. &
		propensity_weighting_factor < 0.0) then
		write(*,*) 'ERROR: invalid MC input.'
		stop
	endif
	if (anneal_stages > nstep_terminal) then
		write(*,*) 'ERROR: ANNEAL_STAGES cannot exceed N_STEPS.'
		stop
	endif
	if (recalcu_switch /= 0 .and. recalcu_switch /= 1) then
		write(*,*) 'ERROR: RESTART must be 0 or 1.'
		stop
	endif
	if (sheetmove_flag /= 0 .and. sheetmove_flag /= 1) then
		write(*,*) 'ERROR: SHEETMOVE must be 0 or 1.'
		stop
	endif
	if (sheetmove_flag == 0) sheet_switch=0.0
	if (ESURF_flag < 0 .or. ESURF_flag > 2) then
		write(*,*) 'ERROR: ESURF_MODE must be 0, 1, or 2.'
		stop
	endif
	if (steps_before_sheetmove < 0) then
		write(*,*) 'ERROR: STEPS_BEFORE_SHEETMOVE must not be negative.'
		stop
	endif
	if (sheetmove_interval <= 0) then
		write(*,*) 'ERROR: STEPS_BETWEEN_SHEETMOVE must be positive.'
		stop
	endif
	if (CONROT_flag /= 0 .and. CONROT_flag /= 1) then
		write(*,*) 'ERROR: CONROT_FLAG must be 0 or 1.'
		stop
	endif
	if (CONROT_angle_limit <= 0.0 .or. CONROT_step <= 0.0 .or. &
		CONROT_closure_tol <= 0.0 .or. CONROT_rmsd_limit < 0.0) then
		write(*,*) 'ERROR: invalid CONROT input.'
		stop
	endif
	if (anneal_stages >= 1) then
		if (ekt_seq_low <= 0.0 .or. ekt_seq_high <= 0.0 .or. &
			ekt_sheet_low <= 0.0 .or. ekt_sheet_high <= 0.0) then
			write(*,*) 'ERROR: annealing kBT values must be positive.'
			stop
		endif
		if (ekt_seq_high < ekt_seq_low .or. &
			ekt_sheet_high < ekt_sheet_low) then
			write(*,*) 'ERROR: annealing kBT high must not be lower than kBT low.'
			stop
		endif
	else
		if (ekt_seq <= 0.0 .or. ekt_sheet <= 0.0) then
			write(*,*) 'ERROR: kBT values must be positive.'
			stop
		endif
	endif

	! N_RESIDUE and N_CHAINS are optional, but they must be used together.
	if (expected_gnum == -1 .and. expected_chains == -1) then
		call analyzepdb
	elseif (expected_gnum > 0 .and. expected_chains > 0) then
		call analyzepdb(expected_gnum,expected_chains)
	else
		write(*,*) 'ERROR: enter both N_RESIDUE and N_CHAINS, or omit both.'
		stop
    endif

	! Determine the hydrogen-bond network before assigning peptide groups or sheets.
	call Hbond_matrix

	! Optional user-defined energy groups.  If omitted, keep the DSSP groups.
	if ((len_trim(group_value(1)) == 0) .neqv. &
		(len_trim(group_value(2)) == 0)) then
		write(*,*) 'ERROR: define both GROUP_1 and GROUP_2, or omit both.'
		stop
	endif

	allocate(chain_counts(repeated_unit))
	chain_counts=0
	if (len_trim(group_value(1)) > 0) then
		do group_index=1,2
			line=group_value(group_index)
			do k=1,len_trim(line)
				if (line(k:k) == ',') line(k:k)=' '
			enddo
			token_count=0
			in_token=.false.
			do k=1,len_trim(line)
				if (line(k:k) == ' ') then
					in_token=.false.
				elseif (.not.in_token) then
					token_count=token_count+1
					in_token=.true.
				endif
			enddo
			if (token_count < 1 .or. token_count > repeated_unit) then
				write(*,'(A,I0,A)') 'ERROR: GROUP_',group_index, &
					' has an invalid number of chain IDs.'
				stop
			endif
			selfassembly(group_index)%num4peptides=token_count
			selfassembly(group_index)%peptideID=0
			read(line,*,iostat=ios) &
				(selfassembly(group_index)%peptideID(i),i=1,token_count)
			if (ios /= 0) then
				write(*,'(A,I0,A)') 'ERROR: GROUP_',group_index, &
					' contains an invalid chain ID.'
				stop
			endif
			do i=1,token_count
				chain_id=selfassembly(group_index)%peptideID(i)
				if (chain_id < 1 .or. chain_id > repeated_unit) then
					write(*,'(A,I0,A)') 'ERROR: GROUP_',group_index, &
						' contains a chain ID outside the PDB range.'
					stop
				endif
				chain_counts(chain_id)=chain_counts(chain_id)+1
			enddo
		enddo
		if (any(chain_counts /= 1)) then
			write(*,*) 'ERROR: GROUP_1 and GROUP_2 must contain every chain exactly once.'
			stop
		endif
		user_groups_used=.true.
	endif

	! Optional user-defined beta sheets.  SHEET numbers must be consecutive,
	! and the complete set must contain every chain exactly once.
	if (max_user_sheet > 0) then
		if (max_user_sheet > repeated_unit) then
			write(*,*) 'ERROR: number of user-defined sheets exceeds number of chains.'
			stop
		endif
		do sheet_number=1,max_user_sheet
			if (.not.sheet_defined(sheet_number)) then
				write(*,'(A,I0,A)') 'ERROR: SHEET_',sheet_number,' is missing.'
				stop
			endif
		enddo

		if (allocated(betasheets)) deallocate(betasheets)
		allocate(betasheets(max_user_sheet))
		chain_counts=0
		do sheet_number=1,max_user_sheet
			allocate(betasheets(sheet_number)%peptideID(repeated_unit))
			betasheets(sheet_number)%peptideID=0
			line=sheet_value(sheet_number)
			do k=1,len_trim(line)
				if (line(k:k) == ',') line(k:k)=' '
			enddo
			token_count=0
			in_token=.false.
			do k=1,len_trim(line)
				if (line(k:k) == ' ') then
					in_token=.false.
				elseif (.not.in_token) then
					token_count=token_count+1
					in_token=.true.
				endif
			enddo
			if (token_count < 1 .or. token_count > repeated_unit) then
				write(*,'(A,I0,A)') 'ERROR: SHEET_',sheet_number, &
					' has an invalid number of chain IDs.'
				stop
			endif
			betasheets(sheet_number)%num4peptides=token_count
			read(line,*,iostat=ios) &
				(betasheets(sheet_number)%peptideID(i),i=1,token_count)
			if (ios /= 0) then
				write(*,'(A,I0,A)') 'ERROR: SHEET_',sheet_number, &
					' contains an invalid chain ID.'
				stop
			endif
			do i=1,token_count
				chain_id=betasheets(sheet_number)%peptideID(i)
				if (chain_id < 1 .or. chain_id > repeated_unit) then
					write(*,'(A,I0,A)') 'ERROR: SHEET_',sheet_number, &
						' contains a chain ID outside the PDB range.'
					stop
				endif
				chain_counts(chain_id)=chain_counts(chain_id)+1
			enddo
		enddo
		if (any(chain_counts /= 1)) then
			write(*,*) 'ERROR: user-defined sheets must contain every chain exactly once.'
			stop
		endif
		num4betasheets=max_user_sheet
		user_sheets_used=.true.
	endif
	deallocate(chain_counts)

	! These three values may be blank or NONE.
	line=single_site_value
	call upcase(line)
	if (trim(adjustl(line)) == 'NONE') single_site_value=''
	line=grouped_pool_value
	call upcase(line)
	if (trim(adjustl(line)) == 'NONE') then
		grouped_pool_value=''
	else
		grouped_pool_value=line
	endif
	line=grouped_site_value
	call upcase(line)
	if (trim(adjustl(line)) == 'NONE') grouped_site_value=''

	call read_void_site_input(trim(single_site_value))
	call read_restraints_input(trim(grouped_pool_value),3)
	call read_nmr_site_input(trim(grouped_site_value),' ')

	! analyzepdb has already placed ACE/NME/NHE sites in void_site.  User-defined
	! single sites are appended to the same list.  A grouped site cannot overlap
	! any entry in that combined void-site list.
	do i=1,void_site_num
		do j=i+1,void_site_num
			if (void_site(i) == void_site(j)) &
				call input_error('duplicate single-site constraint or cap site',0)
		enddo
	enddo
	do i=1,nmr_site_num
		if (void_site_num > 0) then
			if (any(void_site(1:void_site_num) == nmr_site_ID(i))) &
				call input_error('a grouped site overlaps a cap or single-site constraint',0)
		endif
	enddo
	if (nmr_site_num > 0 .and. NMR_pool_size < nmr_site_num) &
		call input_error('grouped-site AA pool is smaller than the grouped-site list',0)

	expected_mutation_sites=gnum-void_site_num-nmr_site_num
	if (expected_mutation_sites < 0) &
		call input_error('positional constraints exceed the number of amino acids',0)

	! Count amino acids fixed by user single-site constraints.  A blank fixed
	! name means retain the residue from the initial PDB structure.
	fixed_aa_counts=0
	cap_site_count=gnum-amino_acid_residue_count
	do i=cap_site_count+1,void_site_num
		fixed_residue_name=void_site_fixed_name(i)
		if (len_trim(fixed_residue_name) == 0) &
			fixed_residue_name=initial_residue_name(void_site(i))
		fixed_residue_name=adjustl(fixed_residue_name)
		call upcase(fixed_residue_name)
		if (len_trim(fixed_residue_name) == 4) then
			if (fixed_residue_name(1:1) == 'N' .or. &
				fixed_residue_name(1:1) == 'C') &
				fixed_residue_name=fixed_residue_name(2:4)
		endif
		if (trim(fixed_residue_name) == 'HIS' .or. &
			trim(fixed_residue_name) == 'HID' .or. &
			trim(fixed_residue_name) == 'HIP') fixed_residue_name='HIE'

		aa_index=0
		do j=1,size(amino_acid_names)
			if (trim(fixed_residue_name) == amino_acid_names(j)) then
				aa_index=j
				exit
			endif
		enddo
		if (aa_index == 0) &
			call input_error('unsupported fixed residue '//trim(fixed_residue_name),0)
		fixed_aa_counts(aa_index)=fixed_aa_counts(aa_index)+1
	enddo

	! Classify fixed single sites and every grouped-pool entry.  Pool entries
	! must be unique because grouped sites select without replacement.
	fixed_category_counts=0
	do aa_index=1,size(amino_acid_names)
		fixed_category_counts(amino_acid_categories(aa_index))= &
			fixed_category_counts(amino_acid_categories(aa_index))+ &
			fixed_aa_counts(aa_index)
	enddo

	grouped_pool_category_counts=0
	do i=1,NMR_pool_size
		aa_index=0
		do j=1,size(amino_acid_names)
			if (trim(NMR_AA_pool(i)) == amino_acid_names(j)) then
				aa_index=j
				exit
			endif
		enddo
		if (aa_index == 0) &
			call input_error('unsupported grouped-pool residue '//trim(NMR_AA_pool(i)),0)
		grouped_pool_category_counts(amino_acid_categories(aa_index))= &
			grouped_pool_category_counts(amino_acid_categories(aa_index))+1
	enddo

	do i=1,4
		grouped_category_minimums(i)=max(0,nmr_site_num- &
			(NMR_pool_size-grouped_pool_category_counts(i)))
		grouped_category_maximums(i)=min(nmr_site_num, &
			grouped_pool_category_counts(i))
	enddo
	positional_category_minimums=fixed_category_counts+grouped_category_minimums
	positional_category_maximums=fixed_category_counts+ &
		grouped_category_maximums+expected_mutation_sites

	! Composition range and feasibility checks.
	if (fpho_min < 0 .or. fpho_min > gnum .or. fpho_max < 0 .or. fpho_max > gnum .or. &
		fpol_min < 0 .or. fpol_min > gnum .or. fpol_max < 0 .or. fpol_max > gnum .or. &
		fchg_min < 0 .or. fchg_min > gnum .or. fchg_max < 0 .or. fchg_max > gnum .or. &
		foth_min < 0 .or. foth_min > gnum .or. foth_max < 0 .or. foth_max > gnum) then
		write(*,*) 'ERROR: invalid composition range.'
		stop
	endif
	if (fpho_min > fpho_max) then
		write(*,*) 'ERROR: N_HYDROPHOBIC_MIN is larger than N_HYDROPHOBIC_MAX.'
		stop
	endif
	if (fpol_min > fpol_max) then
		write(*,*) 'ERROR: N_POLAR_MIN is larger than N_POLAR_MAX.'
		stop
	endif
	if (fchg_min > fchg_max) then
		write(*,*) 'ERROR: N_CHARGED_MIN is larger than N_CHARGED_MAX.'
		stop
	endif
	if (foth_min > foth_max) then
		write(*,*) 'ERROR: N_OTHER_MIN is larger than N_OTHER_MAX.'
		stop
	endif
	if (positional_category_minimums(1) > fpho_max) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: single/grouped constraints require at least ', &
			positional_category_minimums(1), &
			' hydrophobic amino acids, exceeding N_HYDROPHOBIC_MAX (',fpho_max,').'
		stop
	endif
	if (positional_category_minimums(2) > fpol_max) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: single/grouped constraints require at least ', &
			positional_category_minimums(2), &
			' polar amino acids, exceeding N_POLAR_MAX (',fpol_max,').'
		stop
	endif
	if (positional_category_minimums(3) > fchg_max) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: single/grouped constraints require at least ', &
			positional_category_minimums(3), &
			' charged amino acids, exceeding N_CHARGED_MAX (',fchg_max,').'
		stop
	endif
	if (positional_category_minimums(4) > foth_max) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: single/grouped constraints require at least ', &
			positional_category_minimums(4), &
			' other-category amino acids, exceeding N_OTHER_MAX (',foth_max,').'
		stop
	endif

	if (positional_category_maximums(1) < fpho_min) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: positional constraints allow at most ', &
			positional_category_maximums(1), &
			' hydrophobic amino acids, below N_HYDROPHOBIC_MIN (',fpho_min,').'
		stop
	endif
	if (positional_category_maximums(2) < fpol_min) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: positional constraints allow at most ', &
			positional_category_maximums(2), &
			' polar amino acids, below N_POLAR_MIN (',fpol_min,').'
		stop
	endif
	if (positional_category_maximums(3) < fchg_min) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: positional constraints allow at most ', &
			positional_category_maximums(3), &
			' charged amino acids, below N_CHARGED_MIN (',fchg_min,').'
		stop
	endif
	if (positional_category_maximums(4) < foth_min) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: positional constraints allow at most ', &
			positional_category_maximums(4), &
			' other-category amino acids, below N_OTHER_MIN (',foth_min,').'
		stop
	endif
	if (fpho_min+fpol_min+fchg_min+foth_min > amino_acid_residue_count) then
		write(*,*) 'ERROR: composition minimums exceed the number of amino acids.'
		stop
	endif
	if (fpho_max+fpol_max+fchg_max+foth_max < amino_acid_residue_count) then
		write(*,*) 'ERROR: composition maximums is less than the number of amino acids.'
		stop
	endif

	! Lower and upper restrictions are independent arrays. Match them by name.
	do i=1,n_morethan
		if (aa_morethan(i)%min < 0 .or. aa_morethan(i)%min > gnum) then
			write(*,*) 'ERROR: invalid minimum for ',aa_morethan(i)%gtype
			stop
		endif
		aa_morethan(i)%max=gnum
	enddo

	do j=1,n_restrictions
		if (aa_restrictions(j)%max < 0 .or. aa_restrictions(j)%max > gnum) then
			write(*,*) 'ERROR: invalid maximum for ',aa_restrictions(j)%gtype
			stop
		endif
		aa_restrictions(j)%min=0
	enddo

	do i=1,n_morethan
		do j=1,n_restrictions
			if (aa_morethan(i)%gtype == aa_restrictions(j)%gtype) then
				if (aa_morethan(i)%min > aa_restrictions(j)%max) then
					write(*,*) 'ERROR: minimum is larger than maximum for ', &
						aa_morethan(i)%gtype
					stop
				endif
			endif
		enddo
	enddo

	! Fixed single-site residues contribute to the total sequence composition.
	! For each amino acid, the effective minimum is the larger of its user MIN
	! and the number already forced by fixed sites.
	user_aa_minimums=0
	user_aa_maximums=amino_acid_residue_count
	do i=1,n_morethan
		do aa_index=1,size(amino_acid_names)
			if (aa_morethan(i)%gtype == amino_acid_names(aa_index)) then
				user_aa_minimums(aa_index)=aa_morethan(i)%min
				exit
			endif
		enddo
	enddo
	do i=1,n_restrictions
		do aa_index=1,size(amino_acid_names)
			if (aa_restrictions(i)%gtype == amino_acid_names(aa_index)) then
				user_aa_maximums(aa_index)=aa_restrictions(i)%max
				exit
			endif
		enddo
	enddo

	do aa_index=1,size(amino_acid_names)
		if (fixed_aa_counts(aa_index) > user_aa_maximums(aa_index)) then
			write(*,'(3A,I0,A,I0,A)') 'ERROR: fixed ',amino_acid_names(aa_index), &
				' count (',fixed_aa_counts(aa_index),') exceeds its maximum (', &
				user_aa_maximums(aa_index),').' 
			stop
		endif
	enddo

	effective_aa_minimums=max(user_aa_minimums,fixed_aa_counts)
	specific_minimums=0
	do aa_index=1,size(amino_acid_names)
		specific_minimums(amino_acid_categories(aa_index))= &
			specific_minimums(amino_acid_categories(aa_index))+ &
			effective_aa_minimums(aa_index)
	enddo
	if (specific_minimums(1) > fpho_max) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: minimum number of required hydrophobic amino acids (', &
			specific_minimums(1),') exceeds N_HYDROPHOBIC_MAX (',fpho_max,').'
		stop
	endif
	if (specific_minimums(2) > fpol_max) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: minimum number of required polar amino acids (', &
			specific_minimums(2),') exceeds N_POLAR_MAX (',fpol_max,').'
		stop
	endif
	if (specific_minimums(3) > fchg_max) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: minimum number of required charged amino acids (', &
			specific_minimums(3),') exceeds N_CHARGED_MAX (',fchg_max,').'
		stop
	endif
	if (specific_minimums(4) > foth_max) then
		write(*,'(A,I0,A,I0,A)') 'ERROR: minimum number of required other-category amino acids (', &
			specific_minimums(4),') exceeds N_OTHER_MAX (',foth_max,').'
		stop
	endif

	! A category minimum and the individual minima inside that category can
	! impose different requirements.  Their combined requirements must still
	! fit inside the peptide.
	category_minimum_requirements=(/fpho_min,fpol_min,fchg_min,foth_min/)
	category_minimum_requirements=max(category_minimum_requirements,specific_minimums)
	category_minimum_requirements=max(category_minimum_requirements, &
		positional_category_minimums)
	if (sum(category_minimum_requirements) > amino_acid_residue_count) then
		write(*,*) 'ERROR: category and specific amino-acid minima cannot fit in the peptide.'
		stop
	endif

	! Treat an omitted specific maximum as the full peptide size.  The values
	! 9, 5, 4, and 2 are the numbers of amino-acid types in each category.
	specific_maximum_capacity=(/9,5,4,2/)*amino_acid_residue_count
	do j=1,n_restrictions
		select case(trim(aa_restrictions(j)%gtype))
		case('GLY','ALA','VAL','LEU','ILE','MET','PHE','TYR','TRP')
			specific_maximum_capacity(1)=specific_maximum_capacity(1)- &
				amino_acid_residue_count+min(aa_restrictions(j)%max,amino_acid_residue_count)
		case('SER','ASN','GLN','THR','HIE')
			specific_maximum_capacity(2)=specific_maximum_capacity(2)- &
				amino_acid_residue_count+min(aa_restrictions(j)%max,amino_acid_residue_count)
		case('ARG','LYS','GLU','ASP')
			specific_maximum_capacity(3)=specific_maximum_capacity(3)- &
				amino_acid_residue_count+min(aa_restrictions(j)%max,amino_acid_residue_count)
		case('CYS','PRO')
			specific_maximum_capacity(4)=specific_maximum_capacity(4)- &
				amino_acid_residue_count+min(aa_restrictions(j)%max,amino_acid_residue_count)
		end select
	enddo
	if (specific_maximum_capacity(1) < fpho_min) then
		write(*,*) 'ERROR: specific hydrophobic maxima cannot satisfy N_HYDROPHOBIC_MIN.'
		stop
	endif
	if (specific_maximum_capacity(2) < fpol_min) then
		write(*,*) 'ERROR: specific polar maxima cannot satisfy N_POLAR_MIN.'
		stop
	endif
	if (specific_maximum_capacity(3) < fchg_min) then
		write(*,*) 'ERROR: specific charged maxima cannot satisfy N_CHARGED_MIN.'
		stop
	endif
	if (specific_maximum_capacity(4) < foth_min) then
		write(*,*) 'ERROR: specific other-residue maxima cannot satisfy N_OTHER_MIN.'
		stop
	endif

	category_maximum_capacity=(/fpho_max,fpol_max,fchg_max,foth_max/)
	category_maximum_capacity=min(category_maximum_capacity,specific_maximum_capacity)
	category_maximum_capacity=min(category_maximum_capacity, &
		positional_category_maximums)
	if (sum(category_maximum_capacity) < amino_acid_residue_count) then
		write(*,*) 'ERROR: category and specific amino-acid maxima cannot fill the peptide.'
		stop
	endif

	!sheetmove_interval=steps_between_sheetmove

	! Print an aligned summary of the values that will be used.
	write(6,*)
	write(6,'(A)') '========================================================================'
	write(6,'(A)') '                          PepAD input summary'
	write(6,'(A)') '========================================================================'

	write(6,*)
	write(6,'(A)') '[ Initial structure ]'
	write(6,'(1X,A,T38,"= ",A)') 'PDB file name',trim(filename)
	write(6,'(1X,A,T38,"= ",I0)') 'Chain length (including caps)',gnum
	write(6,'(1X,A,T38,"= ",I0)') 'Amino acids (excluding caps)', &
		amino_acid_residue_count
	write(6,'(1X,A,T38,"= ",I0)') 'Total chains',repeated_unit
	if (user_sheets_used) then
		write(6,'(1X,A,T38,"= ",A)') 'Beta-sheet assignment','user-defined'
	else
		write(6,'(1X,A,T38,"= ",A)') 'Beta-sheet assignment','DSSP calculated'
	endif
	if (user_groups_used) then
		write(6,'(1X,A,T38,"= ",A)') 'Energy-group assignment','user-defined'
	else
		write(6,'(1X,A,T38,"= ",A)') 'Energy-group assignment','DSSP calculated'
	endif
	write(6,'(1X,A,T38,"= ",I0)') 'Number of beta sheets',num4betasheets
	do j=1,num4betasheets
		write(sheet_label,'("Beta sheet ",I0," chain IDs")') j
		write(6,'(1X,A,T38,"= ",*(I0,1X))') trim(sheet_label), &
			(betasheets(j)%peptideID(i),i=1,betasheets(j)%num4peptides)
	enddo
	write(6,'(1X,A,T38,"= ",*(I0,1X))') 'Energy group 1 chain IDs', &
		(selfassembly(1)%peptideID(i),i=1,selfassembly(1)%num4peptides)
	write(6,'(1X,A,T38,"= ",*(I0,1X))') 'Energy group 2 chain IDs', &
		(selfassembly(2)%peptideID(i),i=1,selfassembly(2)%num4peptides)

	write(6,*)
	write(6,'(A)') '[ Monte Carlo and energy ]'
	write(6,'(1X,A,T38,"= ",I0)') 'Restart mode',recalcu_switch
	if (seed_switch == 1) then
		write(6,'(1X,A,T38,"= ",I0)') 'Random seed (user-defined)',idum
	else
		write(6,'(1X,A,T38,"= ",I0)') 'Random seed (system clock)',idum
	endif
	write(6,'(1X,A,T38,"= ",I0)') 'Monte Carlo steps',nstep_terminal
	if (anneal_stages >= 1) then
		write(6,'(1X,A,T38,"= ",F5.3," to ",F5.3)') &
			'kBT (sequence), high to low',ekt_seq_high,ekt_seq_low
		write(6,'(1X,A,T38,"= ",F5.3," to ",F5.3)') &
			'kBT (sheet move), high to low',ekt_sheet_high,ekt_sheet_low
	else
		write(6,'(1X,A,T38,"= ",F5.3)') 'kBT (sequence)',ekt_seq
		write(6,'(1X,A,T38,"= ",F5.3)') 'kBT (sheet move)',ekt_sheet
	endif
	write(6,'(1X,A,T38,"= ",F5.3)') 'Maximum RMSD (angstrom)',rmsd_max
	write(6,'(1X,A,T38,"= ",F5.3)') 'Aggregation propensity weight', &
		propensity_weighting_factor
	write(6,'(1X,A,T38,"= ",I0)') 'Annealing stages',anneal_stages
    write(6,'(1X,A,T38,"= ",I0)') 'Sheet-move flag',sheetmove_flag
    if (sheetmove_flag == 1) then
		write(6,'(1X,A,T38,"= ",A)') 'Sheet position perturbation move','enabled'
		!write(6,'(1X,A,T38,"= ",F5.3)') 'Sheet-move proposal probability',sheet_switch
		write(6,'(1X,A,T38,"= ",I0)') 'Steps before sheet move',steps_before_sheetmove
		write(6,'(1X,A,T38,"= ",I0)') 'Steps between sheet moves',sheetmove_interval
    elseif (sheetmove_flag == 0) then 
		write(6,'(1X,A,T38,"= ",A)') 'Sheet position perturbation move','disabled'
    endif
    
    write(6,'(1X,A,T38,"= ",I0)') 'ESURF mode',ESURF_flag
    if (ESURF_flag == 0) then
		write(6,'(1X,A,T38,"= ",A)') 'Non-polar solvation energy','disabled'
    elseif (ESURF_flag == 1) then
		write(6,'(1X,A,T38,"= ",A)') 'Non-polar solvation energy','full calculation'
    elseif (ESURF_flag == 2) then
		write(6,'(1X,A,T38,"= ",A)') 'Non-polar solvation energy','neighbor-list calculation'
	endif
	write(6,*)
	write(6,'(A)') '[ Compositions ]'
	write(6,'(1X,A,T38,"= ",I0," / ",I0)') 'Hydrophobic min/max',fpho_min,fpho_max
	write(6,'(1X,A,T38,"= ",I0," / ",I0)') 'Polar min/max',fpol_min,fpol_max
	write(6,'(1X,A,T38,"= ",I0," / ",I0)') 'Charged min/max',fchg_min,fchg_max
	write(6,'(1X,A,T38,"= ",I0," / ",I0)') 'Other min/max',foth_min,foth_max

	write(6,*)
	write(6,'(A)') '[ Compositional constraints ]'
	if (n_morethan == 0 .and. n_restrictions == 0) then
		write(6,'(1X,A,T38,"= ",A)') 'Restrictions','none'
	else
		! Print every amino acid that has a lower or upper restriction.
		do i=1,n_morethan
			min_value=aa_morethan(i)%min
			max_value=gnum
			do j=1,n_restrictions
				if (aa_morethan(i)%gtype == aa_restrictions(j)%gtype) then
					max_value=aa_restrictions(j)%max
					exit
				endif
			enddo
			restriction_label=trim(aa_morethan(i)%gtype)//' min/max'
			write(6,'(1X,A,T38,"= ",I0," / ",I0)') &
				trim(restriction_label),min_value,max_value
		enddo

		! Print upper-only restrictions with an unrestricted minimum of zero.
		do j=1,n_restrictions
			matching_minimum=.false.
			do i=1,n_morethan
				if (aa_restrictions(j)%gtype == aa_morethan(i)%gtype) then
					matching_minimum=.true.
					exit
				endif
			enddo
			if (.not.matching_minimum) then
				restriction_label=trim(aa_restrictions(j)%gtype)//' min/max'
				write(6,'(1X,A,T38,"= ",I0," / ",I0)') &
					trim(restriction_label),0,aa_restrictions(j)%max
			endif
		enddo
	endif

	write(6,*)
	write(6,'(A)') '[ Positional constraints ]'
	write(6,'(1X,A,T38,"= ",I0)') 'Freely mutable amino-acid sites', &
		expected_mutation_sites
	write(6,'(1X,A,T38,"= ",I0)') 'Number of cap sites', &
		gnum-amino_acid_residue_count
	write(6,'(1X,A,T38,"= ",I0)') 'Number of single-site constraints', &
		void_site_num-(gnum-amino_acid_residue_count)
	if (void_site_num > 0) then
		write(6,'(1X,A,T38,"= ",*(I0,1X))') 'Void-site IDs (caps and single)', &
			(void_site(i),i=1,void_site_num)
	else
		write(6,'(1X,A,T38,"= ",A)') 'Void-site IDs (caps and single)','none'
	endif
	write(6,'(1X,A,T38,"= ",I0)') 'Number of grouped-site constraints',nmr_site_num
	if (nmr_site_num > 0) then
		write(6,'(1X,A,T38,"= ",*(I0,1X))') 'Grouped-site chain positions', &
			(nmr_site_ID(i),i=1,nmr_site_num)
	else
		write(6,'(1X,A,T38,"= ",A)') 'Grouped-site chain positions','none'
	endif
	if (NMR_pool_size > 0) then
		write(6,'(1X,A,T38,"= ",*(A,1X))') 'Grouped-site AA pool', &
			(NMR_AA_pool(i),i=1,NMR_pool_size)
	else
		write(6,'(1X,A,T38,"= ",A)') 'Grouped-site AA pool','none'
	endif
	if (nmr_site_num > 0) then
		write(6,'(1X,A,T38,"= ",I0," / ",I0)') 'Grouped hydrophobic min/max', &
			grouped_category_minimums(1),grouped_category_maximums(1)
		write(6,'(1X,A,T38,"= ",I0," / ",I0)') 'Grouped polar min/max', &
			grouped_category_minimums(2),grouped_category_maximums(2)
		write(6,'(1X,A,T38,"= ",I0," / ",I0)') 'Grouped charged min/max', &
			grouped_category_minimums(3),grouped_category_maximums(3)
		write(6,'(1X,A,T38,"= ",I0," / ",I0)') 'Grouped other min/max', &
			grouped_category_minimums(4),grouped_category_maximums(4)
	endif
	write(6,'(A)') '========================================================================'
	end subroutine inputfile

	subroutine input_error(message,line_number)
		implicit none
		character(len=*), intent(in) :: message
		integer, intent(in) :: line_number

		if (line_number > 0) then
			write(*,'(A,I0,2A)') 'ERROR in input line ',line_number,': ',trim(message)
		else
			write(*,'(2A)') 'ERROR: ',trim(message)
		endif
		stop 1
	end subroutine input_error



    subroutine strip_comment(s)
		character(len=*), intent(inout) :: s
		integer :: p
		p = index(s,'#'); if (p > 0) s = s(:p-1)
		p = index(s,'!'); if (p > 0) s = s(:p-1)
	end subroutine strip_comment

	subroutine upcase(s)
		character(len=*), intent(inout) :: s
		integer :: i, ia
		do i = 1, len_trim(s)
			ia = ichar(s(i:i))
			if (ia.ge.97.and.ia.le.122) s(i:i) = char(ia-32)
		enddo
	end subroutine upcase
    
    subroutine setparameter
    integer :: i
    
    atom_num=60*gnum
    allocate (sitenum4mutation_group(gnum))
    allocate (original_group(repeated_unit, gnum))
	sitenum4mutation=0
	sitenum4mutation_group=0

	hydrationprop%hgtype = "NONE"
	hydrationprop%pgtype = "NONE"
	hydrationprop%cgtype = "NONE"
	hydrationprop%ogtype = "NONE"

    hydrationprop%hnum=9          ! hydrophobic amino acid number
    hydrationprop%pnum=5          ! polar amino acid number
    hydrationprop%cnum=4          ! charged amino acid number
    hydrationprop%onum=2          ! other amino acid number
    ! standard amino acid
    hydrationprop%hgtype(1:9)=(/"ALA", "LEU", "VAL", "ILE", "MET", "PHE", "TYR", "TRP", "GLY"/)
    hydrationprop%pgtype(1:5)=(/"ASN", "GLN", "SER", "THR", "HIE"/)
    hydrationprop%cgtype(1:4)=(/"ARG", "LYS", "GLU", "ASP"/)
    hydrationprop%ogtype(1:2)=(/"PRO", "CYS"/)

    ! detect the amino acid types (hydrophobic, polar, charged) of the restricted ones
    do i = 1, n_restrictions
		do j = 1, hydrationprop%hnum
            if (hydrationprop%hgtype(j) == aa_restrictions(i)%gtype) then
                aa_restrictions(i)%clas = 1									! hydrophobic
                if (aa_restrictions(i)%max == 0) then      
					do k = j, hydrationprop%hnum - 1
						hydrationprop%hgtype(k) = hydrationprop%hgtype(k+1)
					end do
					hydrationprop%hgtype(hydrationprop%hnum) = "NONE"		! Clear last entry if the maximum number is 0
					hydrationprop%hnum = hydrationprop%hnum - 1				! Decrease count
                endif
                exit
            endif
        enddo
        
        do j = 1, hydrationprop%pnum
            if (hydrationprop%pgtype(j) == aa_restrictions(i)%gtype) then
                aa_restrictions(i)%clas = 2									! polar
                if (aa_restrictions(i)%max == 0) then      
					do k = j, hydrationprop%pnum - 1
						hydrationprop%pgtype(k) = hydrationprop%pgtype(k+1)
					end do
					hydrationprop%pgtype(hydrationprop%pnum) = "NONE"		! Clear last entry if the maximum number is 0
					hydrationprop%pnum = hydrationprop%pnum - 1				! Decrease count
                endif
                exit
            endif
        enddo
        
        do j = 1, hydrationprop%cnum
            if (hydrationprop%cgtype(j) == aa_restrictions(i)%gtype) then
                aa_restrictions(i)%clas = 3									! charged
                if (aa_restrictions(i)%max == 0) then      
					do k = j, hydrationprop%cnum - 1
						hydrationprop%cgtype(k) = hydrationprop%cgtype(k+1)
					end do
					hydrationprop%cgtype(hydrationprop%cnum) = "NONE"		! Clear last entry if the maximum number is 0
					hydrationprop%cnum = hydrationprop%cnum - 1				! Decrease count
                endif
                exit
            endif
        enddo
        
        do j = 1, hydrationprop%onum
            if (hydrationprop%ogtype(j) == aa_restrictions(i)%gtype) then
                aa_restrictions(i)%clas = 4									! other
                if (aa_restrictions(i)%max == 0) then      
					do k = j, hydrationprop%onum - 1
						hydrationprop%ogtype(k) = hydrationprop%ogtype(k+1)
					end do
					hydrationprop%ogtype(hydrationprop%onum) = "NONE"		! Clear last entry if the maximum number is 0
					hydrationprop%onum = hydrationprop%onum - 1				! Decrease count
                endif
                exit
            endif
        enddo
    enddo
    
    ! 2026/5/08 isotropic RMSD
    ! max_X : max_Y : max_Z = 1:1:1
    rmsd_max_x = sqrt(rmsd_max*rmsd_max/3.0)
    rmsd_max_y = sqrt(rmsd_max*rmsd_max/3.0)
    rmsd_max_z = sqrt(rmsd_max*rmsd_max/3.0)
    
		end subroutine setparameter

	integer function chain_site_from_amino_acid(amino_acid_site)
		implicit none
		integer, intent(in) :: amino_acid_site
		integer :: chain_site, amino_acid_count, cap_site_count

		chain_site_from_amino_acid=0
		amino_acid_count=0
		cap_site_count=gnum-amino_acid_residue_count

		! analyzepdb stores the cap positions first in void_site.  Count only
		! non-cap positions until the requested amino-acid number is reached.
		do chain_site=1,gnum
			if (cap_site_count > 0) then
				if (any(void_site(1:cap_site_count) == chain_site)) cycle
			endif
			amino_acid_count=amino_acid_count+1
			if (amino_acid_count == amino_acid_site) then
				chain_site_from_amino_acid=chain_site
				return
			endif
		enddo

		call input_error('cannot map amino-acid number to chain position',0)
	end function chain_site_from_amino_acid
	    
	    subroutine read_nmr_site_input(input, delim)
		implicit none
		character(len=*), intent(in) :: input, delim
		character(len=256) :: work,token
		integer :: first,last,dash_at,site_start,site_end,site_id,chain_site,ios,i

		work=adjustl(input)
		nmr_site_num=0
		nmr_site_ID=0
		if (len_trim(work) == 0) return
		do i=1,len_trim(work)
			if (work(i:i) == ',' .or. work(i:i) == achar(9)) work(i:i)=' '
		enddo

		first=1
		do while(first <= len_trim(work))
			do while(first <= len_trim(work) .and. work(first:first) == ' ')
				first=first+1
			enddo
			if (first > len_trim(work)) exit
			last=first
			do while(last <= len_trim(work) .and. work(last:last) /= ' ')
				last=last+1
			enddo
			token=''
			token=work(first:last-1)
			dash_at=index(token,'-')
			if (dash_at > 1) then
				read(token(:dash_at-1),*,iostat=ios) site_start
				if (ios /= 0) call input_error('invalid grouped-site range',0)
				read(token(dash_at+1:),*,iostat=ios) site_end
				if (ios /= 0) call input_error('invalid grouped-site range',0)
			else
				read(token,*,iostat=ios) site_start
				if (ios /= 0) call input_error('invalid grouped-site index',0)
				site_end=site_start
			endif
			if (site_start < 1 .or. site_end > amino_acid_residue_count .or. &
				site_start > site_end) &
				call input_error('grouped-site amino-acid number is outside the peptide',0)
			do site_id=site_start,site_end
				chain_site=chain_site_from_amino_acid(site_id)
				if (nmr_site_num > 0) then
					if (any(nmr_site_ID(1:nmr_site_num) == chain_site)) &
						call input_error('duplicate grouped-site index',0)
				endif
				nmr_site_num=nmr_site_num+1
				if (nmr_site_num > maximum_nmr_site_num) &
					call input_error('too many grouped-site constraints',0)
				nmr_site_ID(nmr_site_num)=chain_site
			enddo
			first=last+1
		enddo
    end subroutine read_nmr_site_input
    
	subroutine read_void_site_input(input)
		implicit none
		character(len=*), intent(in) :: input
		character(len=256) :: s, tok, right_part
		integer :: i, ios, pos, pos_dash, site_start, site_end, chain_site

		s = adjustl(trim(input))
		if (len_trim(s) == 0) then
			return
		end if
		
		! convert commas to spaces
		do i = 1, len_trim(s)
			if (s(i:i) == ',') s(i:i) = ' '
		end do
		do
			tok = ''
			read(s,*,iostat=ios) tok
			if (ios /= 0) exit
			if (len_trim(tok) == 0) exit

			tok = trim(tok)

			pos = index(tok,':')
			pos_dash = index(tok,'-')
			if (pos > 0) then
				read(tok(1:pos-1), *) site_start
				right_part = adjustl(tok(pos+1:))
				if (site_start < 1 .or. site_start > amino_acid_residue_count) then
					stop "ERROR: invalid single positional constraint site."
				endif
				chain_site=chain_site_from_amino_acid(site_start)
				if (len_trim(right_part) < 3 .or. len_trim(right_part) > 4) then
					stop "ERROR: invalid fixed residue name in single positional constraint."
				endif
				void_site_num = void_site_num + 1
				if (void_site_num > maximum_void_site_num) stop "ERROR: too many single positional constraints."
				void_site(void_site_num) = chain_site
				void_site_fixed_name(void_site_num) = adjustl(right_part)
				call normalize_void_site_fixed_name(void_site(void_site_num), void_site_fixed_name(void_site_num))
				write(*,'(A,I0,A,I0,2A)') 'amino acid ',site_start,' (chain site ', &
					chain_site,') will be fixed as ',void_site_fixed_name(void_site_num)
			elseif (pos_dash > 0) then
				read(tok(1:pos_dash-1), *) site_start
				read(tok(pos_dash+1:), *) site_end
				if (site_start < 1 .or. site_end > amino_acid_residue_count .or. &
					site_start > site_end) then
					stop "ERROR: invalid single positional constraint site."
				endif
				do i = site_start, site_end
					chain_site=chain_site_from_amino_acid(i)
					void_site_num = void_site_num + 1
					if (void_site_num > maximum_void_site_num) stop "ERROR: too many single positional constraints."
					void_site(void_site_num) = chain_site
				enddo
			else
				read(tok, *) site_start
				if (site_start < 1 .or. site_start > amino_acid_residue_count) then
					stop "ERROR: invalid single positional constraint site."
				endif
				chain_site=chain_site_from_amino_acid(site_start)
				void_site_num = void_site_num + 1
				if (void_site_num > maximum_void_site_num) stop "ERROR: too many single positional constraints."
				void_site(void_site_num) = chain_site
			endif

			! remove the token we just read from the front of s
			s = adjustl(s(len_trim(tok)+1:))
			if (len_trim(s) == 0) exit
		end do
	end subroutine read_void_site_input

	subroutine normalize_void_site_fixed_name(site_id, residue_name)
		implicit none
		integer, intent(in) :: site_id
		character*4, intent(inout) :: residue_name
		integer :: i, ich

		residue_name = adjustl(residue_name)
		do i = 1, len_trim(residue_name)
			ich = iachar(residue_name(i:i))
			if (ich >= iachar('a') .and. ich <= iachar('z')) residue_name(i:i) = achar(ich - 32)
		end do

		if (residue_name /= "ACE" .and. residue_name /= "NME" .and. residue_name /= "NHE" .and. len_trim(residue_name) == 3) then
			if (site_id == 1) then
				residue_name = "N" // residue_name(1:3)
			elseif (site_id == gnum) then
				residue_name = "C" // residue_name(1:3)
			endif
		endif
	end subroutine normalize_void_site_fixed_name
    
	    subroutine read_restraints_input(input, signal)
		implicit none
		character(len=*), intent(in) :: input
		character(len=100)           :: temp
	        character*4                  :: AA_list(10)
		integer, intent(in)           :: signal
		integer                      :: i, res_num

	        res_num = 0
		AA_list=''
	        temp = input
		call upcase(temp)
		if (signal == 3) then
			NMR_AA_pool=''
			NMR_pool_size=0
		endif
		do i = 1, len_trim(temp)
			if (temp(i:i) == ',') temp(i:i) = ' '
		enddo
        do while (len_trim(temp) > 0)
            i = index(temp, ' ')   
            if (i == 0) then
				if (len_trim(temp) > 0) then                        ! avoid spaces
					if (res_num >= size(AA_list)) &
						call input_error('too many amino acids in grouped-site AA pool',0)
					res_num = res_num + 1
					AA_list(res_num) = trim(adjustl(temp))
                endif
                exit
            else 
	                if (len_trim(temp(:i-1)) > 0) then                  ! avoid spaces
					if (res_num >= size(AA_list)) &
						call input_error('too many amino acids in grouped-site AA pool',0)
					res_num = res_num + 1
					AA_list(res_num) = trim(adjustl(temp(:i-1)))
                endif
                temp = temp(i+1:)
            end if
        end do
        
        select case(signal)
            case (1)
				do i=1, res_num
					NMR_restraint_AA(i) = AA_list(i)
				enddo
				restraint_num = res_num
				write(*,*) "Number of necessary AA: ", restraint_num
				write(*,*) "Necessary AA: ", NMR_restraint_AA
            case (2)
				do i=1, res_num
					exclusion_AA(i) = AA_list(i)
				enddo
				exclusion_num = res_num
				write(*,*) "Number of avoided AA: ", exclusion_num
				write(*,*) "Avoided AA: ", exclusion_AA                
            case (3)
                if (len_trim(input) == 0) then
					NMR_pool_size = 0
					return
				end if
				do i=1, res_num
					NMR_AA_pool(i) = AA_list(i)
                enddo
                NMR_pool_size = res_num
				write(*,*) "Grouped positional constraint AA pool: ", NMR_AA_pool
                
			end select
        return
	end subroutine read_restraints_input
    
    subroutine readdir()
    integer :: realpath
    character(len=1000) :: path, rundir
    logical :: back=.true.
    character(len=100) :: os_name
    character*1 :: sep
    
    call get_environment_variable("OS",os_name)    ! get platform name
    call getcwd(rundir)                            ! path to working directory
	call get_command_argument(0,path)              ! path of the executable file (contains exe file name)
    
    if (index(trim(os_name), "Windows_NT") > 0) then
            sep =  '\'
    else
            sep = '/'
    end if
    
	realpath = scan (path,sep,back)                
	mydir = path(1:realpath)                       ! path to where the executable file located
    
    return
end subroutine readdir
    
subroutine check_help
	implicit none
	integer :: number_of_arguments
	character(len=256) :: argument

	number_of_arguments=command_argument_count()
	if (number_of_arguments == 0) return

	call get_command_argument(1,argument)
	call upcase(argument)
	if (number_of_arguments /= 1 .or. &
		(trim(argument) /= '-H' .and. trim(argument) /= '--HELP')) then
		write(*,'(2A)') 'ERROR: unsupported command-line option: ',trim(argument)
		write(*,'(A)') 'PepAD design parameters belong in input.txt.'
		write(*,'(A)') 'The only command-line options are -h and --help.'
		stop 1
	endif

	write(*,'(A)') 'PepAD v1.42 help'
	write(*,'(A)') '======================================================================'
	write(*,'(A)') 'Usage:'
	write(*,'(A)') '  PepAD                 Run a design using input.txt'
	write(*,'(A)') '  PepAD -h              Display this help and exit'
	write(*,'(A)') '  PepAD --help          Display this help and exit'
	write(*,'(A)') ' '
	write(*,'(A)') 'Input format: FLAG = value'
	write(*,'(A)') 'Flag order does not matter. Comments may start with # or !.'
	write(*,'(A)') 'Each input flag may appear only once.'
	write(*,'(A)') 'Only PDBFILE and N_STEPS are required.'
	write(*,'(A)') ' '
	write(*,'(A)') 'Required parameters:'
	write(*,'(A)') '  PDBFILE                       String; initial PDB structure filename'
	write(*,'(A)') '  N_STEPS                       Integer; number of Monte Carlo steps'
	write(*,'(A)') '-----------------------------------'
	write(*,'(A)') 'The following parameters are optional.'
	write(*,'(A)') 'Initial structure parameters:'
	write(*,'(A)') '  N_RESIDUES                    Integer; chain length including caps'
	write(*,'(A)') '  N_CHAINS                      Integer; total number of peptide chains'
	write(*,'(A)') '  GROUP_1, GROUP_2              Integer lists; chain IDs in each'
	write(*,'(A)') '                                energy-calculation group'
	write(*,'(A)') '  SHEET_1, SHEET_2, ...         Integer lists; chain IDs in each'
	write(*,'(A)') '                                beta-sheet'
	write(*,'(A)') '  RESTART                       Integer; 0 = fresh run, 1 = restart'
	write(*,'(A)') '                                (default 0)'
    write(*,'(A)') '  NOTE:'
	write(*,'(A)') '  If N_RESIDUES and N_CHAINS are omitted, both values are'
	write(*,'(A)') '  determined automatically from the PDB file.'
	write(*,'(A)') '  If either N_RESIDUES or N_CHAINS is provided, both must be provided.'
	write(*,'(A)') '  If GROUP_* and SHEET_* are omitted, PepAD determines the'
	write(*,'(A)') '  energy groups and beta-sheets using DSSP-based analysis.'
	write(*,'(A)') ' '
	write(*,'(A)') 'Monte Carlo and energy parameters:'
	write(*,'(A)') '  RANSEED                       Integer; random seed'
	write(*,'(A)') '                                (default: system clock)'
	write(*,'(A)') '  KBT_SEQ                       Real; sequence kBT used when'
	write(*,'(A)') '                                ANNEAL_STAGES = 0 (default 1.0)'
	write(*,'(A)') '  KBT_SEQ_HIGH                  Real; initial sequence kBT when'
	write(*,'(A)') '                                ANNEAL_STAGES > 0 (default 2.0)'
	write(*,'(A)') '  KBT_SEQ_LOW                   Real; final sequence kBT when'
	write(*,'(A)') '                                ANNEAL_STAGES > 0 (default 0.5)'
	write(*,'(A)') '  KBT_SHEETMOVE                 Real; sheet-move kBT used when'
	write(*,'(A)') '                                ANNEAL_STAGES = 0 (default 0.6)'
	write(*,'(A)') '  KBT_SHEETMOVE_HIGH            Real; initial sheet-move kBT when'
	write(*,'(A)') '                                ANNEAL_STAGES > 0 (default 1.2)'
	write(*,'(A)') '  KBT_SHEETMOVE_LOW             Real; final sheet-move kBT when'
	write(*,'(A)') '                                ANNEAL_STAGES > 0 (default 0.3)'
	write(*,'(A)') '  RMSD_MAX                      Real; maximum RMSD in angstroms'
	write(*,'(A)') '                                (default 3.0)'
	write(*,'(A)') '  PAGG_WEIGHT                   Real; aggregation propensity weight'
	write(*,'(A)') '                                (default 2.0)'
	write(*,'(A)') '  SHEETMOVE                     Integer; 0 = off, 1 = on'
	write(*,'(A)') '                                (default 1)'
	write(*,'(A)') '  ANNEAL_STAGES                 Integer; 0 = no annealing,'
	write(*,'(A)') '                                >0 = number of annealing stages'
	write(*,'(A)') '                                (default 0)'
	write(*,'(A)') '  ESURF_MODE                    Integer; 0 = off,'
	write(*,'(A)') '                                1 = direct full ARVO reference,'
	write(*,'(A)') '                                2 = cached/incremental ARVO (default 0)'
	write(*,'(A)') '  STEPS_BEFORE_SHEETMOVE        Integer; initial waiting steps'
	write(*,'(A)') '                                (default 500)'
	write(*,'(A)') '  STEPS_BETWEEN_SHEETMOVE       Integer; sheet-move interval'
	write(*,'(A)') '                                (default 200)'
	write(*,'(A)') ' '
	write(*,'(A)') 'Composition parameters:'
	write(*,'(A)') '  N_HYDROPHOBIC_MIN             Integer; minimum number of'
	write(*,'(A)') '                                hydrophobic residues'
	write(*,'(A)') '  N_HYDROPHOBIC_MAX             Integer; maximum number of'
	write(*,'(A)') '                                hydrophobic residues'
	write(*,'(A)') '  N_POLAR_MIN                   Integer; minimum number of'
	write(*,'(A)') '                                polar residues'
	write(*,'(A)') '  N_POLAR_MAX                   Integer; maximum number of'
	write(*,'(A)') '                                polar residues'
	write(*,'(A)') '  N_CHARGED_MIN                 Integer; minimum number of'
	write(*,'(A)') '                                charged residues'
	write(*,'(A)') '  N_CHARGED_MAX                 Integer; maximum number of'
	write(*,'(A)') '                                charged residues'
	write(*,'(A)') '  N_OTHER_MIN                   Integer; minimum number of'
	write(*,'(A)') '                                other residues'
	write(*,'(A)') '  N_OTHER_MAX                   Integer; maximum number of'
	write(*,'(A)') '                                other residues'
	write(*,'(A)') ' '
	write(*,'(A)') '  If all category flags are omitted, the initial PDB composition'
	write(*,'(A)') '  is used.'
	write(*,'(A)') '  In custom mode, every category requires MIN, MAX, or both.'
	write(*,'(A)') '  An omitted MIN becomes 0; an omitted MAX becomes N_RESIDUES.'
	write(*,'(A)') '  Equal minimum and maximum values indicate a fixed category number.'
	write(*,'(A)') '  A category MIN or MAX flag may appear only once.'
	write(*,'(A)') '  Composition totals exclude ACE, NME, and NHE cap residues.'
	write(*,'(A)') '  Single-site and grouped-site amino acids are included in composition totals.'
	write(*,'(A)') '  Fixed single-site amino acids contribute to the minimum required AA/category counts.'
	write(*,'(A)') ' '

	write(*,'(A)') 'Compositional constraints:'
	write(*,'(A)') '  N_<AA>_MIN                    Integer; minimum count of AA'
	write(*,'(A)') '  N_<AA>_MAX                    Integer; maximum count of AA'
	write(*,'(A)') '  AA = GLY ALA VAL LEU ILE MET PHE TYR TRP SER ASN GLN THR HIE'
	write(*,'(A)') '       ARG LYS GLU ASP CYS PRO'
	write(*,'(A)') '  Users may provide only a minimum, only a maximum, or both.'
	write(*,'(A)') '  Equal minimum and maximum values impose an exact AA count.'
	write(*,'(A)') '  Omitted amino-acids are unconstrained.'
	write(*,'(A)') '  Each N_<AA>_MIN or N_<AA>_MAX flag may appear only once.'
	write(*,'(A)') '  Specific-AA minima cannot exceed the corresponding category maximum.'
	write(*,'(A)') ' '

	write(*,'(A)') 'Positional constraints:'
	write(*,'(A)') '  ACE, NME, and NHE cap positions are handled automatically.'
	write(*,'(A)') '  Site numbers refer to amino acids only; cap positions are skipped.'
	write(*,'(A)') '  Example: amino acid 3 is chain site 4 when chain site 1 is an ACE cap.'
	write(*,'(A)') '  SINGLE_SITE_CONSTRAINTS       Sites remain unchanged during the design run'
    write(*,'(A)') '                                Usage: '
    write(*,'(A)') '                                "Site number only": retain the amino acid from initial structure'
    write(*,'(A)') '                                Example: 3 5 (3rd and 5 th amino acid are not changed)'
    write(*,'(A)') '                                "Site number:AA": change the site to the specified amino acid'
    write(*,'(A)') '                                and keep it unchanged'
	write(*,'(A)') '                                Example: 3:CYS 12:LEU'
	write(*,'(A)') '                               (Separate entries with spaces, do not place spaces around the colon)'
	write(*,'(A)') ' '
	write(*,'(A)') '  GROUPED_SITE_AA_POOL          Allowed amino acids for grouped sites'
	write(*,'(A)') '                                Example: ASN SER ALA ALA'
	write(*,'(A)') '                                Repeated amino acids are allowed to be selected by the amount entered'
	write(*,'(A)') '                                Each remaining copy has an equal selection probability'
	write(*,'(A)') '                                Number of pool entries must be at least the number'
	write(*,'(A)') '                                of GROUPED_SITE_CONSTRAINTS'
	write(*,'(A)') ' '
	write(*,'(A)') '  GROUPED_SITE_CONSTRAINTS      Sites select from the available grouped-pool copies'
	write(*,'(A)') '                                Example: 4 6 10 12'
	write(*,'(A)') '  PepAD checks grouped-pool category minima/maxima against composition limits.'
	write(*,'(A)') '  Caps and positional-constraint sites are excluded from ordinary MC site selection.'
	write(*,'(A)') ' '
	!write(*,'(A)') 'CONROT parameters:'
	!write(*,'(A)') '  CONROT_FLAG                   Enable CONROT: 0 or 1 (default 0)'
	!write(*,'(A)') '  CONROT_ANGLE_LIMIT            Angle limit (default 2.0)'
	!write(*,'(A)') '  CONROT_STEP                   Angle step (default 0.1)'
	!write(*,'(A)') '  CONROT_CLOSURE_TOL            Closure tolerance (default 0.10)'
	!write(*,'(A)') '  CONROT_RMSD_LIMIT             RMSD limit (default 0.5)'
	write(*,'(A)') '======================================================================'
	stop
	end subroutine check_help
    
subroutine banner()
    implicit none
	integer :: environment_status, command_status, exit_status
	character(len=32) :: os_name
	character(len=128) :: locale_name, command_message
	logical :: unicode_output

	unicode_output=.true.
	os_name=''
	call get_environment_variable('OS',os_name,status=environment_status)
	if (index(os_name,'Windows_NT') > 0) then
		command_message=''
		call execute_command_line('chcp 65001 > nul',wait=.true., &
			exitstat=exit_status,cmdstat=command_status,cmdmsg=command_message)
		if (command_status /= 0 .or. exit_status /= 0) unicode_output=.false.
	else
		! Use ASCII only when Linux explicitly reports a non-UTF-8 locale.
		locale_name=''
		call get_environment_variable('LC_ALL',locale_name,status=environment_status)
		if (len_trim(locale_name) == 0) &
			call get_environment_variable('LC_CTYPE',locale_name,status=environment_status)
		if (len_trim(locale_name) == 0) &
			call get_environment_variable('LANG',locale_name,status=environment_status)
		call upcase(locale_name)
		if (len_trim(locale_name) > 0 .and. index(locale_name,'UTF-8') == 0 &
			.and. index(locale_name,'UTF8') == 0) unicode_output=.false.
	endif

    write(*,'(A)') "======================================================================"
	if (unicode_output) then
		write(*,'(A)') " ██████╗ ███████╗██████╗  █████╗ ██████╗"
		write(*,'(A)') " ██╔══██╗██╔════╝██╔══██╗██╔══██╗██╔══██╗"
		write(*,'(A)') " ██████╔╝█████╗  ██████╔╝███████║██║  ██║"
		write(*,'(A)') " ██╔═══╝ ██╔══╝  ██╔═══╝ ██╔══██║██║  ██║"
		write(*,'(A)') " ██║     ███████╗██║     ██║  ██║██████╔╝"
		write(*,'(A)') " ╚═╝     ╚══════╝╚═╝     ╚═╝  ╚═╝╚═════╝"
    else
        

        write(*,'(A)') '  _____    _______    _____       __      ________  '
		write(*,'(A)') ' |  __ \  | ______|  |  __ \     /  \     | _____ \ '
		write(*,'(A)') ' | |  \ | | |_____   | |  \ |   / /\ \    | |    \ \'
		write(*,'(A)') ' | |__/ / | ______|  | |__/ /  / /__\ \   | |    | |'
        write(*,'(A)') ' |  ___/  | |_____   |  ___/  / _____\ \  | |____/ /'
		write(*,'(A)') ' | |      |_______|  | |     /_/      \_\ |_______/ '
		write(*,'(A)') ' |_|                 |_|                            '

        
		!write(*,'(A)') " PPPPP   EEEEE  PPPPP    AAA    DDDD"
		!write(*,'(A)') " P    P  E      P    P  A   A   D   D"
		!write(*,'(A)') " PPPPP   EEEE   PPPPP   AAAAA   D   D"
		!write(*,'(A)') " P       E      P       A   A   D   D"
		!write(*,'(A)') " P       EEEEE  P       A   A   DDDD"
	endif
    write(*,'(A)') " "
    write(*,'(A)') "  Peptide Assembly Design Algorithm version 1.42"
    write(*,'(A)') "  Authors: Haoyu Wang, Sudeep Sarma, and Carol K. Hall"
    write(*,'(A)') "  Department of Chemical and Biomolecular Engineering"
    write(*,'(A)') "  North Carolina State University"
    write(*,'(A)') "======================================================================"
    write(*,'(A)') " "

end subroutine banner

end module input


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mathfunction

	use constant

	contains
	subroutine normalvector(rsta, rmid, rend, r_nor)
	implicit none
	real						:: rsta(3), rmid(3), rend(3), r_nor(3)
	real						:: a(3), b(3)

	a=rsta-rmid
	b=rend-rmid

	r_nor(1)=a(2)*b(3)-a(3)*b(2)
	r_nor(2)=a(3)*b(1)-a(1)*b(3)
	r_nor(3)=a(1)*b(2)-a(2)*b(1)

	return
	end subroutine normalvector


	subroutine vectorrotation(rsta, rend, m)
	implicit none
	real						:: rsta(3), rend(3), m(3,3)
	real						:: r_cropro(3), a(3), a1(3,3), a2(3,3), a3(3,3)
	real						:: absrsta, absrend, r_dotpro, cos, sin

	absrsta=sqrt(rsta(1)*rsta(1)+rsta(2)*rsta(2)+rsta(3)*rsta(3))
	absrend=sqrt(rend(1)*rend(1)+rend(2)*rend(2)+rend(3)*rend(3))

	r_dotpro=dot_product(rsta, rend)

	r_cropro(1)=rsta(2)*rend(3)-rsta(3)*rend(2)
	r_cropro(2)=rsta(3)*rend(1)-rsta(1)*rend(3)
	r_cropro(3)=rsta(1)*rend(2)-rsta(2)*rend(1)

	cos=r_dotpro/(absrsta*absrend)
	sin=sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))/(absrsta*absrend)

	a(1)=r_cropro(1)/sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))
	a(2)=r_cropro(2)/sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))
	a(3)=r_cropro(3)/sqrt(r_cropro(1)*r_cropro(1)+r_cropro(2)*r_cropro(2)+r_cropro(3)*r_cropro(3))

	a1(1,1)=a(1)*a(1)
	a1(1,2)=a(1)*a(2)
	a1(1,3)=a(1)*a(3)
	a1(2,1)=a(2)*a(1)
	a1(2,2)=a(2)*a(2)
	a1(2,3)=a(2)*a(3)
	a1(3,1)=a(3)*a(1)
	a1(3,2)=a(3)*a(2)
	a1(3,3)=a(3)*a(3)

	a2=-a1
	a2(1,1)=1+a2(1,1)
	a2(2,2)=1+a2(2,2)
	a2(3,3)=1+a2(3,3)
	a2=cos*a2

	a3(1,1)=0.0
	a3(1,2)=-a(3)
	a3(1,3)=a(2)
	a3(2,1)=a(3)
	a3(2,2)=0.0
	a3(2,3)=-a(1)
	a3(3,1)=-a(2)
	a3(3,2)=a(1)
	a3(3,3)=0
	a3=sin*a3

	m=a1+a2+a3
	m=transpose(m)

	return	
	end subroutine vectorrotation

	
	subroutine axisrotation(a, cos, sin, m)
	implicit none
	real					:: cos, sin
	real					:: a(3), a1(3,3), a2(3,3), a3(3,3)
	real					:: m(3,3)

	a1(1,1)=a(1)*a(1)
	a1(1,2)=a(1)*a(2)
	a1(1,3)=a(1)*a(3)
	a1(2,1)=a(2)*a(1)
	a1(2,2)=a(2)*a(2)
	a1(2,3)=a(2)*a(3)
	a1(3,1)=a(3)*a(1)
	a1(3,2)=a(3)*a(2)
	a1(3,3)=a(3)*a(3)

	a2=-a1
	a2(1,1)=1+a2(1,1)
	a2(2,2)=1+a2(2,2)
	a2(3,3)=1+a2(3,3)
	a2=cos*a2

	a3(1,1)=0.0
	a3(1,2)=-a(3)
	a3(1,3)=a(2)
	a3(2,1)=a(3)
	a3(2,2)=0.0
	a3(2,3)=-a(1)
	a3(3,1)=-a(2)
	a3(3,2)=a(1)
	a3(3,3)=0
	a3=sin*a3

	m=a1+a2+a3
	m=transpose(m)

	return	
	end subroutine axisrotation


	subroutine transformatrix(bond_angle, dihedral_angle, T)
	implicit none
	real					:: bond_angle, dihedral_angle
	real					:: T(3,3)

	T(1,1)=cosd(bond_angle)
	T(1,2)=sind(bond_angle)
	T(1,3)=0.0
	T(2,1)=-sind(bond_angle)*cosd(dihedral_angle)
	T(2,2)=cosd(bond_angle)*cosd(dihedral_angle)
	T(2,3)=sind(dihedral_angle)
	T(3,1)=sind(bond_angle)*sind(dihedral_angle)
	T(3,2)=-cosd(bond_angle)*sind(dihedral_angle)
	T(3,3)=cosd(dihedral_angle)

	return
	end subroutine transformatrix
	
	
	subroutine phipsiomg_angle(p1, p2, p3, p4, angle)
	implicit none
	real						:: p1(3), p2(3), p3(3), p4(3)
	real						:: angle, angle_T1, angle_T2
	real						:: rsta(3), rend(3), rmid(3)
	real						:: r_1(3), r_2(3)
	real						:: absrsta, absrend, r_dotpro, value

	call normalvector(p1, p2, p3, rend)
	call normalvector(p2, p3, p4, rsta)

	absrsta=sqrt(rsta(1)*rsta(1)+rsta(2)*rsta(2)+rsta(3)*rsta(3))
	absrend=sqrt(rend(1)*rend(1)+rend(2)*rend(2)+rend(3)*rend(3))

	r_dotpro=dot_product(rsta, rend)
	value=r_dotpro/(absrsta*absrend)
	if(value.lt.-1.00) then
		value=-1.00
	elseif(value.gt.1.00) then	
		value=1.00
	endif
	angle_T1=acosd(value)

	if(abs(180.0-angle_T1).le.(0.5)) then
		angle=180.0
	elseif(abs(angle_T1).le.(0.4)) then
		angle=0.0
	else
		rmid=0.0
		r_2=p3-p2
		call normalvector(rsta, rmid, rend, r_1)

		absrsta=sqrt(r_1(1)*r_1(1)+r_1(2)*r_1(2)+r_1(3)*r_1(3))
		absrend=sqrt(r_2(1)*r_2(1)+r_2(2)*r_2(2)+r_2(3)*r_2(3))

		r_dotpro=dot_product(r_1, r_2)
		if(abs(r_dotpro/(absrsta*absrend)-1.0).le.(0.1)) then
			angle_T2=0.00
		elseif(abs(r_dotpro/(absrsta*absrend)+1.0).le.(0.1)) then
			angle_T2=180
		else
			value=r_dotpro/(absrsta*absrend)
			if(value.lt.-1.00) then
				value=-1.00
			elseif(value.gt.1.00) then	
				value=1.00
			endif	
			angle_T2=acosd(value)
		endif
		if(angle_T2.gt.90) then
			angle=-angle_T1
		else
			angle=angle_T1
		endif
	endif

	return
	end subroutine phipsiomg_angle

	
	subroutine variance_covariance_matrix(matrix,obs,n,det)
	Implicit none

	integer		     					:: obs, n, i, j, k
	real								:: matrix(34,4),a(n,n),avg(n),matrix_dif(obs,n)
	real								:: det

	do i=1,obs
		do j=1, n
			if(matrix(i,j).gt.180) then
					matrix(i,j)=matrix(i,j)-360.0
			elseif(matrix(i,j).lt.-180) then
					matrix(i,j)=matrix(i,j)+360.0
			endif
		enddo
	enddo
	matrix=matrix*0.0174533

	avg=0.0
	do j=1, n
		do i=1, obs
			avg(j)=avg(j)+matrix(i,j)
		enddo
	enddo
	avg=avg/obs

	do j=1, n
		do i=1, obs
			matrix_dif(i,j)=matrix(i,j)-avg(j)
		enddo
	enddo

	a=0.0
	do j=1, n
		do k=1, j
			do i=1, obs
				a(j,k)=a(j,k)+matrix_dif(i,j)*matrix_dif(i,k)
			enddo
			a(j,k)=a(j,k)/(obs-1)
			a(k,j)=a(j,k)
		enddo
	enddo

	call eigenfunction(n,a,det)

	return
	end subroutine variance_covariance_matrix


	subroutine eigenfunction(n,a,det)
	implicit none
	
	integer		:: n,i,j,p,q
	real		:: amax,temp,zemp,coo,sii,co,si,app,aqq,apq,api,aqi
	real		:: rip,riq,det
	real		:: a(n,n),r(n,n)
	character	:: name*12

	do i=1,n
		do j=1,n
			r(i,j)=0.0
		enddo
		r(i,i)=1.0
	enddo

10	amax=abs(a(2,1))
	p=2
	q=1
	do i=2,n
		do j=1,i-1
			if(abs(a(i,j)).gt.amax) then
				amax=abs(a(i,j))
	 			p=i
				q=j
			endif
		enddo
	enddo

	if(amax.le.1.0e-7) then
		goto 20
	endif

	temp=2*a(p,q)/(a(p,p)-a(q,q)+1.0e-30)
	zemp=(a(p,p)-a(q,q))/(2*a(p,q))
	
	if(abs(temp).lt.1.0) then
		coo=(1+temp**2)**(-0.5)
		sii=temp*(1+temp**2)**(-0.5)
	else
		coo=abs(zemp)*(1+zemp**2)**(-0.5)
		sii=sign(1.0,zemp)*(1+zemp**2)**(-0.5)
	endif
	
	co=sqrt(0.5*(1+coo))
	si=sii/(2.0*co)

	do i=1,n
		rip=r(i,p)*co+r(i,q)*si
		riq=-r(i,p)*si+r(i,q)*co
		r(i,p)=rip
		r(i,q)=riq
	enddo

	app=a(p,p)*co**2+a(q,q)*si**2+2*a(p,q)*co*si
	aqq=a(p,p)*si**2+a(q,q)*co**2-2*a(p,q)*co*si
	apq=0.5*(a(q,q)-a(p,p))*sii+a(p,q)*coo
	a(p,p)=app
	a(q,q)=aqq
	a(p,q)=apq
	a(q,p)=a(p,q)
	
	do i=1,n
		if(i.eq.p.or.i.eq.q) then
		else
			api=a(p,i)*co+a(q,i)*si
			aqi=-a(p,i)*si+a(q,i)*co
			a(p,i)=api
			a(q,i)=aqi
			a(i,p)=a(p,i)
			a(i,q)=a(q,i)
		endif
	enddo	
    goto 10

20	continue
	det=a(1,1)
	do i=2,n
		det=det*a(i,i)
	enddo

	return
	end	subroutine eigenfunction

end module mathfunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module database

	use constant
	use randomgenerator
	use mathfunction
	use pdbfile, only: count_amino_acid_categories

	integer, parameter				:: max_library_cache_entries=160
	integer, parameter				:: max_forcefield_atoms=120
	integer, parameter				:: max_atomlink_atoms=120
	integer, parameter				:: max_atomlink_per_atom=4

	type forcefield_cache_entry
		character*4					:: gtype=''
		integer						:: atom_count=0
		character*4					:: igraph(max_forcefield_atoms)
		real						:: charge(max_forcefield_atoms), epsion(max_forcefield_atoms)
		real						:: r(max_forcefield_atoms), rborn(max_forcefield_atoms)
		real						:: fs(max_forcefield_atoms), dielecons(max_forcefield_atoms)
		integer						:: atomid(max_forcefield_atoms)
		logical						:: loaded=.false.
	end type

	type atomlink_cache_entry
		character*4					:: gtype=''
		integer						:: natom=0
		integer						:: linknum(max_atomlink_atoms)
		integer						:: linkindex(max_atomlink_atoms,max_atomlink_per_atom)
		logical						:: loaded=.false.
	end type

	type dihedral_cache_entry
		character*4					:: gtype=''
		integer						:: dihedral_num=0
		type(dihedralparameters)	:: dihedral
		logical						:: loaded=.false.
	end type

	type(forcefield_cache_entry), save	:: forcefield_cache(max_library_cache_entries)
	type(atomlink_cache_entry), save	:: atomlink_cache(max_library_cache_entries)
	type(dihedral_cache_entry), save	:: dihedral_cache(max_library_cache_entries)
	integer, save						:: forcefield_cache_count=0
	integer, save						:: atomlink_cache_count=0
	integer, save						:: dihedral_cache_count=0

	contains
	subroutine pickupsite(length,ic)
	implicit none
	integer							:: length, ic, ic1
	real							:: ran2

	if(length==0) then
		call ran_gen(ran2, 0)
		ic1=int(ran2*sitenum4mutation-1.0e-3)+1
		if(ic1.gt.sitenum4mutation) ic1=sitenum4mutation
		ic=sitenum4mutation_group(ic1)
	else
		call ran_gen(ran2, 0)
		ic=int(ran2*length-1.0e-3)+1
		if(ic.gt.length) ic=length
	endif
	
	return
	end subroutine pickupsite   

	subroutine ensure_forcefield_cache(gtype, cache_id)
    ! 2026/04/16 cache forcefield parameters
	! For every mutation: call energy_parameter -> call ensure_forcefield_cache
    !
    ! ensure_forcefield_cache first checks whether that residue type is already cached:
	!	If cached: it returns immediately, with no open and no close
	!	If not cached yet: it opens the file once, reads it once, closes it once, and stores it in memory
	!
    ! So the new behavior is:
	!	one open/read/close per residue type, not per residue evaluation
	!
    ! That means:
	!	old version: ALA file may be opened hundreds or thousands of times
	!	new version: ALA file should be opened only the first time ALA is encountered
	! functions:
	!	cache check
	!	one-time open/read/close
	!	energy_parameter calling the cache helper in main_v1.37.f90 (line 2139)
    
	implicit none
	integer							:: cache_id, status
	real							:: charge, epsion, r, rborn, fs, dielecons
	character*4						:: gtype, lbres, igraph
	integer							:: atomid

    ! loop all loaded residue check if it is loaded
	do cache_id=1, forcefield_cache_count
		if(forcefield_cache(cache_id)%loaded.and.forcefield_cache(cache_id)%gtype==gtype) return
	enddo

	forcefield_cache_count=forcefield_cache_count+1
	if(forcefield_cache_count.gt.max_library_cache_entries) stop "ERROR: forcefield cache is full."
	cache_id=forcefield_cache_count
	forcefield_cache(cache_id)%gtype=gtype
	forcefield_cache(cache_id)%atom_count=0
	forcefield_cache(cache_id)%igraph='    '
	forcefield_cache(cache_id)%charge=0.0
	forcefield_cache(cache_id)%epsion=0.0
	forcefield_cache(cache_id)%r=0.0
	forcefield_cache(cache_id)%rborn=0.0
	forcefield_cache(cache_id)%fs=0.0
	forcefield_cache(cache_id)%dielecons=0.0
	forcefield_cache(cache_id)%atomid=0

	open(10, file=trim(mydir)//'lib/ForceField/'//trim(adjustl(trim(gtype))), status="old")
	read(10, *)
	do while(.true.)
		read(10, 20, iostat=status) lbres, igraph, charge, epsion, r, rborn, fs, dielecons, atomid
		if(status.ne.0) exit
		forcefield_cache(cache_id)%atom_count=forcefield_cache(cache_id)%atom_count+1
		forcefield_cache(cache_id)%igraph(forcefield_cache(cache_id)%atom_count)=igraph
		forcefield_cache(cache_id)%charge(forcefield_cache(cache_id)%atom_count)=charge
		forcefield_cache(cache_id)%epsion(forcefield_cache(cache_id)%atom_count)=epsion
		forcefield_cache(cache_id)%r(forcefield_cache(cache_id)%atom_count)=r
		forcefield_cache(cache_id)%rborn(forcefield_cache(cache_id)%atom_count)=rborn
		forcefield_cache(cache_id)%fs(forcefield_cache(cache_id)%atom_count)=fs
		forcefield_cache(cache_id)%dielecons(forcefield_cache(cache_id)%atom_count)=dielecons
		forcefield_cache(cache_id)%atomid(forcefield_cache(cache_id)%atom_count)=atomid
	enddo
	close(10)
	forcefield_cache(cache_id)%loaded=.true.

20	format(2a4, 6e16.8, i8)
	return
	end subroutine ensure_forcefield_cache

	subroutine ensure_atomlink_cache(gtype, cache_id)
	implicit none
	integer							:: cache_id, status, id, linknum, k
	integer							:: linkindex(max_atomlink_per_atom)
	character*4						:: gtype

	do cache_id=1, atomlink_cache_count
		if(atomlink_cache(cache_id)%loaded.and.atomlink_cache(cache_id)%gtype==gtype) return
	enddo

	atomlink_cache_count=atomlink_cache_count+1
	if(atomlink_cache_count.gt.max_library_cache_entries) stop "ERROR: atomlink cache is full."
	cache_id=atomlink_cache_count
	atomlink_cache(cache_id)%gtype=gtype
	atomlink_cache(cache_id)%natom=0
	atomlink_cache(cache_id)%linknum=0
	atomlink_cache(cache_id)%linkindex=0

	open(10, file=trim(mydir)//'lib/Atomlink/'//trim(adjustl(trim(gtype))), status="old")
	do while(.true.)
		read(10, 20, iostat=status) id, linknum, (linkindex(k), k=1, linknum)
		if(status.ne.0) exit
		atomlink_cache(cache_id)%linknum(id)=linknum
		do k=1, linknum
			atomlink_cache(cache_id)%linkindex(id,k)=linkindex(k)
		enddo
		atomlink_cache(cache_id)%natom=max(atomlink_cache(cache_id)%natom, id)
	enddo
	close(10)
	atomlink_cache(cache_id)%loaded=.true.

20	format(i6, i7, 4i3)
	return
	end subroutine ensure_atomlink_cache

	subroutine ensure_dihedral_cache(gtype, cache_id)
	implicit none
	integer							:: cache_id, i, j
	character*4						:: gtype

	do cache_id=1, dihedral_cache_count
		if(dihedral_cache(cache_id)%loaded.and.dihedral_cache(cache_id)%gtype==gtype) return
	enddo

	dihedral_cache_count=dihedral_cache_count+1
	if(dihedral_cache_count.gt.max_library_cache_entries) stop "ERROR: dihedral cache is full."
	cache_id=dihedral_cache_count
	dihedral_cache(cache_id)%gtype=gtype
	dihedral_cache(cache_id)%dihedral_num=0

	open(10, file=trim(mydir)//'lib/DihedralAngle/'//trim(gtype), status="old")
		read(10, "(i8)") dihedral_cache(cache_id)%dihedral_num
		do i=1, dihedral_cache(cache_id)%dihedral_num
			read(10,"(5i8)") dihedral_cache(cache_id)%dihedral%iph(i), dihedral_cache(cache_id)%dihedral%jph(i), &
							dihedral_cache(cache_id)%dihedral%kph(i), dihedral_cache(cache_id)%dihedral%lph(i), &
							dihedral_cache(cache_id)%dihedral%multiply(i)
			do j=1, dihedral_cache(cache_id)%dihedral%multiply(i)
				read(10,"(3e16.8)") dihedral_cache(cache_id)%dihedral%pk(i,j), dihedral_cache(cache_id)%dihedral%pn(i,j), &
								 dihedral_cache(cache_id)%dihedral%phase(i,j)
			enddo
		enddo
	close(10)
	dihedral_cache(cache_id)%loaded=.true.

	return
	end subroutine ensure_dihedral_cache

	subroutine rotamerlib
	use iso_fortran_env, only: output_unit
	implicit none
	integer							:: status, grade, rotanum, anum, num, i, j, k
	real							:: x, y, z
	character*4						:: char, atype, name

	! Temporary HPC checkpoints.
	!write(output_unit,*) 'ROTAMER CHECK 1: entered rotamerlib'
	!flush(output_unit)
	open(10, file=trim(mydir)//"lib/rotamer", status="old")
		!write(output_unit,*) 'ROTAMER CHECK 2: opened lib/rotamer'
		!flush(output_unit)
		do while(.true.)
			read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum		
			if(status.ne.0) exit
			!write(output_unit,*) 'ROTAMER CHECK 3: header = ',trim(name), &
			!	', grade = ',grade,', rotamers = ',rotanum
			!flush(output_unit)
			i=0
			if(name=="GLY") then
				i=1
			elseif(name=="LEU") then
				i=2
			elseif(name=="VAL") then
				i=3
			elseif(name=="ILE") then
				i=4
			elseif(name=="MET") then
				i=5
			elseif(name=="PHE") then
				i=6
			elseif(name=="TYR") then
				i=7
			elseif(name=="TRP") then
				i=8
			elseif(name=="ARG") then
				i=9
			elseif(name=="LYS") then
				i=10
			elseif(name=="SER") then
				i=11
			elseif(name=="THR") then
				i=12
			elseif(name=="ASN") then
				i=13
			elseif(name=="GLN") then
				i=14
			elseif(name=="HIE") then
				i=15
			elseif(name=="PRO") then
				i=16
			elseif(name=="CYS") then
				i=17
			elseif(name=="ALA") then
				i=18
			elseif(name=="GLU") then
				i=19
			elseif(name=="ASP") then
				i=20
			elseif(name=="ACE") then
				i=21
			elseif(name=="NME") then
				i=22
 			elseif(name=="NHE") then
				i=23               
			endif
				
			if(i.ne.0) then
				aa_lib(i)%gtype=name; aa_lib(i)%grade=grade; aa_lib(i)%rotanum=rotanum
				if(grade.ne.0) then
					do j=1, rotanum
						read(10,'(4I5)',iostat=status) (aa_lib(i)%dihedralangle(j,k),k=1,grade)
						if (status /= 0) then
							write(output_unit,*) 'ERROR: failed reading rotamer angles for ', &
								trim(name),', rotamer ',j
							flush(output_unit)
							stop
						endif
					enddo
				endif
				!write(output_unit,*) 'ROTAMER CHECK 4: finished entry ',trim(name)
				!flush(output_unit)
			endif
		enddo
	close(10)
	!write(output_unit,*) 'ROTAMER CHECK 5: finished lib/rotamer'
	!flush(output_unit)
	aa_lib%cnum1=0
	aa_lib%cnum2=0
	aa_lib%cnum3=0

	do i=1, 20
		open (10, file=trim(mydir)//'lib/RotamerLibrary/'//trim(aa_lib(i)%gtype), status="old")	
		!write(output_unit,*) 'ROTAMER CHECK 6: opened RotamerLibrary/'//trim(aa_lib(i)%gtype)
		!flush(output_unit)
		do while(.true.)
			read(10, *, iostat=status) char, anum, atype, name, char, num, x, y, z
			if(status.ne.0) exit
			!write(output_unit,*) 'ROTAMER CHECK 7: residue = ',trim(name), &
			!	', atom number = ',anum,', atom = ',trim(atype)
			!flush(output_unit)
			if(atype=="N".or.atype=="H".or.atype=="H1".or.atype=="H2".or.atype=="H3".or.atype=="CA" &
			   .or.atype=="HA".or.atype=="HA2".or.atype=="HA3") then  
				aa_lib(i)%cnum1=aa_lib(i)%cnum1+1
				aa_lib(i)%atype1(aa_lib(i)%cnum1)=atype
				aa_lib(i)%coo1(aa_lib(i)%cnum1,1)=x
				aa_lib(i)%coo1(aa_lib(i)%cnum1,2)=y
				aa_lib(i)%coo1(aa_lib(i)%cnum1,3)=z
			elseif(atype=="C".or.atype=="O".or.atype=="OXT") then
				aa_lib(i)%cnum3=aa_lib(i)%cnum3+1
				aa_lib(i)%atype3(aa_lib(i)%cnum3)=atype
				aa_lib(i)%coo3(aa_lib(i)%cnum3,1)=x
				aa_lib(i)%coo3(aa_lib(i)%cnum3,2)=y
				aa_lib(i)%coo3(aa_lib(i)%cnum3,3)=z
			else
				aa_lib(i)%cnum2=aa_lib(i)%cnum2+1
				aa_lib(i)%atype2(aa_lib(i)%cnum2)=atype
				aa_lib(i)%coo2(aa_lib(i)%cnum2,1)=x
				aa_lib(i)%coo2(aa_lib(i)%cnum2,2)=y
				aa_lib(i)%coo2(aa_lib(i)%cnum2,3)=z
			endif
		enddo
		close(10)
		!write(output_unit,*) 'ROTAMER CHECK 8: finished RotamerLibrary/'//trim(aa_lib(i)%gtype)
		!flush(output_unit)
	enddo
	!write(output_unit,*) 'ROTAMER CHECK 9: finished rotamerlib'
	!flush(output_unit)
	
	return
	end subroutine rotamerlib


	subroutine findrotamer(ic, group, name_original, rotanum, aa_group, ip)
	implicit none
	integer								:: categoryID, status, ic, rotanum, i, ii, j, k, l, n, ip
	integer								:: grade, grade_num(6), monitor(6)
	real								:: nr(3), car(3), cr(3), r_norpep(3)
	real								:: aa_nr(3), aa_car(3), aa_cr(3), r_norrot(3)
	real								:: r_nca(3), aa_r_nca(3), r_trans(3)
	real								:: CA(3), rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3), Tmember(15,3)
	real								:: delta_chi, cos_angle, sin_angle	
	real								:: temp1(20,3), temp2(60,3), temp3(20,3)
	character*4							:: name_original

	type(groupdetails)					:: group(repeated_unit,gnum), aa_group(repeated_unit,40)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6), Tclass(6)
	
	ip=0
	if(name_original=="GLY".or.name_original=="NGLY".or.name_original=="CGLY") then
		ip=1
	elseif(name_original=="LEU".or.name_original=="NLEU".or.name_original=="CLEU") then
		ip=2
	elseif(name_original=="VAL".or.name_original=="NVAL".or.name_original=="CVAL") then
		ip=3
	elseif(name_original=="ILE".or.name_original=="NILE".or.name_original=="CILE") then
		ip=4
	elseif(name_original=="MET".or.name_original=="NMET".or.name_original=="CMET") then
		ip=5
	elseif(name_original=="PHE".or.name_original=="NPHE".or.name_original=="CPHE") then
		ip=6
	elseif(name_original=="TYR".or.name_original=="NTYR".or.name_original=="CTYR".or. &
	       name_original=="TYX".or.name_original=="NTYX".or.name_original=="CTYX") then
		ip=7
	elseif(name_original=="TRP".or.name_original=="NTRP".or.name_original=="CTRP") then
		ip=8
	elseif(name_original=="ARG".or.name_original=="NARG".or.name_original=="CARG".or. &
	       name_original=="ARN".or.name_original=="NARN".or.name_original=="CARN") then
		ip=9
	elseif(name_original=="LYS".or.name_original=="NLYS".or.name_original=="CLYS".or. &
	       name_original=="LYN".or.name_original=="NLYN".or.name_original=="CLYN") then
		ip=10
	elseif(name_original=="SER".or.name_original=="NSER".or.name_original=="CSER") then
		ip=11
	elseif(name_original=="THR".or.name_original=="NTHR".or.name_original=="CTHR") then
		ip=12
	elseif(name_original=="ASN".or.name_original=="NASN".or.name_original=="CASN") then
		ip=13
	elseif(name_original=="GLN".or.name_original=="NGLN".or.name_original=="CGLN") then
		ip=14
	elseif(name_original=="HIE".or.name_original=="NHIE".or.name_original=="CHIE".or. &
	       name_original=="HIP".or.name_original=="NHIP".or.name_original=="CHIP") then
		ip=15
	elseif(name_original=="PRO".or.name_original=="NPRO".or.name_original=="CPRO") then
		ip=16
	elseif(name_original=="CYS".or.name_original=="NCYS".or.name_original=="CCYS".or. &
	       name_original=="CYT".or.name_original=="NCYT".or.name_original=="CCYT") then
		ip=17
	elseif(name_original=="ALA".or.name_original=="NALA".or.name_original=="CALA") then
		ip=18
	elseif(name_original=="GLU".or.name_original=="NGLU".or.name_original=="CGLU".or. &
	       name_original=="GLH".or.name_original=="NGLH".or.name_original=="CGLH") then
		ip=19
	elseif(name_original=="ASP".or.name_original=="NASP".or.name_original=="CASP".or. &
	       name_original=="ASH".or.name_original=="NASH".or.name_original=="CASH") then
		ip=20
	elseif(name_original=="ACE") then
		ip=21
	elseif(name_original=="NME") then
		ip=22
	elseif(name_original=="NHE") then
		ip=23        
	endif

	if(ip==0) then
		open(10, file="error.txt", access="append")
			write(10,*) "findrotamer: unsupported residue name = ", name_original
		close(10)
		write(*,*) "findrotamer error: unsupported residue name =", name_original
		stop
	endif
	
	do ii=1, repeated_unit	
		rotanum=aa_lib(ip)%rotanum
		aa_group(ii,1)%cnum1=aa_lib(ip)%cnum1; aa_group(ii,1)%cnum2=aa_lib(ip)%cnum2; aa_group(ii,1)%cnum3=aa_lib(ip)%cnum3
		aa_group(ii,1)%gtype=name_original
		aa_group(ii,1)%atype1=aa_lib(ip)%atype1; aa_group(ii,1)%atype2=aa_lib(ip)%atype2; aa_group(ii,1)%atype3=aa_lib(ip)%atype3
		aa_group(ii,1)%coo1=aa_lib(ip)%coo1; aa_group(ii,1)%coo2=aa_lib(ip)%coo2; aa_group(ii,1)%coo3=aa_lib(ip)%coo3
	
		do k=1, group(ii,ic)%cnum1
			if(group(ii,ic)%atype1(k)=="N") then
				nr(1)=group(ii,ic)%coo1(k,1)
				nr(2)=group(ii,ic)%coo1(k,2)
				nr(3)=group(ii,ic)%coo1(k,3)
			elseif(group(ii,ic)%atype1(k)=="CA") then
				car(1)=group(ii,ic)%coo1(k,1)
				car(2)=group(ii,ic)%coo1(k,2)
				car(3)=group(ii,ic)%coo1(k,3)
			endif
		enddo
		do k=1, group(ii,ic)%cnum3
			if(group(ii,ic)%atype3(k)=="C") then
				cr(1)=group(ii,ic)%coo3(k,1)
				cr(2)=group(ii,ic)%coo3(k,2)
				cr(3)=group(ii,ic)%coo3(k,3)
			endif
		enddo

		call normalvector(nr, car, cr, r_norpep)

		r_nca(1)=nr(1)-car(1)
		r_nca(2)=nr(2)-car(2)
		r_nca(3)=nr(3)-car(3)

		do k=1, aa_group(ii,1)%cnum1
			if(aa_group(ii,1)%atype1(k)=="N") then
				aa_nr(1)=aa_group(ii,1)%coo1(k,1)
				aa_nr(2)=aa_group(ii,1)%coo1(k,2)
				aa_nr(3)=aa_group(ii,1)%coo1(k,3)
			elseif(aa_group(ii,1)%atype1(k)=="CA") then
				aa_car(1)=aa_group(ii,1)%coo1(k,1)
				aa_car(2)=aa_group(ii,1)%coo1(k,2)
				aa_car(3)=aa_group(ii,1)%coo1(k,3)
			endif
		enddo
		do k=1, aa_group(ii,1)%cnum3
			if(aa_group(ii,1)%atype3(k)=="C") then
				aa_cr(1)=aa_group(ii,1)%coo3(k,1)
				aa_cr(2)=aa_group(ii,1)%coo3(k,2)
				aa_cr(3)=aa_group(ii,1)%coo3(k,3)
			endif
		enddo	
		
		call normalvector(aa_nr, aa_car, aa_cr, r_norrot)

		call vectorrotation(r_norrot, r_norpep, m)
		
		temp1=matmul(aa_group(ii,1)%coo1, m)
		aa_group(ii,1)%coo1=temp1

		temp2=matmul(aa_group(ii,1)%coo2, m)
		aa_group(ii,1)%coo2=temp2

		temp3=matmul(aa_group(ii,1)%coo3, m)
		aa_group(ii,1)%coo3=temp3

		do k=1, aa_group(ii,1)%cnum1
			if(aa_group(ii,1)%atype1(k)=="N") then
				aa_nr(1)=aa_group(ii,1)%coo1(k,1)
				aa_nr(2)=aa_group(ii,1)%coo1(k,2)
				aa_nr(3)=aa_group(ii,1)%coo1(k,3)
			elseif(aa_group(ii,1)%atype1(k)=="CA") then
				aa_car(1)=aa_group(ii,1)%coo1(k,1)
				aa_car(2)=aa_group(ii,1)%coo1(k,2)
				aa_car(3)=aa_group(ii,1)%coo1(k,3)
			endif
		enddo

		aa_r_nca(1)=aa_nr(1)-aa_car(1)
		aa_r_nca(2)=aa_nr(2)-aa_car(2)
		aa_r_nca(3)=aa_nr(3)-aa_car(3)

		call vectorrotation(aa_r_nca, r_nca, m)

		temp1=matmul(aa_group(ii,1)%coo1, m)
		aa_group(ii,1)%coo1=temp1

		temp2=matmul(aa_group(ii,1)%coo2, m)
		aa_group(ii,1)%coo2=temp2

		temp3=matmul(aa_group(ii,1)%coo3, m)
		aa_group(ii,1)%coo3=temp3

		do k=1, aa_group(ii,1)%cnum1
			if(aa_group(ii,1)%atype1(k)=="CA") then
				aa_car(1)=aa_group(ii,1)%coo1(k,1)
				aa_car(2)=aa_group(ii,1)%coo1(k,2)
				aa_car(3)=aa_group(ii,1)%coo1(k,3)
			endif
		enddo

		r_trans(1)=car(1)-aa_car(1)
		r_trans(2)=car(2)-aa_car(2)
		r_trans(3)=car(3)-aa_car(3)

		do k=1, aa_group(ii,1)%cnum1
			aa_group(ii,1)%coo1(k,1)=anint((aa_group(ii,1)%coo1(k,1)+r_trans(1))*1000)/1000
			aa_group(ii,1)%coo1(k,2)=anint((aa_group(ii,1)%coo1(k,2)+r_trans(2))*1000)/1000
			aa_group(ii,1)%coo1(k,3)=anint((aa_group(ii,1)%coo1(k,3)+r_trans(3))*1000)/1000
		enddo
		do k=1, aa_group(ii,1)%cnum2
			aa_group(ii,1)%coo2(k,1)=anint((aa_group(ii,1)%coo2(k,1)+r_trans(1))*1000)/1000
			aa_group(ii,1)%coo2(k,2)=anint((aa_group(ii,1)%coo2(k,2)+r_trans(2))*1000)/1000
			aa_group(ii,1)%coo2(k,3)=anint((aa_group(ii,1)%coo2(k,3)+r_trans(3))*1000)/1000
		enddo
		do k=1, aa_group(ii,1)%cnum3
			aa_group(ii,1)%coo3(k,1)=anint((aa_group(ii,1)%coo3(k,1)+r_trans(1))*1000)/1000
			aa_group(ii,1)%coo3(k,2)=anint((aa_group(ii,1)%coo3(k,2)+r_trans(2))*1000)/1000
			aa_group(ii,1)%coo3(k,3)=anint((aa_group(ii,1)%coo3(k,3)+r_trans(3))*1000)/1000
		enddo	
	
		if(rotanum.le.1) goto 10

		grade_num=0
		if(aa_group(ii,1)%gtype=="VAL".or.aa_group(ii,1)%gtype=="NVAL".or.aa_group(ii,1)%gtype=="CVAL") then
			grade=1
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				else
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="LEU".or.aa_group(ii,1)%gtype=="NLEU".or.aa_group(ii,1)%gtype=="CLEU") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else	
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo		
		elseif(aa_group(ii,1)%gtype=="ILE".or.aa_group(ii,1)%gtype=="NILE".or.aa_group(ii,1)%gtype=="CILE") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2	
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3) 
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB".or.aa_group(ii,1)%atype2(k)=="CG2".or.aa_group(ii,1)%atype2(k)=="HG21".or.aa_group(ii,1)%atype2(k)=="HG22".or. &
					aa_group(ii,1)%atype2(k)=="HG23".or.aa_group(ii,1)%atype2(k)=="CG1") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG1") monitor(2)=grade_num(2)
				else	
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="PHE".or.aa_group(ii,1)%gtype=="NPHE".or.aa_group(ii,1)%gtype=="CPHE") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="TRP".or.aa_group(ii,1)%gtype=="NTRP".or.aa_group(ii,1)%gtype=="CTRP") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else	
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="TYR".or.aa_group(ii,1)%gtype=="NTYR".or.aa_group(ii,1)%gtype=="CTYR".or. &
			   aa_group(ii,1)%gtype=="TYX".or.aa_group(ii,1)%gtype=="NTYX".or.aa_group(ii,1)%gtype=="CTYX") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo						
		elseif(aa_group(ii,1)%gtype=="SER".or.aa_group(ii,1)%gtype=="NSER".or.aa_group(ii,1)%gtype=="CSER") then
			grade=1
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				else
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
				endif
			enddo		
		elseif(aa_group(ii,1)%gtype=="THR".or.aa_group(ii,1)%gtype=="NTHR".or.aa_group(ii,1)%gtype=="CTHR") then
			grade=1
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				else
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="CYS".or.aa_group(ii,1)%gtype=="NCYS".or.aa_group(ii,1)%gtype=="CCYS".or. &
			   aa_group(ii,1)%gtype=="CYT".or.aa_group(ii,1)%gtype=="NCYT".or.aa_group(ii,1)%gtype=="CCYT") then
			grade=1
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				else
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="MET".or.aa_group(ii,1)%gtype=="NMET".or.aa_group(ii,1)%gtype=="CMET") then
			grade=3
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="SD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="SD") monitor(3)=grade_num(3)
				else
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="ASN".or.aa_group(ii,1)%gtype=="NASN".or.aa_group(ii,1)%gtype=="CASN") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="GLN".or.aa_group(ii,1)%gtype=="NGLN".or.aa_group(ii,1)%gtype=="CGLN") then
			grade=3
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="CD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="CD") monitor(3)=grade_num(3)
				else	
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="ASP".or.aa_group(ii,1)%gtype=="NASP".or.aa_group(ii,1)%gtype=="CASP".or. &
			   aa_group(ii,1)%gtype=="ASH".or.aa_group(ii,1)%gtype=="NASH".or.aa_group(ii,1)%gtype=="CASH") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="GLU".or.aa_group(ii,1)%gtype=="NGLU".or.aa_group(ii,1)%gtype=="CGLU".or. &
			   aa_group(ii,1)%gtype=="GLH".or.aa_group(ii,1)%gtype=="NGLH".or.aa_group(ii,1)%gtype=="CGLH") then
			grade=3
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="CD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="CD") monitor(3)=grade_num(3)
				else
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="HIE".or.aa_group(ii,1)%gtype=="NHIE".or.aa_group(ii,1)%gtype=="CHIE".or. &
			   aa_group(ii,1)%gtype=="HIP".or.aa_group(ii,1)%gtype=="NHIP".or.aa_group(ii,1)%gtype=="CHIP") then
			grade=2
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				else
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="LYS".or.aa_group(ii,1)%gtype=="NLYS".or.aa_group(ii,1)%gtype=="CLYS".or. &
			   aa_group(ii,1)%gtype=="LYN".or.aa_group(ii,1)%gtype=="NLYN".or.aa_group(ii,1)%gtype=="CLYN") then
			grade=4
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="CD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="CD") monitor(3)=grade_num(3)
				elseif(aa_group(ii,1)%atype2(k)=="HD2".or.aa_group(ii,1)%atype2(k)=="HD3".or.aa_group(ii,1)%atype2(k)=="CE") then
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
					if(aa_group(ii,1)%atype2(k)=="CE") monitor(4)=grade_num(4)
				else
					grade_num(5)=grade_num(5)+1
					Iclass(5)%member(grade_num(5),1)=aa_group(ii,1)%coo2(k,1); Iclass(5)%member(grade_num(5),2)=aa_group(ii,1)%coo2(k,2); Iclass(5)%member(grade_num(5),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=5; index(k)%member_No=grade_num(5)
				endif
			enddo	
		elseif(aa_group(ii,1)%gtype=="ARG".or.aa_group(ii,1)%gtype=="NARG".or.aa_group(ii,1)%gtype=="CARG".or. &
			   aa_group(ii,1)%gtype=="ARN".or.aa_group(ii,1)%gtype=="NARN".or.aa_group(ii,1)%gtype=="CARN") then
			grade=4
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="CD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="CD") monitor(3)=grade_num(3)
				elseif(aa_group(ii,1)%atype2(k)=="HD2".or.aa_group(ii,1)%atype2(k)=="HD3".or.aa_group(ii,1)%atype2(k)=="NE") then
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
					if(aa_group(ii,1)%atype2(k)=="NE") monitor(4)=grade_num(4)
				else
					grade_num(5)=grade_num(5)+1
					Iclass(5)%member(grade_num(5),1)=aa_group(ii,1)%coo2(k,1); Iclass(5)%member(grade_num(5),2)=aa_group(ii,1)%coo2(k,2); Iclass(5)%member(grade_num(5),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=5; index(k)%member_No=grade_num(5)
				endif
			enddo
		elseif(aa_group(ii,1)%gtype=="PRO".or.aa_group(ii,1)%gtype=="NPRO".or.aa_group(ii,1)%gtype=="CPRO") then
			grade=3
			do k=1, aa_group(ii,1)%cnum2
				if(aa_group(ii,1)%atype2(k)=="CB") then
					grade_num(1)=grade_num(1)+1
					Iclass(1)%member(grade_num(1),1)=aa_group(ii,1)%coo2(k,1); Iclass(1)%member(grade_num(1),2)=aa_group(ii,1)%coo2(k,2); Iclass(1)%member(grade_num(1),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=1; index(k)%member_No=grade_num(1)
					monitor(1)=grade_num(1)
				elseif(aa_group(ii,1)%atype2(k)=="HB2".or.aa_group(ii,1)%atype2(k)=="HB3".or.aa_group(ii,1)%atype2(k)=="CG") then
					grade_num(2)=grade_num(2)+1
					Iclass(2)%member(grade_num(2),1)=aa_group(ii,1)%coo2(k,1); Iclass(2)%member(grade_num(2),2)=aa_group(ii,1)%coo2(k,2); Iclass(2)%member(grade_num(2),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=2; index(k)%member_No=grade_num(2)
					if(aa_group(ii,1)%atype2(k)=="CG") monitor(2)=grade_num(2)
				elseif(aa_group(ii,1)%atype2(k)=="HG2".or.aa_group(ii,1)%atype2(k)=="HG3".or.aa_group(ii,1)%atype2(k)=="CD") then
					grade_num(3)=grade_num(3)+1
					Iclass(3)%member(grade_num(3),1)=aa_group(ii,1)%coo2(k,1); Iclass(3)%member(grade_num(3),2)=aa_group(ii,1)%coo2(k,2); Iclass(3)%member(grade_num(3),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=3; index(k)%member_No=grade_num(3)
					if(aa_group(ii,1)%atype2(k)=="CD") monitor(3)=grade_num(3)
				else
					grade_num(4)=grade_num(4)+1
					Iclass(4)%member(grade_num(4),1)=aa_group(ii,1)%coo2(k,1); Iclass(4)%member(grade_num(4),2)=aa_group(ii,1)%coo2(k,2); Iclass(4)%member(grade_num(4),3)=aa_group(ii,1)%coo2(k,3)
					index(k)%class_No=4; index(k)%member_No=grade_num(4)
				endif
			enddo
		endif
			
		if(grade.ne.aa_lib(ip)%grade) then
			open(10, file="error.txt", access="append")
				write(10,*) aa_lib(ip)%gtype
				write(10,*) "grade=", grade
				write(10,*) "aa_lib(",ip,")%grade=", aa_lib(ip)%grade
				write(10,*) "They are not equal with each other!"
			close(10)
			stop
		endif

		do k=1, aa_group(ii,1)%cnum1
			if(aa_group(ii,1)%atype1(k)=="CA") then
				CA(1)=aa_group(ii,1)%coo1(k,1); CA(2)=aa_group(ii,1)%coo1(k,2); CA(3)=aa_group(ii,1)%coo1(k,3)
			endif
		enddo

		do n=2, rotanum
			aa_group(ii,n)%cnum1=aa_group(ii,1)%cnum1; aa_group(ii,n)%cnum2=aa_group(ii,1)%cnum2; aa_group(ii,n)%cnum3=aa_group(ii,1)%cnum3
			aa_group(ii,n)%gtype=name_original
			aa_group(ii,n)%atype1=aa_group(ii,1)%atype1; aa_group(ii,n)%atype2=aa_group(ii,1)%atype2; aa_group(ii,n)%atype3=aa_group(ii,1)%atype3
			aa_group(ii,n)%coo1=aa_group(ii,1)%coo1; aa_group(ii,n)%coo2=aa_group(ii,1)%coo2; aa_group(ii,n)%coo3=aa_group(ii,1)%coo3

			Tclass=Iclass	
			do j=1, grade
				delta_chi=real(aa_lib(ip)%dihedralangle(n,j)-aa_lib(ip)%dihedralangle(1,j))
				cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)			
				if(j==1) then
					rotaxis_x=Tclass(j)%member(monitor(j),1)-CA(1)
					rotaxis_y=Tclass(j)%member(monitor(j),2)-CA(2)
					rotaxis_z=Tclass(j)%member(monitor(j),3)-CA(3)
					rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				else
					rotaxis_x=Tclass(j)%member(monitor(j),1)-Tclass(j-1)%member(monitor(j-1),1)
					rotaxis_y=Tclass(j)%member(monitor(j),2)-Tclass(j-1)%member(monitor(j-1),2)
					rotaxis_z=Tclass(j)%member(monitor(j),3)-Tclass(j-1)%member(monitor(j-1),3)
					rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				endif

				call axisrotation(rotaxis, cos_angle, sin_angle, m)
				
				do l=(j+1), (grade+1)
					do k=1, grade_num(l)
						Tclass(l)%member(k,1)=Tclass(l)%member(k,1)-Tclass(j)%member(monitor(j),1)
						Tclass(l)%member(k,2)=Tclass(l)%member(k,2)-Tclass(j)%member(monitor(j),2)
						Tclass(l)%member(k,3)=Tclass(l)%member(k,3)-Tclass(j)%member(monitor(j),3)
					enddo
					
					Tmember=matmul(Tclass(l)%member, m)
					Tclass(l)%member=Tmember
					
					do k=1, grade_num(l)
						Tclass(l)%member(k,1)=anint((Tclass(l)%member(k,1)+Tclass(j)%member(monitor(j),1))*1000)/1000
						Tclass(l)%member(k,2)=anint((Tclass(l)%member(k,2)+Tclass(j)%member(monitor(j),2))*1000)/1000				
						Tclass(l)%member(k,3)=anint((Tclass(l)%member(k,3)+Tclass(j)%member(monitor(j),3))*1000)/1000
					enddo
				enddo
			enddo
				
			do k=1, aa_group(ii,n)%cnum2
				aa_group(ii,n)%coo2(k,1)=Tclass(index(k)%class_No)%member(index(k)%member_No,1)
				aa_group(ii,n)%coo2(k,2)=Tclass(index(k)%class_No)%member(index(k)%member_No,2)
				aa_group(ii,n)%coo2(k,3)=Tclass(index(k)%class_No)%member(index(k)%member_No,3)
			enddo
		enddo
10	continue
	enddo

	return
	end subroutine findrotamer

	
	subroutine energy_parameter(group, group_para)
	implicit none
	integer							:: i, j, k, cache_id
	character*4						:: igraph
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(energyparameters)			:: group_para(repeated_unit,gnum)

	do i=1, repeated_unit
		do j=1, gnum
			group_para(i,j)%atomid1=0; group_para(i,j)%atomid2=0; group_para(i,j)%atomid3=0
			group_para(i,j)%charge1=0.0; group_para(i,j)%charge2=0.0; group_para(i,j)%charge3=0.0
			group_para(i,j)%epsion1=0.0; group_para(i,j)%epsion2=0.0; group_para(i,j)%epsion3=0.0
			group_para(i,j)%r1=0.0; group_para(i,j)%r2=0.0; group_para(i,j)%r3=0.0
			group_para(i,j)%rborn1=0.0; group_para(i,j)%rborn2=0.0; group_para(i,j)%rborn3=0.0
			group_para(i,j)%fs1=0.0; group_para(i,j)%fs2=0.0; group_para(i,j)%fs3=0.0
			group_para(i,j)%dielecons1=0.0; group_para(i,j)%dielecons2=0.0; group_para(i,j)%dielecons3=0.0
			call ensure_forcefield_cache(group(i,j)%gtype, cache_id)
			do k=1, forcefield_cache(cache_id)%atom_count
				igraph=forcefield_cache(cache_id)%igraph(k)
				call fill_energy_parameter_from_cache(group(i,j), group_para(i,j), igraph, cache_id, k)
			enddo
		enddo
    enddo

	do i=1, repeated_unit
		do j=1, gnum
			do k=1, group(i,j)%cnum1
				if(group_para(i,j)%dielecons1(k)<=0.1) then
					open(10, file="error.txt", access="append")
						write(10,*) group(i,j)%gtype, i, j, "cluster 1 atom", k, group(i,j)%atype1(k), &
							"has wrong force field parameter", group_para(i,j)%dielecons1(k), "in the LIB!"
						write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
					close(10)
					stop
				endif
			enddo
			do k=1, group(i,j)%cnum2
				if(group_para(i,j)%dielecons2(k)<=0.1) then
					open(10, file="error.txt", access="append")
						write(10,*) group(i,j)%gtype, i, j, "cluster 2 atom", k, group(i,j)%atype2(k), &
							"has wrong force field parameter", group_para(i,j)%dielecons2(k), "in the LIB!"
						write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
					close(10)
					stop
				endif
			enddo
			do k=1, group(i,j)%cnum3
				if(group_para(i,j)%dielecons3(k)<=0.1) then
					open(10, file="error.txt", access="append")
						write(10,*) group(i,j)%gtype, i, j, "cluster 3 atom", k, group(i,j)%atype3(k), &
							"has wrong force field parameter", group_para(i,j)%dielecons3(k), "in the LIB!"
						write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
					close(10)
					stop
				endif
			enddo
		enddo
	enddo

	return
	end subroutine energy_parameter

	subroutine fill_energy_parameter_from_cache(group_residue, para_residue, igraph, cache_id, entry_id)
	implicit none
	integer							:: cache_id, entry_id, k
	character*4						:: igraph
	type(groupdetails)				:: group_residue
	type(energyparameters)			:: para_residue

	do k=1, group_residue%cnum1
		if(group_residue%atype1(k)==igraph) then
			para_residue%charge1(k)=forcefield_cache(cache_id)%charge(entry_id)
			para_residue%epsion1(k)=forcefield_cache(cache_id)%epsion(entry_id)
			para_residue%r1(k)=forcefield_cache(cache_id)%r(entry_id)
			para_residue%rborn1(k)=forcefield_cache(cache_id)%rborn(entry_id)
			para_residue%fs1(k)=forcefield_cache(cache_id)%fs(entry_id)
			para_residue%dielecons1(k)=forcefield_cache(cache_id)%dielecons(entry_id)
			para_residue%atomid1(k)=forcefield_cache(cache_id)%atomid(entry_id)
			return
		endif
	enddo
	do k=1, group_residue%cnum2
		if(group_residue%atype2(k)==igraph) then
			para_residue%charge2(k)=forcefield_cache(cache_id)%charge(entry_id)
			para_residue%epsion2(k)=forcefield_cache(cache_id)%epsion(entry_id)
			para_residue%r2(k)=forcefield_cache(cache_id)%r(entry_id)
			para_residue%rborn2(k)=forcefield_cache(cache_id)%rborn(entry_id)
			para_residue%fs2(k)=forcefield_cache(cache_id)%fs(entry_id)
			para_residue%dielecons2(k)=forcefield_cache(cache_id)%dielecons(entry_id)
			para_residue%atomid2(k)=forcefield_cache(cache_id)%atomid(entry_id)
			return
		endif
	enddo
	do k=1, group_residue%cnum3
		if(group_residue%atype3(k)==igraph) then
			para_residue%charge3(k)=forcefield_cache(cache_id)%charge(entry_id)
			para_residue%epsion3(k)=forcefield_cache(cache_id)%epsion(entry_id)
			para_residue%r3(k)=forcefield_cache(cache_id)%r(entry_id)
			para_residue%rborn3(k)=forcefield_cache(cache_id)%rborn(entry_id)
			para_residue%fs3(k)=forcefield_cache(cache_id)%fs(entry_id)
			para_residue%dielecons3(k)=forcefield_cache(cache_id)%dielecons(entry_id)
			para_residue%atomid3(k)=forcefield_cache(cache_id)%atomid(entry_id)
			return
		endif
	enddo

	return
	end subroutine fill_energy_parameter_from_cache

	
	subroutine atom_links(group, numex, inb, numex4, inb4)
	implicit none
	type atomlink      ! The Data type "atomlink" is used to store the neighboring atoms for each atom.
		integer				:: linknum
		integer				:: linkindex(4)
	end type

	integer							:: categoryID, i, ii, ic, j, k, i1, j1, i2, j2, atomid, natom, cache_id
	integer							:: ipres, numex(repeated_unit*atom_num), inb(repeated_unit*atom_num,20), numex4(repeated_unit*atom_num), inb4(repeated_unit*atom_num,60)

	type(groupdetails)				:: group(repeated_unit,gnum)
	type(atomlink)					:: atom(repeated_unit*atom_num)									

	do i=1, repeated_unit*atom_num
		atom(i)%linknum=0
		atom(i)%linkindex=0
	enddo
	natom=0
	do categoryID=1, num4category
		do i=1, selfassembly(categoryID)%num4peptides ! Number of peptides in this energy group.
			ii=selfassembly(categoryID)%peptideID(i) ! Peptide chain ID.
			do ic=1, gnum
				ipres=natom
				call ensure_atomlink_cache(group(ii,ic)%gtype, cache_id)
				do j=1, atomlink_cache(cache_id)%natom
					atomid=ipres+j
					atom(atomid)%linknum=atomlink_cache(cache_id)%linknum(j)
					do k=1, atom(atomid)%linknum
						atom(atomid)%linkindex(k)=ipres+atomlink_cache(cache_id)%linkindex(j,k)
					enddo
				enddo
				natom=ipres+atomlink_cache(cache_id)%natom
			enddo
		enddo
	enddo

	do i1=1, natom
		numex(i1)=0
		do j1=1, atom(i1)%linknum
			numex(i1)=numex(i1)+1
			inb(i1,numex(i1))=atom(i1)%linkindex(j1)
		enddo

		do j1=1, atom(i1)%linknum
			i2=atom(i1)%linkindex(j1)
			do j2=1, atom(i2)%linknum
				if(atom(i2)%linkindex(j2).eq.i1) goto 40
				do k=1, atom(i1)%linknum
					if(atom(i2)%linkindex(j2).eq.atom(i1)%linkindex(k)) goto 40
				enddo
				numex(i1)=numex(i1)+1
				inb(i1,numex(i1))=atom(i2)%linkindex(j2)
40				continue
			enddo
		enddo
	enddo

	do i1=1, natom
		numex4(i1)=0
		do j1=1, numex(i1)
			i2=inb(i1,j1)
			do j2=1, atom(i2)%linknum
				if(atom(i2)%linkindex(j2).eq.i1) goto 50
				do k=1, numex(i1)
					if(atom(i2)%linkindex(j2).eq.inb(i1,k)) goto 50
				enddo
				numex4(i1)=numex4(i1)+1
				inb4(i1,numex4(i1))=atom(i2)%linkindex(j2)
50			continue
			enddo
		enddo
	enddo

	return
	end subroutine atom_links

	
	subroutine atom_links4sidechain(chainID, ic, group, S_numex, S_inb, S_numex4, S_inb4)
	implicit none
	type atomlink
		integer				:: linknum
		integer				:: linkindex(4)
	end type
	
	integer							:: chainID, i, ii, j, k, i1, j1, i2, j2, natom, cache_id
	integer							:: ic
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(atomlink)					:: atom(50)

	ii=chainID

	do i=1, 50
		atom(i)%linknum=0
		atom(i)%linkindex=0
	enddo
	call ensure_atomlink_cache(group(ii,ic)%gtype, cache_id)
	do i=1, atomlink_cache(cache_id)%natom
		atom(i)%linknum=atomlink_cache(cache_id)%linknum(i)
		do j=1, atom(i)%linknum
			atom(i)%linkindex(j)=atomlink_cache(cache_id)%linkindex(i,j)
		enddo
	enddo
	natom=atomlink_cache(cache_id)%natom

	do i1=1, natom
		S_numex(i1)=0
		do j1=1, atom(i1)%linknum
			S_numex(i1)=S_numex(i1)+1						! Count the bonded atoms for atom i1.
			S_inb(i1,S_numex(i1))=atom(i1)%linkindex(j1)    
        enddo

		do j1=1, atom(i1)%linknum
			i2=atom(i1)%linkindex(j1)

			if(i2.gt.natom) atom(i2)%linknum=0     ! Reset link data beyond the current atom count.
            
			if(i2 > 0) then
				do j2=1, atom(i2)%linknum
					if(atom(i2)%linkindex(j2).eq.i1) goto 40
					do k=1, atom(i1)%linknum
						if(atom(i2)%linkindex(j2).eq.atom(i1)%linkindex(k)) goto 40
					enddo
					S_numex(i1)=S_numex(i1)+1
					S_inb(i1,S_numex(i1))=atom(i2)%linkindex(j2)
40					continue
                enddo
            end if
		enddo
	enddo

	do i1=1, natom
		S_numex4(i1)=0
		do j1=1, S_numex(i1)
			i2=S_inb(i1,j1)
            
            if (i2 > 0) then
				do j2=1, atom(i2)%linknum
					if(atom(i2)%linkindex(j2).eq.i1) goto 50
					do k=1, S_numex(i1)
						if(atom(i2)%linkindex(j2).eq.S_inb(i1,k)) goto 50
					enddo
					S_numex4(i1)=S_numex4(i1)+1
					S_inb4(i1,S_numex4(i1))=atom(i2)%linkindex(j2)
50				continue
                enddo
            end if
		enddo
	enddo

	return
    end subroutine atom_links4sidechain


	subroutine find_rest_aa(group,rest_AA,l)
	implicit none
	integer :: i,j,l
	character(len=4) :: AA
	character(len=4) :: rest_AA(maximum_nmr_site_num)
	logical :: pool_used(maximum_nmr_site_num)
	type(groupdetails) :: group(repeated_unit,gnum)

	rest_AA=''
	pool_used=.false.
    do j = 1, nmr_site_num ! Record the amino-acid identities at all grouped constrained sites.
		AA=adjustl(group(1,nmr_site_ID(j))%gtype)
		if (len_trim(AA) == 4) then
			if (AA(1:1) == 'N' .or. AA(1:1) == 'C') AA=AA(2:4)
		endif
		do i=1,NMR_pool_size
			if (.not.pool_used(i) .and. trim(AA) == trim(NMR_AA_pool(i))) then
				pool_used(i)=.true.
				exit
			endif
		enddo
    enddo

	l=0
	do i=1,NMR_pool_size
		if (pool_used(i)) cycle
		l=l+1
		rest_AA(l)=NMR_AA_pool(i)
	enddo
	end subroutine find_rest_aa
    
    
	subroutine mc_choose_aminoacid(ic,group,aminoacid_name)
	implicit none
	integer, intent(in) :: ic
	type(groupdetails), intent(in) :: group(repeated_unit,gnum)
	character(len=4), intent(out) :: aminoacid_name
	integer :: i,j,k,current_category,candidate_category,choice
	integer :: category_count(4), category_min(4), category_max(4), trial_count(4)
	integer :: specific_count, available_count, candidate_count, available_category(20),rest_count
	character(len=4) :: current_name,current_base,sequence_name,candidate
	character(len=4) :: available_aa(20),candidates(20),rest_AA(maximum_nmr_site_num)
	logical :: blocked,grouped_site
	real :: ran2

	! Default to no change. A legal alternative replaces this value below.
	aminoacid_name=adjustl(group(1,ic)%gtype)

	! Count the complete first-chain composition. Caps are excluded, while
	! fixed and grouped positional sites are included.
	call count_amino_acid_categories(group(1,:)%gtype,category_count(1), &
		category_count(2),category_count(3),category_count(4))
	category_min=(/fpho_min,fpol_min,fchg_min,foth_min/)
	category_max=(/fpho_max,fpol_max,fchg_max,foth_max/)

	! Count every specifically restricted amino acid once for this sequence.
	if (n_restrictions > 0) aa_restrictions(1:n_restrictions)%count=0
	if (n_morethan > 0) aa_morethan(1:n_morethan)%count=0
	do i=1,gnum
		sequence_name=adjustl(group(1,i)%gtype)
		if (trim(sequence_name) == 'ACE' .or. trim(sequence_name) == 'NME' .or. trim(sequence_name) == 'NHE') cycle
		if (len_trim(sequence_name) == 4) then
			if (sequence_name(1:1) == 'N' .or. sequence_name(1:1) == 'C') sequence_name=sequence_name(2:4)
		endif
		do j=1,n_restrictions
			if (trim(sequence_name) == trim(aa_restrictions(j)%gtype)) aa_restrictions(j)%count=aa_restrictions(j)%count+1
		enddo
		do j=1,n_morethan
			if (trim(sequence_name) == trim(aa_morethan(j)%gtype)) aa_morethan(j)%count=aa_morethan(j)%count+1
		enddo
	enddo

	! A grouped site may use only an available entry from its grouped pool.
	! A regular mutable site may use the complete amino-acid library.
	available_aa=''
	available_category=0
	available_count=0
	grouped_site=.false.
	do k=1,nmr_site_num
		if (ic == nmr_site_ID(k)) then
			grouped_site=.true.
			exit
		endif
	enddo

	if (grouped_site) then
		call find_rest_aa(group,rest_AA,rest_count)
		do i=1,rest_count
			candidate=rest_AA(i)
			available_count=available_count+1
			available_aa(available_count)=candidate
			select case(trim(candidate))
			case('ALA','LEU','VAL','ILE','MET','PHE','TYR','TRP','GLY')
				available_category(available_count)=1
			case('ASN','GLN','SER','THR','HIE')
				available_category(available_count)=2
			case('ARG','LYS','GLU','ASP')
				available_category(available_count)=3
			case('PRO','CYS')
				available_category(available_count)=4
			end select
		enddo
	else
		do i=1,hydrationprop%hnum
			available_count=available_count+1
			available_aa(available_count)=hydrationprop%hgtype(i)
			available_category(available_count)=1
		enddo
		do i=1,hydrationprop%pnum
			available_count=available_count+1
			available_aa(available_count)=hydrationprop%pgtype(i)
			available_category(available_count)=2
		enddo
		do i=1,hydrationprop%cnum
			available_count=available_count+1
			available_aa(available_count)=hydrationprop%cgtype(i)
			available_category(available_count)=3
		enddo
		do i=1,hydrationprop%onum
			available_count=available_count+1
			available_aa(available_count)=hydrationprop%ogtype(i)
			available_category(available_count)=4
		enddo
	endif

	current_name=adjustl(group(1,ic)%gtype)
	current_base=current_name
	if (len_trim(current_base) == 4) then
		if (current_base(1:1) == 'N' .or. current_base(1:1) == 'C') current_base=current_base(2:4)
	endif

	select case(trim(current_base))
	case('ALA','LEU','VAL','ILE','MET','PHE','TYR','TRP','GLY')
		current_category=1
	case('ASN','GLN','SER','THR','HIE','HIS','HID','HIP')
		current_category=2
	case('ARG','LYS','GLU','ASP')
		current_category=3
	case('PRO','CYS')
		current_category=4
	case default
		write(*,'(3A)') 'ERROR: cannot mutate unknown residue ',trim(current_base),'.'
		stop
	end select

	! Do not mutate away an amino acid that is already at its specific minimum.
	do j=1,n_morethan
		if (trim(current_base) == trim(aa_morethan(j)%gtype) .and. aa_morethan(j)%count <= aa_morethan(j)%min) return
	enddo

	! Build one combined library. Every entry in this library will have the
	! same probability of being selected.
	candidates=''
	candidate_count=0
	do i=1,available_count
		candidate=available_aa(i)
		if (trim(candidate) == trim(current_base)) cycle

		blocked=.false.
		do j=1,n_restrictions
			if (trim(candidate) == trim(aa_restrictions(j)%gtype)) then
				specific_count=aa_restrictions(j)%count
				if (specific_count >= aa_restrictions(j)%max) blocked=.true.
				exit
			endif
		enddo
		if (blocked) cycle

		candidate_category=available_category(i)
		trial_count=category_count
		trial_count(current_category)=trial_count(current_category)-1
		trial_count(candidate_category)=trial_count(candidate_category)+1
		if (all(trial_count >= category_min) .and. all(trial_count <= category_max)) then
			candidate_count=candidate_count+1
			candidates(candidate_count)=candidate
		endif
	enddo
	if (candidate_count == 0) return

	call ran_gen(ran2,0)
	choice=int(ran2*candidate_count-1.0e-3)+1
	if (choice > candidate_count) choice=candidate_count
	aminoacid_name=candidates(choice)

	! Preserve an N- or C-terminal residue prefix.
	if (len_trim(current_name) == 4) then
		if (current_name(1:1) == 'N' .or. current_name(1:1) == 'C') aminoacid_name=current_name(1:1)//trim(aminoacid_name)
	endif
	end subroutine mc_choose_aminoacid

	subroutine check_restrained_aminoacid(group, ic, aminoacid_name, feedback)
    implicit none
    integer					:: ic, feedback, i, j 
    integer                 :: new_count, old_count                      ! restraint counter
    character*4             :: aminoacid_name, AA, AA_C, AA_N
    character*4             :: new_pep_sequence(gnum), old_pep_sequence(gnum)
    type(groupdetails)      :: group(repeated_unit,gnum)
    
    new_count=0; old_count=0
    if (NMR_restraint_AA(1)  == "None" .or. restraint_num == 0) then
        feedback=1
    else   
    
        do i=1, gnum
            old_pep_sequence(i)=group(1,i)%gtype
        enddo
        
        new_pep_sequence = old_pep_sequence
		new_pep_sequence(ic) = aminoacid_name
    
		do i = 1, restraint_num
            AA   =  NMR_restraint_AA(i)
            AA_C =  "C" // AA
		    AA_N =  "N" // AA
			do j = 1, gnum
				if (new_pep_sequence(j) == AA .or. new_pep_sequence(j) == AA_C .or. new_pep_sequence(j) == AA_N ) then
                    new_count = new_count+1
					exit ! if there is one resraint, exit the loop
                endif
            enddo
            
			do j = 1, gnum
                if (old_pep_sequence(j) == AA .or. old_pep_sequence(j) == AA_C .or. old_pep_sequence(j) == AA_N ) then
                    old_count = old_count+1
					exit ! if there is one resraint, exit the loop
                endif
			enddo
        enddo
    endif
    
	if (new_count == restraint_num) then
		feedback=1        ! all restrained aa were found in the new sequence
	else if (new_count > old_count) then
		feedback=1 
	endif
    
    return
    end subroutine check_restrained_aminoacid
    
	subroutine groupinfo(name, group_name, flag)
	implicit none
	integer					:: i, flag
	character*4				:: name, group_name(3)

	if(name=="GLY".or.name=="NGLY".or.name=="CGLY") then
		group_name(1)="GLY"
		group_name(2)="NGLY"
		group_name(3)="CGLY"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="LEU".or.name=="NLEU".or.name=="CLEU") then
		group_name(1)="LEU"
		group_name(2)="NLEU"
		group_name(3)="CLEU"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="VAL".or.name=="NVAL".or.name=="CVAL") then
		group_name(1)="VAL"
		group_name(2)="NVAL"
		group_name(3)="CVAL"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ILE".or.name=="NILE".or.name=="CILE") then
		group_name(1)="ILE"
		group_name(2)="NILE"
		group_name(3)="CILE"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="MET".or.name=="NMET".or.name=="CMET") then
		group_name(1)="MET"
		group_name(2)="NMET"
		group_name(3)="CMET"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="PHE".or.name=="NPHE".or.name=="CPHE") then
		group_name(1)="PHE"
		group_name(2)="NPHE"
		group_name(3)="CPHE"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="TYR".or.name=="NTYR".or.name=="CTYR") then
		group_name(1)="TYR"
		group_name(2)="NTYR"
		group_name(3)="CTYR"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="TYX".or.name=="NTYX".or.name=="CTYX") then
		group_name(1)="TYX"
		group_name(2)="NTYX"
		group_name(3)="CTYX"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="TRP".or.name=="NTRP".or.name=="CTRP") then
		group_name(1)="TRP"
		group_name(2)="NTRP"
		group_name(3)="CTRP"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ARG".or.name=="NARG".or.name=="CARG") then
		group_name(1)="ARG"
		group_name(2)="NARG"
		group_name(3)="CARG"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ARN".or.name=="NARN".or.name=="CARN") then
		group_name(1)="ARN"
		group_name(2)="NARN"
		group_name(3)="CARN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="LYN".or.name=="NLYN".or.name=="CLYN") then
		group_name(1)="LYN"
		group_name(2)="NLYN"
		group_name(3)="CLYN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="LYS".or.name=="NLYS".or.name=="CLYS") then
		group_name(1)="LYS"
		group_name(2)="NLYS"
		group_name(3)="CLYS"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="SER".or.name=="NSER".or.name=="CSER") then
		group_name(1)="SER"
		group_name(2)="NSER"
		group_name(3)="CSER"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="THR".or.name=="NTHR".or.name=="CTHR") then
		group_name(1)="THR"
		group_name(2)="NTHR"
		group_name(3)="CTHR"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ASN".or.name=="NASN".or.name=="CASN") then
		group_name(1)="ASN"
		group_name(2)="NASN"
		group_name(3)="CASN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="GLN".or.name=="NGLN".or.name=="CGLN") then
		group_name(1)="GLN"
		group_name(2)="NGLN"
		group_name(3)="CGLN"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="HIE".or.name=="NHIE".or.name=="CHIE") then
		group_name(1)="HIE"
		group_name(2)="NHIE"
		group_name(3)="CHIE"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="HIP".or.name=="NHIP".or.name=="CHIP") then
		group_name(1)="HIP"
		group_name(2)="NHIP"
		group_name(3)="CHIP"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="PRO".or.name=="NPRO".or.name=="CPRO") then
		group_name(1)="PRO"
		group_name(2)="NPRO"
		group_name(3)="CPRO"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="CYS".or.name=="NCYS".or.name=="CCYS") then
		group_name(1)="CYS"
		group_name(2)="NCYS"
		group_name(3)="CCYS"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="CYT".or.name=="NCYT".or.name=="CCYT") then
		group_name(1)="CYT"
		group_name(2)="NCYT"
		group_name(3)="CCYT"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ALA".or.name=="NALA".or.name=="CALA") then
		group_name(1)="ALA"
		group_name(2)="NALA"
		group_name(3)="CALA"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="GLH".or.name=="NGLH".or.name=="CGLH") then
		group_name(1)="GLH"
		group_name(2)="NGLH"
		group_name(3)="CGLH"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="GLU".or.name=="NGLU".or.name=="CGLU") then
		group_name(1)="GLU"
		group_name(2)="NGLU"
		group_name(3)="CGLU"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ASH".or.name=="NASH".or.name=="CASH") then
		group_name(1)="ASH"
		group_name(2)="NASH"
		group_name(3)="CASH"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	elseif(name=="ASP".or.name=="NASP".or.name=="CASP") then
		group_name(1)="ASP"
		group_name(2)="NASP"
		group_name(3)="CASP"
		do i=1, 3
			if(name==group_name(i)) then
				flag=i
				goto 5
			endif
		enddo
	endif
5	continue

	return
	end subroutine groupinfo

	
    
    
	subroutine sidechain_category(chainID, ic, group, Iclass, grade, grade_num, index, monitor)
	implicit none
	integer								:: grade, grade_num(6), monitor(6), chainID, i, ii, ic
	type(groupdetails)					:: group(repeated_unit,gnum)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6)

	ii=chainID
	grade_num=0
	if(group(ii,ic)%gtype=="VAL".or.group(ii,ic)%gtype=="NVAL".or.group(ii,ic)%gtype=="CVAL") then
		grade=1
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1)
				Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2)
				Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
			endif
		enddo				
	elseif(group(ii,ic)%gtype=="LEU".or.group(ii,ic)%gtype=="NLEU".or.group(ii,ic)%gtype=="CLEU") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else	
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo		
	elseif(group(ii,ic)%gtype=="ILE".or.group(ii,ic)%gtype=="NILE".or.group(ii,ic)%gtype=="CILE") then
		grade=2
		do i=1, group(ii,ic)%cnum2	
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3) 
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB".or.group(ii,ic)%atype2(i)=="CG2".or.group(ii,ic)%atype2(i)=="HG21".or.group(ii,ic)%atype2(i)=="HG22".or. &
			       group(ii,ic)%atype2(i)=="HG23".or.group(ii,ic)%atype2(i)=="CG1") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG1") monitor(2)=grade_num(2)
			else	
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo		
	elseif(group(ii,ic)%gtype=="PHE".or.group(ii,ic)%gtype=="NPHE".or.group(ii,ic)%gtype=="CPHE") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else	
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="TRP".or.group(ii,ic)%gtype=="NTRP".or.group(ii,ic)%gtype=="CTRP") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else	
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="TYR".or.group(ii,ic)%gtype=="NTYR".or.group(ii,ic)%gtype=="CTYR".or.  &
	       group(ii,ic)%gtype=="TYX".or.group(ii,ic)%gtype=="NTYX".or.group(ii,ic)%gtype=="CTYX") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="SER".or.group(ii,ic)%gtype=="NSER".or.group(ii,ic)%gtype=="CSER") then
		grade=1
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="THR".or.group(ii,ic)%gtype=="NTHR".or.group(ii,ic)%gtype=="CTHR") then
		grade=1
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="CYS".or.group(ii,ic)%gtype=="NCYS".or.group(ii,ic)%gtype=="CCYS".or.  &
	       group(ii,ic)%gtype=="CYT".or.group(ii,ic)%gtype=="NCYT".or.group(ii,ic)%gtype=="CCYT") then
		grade=1
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			else
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
			endif
		enddo
	elseif(group(ii,ic)%gtype=="MET".or.group(ii,ic)%gtype=="NMET".or.group(ii,ic)%gtype=="CMET") then
		grade=3
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ii,ic)%atype2(i)=="HG2".or.group(ii,ic)%atype2(i)=="HG3".or.group(ii,ic)%atype2(i)=="SD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ii,ic)%atype2(i)=="SD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ii,ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ii,ic)%coo2(i,2); Iclass(4)%member(grade_num(4),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="ASN".or.group(ii,ic)%gtype=="NASN".or.group(ii,ic)%gtype=="CASN") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="GLN".or.group(ii,ic)%gtype=="NGLN".or.group(ii,ic)%gtype=="CGLN") then
		grade=3
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ii,ic)%atype2(i)=="HG2".or.group(ii,ic)%atype2(i)=="HG3".or.group(ii,ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ii,ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			else	
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ii,ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ii,ic)%coo2(i,2); Iclass(4)%member(grade_num(4),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo	
	elseif(group(ii,ic)%gtype=="ASP".or.group(ii,ic)%gtype=="NASP".or.group(ii,ic)%gtype=="CASP".or.  &
	       group(ii,ic)%gtype=="ASH".or.group(ii,ic)%gtype=="NASH".or.group(ii,ic)%gtype=="CASH") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ii,ic)%gtype=="GLU".or.group(ii,ic)%gtype=="NGLU".or.group(ii,ic)%gtype=="CGLU".or.  &
	       group(ii,ic)%gtype=="GLH".or.group(ii,ic)%gtype=="NGLH".or.group(ii,ic)%gtype=="CGLH") then
		grade=3
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ii,ic)%atype2(i)=="HG2".or.group(ii,ic)%atype2(i)=="HG3".or.group(ii,ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ii,ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			else
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ii,ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ii,ic)%coo2(i,2); Iclass(4)%member(grade_num(4),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
			endif
		enddo
	elseif(group(ii,ic)%gtype=="HIE".or.group(ii,ic)%gtype=="NHIE".or.group(ii,ic)%gtype=="CHIE".or.  &
	       group(ii,ic)%gtype=="HIP".or.group(ii,ic)%gtype=="NHIP".or.group(ii,ic)%gtype=="CHIP") then
		grade=2
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			else
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
			endif
		enddo
	elseif(group(ii,ic)%gtype=="LYS".or.group(ii,ic)%gtype=="NLYS".or.group(ii,ic)%gtype=="CLYS".or.  &
	       group(ii,ic)%gtype=="LYN".or.group(ii,ic)%gtype=="NLYN".or.group(ii,ic)%gtype=="CLYN") then
		grade=4
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ii,ic)%atype2(i)=="HG2".or.group(ii,ic)%atype2(i)=="HG3".or.group(ii,ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ii,ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			elseif(group(ii,ic)%atype2(i)=="HD2".or.group(ii,ic)%atype2(i)=="HD3".or.group(ii,ic)%atype2(i)=="CE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ii,ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ii,ic)%coo2(i,2); Iclass(4)%member(grade_num(4),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				if(group(ii,ic)%atype2(i)=="CE") monitor(4)=grade_num(4)
			else
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=group(ii,ic)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=group(ii,ic)%coo2(i,2); Iclass(5)%member(grade_num(5),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
			endif
		enddo
	elseif(group(ii,ic)%gtype=="ARG".or.group(ii,ic)%gtype=="NARG".or.group(ii,ic)%gtype=="CARG".or.  &
	       group(ii,ic)%gtype=="ARN".or.group(ii,ic)%gtype=="NARN".or.group(ii,ic)%gtype=="CARN") then
		grade=4
		do i=1, group(ii,ic)%cnum2
			if(group(ii,ic)%atype2(i)=="CB") then
				grade_num(1)=grade_num(1)+1
				Iclass(1)%member(grade_num(1),1)=group(ii,ic)%coo2(i,1); Iclass(1)%member(grade_num(1),2)=group(ii,ic)%coo2(i,2); Iclass(1)%member(grade_num(1),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=1; index(i)%member_No=grade_num(1)
				monitor(1)=grade_num(1)
			elseif(group(ii,ic)%atype2(i)=="HB2".or.group(ii,ic)%atype2(i)=="HB3".or.group(ii,ic)%atype2(i)=="CG") then
				grade_num(2)=grade_num(2)+1
				Iclass(2)%member(grade_num(2),1)=group(ii,ic)%coo2(i,1); Iclass(2)%member(grade_num(2),2)=group(ii,ic)%coo2(i,2); Iclass(2)%member(grade_num(2),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=2; index(i)%member_No=grade_num(2)
				if(group(ii,ic)%atype2(i)=="CG") monitor(2)=grade_num(2)
			elseif(group(ii,ic)%atype2(i)=="HG2".or.group(ii,ic)%atype2(i)=="HG3".or.group(ii,ic)%atype2(i)=="CD") then
				grade_num(3)=grade_num(3)+1
				Iclass(3)%member(grade_num(3),1)=group(ii,ic)%coo2(i,1); Iclass(3)%member(grade_num(3),2)=group(ii,ic)%coo2(i,2); Iclass(3)%member(grade_num(3),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=3; index(i)%member_No=grade_num(3)
				if(group(ii,ic)%atype2(i)=="CD") monitor(3)=grade_num(3)
			elseif(group(ii,ic)%atype2(i)=="HD2".or.group(ii,ic)%atype2(i)=="HD3".or.group(ii,ic)%atype2(i)=="NE") then
				grade_num(4)=grade_num(4)+1
				Iclass(4)%member(grade_num(4),1)=group(ii,ic)%coo2(i,1); Iclass(4)%member(grade_num(4),2)=group(ii,ic)%coo2(i,2); Iclass(4)%member(grade_num(4),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=4; index(i)%member_No=grade_num(4)
				if(group(ii,ic)%atype2(i)=="NE") monitor(4)=grade_num(4)
			else
				grade_num(5)=grade_num(5)+1
				Iclass(5)%member(grade_num(5),1)=group(ii,ic)%coo2(i,1); Iclass(5)%member(grade_num(5),2)=group(ii,ic)%coo2(i,2); Iclass(5)%member(grade_num(5),3)=group(ii,ic)%coo2(i,3)
				index(i)%class_No=5; index(i)%member_No=grade_num(5)
			endif
		enddo
	endif

	return
	end subroutine sidechain_category	

	
	subroutine dihedralangle_reading(gtype, dihedral_num, dihedral)
	implicit none
	integer								:: dihedral_num, cache_id
	character*4							:: gtype
	type(dihedralparameters)			:: dihedral	

	call ensure_dihedral_cache(gtype, cache_id)
	dihedral_num=dihedral_cache(cache_id)%dihedral_num
	dihedral=dihedral_cache(cache_id)%dihedral
	
	return
    end subroutine dihedralangle_reading	

    subroutine convert_AA_name(group, pep_name)
    implicit none
    type(groupdetails)				:: group(repeated_unit,gnum)
    integer							:: i
    character*4					    :: name
    character(len=gnum)				:: pep_name

    do i = 1, gnum
        name = group(1,i)%gtype
		if(name=="GLY".or.name=="NGLY".or.name=="CGLY") then
            pep_name(i:i) = "G"
        elseif(name=="ALA".or.name=="NALA".or.name=="CALA") then
            pep_name(i:i) = "A"
        elseif(name=="LEU".or.name=="NLEU".or.name=="CLEU") then
            pep_name(i:i) = "L"
        elseif(name=="VAL".or.name=="NVAL".or.name=="CVAL") then
            pep_name(i:i) = "V"
        elseif(name=="ILE".or.name=="NILE".or.name=="CILE") then
            pep_name(i:i) = "I"
        elseif(name=="MET".or.name=="NMET".or.name=="CMET") then
            pep_name(i:i) = "M"
        elseif(name=="PHE".or.name=="NPHE".or.name=="CPHE") then
            pep_name(i:i) = "F"
        elseif(name=="TYR".or.name=="NTYR".or.name=="CTYR".or.name=="TYX".or.name=="NTYX".or.name=="CTYX") then
            pep_name(i:i) = "Y"
        elseif(name=="TRP".or.name=="NTRP".or.name=="CTRP") then
            pep_name(i:i) = "W"
        elseif(name=="ARG".or.name=="NARG".or.name=="CARG".or.name=="ARN".or.name=="NARN".or.name=="CARN") then
            pep_name(i:i) = "R"
        elseif(name=="LYS".or.name=="NLYS".or.name=="CLYS".or.name=="LYN".or.name=="NLYN".or.name=="CLYN") then
            pep_name(i:i) = "K"
        elseif(name=="SER".or.name=="NSER".or.name=="CSER") then
            pep_name(i:i) = "S"
        elseif(name=="THR".or.name=="NTHR".or.name=="CTHR") then
            pep_name(i:i) = "T"
        elseif(name=="ASN".or.name=="NASN".or.name=="CASN") then
            pep_name(i:i) = "N"
        elseif(name=="GLN".or.name=="NGLN".or.name=="CGLN") then
            pep_name(i:i) = "Q"
        elseif(name=="HIE".or.name=="NHIE".or.name=="CHIE".or.name=="HIP".or.name=="NHIP".or.name=="CHIP") then
            pep_name(i:i) = "H"
        elseif(name=="GLU".or.name=="NGLU".or.name=="CGLU".or.name=="GLH".or.name=="NGLH".or.name=="CGLH") then
            pep_name(i:i) = "E"
        elseif(name=="ASP".or.name=="NASP".or.name=="CASP".or.name=="ASH".or.name=="NASH".or.name=="CASH") then
            pep_name(i:i) = "D"
        elseif(name=="CYS".or.name=="NCYS".or.name=="CCYS".or.name=="CYT".or.name=="NCYT".or.name=="CCYT") then
            pep_name(i:i) = "C"
        elseif(name=="PRO".or.name=="NPRO".or.name=="CPRO") then
            pep_name(i:i) = "P"
        else 
            pep_name(i:i) = " "
        end if
    end do

    return
    end subroutine convert_AA_name
    
end module database

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module transplant
	
	use constant
	use mathfunction
	use randomgenerator

	contains
	subroutine residue_replace(chainID, ic, group, ip, aa_group, temp_group)
	implicit none
	integer								:: chainID, ic, ip, ii, j, k, l, flag, ha2_pos
	type(groupdetails)					:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), aa_group(repeated_unit,40)

	temp_group=group

	ii=chainID
	if(temp_group(ii,ic)%gtype=="PRO".or.temp_group(ii,ic)%gtype=="NPRO".or.temp_group(ii,ic)%gtype=="CPRO") then
		if(aa_group(ii,ip)%gtype=="PRO".or.aa_group(ii,ip)%gtype=="NPRO".or.aa_group(ii,ip)%gtype=="CPRO") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		elseif(aa_group(ii,ip)%gtype=="GLY".or.aa_group(ii,ip)%gtype=="NGLY".or.aa_group(ii,ip)%gtype=="CGLY") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			flag=0
			do j=1, temp_group(ii,ic)%cnum1
				if(temp_group(ii,ic)%atype1(j)=="H2".or.temp_group(ii,ic)%atype1(j)=="H3") then
					flag=1
				endif
			enddo

			do j=1, aa_group(ii,ip)%cnum1
				if(aa_group(ii,ip)%atype1(j)=="H") then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="N") then
							if(flag==1) then ! Use H1 instead of H for an N-terminal residue.
								temp_group(ii,ic)%atype1(k+1)="H1"
							else
								temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							endif
							goto 10
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo
				elseif(aa_group(ii,ip)%atype1(j)=="HA2") then!!put HA2 coord next to CA in list
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="CA") then
							temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							temp_group(ii,ic)%coo1((k+1),1)=aa_group(ii,ip)%coo1(j,1)
							temp_group(ii,ic)%coo1((k+1),2)=aa_group(ii,ip)%coo1(j,2)
							temp_group(ii,ic)%coo1((k+1),3)=aa_group(ii,ip)%coo1(j,3)
							goto 10
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo	
				elseif(aa_group(ii,ip)%atype1(j)=="HA3") then!!put HA3 coord next to CA in list
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="HA2") then
							temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							temp_group(ii,ic)%coo1((k+1),1)=aa_group(ii,ip)%coo1(j,1)
							temp_group(ii,ic)%coo1((k+1),2)=aa_group(ii,ip)%coo1(j,2)
							temp_group(ii,ic)%coo1((k+1),3)=aa_group(ii,ip)%coo1(j,3)
							goto 10
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo
				endif					
10				continue
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
		else
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			flag=0
			do j=1, temp_group(ii,ic)%cnum1
				if(temp_group(ii,ic)%atype1(j)=="H2".or.temp_group(ii,ic)%atype1(j)=="H3") then
					flag=1
				endif
			enddo

			do j=1, aa_group(ii,ip)%cnum1
				if(aa_group(ii,ip)%atype1(j)=="H") then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="N") then
							if(flag==1) then
								temp_group(ii,ic)%atype1(k+1)="H1"
							else
								temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							endif
							goto 20			
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo
				endif
			enddo
20			continue
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		endif	
	elseif(temp_group(ii,ic)%gtype=="GLY".or.temp_group(ii,ic)%gtype=="NGLY".or.temp_group(ii,ic)%gtype=="CGLY") then
		if(aa_group(ii,ip)%gtype=="PRO".or.aa_group(ii,ip)%gtype=="NPRO".or.aa_group(ii,ip)%gtype=="CPRO") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			flag=0
			do j=1, (temp_group(ii,ic)%cnum1-1)
				if(temp_group(ii,ic)%atype1(j)=="H1".or.temp_group(ii,ic)%atype1(j)=="H".or.flag==1) then
					temp_group(ii,ic)%atype1(j)=temp_group(ii,ic)%atype1(j+1)
					temp_group(ii,ic)%coo1(j,1)=temp_group(ii,ic)%coo1((j+1),1)
					temp_group(ii,ic)%coo1(j,2)=temp_group(ii,ic)%coo1((j+1),2)
					temp_group(ii,ic)%coo1(j,3)=temp_group(ii,ic)%coo1((j+1),3)
					flag=1
				endif
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1

			flag=0
			do j=1, (temp_group(ii,ic)%cnum1-1)
				if(temp_group(ii,ic)%atype1(j)=="HA2") then
					temp_group(ii,ic)%atype1(j)="HA"
				elseif(temp_group(ii,ic)%atype1(j)=="HA3".or.flag==1) then
					temp_group(ii,ic)%atype1(j)=temp_group(ii,ic)%atype1(j+1)
					temp_group(ii,ic)%coo1(j,1)=temp_group(ii,ic)%coo1((j+1),1)
					temp_group(ii,ic)%coo1(j,2)=temp_group(ii,ic)%coo1((j+1),2)
					temp_group(ii,ic)%coo1(j,3)=temp_group(ii,ic)%coo1((j+1),3)
					flag=1
				endif
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		elseif(aa_group(ii,ip)%gtype=="GLY".or.aa_group(ii,ip)%gtype=="NGLY".or.aa_group(ii,ip)%gtype=="CGLY") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype			
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		else
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype			
			flag=0
			do j=1, (temp_group(ii,ic)%cnum1-1)
				if(temp_group(ii,ic)%atype1(j)=="HA2") then
					temp_group(ii,ic)%atype1(j)="HA"
				elseif(temp_group(ii,ic)%atype1(j)=="HA3".or.flag==1) then ! REMOVE HA3
					temp_group(ii,ic)%atype1(j)=temp_group(ii,ic)%atype1(j+1)
					temp_group(ii,ic)%coo1(j,1)=temp_group(ii,ic)%coo1((j+1),1)
					temp_group(ii,ic)%coo1(j,2)=temp_group(ii,ic)%coo1((j+1),2)
					temp_group(ii,ic)%coo1(j,3)=temp_group(ii,ic)%coo1((j+1),3)
					flag=1
				endif
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1											
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		endif			
	else
		if(aa_group(ii,ip)%gtype=="PRO".or.aa_group(ii,ip)%gtype=="NPRO".or.aa_group(ii,ip)%gtype=="CPRO") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			flag=0
			do j=1, (temp_group(ii,ic)%cnum1-1)
				if(temp_group(ii,ic)%atype1(j)=="H1".or.temp_group(ii,ic)%atype1(j)=="H".or.flag==1) then
					temp_group(ii,ic)%atype1(j)=temp_group(ii,ic)%atype1(j+1)
					temp_group(ii,ic)%coo1(j,1)=temp_group(ii,ic)%coo1((j+1),1)
					temp_group(ii,ic)%coo1(j,2)=temp_group(ii,ic)%coo1((j+1),2)
					temp_group(ii,ic)%coo1(j,3)=temp_group(ii,ic)%coo1((j+1),3)
					flag=1
				endif
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo			
		elseif(aa_group(ii,ip)%gtype=="GLY".or.aa_group(ii,ip)%gtype=="NGLY".or.aa_group(ii,ip)%gtype=="CGLY") then
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype
			ha2_pos=0
			do j=1, aa_group(ii,ip)%cnum1
				if(aa_group(ii,ip)%atype1(j)=="HA2") then
					do k=1, temp_group(ii,ic)%cnum1
						if(temp_group(ii,ic)%atype1(k)=="HA") then
							temp_group(ii,ic)%atype1(k)="HA2"
							temp_group(ii,ic)%coo1(k,1)=aa_group(ii,ip)%coo1(j,1)
							temp_group(ii,ic)%coo1(k,2)=aa_group(ii,ip)%coo1(j,2)
							temp_group(ii,ic)%coo1(k,3)=aa_group(ii,ip)%coo1(j,3)
							ha2_pos=k
							goto 30
						endif
					enddo
				endif					
30				continue
			enddo
			do j=1, aa_group(ii,ip)%cnum1
				if(aa_group(ii,ip)%atype1(j)=="HA3".and.ha2_pos.gt.0) then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do l=temp_group(ii,ic)%cnum1, ha2_pos+2, -1
						temp_group(ii,ic)%atype1(l)=temp_group(ii,ic)%atype1(l-1)
						temp_group(ii,ic)%coo1(l,1)=temp_group(ii,ic)%coo1(l-1,1)
						temp_group(ii,ic)%coo1(l,2)=temp_group(ii,ic)%coo1(l-1,2)
						temp_group(ii,ic)%coo1(l,3)=temp_group(ii,ic)%coo1(l-1,3)
					enddo
					temp_group(ii,ic)%atype1(ha2_pos+1)="HA3"
					temp_group(ii,ic)%coo1(ha2_pos+1,1)=aa_group(ii,ip)%coo1(j,1)
					temp_group(ii,ic)%coo1(ha2_pos+1,2)=aa_group(ii,ip)%coo1(j,2)
					temp_group(ii,ic)%coo1(ha2_pos+1,3)=aa_group(ii,ip)%coo1(j,3)
					goto 31
				endif
31				continue
			enddo
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
		else
			temp_group(ii,ic)%gtype=aa_group(ii,ip)%gtype		
			temp_group(ii,ic)%cnum2=aa_group(ii,ip)%cnum2
			do j=1, aa_group(ii,ip)%cnum2
				temp_group(ii,ic)%atype2(j)=aa_group(ii,ip)%atype2(j)
				temp_group(ii,ic)%coo2(j,1)=aa_group(ii,ip)%coo2(j,1)
				temp_group(ii,ic)%coo2(j,2)=aa_group(ii,ip)%coo2(j,2)
				temp_group(ii,ic)%coo2(j,3)=aa_group(ii,ip)%coo2(j,3)
			enddo
		endif
	endif

	return
	end subroutine residue_replace

	
	subroutine backup4sidechain(flag, chainID, ic, group, aa_backup)
	implicit none
	integer							:: flag, chainID, ic, ii, ik
	type(groupdetails)				:: group(repeated_unit,gnum), aa_backup

	ii=chainID
	if(flag==0) then
		aa_backup%gtype=group(ii,ic)%gtype

		aa_backup%cnum1=group(ii,ic)%cnum1
		do ik=1, group(ii,ic)%cnum1
			aa_backup%atype1(ik)=group(ii,ic)%atype1(ik)
			aa_backup%coo1(ik,1)=group(ii,ic)%coo1(ik,1)
			aa_backup%coo1(ik,2)=group(ii,ic)%coo1(ik,2)
			aa_backup%coo1(ik,3)=group(ii,ic)%coo1(ik,3)
		enddo

		aa_backup%cnum2=group(ii,ic)%cnum2
		do ik=1, group(ii,ic)%cnum2
			aa_backup%atype2(ik)=group(ii,ic)%atype2(ik)
			aa_backup%coo2(ik,1)=group(ii,ic)%coo2(ik,1)
			aa_backup%coo2(ik,2)=group(ii,ic)%coo2(ik,2)
			aa_backup%coo2(ik,3)=group(ii,ic)%coo2(ik,3)
		enddo

		aa_backup%cnum3=group(ii,ic)%cnum3
		do ik=1, group(ii,ic)%cnum3
			aa_backup%atype3(ik)=group(ii,ic)%atype3(ik)
			aa_backup%coo3(ik,1)=group(ii,ic)%coo3(ik,1)
			aa_backup%coo3(ik,2)=group(ii,ic)%coo3(ik,2)
			aa_backup%coo3(ik,3)=group(ii,ic)%coo3(ik,3)
		enddo

	elseif(flag==1) then
		group(ii,ic)%gtype=aa_backup%gtype

		group(ii,ic)%cnum1=aa_backup%cnum1
		do ik=1, aa_backup%cnum1
			group(ii,ic)%atype1(ik)=aa_backup%atype1(ik)
			group(ii,ic)%coo1(ik,1)=aa_backup%coo1(ik,1)
			group(ii,ic)%coo1(ik,2)=aa_backup%coo1(ik,2)
			group(ii,ic)%coo1(ik,3)=aa_backup%coo1(ik,3)
		enddo

		group(ii,ic)%cnum2=aa_backup%cnum2
		do ik=1, aa_backup%cnum2
			group(ii,ic)%atype2(ik)=aa_backup%atype2(ik)
			group(ii,ic)%coo2(ik,1)=aa_backup%coo2(ik,1)
			group(ii,ic)%coo2(ik,2)=aa_backup%coo2(ik,2)
			group(ii,ic)%coo2(ik,3)=aa_backup%coo2(ik,3)
		enddo

		group(ii,ic)%cnum3=aa_backup%cnum3
		do ik=1, aa_backup%cnum3
			group(ii,ic)%atype3(ik)=aa_backup%atype3(ik)
			group(ii,ic)%coo3(ik,1)=aa_backup%coo3(ik,1)
			group(ii,ic)%coo3(ik,2)=aa_backup%coo3(ik,2)
			group(ii,ic)%coo3(ik,3)=aa_backup%coo3(ik,3)
		enddo
	endif
	
	return
	end subroutine backup4sidechain

	
	subroutine check_transplant(chainID, ic, group, feedback)
	implicit none
	integer							:: chainID, ii, ic, ik, jj, jc, jk, feedback
	real							:: rij
	type(groupdetails)				:: group(repeated_unit,gnum)

	feedback=1
	ii=chainID
	do ik=1, group(ii,ic)%cnum2
		do jj=1, repeated_unit
			do jc=1, gnum
				if (jj.eq.ii.and.jc.eq.ic) goto 20
				do jk=1, group(jj,jc)%cnum1
					if(jj==ii.and.jc==(ic+1).and.group(ii,ic)%atype2(ik)=="CB".and.group(jj,jc)%atype1(jk)=="N") goto 50
					if(group(ii,ic)%gtype=="PRO".or.group(ii,ic)%gtype=="NPRO".or.group(ii,ic)%gtype=="CPRO") then
						if(jj==ii.and.jc==(ic-1).and.group(ii,ic)%atype2(ik)=="CD".and.group(jj,jc)%atype1(jk)=="CA") goto 50
					endif
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
						(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
					rij=sqrt(rij)
					if(rij<1.55) then
						feedback=0
						goto 10
					endif
50					continue
				enddo
				do jk=1, group(jj,jc)%cnum2
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
						(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
					rij=sqrt(rij)
					if(rij<1.55) then
						feedback=0
						goto 10
					endif
				enddo
				do jk=1, group(jj,jc)%cnum3
					if(jj==ii.and.jc==(ic-1).and.group(ii,ic)%atype2(ik)=="CB".and.group(jj,jc)%atype3(jk)=="C") goto 60
					if(group(ii,ic)%gtype=="PRO".or.group(ii,ic)%gtype=="NPRO".or.group(ii,ic)%gtype=="CPRO") then
						if(jj==ii.and.jc==(ic-1).and.group(ii,ic)%atype2(ik)=="CD") then
							goto 60
						elseif(jj==ii.and.jc==(ic-1).and.(group(ii,ic)%atype2(ik)=="CG".or.group(ii,ic)%atype2(ik)=="HD2".or.group(ii,ic)%atype2(ik)=="HD3")) then
							if(group(jj,jc)%atype3(jk)=="C") goto 60
						endif
					endif
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
						(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
					rij=sqrt(rij)
					if(rij<1.55) then
						feedback=0
						goto 10
					endif
60		    		continue
				enddo
20				continue
			enddo
		enddo
	enddo
10	continue

	return
	end subroutine check_transplant


	subroutine torsionangle4sidechain(group, chainID, ic, grade, dihedral, Dihedral4entropy)
	implicit none
	integer								:: chainID, ic, i, j
	integer								:: natom, grade, ip, jp, kp, lp
	real								:: p1(3), p2(3), p3(3), p4(3), angle
	real							    :: mdcrd(60,3)
	real								:: Dihedral4entropy(4)
	type(groupdetails)					:: group(repeated_unit,gnum)
	type(dihedralparameters)		    :: dihedral

	natom=0
	do i=1, group(chainID,ic)%cnum1
		natom=natom+1
		mdcrd(natom,1)=group(chainID,ic)%coo1(i,1)
		mdcrd(natom,2)=group(chainID,ic)%coo1(i,2)
		mdcrd(natom,3)=group(chainID,ic)%coo1(i,3)
	enddo
	do i=1, group(chainID,ic)%cnum2
		natom=natom+1
		mdcrd(natom,1)=group(chainID,ic)%coo2(i,1)
		mdcrd(natom,2)=group(chainID,ic)%coo2(i,2)
		mdcrd(natom,3)=group(chainID,ic)%coo2(i,3)
	enddo
	do i=1, group(chainID,ic)%cnum3
		natom=natom+1
		mdcrd(natom,1)=group(chainID,ic)%coo3(i,1)
		mdcrd(natom,2)=group(chainID,ic)%coo3(i,2)
		mdcrd(natom,3)=group(chainID,ic)%coo3(i,3)
	enddo

	do j=1, grade
		ip=dihedral%iph(j); jp=dihedral%jph(j); kp=dihedral%kph(j); lp=dihedral%lph(j)
		p1(1)=mdcrd(ip,1); p1(2)=mdcrd(ip,2); p1(3)=mdcrd(ip,3)
		p2(1)=mdcrd(jp,1); p2(2)=mdcrd(jp,2); p2(3)=mdcrd(jp,3)
		p3(1)=mdcrd(kp,1); p3(2)=mdcrd(kp,2); p3(3)=mdcrd(kp,3)
		p4(1)=mdcrd(lp,1); p4(2)=mdcrd(lp,2); p4(3)=mdcrd(lp,3)	
		call phipsiomg_angle(p1,p2,p3,p4,angle)
		Dihedral4entropy(j)=angle
	enddo

	return
	end subroutine torsionangle4sidechain
	
	subroutine sheet_move(ip1, group, Tgroup, ip2)
	implicit none
	integer					:: ip1, ip2, i, j, k, chainID
	real					:: ran2, displacement, ip3
	type(groupdetails)		:: group(repeated_unit,gnum), Tgroup(repeated_unit,gnum)
	
    
	call ran_gen(ran2,0)
	ip2=int(ran2*2-1.0e-3)+1
	call ran_gen(ran2,0)
    
    if (ip1 == 1) then
		ip3=ran2 * displacement_factor_x
    elseif (ip1 == 2) then
		ip3=ran2 * displacement_factor_y
    elseif (ip1 == 3) then
		ip3=ran2 * displacement_factor_z
    endif
        
	if(ip2==1) then
		displacement=ip3
	elseif(ip2==2) then
		displacement=-ip3
	endif
	
	Tgroup=group
	!if (.not.allocated(betasheets)) stop 'ERROR: beta-sheet assignments are not available for sheet move.'
	!if (num4betasheets < 1) stop 'ERROR: beta sheet 1 is not available for sheet move.'

	do i=1,betasheets(1)%num4peptides
		chainID=betasheets(1)%peptideID(i)
        do j=1, gnum
			do k=1,Tgroup(chainID,j)%cnum1
				Tgroup(chainID,j)%coo1(k,ip1)=Tgroup(chainID,j)%coo1(k,ip1)+displacement
            enddo
			
			do k=1,Tgroup(chainID,j)%cnum2
				Tgroup(chainID,j)%coo2(k,ip1)=Tgroup(chainID,j)%coo2(k,ip1)+displacement
            enddo
			
			do k=1,Tgroup(chainID,j)%cnum3
				Tgroup(chainID,j)%coo3(k,ip1)=Tgroup(chainID,j)%coo3(k,ip1)+displacement
            enddo
        enddo
    enddo
	
	return
	end subroutine sheet_move
	
	
	subroutine rmsd_calc(ip1, Tgroup, rmsd)
	implicit none
	integer					:: ip1, i, j, k, feedback
	real					:: rmsd
	real					:: N1(3), N2(3), CA1(3), CA2(3), C1(3), C2(3)
	type(groupdetails)		:: Tgroup(repeated_unit,gnum)
	
	rmsd=0.0
	do i=1, repeated_unit
		do j=1, gnum
			do k=1, Tgroup(i,j)%cnum1
				if(Tgroup(i,j)%atype1(k)=="N") then
					N1(1)=Tgroup(i,j)%coo1(k,1)
					N1(2)=Tgroup(i,j)%coo1(k,2)
					N1(3)=Tgroup(i,j)%coo1(k,3)
				elseif(Tgroup(i,j)%atype1(k)=="CA") then
					CA1(1)=Tgroup(i,j)%coo1(k,1)
					CA1(2)=Tgroup(i,j)%coo1(k,2)
					CA1(3)=Tgroup(i,j)%coo1(k,3)
				endif
			enddo
			do k=1, Tgroup(i,j)%cnum3
				if(Tgroup(i,j)%atype3(k)=="C") then
					C1(1)=Tgroup(i,j)%coo3(k,1)
					C1(2)=Tgroup(i,j)%coo3(k,2)
					C1(3)=Tgroup(i,j)%coo3(k,3)
				endif
			enddo
			
			do k=1, original_group(i,j)%cnum1
				if(original_group(i,j)%atype1(k)=="N") then
					N2(1)=original_group(i,j)%coo1(k,1)
					N2(2)=original_group(i,j)%coo1(k,2)
					N2(3)=original_group(i,j)%coo1(k,3)
				elseif(Tgroup(i,j)%atype1(k)=="CA") then
					CA2(1)=original_group(i,j)%coo1(k,1)
					CA2(2)=original_group(i,j)%coo1(k,2)
					CA2(3)=original_group(i,j)%coo1(k,3)
				endif
			enddo
			do k=1, original_group(i,j)%cnum3
				if(original_group(i,j)%atype3(k)=="C") then
					C2(1)=original_group(i,j)%coo3(k,1)
					C2(2)=original_group(i,j)%coo3(k,2)
					C2(3)=original_group(i,j)%coo3(k,3)
				endif
			enddo
			
			if(ip1==1) then
				rmsd=rmsd+(N1(1)-N2(1))**2+(CA1(1)-CA2(1))**2+(C1(1)-C2(1))**2
			elseif(ip1==2) then
				rmsd=rmsd+(N1(2)-N2(2))**2+(CA1(2)-CA2(2))**2+(C1(2)-C2(2))**2
			elseif(ip1==3) then
				rmsd=rmsd+(N1(3)-N2(3))**2+(CA1(3)-CA2(3))**2+(C1(3)-C2(3))**2
			elseif(ip1==4) then
				rmsd=rmsd+(N1(1)-N2(1))**2+(N1(2)-N2(2))**2+(N1(3)-N2(3))**2+(CA1(1)-CA2(1))**2+(CA1(2)-CA2(2))**2+(CA1(3)-CA2(3))**2+(C1(1)-C2(1))**2+(C1(2)-C2(2))**2+(C1(3)-C2(3))**2
			endif
		enddo
	enddo
	
	rmsd=sqrt(rmsd/(3*repeated_unit*gnum))
	
	return
	end subroutine rmsd_calc
	
	subroutine axis_criteria(ip1, rmsd, feedback)
	implicit none
	real					:: rmsd, rmsd_max_judge
	integer					:: ip1, feedback
	
    if (ip1 == 1) then
        rmsd_max_judge = rmsd_max_x
    elseif (ip1 == 2) then
        rmsd_max_judge = rmsd_max_y
    elseif (ip1 == 3) then
        rmsd_max_judge = rmsd_max_z
    endif
    
	if(rmsd.le.rmsd_max_judge) then
		feedback=1
	else
		feedback=0
	endif

	end subroutine axis_criteria
	
	

end module transplant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module energy_calculation

	use constant
	use mathfunction
	use database

	contains
	subroutine vdwenergy(chainID, ic, group, group_para, energy)
	implicit none
	integer							:: chainID, ii, ic, ik, jj, jc, jk
	real(kind=8)                    :: energy
    real							:: rij, epsion, r0, acoeff, bcoeff, vdw
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(energyparameters)			:: group_para(repeated_unit,gnum)

	energy=0.0
	ii=chainID
	do ik=1, group(ii,ic)%cnum2
		do jj=1, repeated_unit
			do jc=1, gnum
				if (jj==ii.and.(jc==(ic-1).or.jc==ic.or.jc==(ic+1))) goto 10

				do jk=1, group(jj,jc)%cnum1
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
					    (group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
					if(rij>100.0) goto 20
					epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion1(jk))
					r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r1(jk)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
20					continue
				enddo
				do jk=1, group(jj,jc)%cnum2
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
					    (group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
					if(rij>100.0) goto 30
					epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion2(jk))
					r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r2(jk)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
30					continue
				enddo
				do jk=1, group(jj,jc)%cnum3
					rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
					    (group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
					if(rij>100.0) goto 40
					epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion3(jk))
					r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r3(jk)
					acoeff=epsion*(r0**12)
					bcoeff=epsion*2*(r0**6)
					vdw=acoeff/(rij**6)-bcoeff/(rij**3)
					energy=energy+vdw
40					continue
				enddo
10				continue
			enddo
		enddo
	enddo

	return
	end subroutine vdwenergy


	subroutine sidechain_energy(stage, chainID, ic, group, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy)	
	implicit none
	integer								:: chainID, j, k, l, ii, ic, ik, i_id, jj, jc, jk, j_id, flag, stage
	integer								:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	integer								:: natom, dihedral_num, ip, jp, kp, lp
	real								:: rij, epsion, r0, acoeff, bcoeff, dielecons4solute, vdw, ele
    real(kind=16)                       :: energy
	real								:: p1(3), p2(3), p3(3), p4(3), angle, dihedral_energy, potential
	real							    :: mdcrd(60,3)
	type(groupdetails)					:: group(repeated_unit,gnum)
	type(energyparameters)				:: group_para(repeated_unit,gnum)
	type(dihedralparameters)		    :: dihedral

	energy=0.0
	ii=chainID
	natom=0
	do ik=1, group(ii,ic)%cnum1
		natom=natom+1
		mdcrd(natom,1)=group(ii,ic)%coo1(ik,1)
		mdcrd(natom,2)=group(ii,ic)%coo1(ik,2)
		mdcrd(natom,3)=group(ii,ic)%coo1(ik,3)
	enddo
	do ik=1, group(ii,ic)%cnum2
		natom=natom+1
		mdcrd(natom,1)=group(ii,ic)%coo2(ik,1)
		mdcrd(natom,2)=group(ii,ic)%coo2(ik,2)
		mdcrd(natom,3)=group(ii,ic)%coo2(ik,3)
	enddo
	do ik=1, group(ii,ic)%cnum3
		natom=natom+1
		mdcrd(natom,1)=group(ii,ic)%coo3(ik,1)
		mdcrd(natom,2)=group(ii,ic)%coo3(ik,2)
		mdcrd(natom,3)=group(ii,ic)%coo3(ik,3)
	enddo

	dihedral_energy=0.0
	do j=1, dihedral_num
		ip=dihedral%iph(j); jp=dihedral%jph(j); kp=dihedral%kph(j); lp=dihedral%lph(j)
		p1(1)=mdcrd(ip,1); p1(2)=mdcrd(ip,2); p1(3)=mdcrd(ip,3)
		p2(1)=mdcrd(jp,1); p2(2)=mdcrd(jp,2); p2(3)=mdcrd(jp,3)
		p3(1)=mdcrd(kp,1); p3(2)=mdcrd(kp,2); p3(3)=mdcrd(kp,3)
		p4(1)=mdcrd(lp,1); p4(2)=mdcrd(lp,2); p4(3)=mdcrd(lp,3)	
		call phipsiomg_angle(p1,p2,p3,p4,angle)
		do k=1, dihedral%multiply(j)
			potential=dihedral%pk(j,k)*(1+cosd(dihedral%pn(j,k)*angle-dihedral%phase(j,k)))
			dihedral_energy=dihedral_energy+potential
		enddo
	enddo

	do ik=1, group(ii,ic)%cnum2
		if(group(ii,ic)%atype2(ik)=="CB") goto 10
		do jj=1,repeated_unit
			do jc=1, gnum
				if(stage==0) then
					if(jj.eq.ii.and.jc.eq.ic) then
						i_id=group_para(ii,ic)%atomid2(ik)
						do jk=1, group(jj,jc)%cnum1
							j_id=group_para(jj,jc)%atomid1(jk)
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 20
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 30
								endif
							enddo
30							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion1(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r1(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(flag==1) then
								vdw=vdw/vdw14_coeff
							endif
							energy=energy+vdw
20							continue
						enddo				
						do jk=1, group(jj,jc)%cnum2
							j_id=group_para(jj,jc)%atomid2(jk)
							if(i_id.eq.j_id) goto 40					
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 40
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 50
								endif
							enddo
50							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion2(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r2(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(flag==1) then
								vdw=vdw/vdw14_coeff
							endif
							energy=energy+vdw				
40							continue
						enddo
						do jk=1, group(jj,jc)%cnum3
							j_id=group_para(jj,jc)%atomid3(jk)
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 60
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 70
								endif
							enddo
70							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion3(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r3(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(flag==1) then
								vdw=vdw/vdw14_coeff
							endif
							energy=energy+vdw
60							continue
						enddo
					else
						do jk=1, group(jj,jc)%cnum1
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion1(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r1(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							energy=energy+vdw
						enddo
						do jk=1, group(jj,jc)%cnum2
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion2(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r2(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							energy=energy+vdw					
						enddo
						do jk=1, group(jj,jc)%cnum3
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion3(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r3(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							energy=energy+vdw
						enddo
					endif

				elseif(stage==1) then
					if(jj.eq.ii.and.jc.eq.ic) then
						i_id=group_para(ii,ic)%atomid2(ik)
						do jk=1, group(jj,jc)%cnum1
							j_id=group_para(jj,jc)%atomid1(jk)
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 80
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 90
								endif
							enddo
90							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion1(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r1(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons1(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons1(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge1(jk))/(dielecons4solute*sqrt(rij))
							if(flag==1) then
								vdw=vdw/vdw14_coeff
								ele=ele/ele14_coeff
							endif
							energy=energy+vdw+ele
80							continue
						enddo				
						do jk=1, group(jj,jc)%cnum2
							j_id=group_para(jj,jc)%atomid2(jk)
							if(i_id.eq.j_id) goto 100					
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 100
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 110
								endif
							enddo
110							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion2(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r2(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons2(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons2(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge2(jk))/(dielecons4solute*sqrt(rij))
							if(flag==1) then
								vdw=vdw/vdw14_coeff
								ele=ele/ele14_coeff
							endif
							energy=energy+vdw+ele
100							continue
						enddo
						do jk=1, group(jj,jc)%cnum3
							j_id=group_para(jj,jc)%atomid3(jk)
							do l=1, S_numex(j_id)
								if(i_id.eq.S_inb(j_id,l)) goto 120
							enddo
							flag=0
							do l=1, S_numex4(j_id)
								if(i_id.eq.S_inb4(j_id,l)) then
									flag=1
									goto 130
								endif
							enddo
130							continue
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion3(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r3(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons3(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons3(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge3(jk))/(dielecons4solute*sqrt(rij))
							if(flag==1) then
								vdw=vdw/vdw14_coeff
								ele=ele/ele14_coeff
							endif
							energy=energy+vdw+ele
120							continue
						enddo
					else
						do jk=1, group(jj,jc)%cnum1
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo1(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo1(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo1(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion1(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r1(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons1(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons1(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge1(jk))/(dielecons4solute*sqrt(rij))
							energy=energy+vdw+ele
						enddo
						do jk=1, group(jj,jc)%cnum2
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo2(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo2(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo2(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion2(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r2(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons2(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons2(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge2(jk))/(dielecons4solute*sqrt(rij))
							energy=energy+vdw+ele
						enddo
						do jk=1, group(jj,jc)%cnum3
							rij=(group(ii,ic)%coo2(ik,1)-group(jj,jc)%coo3(jk,1))**2+(group(ii,ic)%coo2(ik,2)-group(jj,jc)%coo3(jk,2))**2+ &
								(group(ii,ic)%coo2(ik,3)-group(jj,jc)%coo3(jk,3))**2
							epsion=sqrt(group_para(ii,ic)%epsion2(ik)*group_para(jj,jc)%epsion3(jk))
							r0=group_para(ii,ic)%r2(ik)+group_para(jj,jc)%r3(jk)
							acoeff=epsion*(r0**12)
							bcoeff=epsion*2*(r0**6)
							vdw=acoeff/(rij**6)-bcoeff/(rij**3)
							if(group_para(ii,ic)%dielecons2(ik).ge.group_para(jj,jc)%dielecons3(jk)) then
								dielecons4solute=group_para(ii,ic)%dielecons2(ik)
							else
								dielecons4solute=group_para(jj,jc)%dielecons3(jk)
							endif
							ele=(group_para(ii,ic)%charge2(ik)*group_para(jj,jc)%charge3(jk))/(dielecons4solute*sqrt(rij))
							energy=energy+vdw+ele
						enddo			
					endif
				endif
			enddo
		enddo
10		continue
	enddo
	energy=energy+dihedral_weighting_factor*dihedral_energy
	
	return
	end  subroutine sidechain_energy


	subroutine bindingenergy_noentropy(group, group_para, numex, inb, numex4, inb4, score, binding_energy, &
			binding_vdw, binding_ele, binding_sgb, binding_snp, comp_vdw, score4hydration, Pagg)
	implicit none
	integer							:: ligan_recep_sta(num4category), ligan_recep_end(num4category)	
	integer							:: natom, ipres, atomid(repeated_unit*atom_num)
	integer, allocatable			:: sasa_chain(:),sasa_site(:),sasa_group(:)
	integer    						:: numex(repeated_unit*atom_num), inb(repeated_unit*atom_num,20), numex4(repeated_unit*atom_num), inb4(repeated_unit*atom_num,60)
	real							:: mdcrd(repeated_unit*atom_num,3), charge(repeated_unit*atom_num), epsion(repeated_unit*atom_num)
	real							:: r(repeated_unit*atom_num), rborn(repeated_unit*atom_num), fs(repeated_unit*atom_num), dielecons(repeated_unit*atom_num)
	character*4						:: lbres(repeated_unit*atom_num)
	character(len=4), allocatable	:: sasa_atom_name(:)

	integer							:: categoryID, i, ii, ic, ik, j, k, flag, num4peptides
	real							:: rij, vdw, ele, sgb
	real(kind=8)                    :: comp_vdw, comp_ele, comp_sgb, comp_snp
	real(kind=8)                    :: ligan_recep_vdw(num4category), ligan_recep_ele(num4category), ligan_recep_sgb(num4category), ligan_recep_snp(num4category)
	real(kind=8)                    :: binding_vdw, binding_ele, binding_sgb, binding_snp
    real(kind=8)					:: score, binding_energy
	real							:: score4hydration, Pagg
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(energyparameters)			:: group_para(repeated_unit,gnum)
	real, dimension(:), allocatable :: alpha

	call hydrationscale(group, score4hydration)
    call aggregation_propensity(group, Pagg)
	if (ESURF_flag == 2) allocate(sasa_chain(repeated_unit*atom_num),sasa_site(repeated_unit*atom_num), &
		sasa_group(repeated_unit*atom_num),sasa_atom_name(repeated_unit*atom_num))

	natom=0
	flag=0
	num4peptides=0
	do categoryID=1, num4category
		num4peptides=num4peptides+selfassembly(categoryID)%num4peptides
		do i=1, selfassembly(categoryID)%num4peptides
			ii=selfassembly(categoryID)%peptideID(i)
			do ic=1, gnum			
				if(i==1.and.ic==1.and.flag==0) then
					ligan_recep_sta(categoryID)=natom+1
					flag=1
				endif

				ipres=natom
				do ik=1, group(ii,ic)%cnum1
					natom=natom+1
					if (ESURF_flag == 2) then
						sasa_chain(natom)=ii
						sasa_site(natom)=ic
						sasa_group(natom)=categoryID
						sasa_atom_name(natom)=group(ii,ic)%atype1(ik)
					endif
					charge(natom)=group_para(ii,ic)%charge1(ik)
					epsion(natom)=group_para(ii,ic)%epsion1(ik)
					r(natom)=group_para(ii,ic)%r1(ik)
					rborn(natom)=group_para(ii,ic)%rborn1(ik)
					fs(natom)=group_para(ii,ic)%fs1(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons1(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid1(ik)
					mdcrd(natom,1)=group(ii,ic)%coo1(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo1(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo1(ik,3)
					lbres(natom)=group(ii,ic)%gtype	
				enddo
				
				do ik=1, group(ii,ic)%cnum2
					natom=natom+1
					if (ESURF_flag == 2) then
						sasa_chain(natom)=ii
						sasa_site(natom)=ic
						sasa_group(natom)=categoryID
						sasa_atom_name(natom)=group(ii,ic)%atype2(ik)
					endif
					charge(natom)=group_para(ii,ic)%charge2(ik)
					epsion(natom)=group_para(ii,ic)%epsion2(ik)
					r(natom)=group_para(ii,ic)%r2(ik)
					rborn(natom)=group_para(ii,ic)%rborn2(ik)
					fs(natom)=group_para(ii,ic)%fs2(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons2(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid2(ik)
					mdcrd(natom,1)=group(ii,ic)%coo2(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo2(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo2(ik,3)
					lbres(natom)=group(ii,ic)%gtype
				enddo
	
				do ik=1, group(ii,ic)%cnum3
					natom=natom+1
					if (ESURF_flag == 2) then
						sasa_chain(natom)=ii
						sasa_site(natom)=ic
						sasa_group(natom)=categoryID
						sasa_atom_name(natom)=group(ii,ic)%atype3(ik)
					endif
					charge(natom)=group_para(ii,ic)%charge3(ik)
					epsion(natom)=group_para(ii,ic)%epsion3(ik)
					r(natom)=group_para(ii,ic)%r3(ik)
					rborn(natom)=group_para(ii,ic)%rborn3(ik)
					fs(natom)=group_para(ii,ic)%fs3(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons3(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid3(ik)
					mdcrd(natom,1)=group(ii,ic)%coo3(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo3(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo3(ik,3)
					lbres(natom)=group(ii,ic)%gtype
				enddo
				
				if(i==selfassembly(categoryID)%num4peptides.and.ic==gnum.and.flag==1) then
					ligan_recep_end(categoryID)=natom
					flag=0
				endif
				
			enddo
		enddo
	enddo

	allocate (alpha(natom))
	do i=1, natom
		call effgbradi(1,natom,i,rborn,fs,mdcrd,alpha(i))
	enddo

	comp_vdw=0.0; comp_ele=0.0;	comp_sgb=0.0; comp_snp=0.0
	do i=1, natom
		do j=i, natom
			vdw=0.0; ele=0.0; sgb=0.0
			rij=(mdcrd(i,1)-mdcrd(j,1))**2+(mdcrd(i,2)-mdcrd(j,2))**2+(mdcrd(i,3)-mdcrd(j,3))**2
			if(j.ne.i) then
				do k=1,numex(atomid(i))
					if(atomid(j).eq.inb(atomid(i),k)) goto 40
				enddo
				flag=0
				do k=1,numex4(atomid(i))
					if(atomid(j).eq.inb4(atomid(i),k)) then
						flag=1
						goto 45
					endif
				enddo
45				continue
				call vdwcontri(i,j,rij,epsion,r,vdw)
				call elecontri(i,j,rij,charge,dielecons,ele)
				if(flag==1) then
					vdw=vdw/vdw14_coeff
					ele=ele/ele14_coeff
				endif
			endif
40			continue
			call sgbcontri(i,j,rij,alpha(i),alpha(j),charge,dielecons,sgb)
			if(j.eq.i) sgb=sgb/2.0
			comp_vdw=comp_vdw+vdw
			comp_ele=comp_ele+ele
			comp_sgb=comp_sgb+sgb
		enddo
	enddo
	deallocate (alpha)
	if (ESURF_flag == 1) then
		call snp_energy(1, natom, rborn, mdcrd, comp_snp)
	endif

	do categoryID=1, num4category	
		allocate (alpha(natom))
		do i=ligan_recep_sta(categoryID), ligan_recep_end(categoryID)
			call effgbradi(ligan_recep_sta(categoryID),ligan_recep_end(categoryID),i,rborn,fs,mdcrd,alpha(i))
		enddo

		ligan_recep_vdw(categoryID)=0.0; ligan_recep_ele(categoryID)=0.0; ligan_recep_sgb(categoryID)=0.0; ligan_recep_snp(categoryID)=0.0
		do i=ligan_recep_sta(categoryID), ligan_recep_end(categoryID)
			do j=i, ligan_recep_end(categoryID)
				vdw=0.0; ele=0.0; sgb=0.0
				rij= (mdcrd(i,1)-mdcrd(j,1))**2+(mdcrd(i,2)-mdcrd(j,2))**2+(mdcrd(i,3)-mdcrd(j,3))**2
				if(j.ne.i) then
					do k=1,numex(atomid(i))
						if(atomid(j).eq.inb(atomid(i),k)) goto 50
					enddo
					flag=0
					do k=1,numex4(atomid(i))
						if(atomid(j).eq.inb4(atomid(i),k)) then
							flag=1
							goto 55
						endif
					enddo
55					continue
					call vdwcontri(i,j,rij,epsion,r,vdw)
					call elecontri(i,j,rij,charge,dielecons,ele)
					if(flag==1) then
						vdw=vdw/vdw14_coeff
						ele=ele/ele14_coeff
					endif						
				endif
50				continue
				call sgbcontri(i,j,rij,alpha(i),alpha(j),charge,dielecons,sgb)
				if(j.eq.i) sgb=sgb/2.0		
				ligan_recep_vdw(categoryID)=ligan_recep_vdw(categoryID)+vdw
				ligan_recep_ele(categoryID)=ligan_recep_ele(categoryID)+ele	
 				ligan_recep_sgb(categoryID)=ligan_recep_sgb(categoryID)+sgb
			enddo
		enddo
		if (ESURF_flag == 1) then
			call snp_energy(ligan_recep_sta(categoryID), ligan_recep_end(categoryID), rborn, &
							mdcrd, ligan_recep_snp(categoryID))
		endif
		deallocate (alpha)
		
	enddo
	
	binding_vdw=0; binding_ele=0; binding_sgb=0; binding_snp=0
	do categoryID=1, num4category
		binding_vdw=binding_vdw-ligan_recep_vdw(categoryID)
		binding_ele=binding_ele-ligan_recep_ele(categoryID)
		binding_sgb=binding_sgb-ligan_recep_sgb(categoryID)
		if (ESURF_flag /= 2) binding_snp=binding_snp-ligan_recep_snp(categoryID)
    enddo

    binding_vdw=(binding_vdw+comp_vdw)/real(num4peptides)
    binding_ele=(binding_ele+comp_ele)/real(num4peptides)
    binding_sgb=(binding_sgb+comp_sgb)/real(num4peptides)
	if (ESURF_flag == 2) then
		if (sasa_current%valid) then
			call update_sasa_cache(natom,rborn,mdcrd,sasa_chain,sasa_site,sasa_group, &
				sasa_atom_name,sasa_current,sasa_trial,binding_snp)
		else
			call build_sasa_cache(natom,rborn,mdcrd,sasa_chain,sasa_site,sasa_group, &
				sasa_atom_name,sasa_trial,binding_snp)
		endif
	else
		binding_snp=(binding_snp+comp_snp)/real(num4peptides)
	endif
	if (allocated(sasa_chain)) deallocate(sasa_chain,sasa_site,sasa_group,sasa_atom_name)
    
	binding_energy=binding_vdw + binding_ele + binding_sgb + binding_snp

	score=binding_energy-propensity_weighting_factor*(score4hydration+Pagg)
	return
	end subroutine bindingenergy_noentropy

	
	subroutine entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)
	implicit none
	integer							:: rotanum, obs, grade
	real							:: matrix(34,4)
	real							:: det, entropy
	character*4						:: aminoacid_name

	if(aminoacid_name=="GLY".or.aminoacid_name=="NGLY".or.aminoacid_name=="CGLY".or.aminoacid_name=="PRO".or.aminoacid_name=="NPRO".or.aminoacid_name=="CPRO".or.  &
	   aminoacid_name=="CYS".or.aminoacid_name=="NCYS".or.aminoacid_name=="CCYS".or.aminoacid_name=="ALA".or.aminoacid_name=="NALA".or.aminoacid_name=="CALA".or.  &
	   aminoacid_name=="VAL".or.aminoacid_name=="NVAL".or.aminoacid_name=="CVAL".or.aminoacid_name=="SER".or.aminoacid_name=="NSER".or.aminoacid_name=="CSER".or.  &
	   aminoacid_name=="THR".or.aminoacid_name=="NTHR".or.aminoacid_name=="CTHR".or.aminoacid_name=="ACE".or.aminoacid_name=="NME".or.aminoacid_name=="NHE") then
		entropy=0.0
	elseif(aminoacid_name=="PHE".or.aminoacid_name=="NPHE".or.aminoacid_name=="CPHE".or.aminoacid_name=="TYR".or.aminoacid_name=="NTYR".or.aminoacid_name=="CTYR") then
		if(obs.le.2) then
			entropy=-2.260
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.260
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="ARG".or.aminoacid_name=="NARG".or.aminoacid_name=="CARG") then
		if(obs.le.2) then
			entropy=-4.863
		else	
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-4.863
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="LYS".or.aminoacid_name=="NLYS".or.aminoacid_name=="CLYS") then
		if(obs.le.2) then
			entropy=-4.621
		else	
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-4.621
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="ASN".or.aminoacid_name=="NASN".or.aminoacid_name=="CASN") then
		if(obs.le.2) then
			entropy=-2.056
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.056
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif	
	elseif(aminoacid_name=="ASP".or.aminoacid_name=="NASP".or.aminoacid_name=="CASP") then
		if(obs.le.2) then
			entropy=-1.790
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-1.790
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="GLN".or.aminoacid_name=="NGLN".or.aminoacid_name=="CGLN") then
		if(obs.le.2) then
			entropy=-3.205
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-3.205
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif	
	elseif(aminoacid_name=="GLU".or.aminoacid_name=="NGLU".or.aminoacid_name=="CGLU") then
		if(obs.le.2) then
			entropy=-2.541
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.541
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="HIE".or.aminoacid_name=="NHIE".or.aminoacid_name=="CHIE") then
		if(obs.le.2) then
			entropy=-2.411
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.411
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="ILE".or.aminoacid_name=="NILE".or.aminoacid_name=="CILE") then
		if(obs.le.2) then
			entropy=-2.239
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.239
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif	
	elseif(aminoacid_name=="LEU".or.aminoacid_name=="NLEU".or.aminoacid_name=="CLEU") then
		if(obs.le.2) then
			entropy=-1.946
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-1.946
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="MET".or.aminoacid_name=="NMET".or.aminoacid_name=="CMET") then
		if(obs.le.2) then
			entropy=-3.495
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-3.495
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	elseif(aminoacid_name=="TRP".or.aminoacid_name=="NTRP".or.aminoacid_name=="CTRP") then
		if(obs.le.2) then
			entropy=-2.316
		else
			call variance_covariance_matrix(matrix,obs,grade,det)

			if(det.le.1.0e-14) then
				entropy=-2.316
			else
				entropy=0.3*(real(grade)*(1+log(6.2831852))+log(det))*log(real(obs)/real(rotanum))
			endif
		endif
	endif
	
	return
	end subroutine entropy_calculation

	
	subroutine conformation_entropy(entropy4individual, entropy)
	implicit none
	integer							:: categoryID, chainID, ic, i, j
	real							:: entropy4individual(repeated_unit,gnum), entropy4sort(repeated_unit,gnum)
	real							:: minimum, entropy

	do ic=1, gnum
		chainID=1
		minimum=entropy4individual(chainID,ic)
		do chainID=2, repeated_unit
			if(entropy4individual(chainID,ic).lt.minimum) minimum=entropy4individual(chainID,ic)
		enddo

		do chainID=1, repeated_unit
			entropy4sort(chainID,ic)=minimum
		enddo
	enddo
	
	entropy=0.0
	chainID=1
!	do chainID=1, repeated_unit
		do ic=1,gnum	
			entropy=entropy+entropy4sort(chainID,ic)
		enddo
!	enddo

	return
	end subroutine conformation_entropy	

	
	subroutine bindingenergy(group, group_para, entropy4individual, numex, inb, numex4, inb4, score, &
		binding_energy, binding_vdw, binding_ele, binding_sgb, binding_snp, entropy, score4hydration, Pagg)
	implicit none
	integer							:: ligan_recep_sta(num4category), ligan_recep_end(num4category)	
	integer							:: natom, ipres, atomid(repeated_unit*atom_num)
	integer, allocatable			:: sasa_chain(:),sasa_site(:),sasa_group(:)
	integer    						:: numex(repeated_unit*atom_num), inb(repeated_unit*atom_num,20), numex4(repeated_unit*atom_num), inb4(repeated_unit*atom_num,60)
	real							:: mdcrd(repeated_unit*atom_num,3), charge(repeated_unit*atom_num), epsion(repeated_unit*atom_num)
	real							:: r(repeated_unit*atom_num), rborn(repeated_unit*atom_num), fs(repeated_unit*atom_num), dielecons(repeated_unit*atom_num)
	real							:: entropy4individual(repeated_unit,gnum)
	character*4						:: lbres(repeated_unit*atom_num)
	character(len=4), allocatable	:: sasa_atom_name(:)

	integer							:: categoryID, i, ii, ic, ik, j, k, flag, num4peptides
	real							:: rij, vdw, ele, sgb
	real(kind=8)                    :: comp_vdw, comp_ele, comp_sgb, comp_snp
	real(kind=8)                    :: ligan_recep_vdw(num4category), ligan_recep_ele(num4category), ligan_recep_sgb(num4category), ligan_recep_snp(num4category)
	real(kind=8)                    :: binding_vdw, binding_ele, binding_sgb, binding_snp
	real(kind=8)                    :: score, binding_energy
    real                            :: entropy, score4hydration, Pagg
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(energyparameters)			:: group_para(repeated_unit,gnum)
	real, dimension(:), allocatable :: alpha

	call conformation_entropy(entropy4individual, entropy)	
	call hydrationscale(group, score4hydration)
    call aggregation_propensity(group, Pagg)
	if (ESURF_flag == 2) allocate(sasa_chain(repeated_unit*atom_num),sasa_site(repeated_unit*atom_num), &
		sasa_group(repeated_unit*atom_num),sasa_atom_name(repeated_unit*atom_num))

	natom=0
	flag=0
	num4peptides=0
	do categoryID=1, num4category
		num4peptides=num4peptides+selfassembly(categoryID)%num4peptides
		do i=1, selfassembly(categoryID)%num4peptides
			ii=selfassembly(categoryID)%peptideID(i)
			do ic=1, gnum			
				if(i==1.and.ic==1.and.flag==0) then
					ligan_recep_sta(categoryID)=natom+1
					flag=1
				endif

				ipres=natom
				do ik=1, group(ii,ic)%cnum1
					natom=natom+1
					if (ESURF_flag == 2) then
						sasa_chain(natom)=ii
						sasa_site(natom)=ic
						sasa_group(natom)=categoryID
						sasa_atom_name(natom)=group(ii,ic)%atype1(ik)
					endif
					charge(natom)=group_para(ii,ic)%charge1(ik)
					epsion(natom)=group_para(ii,ic)%epsion1(ik)
					r(natom)=group_para(ii,ic)%r1(ik)
					rborn(natom)=group_para(ii,ic)%rborn1(ik)
					fs(natom)=group_para(ii,ic)%fs1(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons1(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid1(ik)
					mdcrd(natom,1)=group(ii,ic)%coo1(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo1(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo1(ik,3)
					lbres(natom)=group(ii,ic)%gtype	
				enddo

				do ik=1, group(ii,ic)%cnum2
					natom=natom+1
					if (ESURF_flag == 2) then
						sasa_chain(natom)=ii
						sasa_site(natom)=ic
						sasa_group(natom)=categoryID
						sasa_atom_name(natom)=group(ii,ic)%atype2(ik)
					endif
					charge(natom)=group_para(ii,ic)%charge2(ik)
					epsion(natom)=group_para(ii,ic)%epsion2(ik)
					r(natom)=group_para(ii,ic)%r2(ik)
					rborn(natom)=group_para(ii,ic)%rborn2(ik)
					fs(natom)=group_para(ii,ic)%fs2(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons2(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid2(ik)
					mdcrd(natom,1)=group(ii,ic)%coo2(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo2(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo2(ik,3)
					lbres(natom)=group(ii,ic)%gtype
				enddo
	
				do ik=1, group(ii,ic)%cnum3
					natom=natom+1
					if (ESURF_flag == 2) then
						sasa_chain(natom)=ii
						sasa_site(natom)=ic
						sasa_group(natom)=categoryID
						sasa_atom_name(natom)=group(ii,ic)%atype3(ik)
					endif
					charge(natom)=group_para(ii,ic)%charge3(ik)
					epsion(natom)=group_para(ii,ic)%epsion3(ik)
					r(natom)=group_para(ii,ic)%r3(ik)
					rborn(natom)=group_para(ii,ic)%rborn3(ik)
					fs(natom)=group_para(ii,ic)%fs3(ik)
					dielecons(natom)=group_para(ii,ic)%dielecons3(ik)
					atomid(natom)=ipres+group_para(ii,ic)%atomid3(ik)
					mdcrd(natom,1)=group(ii,ic)%coo3(ik,1)
					mdcrd(natom,2)=group(ii,ic)%coo3(ik,2)
					mdcrd(natom,3)=group(ii,ic)%coo3(ik,3)
					lbres(natom)=group(ii,ic)%gtype
				enddo

				if(i==selfassembly(categoryID)%num4peptides.and.ic==gnum.and.flag==1) then
					ligan_recep_end(categoryID)=natom
					flag=0
				endif

			enddo
		enddo
    enddo
    
	allocate (alpha(natom))
	do i=1, natom
		call effgbradi(1,natom,i,rborn,fs,mdcrd,alpha(i))
	enddo

	comp_vdw=0.0; comp_ele=0.0;	comp_sgb=0.0; comp_snp=0.0
	do i=1, natom
		do j=i, natom
			vdw=0.0; ele=0.0; sgb=0.0
			rij=(mdcrd(i,1)-mdcrd(j,1))**2+(mdcrd(i,2)-mdcrd(j,2))**2+(mdcrd(i,3)-mdcrd(j,3))**2
			if(j.ne.i) then
				do k=1,numex(atomid(i))
					if(atomid(j).eq.inb(atomid(i),k)) goto 40
				enddo
				flag=0
				do k=1,numex4(atomid(i))
					if(atomid(j).eq.inb4(atomid(i),k)) then
						flag=1
						goto 45
					endif
				enddo
45				continue
				call vdwcontri(i,j,rij,epsion,r,vdw)
				call elecontri(i,j,rij,charge,dielecons,ele)
				if(flag==1) then
					vdw=vdw/vdw14_coeff
					ele=ele/ele14_coeff
				endif
			endif
40			continue
			call sgbcontri(i,j,rij,alpha(i),alpha(j),charge,dielecons,sgb)
			if(j.eq.i) sgb=sgb/2.0
			comp_vdw=comp_vdw+vdw
			comp_ele=comp_ele+ele
			comp_sgb=comp_sgb+sgb
		enddo
	enddo
	deallocate (alpha)
	if (ESURF_flag == 1) then
		call snp_energy(1, natom, rborn, mdcrd, comp_snp)
	endif

	do categoryID=1, num4category	
		allocate (alpha(natom))
		do i=ligan_recep_sta(categoryID), ligan_recep_end(categoryID)
			call effgbradi(ligan_recep_sta(categoryID),ligan_recep_end(categoryID),i,rborn,fs,mdcrd,alpha(i))
		enddo

		ligan_recep_vdw(categoryID)=0.0; ligan_recep_ele(categoryID)=0.0; ligan_recep_sgb(categoryID)=0.0; ligan_recep_snp(categoryID)=0.0
		do i=ligan_recep_sta(categoryID), ligan_recep_end(categoryID)
			do j=i, ligan_recep_end(categoryID)
				vdw=0.0; ele=0.0; sgb=0.0
				rij= (mdcrd(i,1)-mdcrd(j,1))**2+(mdcrd(i,2)-mdcrd(j,2))**2+(mdcrd(i,3)-mdcrd(j,3))**2
				if(j.ne.i) then
					do k=1,numex(atomid(i))
						if(atomid(j).eq.inb(atomid(i),k)) goto 50
					enddo
					flag=0
					do k=1,numex4(atomid(i))
						if(atomid(j).eq.inb4(atomid(i),k)) then
							flag=1
							goto 55
						endif
					enddo
55					continue
					call vdwcontri(i,j,rij,epsion,r,vdw)
					call elecontri(i,j,rij,charge,dielecons,ele)
					if(flag==1) then
						vdw=vdw/vdw14_coeff
						ele=ele/ele14_coeff
					endif						
				endif
50				continue
				call sgbcontri(i,j,rij,alpha(i),alpha(j),charge,dielecons,sgb)
				if(j.eq.i) sgb=sgb/2.0		
				ligan_recep_vdw(categoryID)=ligan_recep_vdw(categoryID)+vdw
				ligan_recep_ele(categoryID)=ligan_recep_ele(categoryID)+ele	
 				ligan_recep_sgb(categoryID)=ligan_recep_sgb(categoryID)+sgb
			enddo
		enddo
		if (ESURF_flag == 1) then
			call snp_energy(ligan_recep_sta(categoryID), ligan_recep_end(categoryID), rborn, &
							mdcrd, ligan_recep_snp(categoryID))
		endif
		deallocate (alpha)
		
	enddo
	
	binding_vdw=0; binding_ele=0; binding_sgb=0; binding_snp=0
	do categoryID=1, num4category
		binding_vdw=binding_vdw-ligan_recep_vdw(categoryID)
		binding_ele=binding_ele-ligan_recep_ele(categoryID)
		binding_sgb=binding_sgb-ligan_recep_sgb(categoryID)
		if (ESURF_flag /= 2) binding_snp=binding_snp-ligan_recep_snp(categoryID)
    enddo

    binding_vdw=(binding_vdw+comp_vdw)/real(num4peptides)
    binding_ele=(binding_ele+comp_ele)/real(num4peptides)
    binding_sgb=(binding_sgb+comp_sgb)/real(num4peptides)
	if (ESURF_flag == 2) then
		if (sasa_current%valid) then
			call update_sasa_cache(natom,rborn,mdcrd,sasa_chain,sasa_site,sasa_group, &
				sasa_atom_name,sasa_current,sasa_trial,binding_snp)
		else
			call build_sasa_cache(natom,rborn,mdcrd,sasa_chain,sasa_site,sasa_group, &
				sasa_atom_name,sasa_trial,binding_snp)
		endif
		if (.not.sasa_current%valid) sasa_current=sasa_trial
	else
		binding_snp=(binding_snp+comp_snp)/real(num4peptides)
	endif
	if (allocated(sasa_chain)) deallocate(sasa_chain,sasa_site,sasa_group,sasa_atom_name)
    
	binding_energy=binding_vdw + binding_ele + binding_sgb + binding_snp

	score=binding_energy-entropy-propensity_weighting_factor*(score4hydration+Pagg)

	return
	end subroutine bindingenergy


	subroutine hydrationscale(group, score4hydration)
	implicit none
	integer							:: i, j, ic, chainID
	integer							:: group_num(23)
	real							:: gly, leu, val, ile, met, phe, tyr, trp, glu, asp
	real							:: arg, lys, ser, thr, asn, gln, hie, ala, cys, pro
	real							:: avg_pho, avg_pol, avg_oth
	real							:: score_pho, score_pol, score_oth
	real							:: Npho, Npol, Noth
	real							:: Tscore, score4hydration
	type(groupdetails)				:: group(repeated_unit,gnum)

	ala=0.42; leu=2.32; val=1.66; ile=2.46; met=1.68; phe=2.44; tyr=1.31; trp=3.07; gly=0; glu=-0.87
	asp=-1.05; arg=-1.37; lys=-1.35; ser=-0.05; thr=0.35; asn=-0.82; gln=-0.30; hie=0.18; cys=1.34; pro=0.98
	avg_pho=1.58; avg_pol=-0.63; avg_oth=1.10	
	score4hydration=0.0

	chainID=1
	group_num=0.0
	do ic=1,gnum
		if(group(chainID,ic)%gtype=="ALA".or.group(chainID,ic)%gtype=="NALA".or.group(chainID,ic)%gtype=="CALA") then
			group_num(1)=group_num(1)+1.0
		elseif(group(chainID,ic)%gtype=="LEU".or.group(chainID,ic)%gtype=="NLEU".or.group(chainID,ic)%gtype=="CLEU") then
			group_num(2)=group_num(2)+1.0
		elseif(group(chainID,ic)%gtype=="VAL".or.group(chainID,ic)%gtype=="NVAL".or.group(chainID,ic)%gtype=="CVAL") then
			group_num(3)=group_num(3)+1.0
		elseif(group(chainID,ic)%gtype=="ILE".or.group(chainID,ic)%gtype=="NILE".or.group(chainID,ic)%gtype=="CILE") then
			group_num(4)=group_num(4)+1.0
		elseif(group(chainID,ic)%gtype=="MET".or.group(chainID,ic)%gtype=="NMET".or.group(chainID,ic)%gtype=="CMET") then
			group_num(5)=group_num(5)+1.0
		elseif(group(chainID,ic)%gtype=="PHE".or.group(chainID,ic)%gtype=="NPHE".or.group(chainID,ic)%gtype=="CPHE") then
			group_num(6)=group_num(6)+1.0
		elseif(group(chainID,ic)%gtype=="TYR".or.group(chainID,ic)%gtype=="NTYR".or.group(chainID,ic)%gtype=="CTYR") then
			group_num(7)=group_num(7)+1.0
		elseif(group(chainID,ic)%gtype=="TRP".or.group(chainID,ic)%gtype=="NTRP".or.group(chainID,ic)%gtype=="CTRP") then
			group_num(8)=group_num(8)+1.0
		elseif(group(chainID,ic)%gtype=="GLY".or.group(chainID,ic)%gtype=="NGLY".or.group(chainID,ic)%gtype=="CGLY") then
			group_num(9)=group_num(9)+1.0
		elseif(group(chainID,ic)%gtype=="GLU".or.group(chainID,ic)%gtype=="NGLU".or.group(chainID,ic)%gtype=="CGLU") then
			group_num(10)=group_num(10)+1.0
		elseif(group(chainID,ic)%gtype=="ASP".or.group(chainID,ic)%gtype=="NASP".or.group(chainID,ic)%gtype=="CASP") then
			group_num(11)=group_num(11)+1.0
		elseif(group(chainID,ic)%gtype=="ARG".or.group(chainID,ic)%gtype=="NARG".or.group(chainID,ic)%gtype=="CARG") then
			group_num(12)=group_num(12)+1.0
		elseif(group(chainID,ic)%gtype=="LYS".or.group(chainID,ic)%gtype=="NLYS".or.group(chainID,ic)%gtype=="CLYS") then
			group_num(13)=group_num(13)+1.0
		elseif(group(chainID,ic)%gtype=="ASN".or.group(chainID,ic)%gtype=="NASN".or.group(chainID,ic)%gtype=="CASN") then
			group_num(14)=group_num(14)+1.0
		elseif(group(chainID,ic)%gtype=="GLN".or.group(chainID,ic)%gtype=="NGLN".or.group(chainID,ic)%gtype=="CGLN") then
			group_num(15)=group_num(15)+1.0
		elseif(group(chainID,ic)%gtype=="SER".or.group(chainID,ic)%gtype=="NSER".or.group(chainID,ic)%gtype=="CSER") then
			group_num(16)=group_num(16)+1.0
		elseif(group(chainID,ic)%gtype=="THR".or.group(chainID,ic)%gtype=="NTHR".or.group(chainID,ic)%gtype=="CTHR") then
			group_num(17)=group_num(17)+1.0
		elseif(group(chainID,ic)%gtype=="HIE".or.group(chainID,ic)%gtype=="NHIE".or.group(chainID,ic)%gtype=="CHIE") then
			group_num(18)=group_num(18)+1.0
		elseif(group(chainID,ic)%gtype=="CYS".or.group(chainID,ic)%gtype=="NCYS".or.group(chainID,ic)%gtype=="CCYS") then
			group_num(19)=group_num(19)+1.0
		elseif(group(chainID,ic)%gtype=="PRO".or.group(chainID,ic)%gtype=="NPRO".or.group(chainID,ic)%gtype=="CPRO") then
			group_num(20)=group_num(20)+1.0
		elseif(group(chainID,ic)%gtype=="ACE") then
			group_num(21)=group_num(21)+1.0
		elseif(group(chainID,ic)%gtype=="NME") then
			group_num(22)=group_num(22)+1.0
		elseif(group(chainID,ic)%gtype=="NHE") then
			group_num(23)=group_num(23)+1.0            
            
		endif
	enddo

	score_pho=real(group_num(1)*ala+group_num(2)*leu+group_num(3)*val+group_num(4)*ile+group_num(5)*met+group_num(6)*phe+group_num(7)*tyr+group_num(8)*trp+group_num(9)*gly)
	score_pol=real(group_num(10)*glu+group_num(11)*asp+group_num(12)*arg+group_num(13)*lys+group_num(14)*asn+group_num(15)*gln+group_num(16)*ser+group_num(17)*thr+group_num(18)*hie)
	score_oth=real(group_num(19)*cys+group_num(20)*pro)
	
	Npho=group_num(1)+group_num(2)+group_num(3)+group_num(4)+group_num(5)+group_num(6)+group_num(7)+group_num(8)+group_num(9)
	Npol=group_num(10)+group_num(11)+group_num(12)+group_num(13)+group_num(14)+group_num(15)+group_num(16)+group_num(17)+group_num(18)
	Noth=group_num(19)+group_num(20)

	Tscore=abs(score_pho-avg_pho*Npho)+abs(score_pol-avg_pol*Npol)+abs(score_oth-avg_oth*Noth)
!	score4hydration=score4hydration+Tscore
	score4hydration=0-Tscore

	return
	end subroutine hydrationscale

	
	subroutine aggregation_propensity(group,Pagg)
	implicit none	
	integer							:: chainID, ic
	real							:: charge, helix, sheet, pat, pat_old
	real							:: CH, Alpha, Belta, Pattern
	real, parameter					:: coeff_CH=-0.16, coeff_Alpha=-5.7, coeff_Belta=5.0, coeff_Pat=0.39
	real							:: Pagg, p
	type(groupdetails)				:: group(repeated_unit,gnum)

	Pagg=0.0
	chainID=1
	CH=0.0; Alpha=0.0; Belta=0.0; Pattern=0.0; pat_old=0.5
	do ic=1,gnum
		if(group(chainID,ic)%gtype=="ALA".or.group(chainID,ic)%gtype=="NALA".or.group(chainID,ic)%gtype=="CALA") then
			charge=0.0;   helix=0.04;    sheet=0.12;    pat=1.0
		elseif(group(chainID,ic)%gtype=="LEU".or.group(chainID,ic)%gtype=="NLEU".or.group(chainID,ic)%gtype=="CLEU") then
			charge=0.0;   helix=0.38;    sheet=-0.15;   pat=1.0
		elseif(group(chainID,ic)%gtype=="VAL".or.group(chainID,ic)%gtype=="NVAL".or.group(chainID,ic)%gtype=="CVAL") then
			charge=0.0;   helix=0.06;    sheet=0.70;    pat=1.0
		elseif(group(chainID,ic)%gtype=="ILE".or.group(chainID,ic)%gtype=="NILE".or.group(chainID,ic)%gtype=="CILE") then
			charge=0.0;   helix=0.26;    sheet=0.77;    pat=1.0
		elseif(group(chainID,ic)%gtype=="MET".or.group(chainID,ic)%gtype=="NMET".or.group(chainID,ic)%gtype=="CMET") then
			charge=0.0;   helix=0.09;    sheet=0.71;    pat=1.0
		elseif(group(chainID,ic)%gtype=="PHE".or.group(chainID,ic)%gtype=="NPHE".or.group(chainID,ic)%gtype=="CPHE") then
			charge=0.0;   helix=0.01;    sheet=0.67;    pat=1.0
		elseif(group(chainID,ic)%gtype=="TYR".or.group(chainID,ic)%gtype=="NTYR".or.group(chainID,ic)%gtype=="CTYR") then
			charge=0.0;   helix=-0.05;   sheet=0.49;    pat=1.0
		elseif(group(chainID,ic)%gtype=="TRP".or.group(chainID,ic)%gtype=="NTRP".or.group(chainID,ic)%gtype=="CTRP") then
			charge=0.0;   helix=-0.21;   sheet=0.14;    pat=1.0
		elseif(group(chainID,ic)%gtype=="GLY".or.group(chainID,ic)%gtype=="NGLY".or.group(chainID,ic)%gtype=="CGLY") then
			charge=0.0;   helix=-1.24;   sheet=-0.76;   pat=0.0
		elseif(group(chainID,ic)%gtype=="GLU".or.group(chainID,ic)%gtype=="NGLU".or.group(chainID,ic)%gtype=="CGLU") then
			charge=-1.0;  helix=0.33;    sheet=-0.91;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="ASP".or.group(chainID,ic)%gtype=="NASP".or.group(chainID,ic)%gtype=="CASP") then
			charge=-1.0;  helix=-0.27;   sheet=-1.12;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="ARG".or.group(chainID,ic)%gtype=="NARG".or.group(chainID,ic)%gtype=="CARG") then
			charge=1.0;   helix=1.30;    sheet=-1.34;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="LYS".or.group(chainID,ic)%gtype=="NLYS".or.group(chainID,ic)%gtype=="CLYS") then
			charge=1.0;   helix=0.18;    sheet=-0.29;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="ASN".or.group(chainID,ic)%gtype=="NASN".or.group(chainID,ic)%gtype=="CASN") then
			charge=0.0;   helix=-0.25;   sheet=-1.05;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="GLN".or.group(chainID,ic)%gtype=="NGLN".or.group(chainID,ic)%gtype=="CGLN") then
			charge=0.0;   helix=0.02;    sheet=-1.67;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="SER".or.group(chainID,ic)%gtype=="NSER".or.group(chainID,ic)%gtype=="CSER") then
			charge=0.0;   helix=-0.15;   sheet=-1.45;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="THR".or.group(chainID,ic)%gtype=="NTHR".or.group(chainID,ic)%gtype=="CTHR") then
			charge=0.0;   helix=-0.39;   sheet=0.70;    pat=-1.0
		elseif(group(chainID,ic)%gtype=="HIE".or.group(chainID,ic)%gtype=="NHIE".or.group(chainID,ic)%gtype=="CHIE") then
			charge=0.0;   helix=0.11;    sheet=-1.34;   pat=-1.0
		elseif(group(chainID,ic)%gtype=="CYS".or.group(chainID,ic)%gtype=="NCYS".or.group(chainID,ic)%gtype=="CCYS") then
			charge=0.0;   helix=-0.57;   sheet=0.63;    pat=1.0
		elseif(group(chainID,ic)%gtype=="PRO".or.group(chainID,ic)%gtype=="NPRO".or.group(chainID,ic)%gtype=="CPRO") then
			charge=0.0;   helix=0.0;     sheet=0.0;     pat=1.0
		elseif(group(chainID,ic)%gtype=="ACE") then
			goto 10
        elseif(group(chainID,ic)%gtype=="NME") then
			goto 10
        elseif(group(chainID,ic)%gtype=="NHE") then
			goto 10
		endif

		CH=CH+charge
		Alpha=Alpha+helix
		Belta=Belta+sheet
		
		if(pat/=pat_old) then
			Pattern=Pattern+1.0
			pat_old=pat
		endif
10      continue		
	enddo
		
	p=coeff_CH*abs(CH)+(coeff_Alpha*Alpha+coeff_Belta*Belta+coeff_Pat*Pattern)/gnum
	Pagg=Pagg+p
	
	return
	end subroutine aggregation_propensity

	
	subroutine vdwcontri(x,y,rxy,epsion,r,vdw)
	implicit none
	integer					:: x, y
	real					:: rxy, vdw
	real					:: epsion_xy, r_xy
	real					:: acoeff, bcoeff
	real					:: epsion(repeated_unit*atom_num), r(repeated_unit*atom_num)

	epsion_xy=sqrt(epsion(x)*epsion(y))
	r_xy=r(x)+r(y)
	acoeff=epsion_xy*(r_xy**12)
	bcoeff=epsion_xy*2*(r_xy**6)
	vdw=acoeff/(rxy**6)-bcoeff/(rxy**3)

	return
	end subroutine vdwcontri

	
	subroutine	elecontri(x,y,rxy,charge,dielecons,ele)
	implicit none
	integer					:: x, y
	real					:: rxy,ele
	real					:: qx, qy, dielecons4solute
	real					:: charge(repeated_unit*atom_num),dielecons(repeated_unit*atom_num)

	qx=charge(x)
	qy=charge(y)
	if(dielecons(x).ge.dielecons(y)) then
		dielecons4solute=dielecons(x)
	else
		dielecons4solute=dielecons(y)
	endif
	ele=(qx*qy)/(dielecons4solute*sqrt(rxy))
	
	return
	end subroutine elecontri

	
	subroutine	sgbcontri(x,y,rxy,alphax,alphay,charge,dielecons,sgb)
	implicit none
	integer					:: x, y
	real					:: rxy, sgb
	real					:: dielecons4water
	real					:: fgb, alphax, alphay
	real					:: sgb1, sgb2
	real					:: qx, qy, dielecons4solute
	real					:: charge(repeated_unit*atom_num), dielecons(repeated_unit*atom_num)

	dielecons4water=80.0
	qx=charge(x)
	qy=charge(y)
	fgb=sqrt(rxy+alphax*alphay*exp(-rxy/(4*alphax*alphay)))
	if(dielecons(x).ge.dielecons(y)) then
		dielecons4solute=dielecons(x)
	else
		dielecons4solute=dielecons(y)
	endif
	sgb1=(1.0/dielecons4solute)-(1.0/dielecons4water)
	sgb2=qx*qy/fgb
	sgb=-sgb1*sgb2

	return
	end subroutine sgbcontri

	
	subroutine effgbradi(xstart,xend,x,rborn,fs,mdcrd,alphax)
	implicit none
	integer					:: x
	integer					:: xstart, xend
	real					:: alpha, beta, gamma
	real					:: redborn
	real					:: psi, integra, alphax
	real					:: rborn(repeated_unit*atom_num),fs(repeated_unit*atom_num),mdcrd(repeated_unit*atom_num,3)

	alpha=0.8
	beta=0.0
	gamma=2.91
	redborn=rborn(x)-0.09
	call areafract(xstart,xend,x,rborn,fs,mdcrd,integra)
	psi=integra*redborn
	alphax=1.0/(1.0/redborn-tanh(alpha*psi-beta*psi*psi+gamma*psi*psi*psi)/rborn(x))

	return
	end subroutine effgbradi

	
	subroutine areafract(xstart,xend,x,rborn,fs,mdcrd,integra)
	implicit none
	integer					:: x, y
	integer					:: xstart, xend
	real					:: integra, rxy, sum, redborn
	real					:: lxy, uxy
	real					:: rborn(repeated_unit*atom_num), fs(repeated_unit*atom_num), mdcrd(repeated_unit*atom_num,3)

	integra=0.0
	do y=xstart, xend
		if(y.ne.x) then
			rxy=(mdcrd(x,1)-mdcrd(y,1))**2+(mdcrd(x,2)-mdcrd(y,2))**2+(mdcrd(x,3)-mdcrd(y,3))**2
			rxy=sqrt(rxy)
			redborn=fs(y)*(rborn(y)-0.09)
			call evalualij(x,rxy,redborn,rborn,lxy)
			call evaluauij(x,rxy,redborn,rborn,uxy)
			sum=(1.0/lxy)-(1.0/uxy)+(1.0/(uxy*uxy)-1.0/(lxy*lxy))*rxy/4.0+log(lxy/uxy)/(2.0*rxy)+ &
				(1.0/(lxy*lxy)-1.0/(uxy*uxy))*redborn*redborn/(4*rxy)
			integra=integra+sum
		endif
	enddo
	integra=integra/2.0

	return
	end	subroutine areafract

	
	subroutine evalualij(x,rxy,redborn,rborn,lxy)
	implicit none
	integer					:: x
	real					:: rxy, redborn
	real					:: lxy
	real					:: rborn(repeated_unit*atom_num)

	if(rborn(x).le.(rxy-redborn)) then
		lxy=rxy-redborn
	elseif(rborn(x).le.(rxy+redborn)) then
		lxy=rborn(x)-0.09
	else
		lxy=1.0
	endif

	return
	end subroutine evalualij

	
	subroutine evaluauij(x,rxy,redborn,rborn,uxy)
	implicit none
	integer					:: x
	real					:: rxy, redborn, redborn1
	real					:: uxy
	real					:: rborn(repeated_unit*atom_num)

	if(rborn(x).lt.(rxy+redborn)) then
		uxy=rxy+redborn
	else
		uxy=1.0
	endif

	return
	end subroutine evaluauij


	subroutine snp_energy(start_atom,end_atom,rborn,mdcrd,snp)
	implicit none
	integer							:: start_atom, end_atom
	real							:: rborn(:), mdcrd(:,:)
	real(kind=8)					:: volume, surfarea, snp

	call arvo(start_atom,end_atom,rborn,mdcrd,volume,surfarea)
	snp=surftens*surfarea+offsetvalue

	return
	end subroutine snp_energy


	subroutine build_sasa_cache(natom,rborn,mdcrd,atom_chain,atom_site,atom_group, &
		atom_type,cache,binding_snp)
	implicit none
	integer, intent(in) :: natom
	real, intent(in) :: rborn(:),mdcrd(:,:)
	integer, intent(in) :: atom_chain(:),atom_site(:),atom_group(:)
	character(len=4), intent(in) :: atom_type(:)
	type(sasa_cache_type), intent(inout) :: cache
	real(kind=8), intent(out) :: binding_snp
	integer :: i,categoryID,n_group,local_index
	integer, allocatable :: global_index(:)
	real, allocatable :: group_rborn(:),group_coord(:,:)
	real(kind=8), allocatable :: group_atom_area(:)
	real(kind=8) :: volume,area,tolerance

	cache%valid=.false.
	cache%n_atoms=natom
	if (allocated(cache%area_group)) deallocate(cache%area_group)
	if (allocated(cache%atom_area_complex)) deallocate(cache%atom_area_complex)
	if (allocated(cache%atom_area_group)) deallocate(cache%atom_area_group)
	if (allocated(cache%coord)) deallocate(cache%coord)
	if (allocated(cache%sasa_radius)) deallocate(cache%sasa_radius)
	if (allocated(cache%chain_id)) deallocate(cache%chain_id)
	if (allocated(cache%site_id)) deallocate(cache%site_id)
	if (allocated(cache%group_id)) deallocate(cache%group_id)
	if (allocated(cache%atom_name)) deallocate(cache%atom_name)
	allocate(cache%area_group(num4category))
	allocate(cache%atom_area_complex(natom),cache%atom_area_group(natom))
	allocate(cache%coord(natom,3),cache%sasa_radius(natom))
	allocate(cache%chain_id(natom),cache%site_id(natom),cache%group_id(natom))
	allocate(cache%atom_name(natom))

	cache%coord=dble(mdcrd(1:natom,1:3))
	cache%sasa_radius=dble(rborn(1:natom))+1.40d0
	cache%chain_id=atom_chain(1:natom)
	cache%site_id=atom_site(1:natom)
	cache%group_id=atom_group(1:natom)
	cache%atom_name=atom_type(1:natom)
	cache%area_group=0.0d0
	cache%atom_area_group=0.0d0

	call arvo(1,natom,rborn,mdcrd,volume,cache%area_complex,cache%atom_area_complex)
	tolerance=max(1.0d-8,1.0d-10*abs(cache%area_complex))
	if (abs(sum(cache%atom_area_complex)-cache%area_complex) > tolerance) &
		stop 'ERROR: per-atom complex SASA does not reproduce the total SASA.'

	do categoryID=1,num4category
		n_group=count(atom_group(1:natom) == categoryID)
		allocate(group_rborn(n_group),group_coord(n_group,3))
		allocate(group_atom_area(n_group),global_index(n_group))
		local_index=0
		do i=1,natom
			if (atom_group(i) /= categoryID) cycle
			local_index=local_index+1
			global_index(local_index)=i
			group_rborn(local_index)=rborn(i)
			group_coord(local_index,:)=mdcrd(i,:)
		enddo
		call arvo(1,n_group,group_rborn,group_coord,volume,area,group_atom_area)
		cache%area_group(categoryID)=area
		do i=1,n_group
			cache%atom_area_group(global_index(i))=group_atom_area(i)
		enddo
		tolerance=max(1.0d-8,1.0d-10*abs(area))
		if (abs(sum(group_atom_area)-area) > tolerance) &
			stop 'ERROR: per-atom group SASA does not reproduce the total SASA.'
		deallocate(group_rborn,group_coord,group_atom_area,global_index)
	enddo

	binding_snp=surftens*(cache%area_complex-sum(cache%area_group))/ &
		real(repeated_unit,kind=8)
	cache%valid=.true.
	end subroutine build_sasa_cache


	subroutine calculate_atom_sasa(atom,natom,rborn,mdcrd,atom_group,system_id,area)
	implicit none
	integer, intent(in) :: atom,natom,atom_group(:),system_id
	real, intent(in) :: rborn(:),mdcrd(:,:)
	real(kind=8), intent(out) :: area
	integer :: i,n_context
	real(kind=8) :: dx,dy,dz,cutoff2,volume,total_area
	real, allocatable :: local_rborn(:),local_coord(:,:)
	real(kind=8), allocatable :: local_area(:)

	allocate(local_rborn(natom),local_coord(natom,3),local_area(natom))
	n_context=1
	local_rborn(1)=rborn(atom)
	local_coord(1,:)=mdcrd(atom,:)
	! Context-only atoms occlude the updated atom but their cached areas are
	! not replaced. Direct sphere overlap is the exact SASA context criterion.
	do i=1,natom
		if (i == atom) cycle
		if (system_id > 0 .and. atom_group(i) /= system_id) cycle
		dx=dble(mdcrd(atom,1)-mdcrd(i,1))
		dy=dble(mdcrd(atom,2)-mdcrd(i,2))
		dz=dble(mdcrd(atom,3)-mdcrd(i,3))
		cutoff2=(dble(rborn(atom)+rborn(i))+2.80d0)**2
		if (dx*dx+dy*dy+dz*dz >= cutoff2) cycle
		n_context=n_context+1
		local_rborn(n_context)=rborn(i)
		local_coord(n_context,:)=mdcrd(i,:)
	enddo
	call arvo(1,n_context,local_rborn,local_coord,volume,total_area,local_area)
	area=local_area(1)
	deallocate(local_rborn,local_coord,local_area)
	end subroutine calculate_atom_sasa


	subroutine update_sasa_cache(natom,rborn,mdcrd,atom_chain,atom_site,atom_group, &
		atom_type,old_cache,new_cache,binding_snp)
	implicit none
	integer, intent(in) :: natom,atom_chain(:),atom_site(:),atom_group(:)
	real, intent(in) :: rborn(:),mdcrd(:,:)
	character(len=4), intent(in) :: atom_type(:)
	type(sasa_cache_type), intent(in) :: old_cache
	type(sasa_cache_type), intent(inout) :: new_cache
	type(sasa_cache_type) :: full_cache
	real(kind=8), intent(out) :: binding_snp
	integer :: i,j,system_id,n_changed
	integer, allocatable :: old_match(:),new_match(:)
	logical, allocatable :: changed_old(:),changed_new(:),updated_old(:),updated_new(:)
	real(kind=8) :: dx,dy,dz,cutoff2,new_radius,area,full_snp,tolerance

	if (.not.old_cache%valid) then
		call build_sasa_cache(natom,rborn,mdcrd,atom_chain,atom_site,atom_group, &
			atom_type,new_cache,binding_snp)
		return
	endif

	allocate(old_match(old_cache%n_atoms),new_match(natom))
	allocate(changed_old(old_cache%n_atoms),changed_new(natom))
	allocate(updated_old(old_cache%n_atoms),updated_new(natom))
	old_match=0
	new_match=0
	do i=1,old_cache%n_atoms
		do j=1,natom
			if (old_cache%chain_id(i) /= atom_chain(j)) cycle
			if (old_cache%site_id(i) /= atom_site(j)) cycle
			if (old_cache%atom_name(i) /= atom_type(j)) cycle
			old_match(i)=j
			new_match(j)=i
			exit
		enddo
	enddo

	changed_old=.false.
	changed_new=.false.
	do i=1,old_cache%n_atoms
		j=old_match(i)
		if (j == 0) then
			changed_old(i)=.true.
			cycle
		endif
		new_radius=dble(rborn(j))+1.40d0
		dx=old_cache%coord(i,1)-dble(mdcrd(j,1))
		dy=old_cache%coord(i,2)-dble(mdcrd(j,2))
		dz=old_cache%coord(i,3)-dble(mdcrd(j,3))
		if (dx*dx+dy*dy+dz*dz > 1.0d-12 .or. &
			abs(old_cache%sasa_radius(i)-new_radius) > 1.0d-10) then
			changed_old(i)=.true.
			changed_new(j)=.true.
		endif
	enddo
	do j=1,natom
		if (new_match(j) == 0) changed_new(j)=.true.
	enddo
	n_changed=count(changed_old)+count(changed_new)
	if (n_changed > max(40,(old_cache%n_atoms+natom)/10)) then
		call build_sasa_cache(natom,rborn,mdcrd,atom_chain,atom_site,atom_group, &
			atom_type,new_cache,binding_snp)
		deallocate(old_match,new_match,changed_old,changed_new,updated_old,updated_new)
		return
	endif

	new_cache%valid=.false.
	new_cache%n_atoms=natom
	if (allocated(new_cache%area_group)) deallocate(new_cache%area_group)
	if (allocated(new_cache%atom_area_complex)) deallocate(new_cache%atom_area_complex)
	if (allocated(new_cache%atom_area_group)) deallocate(new_cache%atom_area_group)
	if (allocated(new_cache%coord)) deallocate(new_cache%coord)
	if (allocated(new_cache%sasa_radius)) deallocate(new_cache%sasa_radius)
	if (allocated(new_cache%chain_id)) deallocate(new_cache%chain_id)
	if (allocated(new_cache%site_id)) deallocate(new_cache%site_id)
	if (allocated(new_cache%group_id)) deallocate(new_cache%group_id)
	if (allocated(new_cache%atom_name)) deallocate(new_cache%atom_name)
	allocate(new_cache%area_group(num4category))
	allocate(new_cache%atom_area_complex(natom),new_cache%atom_area_group(natom))
	allocate(new_cache%coord(natom,3),new_cache%sasa_radius(natom))
	allocate(new_cache%chain_id(natom),new_cache%site_id(natom),new_cache%group_id(natom))
	allocate(new_cache%atom_name(natom))
	new_cache%coord=dble(mdcrd(1:natom,1:3))
	new_cache%sasa_radius=dble(rborn(1:natom))+1.40d0
	new_cache%chain_id=atom_chain(1:natom)
	new_cache%site_id=atom_site(1:natom)
	new_cache%group_id=atom_group(1:natom)
	new_cache%atom_name=atom_type(1:natom)

	do system_id=0,num4category
		! Build independent update sets for the complete complex and each
		! isolated energy group. The context shell is added inside
		! calculate_atom_sasa when an updated contribution is recalculated.
		updated_old=changed_old
		updated_new=changed_new
		do i=1,old_cache%n_atoms
			if (system_id > 0 .and. old_cache%group_id(i) /= system_id) cycle
			do j=1,old_cache%n_atoms
				if (.not.changed_old(j)) cycle
				if (system_id > 0 .and. old_cache%group_id(j) /= system_id) cycle
				dx=old_cache%coord(i,1)-old_cache%coord(j,1)
				dy=old_cache%coord(i,2)-old_cache%coord(j,2)
				dz=old_cache%coord(i,3)-old_cache%coord(j,3)
				cutoff2=(old_cache%sasa_radius(i)+old_cache%sasa_radius(j))**2
				if (dx*dx+dy*dy+dz*dz < cutoff2) updated_old(i)=.true.
			enddo
		enddo
		do i=1,natom
			if (system_id > 0 .and. atom_group(i) /= system_id) cycle
			do j=1,natom
				if (.not.changed_new(j)) cycle
				if (system_id > 0 .and. atom_group(j) /= system_id) cycle
				dx=dble(mdcrd(i,1)-mdcrd(j,1))
				dy=dble(mdcrd(i,2)-mdcrd(j,2))
				dz=dble(mdcrd(i,3)-mdcrd(j,3))
				cutoff2=(dble(rborn(i)+rborn(j))+2.80d0)**2
				if (dx*dx+dy*dy+dz*dz < cutoff2) updated_new(i)=.true.
			enddo
		enddo
		do i=1,old_cache%n_atoms
			j=old_match(i)
			if (j == 0) cycle
			if (updated_old(i)) updated_new(j)=.true.
			if (updated_new(j)) updated_old(i)=.true.
		enddo

		do i=1,natom
			if (system_id > 0 .and. atom_group(i) /= system_id) cycle
			j=new_match(i)
			if (system_id == 0) then
				if (.not.updated_new(i) .and. j > 0) then
					new_cache%atom_area_complex(i)=old_cache%atom_area_complex(j)
				else
					call calculate_atom_sasa(i,natom,rborn,mdcrd,atom_group,0,area)
					new_cache%atom_area_complex(i)=area
				endif
			else
				if (.not.updated_new(i) .and. j > 0) then
					new_cache%atom_area_group(i)=old_cache%atom_area_group(j)
				else
					call calculate_atom_sasa(i,natom,rborn,mdcrd,atom_group,system_id,area)
					new_cache%atom_area_group(i)=area
				endif
			endif
		enddo
	enddo

	new_cache%area_complex=sum(new_cache%atom_area_complex)
	do system_id=1,num4category
		new_cache%area_group(system_id)=0.0d0
		do i=1,natom
			if (atom_group(i) == system_id) &
				new_cache%area_group(system_id)=new_cache%area_group(system_id)+new_cache%atom_area_group(i)
		enddo
	enddo
	binding_snp=surftens*(new_cache%area_complex-sum(new_cache%area_group))/ &
		real(repeated_unit,kind=8)
	new_cache%valid=.true.
	if (sasa_validation_count < 5) then
		call build_sasa_cache(natom,rborn,mdcrd,atom_chain,atom_site,atom_group, &
			atom_type,full_cache,full_snp)
		tolerance=max(1.0d-6,1.0d-8*abs(full_cache%area_complex))
		if (abs(new_cache%area_complex-full_cache%area_complex) > tolerance .or. &
			any(abs(new_cache%area_group-full_cache%area_group) > tolerance)) then
			stop 'ERROR: incremental SASA cache does not match a full ARVO calculation.'
		endif
		sasa_validation_count=sasa_validation_count+1
	endif
	deallocate(old_match,new_match,changed_old,changed_new,updated_old,updated_new)
	end subroutine update_sasa_cache


	subroutine arvo(start_atom,end_atom,rborn,mdcrd,volume,surfarea,atom_area)
	implicit none
	integer							:: start_atom, end_atom, ks, kl, ka, ki, i, index, np_test, natom
	real							:: rborn(:), mdcrd(:,:)
	real(kind=8)					:: sa, volume, surfarea
	real(kind=8), intent(out), optional :: atom_area(:)
	parameter (ks=5000,kl=500,ka=500,ki=250000)
	integer							:: neighbors_number(ks), index_start(ks), neighbors_indices(ki)
	real(kind=8)					:: spheres(ks,4), av(2)

	natom=end_atom-start_atom+1
	if(natom.ge.ks) then
		open(10, file="error.txt", access="append")
			write(10,*) "ks is too small to perform ESURF calculation!"
			write(10,*) "Increase in function then recompile!"
		close(10)
		stop
	endif

	do i=start_atom,end_atom
		index=i-start_atom+1
		spheres(index,1)=dble(mdcrd(i,1))
		spheres(index,2)=dble(mdcrd(i,2))
		spheres(index,3)=dble(mdcrd(i,3))
		spheres(index,4)=dble(rborn(i))+1.40d0
	enddo

    call make_neighbors(1,natom,spheres,neighbors_number,index_start,neighbors_indices,ks,kl,natom,ki)

	np_test=0
	do while(np_test.eq.0)
		call NPTest(1,natom,spheres,neighbors_number,index_start,neighbors_indices,ks,ki,np_test)
		if(np_test.eq.0) then
			sa=0.324d0
			call spheres_rotation(spheres,ks,natom,sa)
		endif
	enddo

	volume=0d0
	surfarea=0d0
	do i=1,natom
	    call areavolume(i,spheres,neighbors_number,index_start,neighbors_indices,ks,kl,ka,ki,av)
	    volume=volume+av(1)
	    surfarea=surfarea+av(2)
		if (present(atom_area)) atom_area(i)=av(2)
	enddo

	return
	end subroutine arvo


	subroutine make_neighbors(i1,i2,spheres,neighbors_number,index_start,neighbors_indices,ks,kl,ns,ki)
	implicit none
	integer							:: i1, i2, ns, ks, kl, ki, i, j, neighbors
	integer							:: neighbors_number(ks), index_start(ks), neighbors_indices(ki), ind(kl)
	real(kind=8)					:: spheres(ks,4)

    index_start(i1)=1
    do i=i1,i2
		call NEIGHBOR(i,spheres,ind,ks,kl,ns,neighbors)
		neighbors_number(i)=neighbors
		if(neighbors_number(i).gt.kl) then
			open(10, file="error.txt", access="append")
				write(10,*) "Too small of kl!"
			close(10)
			stop
		endif
		if (neighbors_number(i).le.0) then
			index_start(i+1)=index_start(i)
		else
			index_start(i+1)=index_start(i)+neighbors_number(i)
			do j=1,neighbors_number(i)
				neighbors_indices(index_start(i)+j-1)=ind(j)
			enddo
		endif
    enddo

	return
	end subroutine make_neighbors


	subroutine NEIGHBOR(i,spheres,ind,ks,kl,ns,neighbors)
	implicit none
	integer							:: i, k, kl, ks, ns, neighbors
	integer							:: neighbors_num
	real(kind=8)					:: xi, yi, zi, ri, rk, dd
	integer							:: ind(kl)
	real(kind=8)					:: spheres(ks,4)

    neighbors_num=0
    xi=spheres(i,1)
	yi=spheres(i,2)
	zi=spheres(i,3)
	ri=spheres(i,4)
    do k=1,ns
		if (k.ne.i) then
			if(dabs(xi-spheres(k,1)).lt.ri+spheres(k,4)) then
				dd=dsqrt((xi-spheres(k,1))**2+(yi-spheres(k,2))**2+(zi-spheres(k,3))**2)
				rk=spheres(k,4)
				if (dd.lt.ri+rk) then
					if (dd+ri.le.rk) then
						neighbors_num=-1
						exit
					elseif (dd+rk.gt.ri) then
						neighbors_num=neighbors_num+1
						ind(neighbors_num)=k
					endif
				endif
			endif
		endif
      enddo
	neighbors=neighbors_num

	return
	end subroutine NEIGHBOR


	subroutine NPTest(i1,i2,spheres,neighbors_number,index_start,neighbors_indices,ks,ki,np_test)
	implicit none
	integer							:: i1, i2, ks, ki, i, k, ink, npt, np_test
	real(kind=8)					:: eps_north_pole, dmin, d
	integer							:: neighbors_number(ks), index_start(ks), neighbors_indices(ki)
	real(kind=8)					:: spheres(ks,4)

    eps_north_pole=1d-4
	dmin=10000d0
	do i=i1,i2
		do k=1,neighbors_number(i)
			ink=neighbors_indices(index_start(i)+k-1)
			d=dabs(dsqrt((spheres(i,1)-spheres(ink,1))**2+(spheres(i,2)-spheres(ink,2))**2 &
     		  +(spheres(i,3)+spheres(i,4)-spheres(ink,3))**2)-spheres(ink,4))
			if (d.lt.dmin) then
				dmin=d
			endif
		enddo
	enddo

    if (dmin.lt.eps_north_pole) then
		npt=0
	else
		npt=1
	endif
 	np_test=npt

    return
	end subroutine NPTest


    subroutine spheres_rotation(spheres,ks,ns,sa)
	implicit none
	integer							:: ks, ns, i
	real(kind=8)					:: sa, ca, x, z
	real(kind=8)					:: spheres(ks,4)

	ca=dsqrt(1d0-sa*sa)
	do i=1,ns
 		x=spheres(i,1)
		z=spheres(i,3)
		spheres(i,1)=ca*x-sa*z
		spheres(i,3)=sa*x+ca*z
	enddo

	return
	end subroutine spheres_rotation


	subroutine areavolume(i,spheres,neighbors_number,index_start,neighbors_indices,ks,kl,ka,ki,av)
	implicit none
	integer							:: i, ks, kl, ka, ki
	integer							:: circles_to_arcs, j, npos, narcs
	integer							:: nls
	real(kind=8)					:: z1, r1
	integer							:: neighbors_number(ks), index_start(ks), neighbors_indices(ki), ind(kl)
	real(kind=8)					:: spheres(ks,4), av(2), avi(2)
	real(kind=8)					:: circles(kl,4), arcs(ka,3), sphere_local(kl,4)
	real(kind=8), parameter			:: pi4arvo=3.1415926535897932384626433832795d0

	if (neighbors_number(i).lt.0) then
		av(1)=0d0
		av(2)=0d0
    elseif (neighbors_number(i).eq.0) then
		av(1)=4d0*pi4arvo*spheres(i,4)**3/3.d0
		av(2)=4d0*pi4arvo*spheres(i,4)**2
	else
		nls=neighbors_number(i)+1
		ind(1)=i
		do j=1,(nls-1)
			ind(j+1)=neighbors_indices(index_start(i)+j-1)
		enddo
		call local_spheres(spheres,ind,sphere_local,nls,ks,kl)
		av(1)=0d0
		av(2)=0d0
		call make_ts_circles(sphere_local,circles,kl,nls)
		call CirclestoArcs(circles,arcs,kl,nls,ka,circles_to_arcs)
		narcs=circles_to_arcs
		npos=0
		do j=1,(nls-1)
			if (circles(j,4).gt.0) then
				npos=npos+1
			endif
		enddo
		z1=sphere_local(1,3)
		r1=sphere_local(1,4)

		if (npos.gt.0) then
			call avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi)
			av(1)=av(1)+avi(1)
			av(2)=av(2)+avi(2)
         else
			call avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi)
			av(1)=av(1)+avi(1)+4d0*pi4arvo*sphere_local(1,4)**3/3d0
			av(2)=av(2)+avi(2)+4d0*pi4arvo*sphere_local(1,4)**2
         endif
    endif

	return
	end subroutine areavolume


    subroutine local_spheres(spheres,ind,sphere_local,nls,ks,kl)
	implicit none
	integer							:: kl, ks, i, j
	integer							:: nls
	integer							:: ind(kl)
	real(kind=8)					:: spheres(ks,4), sphere_local(kl,4)

	do i=1,nls
		do j=1,4
			sphere_local(i,j)=spheres(ind(i),j)
	    enddo
    enddo

	return
	end subroutine local_spheres


	subroutine make_ts_circles(sphere_local,circles,kl,nls)
	implicit none
	integer							:: kl, k
	integer							:: nls
	real(kind=8)					:: r1, dx, dy, a, b, c, d
	real(kind=8)					:: circles(kl,4), sphere_local(kl,4)

	r1=sphere_local(1,4)
	do k=1,(nls-1)
		dx=sphere_local(1,1)-sphere_local(k+1,1)
		dy=sphere_local(1,2)-sphere_local(k+1,2)
		a=dx*dx+dy*dy+(sphere_local(1,3)+r1-sphere_local(k+1,3))**2-sphere_local(k+1,4)**2
		b=8d0*r1*r1*dx
		c=8d0*r1*r1*dy
		d=4d0*r1*r1*(dx*dx+dy*dy+(sphere_local(1,3)-r1-sphere_local(k+1,3))**2-sphere_local(k+1,4)**2)
		circles(k,1)=-b/(2d0*a)
		circles(k,2)=-c/(2d0*a)
		circles(k,3)=dsqrt((b*b+c*c-4d0*a*d)/(4d0*a*a))
		if (a.gt.0) then
			circles(k,4)=-1
		else
			circles(k,4)=1
		endif
	enddo

	return
	end subroutine make_ts_circles


	subroutine CirclestoArcs(circles,arcs,kl,nls,ka,circles_to_arcs)
	implicit none
	integer							:: kl, ka, number_arc, i, j, k
	integer							:: nna, new_arcs, circles_to_arcs
	integer							:: nls
	real(kind=8)					:: circles(kl,4), arcs(ka,3), arcsnew(ka,3)
	real(kind=8), parameter			:: pi4arvo=3.1415926535897932384626433832795d0

    number_arc=0
	if (nls.eq.2) then
		number_arc=1
		arcs(1,1)=1
        arcs(1,2)=0d0
        arcs(1,3)=2d0*pi4arvo*circles(1,4)
	else
		do i=1,(nls-1)
			call NewArcs(i,circles,arcsnew,kl,ka,nls,new_arcs)
			nna=new_arcs
			if (nna.gt.0) then
				do j=1,nna
					do k=1,3
						arcs(number_arc+j,k)=arcsnew(j,k)
					enddo
				enddo
				number_arc=number_arc+nna
			endif
		enddo
	endif
	circles_to_arcs=number_arc

	return
	end subroutine CirclestoArcs


	subroutine NewArcs(i,circles,arcsnew,kl,ka,nls,new_arcs)
	implicit none
	integer							:: kl, ka, i, j, jj
	integer							:: new_arcs
	integer							:: circle_in_circle, point_in_circle, delete_equal
	integer							:: num_arc, num_angle, number_cond, na
	integer							:: nls
	real(kind=8)					:: ti, si, ri, t, s, r, d
	real(kind=8)					:: a1, a2, b1, b2
	real(kind=8)					:: circles(kl,4), arcsnew(ka,3), angles(ka)
	real(kind=8), parameter			:: pi4arvo=3.1415926535897932384626433832795d0

	num_arc=0
	num_angle=0
	ti=circles(i,1)
	si=circles(i,2)
	ri=circles(i,3)
	do j=1,(nls-1)
		if (j.ne.i) then
			t=circles(j,1)
			s=circles(j,2)
			r=circles(j,3)
			d=dsqrt((ti-t)**2+(si-s)**2)
			if ( (d.lt.r+ri) .AND. (dabs(r-ri).lt.d) ) then
				call circles_intersection(i,j,circles,kl,a1,a2,b1,b2)
				angles(num_angle+1)=a1
                angles(num_angle+2)=a2
                num_angle=num_angle+2
			endif
		endif
	enddo
	if (num_angle.eq.0) then
		number_cond=0
	    do j=1,(nls-1)
			if (j.ne.i) then
				call CircleinCircle(i,j,circles,kl,circle_in_circle)
				number_cond=number_cond+circle_in_circle
			endif
		enddo
		if (number_cond.eq.(nls-2)) then
			num_arc=1
			arcsnew(1,1)=i
			arcsnew(1,2)=0d0
			arcsnew(1,3)=2d0*pi4arvo*circles(i,4)
		endif
	else
		if (circles(i,4).gt.0) then
			call mysort(angles,ka,num_angle)
		else
			call mydsort(angles,ka,num_angle)
		endif
		call DeleteEqual(angles,ka,num_angle,delete_equal)
		na=delete_equal
		num_angle=na
		do j=1,(na-1)
			number_cond=0
			do jj=1,(nls-1)
				if (jj.ne.i) then
					t=ti+ri*dcos((angles(j)+angles(j+1))/2d0)
					s=si+ri*dsin((angles(j)+angles(j+1))/2d0)
					call PointinCircle(t,s,jj,circles,kl,point_in_circle)
					number_cond=number_cond+point_in_circle
				endif
			enddo
			if (number_cond.eq.(nls-2)) then
				num_arc=num_arc+1
				arcsnew(num_arc,1)=i
				arcsnew(num_arc,2)=angles(j)
				arcsnew(num_arc,3)=angles(j+1)-angles(j)
			endif
		enddo
		number_cond=0
		do j=1,(nls-1)
			if (j.ne.i) then
				t=ti+ri*dcos((angles(1)+2d0*pi4arvo+angles(na))/2d0)
				s=si+ri*dsin((angles(1)+2d0*pi4arvo+angles(na))/2d0)
				call PointinCircle(t,s,j,circles,kl,point_in_circle)
				number_cond=number_cond+point_in_circle
			endif
		enddo
		if (number_cond.eq.(nls-2)) then
			num_arc=num_arc+1
			arcsnew(num_arc,1)=i
			arcsnew(num_arc,2)=angles(na)
			arcsnew(num_arc,3)=angles(1)+circles(i,4)*2d0*pi4arvo-angles(na)
		endif
	endif
	new_arcs=num_arc

	return
	end subroutine NewArcs


    subroutine circles_intersection(ic1,ic2,circles,kl,a1,a2,b1,b2)
	implicit none
	integer							:: ic1, ic2, kl
	real(kind=8)					:: a1, a2, b1, b2
	real(kind=8)					:: eps_deltat, t1, s1, r1, t2, s2, r2
	real(kind=8)					:: A, B, C, D
	real(kind=8)					:: circles(kl,4)
	real(kind=8), parameter			:: pi4arvo=3.1415926535897932384626433832795d0

	eps_deltat=1d-12
	t1=circles(ic1,1)
	s1=circles(ic1,2)
	r1=circles(ic1,3)
	t2=circles(ic2,1)
	s2=circles(ic2,2)
	r2=circles(ic2,3)
	if(dabs(t2-t1).lt.eps_deltat) then
		B=((r1*r1-r2*r2)/(s2-s1)-(s2-s1))/2d0
		A=dsqrt(r2*r2-B*B)
		if (B.eq.0) then
			b1=0d0
			b2=pi4arvo
		elseif (B.gt.0) then
			b1=datan(dabs(B/A))
			b2=pi4arvo-b1
		else
			b1=pi4arvo+datan(dabs(B/A))
			b2=3d0*pi4arvo-b1
		endif
		B=B+s2-s1
		if (B.eq.0) then
			a1=0d0
			a2=pi4arvo
		elseif (B.gt.0) then
			a1=datan(dabs(B/A))
			a2=pi4arvo-a1
		else
			a1=pi4arvo+datan(dabs(B/A))
			a2=3d0*pi4arvo-a1
		endif
	else
		C=((r1*r1-r2*r2-(s2-s1)**2)/(t2-t1)-(t2-t1))/2d0
		D=(s1-s2)/(t2-t1)
		B=(-C*D+dsqrt((D*D+1d0)*r2*r2-C*C))/(D*D+1d0)
		A=C+D*B
		if (A.eq.0) then
			if (B.gt.0) then
				b1=pi4arvo/2d0
			else
				b1=-pi4arvo/2d0
			endif
		elseif (A.gt.0) then
			b1=datan(B/A)
		else
			b1=pi4arvo+datan(B/A)
		endif
		B=B+s2-s1
		A=A+t2-t1
		if (A.eq.0) then
			if (B.gt.0) then
				a1=pi4arvo/2d0
			else
				a1=-pi4arvo/2d0
			endif
		elseif (A.gt.0) then
			a1=datan(B/A)
		else
			a1=pi4arvo+datan(B/A)
		endif
		B=(-C*D-dsqrt((D*D+1d0)*r2*r2-C*C))/(D*D+1d0)
		A=C+D*B
		if (A.eq.0) then
			if (B.gt.0) then
				b2=pi4arvo/2d0
			else
				b2=-pi4arvo/2d0
			endif
		elseif (A.gt.0) then
			b2=datan(B/A)
		else
			b2=pi4arvo+datan(B/A)
		endif
		B=B+s2-s1
		A=A+t2-t1
		if (A.eq.0) then
			if (B.gt.0) then
				a2=pi4arvo/2d0
			else
				a2=-pi4arvo/2d0
			endif
		elseif (A.gt.0) then
			a2=datan(B/A)
		else
			a2=pi4arvo+datan(B/A)
		endif
	endif
	if (a1.lt.0) a1=a1+2d0*pi4arvo
	if (a2.lt.0) a2=a2+2d0*pi4arvo
	if (b1.lt.0) b1=b1+2d0*pi4arvo
	if (b2.lt.0) b2=b2+2d0*pi4arvo

    return
	end subroutine circles_intersection


	subroutine CircleinCircle(i,k,circles,kl, circle_in_circle)
	implicit none
	integer							:: i, k, kl
	integer							:: circle_in_circle
	real(kind=8)					:: d
	real(kind=8)					:: circles(kl,4)

	d=dsqrt((circles(i,1)+circles(i,3)-circles(k,1))**2+(circles(i,2)-circles(k,2))**2)
	if (d.lt.circles(k,3)) then
		if (circles(k,4).gt.0) then
			circle_in_circle=1
		else
			circle_in_circle=0
		endif
	elseif (d.gt.circles(k,3)) then
		if (circles(k,4).gt.0) then
			circle_in_circle=0
		else
			circle_in_circle=1
		endif
	else
		d=dsqrt((circles(i,1)-circles(k,1))**2+(circles(i,2)-circles(k,2))**2)
		if (d.lt.circles(k,3)) then
			if (circles(k,4).gt.0) then
				circle_in_circle=1
			else
				circle_in_circle=0
			endif
		else
			if (circles(k,4).gt.0) then
				circle_in_circle=0
			else
				circle_in_circle=1
			endif
		endif
	endif

	return
	end	subroutine CircleinCircle


	subroutine PointinCircle(t,s,k,circles,kl,point_in_circle)
	implicit none
	integer							:: k, kl
	integer							:: point_in_circle
	real(kind=8)					:: t, s, d
	real(kind=8)					:: circles(kl,4)

	d=dsqrt((t-circles(k,1))**2+(s-circles(k,2))**2)
	if (d.lt.circles(k,3)) then
		if (circles(k,4).gt.0) then
			point_in_circle=1
		else
			point_in_circle=0
		endif
	else
		if (circles(k,4).gt.0) then
			point_in_circle=0
		else
			point_in_circle=1
		endif
	endif

	return
	end	subroutine PointinCircle


	subroutine mysort(angles,ka,num_angle)
	implicit none
	integer							:: ka, num_angle, i, ii, j
	real(kind=8)					:: amax
	real(kind=8)					:: angles(ka)

	do i=1,(num_angle-1)
		ii=i
		amax=angles(i)
		do j=i+1,num_angle
			if (amax.gt.angles(j)) then
				ii=j
				amax=angles(j)
			endif
 		enddo
		if (ii.ne.i) then
			angles(ii)=angles(i)
			angles(i)=amax
		endif
	enddo

	return
	end	subroutine mysort


	subroutine mydsort(angles,ka,num_angle)
	implicit none
	integer							:: ka, num_angle, i, ii, j
	real(kind=8)					:: amin
	real(kind=8)					:: angles(ka)

	do i=1,(num_angle-1)
		ii=i
		amin=angles(i)
		do j=i+1,num_angle
			if (amin.lt.angles(j)) then
				ii=j
				amin=angles(j)
			endif
 		enddo
		if (ii.ne.i) then
			angles(ii)=angles(i)
			angles(i)=amin
		endif
	enddo

	return
	end subroutine mydsort


	subroutine DeleteEqual(angles,ka,num_angle,delete_equal)
	implicit none
	integer							:: ka, num_angle, m, i
	integer							:: delete_equal
	real(kind=8)					:: eps_angle, angle
	real(kind=8)					:: angles(ka), anglesnew(ka)

	eps_angle=1d-12
   	m=1
	angle=angles(1)
	anglesnew(1)=angle
	do i=2,num_angle
		if (dabs(angles(i)-angle).gt.eps_angle) then
			angle=angles(i)
			m=m+1
			anglesnew(m)=angle
		endif
	enddo
	delete_equal=m
	do i=1,m
		angles(i)=anglesnew(i)
	enddo

	return
	end	subroutine DeleteEqual


	subroutine avintegral(circles,arcs,kl,ka,narcs,r1,z1,avi)
	implicit none
	integer							:: kl, ka, narcs, k
	real(kind=8)					:: z1, r1, eps_two_pi
	real(kind=8)					:: t, s, r, A, B, C, SS, rr
	real(kind=8)					:: vIone, vItwo, vIthree
	real(kind=8)					:: vJone, vJtwo, vJthree
	real(kind=8)					:: delta_vint, delta_aint
	real(kind=8)					:: be, al, sb, cb, sa, ca
	real(kind=8)					:: a1, a2, b1, b2
	real(kind=8)					:: circles(kl,4), arcs(ka,3), avi(2)
	real(kind=8), parameter			:: pi4arvo=3.1415926535897932384626433832795d0

	eps_two_pi=1d-12
	avi(1)=0d0
	avi(2)=0d0
	do k=1,narcs
	    t=circles(int(arcs(k,1)),1)
	    s=circles(int(arcs(k,1)),2)
	    r=circles(int(arcs(k,1)),3)
	    A=(4d0*r1*r1+t*t+s*s+r*r)/2d0
	    B=t*r
		C=s*r
		SS=dsqrt(A*A-B*B-C*C)
		rr=r*r-A
		if (dabs(dabs(arcs(k,3))-2d0*pi4arvo).lt.eps_two_pi) then
			vIone=2d0*pi4arvo/SS
			vItwo=2d0*pi4arvo*A/(SS**3)
			vIthree=pi4arvo*(2d0*A*A+B*B+C*C)/(SS**5)
			vJone=pi4arvo+rr/2d0*vIone
			vJtwo=(vIone+rr*vItwo)/4d0
			vJthree=(vItwo+rr*vIthree)/8d0
			delta_vint=(128d0*vJthree*r1**7+8d0*vJtwo*r1**5+2d0*vJone*r1**3)/3d0-8d0*r1**4*vJtwo*(z1+r1)
	        delta_aint=2d0*vJone*r1**2
			if (arcs(k,3).lt.0) then
			   delta_vint=-delta_vint
			   delta_aint=-delta_aint
              endif
			avi(1)=avi(1)+delta_vint
			avi(2)=avi(2)+delta_aint
		else
			if (arcs(k,3).lt.0) then
				al=arcs(k,2)+arcs(k,3)
				be=arcs(k,2)
			else
				be=arcs(k,2)+arcs(k,3)
				al=arcs(k,2)
			endif
			vIone=2d0*(pi4arvo/2d0-datan((A*dcos((be-al)/2d0)+B*dcos((al+be)/2d0)+C*dsin((al+be)/2d0))/ &
    			  (SS*dsin((be-al)/2d0))))/SS
			sb=dsin(be)
			cb=dcos(be)
			sa=dsin(al)
			ca=dcos(al)
			call Fract(A,B,C,sa,ca,1,a1)
			call Fract(A,B,C,sa,ca,2,a2)
			call Fract(A,B,C,sb,cb,1,b1)
			call Fract(A,B,C,sb,cb,2,b2)
			vItwo=(b1-a1+A*vIone)/(SS*SS)
			vIthree=(b2-a2+(b1-a1)/A+(2d0*A*A+B*B+C*C)*vItwo/A)/(2d0*s*s)
			vJone=((be-al)+rr*vIone)/2d0
			vJtwo=(vIone+rr*vItwo)/4d0
			vJthree=(vItwo+rr*vIthree)/8d0
			delta_vint=(128d0*vJthree*r1**7+8d0*vJtwo*r1**5+2d0*vJone*r1**3)/3d0-8d0*r1**4*vJtwo*(z1+r1)
	        delta_aint=2d0*vJone*r1**2
			if (arcs(k,3).lt.0) then
			   delta_vint=-delta_vint
			   delta_aint=-delta_aint
            endif
			avi(1)=avi(1)+delta_vint
			avi(2)=avi(2)+delta_aint
		endif
	enddo

	return
	end	subroutine avintegral


	subroutine Fract(A,B,C,sinphi,cosphi,k,fraction)
	implicit none
	integer							:: k
	real(kind=8)					:: A, B, C, sinphi, cosphi
	real(kind=8)					:: fraction

    fraction=(-B*sinphi+C*cosphi)/(A+B*cosphi+C*sinphi)**k

	return
	end	subroutine Fract

end module energy_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module advancedfunction

	use constant
	use randomgenerator
	use input
	use mathfunction
	use database
	use transplant
	use energy_calculation

    contains

    

	subroutine scmf_substitution(group,sub_cycle,tgroup,avoid_site,avoid_aminoacid)
	implicit none
	integer, parameter :: max_sequence_attempts=10000
	integer, parameter :: max_avoid_attempts=10
	character(len=3), parameter :: aa_names(20)=(/ &
		'GLY','ALA','VAL','LEU','ILE','MET','PHE','TYR','TRP','SER', &
		'ASN','GLN','THR','HIE','ARG','LYS','GLU','ASP','CYS','PRO' /)
	integer, parameter :: aa_category(20)=(/ &
		1,1,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4 /)
	type(groupdetails), intent(in) :: group(repeated_unit,gnum)
	type(groupdetails), intent(out) :: tgroup(repeated_unit,gnum)
	integer, intent(out) :: sub_cycle
	integer, intent(in) :: avoid_site
	character(len=4), intent(in) :: avoid_aminoacid
	integer :: attempt,i,j,m,site,chainID,rotanum,ip,feedback
	integer :: aa_index,category,h_count,p_count,c_count,o_count
	integer :: free_site_count,swap_value,choice,changed_count
	integer :: minimum_changes,candidate_changes,tie_count
	integer :: grouped_preserved,free_preserved,lower_sum,upper_sum
	integer :: preserve_gain,remaining_needed,eligible_count
	integer :: repeated_failed_target_count
	integer :: current_aa(gnum),free_site_ID(gnum),free_site_order(gnum),changed_site(gnum)
	integer :: current_free_count(20),fixed_aa_count(20),grouped_aa_count(20)
	integer :: locked_aa_count(20),target_free_count(20),remaining_aa_count(20)
	integer :: trial_aa_count(20),aa_min(20),aa_max(20),aa_lower(20),aa_upper(20)
	integer :: locked_category_count(4),target_category_count(4),trial_category_count(4)
	integer :: free_category_count(4),eligible_aa(20)
	integer :: subset_position(10),best_subset_position(10),selected_aa(10)
	integer :: grouped_site_order(maximum_nmr_site_num),grouped_target_aa(maximum_nmr_site_num)
	logical :: is_free_site(gnum),is_grouped_site(gnum),grouped_site_done(maximum_nmr_site_num)
	logical :: selected_aa_used(10),current_grouped_used(maximum_nmr_site_num)
	logical :: assignment_feasible,structure_built,found_aa
	logical :: final_subset
	real :: ran2
	character(len=4) :: base_name,target_name,target_sequence(gnum)
	type(groupdetails) :: temp_group_1(repeated_unit,gnum)
	type(groupdetails) :: temp_group_2(repeated_unit,gnum)
	type(groupdetails) :: aa_group(repeated_unit,40)

	tgroup=group
	sub_cycle=0
	repeated_failed_target_count=0

	! Specific amino-acid limits which were not entered by the user remain 0..chain length.
	aa_min=0
	aa_max=amino_acid_residue_count
	do i=1,n_morethan
		do aa_index=1,20
			if (trim(aa_morethan(i)%gtype) == aa_names(aa_index)) then
				aa_min(aa_index)=aa_morethan(i)%min
				exit
			endif
		enddo
	enddo
	do i=1,n_restrictions
		do aa_index=1,20
			if (trim(aa_restrictions(i)%gtype) == aa_names(aa_index)) then
				aa_max(aa_index)=aa_restrictions(i)%max
				exit
			endif
		enddo
	enddo

	! Record the current sequence and separate free, grouped, and fixed sites.
	is_grouped_site=.false.
	do i=1,nmr_site_num
		is_grouped_site(nmr_site_ID(i))=.true.
	enddo
	is_free_site=.false.
	free_site_count=0
	do i=1,sitenum4mutation
		site=sitenum4mutation_group(i)
		if (is_grouped_site(site)) cycle
		free_site_count=free_site_count+1
		free_site_ID(free_site_count)=site
		is_free_site(site)=.true.
	enddo

	current_aa=0
	do site=1,gnum
		base_name=adjustl(group(1,site)%gtype)
		if (trim(base_name)=='ACE'.or.trim(base_name)=='NME'.or.trim(base_name)=='NHE') cycle
		if (len_trim(base_name) == 4) then
			if (base_name(1:1) == 'N' .or. base_name(1:1) == 'C') base_name=base_name(2:4)
		endif
		if (trim(base_name)=='HIS'.or.trim(base_name)=='HID'.or.trim(base_name) == 'HIP') base_name='HIE'
		found_aa=.false.
		do aa_index=1,20
			if (trim(base_name) == aa_names(aa_index)) then
				current_aa(site)=aa_index
				found_aa=.true.
				exit
			endif
		enddo
		if (.not.found_aa) then
			write(*,'(3A)') 'ERROR: unsupported residue name ',trim(base_name), &
				' while constructing the initial sequence.'
			stop
		endif
	enddo

	fixed_aa_count=0
	current_free_count=0
	do site=1,gnum
		if (current_aa(site) == 0 .or. is_grouped_site(site)) cycle
		if (is_free_site(site)) then
			current_free_count(current_aa(site))=current_free_count(current_aa(site))+1
		else
			fixed_aa_count(current_aa(site))=fixed_aa_count(current_aa(site))+1
		endif
	enddo
	if (sum(fixed_aa_count)+free_site_count+nmr_site_num /= amino_acid_residue_count) then
		write(*,*) 'ERROR: fixed, grouped, and free SCMF sites do not match the peptide amino-acid count.'
		stop
	endif

	do attempt=1,max_sequence_attempts
		! Examine every possible subset of the grouped pool and every legal
		! category composition. Retain the solution that changes the fewest sites.
		minimum_changes=huge(1)
		tie_count=0
		subset_position=0
		best_subset_position=0
		if (nmr_site_num > 0) then
			do i=1,nmr_site_num
				subset_position(i)=i
			enddo
		endif
		final_subset=.false.
		do
			grouped_aa_count=0
			selected_aa=0
			do i=1,nmr_site_num
				do aa_index=1,20
					if (trim(NMR_AA_pool(subset_position(i))) == aa_names(aa_index)) then
						selected_aa(i)=aa_index
						grouped_aa_count(aa_index)=grouped_aa_count(aa_index)+1
						exit
					endif
				enddo
			enddo

			grouped_preserved=0
			current_grouped_used=.false.
			do i=1,nmr_site_num
				do j=1,nmr_site_num
					if (current_grouped_used(j)) cycle
					if (current_aa(nmr_site_ID(j)) == selected_aa(i)) then
						current_grouped_used(j)=.true.
						grouped_preserved=grouped_preserved+1
						exit
					endif
				enddo
			enddo
			locked_aa_count=fixed_aa_count+grouped_aa_count
			if (.not.any(locked_aa_count > aa_max)) then
				locked_category_count=0
				do aa_index=1,20
					category=aa_category(aa_index)
					locked_category_count(category)=locked_category_count(category)+locked_aa_count(aa_index)
				enddo

				do h_count=fpho_min,fpho_max
					do p_count=fpol_min,fpol_max
						do c_count=fchg_min,fchg_max
							o_count=amino_acid_residue_count-h_count-p_count-c_count
							if (o_count < foth_min .or. o_count > foth_max) cycle
							trial_category_count=(/h_count,p_count,c_count,o_count/)
							free_category_count=trial_category_count-locked_category_count
							if (any(free_category_count < 0)) cycle
							if (sum(free_category_count) /= free_site_count) cycle

							assignment_feasible=.true.
							free_preserved=0
							do category=1,4
								lower_sum=0
								upper_sum=0
								preserve_gain=0
								remaining_needed=0
								do aa_index=1,20
									if (aa_category(aa_index) /= category) cycle
									aa_lower(aa_index)=max(0,aa_min(aa_index)-locked_aa_count(aa_index))
									aa_upper(aa_index)=aa_max(aa_index)-locked_aa_count(aa_index)
									if (aa_upper(aa_index) < aa_lower(aa_index)) assignment_feasible=.false.
									lower_sum=lower_sum+aa_lower(aa_index)
									upper_sum=upper_sum+max(0,aa_upper(aa_index))
									remaining_needed=remaining_needed+min(current_free_count(aa_index), &
										aa_lower(aa_index))
									preserve_gain=preserve_gain+max(0,min(current_free_count(aa_index), &
										aa_upper(aa_index))-aa_lower(aa_index))
								enddo
								if (free_category_count(category) < lower_sum .or. &
									free_category_count(category) > upper_sum) assignment_feasible=.false.
								free_preserved=free_preserved+remaining_needed+ &
									min(max(0,free_category_count(category)-lower_sum),preserve_gain)
							enddo
							if (.not.assignment_feasible) cycle

							candidate_changes=nmr_site_num-grouped_preserved+free_site_count-free_preserved
							if (candidate_changes < minimum_changes) then
								minimum_changes=candidate_changes
								tie_count=1
								target_category_count=trial_category_count
								best_subset_position=subset_position
							elseif (candidate_changes == minimum_changes) then
								tie_count=tie_count+1
								call ran_gen(ran2,0)
								if (ran2 < 1.0/real(tie_count)) then
									target_category_count=trial_category_count
									best_subset_position=subset_position
								endif
							endif
						enddo
					enddo
				enddo
			endif

			if (nmr_site_num == 0) exit
			i=nmr_site_num
			do while(i >= 1)
				if (subset_position(i) < NMR_pool_size-nmr_site_num+i) exit
				i=i-1
			enddo
			if (i == 0) then
				final_subset=.true.
			else
				subset_position(i)=subset_position(i)+1
				do j=i+1,nmr_site_num
					subset_position(j)=subset_position(j-1)+1
				enddo
			endif
			if (final_subset) exit
		enddo
		if (minimum_changes == huge(1)) cycle

		! Build a minimum-change target sequence. Grouped sites first keep any
		! current identity selected from their pool, then receive the unused names.
		target_sequence=group(1,:)%gtype
		selected_aa=0
		do i=1,nmr_site_num
			do aa_index=1,20
				if (trim(NMR_AA_pool(best_subset_position(i))) == aa_names(aa_index)) then
					selected_aa(i)=aa_index
					exit
				endif
			enddo
			grouped_site_order(i)=i
		enddo
		do i=nmr_site_num,2,-1
			call ran_gen(ran2,0)
			j=int(ran2*real(i))+1
			if (j > i) j=i
			swap_value=grouped_site_order(i)
			grouped_site_order(i)=grouped_site_order(j)
			grouped_site_order(j)=swap_value
		enddo
		selected_aa_used=.false.
		grouped_site_done=.false.
		grouped_target_aa=0
		do m=1,nmr_site_num
			i=grouped_site_order(m)
			do j=1,nmr_site_num
				if (.not.selected_aa_used(j) .and. current_aa(nmr_site_ID(i)) == selected_aa(j)) then
					grouped_target_aa(i)=selected_aa(j)
					grouped_site_done(i)=.true.
					selected_aa_used(j)=.true.
					exit
				endif
			enddo
		enddo
		do m=1,nmr_site_num
			i=grouped_site_order(m)
			if (grouped_site_done(i)) cycle
			eligible_count=0
			do j=1,nmr_site_num
				if (selected_aa_used(j)) cycle
				eligible_count=eligible_count+1
				eligible_aa(eligible_count)=j
			enddo
			call ran_gen(ran2,0)
			choice=int(ran2*real(eligible_count))+1
			if (choice > eligible_count) choice=eligible_count
			j=eligible_aa(choice)
			grouped_target_aa(i)=selected_aa(j)
			grouped_site_done(i)=.true.
			selected_aa_used(j)=.true.
		enddo
		do i=1,nmr_site_num
			site=nmr_site_ID(i)
			target_name=aa_names(grouped_target_aa(i))
			base_name=adjustl(group(1,site)%gtype)
			if (len_trim(base_name) == 4) then
				if (base_name(1:1)=='N'.or. base_name(1:1)=='C') target_name=base_name(1:1)//trim(target_name)
			endif
			target_sequence(site)=target_name
		enddo

		grouped_aa_count=0
		do i=1,nmr_site_num
			grouped_aa_count(grouped_target_aa(i))=grouped_aa_count(grouped_target_aa(i))+1
		enddo
		locked_aa_count=fixed_aa_count+grouped_aa_count
		locked_category_count=0
		do aa_index=1,20
			category=aa_category(aa_index)
			locked_category_count(category)=locked_category_count(category)+locked_aa_count(aa_index)
			aa_lower(aa_index)=max(0,aa_min(aa_index)-locked_aa_count(aa_index))
			aa_upper(aa_index)=aa_max(aa_index)-locked_aa_count(aa_index)
		enddo

		! Choose free-site amino-acid counts that retain as many current identities
		! as possible while satisfying all category and specific-AA bounds.
		target_free_count=aa_lower
		assignment_feasible=.true.
		do category=1,4
			remaining_needed=target_category_count(category)-locked_category_count(category)
			do aa_index=1,20
				if (aa_category(aa_index) == category) &
					remaining_needed=remaining_needed-target_free_count(aa_index)
			enddo

			do while(remaining_needed > 0)
				eligible_count=0
				do aa_index=1,20
					if (aa_category(aa_index) /= category) cycle
					if (target_free_count(aa_index) >= min(current_free_count(aa_index), &
						aa_upper(aa_index))) cycle
					eligible_count=eligible_count+1
					eligible_aa(eligible_count)=aa_index
				enddo
				if (eligible_count == 0) exit
				call ran_gen(ran2,0)
				choice=int(ran2*real(eligible_count))+1
				if (choice > eligible_count) choice=eligible_count
				aa_index=eligible_aa(choice)
				target_free_count(aa_index)=target_free_count(aa_index)+1
				remaining_needed=remaining_needed-1
			enddo
			do while(remaining_needed > 0)
				eligible_count=0
				do aa_index=1,20
					if (aa_category(aa_index) /= category) cycle
					if (target_free_count(aa_index) >= aa_upper(aa_index)) cycle
					eligible_count=eligible_count+1
					eligible_aa(eligible_count)=aa_index
				enddo
				if (eligible_count == 0) then
					assignment_feasible=.false.
					exit
				endif
				call ran_gen(ran2,0)
				choice=int(ran2*real(eligible_count))+1
				if (choice > eligible_count) choice=eligible_count
				aa_index=eligible_aa(choice)
				target_free_count(aa_index)=target_free_count(aa_index)+1
				remaining_needed=remaining_needed-1
			enddo
			if (.not.assignment_feasible) exit
		enddo
		if (.not.assignment_feasible) cycle

		! Retain current free-site identities first. Only the remaining sites are
		! assigned new names, randomly among equally minimal choices.
		remaining_aa_count=target_free_count
		do i=1,free_site_count
			free_site_order(i)=i
		enddo
		do i=free_site_count,2,-1
			call ran_gen(ran2,0)
			j=int(ran2*real(i))+1
			if (j > i) j=i
			swap_value=free_site_order(i)
			free_site_order(i)=free_site_order(j)
			free_site_order(j)=swap_value
		enddo
		do i=1,free_site_count
			site=free_site_ID(free_site_order(i))
			aa_index=current_aa(site)
			if (remaining_aa_count(aa_index) <= 0) cycle
			remaining_aa_count(aa_index)=remaining_aa_count(aa_index)-1
			free_site_order(i)=-free_site_order(i)
		enddo
		do i=1,free_site_count
			if (free_site_order(i) < 0) cycle
			site=free_site_ID(free_site_order(i))
			eligible_count=0
			do aa_index=1,20
				if (remaining_aa_count(aa_index) <= 0) cycle
				eligible_count=eligible_count+1
				eligible_aa(eligible_count)=aa_index
			enddo
			if (eligible_count == 0) then
				assignment_feasible=.false.
				exit
			endif
			call ran_gen(ran2,0)
			choice=int(ran2*real(eligible_count))+1
			if (choice > eligible_count) choice=eligible_count
			aa_index=eligible_aa(choice)
			remaining_aa_count(aa_index)=remaining_aa_count(aa_index)-1
			target_name=aa_names(aa_index)
			base_name=adjustl(group(1,site)%gtype)
			if (len_trim(base_name) == 4) then
				if (base_name(1:1) == 'N' .or. base_name(1:1) == 'C') &
					target_name=base_name(1:1)//trim(target_name)
			endif
			target_sequence(site)=target_name
		enddo
		if (.not.assignment_feasible .or. any(remaining_aa_count /= 0)) cycle

		! Verify the complete target before changing coordinates.
		trial_aa_count=0
		do site=1,gnum
			base_name=adjustl(target_sequence(site))
			if (trim(base_name) == 'ACE' .or. trim(base_name) == 'NME' .or. &
				trim(base_name) == 'NHE') cycle
			if (len_trim(base_name) == 4) then
				if (base_name(1:1) == 'N' .or. base_name(1:1) == 'C') base_name=base_name(2:4)
			endif
			do aa_index=1,20
				if (trim(base_name) == aa_names(aa_index)) then
					trial_aa_count(aa_index)=trial_aa_count(aa_index)+1
					exit
				endif
			enddo
		enddo
		trial_category_count=0
		do aa_index=1,20
			category=aa_category(aa_index)
			trial_category_count(category)=trial_category_count(category)+trial_aa_count(aa_index)
		enddo
		if (any(trial_aa_count < aa_min) .or. any(trial_aa_count > aa_max)) cycle
		if (trial_category_count(1) < fpho_min .or. trial_category_count(1) > fpho_max) cycle
		if (trial_category_count(2) < fpol_min .or. trial_category_count(2) > fpol_max) cycle
		if (trial_category_count(3) < fchg_min .or. trial_category_count(3) > fchg_max) cycle
		if (trial_category_count(4) < foth_min .or. trial_category_count(4) > foth_max) cycle

		! After an SCMF packing failure, prefer a new sequence with a different
		! amino acid at that failed chain position. Try up to ten otherwise-legal
		! proposals before allowing the same identity again.
		if (avoid_site >= 1 .and. avoid_site <= gnum) then
			if (trim(target_sequence(avoid_site)) == trim(avoid_aminoacid)) then
				repeated_failed_target_count=repeated_failed_target_count+1
				if (repeated_failed_target_count < max_avoid_attempts) cycle
				write(*,'(A,I0,A,I0,3A)') 'SCMF substitution could not find a different amino acid at site ', &
					avoid_site,' after ',max_avoid_attempts,' proposals; retaining ', &
					trim(target_sequence(avoid_site)),'.'
			endif
		endif

		! Build only the changed sites on chain 1, in random order.
		temp_group_1=group
		structure_built=.true.
		changed_count=0
		chainID=1
		do site=1,gnum
			if (trim(temp_group_1(chainID,site)%gtype) == trim(target_sequence(site))) cycle
			changed_count=changed_count+1
			changed_site(changed_count)=site
		enddo
		do i=changed_count,2,-1
			call ran_gen(ran2,0)
			j=int(ran2*real(i))+1
			if (j > i) j=i
			swap_value=changed_site(i)
			changed_site(i)=changed_site(j)
			changed_site(j)=swap_value
		enddo
		do i=1,changed_count
			site=changed_site(i)
			call findrotamer(site,temp_group_1,target_sequence(site),rotanum,aa_group,ip)
			feedback=0
			do m=1,rotanum
				call residue_replace(chainID,site,temp_group_1,m,aa_group,temp_group_2)
				call check_transplant(chainID,site,temp_group_2,feedback)
				if (feedback == 1) then
					temp_group_1=temp_group_2
					exit
				endif
			enddo
			if (feedback /= 1) then
				structure_built=.false.
				exit
			endif
		enddo

		if (structure_built) then
			sub_cycle=changed_count
			tgroup=temp_group_1
			write(*,'(A)',advance='no') 'Current sequence:  '
			do site=1,gnum
				write(*,'(1X,A)',advance='no') trim(group(1,site)%gtype)
			enddo
			write(*,*)
			write(*,'(A)',advance='no') 'Proposed sequence: '
			do site=1,gnum
				write(*,'(1X,A)',advance='no') trim(target_sequence(site))
			enddo
			write(*,*)
			write(*,'(A,I0)') 'Minimal number of amino acid mutation required: ',changed_count
			return
		endif
	enddo

	write(*,*) 'ERROR: unable to construct a legal minimum-change initial sequence after ', &
		max_sequence_attempts,' attempts.'
	stop
    end subroutine scmf_substitution

   ! subroutine find_nb_residues(ic, group, nbs)
   ! implicit none
   ! real			:: rxi, ryi, rzi, rxj, ryj, rzj
   ! 
   ! 
   ! 
   ! nres
   ! do i, nres
   !     do k, group(i,ic)%cnum1
			!if(group(i,ic)%atype1(k)=="CA") then
			!	rxi = group(ii,ic)%coo1(k,1)
			!	ryi = group(ii,ic)%coo1(k,2)
			!	rzi = group(ii,ic)%coo1(k,3)
   !         endif
   ! !if(group(ii,ic)%atype1(k)=="CA") then
   ! !    car(1)=group(ii,ic)%coo1(k,1)
   ! !    car(2)=group(ii,ic)%coo1(k,2)
   ! !    car(3)=group(ii,ic)%coo1(k,3)
   ! !endif
   ! !group
   ! 
   !     
   ! 
   ! 
   ! end subroutine find_nb_residues

	subroutine sidechain_optimization(stage,chainID,ic,group,group_para,S_numex,S_inb, &
		S_numex4,S_inb4,Dihedral4entropy)
	implicit none
	integer								:: grade, grade_num(6), monitor(6)
	integer								:: chainID, i, j, k, ic, account_num, flag, stage, trial_count
	integer								:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60, 60)
	integer								:: dihedral_num	
	integer, parameter					:: max_sidechain_grade=6
	real								:: delta_chi, cos_angle, sin_angle, error, t
	real								:: CA(3), rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3), Tmember(15,3)
	real(kind=16)                       :: h2_denominator, h3_denominator, Tenergy_min, Tenergy
	real								:: Dihedral4entropy(4)
	type(groupdetails)					:: group(repeated_unit,gnum), Tgroup(repeated_unit,gnum)
	type(energyparameters)				:: group_para(repeated_unit,gnum)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6), Tclass(6), class_old(6), class_new(6), class_min(6)
	type(dihedralparameters)			:: dihedral	
	
	real(kind=16)						:: energy_forward(max_sidechain_grade), energy_backward(max_sidechain_grade)
	real(kind=16)						:: gradient_old(max_sidechain_grade,1), gradient_new(max_sidechain_grade,1)
	real(kind=16)						:: Hessian_old(max_sidechain_grade,max_sidechain_grade), Hessian_new(max_sidechain_grade,max_sidechain_grade)
	real(kind=16)						:: H2(max_sidechain_grade,max_sidechain_grade), H3(max_sidechain_grade,max_sidechain_grade), H31(max_sidechain_grade,1)
	real(kind=16)						:: d(max_sidechain_grade,1), y(max_sidechain_grade,1), s(max_sidechain_grade,1), Tchi(max_sidechain_grade,1)

	call sidechain_category(chainID, ic, group, Iclass, grade, grade_num, index, monitor) ! find AA's side-chain info
	call dihedralangle_reading(group(chainID,ic)%gtype, dihedral_num, dihedral)			  ! find AA's dihedral angle info

	energy_forward=0.0_16
	energy_backward=0.0_16
	gradient_old=0.0_16
	gradient_new=0.0_16
	Hessian_old=0.0_16
	Hessian_new=0.0_16
	H2=0.0_16
	H3=0.0_16
	H31=0.0_16
	d=0.0_16
	y=0.0_16
	s=0.0_16
	Tchi=0.0_16

	Tgroup=group
	do i=1, Tgroup(chainID,ic)%cnum1
		if(Tgroup(chainID,ic)%atype1(i)=="CA") then										! find coord for CA
			CA(1)=Tgroup(chainID,ic)%coo1(i,1); CA(2)=Tgroup(chainID,ic)%coo1(i,2); CA(3)=Tgroup(chainID,ic)%coo1(i,3)
		endif
	enddo
	s=0.0_16
	class_new=Iclass

30	continue
	if(stage==0) then               ! Use a 5-degree step in stage 0 and a 1-degree step in stage 1.
		delta_chi=5
	elseif(stage==1) then
		delta_chi=1
	endif
	cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)	! Calculate the sine and cosine of the trial angle.
	do i=1, grade											! loop all grade (number of heavy atom dihedral axis)
		Tclass=class_new
		if(i==1) then										! The first rotation axis is CA-CB.
			rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)	! Calculate the CA-CB axis vector.
			rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
			rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2) ! Normalize the CA-CB axis vector.
			rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		else												! Calculate each subsequent rotation axis.
			rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
			rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		endif

		call axisrotation(rotaxis, cos_angle, sin_angle, m)
		
		do j=(i+1), (grade+1) ! Translate downstream atoms relative to rotation axis i.
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
				Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
				Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
			enddo
			
			Tmember=matmul(Tclass(j)%member, m)
			Tclass(j)%member=Tmember
			
			do k=1, grade_num(j) ! Return to the original coordinate frame and retain three decimals.
				Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
				Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000
				Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
			enddo
		enddo
		
		do j=1, Tgroup(chainID,ic)%cnum2 ! Store the rotated side-chain coordinates in Tgroup.
			Tgroup(chainID,ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
			Tgroup(chainID,ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
			Tgroup(chainID,ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
		enddo

		call sidechain_energy(stage, chainID, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_forward(i))
	enddo
	
	cos_angle=cosd(-delta_chi); sin_angle=sind(-delta_chi)
	do i=1, grade
		Tclass=class_new
		if(i==1) then
			rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
			rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		else
			rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
			rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		endif

		call axisrotation(rotaxis, cos_angle, sin_angle, m)
		
		do j=(i+1), (grade+1)
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
				Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
				Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
			enddo
			
			Tmember=matmul(Tclass(j)%member, m)
			Tclass(j)%member=Tmember
			
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
				Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000
				Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
			enddo
		enddo
		
		do j=1, Tgroup(chainID,ic)%cnum2
			Tgroup(chainID,ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
			Tgroup(chainID,ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
			Tgroup(chainID,ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
		enddo

		call sidechain_energy(stage, chainID, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_backward(i))
	enddo
	
	do i=1, grade
		gradient_new(i,1)=(energy_forward(i)-energy_backward(i))/(2*delta_chi)	
	enddo

	do i=1, grade
		do j=1, grade
			if(i==j) then
				Hessian_new(i,j)=1
			else
				Hessian_new(i,j)=0
			endif
		enddo
	enddo

	account_num=0
	t=0.0
	do while(.true.)
		Hessian_old=Hessian_new
		gradient_old=gradient_new
		d=-matmul(Hessian_old, gradient_old)
		
		trial_count=0
		do while(.true.)
			class_old=Iclass
			Tchi=t*d
			flag=0
			if(stage==0) then
				do i=1, grade
					if(Tchi(i,1).gt.5.0) then
						Tchi(i,1)=5.0
						flag=1
					elseif(Tchi(i,1).lt.(-5.0)) then
						Tchi(i,1)=-5.0
						flag=1
					endif
				enddo
			elseif(stage==1) then
				do i=1, grade
					if(Tchi(i,1).gt.1.0) then
						Tchi(i,1)=1.0
						flag=1
					elseif(Tchi(i,1).lt.(-1.0)) then
						Tchi(i,1)=-1.0
						flag=1
					endif
				enddo
			endif
			if(flag==1) account_num=account_num+1
			
			s=s+Tchi
			do i=1, grade
				cos_angle=cosd(s(i,1)); sin_angle=sind(s(i,1))
				if(i==1) then
					rotaxis_x=class_old(i)%member(monitor(i),1)-CA(1)
					rotaxis_y=class_old(i)%member(monitor(i),2)-CA(2)
					rotaxis_z=class_old(i)%member(monitor(i),3)-CA(3)
					rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				else
					rotaxis_x=class_old(i)%member(monitor(i),1)-class_old(i-1)%member(monitor(i-1),1)
					rotaxis_y=class_old(i)%member(monitor(i),2)-class_old(i-1)%member(monitor(i-1),2)
					rotaxis_z=class_old(i)%member(monitor(i),3)-class_old(i-1)%member(monitor(i-1),3)
					rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
					rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				endif

				call axisrotation(rotaxis, cos_angle, sin_angle, m)
					
				do j=(i+1), (grade+1)
					do k=1, grade_num(j)
						class_old(j)%member(k,1)=class_old(j)%member(k,1)-class_old(i)%member(monitor(i),1)
						class_old(j)%member(k,2)=class_old(j)%member(k,2)-class_old(i)%member(monitor(i),2)
						class_old(j)%member(k,3)=class_old(j)%member(k,3)-class_old(i)%member(monitor(i),3)
					enddo
						
					Tmember=matmul(class_old(j)%member, m)
					class_old(j)%member=Tmember
						
					do k=1, grade_num(j)
						class_old(j)%member(k,1)=anint((class_old(j)%member(k,1)+class_old(i)%member(monitor(i),1))*1000)/1000
						class_old(j)%member(k,2)=anint((class_old(j)%member(k,2)+class_old(i)%member(monitor(i),2))*1000)/1000				
						class_old(j)%member(k,3)=anint((class_old(j)%member(k,3)+class_old(i)%member(monitor(i),3))*1000)/1000
					enddo
				enddo
			enddo
				
			do j=1, Tgroup(chainID,ic)%cnum2
				Tgroup(chainID,ic)%coo2(j,1)=class_old(index(j)%class_No)%member(index(j)%member_No,1)
				Tgroup(chainID,ic)%coo2(j,2)=class_old(index(j)%class_No)%member(index(j)%member_No,2)
				Tgroup(chainID,ic)%coo2(j,3)=class_old(index(j)%class_No)%member(index(j)%member_No,3)
			enddo

			call sidechain_energy(stage, chainID, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, Tenergy)
			
			if(t==0.0) then
				Tenergy_min=Tenergy
				class_min=class_old
			else
				if(Tenergy.lt.Tenergy_min) then
					Tenergy_min=Tenergy
					class_min=class_old
					trial_count=trial_count+1
				else
					s=s-Tchi
					goto 10
				endif
			endif				
			t=delta_chi
		enddo
10		continue
		if(stage==0) then
			if(account_num.gt.20) goto 20
		elseif(stage==1) then
			if(account_num.gt.40) goto 20
		endif

		if(trial_count==0) goto 20
		error=0.0
		do i=1, grade
			error=error+abs(d(i,1))
		enddo
		if(error.lt.0.01) goto 20

		cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)
		do i=1, grade
			Tclass=class_min
			if(i==1) then
				rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
				rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			else
				rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
				rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			endif

			call axisrotation(rotaxis, cos_angle, sin_angle, m)
			
			do j=(i+1), (grade+1)
				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
					Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
					Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
				enddo
				
				Tmember=matmul(Tclass(j)%member, m)
				Tclass(j)%member=Tmember
				
				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
					Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000				
					Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
				enddo
			enddo
			
			do j=1, Tgroup(chainID,ic)%cnum2
				Tgroup(chainID,ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
				Tgroup(chainID,ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
				Tgroup(chainID,ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
			enddo

			call sidechain_energy(stage, chainID, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_forward(i))
		enddo
		
		cos_angle=cosd(-delta_chi); sin_angle=sind(-delta_chi)
		do i=1, grade
			Tclass=class_min
			if(i==1) then
				rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
				rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			else
				rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
				rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
				rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
				rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
				rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			endif

			call axisrotation(rotaxis, cos_angle, sin_angle, m)
			
			do j=(i+1), (grade+1)
				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
					Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
					Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
				enddo
				
				Tmember=matmul(Tclass(j)%member, m)
				Tclass(j)%member=Tmember
				
				do k=1, grade_num(j)
					Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
					Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000				
					Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
				enddo
			enddo
			
			do j=1, Tgroup(chainID,ic)%cnum2
				Tgroup(chainID,ic)%coo2(j,1)=Tclass(index(j)%class_No)%member(index(j)%member_No,1)
				Tgroup(chainID,ic)%coo2(j,2)=Tclass(index(j)%class_No)%member(index(j)%member_No,2)
				Tgroup(chainID,ic)%coo2(j,3)=Tclass(index(j)%class_No)%member(index(j)%member_No,3)
			enddo

			call sidechain_energy(stage, chainID, ic, Tgroup, group_para, S_numex, S_inb, S_numex4, S_inb4, dihedral_num, dihedral, energy_backward(i))
		enddo
		
		do i=1, grade
			gradient_new(i,1)=(energy_forward(i)-energy_backward(i))/(2*delta_chi)			
		enddo
				
		y=gradient_new-gradient_old

		H2=matmul(s, transpose(s))
		h2_denominator=0.0
		do i=1, grade
			h2_denominator=h2_denominator+s(i,1)*y(i,1)
		enddo

		H3=matmul(Hessian_old, matmul(y, matmul(transpose(y), Hessian_old)))
		H31=matmul(Hessian_old, y)
		h3_denominator=0.0			
		do i=1, grade
			h3_denominator=h3_denominator+y(i,1)*H31(i,1)
		enddo			
		
		Hessian_new=Hessian_old+H2/h2_denominator-H3/h3_denominator
 
	enddo
20	continue
	class_new=class_min
	
	if(stage==0.and.Tenergy_min.lt.100.0) then
		stage=1
		goto 30
	endif

	Iclass=class_min

	do i=1, group(chainID,ic)%cnum2
		group(chainID,ic)%coo2(i,1)=Iclass(index(i)%class_No)%member(index(i)%member_No,1)
		group(chainID,ic)%coo2(i,2)=Iclass(index(i)%class_No)%member(index(i)%member_No,2)
		group(chainID,ic)%coo2(i,3)=Iclass(index(i)%class_No)%member(index(i)%member_No,3)
	enddo

	if(stage==1) then
		call torsionangle4sidechain(group, chainID, ic, grade, dihedral, Dihedral4entropy)
	endif

	return
	end	subroutine sidechain_optimization

end module advancedfunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module optimization_techniques

	use constant
	use pdbfile
	use mathfunction	
	use database
	use transplant
	use energy_calculation
	use advancedfunction	

	type(groupdetails), allocatable, save :: opt_temp_group(:,:),opt_aa_group(:,:)
	type(groupdetails), allocatable, save :: opt_best_group(:,:)
	type(energyparameters), allocatable, save		:: opt_group_para(:,:)
	integer, allocatable, save						:: opt_numex(:), opt_numex4(:), opt_inb(:,:), opt_inb4(:,:)

	contains
	subroutine ensure_optimization_workspace()
	implicit none

	if(.not.allocated(opt_temp_group)) allocate(opt_temp_group(repeated_unit,gnum))
	if(.not.allocated(opt_aa_group)) allocate(opt_aa_group(repeated_unit,40))
	if(.not.allocated(opt_best_group)) allocate(opt_best_group(repeated_unit,gnum))
	if(.not.allocated(opt_group_para)) allocate(opt_group_para(repeated_unit,gnum))
	if(.not.allocated(opt_numex)) allocate(opt_numex(repeated_unit*atom_num))
	if(.not.allocated(opt_numex4)) allocate(opt_numex4(repeated_unit*atom_num))
	if(.not.allocated(opt_inb)) allocate(opt_inb(repeated_unit*atom_num,20))
	if(.not.allocated(opt_inb4)) allocate(opt_inb4(repeated_unit*atom_num,60))

	return
	end subroutine ensure_optimization_workspace

	subroutine replace_single_site_aa(group)
	implicit none
	integer							:: i, ic, feedback
	character*4						:: aminoacid_name, cur_aminoacid_name
	type(groupdetails)				:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum)

	do i=1, void_site_num
		if(len_trim(void_site_fixed_name(i)) == 0) cycle
		ic = void_site(i)
		aminoacid_name = void_site_fixed_name(i)
		cur_aminoacid_name = group(1,ic)%gtype
		if((cur_aminoacid_name==aminoacid_name).or.(cur_aminoacid_name=="ACE").or.(cur_aminoacid_name=="NME")) cycle

		call sequence_mutation_nonthermal(ic, group, aminoacid_name, temp_group, feedback)
		if(feedback == 1) then
			group = temp_group
		else
			open(20, file="error.txt", access="append")
				write(20,*) "Unable to apply single-site constraint at sequence site ", ic, " to ", aminoacid_name
			close(20)
			stop "ERROR: unable to apply single-site constraint at site "
		endif
	enddo

	original_group = group

	return
	end subroutine replace_single_site_aa

	subroutine MC_technique_sequence(score_old, energy_old, binding_vdw_old, binding_ele_old, binding_sgb_old, &
		binding_snp_old, entropy_old, score4hydration_old, Pagg_old, group, entropy4individual, score_new, &
		energy_new, binding_vdw_new, binding_ele_new, binding_sgb_new, binding_snp_new, entropy_new, &
		score4hydration_new, Pagg_new, tgroup, Tentropy4individual, feedback)
	implicit none
	integer							:: feedback
	real							:: ran2
	real(kind=8)                    :: energy_old, energy_new, score_old, score_new, score_change
	real(kind=8)					:: binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old
	real(kind=8)					:: binding_vdw_new, binding_ele_new, binding_sgb_new, binding_snp_new
	real							:: entropy_old, entropy_new
	real							:: score4hydration_old, Pagg_old
	real							:: score4hydration_new, Pagg_new
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum)
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)

	feedback=0
	score_change=score_new-score_old
	if(score_change.le.0) then
		group=tgroup
		if (ESURF_flag == 2 .and. sasa_trial%valid) sasa_current=sasa_trial
		score_old=score_new
		energy_old=energy_new
		binding_vdw_old=binding_vdw_new
		binding_ele_old=binding_ele_new
		binding_sgb_old=binding_sgb_new
		binding_snp_old=binding_snp_new
		entropy_old=entropy_new
		score4hydration_old=score4hydration_new
		Pagg_old=Pagg_new
		entropy4individual=Tentropy4individual
		feedback=1
	else
		call ran_gen(ran2,0)
		if(ran2.le.exp(-score_change/ekt_seq)) then
			group=tgroup
			if (ESURF_flag == 2 .and. sasa_trial%valid) sasa_current=sasa_trial
			score_old=score_new
			energy_old=energy_new
			binding_vdw_old=binding_vdw_new
			binding_ele_old=binding_ele_new
			binding_sgb_old=binding_sgb_new
			binding_snp_old=binding_snp_new
			entropy_old=entropy_new
			score4hydration_old=score4hydration_new
			Pagg_old=Pagg_new
			entropy4individual=Tentropy4individual
			feedback=1
		endif
	endif
	
	return
	end subroutine MC_technique_sequence
	
	subroutine MC_technique_sheet(score_old, energy_old, binding_vdw_old, binding_ele_old, binding_sgb_old, &
		binding_snp_old, entropy_old, score4hydration_old, Pagg_old, group, entropy4individual, score_new, &
		energy_new, binding_vdw_new, binding_ele_new, binding_sgb_new, binding_snp_new, entropy_new, &
		score4hydration_new, Pagg_new, tgroup, Tentropy4individual, feedback)
	implicit none
	integer							:: feedback
    real                            :: ran2
	real(kind=8)                    :: score_old, score_new, score_change
	real(kind=8)                    :: energy_old, energy_new
	real(kind=8)					:: binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old
	real(kind=8)					:: binding_vdw_new, binding_ele_new, binding_sgb_new, binding_snp_new
	real							:: entropy_old, entropy_new
	real							:: score4hydration_old, Pagg_old
	real							:: score4hydration_new, Pagg_new
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum)
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)
	
	feedback=0
	score_change=score_new-score_old
	if(score_change.le.0) then
		group=tgroup
		if (ESURF_flag == 2 .and. sasa_trial%valid) sasa_current=sasa_trial
		score_old=score_new
		energy_old=energy_new
		binding_vdw_old=binding_vdw_new
		binding_ele_old=binding_ele_new
		binding_sgb_old=binding_sgb_new
		binding_snp_old=binding_snp_new
		entropy_old=entropy_new
		score4hydration_old=score4hydration_new
		Pagg_old=Pagg_new
		entropy4individual=Tentropy4individual
		feedback=1
	else
		call ran_gen(ran2,0)
		if(ran2.le.exp(-score_change/ekt_sheet)) then
			group=tgroup
			if (ESURF_flag == 2 .and. sasa_trial%valid) sasa_current=sasa_trial
			score_old=score_new
			energy_old=energy_new
			binding_vdw_old=binding_vdw_new
			binding_ele_old=binding_ele_new
			binding_sgb_old=binding_sgb_new
			binding_snp_old=binding_snp_new
			entropy_old=entropy_new
			score4hydration_old=score4hydration_new
			Pagg_old=Pagg_new
			entropy4individual=Tentropy4individual
			feedback=1
		endif
	endif
	
	return
	end subroutine MC_technique_sheet


	subroutine sequence_mutation_nonthermal(ic, group, aminoacid_name, tgroup, flag)
	implicit none
	integer							:: flag, i, j, m, m_best, rotanum, feedback
	integer							:: chainID, ic, ic_1, ic_2, stage, ip
	integer							:: W_numex(repeated_unit*atom_num), W_inb(repeated_unit*atom_num,20), W_numex4(repeated_unit*atom_num), W_inb4(repeated_unit*atom_num,60)
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	real							:: ran2
    real(kind=8)                    :: score, binding_energy, vdw_energy
	real(kind=8)					:: binding_vdw_tmp, binding_ele_tmp, binding_sgb_tmp, binding_snp_tmp
	real							:: vdw_energy_min, score_min
	real							:: entropy, score4hydration, Pagg
	real							:: Dihedral4entropy(4)
	character*4						:: aminoacid_name

	type(groupdetails)				:: aa_backup
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)

	call ensure_optimization_workspace()
	call findrotamer(ic, group, aminoacid_name, rotanum, opt_aa_group, ip)
	
	flag=0
	if(aminoacid_name=="GLY".or.aminoacid_name=="ALA".or.aminoacid_name=="PRO".or. &
	   aminoacid_name=="NGLY".or.aminoacid_name=="NALA".or.aminoacid_name=="NPRO".or. &
	   aminoacid_name=="CGLY".or.aminoacid_name=="CALA".or.aminoacid_name=="CPRO".or. &
	   aminoacid_name=="NME".or.aminoacid_name=="NHE".or.aminoacid_name=="ACE") then	

		tgroup=group
		do chainID=1, repeated_unit
				
			m_best=0
			vdw_energy_min=500.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, opt_aa_group, opt_temp_group)

				if(m==1) then
					call energy_parameter(opt_temp_group, opt_group_para)
					call atom_links4sidechain(chainID, ic, opt_temp_group, S_numex, S_inb, S_numex4, S_inb4)
				endif

				call check_transplant(chainID, ic, opt_temp_group, feedback)

				if(feedback==1) then
					call vdwenergy(chainID, ic, opt_temp_group, opt_group_para, vdw_energy)
					if(vdw_energy.lt.vdw_energy_min) then
						vdw_energy_min=vdw_energy
						call backup4sidechain(0, chainID, ic, opt_temp_group, aa_backup)
						m_best=m
					endif
				endif
			enddo

			if(m_best==0) goto 10
			call backup4sidechain(1, chainID, ic, tgroup, aa_backup)
		enddo

		flag=1
10		continue

	else
		tgroup=group
        vdw_energy_min=1000000000.0*gnum
        
		do chainID=1, repeated_unit
            
			m_best=0            
			score_min=500.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, opt_aa_group, opt_temp_group)
				if(m==1) then
					call energy_parameter(opt_temp_group, opt_group_para)
					call atom_links4sidechain(chainID, ic, opt_temp_group, S_numex, S_inb, S_numex4, S_inb4)
					call atom_links(opt_temp_group, W_numex, W_inb, W_numex4, W_inb4)
				endif
					
				stage=0
				call sidechain_optimization(stage, chainID, ic, opt_temp_group, opt_group_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)

				if(stage==1) then
					call bindingenergy_noentropy(opt_temp_group, opt_group_para, W_numex, W_inb, W_numex4, &
							W_inb4, score, binding_energy, binding_vdw_tmp, binding_ele_tmp, &
							binding_sgb_tmp, binding_snp_tmp, vdw_energy, score4hydration, Pagg)
					if((score.lt.score_min) .and. (vdw_energy.lt.vdw_energy_min)) then
						score_min=score
						call backup4sidechain(0, chainID, ic, opt_temp_group, aa_backup)
						m_best=m
					endif
				endif
			enddo
				
			if(m_best==0) goto 20
			call backup4sidechain(1, chainID, ic, tgroup, aa_backup)
		enddo

		flag=1
20		continue
	endif
	
	return
	end subroutine sequence_mutation_nonthermal
	
	
	subroutine sequence_mutation(ic, group, entropy4individual, aminoacid_name, tgroup, Tentropy4individual, flag)
	implicit none	
	integer							:: flag, i, j, m, m_best, feedback
	integer							:: categoryID, chainID, ic, ic_1, ic_2, ip, stage
	integer							:: W_numex(repeated_unit*atom_num), W_inb(repeated_unit*atom_num,20), W_numex4(repeated_unit*atom_num), W_inb4(repeated_unit*atom_num,60)
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	integer							:: rotanum, obs, grade
	real							:: ran2
	real							:: vdw_energy_min, score_min
    real(kind=8)					:: score, binding_energy, vdw_energy
	real(kind=8)					:: binding_vdw_tmp, binding_ele_tmp, binding_sgb_tmp, binding_snp_tmp
	real							:: entropy, score4hydration, Pagg
	real							:: Dihedral4entropy(4), matrix(34,4)
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum)
	character*4						:: aminoacid_name

	type(groupdetails)				:: aa_backup
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)

	call ensure_optimization_workspace()
	call findrotamer(ic, group, aminoacid_name, rotanum, opt_aa_group, ip)
	grade=aa_lib(ip)%grade
	
	flag=0
	if(aminoacid_name=="GLY".or.aminoacid_name=="ALA".or.aminoacid_name=="PRO".or. &
	   aminoacid_name=="NGLY".or.aminoacid_name=="NALA".or.aminoacid_name=="NPRO".or. &
	   aminoacid_name=="CGLY".or.aminoacid_name=="CALA".or.aminoacid_name=="CPRO".or. &
	   aminoacid_name=="NME".or.aminoacid_name=="NHE".or.aminoacid_name=="ACE") then	

		tgroup=group
		Tentropy4individual=entropy4individual
		do chainID=1, repeated_unit
				
			m_best=0
			vdw_energy_min=500.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, opt_aa_group, opt_temp_group)

				if(m==1) then
					call energy_parameter(opt_temp_group, opt_group_para)
					call atom_links4sidechain(chainID, ic, opt_temp_group, S_numex, S_inb, S_numex4, S_inb4)
				endif

				call check_transplant(chainID, ic, opt_temp_group, feedback)

				if(feedback==1) then
					call vdwenergy(chainID, ic, opt_temp_group, opt_group_para, vdw_energy)
					if(vdw_energy.lt.vdw_energy_min) then
						vdw_energy_min=vdw_energy
						call backup4sidechain(0, chainID, ic, opt_temp_group, aa_backup)
						m_best=m
					endif
				endif
			enddo

			if(m_best==0) goto 10
			matrix=0.0; obs=4
			call entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)
			Tentropy4individual(chainID,ic)=entropy
			call backup4sidechain(1, chainID, ic, tgroup, aa_backup)
		enddo
		
		flag=1
10		continue
											
	else
		tgroup=group
		Tentropy4individual=entropy4individual
        vdw_energy_min=1000000000.0*gnum
		do chainID=1, repeated_unit

			m_best=0
			score_min=500.0
			obs=0
			matrix=0.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, opt_aa_group, opt_temp_group)
				if(m==1) then
					call energy_parameter(opt_temp_group, opt_group_para)
					call atom_links4sidechain(chainID, ic, opt_temp_group, S_numex, S_inb, S_numex4, S_inb4)
					call atom_links(opt_temp_group, W_numex, W_inb, W_numex4, W_inb4)
				endif
					
				stage=0
				call sidechain_optimization(stage, chainID, ic, opt_temp_group, opt_group_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)
						
				if(stage==1) then
					obs=obs+1
					do j=1, grade
						matrix(obs,j)=Dihedral4entropy(j) ! why not pick the matrix with low
					enddo
						
					call bindingenergy_noentropy(opt_temp_group, opt_group_para, W_numex, W_inb, W_numex4, &
							W_inb4, score, binding_energy, binding_vdw_tmp, binding_ele_tmp, &
							binding_sgb_tmp, binding_snp_tmp, vdw_energy, score4hydration, Pagg)
					if((score.lt.score_min) .and. (vdw_energy.lt.vdw_energy_min)) then !.and. vdw.lt.vdw_min
						score_min=score
						call backup4sidechain(0, chainID, ic, opt_temp_group, aa_backup)
						m_best=m
					endif
				endif
			enddo
										
			if(m_best==0) goto 20			
			call entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)
			Tentropy4individual(chainID,ic)=entropy
			call backup4sidechain(1, chainID, ic, tgroup, aa_backup)
		enddo
		
		flag=1
20		continue
	endif

	return
	end subroutine sequence_mutation

	subroutine sidechain_optimization4scmf(chainID,ic,group,rotanum,aa_group, &
		best_group,best_stage,stage1_count,best_score,best_local_vdw,prefer_local_vdw)
	implicit none
	integer, intent(in) :: chainID,ic,rotanum
	integer, intent(out) :: best_stage,stage1_count
	real(kind=8), intent(out) :: best_score,best_local_vdw
	type(groupdetails) :: group(repeated_unit,gnum)
	type(groupdetails) :: aa_group(repeated_unit,40)
	type(groupdetails), intent(out) :: best_group(repeated_unit,gnum)
	integer :: m,stage
	integer :: W_numex(repeated_unit*atom_num),W_numex4(repeated_unit*atom_num)
	integer :: W_inb(repeated_unit*atom_num,20),W_inb4(repeated_unit*atom_num,60)
	integer :: S_numex(60),S_numex4(60),S_inb(60,20),S_inb4(60,60)
	real :: Dihedral4entropy(4),score4hydration,Pagg
	real(kind=8) :: score,binding_energy,local_vdw,complex_vdw
	real(kind=8) :: binding_vdw,binding_ele,binding_sgb,binding_snp
	logical, intent(in), optional :: prefer_local_vdw
	logical :: candidate_saved,use_local_vdw,better_candidate

	call ensure_optimization_workspace()
	best_group=group
	best_stage=0
	stage1_count=0
	best_score=huge(1.0d0)
	best_local_vdw=huge(1.0d0)
	candidate_saved=.false.
	use_local_vdw=.false.
	if (present(prefer_local_vdw)) use_local_vdw=prefer_local_vdw

	do m=1,rotanum
		call residue_replace(chainID,ic,group,m,aa_group,opt_temp_group)
		if (m == 1) then
			call energy_parameter(opt_temp_group,opt_group_para)
			call atom_links4sidechain(chainID,ic,opt_temp_group,S_numex,S_inb,S_numex4,S_inb4)
			call atom_links(opt_temp_group,W_numex,W_inb,W_numex4,W_inb4)
		endif

		stage=0
		call sidechain_optimization(stage,chainID,ic,opt_temp_group,opt_group_para, &
			S_numex,S_inb,S_numex4,S_inb4,Dihedral4entropy)
		call vdwenergy(chainID,ic,opt_temp_group,opt_group_para,local_vdw)

		if (stage == 1) then
			stage1_count=stage1_count+1
			call bindingenergy_noentropy(opt_temp_group,opt_group_para,W_numex,W_inb, &
				W_numex4,W_inb4,score,binding_energy,binding_vdw,binding_ele, &
				binding_sgb,binding_snp,complex_vdw,score4hydration,Pagg)
			better_candidate=best_stage == 0
			if (best_stage == 1) then
				if (use_local_vdw) then
					better_candidate=local_vdw < best_local_vdw
				else
					better_candidate=score < best_score
				endif
			endif
			if (better_candidate) then
				best_stage=1
				best_score=score
				best_local_vdw=local_vdw
				opt_best_group=opt_temp_group
				candidate_saved=.true.
			endif
		elseif (best_stage == 0 .and. local_vdw < best_local_vdw) then
			best_local_vdw=local_vdw
			opt_best_group=opt_temp_group
			candidate_saved=.true.
		endif
	enddo

	if (candidate_saved) best_group=opt_best_group
	end subroutine sidechain_optimization4scmf
	
	subroutine sequence_mutation_4scmfloop(chainID,ic,group,aminoacid_name,tgroup,flag, &
		rotamer_count,transplant_count,stage1_count,best_attempt_score,best_attempt_vdw, &
		relaxed_energy,prefer_local_vdw)
	implicit none
	integer :: flag,m,m_best,rotanum,feedback,chainID,ic,ip,best_stage
	integer, intent(out) :: rotamer_count,transplant_count,stage1_count
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	real(kind=8) :: vdw_energy,vdw_energy_min
	real(kind=8), intent(out) :: best_attempt_score,best_attempt_vdw
	logical, intent(in), optional :: relaxed_energy
	logical, intent(in), optional :: prefer_local_vdw
	logical :: allow_high_energy,use_local_vdw
	character*4						:: aminoacid_name

	type(groupdetails)				:: aa_backup
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)

	call ensure_optimization_workspace()
	call findrotamer(ic, group, aminoacid_name, rotanum, opt_aa_group, ip)
	allow_high_energy=.false.
	if (present(relaxed_energy)) allow_high_energy=relaxed_energy
	use_local_vdw=.false.
	if (present(prefer_local_vdw)) use_local_vdw=prefer_local_vdw
	rotamer_count=rotanum
	transplant_count=0
	stage1_count=0
	best_attempt_score=huge(1.0d0)
	best_attempt_vdw=huge(1.0d0)
	
	flag=0
	if(aminoacid_name=="GLY".or.aminoacid_name=="ALA".or.aminoacid_name=="PRO".or. &
	   aminoacid_name=="NGLY".or.aminoacid_name=="NALA".or.aminoacid_name=="NPRO".or. &
	   aminoacid_name=="CGLY".or.aminoacid_name=="CALA".or.aminoacid_name=="CPRO".or. &
	   aminoacid_name=="NME".or.aminoacid_name=="NHE".or.aminoacid_name=="ACE") then	

		tgroup=group
		m_best=0
		vdw_energy_min=500.0
		if (allow_high_energy) vdw_energy_min=huge(1.0d0)
		do m=1, rotanum
			call residue_replace(chainID, ic, tgroup, m, opt_aa_group, opt_temp_group)
			if (m == 1) opt_best_group=opt_temp_group

			if(m==1) then
				call energy_parameter(opt_temp_group, opt_group_para)
				call atom_links4sidechain(chainID, ic, opt_temp_group, S_numex, S_inb, S_numex4, S_inb4)
			endif
			call check_transplant(chainID, ic, opt_temp_group, feedback)

			if(feedback==1) then
				transplant_count=transplant_count+1
				call vdwenergy(chainID, ic, opt_temp_group, opt_group_para, vdw_energy)
				if (vdw_energy < best_attempt_vdw) then
					best_attempt_vdw=vdw_energy
					opt_best_group=opt_temp_group
				endif
				if(vdw_energy.lt.vdw_energy_min) then
					vdw_energy_min=vdw_energy
					call backup4sidechain(0, chainID, ic, opt_temp_group, aa_backup)
					m_best=m
				endif
			endif
		enddo

		if(m_best==0) then
			if (rotanum > 0) tgroup=opt_best_group
			goto 10
		endif
		call backup4sidechain(1, chainID, ic, tgroup, aa_backup)

		flag=1
10		continue

	else
		call sidechain_optimization4scmf(chainID,ic,group,rotanum,opt_aa_group, &
			tgroup,best_stage,stage1_count,best_attempt_score,best_attempt_vdw,use_local_vdw)
		if (best_stage == 1) then
			if (allow_high_energy .or. best_attempt_score < 500.0d0) flag=1
		endif
	endif
	
	return
	
	end subroutine sequence_mutation_4scmfloop

	subroutine find_cb_neighbors(group,center_chain,center_site,nb_count, &
		nb_chain,nb_site,nb_dist)
	implicit none
	real, parameter :: ca_cb_neighbor_cutoff=9.0
	type(groupdetails), intent(in) :: group(repeated_unit,gnum)
	integer, intent(in) :: center_chain,center_site
	integer, intent(out) :: nb_count
	integer, intent(out) :: nb_chain(repeated_unit*gnum)
	integer, intent(out) :: nb_site(repeated_unit*gnum)
	real, intent(out) :: nb_dist(repeated_unit*gnum)
	integer :: chainID,site,atom
	real :: center_ca(3),other_cb(3),distance
	logical :: center_found,other_found
	character(len=4) :: residue_name

	nb_count=0
	nb_chain=0
	nb_site=0
	nb_dist=0.0
	call get_backbone_atom(group,center_chain,center_site,'CA',center_ca,center_found)
	if (.not.center_found) then
		write(*,'(A,I0,A,I0,A)') 'WARNING: cannot find CA for chain ',center_chain,', site ',center_site,'. Neighbor repacking skipped.'
		return
	endif

	do chainID=1,repeated_unit
		do site=1,gnum
			if (chainID == center_chain .and. site == center_site) cycle
			residue_name=adjustl(group(chainID,site)%gtype)
			if (trim(residue_name) == 'ACE' .or. trim(residue_name) == 'NME' .or. &
				trim(residue_name) == 'NHE') cycle
			other_found=.false.
			other_cb=0.0
			do atom=1,group(chainID,site)%cnum2
				if (group(chainID,site)%atype2(atom) == 'CB') then
					other_cb=group(chainID,site)%coo2(atom,:)
					other_found=.true.
					exit
				endif
			enddo
			if (.not.other_found) cycle
			distance=sqrt(sum((center_ca-other_cb)**2))
			if (distance <= ca_cb_neighbor_cutoff) then
				nb_count=nb_count+1
				nb_chain(nb_count)=chainID
				nb_site(nb_count)=site
				nb_dist(nb_count)=distance
			endif
		enddo
	enddo

	end subroutine find_cb_neighbors

	subroutine scmf_loop(group,scmf_success,failed_site,failed_aminoacid)
	implicit none
	integer, parameter :: max_repack_sweeps=5
	integer :: i,ic,chainID,feedback,nbi,nb_count
	integer :: repack_flag,n_repacked,placeholder_count,repack_sweep
	integer :: rotamer_count,transplant_count,stage1_count
	integer :: nb_chain(repeated_unit*gnum),nb_site(repeated_unit*gnum)
	real :: nb_dist(repeated_unit*gnum)
	real(kind=8) :: best_attempt_score,best_attempt_vdw
	real(kind=8) :: recovery_best_vdw,sweep_start_vdw
	character(len=4) :: aminoacid_name,nb_name,nb_base
	character(len=4) :: original_name(repeated_unit)
	character(len=4), intent(out) :: failed_aminoacid
	integer, intent(out) :: scmf_success,failed_site
	type(groupdetails), intent(inout) :: group(repeated_unit,gnum)
	type(groupdetails) :: tgroup(repeated_unit,gnum),temp_group(repeated_unit,gnum)
	type(groupdetails), allocatable :: work_group(:,:),nb_group(:,:),best_recovery_group(:,:)
	logical :: simple_target,changed_chain(repeated_unit)
	
	scmf_success=0
	failed_site=0
	failed_aminoacid=''
	tgroup=group

	! Work one sequence position at a time. First place the chain-1 target
	! identity on every changed chain, retaining the best stage-0 candidate as
	! a placeholder. Only after every placeholder at this position is present
	! do we optimize each changed chain and, when needed, repack its neighbors.
	do ic=1,gnum
		aminoacid_name=tgroup(1,ic)%gtype
		changed_chain=.false.
		original_name=''
		placeholder_count=0

		! Pass 1: prepare the target-residue volume on all peptide chains.
		do chainID=2,repeated_unit
			original_name(chainID)=tgroup(chainID,ic)%gtype
			if (trim(original_name(chainID)) == trim(aminoacid_name)) cycle
			changed_chain(chainID)=.true.
			call sequence_mutation_4scmfloop(chainID,ic,tgroup,aminoacid_name,temp_group,feedback, &
				rotamer_count,transplant_count,stage1_count,best_attempt_score,best_attempt_vdw)
			if (trim(temp_group(chainID,ic)%gtype) /= trim(aminoacid_name)) then
				write(*,'(A,I0,A,I0,3A)') 'SCMF attempt could not prepare a temporary residue on chain ', &
					chainID,', site ',ic,', for ',trim(aminoacid_name),'.'
				failed_site=ic
				failed_aminoacid=aminoacid_name
				return
			endif
			tgroup=temp_group
			placeholder_count=placeholder_count+1
		enddo
		if (placeholder_count > 0) then
			write(*,'(A,I0,A,I0,A)') 'Prepared ',placeholder_count,' temporary residue at site ',ic,'.'
		endif

		! Pass 2: optimize every changed target in the complete placeholder
		! environment. A stage-0 target remains present while its neighbors are
		! repacked, and the target is then optimized once more.
		optimization_chain_loop: do chainID=2,repeated_unit
			if (.not.changed_chain(chainID)) cycle
			call sequence_mutation_4scmfloop(chainID,ic,tgroup,aminoacid_name,temp_group,feedback, &
				rotamer_count,transplant_count,stage1_count,best_attempt_score,best_attempt_vdw)
			if (feedback == 1) then
				tgroup=temp_group
				cycle optimization_chain_loop
			endif

			allocate(work_group(repeated_unit,gnum),nb_group(repeated_unit,gnum),best_recovery_group(repeated_unit,gnum))
			work_group=temp_group
			best_recovery_group=temp_group
			recovery_best_vdw=best_attempt_vdw
			call find_cb_neighbors(work_group,chainID,ic,nb_count,nb_chain,nb_site,nb_dist)
			write(*,'(5A,I0,A,I0,A)') 'SCMF loop: replacing ',trim(original_name(chainID)), &
				' with ',trim(aminoacid_name),' on chain ',chainID,', site ',ic,'.'
			write(*,'(3A,I0,A)') 'SCMF loop: ',trim(aminoacid_name),' does not fit. Repacking ',nb_count,' neighboring side chains.'
			!do nbi=1,nb_count
			!	write(*,'(A,I0,A,I0,A,F7.3)') '  neighbor chain ',nb_chain(nbi), &
			!		', site ',nb_site(nbi),', target-CA to neighbor-CB distance ',nb_dist(nbi)
			!enddo

			do repack_sweep=1,max_repack_sweeps
				sweep_start_vdw=recovery_best_vdw
				n_repacked=0
				do nbi=1,nb_count
					nb_name=adjustl(work_group(nb_chain(nbi),nb_site(nbi))%gtype)
					nb_base=nb_name
					if (len_trim(nb_base) == 4) then
						if (nb_base(1:1) == 'N' .or. nb_base(1:1) == 'C') nb_base=nb_base(2:4)
					endif
					if (trim(nb_base)=='GLY'.or.trim(nb_base)=='ALA'.or.trim(nb_base)=='PRO') cycle
					call sequence_mutation_4scmfloop(nb_chain(nbi),nb_site(nbi),work_group, &
						nb_name,nb_group,repack_flag,rotamer_count,transplant_count,stage1_count, &
						best_attempt_score,best_attempt_vdw,.true.,.true.)
					if (repack_flag == 1) then
						work_group=nb_group
						n_repacked=n_repacked+1
						call sequence_mutation_4scmfloop(chainID,ic,work_group,aminoacid_name, &
							temp_group,feedback,rotamer_count,transplant_count,stage1_count, &
							best_attempt_score,best_attempt_vdw)
						if (feedback == 1) then
							tgroup=temp_group
							write(*,'(5A,I0,A,I0,A)') 'SCMF: replaced ',trim(original_name(chainID)), &
								' with ',trim(aminoacid_name),' on chain ',chainID,', site ',ic, &
								' after repacking neighboring side chains.'
							deallocate(work_group,nb_group,best_recovery_group)
							cycle optimization_chain_loop
						endif
						work_group=temp_group
						if (best_attempt_vdw < recovery_best_vdw) then
							recovery_best_vdw=best_attempt_vdw
							best_recovery_group=temp_group
						endif
					endif
				enddo

				call sequence_mutation_4scmfloop(chainID,ic,work_group,aminoacid_name, &
					temp_group,feedback,rotamer_count,transplant_count,stage1_count, &
					best_attempt_score,best_attempt_vdw)
				if (feedback == 1) then
					tgroup=temp_group
					write(*,'(5A,I0,A,I0,A)') 'SCMF: replaced ',trim(original_name(chainID)), &
						' with ',trim(aminoacid_name),' on chain ',chainID,', site ',ic, &
						' after repacking neighboring side chains.'
					deallocate(work_group,nb_group,best_recovery_group)
					cycle optimization_chain_loop
				endif
				write(*,'(5A,I0,A,I0,A,I0,A)') 'SCMF: replacing ',trim(original_name(chainID)), &
					' with ',trim(aminoacid_name),' on chain ',chainID,', site ',ic, &
					' still does not fit after repacking ',n_repacked,' neighboring side chains.'

				if (best_attempt_vdw < recovery_best_vdw) then
					recovery_best_vdw=best_attempt_vdw
					best_recovery_group=temp_group
				endif
				if (recovery_best_vdw >= sweep_start_vdw-1.0d-6) exit
				work_group=best_recovery_group
			enddo

			temp_group=best_recovery_group
			best_attempt_vdw=recovery_best_vdw

			call generatepdb_debug('scmf_failure_before',ic,chainID,tgroup)
			call generatepdb_debug('scmf_failure_candidate',ic,chainID,temp_group)
			write(*,*)
			write(*,'(5A,I0,A,I0)') 'SCMF attempt could not replace ', &
				trim(original_name(chainID)),' with ',trim(aminoacid_name), &
				' on chain ',chainID,' at chain site ',ic
			write(*,'(A)') 'Chain 1 target sequence (chain-site:residue):'
			do i=1,gnum
				write(*,'(1X,I0,":",A)',advance='no') i,trim(tgroup(1,i)%gtype)
			enddo
			write(*,*)
			write(*,'(A,I0)') 'Rotamers tested                  = ',rotamer_count
			simple_target=aminoacid_name=='GLY'.or.aminoacid_name=='ALA'.or. &
				aminoacid_name=='PRO'.or.aminoacid_name=='NGLY'.or. &
				aminoacid_name=='NALA'.or.aminoacid_name=='NPRO'.or. &
				aminoacid_name=='CGLY'.or.aminoacid_name=='CALA'.or. &
				aminoacid_name=='CPRO'.or.aminoacid_name=='NME'.or. &
				aminoacid_name=='NHE'.or.aminoacid_name=='ACE'
			if (simple_target) then
				write(*,'(A,I0)') 'Rotamers passing transplant      = ',transplant_count
				if (transplant_count == 0) then
					write(*,'(A)') 'Failure gate: every rotamer failed check_transplant.'
				else
					write(*,'(A,ES14.6)') 'Best attempted VDW energy        = ',best_attempt_vdw
					if (best_attempt_vdw >= 500.0d0) &
						write(*,'(A)') 'Failure gate: every valid rotamer had VDW energy >= 500.'
				endif
			else
				write(*,'(A,I0)') 'Rotamers reaching stage 1        = ',stage1_count
				if (stage1_count == 0) then
					write(*,'(A,ES14.6)') 'Lowest stage-0 local VDW         = ',best_attempt_vdw
					write(*,'(A)') 'Failure gate: no rotamer reached side-chain optimization stage 1.'
				else
					write(*,'(A,ES14.6)') 'Lowest stage-1 binding score     = ',best_attempt_score
					write(*,'(A,ES14.6)') 'Its local VDW energy             = ',best_attempt_vdw
					if (best_attempt_score >= 500.0d0) &
						write(*,'(A)') 'Failure gate: every optimized rotamer had score >= 500.'
				endif
			endif
			write(*,'(A)') 'Failure structures written to:'
			write(*,'(A)') '  pdbfiles/scmf_failure_before.pdb'
			write(*,'(A)') '  pdbfiles/scmf_failure_candidate.pdb'
			failed_site=ic
			failed_aminoacid=aminoacid_name
			deallocate(work_group,nb_group,best_recovery_group)
			return
		enddo optimization_chain_loop
	enddo
	group=tgroup
	scmf_success=1

	return
	end subroutine scmf_loop


	subroutine initial_entropy_individual(group, entropy4individual)
	implicit none
	integer							:: i, j, ii, ip, m, stage
	integer							:: ic, rotanum, obs, grade
	integer							:: W_numex(repeated_unit*atom_num), W_inb(repeated_unit*atom_num,20), W_numex4(repeated_unit*atom_num), W_inb4(repeated_unit*atom_num,60)
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	real							:: Dihedral4entropy(4), matrix(34,4)
	real							:: entropy,entropy4individual(repeated_unit,gnum)
	character*4						:: aminoacid_name

	type(groupdetails)				:: group(repeated_unit,gnum)

	call ensure_optimization_workspace()
	entropy4individual=0.0
	do ii=1, repeated_unit
		do ic=1,gnum
			aminoacid_name=group(ii,ic)%gtype
			
			if(aminoacid_name=="GLY".or.aminoacid_name=="NGLY".or.aminoacid_name=="CGLY".or.aminoacid_name=="PRO".or.aminoacid_name=="NPRO".or.aminoacid_name=="CPRO".or.  &
			   aminoacid_name=="CYS".or.aminoacid_name=="NCYS".or.aminoacid_name=="CCYS".or.aminoacid_name=="ALA".or.aminoacid_name=="NALA".or.aminoacid_name=="CALA".or.  &
			   aminoacid_name=="VAL".or.aminoacid_name=="NVAL".or.aminoacid_name=="CVAL".or.aminoacid_name=="SER".or.aminoacid_name=="NSER".or.aminoacid_name=="CSER".or.  &
			   aminoacid_name=="THR".or.aminoacid_name=="NTHR".or.aminoacid_name=="CTHR".or.aminoacid_name=="NME".or.aminoacid_name=="NHE".or.aminoacid_name=="ACE") then

				rotanum=4; matrix=0.0; obs=4; grade=4
				call entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)
				entropy4individual(ii,ic)=entropy
			else
				call findrotamer(ic, group, aminoacid_name, rotanum, opt_aa_group, ip)
				grade=aa_lib(ip)%grade

				obs=0
				matrix=0.0
				do m=1, rotanum
					call residue_replace(ii, ic, group, m, opt_aa_group, opt_temp_group)
					if(m==1) then
						call energy_parameter(opt_temp_group, opt_group_para)
						call atom_links4sidechain(ii, ic, opt_temp_group, S_numex, S_inb, S_numex4, S_inb4)
						call atom_links(opt_temp_group, W_numex, W_inb, W_numex4, W_inb4)
					endif

					stage=0
					call sidechain_optimization(stage, ii, ic, opt_temp_group, opt_group_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)

					if(stage==1) then
						obs=obs+1
						do j=1, grade
							matrix(obs,j)=Dihedral4entropy(j)								
						enddo
					endif
				enddo

				call entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)
				entropy4individual(ii,ic)=entropy
			endif		
		enddo
	enddo

	return
	end subroutine initial_entropy_individual
	
	
	subroutine sequence_optimization_nonthermal(group)
	implicit none
	integer							:: attempt, ic_1, ic_2
	integer							:: i, j, l, ii, feedback_1, feedback_2, feedback_4, flag1, flag2
	real							:: ran2
    logical							:: ic_1_NMR, ic_2_NMR, mustAA
	character*4						:: aminoacid_name_1, aminoacid_name_2
	character*4						:: group_name_1(3), group_name_2(3)
	character*4						:: rest_AA(maximum_nmr_site_num)

	type(groupdetails)				:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), tgroup(repeated_unit,gnum)

	do attempt=1, gnum
		call ran_gen(ran2,0)
		if(ran2.le.scmfswitch1) then ! mutate one AA
            feedback_4=0
            do while(feedback_4 == 0)
				call pickupsite(0, ic_1)
				call mc_choose_aminoacid(ic_1,group,aminoacid_name_1)
                call check_restrained_aminoacid(group, ic_1, aminoacid_name_1, feedback_4)
            enddo
            
			call sequence_mutation_nonthermal(ic_1, group, aminoacid_name_1, tgroup, feedback_1)

			if(feedback_1==1) then
				group=tgroup
			endif

        else ! exchange two AA. For NMR sites, there are 3 situation
19			continue            
			call pickupsite(0, ic_1)
			do while(.true.)
				call pickupsite(0, ic_2)
				if(ic_1.ne.ic_2) goto 20
			enddo
20          continue

            ic_1_NMR = .false. ! assume both ic1 and ic2 are not nmr sites
            ic_2_NMR = .false.
            do i=1, nmr_site_num
                if (ic_1 == nmr_site_ID(i)) then
					ic_1_NMR = .true.
                    exit
                endif
            enddo
            
            do i=1, nmr_site_num
                if (ic_2 == nmr_site_ID(i)) then
					ic_2_NMR = .true.
                    exit
                endif
            enddo
            
			call groupinfo(group(1,ic_1)%gtype, group_name_1, flag1)
			call groupinfo(group(1,ic_2)%gtype, group_name_2, flag2)
			if (group_name_1(1) == group_name_2(1)) goto 19 ! if the two exchanging AA are same, choose again
            if (ic_1_NMR .and. .not. ic_2_NMR) then ! one site is NMR and the other is not
                call find_rest_aa(group, rest_AA, l)
                do i=1, l
					if (group_name_2(1) == rest_AA(i)) goto 96 ! Compare base residue names; exchange if the unconstrained residue is in the remaining grouped-site pool.
                enddo
                goto 19 ! Otherwise, select another residue pair.
            else if(ic_2_NMR .and. .not. ic_1_NMR) then
                call find_rest_aa(group, rest_AA, l)
                do i=1, l
					if (group_name_1(1) == rest_AA(i)) goto 96 ! Exchange if the unconstrained residue is in the remaining grouped-site pool.
                enddo
				goto 19 ! Otherwise, select another residue pair.
            endif
                
96          continue
            aminoacid_name_1=group_name_2(flag1)
			aminoacid_name_2=group_name_1(flag2)
			call sequence_mutation_nonthermal(ic_1, group, aminoacid_name_1, temp_group, feedback_1)
	
			if(feedback_1==1) then
				call sequence_mutation_nonthermal(ic_2, temp_group, aminoacid_name_2, tgroup, feedback_2)
				if(feedback_2==1) then
					group=tgroup
				endif
			endif			
	
		endif
	enddo
	
	return
	end subroutine sequence_optimization_nonthermal	


	subroutine sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, &
		binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
		score4hydration_old, Pagg_old)
	implicit none
	integer							:: step, attempt, Num4attempt, ic_1, ic_2, categoryID_1, categoryID_2, chainID_1, chainID_2 ! ic_1 and ic_2 are randomly selected residue positions.
	integer							:: i, ii, l, j, feedback_1, feedback_2, feedback_3, feedback_4, flag1, flag2
	real							:: ran2
	real(kind=8)                    :: score_old, score_new, binding_energy_old, binding_energy_new
	real(kind=8)					:: binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old
	real(kind=8)					:: binding_vdw_new, binding_ele_new, binding_sgb_new, binding_snp_new
	real(kind=8)					:: delta_E_vdw,delta_E_ele,delta_E_sgb,delta_E_sur
	real(kind=8)					:: delta_TS,lambda_delta_Pagg,total_abs_score_change
	real							:: entropy_old, entropy_new
	real							:: score4hydration_old, score4hydration_new
	real							:: Pagg_old, Pagg_new
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum), temp_entropy4individual(repeated_unit,gnum)
	character*4						:: aminoacid_name_1, aminoacid_name_2
    character*20					:: accpt
	character*4						:: group_name_1(3), group_name_2(3)
	character*4						:: rest_AA(maximum_nmr_site_num)
	logical							:: ic_1_NMR, ic_2_NMR, trial_energy_available
    character(len=:), allocatable   :: pep_name
	type(groupdetails)				:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), tgroup(repeated_unit,gnum)
	type(energyparameters), dimension(:,:), allocatable &
									:: group_para
	integer, dimension(:), allocatable &
									:: W_numex, W_numex4
	integer, dimension(:,:), allocatable &
									:: W_inb, W_inb4
	
	call ensure_optimization_workspace()
    feedback_2=0
    feedback_3=0
    !penalty=0
	do attempt=1, 7
		feedback_1=0
		feedback_2=0
		feedback_3=0
		trial_energy_available=.false.
		call ran_gen(ran2,0)
		if(ran2.le.scmfswitch1) then
            feedback_4=0
            do while(feedback_4 == 0)
				call pickupsite(0, ic_1)
				call mc_choose_aminoacid(ic_1,group,aminoacid_name_1)
                call check_restrained_aminoacid(group, ic_1, aminoacid_name_1, feedback_4)
            enddo

			call sequence_mutation(ic_1, group, entropy4individual, aminoacid_name_1, tgroup, Tentropy4individual, feedback_1)

			if(feedback_1==1) then
				call atom_links(tgroup, opt_numex, opt_inb, opt_numex4, opt_inb4)
				call energy_parameter(tgroup, opt_group_para)

				call bindingenergy(tgroup, opt_group_para, Tentropy4individual, opt_numex, opt_inb, &
						opt_numex4, opt_inb4, score_new, binding_energy_new, binding_vdw_new, &
						binding_ele_new, binding_sgb_new, binding_snp_new, entropy_new, &
						score4hydration_new, Pagg_new)
				trial_energy_available=.true.
				delta_E_vdw=binding_vdw_new-binding_vdw_old
				delta_E_ele=binding_ele_new-binding_ele_old
				delta_E_sgb=binding_sgb_new-binding_sgb_old
				delta_E_sur=binding_snp_new-binding_snp_old
				delta_TS=real(entropy_new-entropy_old,kind=8)
				lambda_delta_Pagg=real(propensity_weighting_factor,kind=8)*real((score4hydration_new+Pagg_new)-(score4hydration_old+Pagg_old),kind=8)
				total_abs_score_change=abs(delta_E_vdw)+abs(delta_E_ele)+abs(delta_E_sgb)+abs(delta_E_sur)+abs(delta_TS)+abs(lambda_delta_Pagg)
                !Using MC-Metropolis algorithm to accept or reject sequence
				call MC_technique_sequence(score_old, binding_energy_old, binding_vdw_old, binding_ele_old, &
					binding_sgb_old, binding_snp_old, entropy_old, score4hydration_old, Pagg_old, group, &
					entropy4individual, score_new, binding_energy_new, binding_vdw_new, binding_ele_new, &
					binding_sgb_new, binding_snp_new, entropy_new, score4hydration_new, Pagg_new, tgroup, &
					Tentropy4individual, feedback_3)
			endif
			
			if(feedback_3==1) then
				accpt="Accept"
			elseif(feedback_1==1 .and. feedback_3==0) then
				accpt="Reject-MC"
			!PENALTY = PENALTY + 1                
			elseif(feedback_1==0) then
				accpt="Reject-Rotamer"
			endif

			if (trial_energy_available) then
				open(3, file="energydetails.txt",  access="append")
					write(3,4) step,attempt,score_new,binding_vdw_new,binding_ele_new,binding_sgb_new, &
						binding_snp_new,entropy_new, &
						(score4hydration_new+Pagg_new)*propensity_weighting_factor,delta_E_vdw, &
						delta_E_ele,delta_E_sgb,delta_E_sur,delta_TS,lambda_delta_Pagg, &
						total_abs_score_change,"Mutation ",accpt
					write(3,'(*(A5))') (tgroup(1,j)%gtype,j=1,gnum)
				close(3)
			endif
        else ! Exchange the amino-acid identities at two randomly selected sites.
19			continue            
			call pickupsite(0, ic_1)
			do while(.true.)
				call pickupsite(0, ic_2)
				if(ic_1.ne.ic_2) goto 20
			enddo
20          continue

            ic_1_NMR = .false. ! assume both ic1 and ic2 are not nmr sites
            ic_2_NMR = .false.
            do i=1, nmr_site_num
                if (ic_1 == nmr_site_ID(i)) ic_1_NMR = .true.
            enddo
            
            do i=1, nmr_site_num
                if (ic_2 == nmr_site_ID(i)) ic_2_NMR = .true.
            enddo
            
			call groupinfo(group(1,ic_1)%gtype, group_name_1, flag1)
			call groupinfo(group(1,ic_2)%gtype, group_name_2, flag2)
			if (group_name_1(1) == group_name_2(1)) goto 19 ! if the two exchanging AA are same, choose again
            if (ic_1_NMR .and. .not. ic_2_NMR) then ! one site is NMR and the other is not
                call find_rest_aa(group, rest_AA, l)
                do i=1, l
					if (group_name_2(1) == rest_AA(i)) goto 96 ! Compare base residue names; exchange if the unconstrained residue is in the remaining grouped-site pool.
                enddo
                goto 19 ! Otherwise, select another residue pair.
            else if(ic_2_NMR .and. .not. ic_1_NMR) then
                call find_rest_aa(group, rest_AA, l)
                do i=1, l
					if (group_name_1(1) == rest_AA(i)) goto 96 ! Exchange if the unconstrained residue is in the remaining grouped-site pool.
                enddo
				goto 19 ! Otherwise, select another residue pair.
            endif               
96          continue
            
            aminoacid_name_1=group_name_2(flag1)
			aminoacid_name_2=group_name_1(flag2)
                 
			call sequence_mutation(ic_1, group, entropy4individual, aminoacid_name_1, temp_group, temp_entropy4individual, feedback_1)
			
			if(feedback_1==1) then

				call sequence_mutation(ic_2, temp_group, temp_entropy4individual, aminoacid_name_2, tgroup, Tentropy4individual, feedback_2)

				if(feedback_2==1) then
					call atom_links(tgroup, opt_numex, opt_inb, opt_numex4, opt_inb4)
					call energy_parameter(tgroup, opt_group_para)
		
					call bindingenergy(tgroup, opt_group_para, Tentropy4individual, opt_numex, opt_inb, &
							opt_numex4, opt_inb4, score_new, binding_energy_new, binding_vdw_new, &
							binding_ele_new, binding_sgb_new, binding_snp_new, entropy_new, &
							score4hydration_new, Pagg_new)
					trial_energy_available=.true.
					delta_E_vdw=binding_vdw_new-binding_vdw_old
					delta_E_ele=binding_ele_new-binding_ele_old
					delta_E_sgb=binding_sgb_new-binding_sgb_old
					delta_E_sur=binding_snp_new-binding_snp_old
					delta_TS=real(entropy_new-entropy_old,kind=8)
					lambda_delta_Pagg=real(propensity_weighting_factor,kind=8)*real((score4hydration_new+Pagg_new)-(score4hydration_old+Pagg_old),kind=8)
					total_abs_score_change=abs(delta_E_vdw)+abs(delta_E_ele)+abs(delta_E_sgb)+abs(delta_E_sur)+abs(delta_TS)+abs(lambda_delta_Pagg)
					call MC_technique_sequence(score_old, binding_energy_old, binding_vdw_old, binding_ele_old, &
						binding_sgb_old, binding_snp_old, entropy_old, score4hydration_old, Pagg_old, group, &
						entropy4individual, score_new, binding_energy_new, binding_vdw_new, binding_ele_new, &
						binding_sgb_new, binding_snp_new, entropy_new, score4hydration_new, Pagg_new, tgroup, &
						Tentropy4individual, feedback_3)
				endif
			endif

			if(feedback_3==1) then
				accpt="Accept" ! Accepted trial: report the new peptide and score.
			elseif(feedback_2==1 .and. feedback_3==0) then
				accpt="Reject-MC" ! Rejected trial: report the previous peptide and score.
			elseif(feedback_2==0) then
				accpt="Reject-Rotamer"
			endif

			if (trial_energy_available) then
				open(3, file="energydetails.txt",  access="append")
					write(3,4) step,attempt,score_new,binding_vdw_new,binding_ele_new,binding_sgb_new, &
						binding_snp_new,entropy_new, &
						(score4hydration_new+Pagg_new)*propensity_weighting_factor,delta_E_vdw, &
						delta_E_ele,delta_E_sgb,delta_E_sur,delta_TS,lambda_delta_Pagg, &
						total_abs_score_change,"Exchange ",accpt
					write(3,'(*(A5))') (tgroup(1,j)%gtype,j=1,gnum)
				close(3)
			endif
		endif
		
		if(score_old.lt.energy_min(1)) then
			do i=(num4pdbsave-1),1,-1
				energy_min(i+1)=energy_min(i)
			enddo
			energy_min(1)=score_old
            allocate(character(len=gnum) :: pep_name)
            call convert_AA_name(group, pep_name)
			open(2, file="minimum_energy.txt", access="append")
				write(2,*) step, attempt, score_old, pep_name
			close(2)
			call generatepdb(step, attempt, group)
            deallocate(pep_name)
		elseif(score_old.lt.energy_min(num4pdbsave)) then
			do i=1,num4pdbsave
				if(score_old.eq.energy_min(i)) then
					goto 50
				elseif(score_old.lt.energy_min(i)) then
					do j=(num4pdbsave-1),i,-1
						energy_min(j+1)=energy_min(j)
					enddo
					energy_min(i)=score_old
					goto 60
				endif
			enddo
60			continue
			call generatepdb(step, attempt, group)
50			continue
        endif
        !if (accpt=="Accept") exit
	enddo
4	format(i7,i7,14f15.4,a15,a15)

	return
	end subroutine sequence_optimization
	

	real function angle_difference(angle_new, angle_old)
	implicit none
	real					:: angle_new, angle_old

	angle_difference=angle_new-angle_old
	do while(angle_difference.gt.180.0)
		angle_difference=angle_difference-360.0
	enddo
	do while(angle_difference.le.-180.0)
		angle_difference=angle_difference+360.0
	enddo

	return
	end function angle_difference


	subroutine get_backbone_atom(group, chainID, site, atom_name, coord, found)
	implicit none
	integer					:: chainID, site, i
	logical					:: found
	real					:: coord(3)
	character(len=*)		:: atom_name
	type(groupdetails)		:: group(repeated_unit,gnum)

	found=.false.
	coord=0.0
	do i=1, group(chainID,site)%cnum1
		if(group(chainID,site)%atype1(i)==atom_name) then
			coord=group(chainID,site)%coo1(i,:)
			found=.true.
			return
		endif
	enddo
	do i=1, group(chainID,site)%cnum3
		if(group(chainID,site)%atype3(i)==atom_name) then
			coord=group(chainID,site)%coo3(i,:)
			found=.true.
			return
		endif
	enddo

	return
	end subroutine get_backbone_atom


	subroutine rotate_point(point, origin, m)
	implicit none
	real					:: point(3), origin(3), m(3,3), temp(3)

	temp=point-origin
	temp=matmul(temp,m)
	point=anint((temp+origin)*1000.0)/1000.0

	return
	end subroutine rotate_point


	subroutine rotate_phi_fragment(Tgroup, chainID, site, last_site, phi_new)
	implicit none
	integer					:: chainID, site, last_site, j, k
	logical					:: found
	real					:: phi_new, phi_old, delta_phi, cos_angle, sin_angle
	real					:: C0(3), N(3), CA(3), C(3), rotaxis(3), m(3,3), norm
	type(groupdetails)		:: Tgroup(repeated_unit,gnum)

	call get_backbone_atom(Tgroup, chainID, site-1, "C", C0, found); if(.not.found) return
	call get_backbone_atom(Tgroup, chainID, site, "N", N, found); if(.not.found) return
	call get_backbone_atom(Tgroup, chainID, site, "CA", CA, found); if(.not.found) return
	call get_backbone_atom(Tgroup, chainID, site, "C", C, found); if(.not.found) return
	call phipsiomg_angle(C0, N, CA, C, phi_old)

	delta_phi=-(phi_new-phi_old)
	cos_angle=cosd(delta_phi)
	sin_angle=sind(delta_phi)
	rotaxis=CA-N
	norm=sqrt(dot_product(rotaxis,rotaxis))
	if(norm.le.1.0e-6) return
	rotaxis=rotaxis/norm
	call axisrotation(rotaxis, cos_angle, sin_angle, m)

	do k=1, Tgroup(chainID,site)%cnum1
		if(Tgroup(chainID,site)%atype1(k)=="HA".or.Tgroup(chainID,site)%atype1(k)=="HA2".or. &
			Tgroup(chainID,site)%atype1(k)=="HA3") then
			call rotate_point(Tgroup(chainID,site)%coo1(k,:), N, m)
		endif
	enddo
	do k=1, Tgroup(chainID,site)%cnum2
		call rotate_point(Tgroup(chainID,site)%coo2(k,:), N, m)
	enddo
	do k=1, Tgroup(chainID,site)%cnum3
		call rotate_point(Tgroup(chainID,site)%coo3(k,:), N, m)
	enddo

	do j=site+1, last_site
		do k=1, Tgroup(chainID,j)%cnum1
			call rotate_point(Tgroup(chainID,j)%coo1(k,:), N, m)
		enddo
		do k=1, Tgroup(chainID,j)%cnum2
			call rotate_point(Tgroup(chainID,j)%coo2(k,:), N, m)
		enddo
		if(j.ne.last_site) then
			do k=1, Tgroup(chainID,j)%cnum3
				call rotate_point(Tgroup(chainID,j)%coo3(k,:), N, m)
			enddo
		endif
	enddo

	return
	end subroutine rotate_phi_fragment


	subroutine rotate_psi_fragment(Tgroup, chainID, site, last_site, psi_new)
	implicit none
	integer					:: chainID, site, last_site, j, k
	logical					:: found
	real					:: psi_new, psi_old, delta_psi, cos_angle, sin_angle
	real					:: Nnext(3), N(3), CA(3), C(3), rotaxis(3), m(3,3), norm
	type(groupdetails)		:: Tgroup(repeated_unit,gnum)

	call get_backbone_atom(Tgroup, chainID, site, "N", N, found); if(.not.found) return
	call get_backbone_atom(Tgroup, chainID, site, "CA", CA, found); if(.not.found) return
	call get_backbone_atom(Tgroup, chainID, site, "C", C, found); if(.not.found) return
	call get_backbone_atom(Tgroup, chainID, site+1, "N", Nnext, found); if(.not.found) return
	call phipsiomg_angle(N, CA, C, Nnext, psi_old)

	delta_psi=-(psi_new-psi_old)
	cos_angle=cosd(delta_psi)
	sin_angle=sind(delta_psi)
	rotaxis=C-CA
	norm=sqrt(dot_product(rotaxis,rotaxis))
	if(norm.le.1.0e-6) return
	rotaxis=rotaxis/norm
	call axisrotation(rotaxis, cos_angle, sin_angle, m)

	do k=1, Tgroup(chainID,site)%cnum3
		call rotate_point(Tgroup(chainID,site)%coo3(k,:), CA, m)
	enddo
	do j=site+1, last_site
		do k=1, Tgroup(chainID,j)%cnum1
			call rotate_point(Tgroup(chainID,j)%coo1(k,:), CA, m)
		enddo
		do k=1, Tgroup(chainID,j)%cnum2
			call rotate_point(Tgroup(chainID,j)%coo2(k,:), CA, m)
		enddo
		if(j.ne.last_site) then
			do k=1, Tgroup(chainID,j)%cnum3
				call rotate_point(Tgroup(chainID,j)%coo3(k,:), CA, m)
			enddo
		endif
	enddo

	return
	end subroutine rotate_psi_fragment


	subroutine backbone_rotation_center_chain(group, chainID, ran_resi, phi1, psi1, phi2, psi2, Tgroup)
	implicit none
	integer					:: chainID, ran_resi(3)
	real					:: phi1, psi1, phi2, psi2
	type(groupdetails)		:: group(repeated_unit,gnum), Tgroup(repeated_unit,gnum)

	Tgroup=group
	call rotate_phi_fragment(Tgroup, chainID, ran_resi(1), ran_resi(3), phi1)
	call rotate_psi_fragment(Tgroup, chainID, ran_resi(1), ran_resi(3), psi1)
	call rotate_phi_fragment(Tgroup, chainID, ran_resi(2), ran_resi(3), phi2)
	call rotate_psi_fragment(Tgroup, chainID, ran_resi(2), ran_resi(3), psi2)

	return
	end subroutine backbone_rotation_center_chain


	subroutine chain_backbone_rmsd_between(group_old, group_new, chainID, first_site, last_site, rmsd)
	implicit none
	integer					:: chainID, first_site, last_site, site, count
	logical					:: found1, found2
	real					:: rmsd, p1(3), p2(3)
	type(groupdetails)		:: group_old(repeated_unit,gnum), group_new(repeated_unit,gnum)

	rmsd=0.0
	count=0
	do site=first_site, last_site
		call get_backbone_atom(group_old, chainID, site, "N", p1, found1)
		call get_backbone_atom(group_new, chainID, site, "N", p2, found2)
		if(found1.and.found2) then
			rmsd=rmsd+sum((p1-p2)**2)
			count=count+1
		endif
		call get_backbone_atom(group_old, chainID, site, "CA", p1, found1)
		call get_backbone_atom(group_new, chainID, site, "CA", p2, found2)
		if(found1.and.found2) then
			rmsd=rmsd+sum((p1-p2)**2)
			count=count+1
		endif
		call get_backbone_atom(group_old, chainID, site, "C", p1, found1)
		call get_backbone_atom(group_new, chainID, site, "C", p2, found2)
		if(found1.and.found2) then
			rmsd=rmsd+sum((p1-p2)**2)
			count=count+1
		endif
	enddo
	if(count.gt.0) rmsd=sqrt(rmsd/count)

	return
	end subroutine chain_backbone_rmsd_between


	subroutine pick_conrot_sites(ran_resi, feedback)
	implicit none
	integer					:: ran_resi(3), feedback, i, j, temp
	logical					:: duplicate
	real					:: ran2

	feedback=0
	if(gnum.lt.5) return

	do i=1, 3
		do
			call ran_gen(ran2,0)
			ran_resi(i)=int(ran2*(gnum-2)-1.0e-3)+2
			if(ran_resi(i).gt.gnum-1) ran_resi(i)=gnum-1
			duplicate=.false.
			do j=1, i-1
				if(ran_resi(i).eq.ran_resi(j)) duplicate=.true.
			enddo
			if(.not.duplicate) exit
		enddo
	enddo

	do i=1, 2
		do j=i+1, 3
			if(ran_resi(i).gt.ran_resi(j)) then
				temp=ran_resi(i)
				ran_resi(i)=ran_resi(j)
				ran_resi(j)=temp
			endif
		enddo
	enddo
	feedback=1

	return
	end subroutine pick_conrot_sites


	subroutine concertedrotation_center_chain(group, chainID, ran_resi, Tgroup_best, feedback)
	implicit none
	integer					:: chainID, ran_resi(3), feedback
	integer					:: scan_i, scan_num, i, j, k, H, psi_num, flag
	logical					:: found
	real					:: C0(3), N1(3), CA1(3), C1(3), N2(3), CA2(3), C2(3)
	real					:: N3(3), CA3(3), C3(3), N4(3), CA3_target(3), CA3_temp(3)
	real					:: c11, c22, c33, c44, dist, rmsd
	real					:: omg1, omg2, theta_omg0, theta_phi1, theta_psi1
	real					:: theta_omg1, theta_phi2, theta_psi2, theta_omg2
	real					:: l_omg0, l_phi1, l_psi1, l_omg1, l_phi2, l_psi2, l_omg2, l_phi3
	real					:: U_phi(3), U_omg(3), U_omg0(3), U_phi1(3), U_psi1(3)
	real					:: U_omg1(3), U_phi2(3), U_psi2(3), U_omg2(3), U_phi3(3)
	real					:: Lomg1(3), Lomg2(3), Q_omg1(3), Q_omg2(3), e(3), t(3)
	real					:: t_prime(3), q2(3), q2_prime(3), mvec(3), nvec(3)
	real					:: Tlab(3,3), Tphi1(3,3), Tpsi1(3,3), Tomg1(3,3)
	real					:: r_CA_N(3), r1(3), delta_phi, phi1_old, psi1_old, phi2_old, psi2_old
	real					:: phi1, phi2, phi(2), psi(2), psi1_sol(2)
	real					:: phi1_temp, psi1_temp, phi2_temp, psi2_temp
	type(groupdetails)		:: group(repeated_unit,gnum), Tgroup_best(repeated_unit,gnum), Tgroup_new(repeated_unit,gnum)

	feedback=0
	Tgroup_best=group

	call get_backbone_atom(group, chainID, ran_resi(1)-1, "C", C0, found); if(.not.found) return
	call get_backbone_atom(group, chainID, ran_resi(1), "N", N1, found); if(.not.found) return
	call get_backbone_atom(group, chainID, ran_resi(1), "CA", CA1, found); if(.not.found) return
	call get_backbone_atom(group, chainID, ran_resi(1), "C", C1, found); if(.not.found) return
	call get_backbone_atom(group, chainID, ran_resi(2), "N", N2, found); if(.not.found) return
	call get_backbone_atom(group, chainID, ran_resi(2), "CA", CA2, found); if(.not.found) return
	call get_backbone_atom(group, chainID, ran_resi(2), "C", C2, found); if(.not.found) return
	call get_backbone_atom(group, chainID, ran_resi(3), "N", N3, found); if(.not.found) return
	call get_backbone_atom(group, chainID, ran_resi(3), "CA", CA3, found); if(.not.found) return
	call get_backbone_atom(group, chainID, ran_resi(3), "C", C3, found); if(.not.found) return
	call get_backbone_atom(group, chainID, ran_resi(3)+1, "N", N4, found); if(.not.found) return

	U_omg0=N1-C0; l_omg0=sqrt(dot_product(U_omg0,U_omg0))
	U_phi1=CA1-N1; U_psi1=C1-CA1; U_omg1=N2-C1
	l_phi1=sqrt(dot_product(U_phi1,U_phi1)); l_psi1=sqrt(dot_product(U_psi1,U_psi1))
	l_omg1=sqrt(dot_product(U_omg1,U_omg1))
	U_phi2=CA2-N2; U_psi2=C2-CA2; U_omg2=N3-C2
	l_phi2=sqrt(dot_product(U_phi2,U_phi2)); l_psi2=sqrt(dot_product(U_psi2,U_psi2))
	l_omg2=sqrt(dot_product(U_omg2,U_omg2))
	U_phi3=CA3-N3; l_phi3=sqrt(dot_product(U_phi3,U_phi3))
	if(min(l_omg0,l_phi1,l_psi1,l_omg1,l_phi2,l_psi2,l_omg2,l_phi3).le.1.0e-6) return

	call phipsiomg_angle(C0, N1, CA1, C1, phi1_old)
	call phipsiomg_angle(N1, CA1, C1, N2, psi1_old)
	call phipsiomg_angle(C1, N2, CA2, C2, phi2_old)
	call phipsiomg_angle(N2, CA2, C2, N3, psi2_old)
	call phipsiomg_angle(CA1, C1, N2, CA2, omg1)
	call phipsiomg_angle(CA2, C2, N3, CA3, omg2)

	theta_omg0=acosd(dot_product(U_omg0,U_phi1)/(l_omg0*l_phi1))
	theta_phi1=acosd(dot_product(U_phi1,U_psi1)/(l_phi1*l_psi1))
	theta_psi1=acosd(dot_product(U_psi1,U_omg1)/(l_psi1*l_omg1))
	theta_omg1=acosd(dot_product(U_omg1,U_phi2)/(l_omg1*l_phi2))
	theta_phi2=acosd(dot_product(U_phi2,U_psi2)/(l_phi2*l_psi2))
	theta_psi2=acosd(dot_product(U_psi2,U_omg2)/(l_psi2*l_omg2))
	theta_omg2=acosd(dot_product(U_omg2,U_phi3)/(l_omg2*l_phi3))

	Tlab(1,:)=U_phi1/l_phi1
	U_phi=Tlab(1,:)
	U_omg=U_omg0/l_omg0
	if(abs(sind(theta_omg0)).le.1.0e-6) return
	Tlab(3,1)=(U_phi(2)*U_omg(3)-U_phi(3)*U_omg(2))/sind(theta_omg0)
	Tlab(3,2)=(U_phi(3)*U_omg(1)-U_phi(1)*U_omg(3))/sind(theta_omg0)
	Tlab(3,3)=(U_phi(1)*U_omg(2)-U_phi(2)*U_omg(1))/sind(theta_omg0)
	Tlab(2,1)=Tlab(3,2)*Tlab(1,3)-Tlab(3,3)*Tlab(1,2)
	Tlab(2,2)=Tlab(3,3)*Tlab(1,1)-Tlab(3,1)*Tlab(1,3)
	Tlab(2,3)=Tlab(3,1)*Tlab(1,2)-Tlab(3,2)*Tlab(1,1)

	CA3_target=CA3
	r_CA_N=CA3_target-N1
	r1=matmul(Tlab,r_CA_N)
	r1(1)=r1(1)-l_phi1

	Lomg1(1)=l_omg1+l_phi2*cosd(theta_omg1)
	Lomg1(2)=-l_phi2*sind(theta_omg1)*cosd(omg1)
	Lomg1(3)=l_phi2*sind(theta_omg1)*sind(omg1)
	Lomg2(1)=l_omg2+l_phi3*cosd(theta_omg2)
	Lomg2(2)=-l_phi3*sind(theta_omg2)*cosd(omg2)
	Lomg2(3)=l_phi3*sind(theta_omg2)*sind(omg2)
	Q_omg1(1)=l_psi1*cosd(theta_psi1)+Lomg1(1)
	Q_omg1(2)=l_psi1*sind(theta_psi1)+Lomg1(2)
	Q_omg1(3)=Lomg1(3)
	Q_omg2(1)=l_psi2*cosd(theta_psi2)+Lomg2(1)
	Q_omg2(2)=l_psi2*sind(theta_psi2)+Lomg2(2)
	Q_omg2(3)=Lomg2(3)
	e=(/1.0,0.0,0.0/)

	scan_num=int(2.0*CONROT_angle_limit/CONROT_step+0.5)
	do scan_i=0, scan_num
		delta_phi=-CONROT_angle_limit+CONROT_step*scan_i
		if(abs(delta_phi).le.1.0e-6) cycle
		phi1=phi1_old+delta_phi
		call transformatrix(theta_phi1, phi1, Tphi1)
		t=matmul(transpose(Tphi1), r1)

		t_prime(1)=sum(t**2)+sum(Q_omg1**2)-sum(Q_omg2**2)- &
			2.0*t(1)*(cosd(theta_psi1)*Q_omg1(1)+sind(theta_psi1)*Q_omg1(2))
		t_prime(2)=2.0*t(2)*(cosd(theta_psi1)*Q_omg1(2)-sind(theta_psi1)*Q_omg1(1))+2.0*t(3)*Q_omg1(3)
		t_prime(3)=2.0*t(3)*(cosd(theta_psi1)*Q_omg1(2)-sind(theta_psi1)*Q_omg1(1))-2.0*t(2)*Q_omg1(3)

		if(sqrt(t_prime(2)**2+t_prime(3)**2).le.1.0e-6) cycle
		c11=t_prime(1)/sqrt(t_prime(2)**2+t_prime(3)**2)
		if(abs(c11).gt.1.0) cycle
		if(abs(t_prime(3)).le.1.0e-6) cycle
		if(t_prime(3).gt.0.0) then
			H=0
		else
			H=180
		endif
		psi(1)=atand(t_prime(2)/t_prime(3))-asind(c11)+H
		psi(2)=atand(t_prime(2)/t_prime(3))+asind(c11)-180.0+H
		do j=1, 2
			do while(psi(j).le.-179.5.or.psi(j).gt.180.5)
				if(psi(j).le.-179.5) psi(j)=psi(j)+360.0
				if(psi(j).gt.180.5) psi(j)=psi(j)-360.0
			enddo
		enddo
		psi_num=2
		psi1_sol=psi

		do i=1, psi_num
			! CONROT debug: do not reject by dependent-angle size while inspecting trial structures.
			!if(abs(angle_difference(psi1_sol(i),psi1_old)).gt.CONROT_angle_limit+0.001) cycle
			call transformatrix(theta_psi1, psi1_sol(i), Tpsi1)
			q2=matmul(transpose(Tpsi1), t)
			q2=q2-Q_omg1
			call transformatrix(theta_omg1, omg1, Tomg1)
			q2_prime=matmul(transpose(Tomg1), q2)

			mvec(1)=sind(theta_phi2)*(sind(theta_psi2)*Q_omg2(1)-cosd(theta_psi2)*Q_omg2(2))
			mvec(2)=sind(theta_phi2)*Q_omg2(3)
			mvec(3)=cosd(theta_phi2)*(cosd(theta_psi2)*Q_omg2(1)+sind(theta_psi2)*Q_omg2(2))-q2_prime(1)
			if(sqrt(mvec(1)**2+mvec(2)**2).le.1.0e-6) cycle
			c22=mvec(3)/sqrt(mvec(1)**2+mvec(2)**2)
			if(abs(c22).gt.1.0) cycle
			if(abs(mvec(2)).le.1.0e-6) cycle
			if(mvec(2).gt.0.0) then
				H=0
			else
				H=180
			endif
			psi(1)=atand(mvec(1)/mvec(2))-asind(c22)+H
			psi(2)=atand(mvec(1)/mvec(2))+asind(c22)-180.0+H
			do j=1, 2
				do while(psi(j).le.-179.5.or.psi(j).gt.180.5)
					if(psi(j).le.-179.5) psi(j)=psi(j)+360.0
					if(psi(j).gt.180.5) psi(j)=psi(j)-360.0
				enddo
			enddo

			do j=1, 2
				! CONROT debug: do not reject by dependent-angle size while inspecting trial structures.
				!if(abs(angle_difference(psi(j),psi2_old)).gt.CONROT_angle_limit+0.001) cycle
				nvec(1)=(sind(theta_phi2)*cosd(theta_psi2)+cosd(theta_phi2)*sind(theta_psi2)*cosd(psi(j)))*Q_omg2(1)
				nvec(1)=nvec(1)+(sind(theta_phi2)*sind(theta_psi2)-cosd(theta_phi2)*cosd(theta_psi2)* &
					cosd(psi(j)))*Q_omg2(2)
				nvec(1)=nvec(1)-cosd(theta_phi2)*sind(psi(j))*Q_omg2(3)
				nvec(2)=sind(theta_psi2)*sind(psi(j))*Q_omg2(1)-cosd(theta_psi2)*sind(psi(j))*Q_omg2(2)+ &
					cosd(psi(j))*Q_omg2(3)
				nvec(3)=-q2_prime(2)
				c44=q2_prime(3)
				if(sqrt(nvec(1)**2+nvec(2)**2).le.1.0e-6) cycle
				c33=nvec(3)/sqrt(nvec(1)**2+nvec(2)**2)
				if(abs(c33).gt.1.0) cycle
				if(abs(nvec(2)).le.1.0e-6) cycle
				if(nvec(2).gt.0.0) then
					H=0
				else
					H=180
				endif
				phi(1)=atand(nvec(1)/nvec(2))-asind(c33)+H
				phi(2)=atand(nvec(1)/nvec(2))+asind(c33)-180.0+H
				do k=1, 2
					do while(phi(k).le.-179.5.or.phi(k).gt.180.5)
						if(phi(k).le.-179.5) phi(k)=phi(k)+360.0
						if(phi(k).gt.180.5) phi(k)=phi(k)-360.0
					enddo
				enddo

				flag=0
				do k=1, 2
					if(abs((nvec(1)*sind(phi(k))+nvec(2)*cosd(phi(k)))-c44).le.0.005) then
						phi2=phi(k)
						flag=1
					endif
				enddo
				if(flag.eq.0) cycle
				! CONROT debug: do not reject by dependent-angle size while inspecting trial structures.
				!if(abs(angle_difference(phi2,phi2_old)).gt.CONROT_angle_limit+0.001) cycle

				call backbone_rotation_center_chain(group, chainID, ran_resi, phi1, psi1_sol(i), phi2, psi(j), Tgroup_new)
				call get_backbone_atom(Tgroup_new, chainID, ran_resi(3), "CA", CA3_temp, found); if(.not.found) cycle
				dist=sqrt(sum((CA3_temp-CA3_target)**2))
				! CONROT debug: write trial structures even if the fragment does not close tightly.
				!if(dist.gt.CONROT_closure_tol) cycle

				call get_backbone_atom(Tgroup_new, chainID, ran_resi(1)-1, "C", C0, found); if(.not.found) cycle
				call get_backbone_atom(Tgroup_new, chainID, ran_resi(1), "N", N1, found); if(.not.found) cycle
				call get_backbone_atom(Tgroup_new, chainID, ran_resi(1), "CA", CA1, found); if(.not.found) cycle
				call get_backbone_atom(Tgroup_new, chainID, ran_resi(1), "C", C1, found); if(.not.found) cycle
				call get_backbone_atom(Tgroup_new, chainID, ran_resi(2), "N", N2, found); if(.not.found) cycle
				call get_backbone_atom(Tgroup_new, chainID, ran_resi(2), "CA", CA2, found); if(.not.found) cycle
				call get_backbone_atom(Tgroup_new, chainID, ran_resi(2), "C", C2, found); if(.not.found) cycle
				call get_backbone_atom(Tgroup_new, chainID, ran_resi(3), "N", N3, found); if(.not.found) cycle
				call phipsiomg_angle(C0, N1, CA1, C1, phi1_temp)
				call phipsiomg_angle(N1, CA1, C1, N2, psi1_temp)
				call phipsiomg_angle(C1, N2, CA2, C2, phi2_temp)
				call phipsiomg_angle(N2, CA2, C2, N3, psi2_temp)
				! CONROT debug: do not reject by final angle re-checks while inspecting trial structures.
				!if(abs(angle_difference(phi1_temp,phi1_old)).gt.CONROT_angle_limit+0.05) cycle
				!if(abs(angle_difference(psi1_temp,psi1_old)).gt.CONROT_angle_limit+0.05) cycle
				!if(abs(angle_difference(phi2_temp,phi2_old)).gt.CONROT_angle_limit+0.05) cycle
				!if(abs(angle_difference(psi2_temp,psi2_old)).gt.CONROT_angle_limit+0.05) cycle

				! CONROT debug: do not reject by local backbone RMSD while inspecting trial structures.
				!call chain_backbone_rmsd_between(group, Tgroup_new, chainID, ran_resi(1), ran_resi(3), rmsd)
				!if(rmsd.gt.CONROT_rmsd_limit) cycle

				Tgroup_best=Tgroup_new
				feedback=1
				return
			enddo
		enddo
	enddo

	return
	end subroutine concertedrotation_center_chain


	subroutine conrot_perturbation(group, Tgroup, feedback)
	implicit none
	integer					:: chainID, ran_resi(3), feedback, feedback_site, feedback_chain
	type(groupdetails)		:: group(repeated_unit,gnum), Tgroup(repeated_unit,gnum), Tchain_group(repeated_unit,gnum)

	feedback=0
	Tgroup=group
	do chainID=1, repeated_unit
		call pick_conrot_sites(ran_resi, feedback_site)
		if(feedback_site.eq.0) return
		call concertedrotation_center_chain(Tgroup, chainID, ran_resi, Tchain_group, feedback_chain)
		if(feedback_chain.eq.0) return
		Tgroup=Tchain_group
	enddo
	feedback=1

	return
	end subroutine conrot_perturbation


	subroutine backbone_optimization_conrot(group, step, entropy4individual, score_old, binding_energy_old, &
		binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
		score4hydration_old, Pagg_old)
	implicit none
	integer						:: step, attempt, ic_1, feedback_1, feedback_2, feedback_3
	integer						:: flag1, i, j
	real(kind=8)                :: score_old, score_new
	real(kind=8)                :: binding_energy_old, binding_energy_new
	real(kind=8)				:: binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old
	real(kind=8)				:: binding_vdw_new, binding_ele_new, binding_sgb_new, binding_snp_new
	real						:: entropy_old, entropy_new
	real						:: score4hydration_old, score4hydration_new
	real						:: Pagg_old, Pagg_new
	real						:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum)
	real						:: temp_entropy4individual(repeated_unit,gnum)
	character*20				:: backbone_change, accpt
	character*4					:: aminoacid_name_1
	character*4					:: group_name_1(3)
	character(len=:), allocatable   :: pep_name
	type(groupdetails)			:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), Tgroup(repeated_unit,gnum)
	type(groupdetails)			:: repack_group(repeated_unit,gnum)

	call ensure_optimization_workspace()
	backbone_change="CONROT"

	do attempt=1, 7
		feedback_1=0
		feedback_2=0
		call conrot_perturbation(group, Tgroup, feedback_3)
		if(feedback_3.eq.0) then
			!call generatepdb_debug("conrot_failed", step, attempt, Tgroup)
			accpt="Reject-CONROT"
			!goto 40
        endif
		! Debug PDB output disabled; normal accepted/minimum-energy PDB output remains active.
		!call generatepdb_debug("conrot_failed", step, attempt, Tgroup)
		!call generatepdb_debug("conrot_raw", step, attempt, Tgroup)

		temp_group=Tgroup
		temp_entropy4individual=entropy4individual
		do ic_1=1, gnum
			if(group(1,ic_1)%gtype=="ACE".or.group(1,ic_1)%gtype=="NME".or. &
				group(1,ic_1)%gtype=="NHE") cycle
			call groupinfo(group(1,ic_1)%gtype, group_name_1, flag1)
			aminoacid_name_1=group_name_1(flag1)
			call sequence_mutation(ic_1, temp_group, temp_entropy4individual, aminoacid_name_1, &
				repack_group, Tentropy4individual, feedback_1)
			if(feedback_1.eq.0) then
				goto 20
			endif
			temp_group=repack_group
			temp_entropy4individual=Tentropy4individual
		enddo

		call atom_links(temp_group, opt_numex, opt_inb, opt_numex4, opt_inb4)
		call energy_parameter(temp_group, opt_group_para)
		call bindingenergy(temp_group, opt_group_para, temp_entropy4individual, opt_numex, opt_inb, &
				opt_numex4, opt_inb4, score_new, binding_energy_new, binding_vdw_new, &
				binding_ele_new, binding_sgb_new, binding_snp_new, entropy_new, &
				score4hydration_new, Pagg_new)
		call MC_technique_sheet(score_old, binding_energy_old, binding_vdw_old, binding_ele_old, &
			binding_sgb_old, binding_snp_old, entropy_old, score4hydration_old, Pagg_old, group, &
			entropy4individual, score_new, binding_energy_new, binding_vdw_new, binding_ele_new, &
			binding_sgb_new, binding_snp_new, entropy_new, score4hydration_new, Pagg_new, temp_group, &
			temp_entropy4individual, feedback_2)
		feedback_2=1 ! debug
		if(feedback_2.eq.1) flag4sheet=1

20 		continue
		if(feedback_2.eq.1) then
			accpt="Accept"
		elseif(feedback_1.eq.1.and.feedback_2.eq.0) then
			accpt="Reject-MC"
		elseif(feedback_1.eq.0) then
			accpt="Reject-Rotamer"
		endif

40 		continue
		open(3, file="energydetails.txt",  access="append")
		write(3,4) step, attempt, score_old, binding_vdw_old, binding_ele_old, &
			binding_sgb_old, binding_snp_old, entropy_old, &
			(score4hydration_old+Pagg_old)*propensity_weighting_factor, backbone_change, accpt
		write(3,'(*(A5))') (group(1,j)%gtype, j=1, gnum)
		close(3)

		if(score_old.lt.energy_min(1)) then
			do i=(num4pdbsave-1),1,-1
				energy_min(i+1)=energy_min(i)
			enddo
			energy_min(1)=score_old
			allocate(character(len=gnum) :: pep_name)
			call convert_AA_name(group, pep_name)
			open(2, file="minimum_energy.txt", access="append")
				write(2,*) step, attempt, score_old, pep_name
			close(2)
			call generatepdb(step, attempt, group)
			deallocate(pep_name)
		elseif(score_old.lt.energy_min(num4pdbsave)) then
			do i=1,num4pdbsave
				if(score_old.eq.energy_min(i)) then
					goto 50
				elseif(score_old.lt.energy_min(i)) then
					do j=(num4pdbsave-1),i,-1
						energy_min(j+1)=energy_min(j)
					enddo
					energy_min(i)=score_old
					goto 60
				endif
			enddo
60			continue
			call generatepdb(step, attempt, group)
50			continue
		endif
	enddo

4	format(i7,i7,7f15.4,a15,a15)
	return
	end subroutine backbone_optimization_conrot
	
	
	subroutine sheet_optimization(group, step, entropy4individual, score_old, binding_energy_old, &
		binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
		score4hydration_old, Pagg_old)
	implicit none
	integer						:: step, attempt, ip1, ip2, ic_1, feedback_1, feedback_2, feedback_3
	integer						:: flag1, i, j, flag4conformer  
	real						:: ran2, rmsd
	real(kind=8)                :: score_old, score_new
	real(kind=8)                :: binding_energy_old, binding_energy_new
	real(kind=8)				:: binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old
	real(kind=8)				:: binding_vdw_new, binding_ele_new, binding_sgb_new, binding_snp_new
	real(kind=8)				:: delta_E_vdw,delta_E_ele,delta_E_sgb,delta_E_sur
	real(kind=8)				:: delta_TS,lambda_delta_Pagg,total_abs_score_change
	real						:: entropy_old, entropy_new
	real						:: score4hydration_old, score4hydration_new
	real						:: Pagg_old, Pagg_new
	real						:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum)
	real						:: temp_entropy4individual(repeated_unit,gnum)
	character*20				:: sheet_change, accpt
	character*4					:: aminoacid_name_1
	character*4					:: group_name_1(3)
	character(len=:), allocatable   :: pep_name
	logical						:: trial_energy_available
	type(groupdetails)			:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), Tgroup(repeated_unit,gnum)
	type(groupdetails)			:: repack_group(repeated_unit,gnum)
	
    type(energyparameters), dimension(:,:), allocatable :: group_para
	integer, dimension(:), allocatable :: W_numex, W_numex4
	integer, dimension(:,:), allocatable :: W_inb, W_inb4

	call ensure_optimization_workspace()

!*********************************************************************************
!Here ip1 decides whether the sheet moves in +x, -x, +y, -y, +z and -z direction 
!ip1=1 moves sheet in x and -x
!ip1=2 moves sheet in y and -y
!ip1=3 moves sheet in z and -z
!data passing: group -> (move sheet) -> Tgroup -> (repacking) -> temp_group  
!*********************************************************************************
	flag4conformer=0
	do attempt=1, 7
		feedback_1=0
		feedback_2=0
		feedback_3=0
		trial_energy_available=.false.
		if(flag4conformer==0) then
			call ran_gen(ran2,0)
			ip1=int(ran2*3-1.0e-3)+1
			call sheet_move(ip1, group, Tgroup, ip2)
			flag4conformer=1
		elseif(flag4conformer==1) then
			call sheet_move(ip1, group, Tgroup, ip2)
		endif

		if(ip1==1) then
			if(ip2.eq.1) then
				sheet_change="   +X-axis "
			elseif(ip2.eq.2) then
				sheet_change="   -X-axis "
			endif
		elseif(ip1==2) then
			if(ip2.eq.1) then
				sheet_change="   +Y-axis "
			elseif(ip2.eq.2) then
				sheet_change="   -Y-axis "
			endif
		elseif(ip1==3) then
			if(ip2.eq.1) then
				sheet_change="   +Z-axis "
			elseif(ip2.eq.2) then
				sheet_change="   -Z-axis "
			endif
		endif
		
		call rmsd_calc(ip1, Tgroup, rmsd)
		call axis_criteria(ip1, rmsd, feedback_3)
		
		if(feedback_3==1) then 
			goto 30
		elseif(feedback_3==0) then
			accpt="Reject-rmsd"
			goto 40
		endif
30      continue

!*******************************************		
!Rotamer repacking
		temp_group=Tgroup
		temp_entropy4individual=entropy4individual
		do ic_1=1, gnum
			if(group(1,ic_1)%gtype=="ACE".or.group(1,ic_1)%gtype=="NME".or. group(1,ic_1)%gtype=="NHE") cycle
			call groupinfo(group(1,ic_1)%gtype, group_name_1, flag1)
			aminoacid_name_1=group_name_1(flag1)
			call sequence_mutation(ic_1, temp_group, temp_entropy4individual, aminoacid_name_1, &
				repack_group, Tentropy4individual, feedback_1)
			if(feedback_1==1) then 
				temp_group=repack_group
				temp_entropy4individual=Tentropy4individual
				goto 10
10				continue
			else
				goto 20
			endif
		enddo
!********************************************				
		call atom_links(temp_group, opt_numex, opt_inb, opt_numex4, opt_inb4)
		call energy_parameter(temp_group, opt_group_para)
		
        call bindingenergy(temp_group, opt_group_para, temp_entropy4individual, opt_numex, opt_inb, &
				opt_numex4, opt_inb4, score_new, binding_energy_new, binding_vdw_new, &
				binding_ele_new, binding_sgb_new, binding_snp_new, entropy_new, &
				score4hydration_new, Pagg_new)
		trial_energy_available=.true.
		delta_E_vdw=binding_vdw_new-binding_vdw_old
		delta_E_ele=binding_ele_new-binding_ele_old
		delta_E_sgb=binding_sgb_new-binding_sgb_old
		delta_E_sur=binding_snp_new-binding_snp_old
		delta_TS=real(entropy_new-entropy_old,kind=8)
		lambda_delta_Pagg=real(propensity_weighting_factor,kind=8)* &
			real((score4hydration_new+Pagg_new)-(score4hydration_old+Pagg_old),kind=8)
		total_abs_score_change=abs(delta_E_vdw)+abs(delta_E_ele)+abs(delta_E_sgb)+ &
			abs(delta_E_sur)+abs(delta_TS)+abs(lambda_delta_Pagg)
		call MC_technique_sheet(score_old, binding_energy_old, binding_vdw_old, binding_ele_old, &
			binding_sgb_old, binding_snp_old, entropy_old, score4hydration_old, Pagg_old, group, &
			entropy4individual, score_new, binding_energy_new, binding_vdw_new, binding_ele_new, &
			binding_sgb_new, binding_snp_new, entropy_new, score4hydration_new, Pagg_new, temp_group, &
			temp_entropy4individual, feedback_2)

		if(feedback_2==1) then
			flag4sheet=1
		endif
			
!If no rotamer combinations were found print previous peptide information otherwise print new peptide information	
20 		continue
		
		if(feedback_2==1) then
			accpt="Accept"
		elseif(feedback_1==1 .and. feedback_2==0) then
			accpt="Reject-MC"
		elseif(feedback_1==0) then
			accpt="Reject-Rotamer"
		endif

40 		continue
			
		if (trial_energy_available) then
			open(3, file="energydetails.txt",  access="append")
				write(3,4) step,attempt,score_new,binding_vdw_new,binding_ele_new,binding_sgb_new, &
					binding_snp_new,entropy_new, &
					(score4hydration_new+Pagg_new)*propensity_weighting_factor,delta_E_vdw, &
					delta_E_ele,delta_E_sgb,delta_E_sur,delta_TS,lambda_delta_Pagg, &
					total_abs_score_change,sheet_change,accpt
				write(3,'(*(A5))') (temp_group(1,j)%gtype,j=1,gnum)
			close(3)
		endif
	
		if(score_old.lt.energy_min(1)) then
			do i=(num4pdbsave-1),1,-1
				energy_min(i+1)=energy_min(i)
			enddo
			energy_min(1)=score_old
            
            allocate(character(len=gnum) :: pep_name)
            call convert_AA_name(group, pep_name)
			open(2, file="minimum_energy.txt", access="append")
				write(2,*) step, attempt, score_old, pep_name
			close(2)
                    deallocate(pep_name)
                    call generatepdb(step, attempt, group)
		elseif(score_old.lt.energy_min(num4pdbsave)) then
			do i=1,num4pdbsave
				if(score_old.eq.energy_min(i)) then
					goto 50
				elseif(score_old.lt.energy_min(i)) then
					do j=(num4pdbsave-1),i,-1
						energy_min(j+1)=energy_min(j)
					enddo
					energy_min(i)=score_old
					goto 60
				endif
			enddo
60			continue
			call generatepdb(step, attempt, group)
50			continue
        endif
        !if (accpt=="Accept") exit
	enddo
	
4	format(i7,i7,14f15.4,a15,a15)
	return
	end subroutine sheet_optimization


end module optimization_techniques							

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program ProteinDesign

	use constant
	use randomgenerator
	use input
	use pdbfile
	use mathfunction
	use database
	use energy_calculation
	use advancedfunction
	use transplant
	use optimization_techniques

	implicit none
	integer, parameter				:: max_scmf_attempts=3
	integer							:: step, sub_circle, i, j, k
	integer							:: stage_i, previous_stage_steps, local_stage_step
	integer							:: mc_step_count, scmf_attempt, scmf_feedback, failed_scmf_site
	integer, allocatable			:: stage_steps(:)
	integer(kind=8)					:: clk_start, clk_end, clk_rate
	integer(kind=8)					:: init_elapsed, mc_elapsed
	real							:: ran2, rmsd_x, rmsd_y, rmsd_z, rmsd
	real(kind=8)                    :: score_old, binding_energy_old
	real(kind=8)					:: binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old
    real							:: entropy_old, score4hydration_old, Pagg_old
	real							:: anneal_factor, dekt_seq, dekt_sheet
	real(kind=8)					:: init_seconds, mc_seconds
    !real							:: entropy4individual(repeated_unit,gnum)
	!type(groupdetails)				:: group(repeated_unit,gnum)
    character(len=:), allocatable   :: pep_name
    real, dimension(:,:), allocatable :: entropy4individual
	type(groupdetails), dimension(:,:), allocatable :: group
	type(groupdetails), dimension(:,:), allocatable :: temp_group
	type(groupdetails), dimension(:,:), allocatable :: scmf_start_group
	type(energyparameters), dimension(:,:), allocatable :: group_para
	integer, dimension(:), allocatable :: W_numex, W_numex4
	integer, dimension(:,:), allocatable :: W_inb, W_inb4
	character*4						:: aminoacid_name
	character(len=4)				:: failed_scmf_aminoacid
    character(len=100)				:: os_name
	character(len=32)				:: init_time,mc_time
	logical :: pdbfiles_exists
    
	call system_clock(count_rate=clk_rate)
    call system_clock(count=clk_start)
	if (clk_rate <= 0) stop 'ERROR: system clock is unavailable.'
	init_seconds=0.0d0
	mc_seconds=0.0d0
    
	
    call check_help
    call banner
	call inputfile
    call setparameter
    call readdir
	
    allocate(character(len=gnum) :: pep_name)
    allocate(entropy4individual(repeated_unit,gnum))
    allocate(group(repeated_unit,gnum))
    call init_random_seed()
    ! Create a pdbfile folder
    call get_environment_variable("OS",os_name)
    if (index(trim(os_name), "Windows_NT") > 0) then
        inquire(file='pdbfiles',exist=pdbfiles_exists)
		if (.not.pdbfiles_exists) call execute_command_line('mkdir pdbfiles')
    else
        call execute_command_line('mkdir -p pdbfiles')
	endif
	call readpdb(group) ! load PDB files
    call rotamerlib

	if(recalcu_switch==0) then
		allocate(temp_group(repeated_unit,gnum))
		allocate(scmf_start_group(repeated_unit,gnum))
		call replace_single_site_aa(group)
		scmf_start_group=group
		failed_scmf_site=0
		failed_scmf_aminoacid=''
		scmf_feedback=0
		do scmf_attempt=1,max_scmf_attempts
			write(*,'(A,I0,A,I0,A)') 'SCMF initialization attempt ',scmf_attempt, &
				' of ',max_scmf_attempts,'.'
			call scmf_substitution(scmf_start_group,sub_circle,temp_group, &
				failed_scmf_site,failed_scmf_aminoacid)
			call scmf_loop(temp_group,scmf_feedback,failed_scmf_site,failed_scmf_aminoacid)
			if (scmf_feedback == 1) exit
			if (scmf_attempt < max_scmf_attempts) then
				write(*,'(A,I0,3A)') 'SCMF will propose another sequence that preferably avoids ', &
					failed_scmf_site,':',trim(failed_scmf_aminoacid),'.'
			endif
		enddo
		if (scmf_feedback /= 1) then
			write(*,'(A,I0,A)') 'ERROR: SCMF failed to replace amino acids after ',max_scmf_attempts,' attempts.'
			write(*,'(A)') 'Suggestion: use an MD simulation to remove unfavorable atomic contacts.'
			stop
		endif
		group=temp_group
		deallocate(temp_group,scmf_start_group)
		!write(*,*) 'CHECK: scmf_loop finished'
		!flush(6)
        
		do step=1, 5
			!write(*,'(A,I0)') 'CHECK: calling nonthermal optimization ',step
			!flush(6)
			call sequence_optimization_nonthermal(group)
			!write(*,'(A,I0)') 'CHECK: nonthermal optimization finished ',step
			!flush(6)
		enddo

		!write(*,*) 'CHECK: calling initial generatepdb'
		!flush(6)
		call generatepdb(0, 0, group)
		!write(*,*) 'CHECK: initial generatepdb finished'
		!flush(6)
		energy_min=0.0
		
		allocate(group_para(repeated_unit,gnum))
		allocate(W_numex(repeated_unit*atom_num)); allocate(W_numex4(repeated_unit*atom_num))
		allocate(W_inb(repeated_unit*atom_num,20)); allocate(W_inb4(repeated_unit*atom_num,60))

		call energy_parameter(group, group_para)
		call atom_links(group, W_numex, W_inb, W_numex4, W_inb4)
	
		call initial_entropy_individual(group, entropy4individual)
		call bindingenergy(group, group_para, entropy4individual, W_numex, W_inb, W_numex4, &
				W_inb4, score_old, binding_energy_old, binding_vdw_old, binding_ele_old, &
				binding_sgb_old, binding_snp_old, entropy_old, score4hydration_old, Pagg_old)
		deallocate(W_inb); deallocate(W_inb4)
		deallocate(W_numex); deallocate(W_numex4)
		deallocate(group_para)
	
        call convert_AA_name(group, pep_name)
		open(5, file="energyprofile.txt", access="append")
            write(5,6) 0," ", pep_name, score_old,(binding_energy_old-entropy_old), &
					(score4hydration_old+Pagg_old)*propensity_weighting_factor, 0.0
		close(5)
        
		if(score_old.lt.0) energy_min(1)=score_old		
        
        flag4sheet=0 ! sheet move flag, 0 = enters judgement of sheet move, >0, need to reach a setpoint to reset
	else
		open(5, file="backup4backbone.txt", status="old")
			do i=1, num4pdbsave
				read(5,*) energy_min(i)
			enddo
			read(5,*) score_old
			read(5,*) binding_energy_old
			read(5,*) entropy_old
			read(5,*) score4hydration_old
			read(5,*) Pagg_old
            read(5,*) flag4sheet
			do i=1, repeated_unit
				do j=1, gnum
					read(5,*) entropy4individual(i,j)
				enddo
            enddo
            read(5,*) nstep_start                  ! read in the last finished step 
            nstep_start = nstep_start + 1           ! the next step is the starting step
            write(*,*) "restart at", nstep_start, "step"
		close(5)

		allocate(group_para(repeated_unit,gnum))
		allocate(W_numex(repeated_unit*atom_num)); allocate(W_numex4(repeated_unit*atom_num))
		allocate(W_inb(repeated_unit*atom_num,20)); allocate(W_inb4(repeated_unit*atom_num,60))
		call energy_parameter(group, group_para)
		call atom_links(group, W_numex, W_inb, W_numex4, W_inb4)
		call bindingenergy(group, group_para, entropy4individual, W_numex, W_inb, W_numex4, &
				W_inb4, score_old, binding_energy_old, binding_vdw_old, binding_ele_old, &
				binding_sgb_old, binding_snp_old, entropy_old, score4hydration_old, Pagg_old)
		deallocate(W_inb); deallocate(W_inb4)
		deallocate(W_numex); deallocate(W_numex4)
		deallocate(group_para)
    endif

    call system_clock(count=clk_end)
	init_seconds=real(clk_end-clk_start,kind=8)/real(clk_rate,kind=8)
	init_elapsed=nint(init_seconds,kind=8)
	write(init_time,'(I20.2,":",I2.2,":",I2.2)') &
		init_elapsed/3600,mod(init_elapsed,3600)/60,mod(init_elapsed,60)
	init_time=adjustl(init_time)
	write(*,'(1X,A,T45,"= ",A)') 'Initialization time spent',trim(init_time)
    
    !!!! 2026/07/15/ Prepare for annealing stages
	if (anneal_stages > 0) then
		allocate(stage_steps(anneal_stages))
		stage_steps=nstep_terminal/anneal_stages
		stage_steps(anneal_stages)=stage_steps(anneal_stages)+mod(nstep_terminal,anneal_stages)
		dekt_seq=ekt_seq_high-ekt_seq_low
		dekt_sheet=ekt_sheet_high-ekt_sheet_low
	endif
        
	!!!!Main loop here: Sequence_optimization and Sheet_optimization
	call system_clock(count=clk_start)
	do step=nstep_start, nstep_terminal
		if (anneal_stages > 0) then
			previous_stage_steps=0
			do stage_i=1,anneal_stages ! find current stage
				if (step <= previous_stage_steps+stage_steps(stage_i)) exit
				previous_stage_steps=previous_stage_steps+stage_steps(stage_i)
			enddo
			local_stage_step=step-previous_stage_steps ! find local stages
			if (stage_steps(stage_i) > 1) then
				anneal_factor=real(stage_steps(stage_i)-local_stage_step)/real(stage_steps(stage_i)-1)
			else
				anneal_factor=0.0
			endif
			ekt_seq=ekt_seq_low+anneal_factor*dekt_seq
			ekt_sheet=ekt_sheet_low+anneal_factor*dekt_sheet
		endif

		if(step.le.steps_before_sheetmove) then
			call sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, &
					binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
					score4hydration_old, Pagg_old)
		else
			if(flag4sheet==0) then
				call ran_gen(ran2,0)
				if(ran2.ge.sheet_switch) then
					call sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, &
							binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
							score4hydration_old, Pagg_old)
				else
					! The CONROT backbone-angle move is disabled.
					call sheet_optimization(group, step, entropy4individual, score_old, binding_energy_old, &
							binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
							score4hydration_old, Pagg_old)
				endif
			else
				call sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, &
						binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
						score4hydration_old, Pagg_old)
				flag4sheet=flag4sheet+1
				if(flag4sheet.gt.sheetmove_interval) then
					flag4sheet=0
				endif
			endif

			! Previous CONROT test driver retained below as comments only.
			!if(flag4sheet==0) then
				!call ran_gen(ran2,0)
				!if(ran2.ge.sheet_switch) then
    !                ekt=ekt_seq
				!	call sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, &
				!			binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
				!			score4hydration_old, Pagg_old)
    !            else ! sheet_optimization sets flag4sheet to 1.
					! Backbone-angle CONROT perturbation is disabled.  The v1.41 implementation is
					! retained in module optimization_techniques for later debugging, but its driver
					! is commented out here so it cannot modify the peptide backbone.
                    !call ran_gen(ran2,0)
                    !CONROT_flag=1
					!if(CONROT_flag.eq.1.and.ran2.le.CONROT_switch) then
					!	call backbone_optimization_conrot(group, step, entropy4individual, score_old, binding_energy_old, &
					!			binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
					!			score4hydration_old, Pagg_old)
					!else
					!	call sheet_optimization(group, step, entropy4individual, score_old, binding_energy_old, &
					!			binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
					!			score4hydration_old, Pagg_old)
					!endif
				!endif
    !        else
    !            !ekt=ekt_seq * (sheetmove_interval - flag4sheet) / sheetmove_interval * 1.5 + 0.5
    !            ekt=ekt_seq
				!call sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, &
				!		binding_vdw_old, binding_ele_old, binding_sgb_old, binding_snp_old, entropy_old, &
				!		score4hydration_old, Pagg_old)
				!flag4sheet=flag4sheet+1
				!if(flag4sheet.gt. sheetmove_interval) then
				!	flag4sheet=0
    !            end if
    !        endif				
        endif
        
        call rmsd_calc(4, group, rmsd)
        
        call convert_AA_name(group, pep_name)
		open(5, file="energyprofile.txt", access="append")
            write(5,6) step," ",pep_name, score_old,(binding_energy_old-entropy_old), &
					(score4hydration_old+Pagg_old)*propensity_weighting_factor, rmsd
		close(5)
        

		call generatebackup4pdb(group) ! The start pdb for next step or final pdb for current step
		call ran_gen(ran2,1)	
		open(5, file="backup4backbone.txt", status="replace")
			do i=1, num4pdbsave
				write(5,*) energy_min(i)
			enddo
			write(5,*) score_old
			write(5,*) binding_energy_old
			write(5,*) entropy_old
			write(5,*) score4hydration_old
			write(5,*) pagg_old
            write(5,*) flag4sheet
			do i=1, repeated_unit
				do j=1, gnum
					write(5,*) entropy4individual(i,j)
				enddo
			enddo
			write(5,*) step
		close(5)

    enddo
    
    write(*,'(A)') '========================================================================'
	call system_clock(count=clk_end)
	mc_seconds=real(clk_end-clk_start,kind=8)/real(clk_rate,kind=8)
	mc_step_count=max(0,nstep_terminal-nstep_start+1)
	mc_elapsed=nint(mc_seconds,kind=8)
	write(mc_time,'(I20.2,":",I2.2,":",I2.2)') &
		mc_elapsed/3600,mod(mc_elapsed,3600)/60,mod(mc_elapsed,60)
	mc_time=adjustl(mc_time)
	write(*,'(1X,A,I0,A,T45,"= ",A)') 'Search time spent on ', &
		mc_step_count,' MC steps ',trim(mc_time)

	if (mc_step_count > 0) then
		mc_elapsed=nint(mc_seconds/real(mc_step_count,kind=8),kind=8)
		write(mc_time,'(I20.2,":",I2.2,":",I2.2)') &
			mc_elapsed/3600,mod(mc_elapsed,3600)/60,mod(mc_elapsed,60)
		mc_time=adjustl(mc_time)
		write(*,'(1X,A,T45,"= ",A)') 'Average time spent on each step',trim(mc_time)
	else
		write(*,'(1X,A,T45,"= ",A)') 'Average time spent on each step','no MC steps performed'
	endif
    
	!write(*,'(1X,A,T45,"= ",F12.3,A)') 'Total measured time',initialization_seconds+mc_seconds,' s'
	write(*,'(A)') '========================================================================'
6	format(i5,2a,4f10.4)
7	format(i5,2a,5f10.4)

end
	
