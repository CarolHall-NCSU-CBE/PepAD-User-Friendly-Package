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
		integer			    :: num4peptides, peptideID(8) ! need to use dynamic array for this, or assign big value
    end type
	
    type hydration_types      ! new data type, contains hydration properties
        integer				:: hnum, pnum, cnum, onum
        character*4			:: hgtype(10), pgtype(10), cgtype(10), ogtype(10)
    end type
    
    type :: AminoRestriction
		character*4 :: gtype
		integer		:: max
        integer		:: count
        integer		:: clas
	end type AminoRestriction

	type(AminoRestriction)  :: aa_restrictions(10)
	integer					:: n_restrictions		! number of restriction of AA
    
    integer					:: gnum, repeated_unit, num4category				
	integer					:: atom_num

	integer, parameter		:: num4pdbsave=10
	integer, parameter		:: maximum_nmr_site_num=5
	integer, parameter		:: maximum_void_site_num=30
	integer					:: void_site_start(maximum_void_site_num), void_site_end(maximum_void_site_num)
	integer					:: void_site_num
	integer					:: nmr_site_num, nmr_site_ID(maximum_nmr_site_num)
    integer					:: restraint_num, exclusion_num, NMR_pool_size
	integer					:: nstep_start, nstep_terminal, idum, recalcu_switch, flag4sheet
	integer					:: fpho,fpol,fchg,foth
	integer					:: sitenum4mutation
    integer, allocatable    :: sitenum4mutation_group(:)
	integer					:: sheetmove_interval
    integer					:: seed_switch, generate_traj_switch
    integer					:: ASN_limit, ASN_count

	real, parameter			:: scmfswitch1=0.6
!	real, parameter			:: scmfswitch2=0.8
	real, parameter			:: dihedral_weighting_factor=0.20
	!real, parameter			:: propensity_weighting_factor=3.0 !数值为固定值
    real					:: propensity_weighting_factor
	real, parameter			:: vdw14_coeff=2.0  !1-4 VDW scale factor
	real, parameter			:: ele14_coeff=1.2  !1-4 ELE scale factor

	real					:: ekt, ekt_sequence
	real					:: energy_min(num4pdbsave)
	real					:: sheet_switch, ekt_sheet, rmsd_max, rmsd_max_x, rmsd_max_y, rmsd_max_z
	real					:: displacement_factor, displacement_factor_x, displacement_factor_y, displacement_factor_z
	!type RES4chain
	!	integer				:: num
	!	integer				:: IDs(gnum)
	!end type

	type(lib4aminoacid)		:: aa_lib(23)
	type(hydration_types)   :: hydrationprop
	type(peptideassembly), allocatable 	:: selfassembly(:)
	character*50			:: filename							! The input PDB file
	character*4             :: NMR_restraint_AA(10), exclusion_AA(10), NMR_AA_pool(10)
    character(len=:), allocatable :: mydir						! The exe file located directory
    type(groupdetails),  allocatable  :: original_group(:,:)

end module constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module randomgenerator
	use constant
    contains
! Generate the random number
    
	subroutine init_random_seed()
		integer ::  iiii, kkk, clock
		integer, dimension(:), allocatable :: seed
		real :: ran
   
		call random_seed(size = kkk)
		allocate(seed(kkk))   
		call system_clock(count=clock)
        seed = clock + 37 * (/ (iiii - 1, iiii = 1, kkk) /)
    
        if(seed_switch == 1) then
			seed(1) = 4390990 + idum
			seed(2) = 5570636 + 0.5 * idum
		endif
        
		call random_seed(put = seed)
		deallocate(seed)
		call random_number(ran) !discard first random value
    
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
module input

	use constant

	contains
	subroutine inputfile
	implicit none
	integer						:: i,j,pos,iostat
    character(len=100)          :: void_site_input, nmr_site_input
    character(len=100)          :: restraints_input, exception_input, restriction_inputs, NMR_res_pool_input
    n_restrictions=0			! set number of restrictions, this limit maximum number of certain AA
	open(10, file="input.txt", status="old")
		read(10,*)
        read(10,*) gnum, repeated_unit, num4category
		read(10,*)
		read(10,*) filename, recalcu_switch
		read(10,*)
		read(10,*) nstep_start, nstep_terminal
		read(10,*)
		read(10,*) seed_switch, idum, ekt_sequence
		read(10,*)
		read(10,*) sheet_switch, ekt_sheet, sheetmove_interval
        read(10,*)
		read(10,*) rmsd_max_x, rmsd_max_y, rmsd_max_z, displacement_factor_x, displacement_factor_y, displacement_factor_z
        read(10,*)
        read(10,*) propensity_weighting_factor
		read(10,*)
		read(10,*) fpho,fpol,fchg,foth
		read(10,*)
		read(10,*)
        
        allocate (selfassembly(num4category))
        
		do i=1, num4category
			read(10,*) selfassembly(i)%num4peptides
			read(10,*) (selfassembly(i)%peptideID(j),j=1,selfassembly(i)%num4peptides)
        enddo
        
 		read(10,*) 
        read(10,'(A)') nmr_site_input
        call read_nmr_site_input(trim(nmr_site_input), ',')
        
        !*******************output the info*********
		write(*, "(A,1X,I0)") "Number of NMR sites:", nmr_site_num
		write(*, "(A, *(I0, 1X))") "Residue ID of the NMR sites:", (nmr_site_ID(i), i = 1, nmr_site_num)
		!*******************output the info*********
        
  !      read(10,*) nmr_site_num
		!read(10,*) (nmr_site_ID(j),j=1,nmr_site_num)
		read(10,*)
        read(10,'(A)') void_site_input							! add "(a)" in the read to read entire line
        call read_void_site_input(trim(void_site_input), ',')
        
		read(10,*)
        read(10,'(A)') restraints_input        
        call read_restraints_input(trim(restraints_input), 1)
        
 		!read(10,*)
   !     read(10,'(A)') exception_input
   !     call read_restraints_input(trim(exception_input), 2)       
        
 		read(10,*)
        read(10,'(A)') NMR_res_pool_input
        call read_restraints_input(trim(NMR_res_pool_input), 3)          
        
 		read(10,*)
   !     read(10,*) ASN_limit
        
        do
			read(10, '(A)', iostat=iostat) restriction_inputs
			if (iostat /= 0) exit  ! End of file
			if (trim(restriction_inputs) == '') exit   ! Blank line ends block
			if (trim(restriction_inputs)=='NONE'.or.trim(restriction_inputs)=='None'.or.trim(restriction_inputs)=='none') exit   ! "None" ends block
			if (restriction_inputs(1:2) == '**') exit  ! New section starts
			n_restrictions = n_restrictions + 1
            read(restriction_inputs, *) aa_restrictions(n_restrictions)%gtype, aa_restrictions(n_restrictions)%max
        end do
        
        !*******************output restrictions*********
		write(*, "(A,1X,I0)") "Number of restrictions:", n_restrictions
		do i = 1, n_restrictions
			write(*, "(A3, ': ', I0)") aa_restrictions(i)%gtype, aa_restrictions(i)%max
		end do
		!***********************************************
        
	close(10)
	return
    end subroutine inputfile

    subroutine setparameter
    integer :: i
    
    atom_num=60*gnum
    allocate (sitenum4mutation_group(gnum))
    allocate (original_group(repeated_unit, gnum))

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
    
    !do i = 1, exclusion_num
    !    do j = 1, hydrationprop%hnum
    !        if (hydrationprop%hgtype(j) == exclusion_AA(i)) then
    !            ! Shift left to remove the item
    !            do k = j, hydrationprop%hnum - 1
    !                hydrationprop%hgtype(k) = hydrationprop%hgtype(k+1)
    !            end do
    !            hydrationprop%hgtype(hydrationprop%hnum) = "NONE"  ! Clear last entry
    !            hydrationprop%hnum = hydrationprop%hnum - 1  ! Decrease count
    !            exit
    !        end if
    !    end do
    !end do
    
    !do i = 1, exclusion_num
    !    do j = 1, hydrationprop%hnum
    !        if (hydrationprop%hgtype(j) == exclusion_AA(i)) then
    !            ! Shift left to remove the item
    !            do k = j, hydrationprop%hnum - 1
    !                hydrationprop%hgtype(k) = hydrationprop%hgtype(k+1)
    !            end do
    !            hydrationprop%hgtype(hydrationprop%hnum) = "NONE"  ! Clear last entry
    !            hydrationprop%hnum = hydrationprop%hnum - 1  ! Decrease count
    !            exit
    !        end if
    !    end do
    !end do
    !
    !! Remove from pgtype
    !do i = 1, exclusion_num
    !    do j = 1, hydrationprop%pnum
    !        if (hydrationprop%pgtype(j) == exclusion_AA(i)) then
    !            do k = j, hydrationprop%pnum - 1
    !                hydrationprop%pgtype(k) = hydrationprop%pgtype(k+1)
    !            end do
    !            hydrationprop%pgtype(hydrationprop%pnum) = "NONE"
    !            hydrationprop%pnum = hydrationprop%pnum - 1
    !            exit
    !        end if
    !    end do
    !end do
    !
    !! Remove from cgtype
    !do i = 1, exclusion_num
    !    do j = 1, hydrationprop%cnum
    !        if (hydrationprop%cgtype(j) == exclusion_AA(i)) then
    !            do k = j, hydrationprop%cnum - 1
    !                hydrationprop%cgtype(k) = hydrationprop%cgtype(k+1)
    !            end do
    !            hydrationprop%cgtype(hydrationprop%cnum) = "NONE"
    !            hydrationprop%cnum = hydrationprop%cnum - 1
    !            exit
    !        end if
    !    end do
    !end do
    !
    !! Remove from ogtype
    !do i = 1, exclusion_num
    !    do j = 1, hydrationprop%onum
    !        if (hydrationprop%ogtype(j) == exclusion_AA(i)) then
    !            do k = j, hydrationprop%onum - 1
    !                hydrationprop%ogtype(k) = hydrationprop%ogtype(k+1)
    !            end do
    !            hydrationprop%ogtype(hydrationprop%onum) = "NONE"
    !            hydrationprop%onum = hydrationprop%onum - 1
    !            exit
    !        end if
    !    end do
    !end do
    !
    !write(*,*) hydrationprop
    
    
    end subroutine
    
    subroutine read_nmr_site_input(input, delim)
		character(len=*), intent(in) :: input, delim
		character(len=20), dimension(:), allocatable :: segments
		character(len=100) :: temp
		integer :: start, ending, len_input, len_delim, seg_index, site_num, i, j
		character(len=20), allocatable :: temp_segments(:)
        character(len=20) :: segment

		len_input = len_trim(input)
		len_delim = len_trim(delim)
		allocate(segments(len_input / 2 + 1))  ! Over-allocate initially

		temp = trim(input)
		start = 1
		site_num = 0	 ! residue "1, 4-6, 14" = three sites
		nmr_site_num = 0 ! one residue = one NMR site
        
		do while (start <= len_input)
			ending = index(temp(start:), delim)
			if (ending == 0) then
				site_num = site_num + 1
				segments(site_num) = temp(start:)
				exit
			else
				site_num = site_num + 1
				segments(site_num) = temp(start:start + ending - 2)
				start = start + ending - 1 + len_delim
			end if
        end do

		! Use temporary allocatable array
		temp_segments = segments(1:site_num)
		call move_alloc(temp_segments, segments)		! Reassign resized array
        
        do i = 1, site_num
			segment = trim(segments(i))
			pos = index(segment, '-')
			if (pos > 0) then
				read(segment(1:pos-1), *) site_start
				read(segment(pos+1:), *) site_ending
                do j = site_start, site_ending
                    nmr_site_num = nmr_site_num + 1
                    nmr_site_ID(nmr_site_num) = j		! j records residues that are NMR sites
                enddo
            else
                nmr_site_num = nmr_site_num + 1
				read(segment, *)  nmr_site_ID(nmr_site_num)
			end if
        enddo
        return
    end subroutine read_nmr_site_input
    
    subroutine read_void_site_input(input, delim)
		character(len=*), intent(in) :: input, delim
		character(len=20), dimension(:), allocatable :: segments
		character(len=100) :: temp
		integer :: start, ending, len_input, len_delim, seg_index, i
		character(len=20), allocatable :: temp_segments(:)
        character(len=20) :: segment

		len_input = len_trim(input)
		len_delim = len_trim(delim)
		allocate(segments(len_input / 2 + 1))  ! Over-allocate initially

		temp = trim(input)
		start = 1
		void_site_num = 0

		do while (start <= len_input)
			ending = index(temp(start:), delim)
			if (ending == 0) then
				void_site_num = void_site_num + 1
				segments(void_site_num) = temp(start:)
				exit
			else
				void_site_num = void_site_num + 1
				segments(void_site_num) = temp(start:start + ending - 2)
				start = start + ending - 1 + len_delim
			end if
        end do

        allocate(temp_segments(void_site_num))
		temp_segments = segments(1:void_site_num)
		call move_alloc(temp_segments, segments)  ! Reassign resized array
        
        do i = 1, void_site_num
			segment = trim(segments(i))
			pos = index(segment, '-')
			if (pos > 0) then
				read(segment(1:pos-1), *) void_site_start(i)
				read(segment(pos+1:), *) void_site_end(i)
			else
				read(segment, *)  void_site_start(i)
				void_site_end(i) = void_site_start(i)
			end if
        enddo
        return
    end subroutine read_void_site_input
    
    subroutine read_restraints_input(input, signal)
		character(len=*), intent(in) :: input
		character(len=100)           :: temp
        character*4                  :: AA_list(10)
		integer                      :: i, res_num, signal

        res_num = 0
        temp = input
           
        do while (len_trim(temp) > 0)
            i = index(temp, ',')   
            if (i == 0) then
				if (len_trim(temp) > 0) then                        ! avoid spaces
					res_num = res_num + 1
					AA_list(res_num) = trim(adjustl(temp))
                endif
                exit
            else 
                if (len_trim(temp(:i-1)) > 0) then                  ! avoid spaces
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
				do i=1, res_num
					NMR_AA_pool(i) = AA_list(i)
                enddo
                NMR_pool_size = res_num
				write(*,*) "NMR constrained AA pool: ", NMR_AA_pool
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
end module input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pdbfile

	use	constant

	contains
! Read an initial file.
	subroutine readpdb(group)
	implicit none
	integer						:: anum, status, chainID, ic, i, flag
	real						:: x, y, z
	character*4					:: char, atype, name
	type(groupdetails)			:: group(repeated_unit,gnum)
	character*150               :: line
    
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
		if(ic.gt.(chainID*gnum)) chainID=chainID+1 !判断peptide 序号
		ic=ic-(chainID-1)*gnum						!判断当前aa 在peptide中的序号
        
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

		if(chainID==1) then
			if(ic.ne.flag) then
				flag=ic
				do i=1, void_site_num
					if(ic.ge.void_site_start(i).and.ic.le.void_site_end(i)) then
						goto 10
					endif
				enddo
				sitenum4mutation=sitenum4mutation+1				! record how many sites for mutation 
				sitenum4mutation_group(sitenum4mutation)=ic     ! record which AA (第几个) in the peptide for mutation 
			endif
		endif
10		continue
		!write(*,*) chainID, ic, group(chainID,ic)%gtype
	enddo
	
	close(10)
	
	original_group=group
	
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
	open(10,file='pdbfiles/'//trim(adjustl(stepchar))//'_'//trim(adjustl(attemptchar))//'.pdb', access="append")
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
    
    

end module pdbfile


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
	
	subroutine rotamerlib
	implicit none
	integer							:: status, grade, rotanum, anum, num, i, j, k
	real							:: x, y, z
	character*4						:: char, atype, name

	open(10, file=trim(mydir)//"lib/rotamer", status="old")
		do while(.true.)
			read(10, "(a,i3,i3)", iostat=status) name, grade, rotanum		
			if(status.ne.0) exit
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
						read(10,"(<grade>(i5))") (aa_lib(i)%dihedralangle(j,k), k=1, grade)
					enddo
				endif
			endif
		enddo
	close(10)

	aa_lib%cnum1=0
	aa_lib%cnum2=0
	aa_lib%cnum3=0

	do i=1, 20
		open (10, file=trim(mydir)//'lib/RotamerLibrary/'//trim(aa_lib(i)%gtype), status="old")	
		do while(.true.)
			read(10, *, iostat=status) char, anum, atype, name, char, num, x, y, z
			if(status.ne.0) exit
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
	enddo
	
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
	integer							:: i, j, k, status, atomid
	real							:: charge, epsion, r, rborn, fs, dielecons
	character*4						:: lbres, igraph
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(energyparameters)			:: group_para(repeated_unit,gnum)
    
	! declaring Tgroup_para as a two-dimensional allocatable array of type (energy parameters). T指临时
	type(energyparameters), dimension(:,:), allocatable :: Tgroup_para	

	allocate(Tgroup_para(repeated_unit,gnum))
    
	do i=1, repeated_unit
		do j=1, gnum
			open(10, file=trim(mydir)//'lib/ForceField/'//trim(adjustl(trim(group(i,j)%gtype))), status="old")
			read(10, *)
			do while(.true.)
				read(10, 20, iostat=status) lbres, igraph, charge, epsion, r, rborn, fs, dielecons, atomid
				if(status.ne.0) goto 30
				do k=1, group(i,j)%cnum1
					if(group(i,j)%atype1(k)==igraph) then
						Tgroup_para(i,j)%charge1(k)=charge
						Tgroup_para(i,j)%epsion1(k)=epsion
						Tgroup_para(i,j)%r1(k)=r
						Tgroup_para(i,j)%rborn1(k)=rborn
						Tgroup_para(i,j)%fs1(k)=fs
						Tgroup_para(i,j)%dielecons1(k)=dielecons
						Tgroup_para(i,j)%atomid1(k)=atomid                    
                        
						goto 40
					endif
				enddo
				do k=1, group(i,j)%cnum2
					if(group(i,j)%atype2(k)==igraph) then
						Tgroup_para(i,j)%charge2(k)=charge
						Tgroup_para(i,j)%epsion2(k)=epsion
						Tgroup_para(i,j)%r2(k)=r
						Tgroup_para(i,j)%rborn2(k)=rborn
						Tgroup_para(i,j)%fs2(k)=fs
						Tgroup_para(i,j)%dielecons2(k)=dielecons
						Tgroup_para(i,j)%atomid2(k)=atomid
						goto 40
					endif
				enddo
				do k=1, group(i,j)%cnum3
					if(group(i,j)%atype3(k)==igraph) then
						Tgroup_para(i,j)%charge3(k)=charge
						Tgroup_para(i,j)%epsion3(k)=epsion
						Tgroup_para(i,j)%r3(k)=r
						Tgroup_para(i,j)%rborn3(k)=rborn
						Tgroup_para(i,j)%fs3(k)=fs
						Tgroup_para(i,j)%dielecons3(k)=dielecons
						Tgroup_para(i,j)%atomid3(k)=atomid
						goto 40
					endif
				enddo
40				continue
			enddo
30			continue
			close(10)

20	format(2a4, 6e16.8, i8)
		enddo
    enddo
    
    


            
	do i=1, repeated_unit
		do j=1, gnum
			do k=1, group(i,j)%cnum1
                
				if(Tgroup_para(i,j)%dielecons1(k)<=0.1) then
					open(10, file="error.txt", access="append")
						write(10,*) group(i,j)%gtype, i, j, "cluster 1 has wrong force field parameter", Tgroup_para(i,j)%dielecons1(k), "in the LIB!"
						write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
					close(10)
					stop
				endif
			enddo
			do k=1, group(i,j)%cnum2
				if(Tgroup_para(i,j)%dielecons2(k)<=0.1) then
					open(10, file="error.txt", access="append")
						write(10,*) group(i,j)%gtype, i, j, "cluster 2 has wrong force field parameter", Tgroup_para(i,j)%dielecons1(k), "in the LIB!"
						write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
					close(10)
					stop
				endif
			enddo
			do k=1, group(i,j)%cnum3
				if(Tgroup_para(i,j)%dielecons3(k)<=0.1) then
					open(10, file="error.txt", access="append")
						write(10,*) group(i,j)%gtype, i, j, "cluster 3 has wrong force field parameter", Tgroup_para(i,j)%dielecons1(k), "in the LIB!"
						write(10,*) "Please check whether the atom type of PDB file matches the atom type of Force Field LIB or not!"
					close(10)
					stop
				endif
			enddo
		enddo
	enddo	
	group_para=Tgroup_para
	deallocate(Tgroup_para)
	
	return
	end subroutine energy_parameter

	
	subroutine atom_links(group, numex, inb, numex4, inb4)
	implicit none
	type atomlink      ! The Data type "atomlink" is used to store the neighboring atoms for each atom.
		integer				:: linknum
		integer				:: linkindex(4)
	end type

	integer							:: categoryID, i, ii, ic, j, k, status, i1, j1, i2, j2, atomid, natom
	integer							:: id, linknum, linkindex(4)
	integer							:: ipres, numex(repeated_unit*atom_num), inb(repeated_unit*atom_num,20), numex4(repeated_unit*atom_num), inb4(repeated_unit*atom_num,60)

	type(groupdetails)				:: group(repeated_unit,gnum)
	type(atomlink)					:: atom(repeated_unit*atom_num)									

	natom=0
	do categoryID=1, num4category
		do i=1, selfassembly(categoryID)%num4peptides ! sheet中peptides个数
			ii=selfassembly(categoryID)%peptideID(i) !peptide编号
			do ic=1, gnum
				ipres=natom
                
				open(10, file=trim(mydir)//'lib/Atomlink/'//trim(adjustl(trim(group(ii,ic)%gtype))), status="old")
                
				do while(.true.)
					read(10, 20, iostat=status) id, linknum, (linkindex(k), k=1, linknum)
					if(status.ne.0) goto 30
					atomid=ipres+id
					atom(atomid)%linknum=linknum
					do k=1, linknum
						atom(atomid)%linkindex(k)=ipres+linkindex(k)
					enddo
				enddo			
30				continue
				close(10)
				natom=atomid
			enddo
		enddo
	enddo
20  format(i6, i7, 4i3)	

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
	
	integer							:: chainID, i, ii, j, k, status, i1, j1, i2, j2, atomid, natom
	integer							:: ic, id, linknum, linkindex(4)
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	type(groupdetails)				:: group(repeated_unit,gnum)
	type(atomlink)					:: atom(50)

	ii=chainID
    
	open(10, file=trim(mydir)//'lib/Atomlink/'//trim(adjustl(trim(group(ii,ic)%gtype))), status="old")
    
	do while(.true.)
		read(10, 20, iostat=status) id, linknum, (linkindex(j), j=1, linknum)
		if(status.ne.0) goto 30
		atomid=id
		atom(atomid)%linknum=linknum
		do j=1, linknum
			atom(atomid)%linkindex(j)=linkindex(j)
		enddo
	enddo
30	continue
	close(10)
	natom=atomid
20  format(i6, i7, 4i3)

	do i1=1, natom
		S_numex(i1)=0
		do j1=1, atom(i1)%linknum
			S_numex(i1)=S_numex(i1)+1						!S_numex(i1)最后还是等于atom(i1)%linknum
			S_inb(i1,S_numex(i1))=atom(i1)%linkindex(j1)    
        enddo

		do j1=1, atom(i1)%linknum
			i2=atom(i1)%linkindex(j1)

			if(i2.gt.natom) atom(i2)%linknum=0     ! v1.20, array里面不能放-1 和0               
            
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


    subroutine find_rest_aa(group, rest_AA, l)
    implicit none
    integer						:: ii, j, k, l, i
    character*4					:: NMR_AA(6), rest_AA(6), AA
    type(groupdetails)			:: group(repeated_unit,gnum)

    ii=1
    do j = 1, nmr_site_num ! 如果当前选中的site是nmr site, 则loop所有的nmr site, 记录他们的amino acid type
        AA = group(ii,nmr_site_ID(j))%gtype
        do i=1, NMR_pool_size
			if (AA == NMR_AA_pool(i) .or. AA == "N" // NMR_AA_pool(i) .or. AA == "C" // NMR_AA_pool(i) ) then
				NMR_AA(j) = NMR_AA_pool(i)
                exit
            endif
        enddo
    enddo

    l=0
    do j = 1, NMR_pool_size
        do k = 1, nmr_site_num
            if (NMR_AA_pool(j) == NMR_AA(k)) goto 98 ! loop all NMR aa pool, if the one of the NMR AA is not in the NRM pool, pick it out
        end do
        l = l + 1
        rest_AA(l) = NMR_AA_pool(j)      
98      continue
    end do
    
    return ! only return rest_AA
    end subroutine find_rest_aa
    
    
	subroutine mc_choose_aminoacid(ic, group, aminoacid_name)
    implicit none
	integer						:: ii, ic, ip, ip1, i, j, k, l
    integer						:: res_h_avail_num, res_p_avail_num, res_c_avail_num, res_o_avail_num, clas
	real						:: ran2
	character*4					:: restrict_aa, aminoacid_name, ip2, rest_AA(6), current_AA_name 
    character*4					:: res_h_avail(10), res_p_avail(10), res_c_avail(10), res_o_avail(10)
	type(groupdetails)			:: group(repeated_unit,gnum)
	
    ii=1
    ! count current peptide sequence if they are reach limit of restricted amino acids
    do j=1, n_restrictions
		aa_restrictions(j)%count = 0
		do i=1, gnum
            restrict_aa = aa_restrictions(j)%gtype
			if (group(ii,i)%gtype=="N"//restrict_aa.or.group(ii,i)%gtype==restrict_aa.or.group(ii,i)%gtype =="C"//restrict_aa) then
				aa_restrictions(j)%count = aa_restrictions(j)%count + 1
			endif
        enddo
    enddo    
          
    do i = 1, nmr_site_num
        if (nmr_site_ID(i) == ic) then  ! looping every nmr site ID to see if it is equal to current random picked site
            call find_rest_aa(group, rest_AA, l) ! rest_aa contains nmr_AA that is not appear in the pep, "l" is the list length 
            call ran_gen(ran2,0)
            ip1=int(ran2 * l -1.0e-3)+1    ! 选择剩余aa中的第几个, l=how many AA in the rest_AA list
            ip2=rest_AA(ip1)
            
            if (ic == 1) then
                aminoacid_name = "N" // ip2
            else if (ic == gnum) then
                aminoacid_name = "C" // ip2
            else
                aminoacid_name = ip2
            end if
            goto 99						! Exit the loop earlier if a match is found
        end if
    end do
    
    current_AA_name = group(ii,ic)%gtype
    !!*** need to revise this part
    if(current_AA_name=="ALA".or.current_AA_name=="LEU".or.current_AA_name=="VAL".or.current_AA_name=="ILE".or.current_AA_name=="MET".or. &
	   current_AA_name=="PHE".or.current_AA_name=="TYR".or.current_AA_name=="TRP".or.current_AA_name=="GLY") then
		ip=10
	elseif(current_AA_name=="ASN".or.current_AA_name=="GLN".or.current_AA_name=="SER".or.current_AA_name=="THR".or.current_AA_name=="HIE") then
		ip=20
	elseif(current_AA_name=="ARG".or.current_AA_name=="LYS".or.current_AA_name=="GLU".or.current_AA_name=="ASP") then
		ip=30
	elseif(current_AA_name=="PRO".or.current_AA_name=="CYS") then
		ip=40
	elseif(current_AA_name=="NALA".or.current_AA_name=="NLEU".or.current_AA_name=="NVAL".or.current_AA_name=="NILE".or.current_AA_name=="NMET".or. &
		   current_AA_name=="NPHE".or.current_AA_name=="NTYR".or.current_AA_name=="NTRP".or.current_AA_name=="NGLY") then
		ip=11
	elseif(current_AA_name=="NASN".or.current_AA_name=="NGLN".or.current_AA_name=="NSER".or.current_AA_name=="NTHR".or.current_AA_name=="NHIE") then
		ip=21
	elseif(current_AA_name=="NARG".or.current_AA_name=="NLYS".or.current_AA_name=="NGLU".or.current_AA_name=="NASP") then
		ip=31
	elseif(current_AA_name=="NPRO".or.current_AA_name=="NCYS") then
		ip=41
	elseif(current_AA_name=="CALA".or.current_AA_name=="CLEU".or.current_AA_name=="CVAL".or.current_AA_name=="CILE".or.current_AA_name=="CMET".or. &
		   current_AA_name=="CPHE".or.current_AA_name=="CTYR".or.current_AA_name=="CTRP".or.current_AA_name=="CGLY") then
		ip=12
	elseif(current_AA_name=="CASN".or.current_AA_name=="CGLN".or.current_AA_name=="CSER".or.current_AA_name=="CTHR".or.current_AA_name=="CHIE") then
		ip=22
	elseif(current_AA_name=="CARG".or.current_AA_name=="CLYS".or.current_AA_name=="CGLU".or.current_AA_name=="CASP") then
		ip=32
	elseif(current_AA_name=="CPRO".or.current_AA_name=="CCYS") then
		ip=42
	endif

40  continue
    
    !!!!!!!!!!! Modify available AA !!!!!!!!!!!!!!!!
    ! copy the original hydration property
    res_h_avail = hydrationprop%hgtype; res_p_avail = hydrationprop%pgtype
    res_c_avail = hydrationprop%cgtype; res_o_avail = hydrationprop%ogtype
	res_h_avail_num = hydrationprop%hnum; res_p_avail_num = hydrationprop%pnum
    res_c_avail_num = hydrationprop%cnum; res_o_avail_num = hydrationprop%onum
    do i = 1, n_restrictions
        clas = aa_restrictions(i)%clas
        restrict_aa = aa_restrictions(i)%gtype
        
        if(aa_restrictions(i)%count >= aa_restrictions(i)%max) then
        select case(clas)
            case (1)
				do j = 1, res_h_avail_num
					if (res_h_avail(j) == restrict_aa) then         ! drop the restricted aa in the list
						do k = j, res_h_avail_num - 1				! Shift left to remove the item
							res_h_avail(k) = res_h_avail(k+1)
						end do
						res_h_avail(res_h_avail_num) = "NONE"		! Clear last entry
						res_h_avail_num = res_h_avail_num - 1		! Decrease count
						exit
					end if
				end do
            case (2)
                do j = 1, res_p_avail_num
					if (res_p_avail(j) == restrict_aa) then         ! drop the restricted aa in the list
						do k = j, res_p_avail_num - 1               ! Shift left to remove the item
							res_p_avail(k) = res_p_avail(k+1)
						end do
						res_p_avail(res_p_avail_num) = "NONE"       ! Clear last entry
						res_p_avail_num = res_p_avail_num - 1		! Decrease count
						exit
					end if
				end do
            case (3)
                do j = 1, res_c_avail_num
					if (res_c_avail(j) == restrict_aa) then         ! drop the restricted aa in the list
						do k = j, res_c_avail_num - 1               ! Shift left to remove the item
							res_c_avail(k) = res_c_avail(k+1)
						end do
						res_c_avail(res_c_avail_num) = "NONE"       ! Clear last entry
						res_c_avail_num = res_c_avail_num - 1		! Decrease count
						exit
					end if
				end do
            case (4)
                do j = 1, res_o_avail_num
					if (res_o_avail(j) == restrict_aa) then         ! drop the restricted aa in the list
						do k = j, res_o_avail_num - 1					! Shift left to remove the item
							res_o_avail(k) = res_o_avail(k+1)
						end do
						res_o_avail(res_o_avail_num) = "NONE"       ! Clear last entry
						res_o_avail_num = res_o_avail_num - 1		! Decrease count
						exit
					end if
                end do 
            case default
				open(20, file="error.txt", access="append")
					write(20,*) "Warning: unknown class ", clas, " for amino acid ", restrict_aa
                close(20)
        end select
        endif
    enddo
    
	if(ip.eq.10) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_h_avail_num - 1.0e-3)+1   ! hydrationprop%hnum cannot be zero,不能排除一个类型所有的氨基酸
        aminoacid_name=res_h_avail(ip1) 
        
    elseif(ip.eq.20) then
        call ran_gen(ran2,0)
        ip1=int(ran2*res_p_avail_num - 1.0e-3)+1
        aminoacid_name=res_p_avail(ip1)
      
	elseif(ip.eq.30) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_c_avail_num - 1.0e-3)+1
        aminoacid_name=res_c_avail(ip1)
	elseif(ip.eq.40) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_o_avail_num - 1.0e-3)+1
        aminoacid_name=res_o_avail(ip1)

	elseif(ip.eq.11) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_h_avail_num - 1.0e-3)+1
        aminoacid_name="N"//hydrationprop%hgtype(ip1)

	elseif(ip.eq.21) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_p_avail_num - 1.0e-3)+1
        aminoacid_name="N"//res_p_avail(ip1)
        
	elseif(ip.eq.31) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_c_avail_num - 1.0e-3)+1
        aminoacid_name="N"//res_c_avail(ip1)

	elseif(ip.eq.41) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_o_avail_num - 1.0e-3)+1
        aminoacid_name="N"//res_o_avail(ip1)

	elseif(ip.eq.12) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_h_avail_num - 1.0e-3)+1
        aminoacid_name="C"//res_h_avail(ip1)        

	elseif(ip.eq.22) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_p_avail_num - 1.0e-3)+1
        aminoacid_name="C"//res_p_avail(ip1)			

	elseif(ip.eq.32) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_c_avail_num - 1.0e-3)+1
        aminoacid_name="C"//res_c_avail(ip1)    

	elseif(ip.eq.42) then
		call ran_gen(ran2,0)
        ip1=int(ran2*res_o_avail_num - 1.0e-3)+1
        aminoacid_name="C"//res_o_avail(ip1)            
    endif
    
    if (current_AA_name == aminoacid_name) goto 40
    
99  continue
    
	return
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

	
	subroutine scmf_choose_aminoacid(ip, aminoacid_number, aminoacid_name)
	implicit none
	integer					:: aminoacid_number
	integer					:: ip, ip1, ip2, i
	real					:: ran2
	character*4				:: char, aminoacid_name(10)

	if(ip==1) then
		aminoacid_number=hydrationprop%hnum
        aminoacid_name=hydrationprop%hgtype
		!aminoacid_name(1)="ALA"
		!aminoacid_name(2)="LEU"
		!aminoacid_name(3)="VAL"
		!aminoacid_name(4)="ILE"
		!aminoacid_name(5)="MET"
		!aminoacid_name(6)="PHE"
		!aminoacid_name(7)="TYR"
		!aminoacid_name(8)="TRP"
		!aminoacid_name(9)="GLY"	
	elseif(ip==2) then
		aminoacid_number=hydrationprop%pnum
        aminoacid_name=hydrationprop%pgtype
        
		!aminoacid_name(1)="ASN"
		!aminoacid_name(2)="GLN"
		!aminoacid_name(3)="SER"
		!aminoacid_name(4)="THR"
		!aminoacid_name(5)="HIE"
	elseif(ip==3) then
		aminoacid_number=hydrationprop%cnum
        aminoacid_name=hydrationprop%cgtype
		!aminoacid_name(1)="ARG"
		!aminoacid_name(2)="LYS"
		!aminoacid_name(3)="GLU"
		!aminoacid_name(4)="ASP"
	elseif(ip==4) then
		aminoacid_number=hydrationprop%onum
        aminoacid_name=hydrationprop%ogtype
		!aminoacid_name(1)="PRO"
		!aminoacid_name(2)="CYS"
	endif

	do i=1, (aminoacid_number-1)
		call ran_gen(ran2,0)
		ip1=int(ran2*aminoacid_number-1.0e-3)+1
		if(ip1.gt.aminoacid_number) ip1=aminoacid_number
		call ran_gen(ran2,0)
		ip2=int(ran2*aminoacid_number-1.0e-3)+1
		if(ip2.gt.aminoacid_number) ip2=aminoacid_number

		char=aminoacid_name(ip1)
		aminoacid_name(ip2)=aminoacid_name(ip1)
		aminoacid_name(ip1)=char
	enddo

	return
    end subroutine scmf_choose_aminoacid
	
    
    subroutine scmf_nmr_choose_aminoacid(aminoacid_name)
	implicit none
	integer					:: aminoacid_number
	integer					:: ip, i
	real					:: ran2
	character*4				:: aminoacid_name(10), k
    
    aminoacid_number=NMR_pool_size
    
    do i=1, NMR_pool_size
         aminoacid_name(i)=NMR_AA_pool(i)
    enddo
   
    do i = aminoacid_number, (aminoacid_number - nmr_site_num), -1
        call ran_gen(ran2,0)
		ip = int(ran2*aminoacid_number-1.0e-3)+1
        if(ip.gt.aminoacid_number) ip=aminoacid_number
        
        k = aminoacid_name(i)
        aminoacid_name(i) = aminoacid_name(ip)
        aminoacid_name(ip) = k
    enddo

	return
    end subroutine scmf_nmr_choose_aminoacid

    
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
	integer								:: dihedral_num, i, j
	character*4							:: gtype
	type(dihedralparameters)			:: dihedral	

	open(10, file=trim(mydir)//'lib/DihedralAngle/'//trim(gtype), status="old")
		read(10, "(i8)") dihedral_num
		do i=1, dihedral_num
			read(10,"(5i8)") dihedral%iph(i), dihedral%jph(i), dihedral%kph(i), dihedral%lph(i), dihedral%multiply(i)
			do j=1, dihedral%multiply(i)
				read(10,"(3e16.8)") dihedral%pk(i,j), dihedral%pn(i,j), dihedral%phase(i,j)
			enddo
		enddo
	close(10)
	
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

!	subroutine groupRES_1(ic, group, g_residue)
!	implicit none
!	integer							:: ic, chainID, i

!	type(groupdetails)				:: group(repeated_unit,gnum)
!	type(RES4chain)					:: g_residue	

!	chainID=1
!	g_residue%num=0
!	if(group(chainID,ic)%gtype=="ALA".or.group(chainID,ic)%gtype=="LEU".or.group(chainID,ic)%gtype=="VAL".or.group(chainID,ic)%gtype=="ILE".or.group(chainID,ic)%gtype=="MET".or. &
!	   group(chainID,ic)%gtype=="PHE".or.group(chainID,ic)%gtype=="TYR".or.group(chainID,ic)%gtype=="TRP".or.group(chainID,ic)%gtype=="GLY".or.group(chainID,ic)%gtype=="NALA".or.group(chainID,ic)%gtype=="NLEU".or. &
!	   group(chainID,ic)%gtype=="NVAL".or.group(chainID,ic)%gtype=="NILE".or.group(chainID,ic)%gtype=="NMET".or.group(chainID,ic)%gtype=="NPHE".or.group(chainID,ic)%gtype=="NTYR".or. &
!	   group(chainID,ic)%gtype=="NTRP".or.group(chainID,ic)%gtype=="NGLY".or.group(chainID,ic)%gtype=="CALA".or.group(chainID,ic)%gtype=="CLEU".or.group(chainID,ic)%gtype=="CVAL".or.group(chainID,ic)%gtype=="CILE".or. &
!	   group(chainID,ic)%gtype=="CMET".or.group(chainID,ic)%gtype=="CPHE".or.group(chainID,ic)%gtype=="CTYR".or.group(chainID,ic)%gtype=="CTRP".or.group(chainID,ic)%gtype=="CGLY") then
!	   	do i=1,gnum
!			if(group(chainID,i)%gtype=="ALA".or.group(chainID,i)%gtype=="LEU".or.group(chainID,i)%gtype=="VAL".or.group(chainID,i)%gtype=="ILE".or.group(chainID,i)%gtype=="MET".or. &
!			   group(chainID,i)%gtype=="PHE".or.group(chainID,i)%gtype=="TYR".or.group(chainID,i)%gtype=="TRP".or.group(chainID,i)%gtype=="GLY".or.group(chainID,i)%gtype=="NALA".or.group(chainID,i)%gtype=="NLEU".or. &
!			   group(chainID,i)%gtype=="NVAL".or.group(chainID,i)%gtype=="NILE".or.group(chainID,i)%gtype=="NMET".or.group(chainID,i)%gtype=="NPHE".or.group(chainID,i)%gtype=="NTYR".or. &
!			   group(chainID,i)%gtype=="NTRP".or.group(chainID,i)%gtype=="NGLY".or.group(chainID,i)%gtype=="CALA".or.group(chainID,i)%gtype=="CLEU".or.group(chainID,i)%gtype=="CVAL".or.group(chainID,i)%gtype=="CILE".or. &
!			   group(chainID,i)%gtype=="CMET".or.group(chainID,i)%gtype=="CPHE".or.group(chainID,i)%gtype=="CTYR".or.group(chainID,i)%gtype=="CTRP".or.group(chainID,i)%gtype=="CGLY") then
!				g_residue%num=g_residue%num+1
!				g_residue%IDs(g_residue%num)=i
!			endif
!		enddo
!	elseif(group(chainID,ic)%gtype=="GLU".or.group(chainID,ic)%gtype=="ASP".or.group(chainID,ic)%gtype=="ARG".or.group(chainID,ic)%gtype=="LYS".or.group(chainID,ic)%gtype=="ASN".or. &
!		   group(chainID,ic)%gtype=="GLN".or.group(chainID,ic)%gtype=="SER".or.group(chainID,ic)%gtype=="THR".or.group(chainID,ic)%gtype=="HIE".or.group(chainID,ic)%gtype=="NGLU".or. &
!		   group(chainID,ic)%gtype=="NASP".or.group(chainID,ic)%gtype=="NARG".or.group(chainID,ic)%gtype=="NLYS".or.group(chainID,ic)%gtype=="NASN".or.group(chainID,ic)%gtype=="NGLN".or. &
!		   group(chainID,ic)%gtype=="NSER".or.group(chainID,ic)%gtype=="NTHR".or.group(chainID,ic)%gtype=="NHIE".or.group(chainID,ic)%gtype=="CGLU".or.group(chainID,ic)%gtype=="CASP".or. &
!		   group(chainID,ic)%gtype=="CARG".or.group(chainID,ic)%gtype=="CLYS".or.group(chainID,ic)%gtype=="CASN".or.group(chainID,ic)%gtype=="CGLN".or.group(chainID,ic)%gtype=="CSER".or. &
!		   group(chainID,ic)%gtype=="CTHR".or.group(chainID,ic)%gtype=="CHIE") then
!		do i=1,gnum
!			if(group(chainID,i)%gtype=="GLU".or.group(chainID,i)%gtype=="ASP".or.group(chainID,i)%gtype=="ARG".or.group(chainID,i)%gtype=="LYS".or.group(chainID,i)%gtype=="ASN".or. &
!			   group(chainID,i)%gtype=="GLN".or.group(chainID,i)%gtype=="SER".or.group(chainID,i)%gtype=="THR".or.group(chainID,i)%gtype=="HIE".or.group(chainID,i)%gtype=="NGLU".or. &
!			   group(chainID,i)%gtype=="NASP".or.group(chainID,i)%gtype=="NARG".or.group(chainID,i)%gtype=="NLYS".or.group(chainID,i)%gtype=="NASN".or.group(chainID,i)%gtype=="NGLN".or. &
!			   group(chainID,i)%gtype=="NSER".or.group(chainID,i)%gtype=="NTHR".or.group(chainID,i)%gtype=="NHIE".or.group(chainID,i)%gtype=="CGLU".or.group(chainID,i)%gtype=="CASP".or. &
!			   group(chainID,i)%gtype=="CARG".or.group(chainID,i)%gtype=="CLYS".or.group(chainID,i)%gtype=="CASN".or.group(chainID,i)%gtype=="CGLN".or.group(chainID,i)%gtype=="CSER".or. &
!			   group(chainID,i)%gtype=="CTHR".or.group(chainID,i)%gtype=="CHIE") then
!				g_residue%num=g_residue%num+1
!				g_residue%IDs(g_residue%num)=i
!			endif
!		enddo
!	endif
	
!	return
!	end subroutine groupRES_1

	
end module database

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module transplant
	
	use constant
	use mathfunction
	use randomgenerator

	contains
	subroutine residue_replace(chainID, ic, group, ip, aa_group, temp_group)
	implicit none
	integer								:: chainID, ic, ip, ii, j, k, flag
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
							if(flag==1) then !! 如果是N-terminus AA, 将N的H替换为H1, 不是的话，则为H
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
			do j=1, aa_group(ii,ip)%cnum1
				if(aa_group(ii,ip)%atype1(j)=="HA2") then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="CA") then
							temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							temp_group(ii,ic)%coo1((k+1),1)=aa_group(ii,ip)%coo1(j,1)
							temp_group(ii,ic)%coo1((k+1),2)=aa_group(ii,ip)%coo1(j,2)
							temp_group(ii,ic)%coo1((k+1),3)=aa_group(ii,ip)%coo1(j,3)
							goto 30
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo										
				elseif(aa_group(ii,ip)%atype1(j)=="HA3") then
					temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1+1
					do k=(temp_group(ii,ic)%cnum1-1), 1, -1
						if(temp_group(ii,ic)%atype1(k)=="HA2") then
							temp_group(ii,ic)%atype1(k+1)=aa_group(ii,ip)%atype1(j)
							temp_group(ii,ic)%coo1((k+1),1)=aa_group(ii,ip)%coo1(j,1)
							temp_group(ii,ic)%coo1((k+1),2)=aa_group(ii,ip)%coo1(j,2)
							temp_group(ii,ic)%coo1((k+1),3)=aa_group(ii,ip)%coo1(j,3)
							goto 30
						else
							temp_group(ii,ic)%atype1(k+1)=temp_group(ii,ic)%atype1(k)
							temp_group(ii,ic)%coo1((k+1),1)=temp_group(ii,ic)%coo1(k,1)
							temp_group(ii,ic)%coo1((k+1),2)=temp_group(ii,ic)%coo1(k,2)
							temp_group(ii,ic)%coo1((k+1),3)=temp_group(ii,ic)%coo1(k,3)	
						endif
					enddo
				endif					
30				continue
			enddo
			temp_group(ii,ic)%cnum1=temp_group(ii,ic)%cnum1-1
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
	integer					:: ip1, ip2, i, j, k, chain
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
	
	chain=repeated_unit/2
	Tgroup=group
    
	do i=1, chain
        do j=1, gnum
            do k=1, Tgroup(i,j)%cnum1
                Tgroup(i,j)%coo1(k,ip1)=Tgroup(i,j)%coo1(k,ip1) + displacement
            enddo
			
            do k=1, Tgroup(i,j)%cnum2
                Tgroup(i,j)%coo2(k,ip1)=Tgroup(i,j)%coo2(k,ip1) + displacement
            enddo
			
            do k=1, Tgroup(i,j)%cnum3
                Tgroup(i,j)%coo3(k,ip1)=Tgroup(i,j)%coo3(k,ip1) + displacement
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


	subroutine bindingenergy_noentropy(group, group_para, numex, inb, numex4, inb4, score, binding_energy, comp_vdw, score4hydration, Pagg)
	implicit none
	integer							:: ligan_recep_sta(num4category), ligan_recep_end(num4category)	
	integer							:: natom, ipres, atomid(repeated_unit*atom_num)
	integer    						:: numex(repeated_unit*atom_num), inb(repeated_unit*atom_num,20), numex4(repeated_unit*atom_num), inb4(repeated_unit*atom_num,60)
	real							:: mdcrd(repeated_unit*atom_num,3), charge(repeated_unit*atom_num), epsion(repeated_unit*atom_num)
	real							:: r(repeated_unit*atom_num), rborn(repeated_unit*atom_num), fs(repeated_unit*atom_num), dielecons(repeated_unit*atom_num)
	character*4						:: lbres(repeated_unit*atom_num)

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
		deallocate (alpha)
		
	enddo
	
	binding_vdw=0; binding_ele=0; binding_sgb=0; binding_snp=0
	do categoryID=1, num4category
		binding_vdw=binding_vdw+ligan_recep_vdw(categoryID)
		binding_ele=binding_ele+ligan_recep_ele(categoryID)
		binding_sgb=binding_sgb+ligan_recep_sgb(categoryID)
	enddo

    binding_energy=((comp_vdw+comp_ele+comp_sgb+comp_snp)-(binding_vdw+binding_ele+binding_sgb+binding_snp))/real(num4peptides)

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

	
	subroutine bindingenergy(group, group_para, entropy4individual, numex, inb, numex4, inb4, score, binding_energy, entropy, score4hydration, Pagg)
	implicit none
	integer							:: ligan_recep_sta(num4category), ligan_recep_end(num4category)	
	integer							:: natom, ipres, atomid(repeated_unit*atom_num)
	integer    						:: numex(repeated_unit*atom_num), inb(repeated_unit*atom_num,20), numex4(repeated_unit*atom_num), inb4(repeated_unit*atom_num,60)
	real							:: mdcrd(repeated_unit*atom_num,3), charge(repeated_unit*atom_num), epsion(repeated_unit*atom_num)
	real							:: r(repeated_unit*atom_num), rborn(repeated_unit*atom_num), fs(repeated_unit*atom_num), dielecons(repeated_unit*atom_num)
	real							:: entropy4individual(repeated_unit,gnum)
	character*4						:: lbres(repeated_unit*atom_num)

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
		deallocate (alpha)
		
	enddo
	
	binding_vdw=0; binding_ele=0; binding_sgb=0; binding_snp=0
	do categoryID=1, num4category
		binding_vdw=binding_vdw+ligan_recep_vdw(categoryID)
		binding_ele=binding_ele+ligan_recep_ele(categoryID)
		binding_sgb=binding_sgb+ligan_recep_sgb(categoryID)
!		binding_snp=binding_snp+ligan_recep_snp(categoryID)
	enddo

    binding_energy=((comp_vdw+comp_ele+comp_sgb+comp_snp)-(binding_vdw+binding_ele+binding_sgb+binding_snp))/real(num4peptides)  
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
    
    
    subroutine replace_extra_restricted_aa(group, temp_group)
	implicit none
	integer 			:: restrict_count, non_restrict_count, restrict_locs(40), non_restrict_locs(40), chainID, index
	integer 			:: iii, i, j, k, num_over, index_non_asn, flag, flag1, flag2, rotanum, feedback, ip, ic
    type(groupdetails)  :: group(repeated_unit,gnum), temp_group(repeated_unit,gnum),  aa_group(repeated_unit,40)
	real				:: ran2
	character*4			:: group_name_1(3), group_name_2(3), aminoacid_name_1, aminoacid_name_2, restrict_aa
	
    chainID = 1
    
    do iii=1, n_restrictions
        restrict_locs = 0
        non_restrict_locs = 0						! initial position of non-restricted AA
        restrict_count = 0							! initialize how many this restricted AA
        non_restrict_count = 0						! initialize AA other than restricted AA
        restrict_aa = aa_restrictions(iii)%gtype	! get AA name of restricted AA
		index=0
		do i=1, gnum
			if ((group(chainID,i)%gtype=="N"//restrict_aa).or.(group(chainID,i)%gtype==restrict_aa).or.(group(chainID,i)%gtype=="C"//restrict_aa)) then
				restrict_locs(restrict_count+1) = i
				restrict_count = restrict_count + 1
			else
				non_restrict_locs(non_restrict_count+1) = i
				non_restrict_count = non_restrict_count + 1
			endif
        enddo
        
        num_over = restrict_count -  aa_restrictions(iii)%max
        
		temp_group = group
		do i=1,num_over
			!pick a randomly restricted site, and do MC replacement
			call ran_gen(ran2,0)
			index=int(ran2*restrict_count-1.0e-3)+1
			if (index.gt.restrict_count) index = restrict_count ! index is the index of restricted AA in its array
		
			ic = restrict_locs(index)							! actual position in peptide sequence
            
			call mc_choose_aminoacid(ic, group, aminoacid_name_1)
			call findrotamer(ic, temp_group, aminoacid_name_1, rotanum, aa_group, ip)

			do j=1, rotanum
            
				call residue_replace(chainID, ic, group, j, aa_group, temp_group)

				call check_transplant(chainID, ic, temp_group, feedback)
			
				if(feedback==1) then
					group=temp_group
					non_restrict_locs(non_restrict_count+1) = restrict_locs(index)
					do k=index, restrict_count-1
						restrict_locs(k) = restrict_locs(k+1)
					enddo
					restrict_count = restrict_count - 1
					non_restrict_count = non_restrict_count + 1
					goto 10
				endif
			enddo
			open(20, file="error.txt", access="append")
				write(20,*) "Unable to replace ", restrict_aa, " at site ", ic, ". Terminating program."
				close(20)
			stop
			!if we reach this part, then we couldn't replace the ASN with SER
			!for now, will have program exit since I don't expect this will happen often.
			!if this is a problem, then will add a solution         
10			continue
		enddo
	enddo
	return
	end subroutine replace_extra_restricted_aa
    
    
	subroutine scmf_substitution(group, sub_cycle, tgroup)
	implicit none
	integer							:: sub_cycle, i, j, k, l, m, n, ic, ic1, ic2, flag, flag1, flag2
	integer							:: ii, chainID, ip
	integer							:: Npho, Npol, Nchg, Noth
	integer							:: positive_number, negative_number, aminoacid_number, rotanum, feedback
	integer							:: icpointer(4,gnum), pos_type(gnum), neg_type(gnum)
	integer							:: group_num(4), num(4)
	integer							:: num4peptides, peptideID(8)
	real							:: ran2
	character*4						:: aminoacid_name(10), aminoacid_name_1, group_name_1(3), group_name_2(3)
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum), temp_group_1(repeated_unit,gnum), temp_group_2(repeated_unit,gnum), aa_group(repeated_unit,40)

	tgroup=group
	group_num=0
	ii=1
	do ic=1, gnum
		do j=1,void_site_num
			if(ic.ge.void_site_start(j).and.ic.le.void_site_end(j)) then 
				goto 15
			endif
        enddo
        
        do j=1, nmr_site_num
            if(ic == nmr_site_ID(j)) then 
				goto 15
            endif
        enddo
        
		if(tgroup(ii,ic)%gtype=="ALA".or.tgroup(ii,ic)%gtype=="LEU".or.tgroup(ii,ic)%gtype=="VAL".or.tgroup(ii,ic)%gtype=="ILE".or.tgroup(ii,ic)%gtype=="MET".or.tgroup(ii,ic)%gtype=="PHE".or.tgroup(ii,ic)%gtype=="TYR".or.tgroup(ii,ic)%gtype=="TRP".or.tgroup(ii,ic)%gtype=="GLY"  &
			.or.tgroup(ii,ic)%gtype=="NALA".or.tgroup(ii,ic)%gtype=="NLEU".or.tgroup(ii,ic)%gtype=="NVAL".or.tgroup(ii,ic)%gtype=="NILE".or.tgroup(ii,ic)%gtype=="NMET".or.tgroup(ii,ic)%gtype=="NPHE".or.tgroup(ii,ic)%gtype=="NTYR".or.tgroup(ii,ic)%gtype=="NTRP".or.tgroup(ii,ic)%gtype=="NGLY"  &
			.or.tgroup(ii,ic)%gtype=="CALA".or.tgroup(ii,ic)%gtype=="CLEU".or.tgroup(ii,ic)%gtype=="CVAL".or.tgroup(ii,ic)%gtype=="CILE".or.tgroup(ii,ic)%gtype=="CMET".or.tgroup(ii,ic)%gtype=="CPHE".or.tgroup(ii,ic)%gtype=="CTYR".or.tgroup(ii,ic)%gtype=="CTRP".or.tgroup(ii,ic)%gtype=="CGLY") then			   			   			   
			group_num(1)=group_num(1)+1
			icpointer(1, group_num(1))=ic
		elseif(tgroup(ii,ic)%gtype=="ASN".or.tgroup(ii,ic)%gtype=="GLN".or.tgroup(ii,ic)%gtype=="SER".or.tgroup(ii,ic)%gtype=="THR".or.tgroup(ii,ic)%gtype=="HIE" &
			.or.tgroup(ii,ic)%gtype=="NASN".or.tgroup(ii,ic)%gtype=="NGLN".or.tgroup(ii,ic)%gtype=="NSER".or.tgroup(ii,ic)%gtype=="NTHR".or.tgroup(ii,ic)%gtype=="NHIE" &
			.or.tgroup(ii,ic)%gtype=="CASN".or.tgroup(ii,ic)%gtype=="CGLN".or.tgroup(ii,ic)%gtype=="CSER".or.tgroup(ii,ic)%gtype=="CTHR".or.tgroup(ii,ic)%gtype=="CHIE") then			
			group_num(2)=group_num(2)+1
			icpointer(2, group_num(2))=ic		
		elseif(tgroup(ii,ic)%gtype=="ARG".or.tgroup(ii,ic)%gtype=="LYS".or.tgroup(ii,ic)%gtype=="GLU".or.tgroup(ii,ic)%gtype=="ASP" &
			.or.tgroup(ii,ic)%gtype=="NARG".or.tgroup(ii,ic)%gtype=="NLYS".or.tgroup(ii,ic)%gtype=="NGLU".or.tgroup(ii,ic)%gtype=="NASP" &
			.or.tgroup(ii,ic)%gtype=="CARG".or.tgroup(ii,ic)%gtype=="CLYS".or.tgroup(ii,ic)%gtype=="CGLU".or.tgroup(ii,ic)%gtype=="CASP") then			
			group_num(3)=group_num(3)+1
			icpointer(3, group_num(3))=ic
		elseif(tgroup(ii,ic)%gtype=="PRO".or.tgroup(ii,ic)%gtype=="CYS"  &
		    .or.tgroup(ii,ic)%gtype=="NPRO".or.tgroup(ii,ic)%gtype=="NCYS"  &
			.or.tgroup(ii,ic)%gtype=="CPRO".or.tgroup(ii,ic)%gtype=="CCYS") then
			group_num(4)=group_num(4)+1
			icpointer(4, group_num(4))=ic
		endif
15		continue
	enddo
	
	Npho=fpho; Npol=fpol; Nchg=fchg; Noth=foth
	
	if((Npho+Npol+Nchg+Noth).ne.(group_num(1)+group_num(2)+group_num(3)+group_num(4))) then
		open(20, file="error.txt", access="append")
			write(20,*) "Input setting:       ", "Npho=", Npho, "Npol=", Npol, "Nchg=", Nchg, "Noth=", Noth
			write(20,*) "Program calculation: ", "Npho=", group_num(1), "Npol=", group_num(2), "Nchg=", group_num(3), "Noth=", group_num(4)
			write(20,*) "Please adjust the number of residue in each category of residue types in the Input file!"
		close(20)
		stop	
	endif

	num(1)=group_num(1)-Npho
	num(2)=group_num(2)-Npol
	num(3)=group_num(3)-Nchg
	num(4)=group_num(4)-Noth

	positive_number=0
	negative_number=0
	do i=1, 4
		if(num(i)>0) then
			do j=1, (group_num(i)-1)
				call ran_gen(ran2,0)
				ic1=int(ran2*group_num(i)-1.0e-3)+1
				if(ic1.gt.group_num(i)) ic1=group_num(i)
				call ran_gen(ran2,0)
				ic2=int(ran2*group_num(i)-1.0e-3)+1
				if(ic2.gt.group_num(i)) ic2=group_num(i)

				k=icpointer(i,ic1)
				icpointer(i,ic1)=icpointer(i,ic2)
				icpointer(i,ic2)=k
			enddo
			do j=1, num(i)
				pos_type(positive_number+j)=i
			enddo
			positive_number=positive_number+num(i)		
		elseif(num(i)<0) then
			do j=1, abs(num(i))
				neg_type(negative_number+j)=i
			enddo
			negative_number=negative_number+abs(num(i))
		endif
	enddo

	do i=1, (positive_number-1)
		call ran_gen(ran2,0)
		ic1=int(ran2*positive_number-1.0e-3)+1
		if(ic1.gt.positive_number) ic1=positive_number
		call ran_gen(ran2,0)
		ic2=int(ran2*positive_number-1.0e-3)+1
		if(ic2.gt.positive_number) ic2=positive_number

		k=pos_type(ic1)
		pos_type(ic1)=pos_type(ic2)
		pos_type(ic2)=k
	enddo
	do i=1, (negative_number-1)
		call ran_gen(ran2,0)
		ic1=int(ran2*negative_number-1.0e-3)+1
		if(ic1.gt.negative_number) ic1=negative_number
		call ran_gen(ran2,0)
		ic2=int(ran2*negative_number-1.0e-3)+1
		if(ic2.gt.negative_number) ic2=negative_number

		k=neg_type(ic1)
		neg_type(ic1)=neg_type(ic2)
		neg_type(ic2)=k
    enddo

	temp_group_1=tgroup
	sub_cycle=min(positive_number, negative_number)
	do i=1, sub_cycle
		call scmf_choose_aminoacid(neg_type(i), aminoacid_number, aminoacid_name)
		l=1
		do while(l<=group_num(pos_type(i)))
			do j=1, aminoacid_number
				call groupinfo(temp_group_1(ii,icpointer(pos_type(i),l))%gtype, group_name_1, flag1)
				call groupinfo(aminoacid_name(j), group_name_2, flag2)

				aminoacid_name_1=group_name_2(flag1)

				call findrotamer(icpointer(pos_type(i),l), temp_group_1, aminoacid_name_1, rotanum, aa_group, ip)				

				flag=0
				chainID=1
				do m=1, rotanum
					call residue_replace(chainID, icpointer(pos_type(i),l), temp_group_1, m, aa_group, temp_group_2)

					call check_transplant(chainID, icpointer(pos_type(i),l), temp_group_2, feedback)
			
					if(feedback==1) then
						temp_group_1=temp_group_2
						flag=flag+1
						goto 10
					endif
				enddo
10				continue
				!enddo
				if(flag==1) goto 20				
			enddo
			l=l+1		
		enddo
20		continue

		if(flag==1) then
			do k=1, (group_num(pos_type(i))+l)
				if(k<=l) then
					icpointer(pos_type(i),(group_num(pos_type(i))+k))=icpointer(pos_type(i),k)
				else
					icpointer(pos_type(i),(k-l))=icpointer(pos_type(i),k)
				endif
			enddo
			group_num(pos_type(i))=group_num(pos_type(i))-1
		else
			open(20, file="error.txt", access="append")
				write(20,*) "scmf_substitution is wrong!"
				write(20,*) "Please check whether the code in the module MAINTENACE is right or not!"
			close(20)
			stop
		endif
    enddo

    call scmf_nmr_choose_aminoacid(aminoacid_name)
    do i=1, nmr_site_num
        
        call groupinfo(temp_group_1(ii,nmr_site_ID(i))%gtype, group_name_1, flag1) ! 返回nmr site 的AA类型(flag1)
        call groupinfo(aminoacid_name(i), group_name_2, flag2)  ! aminoacid_name 包含要替换的4个 NMR AA， 返回group_name_2 (AA的3种类型)
        aminoacid_name_1=group_name_2(flag1)
        
        call findrotamer(nmr_site_ID(i), temp_group_1, aminoacid_name_1, rotanum, aa_group, ip)
        flag=0
        chainID=1
        do m=1, rotanum
            call residue_replace(chainID, nmr_site_ID(i), temp_group_1, m, aa_group, temp_group_2)

            call check_transplant(chainID, nmr_site_ID(i), temp_group_2, feedback)
			
            if(feedback==1) then
                temp_group_1=temp_group_2
                flag=flag+1
                goto 97
            endif
        enddo
97      continue
    enddo
    
	tgroup=temp_group_1
	
	return 
	end subroutine scmf_substitution


	subroutine sidechain_optimization(stage, chainID, ic, group, group_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)
	implicit none
	integer								:: grade, grade_num(6), monitor(6)
	integer								:: chainID, i, j, k, ic, account_num, flag, stage, trial_count
	integer								:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60, 60)
	integer								:: dihedral_num	
	real								:: delta_chi, cos_angle, sin_angle, error, t
	real								:: CA(3), rotaxis_x, rotaxis_y, rotaxis_z, rotaxis(3), m(3,3), Tmember(15,3)
	real(kind=16)                       :: h2_denominator, h3_denominator, Tenergy_min, Tenergy
	real								:: Dihedral4entropy(4)
	type(groupdetails)					:: group(repeated_unit,gnum), Tgroup(repeated_unit,gnum)
	type(energyparameters)				:: group_para(repeated_unit,gnum)
	type(index4sidechain)				:: index(60)
	type(conformer4sidechain)			:: Iclass(6), Tclass(6), class_old(6), class_new(6), class_min(6)
	type(dihedralparameters)			:: dihedral	
	
	real(kind=16), dimension(:), allocatable 	:: energy_forward, energy_backward
	real(kind=16), dimension(:,:), allocatable 	:: gradient_old, gradient_new, Hessian_old, Hessian_new
	real(kind=16), dimension(:,:), allocatable   :: H2, H3, H31
	real(kind=16), dimension(:,:), allocatable	:: d, y, s, Tchi

	call sidechain_category(chainID, ic, group, Iclass, grade, grade_num, index, monitor) ! find AA's side-chain info
	call dihedralangle_reading(group(chainID,ic)%gtype, dihedral_num, dihedral)			  ! find AA's dihedral angle info

	allocate(energy_forward(grade)); allocate(energy_backward(grade))
	allocate(gradient_old(grade,1)); allocate(gradient_new(grade,1))
	allocate(Hessian_old(grade,grade)); allocate(Hessian_new(grade,grade))
	allocate(H2(grade,grade)); allocate(H3(grade,grade)); allocate(H31(grade,1))
	allocate(d(grade,1)); allocate(y(grade,1)); allocate(s(grade,1))

	Tgroup=group
	do i=1, Tgroup(chainID,ic)%cnum1
		if(Tgroup(chainID,ic)%atype1(i)=="CA") then										! find coord for CA
			CA(1)=Tgroup(chainID,ic)%coo1(i,1); CA(2)=Tgroup(chainID,ic)%coo1(i,2); CA(3)=Tgroup(chainID,ic)%coo1(i,3)
		endif
	enddo
	s=0.0
	class_new=Iclass

30	continue
	if(stage==0) then               ! stage=0,dk=5° (微小变化角度), stage=1,dk=1°
		delta_chi=5
	elseif(stage==1) then
		delta_chi=1
	endif
	cos_angle=cosd(delta_chi); sin_angle=sind(delta_chi)	! dk 微小角度的cos值和sin值
	do i=1, grade											! loop all grade (number of heavy atom dihedral axis)
		Tclass=class_new
		if(i==1) then										! 第一个axis是 CA-CB 轴
			rotaxis_x=Tclass(i)%member(monitor(i),1)-CA(1)	! 计算 CA-CB 轴向量
			rotaxis_y=Tclass(i)%member(monitor(i),2)-CA(2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-CA(3)
			rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2) ! 计算 CA-CB 轴单位向量
			rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		else												! 计算 i>1 轴向量
			rotaxis_x=Tclass(i)%member(monitor(i),1)-Tclass(i-1)%member(monitor(i-1),1)
			rotaxis_y=Tclass(i)%member(monitor(i),2)-Tclass(i-1)%member(monitor(i-1),2)
			rotaxis_z=Tclass(i)%member(monitor(i),3)-Tclass(i-1)%member(monitor(i-1),3)
			rotaxis(1)=rotaxis_x/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(2)=rotaxis_y/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
			rotaxis(3)=rotaxis_z/sqrt(rotaxis_x**2+rotaxis_y**2+rotaxis_z**2)
		endif

		call axisrotation(rotaxis, cos_angle, sin_angle, m)
		
		do j=(i+1), (grade+1) ! 获得所有原子相对于旋转轴i的一端重原子的坐标
			do k=1, grade_num(j)
				Tclass(j)%member(k,1)=Tclass(j)%member(k,1)-Tclass(i)%member(monitor(i),1)
				Tclass(j)%member(k,2)=Tclass(j)%member(k,2)-Tclass(i)%member(monitor(i),2)
				Tclass(j)%member(k,3)=Tclass(j)%member(k,3)-Tclass(i)%member(monitor(i),3)
			enddo
			
			Tmember=matmul(Tclass(j)%member, m)
			Tclass(j)%member=Tmember
			
			do k=1, grade_num(j) !转换为原坐标，保留三位小数
				Tclass(j)%member(k,1)=anint((Tclass(j)%member(k,1)+Tclass(i)%member(monitor(i),1))*1000)/1000
				Tclass(j)%member(k,2)=anint((Tclass(j)%member(k,2)+Tclass(i)%member(monitor(i),2))*1000)/1000
				Tclass(j)%member(k,3)=anint((Tclass(j)%member(k,3)+Tclass(i)%member(monitor(i),3))*1000)/1000
			enddo
		enddo
		
		do j=1, Tgroup(chainID,ic)%cnum2 ! 替换Tgroup中侧链上的坐标为旋转后的坐标
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
		
		allocate(Tchi(grade,1))
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
		deallocate(Tchi)
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

	deallocate(energy_forward); deallocate(energy_backward)
	deallocate(gradient_old); deallocate(gradient_new)
	deallocate(Hessian_old); deallocate(Hessian_new)
	deallocate(H2); deallocate(H3); deallocate(H31)	
	deallocate(d); deallocate(y); deallocate(s)
	
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

	contains
	subroutine MC_technique_sequence(score_old, energy_old, entropy_old, score4hydration_old, Pagg_old, group, entropy4individual, score_new, energy_new, entropy_new, score4hydration_new, Pagg_new, tgroup, Tentropy4individual, feedback)
	implicit none
	integer							:: feedback
	real							:: ran2
	real(kind=8)                    :: energy_old, energy_new, score_old, score_new, score_change
	real							:: entropy_old, entropy_new
	real							:: score4hydration_old, Pagg_old
	real							:: score4hydration_new, Pagg_new
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum)
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)

	feedback=0
	score_change=score_new-score_old
	if(score_change.le.0) then
		group=tgroup
		score_old=score_new
		energy_old=energy_new
		entropy_old=entropy_new
		score4hydration_old=score4hydration_new
		Pagg_old=Pagg_new
		entropy4individual=Tentropy4individual
		feedback=1
	else
		call ran_gen(ran2,0)
		if(ran2.le.exp(-score_change/ekt)) then
			group=tgroup
			score_old=score_new
			energy_old=energy_new
			entropy_old=entropy_new
			score4hydration_old=score4hydration_new
			Pagg_old=Pagg_new
			entropy4individual=Tentropy4individual
			feedback=1
		endif
	endif
	
	return
	end subroutine MC_technique_sequence
	
	subroutine MC_technique_sheet(score_old, energy_old, entropy_old, score4hydration_old, Pagg_old, group, entropy4individual, score_new, energy_new, entropy_new, score4hydration_new, Pagg_new, tgroup, Tentropy4individual, feedback)
	implicit none
	integer							:: feedback
    real                            :: ran2
	real(kind=8)                    :: score_old, score_new, score_change
	real(kind=8)                    :: energy_old, energy_new
	real							:: entropy_old, entropy_new
	real							:: score4hydration_old, Pagg_old
	real							:: score4hydration_new, Pagg_new
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum)
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)
	
	feedback=0
	score_change=score_new-score_old
	if(score_change.le.0) then
		group=tgroup
		score_old=score_new
		energy_old=energy_new
		entropy_old=entropy_new
		score4hydration_old=score4hydration_new
		Pagg_old=Pagg_new
		entropy4individual=Tentropy4individual
		feedback=1
	else
		call ran_gen(ran2,0)
		if(ran2.le.exp(-score_change/ekt_sheet)) then
			group=tgroup
			score_old=score_new
			energy_old=energy_new
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
	real							:: vdw_energy_min, score_min
	real							:: entropy, score4hydration, Pagg
	real							:: Dihedral4entropy(4)
	character*4						:: aminoacid_name

	type(groupdetails)				:: aa_backup
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)
	type(groupdetails), dimension(:,:), allocatable &
									:: temp_group, aa_group
	type(energyparameters), dimension(:,:), allocatable &
									:: tgroup_para


	allocate(aa_group(repeated_unit,40))
	call findrotamer(ic, group, aminoacid_name, rotanum, aa_group, ip)
	
	flag=0
	if(aminoacid_name=="GLY".or.aminoacid_name=="ALA".or.aminoacid_name=="PRO".or. &
	   aminoacid_name=="NGLY".or.aminoacid_name=="NALA".or.aminoacid_name=="NPRO".or. &
	   aminoacid_name=="CGLY".or.aminoacid_name=="CALA".or.aminoacid_name=="CPRO".or. &
	   aminoacid_name=="NME".or.aminoacid_name=="NHE".or.aminoacid_name=="ACE") then	

		allocate(temp_group(repeated_unit,gnum))
		allocate(tgroup_para(repeated_unit,gnum))

		tgroup=group
		do chainID=1, repeated_unit
				
			m_best=0
			vdw_energy_min=500.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)

				if(m==1) then
					call energy_parameter(temp_group, tgroup_para)
					call atom_links4sidechain(chainID, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
				endif

				call check_transplant(chainID, ic, temp_group, feedback)

				if(feedback==1) then
					call vdwenergy(chainID, ic, temp_group, tgroup_para, vdw_energy)
					if(vdw_energy.lt.vdw_energy_min) then
						vdw_energy_min=vdw_energy
						call backup4sidechain(0, chainID, ic, temp_group, aa_backup)
						m_best=m
					endif
				endif
			enddo

			if(m_best==0) goto 10
			call backup4sidechain(1, chainID, ic, tgroup, aa_backup)
		enddo

		flag=1
10		continue
		deallocate(tgroup_para)
		deallocate(temp_group)

	else
		allocate(temp_group(repeated_unit,gnum))
		allocate(tgroup_para(repeated_unit,gnum))

		tgroup=group
        vdw_energy_min=1000000000.0*gnum
        
		do chainID=1, repeated_unit
            
			m_best=0            
			score_min=500.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)
				if(m==1) then
					call energy_parameter(temp_group, tgroup_para)
					call atom_links4sidechain(chainID, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
					call atom_links(temp_group, W_numex, W_inb, W_numex4, W_inb4)
				endif
					
				stage=0
				call sidechain_optimization(stage, chainID, ic, temp_group, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)

				if(stage==1) then
					call bindingenergy_noentropy(temp_group, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, binding_energy, vdw_energy, score4hydration, Pagg)
					if((score.lt.score_min) .and. (vdw_energy.lt.vdw_energy_min)) then
						score_min=score
						call backup4sidechain(0, chainID, ic, temp_group, aa_backup)
						m_best=m
					endif
				endif
			enddo
				
			if(m_best==0) goto 20
			call backup4sidechain(1, chainID, ic, tgroup, aa_backup)
		enddo

		flag=1
20		continue
		deallocate(tgroup_para)
		deallocate(temp_group)
	endif
	deallocate(aa_group)
	
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
	real							:: entropy, score4hydration, Pagg
	real							:: Dihedral4entropy(4), matrix(34,4)
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum)
	character*4						:: aminoacid_name

	type(groupdetails)				:: aa_backup
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)
	type(groupdetails), dimension(:,:), allocatable &
									:: temp_group, aa_group
	type(energyparameters), dimension(:,:), allocatable &
									:: tgroup_para

	allocate(aa_group(repeated_unit,40))
	call findrotamer(ic, group, aminoacid_name, rotanum, aa_group, ip)
	grade=aa_lib(ip)%grade
	
	flag=0
	if(aminoacid_name=="GLY".or.aminoacid_name=="ALA".or.aminoacid_name=="PRO".or. &
	   aminoacid_name=="NGLY".or.aminoacid_name=="NALA".or.aminoacid_name=="NPRO".or. &
	   aminoacid_name=="CGLY".or.aminoacid_name=="CALA".or.aminoacid_name=="CPRO".or. &
	   aminoacid_name=="NME".or.aminoacid_name=="NHE".or.aminoacid_name=="ACE") then	

		allocate(temp_group(repeated_unit,gnum))
		allocate(tgroup_para(repeated_unit,gnum))

		tgroup=group
		Tentropy4individual=entropy4individual
		do chainID=1, repeated_unit
				
			m_best=0
			vdw_energy_min=500.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)

				if(m==1) then
					call energy_parameter(temp_group, tgroup_para)
					call atom_links4sidechain(chainID, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
				endif

				call check_transplant(chainID, ic, temp_group, feedback)

				if(feedback==1) then
					call vdwenergy(chainID, ic, temp_group, tgroup_para, vdw_energy)
					if(vdw_energy.lt.vdw_energy_min) then
						vdw_energy_min=vdw_energy
						call backup4sidechain(0, chainID, ic, temp_group, aa_backup)
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
		deallocate(tgroup_para)
		deallocate(temp_group)
											
	else
		allocate(temp_group(repeated_unit,gnum))
		allocate(tgroup_para(repeated_unit,gnum))

		tgroup=group
		Tentropy4individual=entropy4individual
        vdw_energy_min=1000000000.0*gnum
		do chainID=1, repeated_unit

			m_best=0
			score_min=500.0
			obs=0
			matrix=0.0
			do m=1, rotanum
				call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)
				if(m==1) then
					call energy_parameter(temp_group, tgroup_para)
					call atom_links4sidechain(chainID, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
					call atom_links(temp_group, W_numex, W_inb, W_numex4, W_inb4)
				endif
					
				stage=0
				call sidechain_optimization(stage, chainID, ic, temp_group, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)
						
				if(stage==1) then
					obs=obs+1
					do j=1, grade
						matrix(obs,j)=Dihedral4entropy(j) ! why not pick the matrix with low
					enddo
						
					call bindingenergy_noentropy(temp_group, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, binding_energy, vdw_energy, score4hydration, Pagg)
					if((score.lt.score_min) .and. (vdw_energy.lt.vdw_energy_min)) then !.and. vdw.lt.vdw_min
						score_min=score
						call backup4sidechain(0, chainID, ic, temp_group, aa_backup)
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
		deallocate(tgroup_para)
		deallocate(temp_group)
	endif
	deallocate(aa_group)

	return
	end subroutine sequence_mutation
	
	subroutine sequence_mutation_4scmfloop(chainID, ic, group, aminoacid_name, tgroup, flag)
	implicit none
	integer							:: flag, i, j, m, m_best, rotanum, feedback
	integer							:: chainID, ic, ic_1, ic_2, stage, ip
	integer							:: W_numex(repeated_unit*atom_num), W_inb(repeated_unit*atom_num,20), W_numex4(repeated_unit*atom_num), W_inb4(repeated_unit*atom_num,60)
	integer							:: S_numex(60), S_inb(60,20), S_numex4(60), S_inb4(60,60)
	real							:: ran2
    real(kind=8)					:: score, binding_energy, vdw_energy
	real							:: vdw_energy_min, score_min
	real							:: entropy, score4hydration, Pagg
	real							:: Dihedral4entropy(4)
	character*4						:: aminoacid_name

	type(groupdetails)				:: aa_backup
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum)
	type(groupdetails), dimension(:,:), allocatable &
									:: temp_group, aa_group
	type(energyparameters), dimension(:,:), allocatable &
									:: tgroup_para


	allocate(aa_group(repeated_unit,40))
	call findrotamer(ic, group, aminoacid_name, rotanum, aa_group, ip)
	
	flag=0
	if(aminoacid_name=="GLY".or.aminoacid_name=="ALA".or.aminoacid_name=="PRO".or. &
	   aminoacid_name=="NGLY".or.aminoacid_name=="NALA".or.aminoacid_name=="NPRO".or. &
	   aminoacid_name=="CGLY".or.aminoacid_name=="CALA".or.aminoacid_name=="CPRO".or. &
	   aminoacid_name=="NME".or.aminoacid_name=="NHE".or.aminoacid_name=="ACE") then	

		allocate(temp_group(repeated_unit,gnum))
		allocate(tgroup_para(repeated_unit,gnum))

		tgroup=group
		m_best=0
		vdw_energy_min=500.0
		do m=1, rotanum
			call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)

			if(m==1) then
				call energy_parameter(temp_group, tgroup_para)
				call atom_links4sidechain(chainID, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
			endif
			call check_transplant(chainID, ic, temp_group, feedback)

			if(feedback==1) then
				call vdwenergy(chainID, ic, temp_group, tgroup_para, vdw_energy)
				if(vdw_energy.lt.vdw_energy_min) then
					vdw_energy_min=vdw_energy
					call backup4sidechain(0, chainID, ic, temp_group, aa_backup)
					m_best=m
				endif
			endif
		enddo

		if(m_best==0) goto 10
		call backup4sidechain(1, chainID, ic, tgroup, aa_backup)

		flag=1
10		continue
		deallocate(tgroup_para)
		deallocate(temp_group)

	else
		allocate(temp_group(repeated_unit,gnum))
		allocate(tgroup_para(repeated_unit,gnum))

		tgroup=group
		m_best=0
        vdw_energy_min=1000000000.0*gnum
		score_min=500.0
		do m=1, rotanum
			call residue_replace(chainID, ic, tgroup, m, aa_group, temp_group)
			if(m==1) then
				call energy_parameter(temp_group, tgroup_para)
				call atom_links4sidechain(chainID, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
				call atom_links(temp_group, W_numex, W_inb, W_numex4, W_inb4)
			endif
					
			stage=0
			call sidechain_optimization(stage, chainID, ic, temp_group, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)

			if(stage==1) then
				call bindingenergy_noentropy(temp_group, tgroup_para, W_numex, W_inb, W_numex4, W_inb4, score, binding_energy, vdw_energy, score4hydration, Pagg)
                if(score.lt.score_min .and. vdw_energy.lt.vdw_energy_min) then
					score_min=score
					call backup4sidechain(0, chainID, ic, temp_group, aa_backup)
					m_best=m
				endif
			endif
		enddo
				
		if(m_best==0) goto 20
		call backup4sidechain(1, chainID, ic, tgroup, aa_backup)

		flag=1
20		continue
		deallocate(tgroup_para)
		deallocate(temp_group)
	endif
	deallocate(aa_group)
	
	return
	
	end subroutine sequence_mutation_4scmfloop

	subroutine scmf_loop(group)
	implicit none
	integer							:: ic, m, flag, flag1, flag2
	integer							:: rotanum, chainID, ip, feedback
	character*4						:: aminoacid_name
	type(groupdetails)				:: group(repeated_unit,gnum), tgroup(repeated_unit,gnum), temp_group(repeated_unit,gnum), aa_group(repeated_unit,40)
	character*4						:: group_name(3)
	
	tgroup=group
    
	do ic=1, gnum	
		flag=0
		do chainID=2, repeated_unit
			if(tgroup(chainID,ic)%gtype==tgroup(1,ic)%gtype) then
				flag=flag+1
				goto 10
			else
				call groupinfo(tgroup(1,ic)%gtype, group_name, flag1)
				aminoacid_name=group_name(flag1)
				call sequence_mutation_4scmfloop(chainID, ic, tgroup, aminoacid_name, temp_group, feedback)
				if(feedback==1) then
					tgroup=temp_group
					flag=flag+1
					goto 20
				endif
			endif
20			continue
		tgroup=temp_group

		enddo
10		continue
		if(flag==0) then
			open(5, file="error.txt", access="append")
				write(5,*) "scmf_loop can't find all rotamers"
			close(5)
			stop
		endif
	enddo
	group=tgroup

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
	type(groupdetails), dimension(:,:), allocatable &
									:: temp_group, aa_group
	type(energyparameters), dimension(:,:), allocatable &
									:: tgroup_para

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
				allocate(aa_group(repeated_unit,40))
				call findrotamer(ic, group, aminoacid_name, rotanum, aa_group, ip)
				grade=aa_lib(ip)%grade

				allocate(temp_group(repeated_unit,gnum))
				allocate(tgroup_para(repeated_unit,gnum))

				obs=0
				matrix=0.0
				do m=1, rotanum
					call residue_replace(ii, ic, group, m, aa_group, temp_group)
					if(m==1) then
						call energy_parameter(temp_group, tgroup_para)
						call atom_links4sidechain(ii, ic, temp_group, S_numex, S_inb, S_numex4, S_inb4)
						call atom_links(temp_group, W_numex, W_inb, W_numex4, W_inb4)
					endif

					stage=0
					call sidechain_optimization(stage, ii, ic, temp_group, tgroup_para, S_numex, S_inb, S_numex4, S_inb4, Dihedral4entropy)

					if(stage==1) then
						obs=obs+1
						do j=1, grade
							matrix(obs,j)=Dihedral4entropy(j)								
						enddo
					endif
				enddo

				call entropy_calculation(aminoacid_name,rotanum,matrix,obs,grade,entropy)
				entropy4individual(ii,ic)=entropy
		
				deallocate(tgroup_para)
				deallocate(temp_group)
				deallocate(aa_group)
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
	character*4						:: group_name_1(3), group_name_2(3), rest_AA(6)

	type(groupdetails)				:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), tgroup(repeated_unit,gnum)

	do attempt=1, gnum
		call ran_gen(ran2,0)
		if(ran2.le.scmfswitch1) then ! mutate one AA
            feedback_4=0
            do while(feedback_4 == 0)
				call pickupsite(0, ic_1)
				call mc_choose_aminoacid(ic_1, group, aminoacid_name_1)
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
					if (group_name_2(1) == rest_AA(i)) goto 96 ! 用AA不含N C端的名称比较，如果挑选的非nmr AA是剩余list中的一个, 则直接交换
                enddo
                goto 19 ! 如果不是则重新选择aa
            else if(ic_2_NMR .and. .not. ic_1_NMR) then
                call find_rest_aa(group, rest_AA, l)
                do i=1, l
					if (group_name_1(1) == rest_AA(i)) goto 96 ! 如果挑选的非nmr AA是剩余list中的一个, 则直接交换
                enddo
				goto 19 ! 如果不是则重新选择aa
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


	subroutine sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old)
	implicit none
	integer							:: step, attempt, Num4attempt, ic_1, ic_2, categoryID_1, categoryID_2, chainID_1, chainID_2 ! ic_1, 2 是随机取的sites
	integer							:: i, ii, l, j, feedback_1, feedback_2, feedback_3, feedback_4, flag1, flag2
	real							:: ran2
	real(kind=8)                    :: score_old, score_new, binding_energy_old, binding_energy_new
	real							:: entropy_old, entropy_new
	real							:: score4hydration_old, score4hydration_new
	real							:: Pagg_old, Pagg_new
	real							:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum), temp_entropy4individual(repeated_unit,gnum)
	character*4						:: aminoacid_name_1, aminoacid_name_2
    character*20					:: accpt
	character*4						:: group_name_1(3), group_name_2(3), rest_AA(6)
	logical							:: ic_1_NMR, ic_2_NMR
    character(len=:), allocatable   :: pep_name
	type(groupdetails)				:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), tgroup(repeated_unit,gnum)
	type(energyparameters), dimension(:,:), allocatable &
									:: group_para
	integer, dimension(:), allocatable &
									:: W_numex, W_numex4
	integer, dimension(:,:), allocatable &
									:: W_inb, W_inb4
	
    feedback_2=0
    feedback_3=0
    !penalty=0
	do attempt=1, 7
		call ran_gen(ran2,0)
		if(ran2.le.scmfswitch1) then
            feedback_4=0
            do while(feedback_4 == 0)
				call pickupsite(0, ic_1)
				call mc_choose_aminoacid(ic_1, group, aminoacid_name_1)
                call check_restrained_aminoacid(group, ic_1, aminoacid_name_1, feedback_4)
            enddo

			call sequence_mutation(ic_1, group, entropy4individual, aminoacid_name_1, tgroup, Tentropy4individual, feedback_1)

			if(feedback_1==1) then
				allocate(group_para(repeated_unit,gnum))
				allocate(W_numex(repeated_unit*atom_num)); allocate(W_numex4(repeated_unit*atom_num))
				allocate(W_inb(repeated_unit*atom_num,20)); allocate(W_inb4(repeated_unit*atom_num,60))
				
				call atom_links(tgroup, W_numex, W_inb, W_numex4, W_inb4)
				call energy_parameter(tgroup, group_para)

				call bindingenergy(tgroup, group_para, Tentropy4individual, W_numex, W_inb, W_numex4, W_inb4, score_new, binding_energy_new, entropy_new, score4hydration_new, Pagg_new)
                !MC 方法 判断是否接受新的变化 如果分数更低就接受，如果分数变高，取随机数和exp(-E/kT)比较，是否接受
				call MC_technique_sequence(score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old, group, entropy4individual, score_new, binding_energy_new, entropy_new, score4hydration_new, Pagg_new, tgroup, Tentropy4individual, feedback_3)

				deallocate(W_inb); deallocate(W_inb4)
				deallocate(W_numex); deallocate(W_numex4)
				deallocate(group_para)
			endif
			
			if(feedback_3==1) then
				accpt="Accept"
			elseif(feedback_1==1 .and. feedback_3==0) then
				accpt="Reject-MC"
			!PENALTY = PENALTY + 1                
			elseif(feedback_1==0) then
				accpt="Reject-Rotamer"
			endif

			open(3, file="energydetails.txt",  access="append")
				write(3,4) step, attempt, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old, "Not-SCMF ", accpt
				write(3,"(<gnum>(a5))") (group(1,j)%gtype, j=1, gnum)
				!write(3,*) "*******************************"
			close(3)
        else !交换amino acid，随机取两个sites，site 1的amino acid mutate到site 2上，site 2到1上
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
					if (group_name_2(1) == rest_AA(i)) goto 96 ! 用AA不含N C端的名称比较，如果挑选的非nmr AA是剩余list中的一个, 则直接交换
                enddo
                goto 19 ! 如果不是则重新选择aa
            else if(ic_2_NMR .and. .not. ic_1_NMR) then
                call find_rest_aa(group, rest_AA, l)
                do i=1, l
					if (group_name_1(1) == rest_AA(i)) goto 96 ! 如果挑选的非nmr AA是剩余list中的一个, 则直接交换
                enddo
				goto 19 ! 如果不是则重新选择aa
            endif               
96          continue
            
            aminoacid_name_1=group_name_2(flag1)
			aminoacid_name_2=group_name_1(flag2)
                 
			call sequence_mutation(ic_1, group, entropy4individual, aminoacid_name_1, temp_group, temp_entropy4individual, feedback_1)
			
			if(feedback_1==1) then

				call sequence_mutation(ic_2, temp_group, temp_entropy4individual, aminoacid_name_2, tgroup, Tentropy4individual, feedback_2)

				if(feedback_2==1) then
					allocate(group_para(repeated_unit,gnum))
					allocate(W_numex(repeated_unit*atom_num)); allocate(W_numex4(repeated_unit*atom_num))
					allocate(W_inb(repeated_unit*atom_num,20)); allocate(W_inb4(repeated_unit*atom_num,60))
					
					call atom_links(tgroup, W_numex, W_inb, W_numex4, W_inb4)
					call energy_parameter(tgroup, group_para)
		
					call bindingenergy(tgroup, group_para, Tentropy4individual, W_numex, W_inb, W_numex4, W_inb4, score_new, binding_energy_new, entropy_new, score4hydration_new, Pagg_new)
					call MC_technique_sequence(score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old, group, entropy4individual, score_new, binding_energy_new, entropy_new, score4hydration_new, Pagg_new, tgroup, Tentropy4individual, feedback_3)

					deallocate(W_inb); deallocate(W_inb4)
					deallocate(W_numex); deallocate(W_numex4)
					deallocate(group_para)
				endif
			endif

			if(feedback_3==1) then
				accpt="Accept" ! accpt,输出新的peptide和score
			elseif(feedback_2==1 .and. feedback_3==0) then
				accpt="Reject-MC" ! reject,输出旧的peptide和score
			elseif(feedback_2==0) then
				accpt="Reject-Rotamer"
			endif

			open(3, file="energydetails.txt",  access="append")
				write(3,4) step, attempt, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old, "intra-SCMF ", accpt
				write(3,"(<gnum>(a5))") (group(1,j)%gtype, j=1, gnum)
				!write(3,*) "*******************************"
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
4	format(i7,i7,5f20.13,a15,a15)

	return
	end subroutine sequence_optimization
	
	
	subroutine sheet_optimization(group, step, entropy4individual, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old)
	implicit none
	integer						:: step, attempt, ip1, ip2, ic_1, feedback_1, feedback_2, feedback_3
	integer						:: flag1, i, j, flag4conformer  
	real						:: ran2, rmsd
	real(kind=8)                :: score_old, score_new
	real(kind=8)                :: binding_energy_old, binding_energy_new
	real						:: entropy_old, entropy_new
	real						:: score4hydration_old, score4hydration_new
	real						:: Pagg_old, Pagg_new
	real						:: entropy4individual(repeated_unit,gnum), Tentropy4individual(repeated_unit,gnum), temp_entropy4individual(repeated_unit,gnum)
	character*20				:: sheet_change, accpt
	character*4					:: aminoacid_name_1
	character*4					:: group_name_1(3)
	character(len=:), allocatable   :: pep_name
	type(groupdetails)			:: group(repeated_unit,gnum), temp_group(repeated_unit,gnum), Tgroup(repeated_unit,gnum)
	
    type(energyparameters), dimension(:,:), allocatable :: group_para
	integer, dimension(:), allocatable :: W_numex, W_numex4
	integer, dimension(:,:), allocatable :: W_inb, W_inb4

!*********************************************************************************
!Here ip1 decides whether the sheet moves in +x, -x, +y, -y, +z and -z direction 
!ip1=1 moves sheet in x and -x
!ip1=2 moves sheet in y and -y
!ip1=3 moves sheet in z and -z
!data passing: group -> (move sheet) -> Tgroup -> (repacking) -> temp_group  
!*********************************************************************************
	flag4conformer=0
	do attempt=1, 7
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
		do ic_1=1, gnum
			do j=1, void_site_num
				if(ic_1.ge.void_site_start(j).and.ic_1.le.void_site_end(j)) then 
				goto 15
			endif
			enddo
			call groupinfo(group(1,ic_1)%gtype, group_name_1, flag1)
			aminoacid_name_1=group_name_1(flag1)
			call sequence_mutation(ic_1, Tgroup, entropy4individual, aminoacid_name_1, temp_group, temp_entropy4individual, feedback_1)
			if(feedback_1==1) then 
				goto 10
10				continue
			else
				goto 20
			endif
15			continue
		enddo
!********************************************				
		allocate(group_para(repeated_unit,gnum))
		allocate(W_numex(repeated_unit*atom_num)); allocate(W_numex4(repeated_unit*atom_num))
		allocate(W_inb(repeated_unit*atom_num,20)); allocate(W_inb4(repeated_unit*atom_num,60))
				
		call atom_links(temp_group, W_numex, W_inb, W_numex4, W_inb4)
		call energy_parameter(temp_group, group_para)
		
        call bindingenergy(temp_group, group_para, temp_entropy4individual, W_numex, W_inb, W_numex4, W_inb4, score_new, binding_energy_new, entropy_new, score4hydration_new, Pagg_new)
		call MC_technique_sheet(score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old, Tgroup, entropy4individual, score_new, binding_energy_new, entropy_new, score4hydration_new, Pagg_new, temp_group, temp_entropy4individual, feedback_2)

		if(feedback_2==1) then
			flag4sheet=1
		endif
			
		deallocate(W_inb); deallocate(W_inb4)
		deallocate(W_numex); deallocate(W_numex4)
		deallocate(group_para)
	
		group=Tgroup
			
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
			
		open(3, file="energydetails.txt",  access="append")
		write(3,4) step, attempt, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old, sheet_change, ip2, accpt
		write(3,"(<gnum>(a5))") (group(1,j)%gtype, j=1, gnum)
		!write(3,*) "*******************************"
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
	
4	format(i7,i7,5f20.13,a15,i7,a15)
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
	integer							:: step, sub_circle, i, j, k
	real							:: ran2, rmsd_x, rmsd_y, rmsd_z, rmsd
	real(kind=8)                    :: score_old, binding_energy_old
    real							:: entropy_old, score4hydration_old, Pagg_old
    !real							:: entropy4individual(repeated_unit,gnum)
	!type(groupdetails)				:: group(repeated_unit,gnum)
    character(len=:), allocatable   :: pep_name
    real, dimension(:,:), allocatable :: entropy4individual
	type(groupdetails), dimension(:,:), allocatable :: group
	type(groupdetails), dimension(:,:), allocatable :: temp_group
	type(energyparameters), dimension(:,:), allocatable :: group_para
	integer, dimension(:), allocatable :: W_numex, W_numex4
	integer, dimension(:,:), allocatable :: W_inb, W_inb4
	character*4						:: aminoacid_name
    
	call inputfile
    call setparameter
    call readdir
    
    allocate(character(len=gnum) :: pep_name)
    allocate(entropy4individual(repeated_unit,gnum))
    allocate(group(repeated_unit,gnum))
    
    call init_random_seed()
	call readpdb(group) ! 读取peptide的坐标，种类等等
	call rotamerlib  
    
	if(recalcu_switch==0) then
		allocate(temp_group(repeated_unit,gnum))
		call scmf_substitution(group, sub_circle, temp_group)
		if(sub_circle.ne.0) then
			group=temp_group
        endif
        
        call replace_extra_restricted_aa(group, temp_group)
        group=temp_group
        deallocate(temp_group)
        
		call scmf_loop(group)
        
		do step=1, 5
			call sequence_optimization_nonthermal(group)
        enddo

		call generatepdb(0, 0, group)
		energy_min=0.0
		
		allocate(group_para(repeated_unit,gnum))
		allocate(W_numex(repeated_unit*atom_num)); allocate(W_numex4(repeated_unit*atom_num))
		allocate(W_inb(repeated_unit*atom_num,20)); allocate(W_inb4(repeated_unit*atom_num,60))

		call energy_parameter(group, group_para)
		call atom_links(group, W_numex, W_inb, W_numex4, W_inb4)
	
		call initial_entropy_individual(group, entropy4individual)
		call bindingenergy(group, group_para, entropy4individual, W_numex, W_inb, W_numex4, W_inb4, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old)

	
		deallocate(W_inb); deallocate(W_inb4)
		deallocate(W_numex); deallocate(W_numex4)
		deallocate(group_para)
	
        call convert_AA_name(group, pep_name)
		open(5, file="energyprofile.txt", access="append")
			write(5,6) 0,  " ", pep_name, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old
		close(5)
        
        open(5,file="rmsd.txt",access="append")	! add rmsd header
			write(5,*) "step", " rmsd_x", " rmsd_y", " rmsd_z", " rmsd"
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
	endif

    !!!!这里是主体循环，进行sequence_optimization and sheet_optimization
	do step=nstep_start, nstep_terminal
		if(step.le.500) then
            ekt=ekt_sequence
			call sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old)
		else
			if(flag4sheet==0) then
				call ran_gen(ran2,0)
				if(ran2.ge.sheet_switch) then
                    ekt=ekt_sequence
					call sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old)
                else ! sheet_opt内会将flag4shee设为1  
					call sheet_optimization(group, step, entropy4individual, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old)
				endif
            else
                !ekt=ekt_sequence * (sheetmove_interval - flag4sheet) / sheetmove_interval * 1.5 + 0.5
                ekt=ekt_sequence
				call sequence_optimization(group, step, entropy4individual, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old)
				flag4sheet=flag4sheet+1
				if(flag4sheet.gt. sheetmove_interval) then
					flag4sheet=0
                end if
            endif				
        endif
        
        open(5,file="rmsd.txt",access="append")
            call rmsd_calc(1, group, rmsd_x)
            call rmsd_calc(2, group, rmsd_y)
            call rmsd_calc(3, group, rmsd_z)
            call rmsd_calc(4, group, rmsd)
            write(5,*) step, rmsd_x, rmsd_y, rmsd_z, rmsd
        close(5)
        
        call convert_AA_name(group, pep_name)
		open(5, file="energyprofile.txt", access="append")
			write(5,6) step,  " ", pep_name, score_old, binding_energy_old, entropy_old, score4hydration_old, Pagg_old
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
6	format(i5,2a,5f10.4)

end
	