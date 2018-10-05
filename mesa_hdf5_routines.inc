!======================================================================================================================!


      subroutine hdf5_startup(s, restart,                                    &
           hdf5_codev, hdf5_modname, hdf5_prefix, hdf5_num_mod_output,       &
           ierr)
        
        type (star_info), pointer :: s
        logical, intent(in) :: restart
        character (len=256), intent(in) :: hdf5_codev, hdf5_modname, hdf5_prefix
        integer, intent(in) :: hdf5_num_mod_output
        integer, intent(out) :: ierr

        ierr = 0
        firstmodel = s% model_number + 1

        call start_new_hdf5_output_file(s, hdf5_codev, hdf5_modname, hdf5_prefix, hdf5_num_mod_output, &
             firstmodel,ierr)
        if (ierr /= 0) return

      end subroutine hdf5_startup


!======================================================================================================================!


      subroutine hdf5_finish_step(s, change_names,                              &
           hdf5_codev, hdf5_modname, hdf5_prefix, hdf5_num_mod_output,          &
           how_many_extra_history_columns, data_for_extra_history_columns,      &
           how_many_extra_profile_columns, data_for_extra_profile_columns,      &
           history_names, history_vals, history_is_int, num_history_columns,    &
           profile_names, profile_vals, profile_is_int, num_profile_columns,    &
           ierr)
        
        type (star_info), pointer :: s
        logical, intent(in) :: change_names
        character (len=256), intent(in) :: hdf5_codev, hdf5_modname, hdf5_prefix
        integer, intent(in) :: hdf5_num_mod_output
        character (len=maxlen_history_column_name), intent(in), pointer :: history_names(:) ! (num_history_columns)
        double precision, intent(in), pointer :: history_vals(:) ! (num_history_columns)
        logical, intent(in), pointer :: history_is_int(:) ! true if the values in the profile column are integers
        integer, intent(in) :: num_history_columns
        character (len=maxlen_profile_column_name), intent(in), pointer :: profile_names(:) ! num_profiles_columns
        double precision, intent(in), pointer :: profile_vals(:,:) ! (nz,num_profile_columns)
        logical, intent(in), pointer :: profile_is_int(:) ! true if the values in the profile column are integers
        integer, intent(in) :: num_profile_columns
        integer, intent(out) :: ierr

        interface
           include 'extra_history_cols.inc'
           include 'extra_profile_cols.inc'
        end interface

        integer :: modini
        
        ierr = 0
        modini = s% model_number

        if (modulo(modini-firstmodel,hdf5_num_mod_output) == 0) then
           call start_new_hdf5_output_file(s, hdf5_codev, hdf5_modname, hdf5_prefix, hdf5_num_mod_output,  &
                modini, ierr)
           if (ierr /= 0) return
        end if
    
        ! initialize FORTRAN interface.
        call h5open_f(ierr)
        if (ierr /= 0) return

        ! Open existing file
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr)
        if (ierr /= 0) return

        call hdf5_writecycle(file_id, s,s% model_number, change_names,            &
             history_names, history_vals, history_is_int, num_history_columns,    &
             profile_names, profile_vals, profile_is_int, num_profile_columns,    &
             ierr)
        if (ierr /= 0) return

        ! Terminate access to the file.
        call h5fclose_f(file_id, ierr)
        if (ierr /= 0) return

        ! Close FORTRAN interface.
        call h5close_f(ierr)
        if (ierr /= 0) return

      end subroutine hdf5_finish_step


!======================================================================================================================!


      subroutine start_new_hdf5_output_file(s, hdf5_codev, hdf5_modname, hdf5_prefix, &
           hdf5_num_mod_output, model, ierr)

        type (star_info), pointer :: s
        character (len=256), intent(in) :: hdf5_codev, hdf5_modname, hdf5_prefix
        integer, intent(in) :: hdf5_num_mod_output
        integer, intent(in) :: model
        integer, intent(out) :: ierr
        
        character (len=80) :: char_modname
    
        ierr = 0

        write(*,*) "initializing HDF5 output ..."
        write(*,*) "Starting new HDF5 output file for model ", model

        write(char_modname,format_file) model
        prefix=hdf5_modname(1:len_trim(hdf5_modname))//'/'//  &
             hdf5_prefix(1:len_trim(hdf5_prefix))//'.'//char_modname(1:len_trim(char_modname))//'.h5'
        filename = prefix(1:len_trim(prefix))

        call h5open_f(ierr)
        if (ierr /= 0) return

        ! Create a new file using default properties.
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)
        if (ierr /= 0) return
        
        ! Write header
        call hdf5_writefhead(file_id, s, model, hdf5_codev, hdf5_modname, ierr)
        if (ierr /= 0) return
    
        ! Terminate access to the file.
        call h5fclose_f(file_id, ierr)
        if (ierr /= 0) return
    
        ! Close FORTRAN interface.
        call h5close_f(ierr)
        if (ierr /= 0) return

      end subroutine start_new_hdf5_output_file


!======================================================================================================================!


      subroutine hdf5_writefhead(file_id, s, model, codev, modname, ierr)

        type (star_info), pointer        :: s
        integer(HID_T),   intent(in)     :: file_id           ! File identifier
        integer,          intent(in)     :: model
        character(len=*), intent(in)     :: modname, codev
        integer, intent(out) :: ierr

        integer,           dimension(1) :: value_int  ! Attribute int
        character(len=80), dimension(1) :: value_chr  ! Attribute character
        double precision,  dimension(1) :: value_dbl  ! Attribute double

        character(len=80),          pointer :: prof_name(:)   ! Memory datatype names
        double precision,           pointer :: prof_vals(:,:) ! Memory datatype types

        integer                             :: number_of_species, myiso
        integer, pointer                    :: isotoprint(:)
        logical, pointer                    :: is_integer(:)
        double precision, pointer           :: iso_charge(:), iso_mass(:), iso_meric_state(:)
      
        integer i,majnum, minnum, relnum
        character(len=80) :: hdf5version,time
        character(len=16) :: mesaversion,format_hdf5

        ierr = 0
    
        number_of_species = s% species
        allocate(isotoprint(number_of_species))
        allocate(iso_charge(number_of_species))
        allocate(iso_mass(number_of_species))
        allocate(iso_meric_state(number_of_species))
        
        do i=1, number_of_species
           isotoprint(i) = s% chem_id(i)
           myiso = isotoprint(i)
           iso_charge(i)      = chem_isos% Z(myiso)
           iso_mass(i)        = chem_isos% Z_plus_N(myiso)
           iso_meric_state(i) = 1
        end do
        
        !! ATTRIBUTES
        ! Add integer attributes
        call scalar_to_rank1_int(hdf5_num_mod_output, value_int)
        call add_attribute_int(file_id, 1, "icyclenb", value_int)
        call scalar_to_rank1_int(model, value_int)
        call add_attribute_int(file_id, 1, "firstcycle", value_int)
	call scalar_to_rank1_int(utils_OMP_GET_MAX_THREADS(), value_int)		
	call add_attribute_int(file_id, 1, "OMP_NUM_THREADS", value_int)
	call scalar_to_rank1_int(version_number, value_int)		
	call add_attribute_int(file_id, 1, "MESA_version", value_int)

        ! Add character attributes
        call scalar_to_rank1_chr(codev, value_chr)
        call add_attribute_chr(file_id, 1, "codev", value_chr)
        call scalar_to_rank1_chr(modname, value_chr)
        call add_attribute_chr(file_id, 1, "modname", value_chr)
	call scalar_to_rank1_chr(mesa_dir, value_chr)	
	call add_attribute_chr(file_id, 1, "MESA_dir", value_chr)

        ! HDF5 version
        call h5get_libversion_f(majnum, minnum, relnum, error)
        format_hdf5 = '(i1,a1,i1,a1,i2)'
        if (minnum > 9) then
           format_hdf5 = '(i1,a1,i2,a1,i2)'
        end if
        write(hdf5version,format_hdf5) majnum,'.',minnum,'.',relnum
        call scalar_to_rank1_chr(hdf5version, value_chr)
        call add_attribute_chr(file_id, 1, "HDF5_version", value_chr)

        ! Time
        call when(time)
        call scalar_to_rank1_chr(time, value_chr)
        call add_attribute_chr(file_id, 1, "Creation", value_chr)

        ! Add double attributes
        call scalar_to_rank1_dbl(s% initial_mass, value_dbl)
        call add_attribute_dbl(file_id, 1, "mini", value_dbl)
        call scalar_to_rank1_dbl(s% initial_z, value_dbl)
        call add_attribute_dbl(file_id, 1, "zini", value_dbl)
        call scalar_to_rank1_dbl(0.d0, value_dbl)
        call add_attribute_dbl(file_id, 1, "rotini", value_dbl)
        call scalar_to_rank1_dbl(s% overshoot_f_above_burn_h_core, value_dbl)
        call add_attribute_dbl(file_id, 1, "overini", value_dbl)
        call scalar_to_rank1_dbl(secyer, value_dbl)
        call add_attribute_dbl(file_id, 1, "age_unit", value_dbl)
        call scalar_to_rank1_dbl(secyer, value_dbl)
        call add_attribute_dbl(file_id, 1, "one_year", value_dbl)
        call scalar_to_rank1_dbl(msol, value_dbl)
        call add_attribute_dbl(file_id, 1, "mass_unit", value_dbl)
        call scalar_to_rank1_dbl(rsol, value_dbl)
        call add_attribute_dbl(file_id, 1, "radius_unit", value_dbl)
        call scalar_to_rank1_dbl(1.d0, value_dbl)
        call add_attribute_dbl(file_id, 1, "rho_unit", value_dbl)
        call scalar_to_rank1_dbl(1.d0, value_dbl)
        call add_attribute_dbl(file_id, 1, "temperature_unit", value_dbl)
        call scalar_to_rank1_dbl(1.d0, value_dbl)
        call add_attribute_dbl(file_id, 1, "dcoeff_unit", value_dbl)
        
        !! DATASETS
        allocate(is_integer(1))
        allocate(prof_name(1))
        allocate(prof_vals(number_of_species,1))
        prof_name = 'data'

        ! isomeric_state (int)
        is_integer = .true.
        do i=1,number_of_species
           prof_vals(i,1) = iso_meric_state(i)
        enddo
        call add_dataset(file_id, "isomeric_state", number_of_species, &
             is_integer, prof_name, prof_vals, 1)
        
        ! A (double)
        is_integer = .false.
        do i=1,number_of_species
           prof_vals(i,1) = iso_mass(i)
        enddo
        call add_dataset(file_id, "A", number_of_species, &
             is_integer, prof_name, prof_vals, 1)
        
        ! Z (double)
        is_integer = .false.
        do i=1,number_of_species
           prof_vals(i,1) = iso_charge(i)
        enddo
        call add_dataset(file_id, "Z", number_of_species, &
             is_integer, prof_name, prof_vals, 1)
        
        deallocate(is_integer)
        deallocate(prof_name)
        deallocate(prof_vals)
        deallocate(isotoprint)
        deallocate(iso_charge)
        deallocate(iso_mass)
        deallocate(iso_meric_state)
        
      end subroutine hdf5_writefhead


!======================================================================================================================!
      

      subroutine hdf5_writecycle(file_id,s,icounthdf, change_names,            &
           history_names, history_vals, history_is_int, num_history_columns,   &
           profile_names, profile_vals, profile_is_int, num_profile_columns,   &
           ierr)
        
        type (star_info), pointer       :: s
        integer, intent(in)             :: icounthdf        ! Model number
        integer(HID_T), intent(in)      :: file_id          ! File identifier
        logical, intent(in) :: change_names
        character (len=maxlen_history_column_name), pointer, intent(in) :: history_names(:) ! (num_columns)
        double precision, intent(in), pointer:: history_vals(:) ! (num_columns)
        logical, intent(in), pointer :: history_is_int(:) ! true if the values in the profile column are integers
        integer :: num_history_columns
        character (len=maxlen_profile_column_name), intent(in), pointer :: profile_names(:) ! num_profiles_columns
        double precision, intent(in), pointer :: profile_vals(:,:) ! (nz,num_profile_columns)
        logical, intent(in), pointer :: profile_is_int(:) ! true if the values in the profile column are integers
        integer, intent(in) :: num_profile_columns
        integer, intent(out) :: ierr

        character (len=maxlen_profile_column_name), pointer :: iso_names(:) ! number_of_species
        logical, pointer :: iso_is_int(:) ! true if the values in the profile column are integers
        double precision, pointer :: iso_massf(:,:) !(nz, number_of_species)
        character(len=80) :: iso_str
        
        character(len=80)               :: groupname        ! Group name
        integer(HID_T)                  :: group_id         ! Group identifier
        character(len=80)               :: str_cycle        ! Cycle in string format
        integer                         :: str_len_cycle    ! length
        integer,           dimension(1) :: value_int        ! Attribute int
        double precision,  dimension(1) :: value_dbl        ! Attribute double

        integer           :: number_of_species
        integer i,j,k

        integer, allocatable :: mask(:) ! for masking profile columns
        logical, allocatable :: mask_profile_is_int(mask)
        character (len=maxlen_profile_column_name), allocatable :: &
           mask_profile_names(:)
        double precision, allocatable :: mask_profile_vals(:,:)


        ierr = 0

        ! Name of the group (cycle)
        write(str_cycle,format_cycle) icounthdf
        groupname = head_cycle_name//str_cycle
        
        ! Replace blanks with zeros
        str_len_cycle = len_trim(groupname)
        do i=1, str_len_cycle
           if (groupname(i:i) == ' ') then
              groupname(i:i) = '0'
           endif
        enddo
        
        call h5gcreate_f(file_id, groupname, group_id, ierr)
        if (ierr /= 0) return

        ! Add attributes
        do i=1, num_history_columns
           if (change_names) then
              call change_history_names(history_names(i))
           endif
           
           if (history_is_int(i)) then
              call scalar_to_rank1_int(int(history_vals(i)), value_int)
              call add_attribute_int(group_id, 1, history_names(i), value_int)
           else
              call scalar_to_rank1_dbl(history_vals(i), value_dbl)
              call add_attribute_dbl(group_id, 1, history_names(i), value_dbl)
           endif
        enddo
        
        ! Add dataset
        ! First do SE_DATASET
        do i=1, num_profile_columns
           if (change_names) then
              call change_profile_names(profile_names(i))
           endif
        enddo
        
        !call add_dataset(group_id, dataset_cycle_name, s% nz, &
        !     profile_is_int, profile_names, profile_vals, num_profile_columns)

        ! -_- only add what is necessary to the HDF5 files for post-processing
        !
        ! this is (6 columns): dcoeff, delta_mass, temperature, radius, mass, rho
        allocate( &
           mask(6), &
           mask_profile_names(6), &
           mask_profile_is_int(6), &
           mask_profile_vals(s% nz, 6) &
           )

        j = 0
        do i = 1, num_profile_columns
           select case(profile_names(i))
              case ("dcoeff","delta_mass","temperature","radius","mass","rho")
                 j += 1
                 mask(j) = i ! we need this column
              case default
                 continue ! we don't need this column for post-processing
           end select
        end do

        mask_profile_names  = profile_names(mask)
        mask_profile_is_int = profile_is_int(mask)
        mask_profile_vals   = profile_vals(:,mask)
  
        call add_dataset(group_id, dataset_cycle_name, s% nz, &
           mask_profile_is_int, mask_profile_names, &
           mask_profile_vals, 6)

        deallocate(mask, mask_profile_names, mask_profile_is_int, &
           mask_profile_vals)
  
        ! -_- erase things for writing iso_massf - we don't need those for
        ! post-processing

        number_of_species = s% species
        
        ! Close group
        call h5gclose_f(group_id, ierr)
        if (ierr /= 0) return

      end subroutine hdf5_writecycle        


!======================================================================================================================!


      subroutine scalar_to_rank1_int(scalar,rank1)
        integer, intent(in)                :: scalar ! Scalar to be transformed in an array
        integer, intent(out), dimension(1) :: rank1  ! Result
        rank1(1) = scalar
      end subroutine scalar_to_rank1_int


!======================================================================================================================!


      subroutine scalar_to_rank1_dbl(scalar,rank1)
        double precision, intent(in)                :: scalar ! Scalar to be transformed in an array
        double precision, intent(out), dimension(1) :: rank1  ! Result        
        rank1(1) = scalar
      end subroutine scalar_to_rank1_dbl


!======================================================================================================================!


      subroutine scalar_to_rank1_chr(scalar,rank1)
        character(len=80), intent(in)                :: scalar ! Scalar to be transformed in an array
        character(len=80), intent(out), dimension(1) :: rank1  ! Result
        rank1(1) = scalar
      end subroutine scalar_to_rank1_chr


!======================================================================================================================!


      subroutine add_attribute_int(target_id,size,name,value)
        
        integer(HID_T),    intent(in)        :: target_id      ! Target (file or group) identifier
        integer,           intent(in)        :: size           ! Attribute size
        character(len=*),  intent(in)        :: name           ! Attribute name
        integer, intent(in), dimension(size) :: value          ! Attribute value

        integer                              :: rank = 1       ! Attribure rank
        integer(HSIZE_T), dimension(1)       :: dims = (/1/)   ! Attribute dimension
        integer(HID_T)                       :: space_id       ! Attribute Dataspace identifier
        integer(HID_T)                       :: type_id        ! Attribute Dataspace identifier
        integer(HID_T)                       :: id             ! Attribute identifier

        dims(1) = size
        ! Create scalar data space for the attribute.
        call h5screate_simple_f(rank, dims, space_id, error)

        ! Create datatype for the attribute.
        call h5tcopy_f(H5T_NATIVE_integer, type_id, error)
        
        ! Create dataset attribute.
        call h5acreate_f(target_id, name, type_id, space_id, id, error)
        
        ! Write the attribute data.
        call h5awrite_f(id, type_id, value, dims, error)
        
        ! Close the attribute.
        call h5aclose_f(id, error)
        
        ! Terminate access to the data space.
        call h5sclose_f(space_id, error)

      end subroutine add_attribute_int


!======================================================================================================================!


      subroutine add_attribute_dbl(target_id,size,name,value)

        integer(HID_T),    intent(in)        :: target_id      ! Target (file or group) identifier
        integer,           intent(in)        :: size           ! Attribute size
        character(len=*),  intent(in)        :: name           ! Attribute name
        double precision,  intent(in), dimension(size) :: value          ! Attribute value

        integer                              :: rank = 1       ! Attribure rank
        integer(HSIZE_T), dimension(1)       :: dims = (/1/)   ! Attribute dimension
        integer(HID_T)                       :: space_id       ! Attribute Dataspace identifier
        integer(HID_T)                       :: type_id        ! Attribute Dataspace identifier
        integer(HID_T)                       :: id             ! Attribute identifier

        dims(1) = size
        ! Create scalar data space for the attribute.
        call h5screate_simple_f(rank, dims, space_id, error)

        ! Create datatype for the attribute.
        call h5tcopy_f(H5T_NATIVE_DOUBLE, type_id, error)

        ! Create dataset attribute.
        call h5acreate_f(target_id, name, type_id, space_id, id, error)

        ! Write the attribute data.
        call h5awrite_f(id, type_id, value, dims, error)

        ! Close the attribute.
        call h5aclose_f(id, error)

        ! Terminate access to the data space.
        call h5sclose_f(space_id, error)
        
      end subroutine add_attribute_dbl
        

!======================================================================================================================!


      subroutine add_attribute_chr(target_id,size,name,value)
       
        integer(HID_T),    intent(in)        :: target_id      ! Target (file or group) identifier
        integer,           intent(in)        :: size           ! Attribute size
        character(len=*),  intent(in)        :: name           ! Attribute name
        character(len=*),  intent(in), dimension(size) :: value          ! Attribute value
        
        integer                              :: rank = 1       ! Attribure rank
        integer(HSIZE_T), dimension(1)       :: dims = (/1/)   ! Attribute dimension
        integer(HID_T)                       :: space_id       ! Attribute Dataspace identifier
        integer(HID_T)                       :: type_id        ! Attribute Dataspace identifier
        integer(HID_T)                       :: id             ! Attribute identifier
        
        integer(SIZE_T), PARAMETER           :: len = 80       ! length of the attribute
        
        dims(1) = size
        ! Create scalar data space for the attribute.
        call h5screate_simple_f(rank, dims, space_id, error)
        
        ! Create datatype for the attribute.
        call h5tcopy_f(H5T_NATIVE_character, type_id, error)
        call h5tset_size_f(type_id, len, error)
        
        ! Create dataset attribute.
        call h5acreate_f(target_id, name, type_id, space_id, id, error)
        
        ! Write the attribute data.
        call h5awrite_f(id, type_id, value, dims, error)
        
        ! Close the attribute.
        call h5aclose_f(id, error)
        
       ! Terminate access to the data space.
        call h5sclose_f(space_id, error)
        
      end subroutine add_attribute_chr


!======================================================================================================================!

      
      subroutine add_dataset(target_id, dset_name, n_points, &
           profile_is_int, profile_names, profile_vals, num_profile_columns)
        
        integer(HID_T),    intent(in) :: target_id        ! Target (file or group) identifier
        character(len=*),  intent(in) :: dset_name        ! Dataset name
        integer,           intent(in) :: n_points         ! Number of mesh points (size of the dataset)
        
        character (len=maxlen_profile_column_name), intent(in), pointer :: profile_names(:) ! num_profiles_columns
        double precision, intent(in), pointer :: profile_vals(:,:) ! (nz,num_profile_columns)
        logical, intent(in), pointer :: profile_is_int(:) ! true if the values in the profile column are integers
        integer :: num_profile_columns
        integer :: number_of_species
        integer, pointer :: array_int(:)

        integer(HID_T) :: id          ! Dataset identifier
        integer(HID_T) :: space_id    ! Dataspace identifier
        integer(HID_T) :: plist_id    ! Dataset transfer property
        integer(HID_T) :: type_id     ! Compound datatype identifier
        integer        :: rank = 1    ! Dataset rank

        integer(SIZE_T) :: type_size_total  ! Size of the entire datatype
        integer(SIZE_T) :: type_sized       ! Size of the double precision datatype
        integer(SIZE_T) :: type_sizei       ! Size of the integer datatype
        integer(SIZE_T) :: type_sizec       ! Size of the character datatype
        integer(SIZE_T) :: offset           ! Member's offset

        integer(SIZE_T), dimension(num_profile_columns+1) :: profile_size     ! Memory datatype sizes
        integer(SIZE_T), dimension(num_profile_columns+1) :: profile_offset   ! Memory datatype offsets
        integer(HID_T),  dimension(num_profile_columns+1) :: profile_id       ! Memory datatype identifier

        integer(HSIZE_T), dimension(1) :: dims
        integer(HSIZE_T), dimension(1) :: data_dims
        
        integer :: i

        ! initialization
        ! First calculate total size by calculating sizes of each member
        call h5tget_size_f(H5T_NATIVE_DOUBLE, type_sized, error)
        call h5tget_size_f(H5T_NATIVE_integer, type_sizei, error)
        call h5tget_size_f(H5T_NATIVE_character, type_sizec, error)

        type_size_total = 0
        profile_offset(1) = 0

        do i=1, num_profile_columns
           if (profile_is_int(i)) then
              profile_size(i) = type_sizei
           else
              profile_size(i) = type_sized
           endif
           if (i>1) then
              profile_offset(i) = profile_offset(i-1) + profile_size(i-1)
           endif
           type_size_total = type_size_total + profile_size(i)
        enddo
        
        ! Dims
        dims(1) = n_points
        data_dims(1) = n_points

        ! Set dataset transfer property to preserve partially initialized fields
        ! during write/read to/from dataset with compound datatype.
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_preserve_f(plist_id, .TRUE., error)

        ! Create the dataspace
        call h5screate_simple_f(rank, dims, space_id, error)

        ! Create compound datatype.
        call h5tcreate_f(H5T_COMPOUND_F, type_size_total, type_id, error)

        ! insert members
        do i=1, num_profile_columns
           if (profile_is_int(i)) then
              call h5tinsert_f(type_id, profile_names(i), profile_offset(i), H5T_NATIVE_integer, error)
           else
              call h5tinsert_f(type_id, profile_names(i), profile_offset(i), H5T_NATIVE_DOUBLE, error)
           endif
        enddo

        ! Create the dataset with compound datatype.
        call h5dcreate_f(target_id, dset_name, type_id, space_id, id, error)

        ! Create memory types. We have to create a compound datatype
        ! for each member we want to write.
        offset = 0
        ! integers
        do i=1, num_profile_columns
           call h5tcreate_f(H5T_COMPOUND_F, profile_size(i), profile_id(i), error)
           if (profile_is_int(i)) then
              call h5tinsert_f(profile_id(i), profile_names(i), offset, H5T_NATIVE_integer, error)
           else
              call h5tinsert_f(profile_id(i), profile_names(i), offset, H5T_NATIVE_DOUBLE, error)
           endif
        enddo
    
        ! Write data by fields in the datatype. Fields order is not important.
        ! integers
        do i=1, num_profile_columns
           if (profile_is_int(i)) then
              call array_double2int(profile_vals(:,i), n_points, array_int)
              call h5dwrite_f(id, profile_id(i), array_int, data_dims, error, xfer_prp = plist_id)
           else
              call h5dwrite_f(id, profile_id(i), profile_vals(:,i), data_dims, error, xfer_prp = plist_id)
           endif
        enddo

        ! End access to the dataset and release resources used by it.
        call h5dclose_f(id, error)
        
        ! Terminate access to the data space.
        call h5sclose_f(space_id, error)
        
        ! Terminate access to the datatype
        call h5tclose_f(type_id, error)
        do i=1, num_profile_columns
           call h5tclose_f(profile_id(i), error)
        enddo

      end subroutine add_dataset


!======================================================================================================================!


      subroutine array_double2int(array_dbl,n,array_int)

        integer, intent(in) :: n
        double precision, intent(in), dimension(n) :: array_dbl
        integer, intent(out), pointer :: array_int(:)
        integer :: i
        
        nullify(array_int)
        allocate(array_int(n))

        do i=1, n
           array_int(i) = int(array_dbl(i))
        enddo
        
      end subroutine array_double2int        


!======================================================================================================================!


      subroutine change_history_names(str)

        ! Change the name of a cycle attribute
	! For instance, to change 'star_mass' into 'stellar_mass':
	! if (str(1:len_trim(str)) .eq. 'star_mass') str = 'stellar_mass'

        character(len=maxlen_history_column_name) :: str

        if (str(1:len_trim(str)) .eq. 'radius')        str = 'R_sol'
        if (str(1:len_trim(str)) .eq. 'photosphere_L') str = 'L_photosphere_Lsun'
        if (str(1:len_trim(str)) .eq. 'effective_T')   str = 'Teff'
        if (str(1:len_trim(str)) .eq. 'surface_h1')    str = 'X_surface_h1'
        if (str(1:len_trim(str)) .eq. 'surface_he4')   str = 'X_surface_he4'
        if (str(1:len_trim(str)) .eq. 'surface_c12')   str = 'X_surface_c12'
        if (str(1:len_trim(str)) .eq. 'surface_o16')   str = 'X_surface_o16'
        if (str(1:len_trim(str)) .eq. 'star_age')      str = 'age'
        if (str(1:len_trim(str)) .eq. 'time_step')     str = 'deltat'
        if (str(1:len_trim(str)) .eq. 'max_eps_h_m')   str = 'eps_h_max_m'
        if (str(1:len_trim(str)) .eq. 'max_eps_he_m')  str = 'eps_he_max_m'
        if (str(1:len_trim(str)) .eq. 'log_L')         str = 'logL'
        if (str(1:len_trim(str)) .eq. 'log_Teff')      str = 'logTeff'
        if (str(1:len_trim(str)) .eq. 'num_zones')     str = 'shellnb'
        if (str(1:len_trim(str)) .eq. 'star_mass')     str = 'total_mass'
        
      end subroutine change_history_names


!======================================================================================================================!


      subroutine change_profile_names(str)        

        ! Change the name of a cycle profile
	! For instance, to change 'mixing_type' into 'convection_indicator':
        ! if (str(1:len_trim(str)) .eq. 'mixing_type ') str = 'convection_indicator'

        character(len=maxlen_profile_column_name) :: str
        
      end subroutine change_profile_names


!======================================================================================================================!

        
      subroutine when(time)

        character(len=80), intent(out) :: time
        integer today(3), now(3)

        call idate(today)   ! today(1)=day, (2)=month, (3)=year
        call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
        write (time, 1000)  today(2), today(1), today(3), now

1000    format ( i2.2, '/', i2.2, '/', i4.4, ' at ',  &
             i2.2, ':', i2.2, ':', i2.2 )

      end subroutine when


!======================================================================================================================!

