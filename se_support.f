! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHout ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module se_support

  use star_lib
  use star_def
  use const_def
  use chem_def
  use HDF5 ! This module contains all necessary modules

  implicit none

  character(len=*), parameter :: head_cycle_name = 'cycle'          ! head of a cycle group name
  character(len=*), parameter :: dataset_cycle_name  = 'SE_DATASET'  ! Name of the dataset within a cycle
  character(len=*), parameter :: dataset_cycle_name2 = 'ISO_MASSF'   ! Name of the dataset within a cycle
  character(len=*), parameter :: format_file = "(I7.7)"          ! Format for the file  name
  character(len=*), parameter :: format_cycle = "(I10)"          ! Format for the cycle name
  character(len=80), parameter :: seversion = "2.0"              ! SE Version number

  character(len=256) :: filename, prefix
  integer            :: num_mod_seoutput        
  integer(HID_T)     :: file_id                   ! File identifier 
  integer            :: error                     ! Error flag
  integer            :: firstmodel



contains
    
  
      
  subroutine se_startup(s, restart, use_hdf5_output,                     &
       se_codev, se_modname, se_prefix, se_num_mod_output,               &
       log_names, log_vals, log_is_int, num_log_columns,                 &
       profile_names, profile_vals, profile_is_int, num_profile_columns, &
       ierr)

    type (star_info), pointer :: s
    logical, intent(in) :: restart, use_hdf5_output
    character (len=256), intent(in) :: se_codev, se_modname, se_prefix
    integer, intent(in) :: se_num_mod_output
    character (len=maxlen_history_column_name), intent(in), pointer :: log_names(:) ! (num_log_columns)
    double precision, intent(in), pointer :: log_vals(:) ! (num_log_columns)
    logical, intent(in), pointer :: log_is_int(:) ! true if the values in the log column are integers
    integer, intent(in) :: num_log_columns
    character (len=maxlen_profile_column_name), intent(in), pointer :: profile_names(:) ! num_profiles_columns
    double precision, intent(in), pointer :: profile_vals(:,:) ! (nz,num_profile_columns)
    logical, intent(in), pointer :: profile_is_int(:) ! true if the values in the profile column are integers
    integer, intent(in)  :: num_profile_columns
    integer, intent(out) :: ierr
  
    

    ierr = 0
    firstmodel = s% model_number + 1 
    
    if (use_hdf5_output) then
       call start_new_se_output_file(s,se_codev,se_modname,se_prefix,se_num_mod_output,    &
            log_names,log_vals,log_is_int,num_log_columns,                                 &
            profile_names,profile_vals,profile_is_int,num_profile_columns,firstmodel,ierr)
    endif
  end subroutine se_startup
      
      

  subroutine start_new_se_output_file(s,se_codev,se_modname,se_prefix,se_num_mod_output,  &
       log_names,log_vals,log_is_int,num_log_columns,                                     &
       profile_names,profile_vals,profile_is_int,num_profile_columns,model,ierr)
    
    type (star_info), pointer :: s
    character (len=256), intent(in) :: se_codev, se_modname, se_prefix
    integer, intent(in) :: se_num_mod_output
    character (len=maxlen_history_column_name), intent(in), pointer :: log_names(:)
    double precision, intent(in), pointer :: log_vals(:)
    logical, intent(in), pointer :: log_is_int(:)
    integer, intent(in) :: num_log_columns
    character (len=maxlen_profile_column_name), intent(in), pointer :: profile_names(:)
    double precision, intent(in), pointer :: profile_vals(:,:)
    logical, intent(in), pointer :: profile_is_int(:)
    integer, intent(in) :: num_profile_columns, model
    integer, intent(out) :: ierr

    character (len=80) :: char_modname
    
    ierr = 0

    write(*,*) "initializing se output ..."
    write(*,*) "Starting new se output file for model ", model
    write(char_modname,format_file) model
    prefix=se_modname(1:len_trim(se_modname))//'/'//                                            &
         se_prefix(1:len_trim(se_prefix))//'.'//char_modname(1:len_trim(char_modname))//'.se.h5'
    filename = prefix(1:len_trim(prefix))
    num_mod_seoutput = se_num_mod_output


    ! initialize FORTRAN interface.
    call h5open_f(error)
    
    ! Create a new file using default properties.
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

    ! Write header
    call sewritefhead(file_id,s,model,se_codev,se_modname)

    ! Terminate access to the file.
    call h5fclose_f(file_id, error)
    
    ! Close FORTRAN interface.
    call h5close_f(error)
    
  end subroutine start_new_se_output_file
      
      
      
  integer function se_finish_step(s, use_hdf5_output, use_se_names,    &
       se_codev, se_modname, se_prefix, se_num_mod_output,             &
       how_many_extra_history_columns, data_for_extra_history_columns, &
       how_many_extra_profile_columns, data_for_extra_profile_columns, &
       log_names, log_vals, log_is_int, num_log_columns,               &
       profile_names, profile_vals, profile_is_int, num_profile_columns)

    type (star_info), pointer :: s
    logical, intent(in) :: use_hdf5_output, use_se_names
    character (len=256), intent(in) :: se_codev, se_modname, se_prefix
    integer, intent(in) :: se_num_mod_output
    character (len=maxlen_history_column_name), intent(in), pointer :: log_names(:) ! (num_log_columns)
    double precision, intent(in), pointer :: log_vals(:) ! (num_log_columns)
    logical, intent(in), pointer :: log_is_int(:) ! true if the values in the profile column are integers
    integer, intent(in) :: num_log_columns
    character (len=maxlen_profile_column_name), intent(in), pointer :: profile_names(:) ! num_profiles_columns
    double precision, intent(in), pointer :: profile_vals(:,:) ! (nz,num_profile_columns)
    logical, intent(in), pointer :: profile_is_int(:) ! true if the values in the profile column are integers
    integer, intent(in) :: num_profile_columns

    interface
       include 'extra_history_cols.inc'
       include 'extra_profile_cols.inc'
    end interface
         
    integer :: modini
    integer :: ierr

    se_finish_step = keep_going
         
    if (.not. use_hdf5_output) return
    
    modini = s% model_number
    if (modulo(modini-firstmodel,num_mod_seoutput) == 0) then
       call start_new_se_output_file(s,se_codev,se_modname,se_prefix,se_num_mod_output, &
            log_names,log_vals,log_is_int,num_log_columns,                              &
            profile_names,profile_vals,profile_is_int,num_profile_columns,modini,ierr)
       if (ierr /= 0) then
          write(*,*) 'extras_finish_step:start_new_se_output_file' 
          stop ' could not open new se file' 
       end if
    end if
    
    ! initialize FORTRAN interface.
    call h5open_f(error)         
    
    ! Open existing file
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
         
    call sewritecycle(file_id, s,s% model_number, use_se_names,          &
         log_names, log_vals, log_is_int, num_log_columns,               &
         profile_names, profile_vals, profile_is_int, num_profile_columns)
    
    ! Terminate access to the file.
    call h5fclose_f(file_id, error)
    
    ! Close FORTRAN interface.
    call h5close_f(error)
    
  end function se_finish_step



  subroutine sewritefhead(file_id,s,model,codev,modname)

    type (star_info), pointer        :: s
    integer(HID_T),   intent(in)     :: file_id           ! File identifier 
    integer,          intent(in)     :: model
    character(len=*), intent(in)     :: modname, codev
    
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
    character(len=80) :: hdf5version
    character(len=16) :: format_hdf5

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
    call scalar_to_rank1_int(num_mod_seoutput, value_int)
    call add_attribute_int(file_id, 1, "icyclenb", value_int)
    call scalar_to_rank1_int(model, value_int)
    call add_attribute_int(file_id, 1, "firstcycle", value_int)  
    ! Add character attributes
    call scalar_to_rank1_chr(codev, value_chr)
    call add_attribute_chr(file_id, 1, "codev", value_chr)
    call scalar_to_rank1_chr(modname, value_chr)
    call add_attribute_chr(file_id, 1, "modname", value_chr)
    ! Version numbers
    ! HDF5
    call h5get_libversion_f(majnum, minnum, relnum, error)
    format_hdf5 = '(i1,a1,i1,a1,i2)'
    if (minnum > 9) then
       format_hdf5 = '(i1,a1,i2,a1,i2)'
    end if
    write(hdf5version,format_hdf5) majnum,'.',minnum,'.',relnum
    call scalar_to_rank1_chr(hdf5version, value_chr)
    call add_attribute_chr(file_id, 1, "HDF5_version", value_chr)
    ! SE
    call scalar_to_rank1_chr(seversion, value_chr)
    call add_attribute_chr(file_id, 1, "SE_version", value_chr)
    ! Add double attributes
    call scalar_to_rank1_dbl(s% initial_mass, value_dbl)
    call add_attribute_dbl(file_id, 1, "mini", value_dbl)
    call scalar_to_rank1_dbl(s% initial_z, value_dbl)
    call add_attribute_dbl(file_id, 1, "zini", value_dbl)
    call scalar_to_rank1_dbl(0.d0, value_dbl)
    call add_attribute_dbl(file_id, 1, "rotini", value_dbl)
    call scalar_to_rank1_dbl(s% overshoot_f_above_burn_h, value_dbl)
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

  end subroutine sewritefhead

  
      
  subroutine sewritecycle(file_id,s,icounthdf, use_se_names,            &
       log_names, log_vals, log_is_int, num_log_columns,                &
       profile_names, profile_vals, profile_is_int, num_profile_columns)

    type (star_info), pointer       :: s
    integer, intent(in)             :: icounthdf        ! Model number
    integer(HID_T), intent(in)      :: file_id          ! File identifier
    logical, intent(in) :: use_se_names
    character (len=maxlen_history_column_name), pointer, intent(in) :: log_names(:) ! (num_columns)
    double precision, intent(in), pointer:: log_vals(:) ! (num_columns)
    logical, intent(in), pointer :: log_is_int(:) ! true if the values in the profile column are integers
    integer :: num_log_columns
    character (len=maxlen_profile_column_name), intent(in), pointer :: profile_names(:) ! num_profiles_columns
    double precision, intent(in), pointer :: profile_vals(:,:) ! (nz,num_profile_columns)
    logical, intent(in), pointer :: profile_is_int(:) ! true if the values in the profile column are integers
    integer, intent(in) :: num_profile_columns

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
    call h5gcreate_f(file_id, groupname, group_id, error)

    ! Add attributes
    do i=1, num_log_columns       
       if (use_se_names) then
          call mesa2se_history(log_names(i))
       endif

       if (log_is_int(i)) then
          call scalar_to_rank1_int(int(log_vals(i)), value_int)
          call add_attribute_int(group_id, 1, log_names(i), value_int)
       else
          call scalar_to_rank1_dbl(log_vals(i), value_dbl)
          call add_attribute_dbl(group_id, 1, log_names(i), value_dbl)
       endif
    enddo

    ! Add dataset
    ! First do SE_DATASET
    do i=1, num_profile_columns
       if (use_se_names) then
          call mesa2se_profile(profile_names(i))
       endif
    enddo

    call add_dataset(group_id, dataset_cycle_name, s% nz, &
         profile_is_int, profile_names, profile_vals, num_profile_columns)

    ! Now do iso_massf
    number_of_species = s% species

    allocate(iso_massf(s% nz,number_of_species))
    allocate(iso_is_int(number_of_species))
    allocate(iso_names(number_of_species))

    do i=1, number_of_species
       iso_is_int(i) = .false.
       write(iso_str,*) i
       iso_names(i) = iso_str
       k = s% net_iso(s% chem_id(i))
       do j=1, s% nz
          iso_massf(j,i) = s% xa(k,j)
       enddo
    enddo

    call add_dataset(group_id, dataset_cycle_name2, s% nz, &
         iso_is_int, iso_names, iso_massf, number_of_species)

    deallocate(iso_massf)
    deallocate(iso_is_int)
    deallocate(iso_names)
   
    ! Close group
    call h5gclose_f(group_id, error)
    
  end subroutine sewritecycle
      

      
  subroutine se_after_evolve(s, ierr)
    type (star_info), pointer :: s
    integer, intent(out) :: ierr
    ierr = 0
  end subroutine se_after_evolve



  subroutine scalar_to_rank1_int(scalar,rank1)

    integer, intent(in)                :: scalar ! Scalar to be transformed in an array
    integer, intent(out), dimension(1) :: rank1  ! Result
         
    rank1(1) = scalar
    
  end subroutine scalar_to_rank1_int



  subroutine scalar_to_rank1_dbl(scalar,rank1)

    double precision, intent(in)                :: scalar ! Scalar to be transformed in an array
    double precision, intent(out), dimension(1) :: rank1  ! Result
         
    rank1(1) = scalar

  end subroutine scalar_to_rank1_dbl



  subroutine scalar_to_rank1_chr(scalar,rank1)
    
    character(len=80), intent(in)                :: scalar ! Scalar to be transformed in an array
    character(len=80), intent(out), dimension(1) :: rank1  ! Result
    
    rank1(1) = scalar
    
  end subroutine scalar_to_rank1_chr



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

  
  subroutine mesa2se_history(str)

    ! Change the name of a cycle attribute
    character(len=maxlen_history_column_name) :: str
    
    if (str(1:len_trim(str)) .eq. 'radius') then
       str = 'R_sol'
    elseif (str(1:len_trim(str)) .eq. 'photosphere_L') then
       str = 'L_photosphere_Lsun'
    elseif (str(1:len_trim(str)) .eq. 'effective_T') then
       str = 'Teff'
    elseif (str(1:len_trim(str)) .eq. 'surface_h1') then
       str = 'X_surface_h1'
    elseif (str(1:len_trim(str)) .eq. 'surface_he4') then
       str = 'X_surface_he4'
    elseif (str(1:len_trim(str)) .eq. 'surface_c12') then
       str = 'X_surface_c12'
    elseif (str(1:len_trim(str)) .eq. 'surface_o16') then
       str = 'X_surface_o16'
    elseif (str(1:len_trim(str)) .eq. 'star_age') then
       str = 'age'
    elseif (str(1:len_trim(str)) .eq. 'time_step') then
       str = 'deltat'
    elseif (str(1:len_trim(str)) .eq. 'max_eps_h_m') then
       str = 'eps_h_max_m'
    elseif (str(1:len_trim(str)) .eq. 'max_eps_he_m') then
       str = 'eps_he_max_m'
    elseif (str(1:len_trim(str)) .eq. 'log_L') then
       str = 'logL'
    elseif (str(1:len_trim(str)) .eq. 'log_Teff') then
       str = 'logTeff'
    elseif (str(1:len_trim(str)) .eq. 'num_zones') then
       str = 'shellnb'
    elseif (str(1:len_trim(str)) .eq. 'star_mass') then
       str = 'total_mass'
    endif
       
  end subroutine mesa2se_history


  subroutine mesa2se_profile(str)

    ! Change the name of a cycle profile
    character(len=maxlen_profile_column_name) :: str
    
    !if (str(1:len_trim(str)) .eq. 'mixing_type ') then
    !   str = 'convection_indicator'
    !endif
       
  end subroutine mesa2se_profile


end module se_support
