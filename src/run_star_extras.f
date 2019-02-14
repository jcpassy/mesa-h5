! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib
      use chem_def
      use run_star_support, only: failed
      use utils_lib, only: utils_OMP_GET_MAX_THREADS
      use HDF5

      implicit none

      include 'mesa_hdf5_params.inc'

      character(len=*), parameter :: format_file  = "(I7.7)"
      character(len=*), parameter :: format_cycle = "(I10)"

      integer :: num_history_columns, num_profile_columns
      character (len=maxlen_profile_column_name), pointer :: profile_names(:) ! num_profiles_columns
      double precision, pointer :: profile_vals(:,:) ! (nz,num_profile_columns)
      logical, pointer :: profile_is_int(:) ! true if the values in the profile column are integers

      ! These variables are just for the logs to check the profile names for the HDF5 file
      integer :: num_profile_columns_logs
      character (len=maxlen_profile_column_name), pointer :: profile_names_logs(:) ! num_profiles_columns_logs
      double precision, pointer :: profile_vals_logs(:,:) ! (nz,num_profile_columns)
      logical, pointer :: profile_is_int_logs(:) ! true if the values in the profile column are integers

      character(len=256) :: filename, prefix
      integer(HID_T)     :: file_id                   ! File identifier
      integer            :: error                     ! Error flag
      integer            :: firstmodel

    ! these routines are called by the standard run_star check_model
    contains

      ! HDF5 subroutines
      include 'mesa_hdf5_routines.inc'

      subroutine extras_controls(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        ! this is the place to set any procedure pointers you want to change
        ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

        ! Uncomment these lines if you wish to use the functions in this file,
        ! otherwise we use a null_ version which does nothing.
        s% extras_startup => extras_startup
        !s% extras_check_model => extras_check_model
        s% extras_finish_step => extras_finish_step
        !s% extras_after_evolve => extras_after_evolve
        s% how_many_extra_history_columns => how_many_extra_history_columns
        s% data_for_extra_history_columns => data_for_extra_history_columns
        s% how_many_extra_profile_columns => how_many_extra_profile_columns
        s% data_for_extra_profile_columns => data_for_extra_profile_columns

        ! Once you have set the function pointers you want,
        ! then uncomment this (or set it in your star_job inlist)
        ! to disable the printed warning message,
        s% job% warn_run_star_extras =.false.

      end subroutine extras_controls

      ! None of the following functions are called unless you set their
      ! function point in extras_control.

      integer function extras_startup(id, restart, ierr)
        integer, intent(in) :: id
        logical, intent(in) :: restart
        integer, intent(out) :: ierr

        integer :: i,id_extra, revfile
        character(len=256) :: command, rev, mass_str, met_str, tmp_str
        double precision :: mass_f

        integer :: ios, line, index
        character(len=256) :: buffer, name

        type (star_info), pointer :: s
        ierr = 0
        revfile = 41

        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        extras_startup = 0
        if (.not. restart) then
           call alloc_extra_info(s)
        else ! it is a restart
           call unpack_extra_info(s)
        end if

        ! Set up global parameters if needed
        ! hdf5_codev
        if (len_trim(hdf5_codev) == 0) then
           command = 'svn info '//mesa_dir(1:len_trim(mesa_dir))//' --show-item revision > rev.txt'
           ! Tidious way to get the revision number
           call system(command)
           open(unit=revfile, file='rev.txt', status="old")
           read(revfile,*,iostat=ierr) rev
           hdf5_codev = "mesa rev "//rev(1:len_trim(rev))
           close(revfile)
           call system('rm rev.txt')
        end if

        ! hdf5_prefix
        if (len_trim(hdf5_prefix) == 0) then
           ! These formats are REALLY such a pain...
           ! May be there is a better way to do it?
           write(mass_str,"(1pd26.5)") s% initial_mass
           call str_to_double(mass_str, mass_f, ierr)
           write(mass_str,"(f10.5)") mass_f
           write(met_str, "(f4.3)") s% initial_z
           tmp_str = "M"//mass_str(1:len_trim(mass_str))//"Z0"//met_str(1:len_trim(met_str))
           call remove_white_spaces(tmp_str, hdf5_prefix)
        end if

        ! Create HDF5 folder
        call system("mkdir -p "//hdf5_modname)

        ! Get names and values for history
        call get_data_for_history_columns(s, id_extra, ierr)
        num_history_columns = s% number_of_history_columns

        ! Get names and values for profiles
        ! We do not use the same for the profiles files and the HDF5
        ! HDF5 profiles are a subset of the profiles files
        num_profile_columns_logs = num_standard_profile_columns(s) + how_many_extra_profile_columns(id, id_extra)
        allocate(                                               &
             profile_names_logs(num_profile_columns_logs),      &
             profile_vals_logs(s% nz,num_profile_columns_logs), &
             profile_is_int_logs(num_profile_columns_logs),     &
             stat=ierr)

        call get_data_for_profile_columns(s, id_extra, s% nz, &
         profile_names_logs, profile_vals_logs, profile_is_int_logs,ierr)

        ! Parse the hdf5_profile_columns.list to get number of HDF5 profiles
        open(41, file=hdf5_profile_list, iostat=ios, status='old')
        line = 0
        ! Read once to get number of lines
        do while (ios == 0)
          read(41, '(A)', iostat=ios) buffer
          call extract_profile_name(buffer, name)
          if ((ios == 0) .and. (len_trim(name) .gt. 0)) then
            line = line + 1
          end if
        end do
        close(41)
        num_profile_columns = line

        ! Parse the hdf5_profile_columns.list to extract the names
        allocate(                                     &
             profile_names(num_profile_columns),      &
             profile_vals(s% nz,num_profile_columns), &
             profile_is_int(num_profile_columns),     &
             stat=ierr)
        open(41, file=hdf5_profile_list, iostat=ios, status='old')
        line = 0
        do while (ios == 0)
          read(41, '(A)', iostat=ios) buffer
          call extract_profile_name(buffer, name)
          if ((ios == 0) .and. (len_trim(name) .gt. 0)) then
            line = line + 1
            profile_names(line) = name
          end if
        end do
        close(41)

        ! Copy names and get is_int
        do line = 1, num_profile_columns
          call get_profile_index(line, index, ierr)
          if (ierr /= 0) return
          profile_is_int(line) = profile_is_int_logs(index)
        end do

        ! Call HDF5 routine
        call hdf5_startup(s, restart,                                          &
             hdf5_codev, hdf5_modname, hdf5_prefix, hdf5_num_mod_output,       &
             ierr)
        if (failed('hdf5_startup', ierr)) return

        ! Deallocate only profile_vals (because nz changes during the evolution)
        deallocate(profile_vals)

      end function extras_startup

      subroutine get_profile_index(index, index_logs, ierr)

        integer, intent(in) :: index
        integer, intent(out) :: index_logs, ierr

        ierr = 0

        ! look for index in profile_names_logs corresponding profile_names(i)
        do index_logs=1, num_profile_columns_logs
          if (profile_names(index) .eq. profile_names_logs(index_logs)) then
            goto 123
          end if
        end do
        write(*,*) 'Error: the following profile has not been not found in profile_names_logs:'
        write(*,*) profile_names(index)
        ierr = 412
123     return

      end subroutine get_profile_index

      subroutine get_data_for_profile_columns_hfd5(s, id_extra, nz, ierr)

        type (star_info), pointer :: s
        integer, intent(in) :: id_extra, nz
        integer, intent(out) :: ierr

        integer i, j, index

        ierr = 0

        ! Call for the logs
        call get_data_for_profile_columns(s, id_extra, nz, &
            profile_names_logs, profile_vals_logs, profile_is_int_logs, ierr)
        if (ierr /= 0) return

        ! Get corresponding index
        do i=1, num_profile_columns
          call get_profile_index(i, index, ierr)
          if (ierr /= 0) return

          ! Copy
          do j=1, nz
            profile_vals(j,i) = profile_vals_logs(j, index)
          end do
        end do

      end subroutine get_data_for_profile_columns_hfd5

      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
        integer, intent(in) :: id, id_extra
        integer :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        extras_check_model = keep_going
        if (.false. .and. s% star_mass_h1 < 0.35d0) then
           ! stop when star hydrogen mass drops to specified level
           extras_check_model = terminate
           write(*, *) 'have reached desired hydrogen mass'
           return
        end if

        ! if you want to check multiple conditions, it can be useful
        ! to set a different termination code depending on which
        ! condition was triggered.  MESA provides 9 customizeable
        ! termination codes, named t_xtra1 .. t_xtra9.  You can
        ! customize the messages that will be printed upon exit by
        ! setting the corresponding termination_code_str value.
        ! termination_code_str(t_xtra1) = 'my termination condition'

        ! by default, indicate where (in the code) MESA terminated
        if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
        integer, intent(in) :: id, id_extra
        integer :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        how_many_extra_history_columns = 2
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
        integer, intent(in) :: id, id_extra, n
        character (len=maxlen_history_column_name) :: names(n)
        real(dp) :: vals(n)
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        !note: do NOT add the extras names to history_columns.list
        ! the history_columns.list is only for the built-in log column options.
        ! it must not include the new column names you are adding here.

        names(1) = 'L_spec'
        names(2) = 'diff_L'

        vals(1) = ((s% Teff)**4.0 / s% grav(1)) / Lsun
        vals(2) = vals(1) / (s% photosphere_L)

      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id, id_extra)
        use star_def, only: star_info
        integer, intent(in) :: id, id_extra
        integer :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        how_many_extra_profile_columns = 3
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
        use star_def, only: star_info, maxlen_profile_column_name
        use const_def, only: dp
        integer, intent(in) :: id, id_extra, n, nz
        character (len=maxlen_profile_column_name) :: names(n)
        real(dp) :: vals(nz,n)
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        integer :: k
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        !note: do NOT add the extra names to profile_columns.list
        ! the profile_columns.list is only for the built-in profile column options.
        ! it must not include the new column names you are adding here.

        ! here is an example for adding a profile column
        !if (n /= 1) stop 'data_for_extra_profile_columns'
        !names(1) = 'beta'
        !do k = 1, nz
        !   vals(k,1) = s% Pgas(k)/s% P(k)
        !end do

        names(1) = 'delta_mass'
        names(2) = 'rho'
        names(3) = 'dcoeff'
        do k = 1, nz
          vals(k,1) = s% mstar * s% dq(k) / msol
          vals(k,2) = s% rho(k)
          vals(k,3) = s% D_mix(k)
        end do

      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
        integer, intent(in) :: id, id_extra
        integer :: ierr
        integer :: ios, line
        character(len=100) :: buffer

        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return
        extras_finish_step = keep_going
        call store_extra_info(s)

        ! to save a profile,
           ! s% need_to_save_profiles_now = .true.
        ! to update the star log,
           ! s% need_to_update_history_now = .true.

        ! Get names and values for history
        call get_data_for_history_columns(s, id_extra, ierr)

        allocate(profile_vals_logs(s% nz,num_profile_columns_logs), &
                 profile_vals(s% nz,num_profile_columns),           &
                 stat=ierr)
        call get_data_for_profile_columns_hfd5(s, id_extra, s% nz, ierr)

        ! Main function call
        call hdf5_finish_step(s, change_names,                                                       &
             hdf5_codev, hdf5_modname, hdf5_prefix, hdf5_num_mod_output,                             &
             how_many_extra_history_columns, data_for_extra_history_columns,                         &
             how_many_extra_profile_columns, data_for_extra_profile_columns,                         &
             s% history_names, s% history_values, s% history_value_is_integer, num_history_columns,  &
             profile_names, profile_vals, profile_is_int, num_profile_columns,                       &
             ierr)

        deallocate(profile_vals_logs, profile_vals)
        if (failed('hdf5_finish_step', ierr)) return

        ! see extras_check_model for information about custom termination codes
        ! by default, indicate where (in the code) MESA terminated
        if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step


      subroutine extras_after_evolve(id, id_extra, ierr)
        integer, intent(in) :: id, id_extra
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        ierr = 0
        call star_ptr(id, s, ierr)

        ! Deallocate everything but profile_vals
        deallocate(profile_names_logs, profile_vals_logs, profile_is_int_logs)
        deallocate(profile_names, profile_is_int)
        if (ierr /= 0) return

      end subroutine extras_after_evolve


      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3


      subroutine alloc_extra_info(s)
        integer, parameter :: extra_info_alloc = 1
        type (star_info), pointer :: s
        call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
        integer, parameter :: extra_info_get = 2
        type (star_info), pointer :: s
        call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
        integer, parameter :: extra_info_put = 3
        type (star_info), pointer :: s
        call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
        integer, parameter :: extra_info_alloc = 1
        integer, parameter :: extra_info_get = 2
        integer, parameter :: extra_info_put = 3
        type (star_info), pointer :: s
        integer, intent(in) :: op

        integer :: i, j, num_ints, num_dbls, ierr

        i = 0
        ! call move_int or move_flg
        num_ints = i

        i = 0
        ! call move_dbl

        num_dbls = i

        if (op /= extra_info_alloc) return
        if (num_ints == 0 .and. num_dbls == 0) return

        ierr = 0
        call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
        if (ierr /= 0) then
           write(*,*) 'failed in star_alloc_extras'
           write(*,*) 'alloc_extras num_ints', num_ints
           write(*,*) 'alloc_extras num_dbls', num_dbls
           stop 1
        end if

      contains

        subroutine move_dbl(dbl)
          real(dp) :: dbl
          i = i+1
          select case (op)
          case (extra_info_get)
             dbl = s% extra_work(i)
          case (extra_info_put)
             s% extra_work(i) = dbl
          end select
        end subroutine move_dbl

        subroutine move_int(int)
          integer :: int
          i = i+1
          select case (op)
          case (extra_info_get)
             int = s% extra_iwork(i)
          case (extra_info_put)
             s% extra_iwork(i) = int
          end select
        end subroutine move_int

        subroutine move_flg(flg)
          logical :: flg
          i = i+1
          select case (op)
          case (extra_info_get)
             flg = (s% extra_iwork(i) /= 0)
          case (extra_info_put)
             if (flg) then
                s% extra_iwork(i) = 1
             else
                s% extra_iwork(i) = 0
             end if
          end select
        end subroutine move_flg

      end subroutine move_extra_info

    end module run_star_extras
