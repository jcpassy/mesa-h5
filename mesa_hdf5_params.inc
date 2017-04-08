!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! These are the parameters that define the location/names/size of the HDF5 file !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! hdf5_prefix         -> basename for each HDF5 file
!!!!!!!                        default: M3.00000Z0.020 for a 3 Msun star with Z = 0.02
!!!!!!! hdf5_codev          -> some extra string attribute for each HDF5 file
!!!!!!!                        default: mesa rev `revision number`
!!!!!!! hdf5_modname        -> folder in which to write out the HDF5 files
!!!!!!! head_cycle_name     -> basename for each cycle
!!!!!!! dataset_cycle_name  -> name of the dataset containing the profiles
!!!!!!! dataset_cycle_name2 -> name of the dataset containing the abundances
!!!!!!! format_file         -> suffix format in the HDF5 name
!!!!!!! format_cycle        -> suffix format in the cycle name
!!!!!!! hdf5_num_mod_output -> number of cycles per HDF5 file
!!!!!!! change_names        -> whether you want to change the name of some quantities,
!!!!!!!			       for instance 'star_mass' into 'stellar_mass'.
!!!!!!!			       See subroutines change_history_names and change_profile_names
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! If left as empty strings, these parameters will take default values 
      ! determined at runtime (see header).
      character (len=256) :: hdf5_prefix  = ''
      character (len=256) :: hdf5_codev   = ''
      character(len=256), parameter :: hdf5_modname      = 'HDF5'

      character(len=*), parameter :: head_cycle_name     = 'cycle' 
      character(len=*), parameter :: dataset_cycle_name  = 'SE_DATASET'
      character(len=*), parameter :: dataset_cycle_name2 = 'ISO_MASSF'

      character(len=*), parameter :: format_file  = "(I7.7)"
      character(len=*), parameter :: format_cycle = "(I10)"

      integer, parameter :: hdf5_num_mod_output = 50

      logical, parameter :: change_names = .true.
