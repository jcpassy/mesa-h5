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
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************
 
      module run_star_support

      use star_lib
      use star_def
      use chem_def
      use chem_lib
      use eos_lib
      use net_lib
      use atm_lib
      use rates_lib, only: set_which_rate_1212
      use rates_def, only: maxlen_reaction_Name

      implicit none
      
      character (len=1000) :: mesa_dir, &
         eosDT_cache_dir, &
         eosPT_cache_dir, &
         ionization_cache_dir, &
         kap_cache_dir, &
         rates_cache_dir, &
         weaklib_cache_dir

      logical :: pause_before_terminate

      character (len=1000) :: profile_columns_file, history_columns_file
      
      logical :: show_log_description_at_start
      logical :: show_net_reactions_info, show_net_species_info
      
      integer :: first_model_for_timing
      
      logical :: show_eqns_and_vars_names
      
      logical :: pgstar_flag
      
      logical :: use_se_output, use_se_names
      character (len=256) :: se_codev
      character (len=256) :: se_modname
      character (len=256) :: se_prefix
      integer :: se_num_mod_output

      logical :: save_photo_when_terminate
      
      logical :: create_pre_main_sequence_model
      integer :: pre_ms_relax_num_steps
      real(dp) :: pre_ms_T_c, pre_ms_guess_rho_c, pre_ms_d_log10_P
      real(dp) :: pre_ms_logT_surf_limit, pre_ms_logP_surf_limit
      
      logical :: create_initial_model
      integer :: initial_model_relax_num_steps
      real(dp) :: radius_in_cm_for_create_initial_model
      real(dp) :: mass_in_gm_for_create_initial_model
      real(dp) :: initial_model_eps
      
      logical :: save_star_job_namelist
      character (len=256) :: star_job_namelist_name 
      
      logical :: load_saved_model
      character (len=256) :: saved_model_name 

      logical :: set_max_dt_to_frac_lifetime
      real(dp) :: max_frac_of_lifetime_per_step

      logical :: astero_just_call_my_extras_check_model

      logical :: relax_mass
      logical :: relax_initial_mass
      real(dp) :: new_mass
      real(dp) :: lg_max_abs_mdot
      logical :: relax_mass_scale
      real(dp) :: dlgm_per_step, relax_M_center_dt, change_mass_years_for_dt
      logical :: relax_M_center
      
      logical :: relax_core
      real(dp) :: new_core_mass, dlg_core_mass_per_step, &
         relax_core_years_for_dt, core_avg_rho, core_avg_eps
      
      logical :: relax_R_center
      real(dp) :: new_R_center
      real(dp) :: dlgR_per_step, relax_R_center_dt
      
      logical :: relax_L_center
      real(dp) :: new_L_center
      real(dp) :: dlgL_per_step, relax_L_center_dt
      
      integer :: remove_center_at_cell_k
      real(dp) :: remove_inner_fraction_q
      logical :: report_mass_not_fe56

      logical :: relax_dxdt_nuc_factor
      real(dp) :: new_dxdt_nuc_factor
      real(dp) :: dxdt_nuc_factor_multiplier

      logical :: relax_eps_nuc_factor
      real(dp) :: new_eps_nuc_factor
      real(dp) :: eps_nuc_factor_multiplier

      logical :: relax_opacity_max
      real(dp) :: new_opacity_max
      real(dp) :: opacity_max_multiplier

      logical :: relax_max_surf_dq
      real(dp) :: new_max_surf_dq
      real(dp) :: max_surf_dq_multiplier

      logical :: relax_surf_bc_offset_factor
      real(dp) :: new_bc_offset_factor
      real(dp) :: dlnbc_offset_factor

      logical :: relax_tau_factor, set_tau_factor
      real(dp) :: relax_to_this_tau_factor, set_to_this_tau_factor
      real(dp) :: dlogtau_factor, &
         set_tau_factor_after_core_He_burn, set_tau_factor_after_core_C_burn, &
         relax_tau_factor_after_core_He_burn, relax_tau_factor_after_core_C_burn

      logical :: relax_irradiation, set_irradiation
      integer :: relax_irradiation_min_steps
      real(dp) :: relax_to_this_irrad_flux, set_to_this_irrad_flux, &
         irrad_col_depth, relax_irradiation_max_yrs_dt

      logical :: relax_mass_change
      integer :: relax_mass_change_min_steps
      real(dp) :: relax_mass_change_max_yrs_dt, &
         relax_mass_change_init_mdot, relax_mass_change_final_mdot
      
      logical :: change_lnPgas_flag, new_lnPgas_flag
      
      logical :: change_v_flag, new_v_flag
      real(dp) :: center_ye_limit_for_v_flag
      real(dp) :: gamma1_integral_for_v_flag
      
      logical :: change_rotation_flag, new_rotation_flag
      
      logical :: set_omega, set_initial_omega
      integer :: set_omega_step_limit, set_near_zams_omega_steps
      real(dp) :: new_omega
      
      logical :: set_omega_div_omega_crit, set_initial_omega_div_omega_crit
      integer :: set_omega_div_omega_crit_step_limit, set_near_zams_omega_div_omega_crit_steps
      real(dp) :: new_omega_div_omega_crit
      
      logical :: set_surface_rotation_v, set_initial_surface_rotation_v
      integer :: set_surf_rotation_v_step_limit, set_near_zams_surface_rotation_v_steps
      real(dp) :: new_surface_rotation_v
      
      logical :: &
         relax_omega, &
         relax_initial_omega, &
         near_zams_relax_omega, &
         relax_omega_div_omega_crit, &
         relax_initial_omega_div_omega_crit, &
         near_zams_relax_omega_div_omega_crit, &
         relax_surface_rotation_v, &
         relax_initial_surface_rotation_v, &
         near_zams_relax_initial_surface_rotation_v
      integer :: num_steps_to_relax_rotation
      
      logical :: set_uniform_initial_composition
      real(dp) :: initial_h1, initial_h2, initial_he3, initial_he4
      integer :: initial_zfracs
      
      logical :: relax_initial_composition, relax_initial_to_xaccrete
      character (len=256) :: relax_composition_filename
      integer :: num_steps_to_relax_composition
      
      real(dp) :: report_cell_for_xm
      logical :: set_to_xa_for_accretion
      integer :: set_nzlo, set_nzhi
      
      logical :: set_kap_base_CNO_Z_fracs
      real(dp) :: kap_base_fC
      real(dp) :: kap_base_fN
      real(dp) :: kap_base_fO
            
      logical :: change_Y
      logical :: change_initial_Y
      logical :: relax_Y
      real(dp) :: new_Y
      
      logical :: change_Z
      logical :: change_initial_Z
      logical :: relax_Z
      real(dp) :: new_Z
      
      integer :: steps_to_take_before_terminate
      character (len=256) :: stop_if_this_file_exists

      logical :: set_initial_age
      real(dp) :: initial_age
      
      logical :: set_initial_model_number
      integer :: initial_model_number

      logical :: set_initial_dt, limit_initial_dt
      real(dp) :: years_for_initial_dt

      logical :: set_Z_all_HELM
      real(dp) :: Z_all_HELM

      logical :: set_logRho_OPAL_SCVH_limits
      real(dp) :: logRho1_OPAL_SCVH_limit, logRho2_OPAL_SCVH_limit

      logical :: set_HELM_OPAL_lgTs
      real(dp) :: logT_all_HELM, logT_all_OPAL

      logical :: set_HELM_SCVH_lgTs
      real(dp) :: logT_low_all_HELM, logT_low_all_SCVH

      logical :: set_eos_PC_parameters
      real(dp) :: mass_fraction_limit_for_PC, logRho1_PC_limit, logRho2_PC_limit
      real(dp) :: log_Gamma_all_HELM, log_Gamma_all_PC, PC_min_Z

      logical :: change_net
      character (len=256) :: new_net_name, h_he_net, co_net, adv_net
      
      logical :: set_rates_preference
      integer :: new_rates_preference
      character (len=16) :: set_rate_c12ag, set_rate_n14pg, set_rate_3a, &
         set_rate_1212_to_ne20, set_rate_1212

      logical :: set_uniform_xa_from_file, set_uniform_initial_xa_from_file
      character (len=256) :: file_for_uniform_xa

      logical :: turn_on_T_limits, turn_off_T_limits
      real(dp) :: max_factor_2n, &
         T_lo_2n, T_hi_2n, min_factor_2n, &
         T_lo_prot, T_hi_prot, min_factor_prot, & 
         T_lo_combo_a_capture, T_hi_combo_a_capture, min_factor_combo_a_capture, &
         T_lo_a_cap_high_mass, T_hi_a_cap_high_mass, min_factor_a_cap_high_mass, &  
         T_lo_a_cap_intermediate, T_hi_a_cap_intermediate, min_factor_a_cap_intermediate
      
      real(dp) :: mix_envelope_down_to_T, mix_initial_envelope_down_to_T
      
      logical :: auto_extend_net
      
      integer :: save_model_number
      character (len=256) :: save_model_filename
      logical :: save_prev_along_with_current, save_model_when_terminate
      
      logical :: profile_starting_model
      integer :: profile_model_number
      
      integer :: internals_num
      
      logical :: report_retries
      logical :: report_backups
      
      logical :: do_abundance_smoothing
      integer :: abundance_smoothing_cnt, abundance_smoothing_nzlo, abundance_smoothing_nzhi
      
      character (len=256) :: net_reaction_filename
      
      character (len=256) :: rate_tables_dir, rate_cache_suffix
      
      logical :: read_extra_star_job_inlist1
      character (len=256) :: extra_star_job_inlist1_name 
      
      logical :: read_extra_star_job_inlist2
      character (len=256) :: extra_star_job_inlist2_name 
      
      logical :: read_extra_star_job_inlist3
      character (len=256) :: extra_star_job_inlist3_name 
      
      logical :: read_extra_star_job_inlist4
      character (len=256) :: extra_star_job_inlist4_name 
      
      logical :: read_extra_star_job_inlist5
      character (len=256) :: extra_star_job_inlist5_name 

      integer :: nzlo, nzhi
      logical :: set_abundance
      integer :: chem_id
      character (len=5) :: chem_name
      real(dp) :: new_frac
      
      logical :: replace_element
      integer :: chem_id1, chem_id2
      character (len=5) :: chem_name1, chem_name2
      
      logical :: do_special_test
      
      integer :: save_pulsation_info_for_model_number
      logical :: save_pulsation_info_when_terminate
      character (len=256) :: save_pulsation_info_filename
      
      integer :: save_short_format_for_model_number
      character (len=256) :: save_short_format_filename 
      
      character (len=256) :: chem_isotopes_filename
      
      character(len=256) :: kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix, &
         other_kappa_file_prefix
      real(dp) :: kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy, &
         kappa_type2_logT_lower_bdy
      character(len=256) :: eos_file_prefix, other_eos_file_prefix
      character(len=256) :: ionization_file_prefix, ionization_Z1_suffix
      character(len=256) :: eosDT_Z1_suffix, eosPT_Z1_suffix
      
      integer, parameter :: max_extras_params = 20, max_extras_cpar_len = 256
      integer :: extras_lipar, extras_lrpar, extras_lcpar, extras_llpar
      integer :: extras_ipar(max_extras_params)
      real(dp) :: extras_rpar(max_extras_params)
      character(len=max_extras_cpar_len) :: extras_cpar(max_extras_params)
      logical :: extras_lpar(max_extras_params)
      
      integer, parameter :: max_num_special_rate_factors = 20
      integer :: num_special_rate_factors
      real(dp) :: special_rate_factor(max_num_special_rate_factors)
      character(len=maxlen_reaction_Name) :: &
         reaction_for_special_factor(max_num_special_rate_factors)


      namelist /star_job/ &
         pause_before_terminate, &
         mesa_dir, eosDT_cache_dir, eosPT_cache_dir, &
         ionization_cache_dir, kap_cache_dir, rates_cache_dir, weaklib_cache_dir, &
         profile_columns_file, history_columns_file, save_photo_when_terminate, &
         create_pre_main_sequence_model, pre_ms_relax_num_steps, load_saved_model, saved_model_name, &
         create_initial_model, initial_model_relax_num_steps, initial_model_eps, &
         radius_in_cm_for_create_initial_model, mass_in_gm_for_create_initial_model, &
         steps_to_take_before_terminate, stop_if_this_file_exists, &
         set_initial_age, initial_age, set_initial_dt, limit_initial_dt, years_for_initial_dt, &
         save_star_job_namelist, star_job_namelist_name, &
         set_logRho_OPAL_SCVH_limits, logRho1_OPAL_SCVH_limit, logRho2_OPAL_SCVH_limit, &
         set_HELM_SCVH_lgTs, logT_low_all_HELM, logT_low_all_SCVH, &
         set_Z_all_HELM, Z_all_HELM, &
         set_HELM_OPAL_lgTs, logT_all_HELM, logT_all_OPAL, &
         set_eos_PC_parameters, mass_fraction_limit_for_PC, &
         logRho1_PC_limit, logRho2_PC_limit, log_Gamma_all_HELM, log_Gamma_all_PC, PC_min_Z, &
         change_net, new_net_name, auto_extend_net, h_he_net, co_net, adv_net, &
         set_rates_preference, new_rates_preference, &
         turn_on_T_limits, turn_off_T_limits, max_factor_2n, &
         T_lo_2n, T_hi_2n, min_factor_2n, &
         T_lo_prot, T_hi_prot, min_factor_prot, & 
         T_lo_combo_a_capture, T_hi_combo_a_capture, min_factor_combo_a_capture, &
         T_lo_a_cap_high_mass, T_hi_a_cap_high_mass, min_factor_a_cap_high_mass, &  
         T_lo_a_cap_intermediate, T_hi_a_cap_intermediate, min_factor_a_cap_intermediate, &     
         set_rate_c12ag, set_rate_n14pg, set_rate_3a, set_rate_1212, set_rate_1212_to_ne20, &
         num_special_rate_factors, reaction_for_special_factor, special_rate_factor, &
         save_pulsation_info_for_model_number, save_pulsation_info_when_terminate, &
         save_pulsation_info_filename, &
         save_short_format_for_model_number, save_short_format_filename, &
         set_uniform_xa_from_file, set_uniform_initial_xa_from_file, file_for_uniform_xa, &
         mix_envelope_down_to_T, mix_initial_envelope_down_to_T, &
         relax_mass, relax_initial_mass, new_mass, lg_max_abs_mdot, &
         relax_mass_scale, dlgm_per_step, change_mass_years_for_dt, relax_M_center, relax_M_center_dt, &
         relax_core, new_core_mass, dlg_core_mass_per_step, &
         relax_core_years_for_dt, core_avg_rho, core_avg_eps, &
         relax_L_center, new_L_center, dlgL_per_step, relax_L_center_dt, &
         relax_R_center, new_R_center, dlgR_per_step, relax_R_center_dt, &
         remove_center_at_cell_k, remove_inner_fraction_q, report_mass_not_fe56, &
         relax_dxdt_nuc_factor, new_dxdt_nuc_factor, dxdt_nuc_factor_multiplier, &
         relax_eps_nuc_factor, new_eps_nuc_factor, eps_nuc_factor_multiplier, &
         relax_opacity_max, new_opacity_max, opacity_max_multiplier, &
         relax_max_surf_dq, new_max_surf_dq, max_surf_dq_multiplier, &
         relax_surf_bc_offset_factor, new_bc_offset_factor, dlnbc_offset_factor, &
         relax_tau_factor, relax_to_this_tau_factor, dlogtau_factor, &
         relax_tau_factor_after_core_He_burn, relax_tau_factor_after_core_C_burn, &
         set_tau_factor, set_tau_factor_after_core_He_burn, &
         set_tau_factor_after_core_C_burn, set_to_this_tau_factor, &
         set_irradiation, relax_irradiation, relax_irradiation_min_steps, relax_irradiation_max_yrs_dt, &
         relax_to_this_irrad_flux, irrad_col_depth, set_to_this_irrad_flux, &
         relax_mass_change, relax_mass_change_min_steps, relax_mass_change_max_yrs_dt, &
         relax_mass_change_init_mdot, relax_mass_change_final_mdot, &
         change_Y, change_initial_Y, relax_Y, new_Y, relax_Z, change_Z, change_initial_Z, new_Z, &
         change_lnPgas_flag, new_lnPgas_flag, &
         change_v_flag, new_v_flag, center_ye_limit_for_v_flag, gamma1_integral_for_v_flag, &
         change_rotation_flag, new_rotation_flag, &
         
         set_omega, set_initial_omega, set_near_zams_omega_steps, set_omega_step_limit, new_omega, &
            
         set_omega_div_omega_crit, set_initial_omega_div_omega_crit, set_near_zams_omega_div_omega_crit_steps, &
            set_omega_div_omega_crit_step_limit, new_omega_div_omega_crit, &
         
         set_surface_rotation_v, set_initial_surface_rotation_v, set_near_zams_surface_rotation_v_steps, &
            set_surf_rotation_v_step_limit, new_surface_rotation_v, &
      
         relax_omega, &
         relax_initial_omega, &
         near_zams_relax_omega, &
         relax_omega_div_omega_crit, &
         relax_initial_omega_div_omega_crit, &
         near_zams_relax_omega_div_omega_crit, &
         relax_surface_rotation_v, &
         relax_initial_surface_rotation_v, &
         near_zams_relax_initial_surface_rotation_v, &
         num_steps_to_relax_rotation, &
            
         set_uniform_initial_composition, set_kap_base_CNO_Z_fracs, kap_base_fC, kap_base_fN, kap_base_fO, &
         relax_initial_composition, relax_initial_to_xaccrete, &
         num_steps_to_relax_composition, relax_composition_filename, &
         initial_h1, initial_h2, initial_he3, initial_he4, initial_zfracs, &
         save_model_when_terminate, save_model_number, &
         set_to_xa_for_accretion, set_nzlo, set_nzhi, report_cell_for_xm, &
         pre_ms_T_c, pre_ms_guess_rho_c, pre_ms_d_log10_P, pre_ms_logT_surf_limit, pre_ms_logP_surf_limit, &
         save_model_filename, save_prev_along_with_current, &
         profile_starting_model, profile_model_number, internals_num, report_retries, report_backups, &
         do_abundance_smoothing, abundance_smoothing_cnt, &
         abundance_smoothing_nzlo, abundance_smoothing_nzhi, &
         read_extra_star_job_inlist1, extra_star_job_inlist1_name, &
         read_extra_star_job_inlist2, extra_star_job_inlist2_name, &
         read_extra_star_job_inlist3, extra_star_job_inlist3_name, &
         read_extra_star_job_inlist4, extra_star_job_inlist4_name, &
         read_extra_star_job_inlist5, extra_star_job_inlist5_name, &
         nzlo, nzhi, set_abundance, chem_name, new_frac, &
         replace_element, chem_name1, chem_name2, &
         net_reaction_filename, rate_tables_dir, rate_cache_suffix, first_model_for_timing, &
         show_log_description_at_start, show_net_reactions_info, pgstar_flag, use_se_output, &
         use_se_names, se_codev, se_modname, se_prefix, se_num_mod_output, &
         show_net_species_info, show_eqns_and_vars_names, &
         set_max_dt_to_frac_lifetime, max_frac_of_lifetime_per_step, &
         chem_isotopes_filename, astero_just_call_my_extras_check_model, &
         do_special_test, kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix, &
         kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy, &
         kappa_type2_logT_lower_bdy, eos_file_prefix, &
         other_eos_file_prefix, other_kappa_file_prefix, &
         ionization_file_prefix, ionization_Z1_suffix, &
         eosDT_Z1_suffix, eosPT_Z1_suffix, &
         set_initial_model_number, initial_model_number, &
         extras_lipar, extras_lrpar, extras_lcpar, extras_llpar, &
         extras_ipar, extras_rpar, extras_cpar, extras_lpar
      
      real(dp) :: &
         step_loop_timing, after_step_timing, before_step_timing, &
         check_time_start, check_time_end, &
         check_step_loop_timing, check_after_step_timing, check_before_step_timing
            
            
      contains
      
      
      subroutine do_load1_star(id, s, restart, restart_filename, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         character (len=*), intent(in) :: restart_filename
         integer, intent(out) :: ierr
      
         if (restart) then
            call star_load_restart_photo(id, restart_filename, ierr)
            if (failed('star_load_restart_photo')) return
         else if (load_saved_model) then
            if (create_pre_main_sequence_model) then
               write(*,*) 'you have both load_saved_model and create_pre_main_sequence_model set true'
               write(*,*) 'please pick one and try again'
               stop 1
            end if
            if (create_initial_model) then
               write(*,*) 'you have both load_saved_model and create_initial_model set true'
               write(*,*) 'please pick one and try again'
               stop 1
            end if
            write(*,'(a)') 'load saved model ' // trim(saved_model_name)
            write(*,*)
            call star_read_model(id, saved_model_name, ierr)
            if (failed('star_read_model')) return
         else if (create_pre_main_sequence_model) then
            if (.not. restart) write(*, *) 'create pre-main-sequence model'
            if (create_initial_model) then
               write(*,*) 'you have both create_pre_main_sequence_model and create_initial_model set true'
               write(*,*) 'please pick one and try again'
               stop 1
            end if
            call star_create_pre_ms_model(id, pre_ms_T_c, pre_ms_guess_rho_c, &
                  pre_ms_d_log10_P, pre_ms_logT_surf_limit, pre_ms_logP_surf_limit, &
                  initial_zfracs, pre_ms_relax_num_steps, ierr)
            if (failed('star_create_pre_ms_model')) return
         else if (create_initial_model) then
            if (.not. restart) write(*, *) 'create initial model'
            if (create_pre_main_sequence_model) then
               write(*,*) 'you have both create_initial_model and create_pre_main_sequence_model set true'
               write(*,*) 'please pick one and try again'
               stop 1
            end if
            call star_create_initial_model(id, &
                  radius_in_cm_for_create_initial_model, mass_in_gm_for_create_initial_model, &
                  initial_zfracs, initial_model_relax_num_steps, initial_model_eps, ierr)
            if (failed('star_create_initial_model')) return
         else ! use mesa starting models
            call star_load_zams(id, ierr)
            if (failed('star_load_zams')) return
         end if
         
         contains
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed

      end subroutine do_load1_star
      
      
      subroutine set_which_rates(id, ierr)
         use rates_def
         use rates_lib
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: which_rate
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (set_rates_preference) then
            write(*,*) 'change rates preference to', new_rates_preference
            s% which_rates(:) = new_rates_preference
         else
            s% which_rates(:) = rates_NACRE_if_available
         end if
         
         if (len_trim(set_rate_c12ag) > 0) then
            if (set_rate_c12ag == 'NACRE') then
               which_rate = use_rate_c12ag_NACRE
            else if (set_rate_c12ag == 'Buchmann') then
               which_rate = use_rate_c12ag_JR
            else if (set_rate_c12ag == 'Kunz') then
               which_rate = use_rate_c12ag_Kunz
            else if (set_rate_c12ag == 'CF88') then
               which_rate = use_rate_c12ag_CF88
            else
               write(*,*) 'invalid string for set_rate_c12ag ' // trim(set_rate_c12ag)
               ierr = -1
               return
            end if
            call set_which_rate_c12ag(s% which_rates, which_rate)
         end if
         
         if (len_trim(set_rate_n14pg) > 0) then
            if (set_rate_n14pg == 'NACRE') then
               which_rate = use_rate_n14pg_NACRE
            else if (set_rate_n14pg == 'Imbriani') then
               which_rate = use_rate_n14pg_JR
            else if (set_rate_n14pg == 'CF88') then
               which_rate = use_rate_n14pg_CF88
            else
               write(*,*) 'invalid string for set_rate_n14pg ' // trim(set_rate_n14pg)
               ierr = -1
               return
            end if
            call set_which_rate_n14pg(s% which_rates, which_rate)
         end if
         
         if (len_trim(set_rate_3a) > 0) then
            if (set_rate_3a == 'NACRE') then
               which_rate = use_rate_3a_NACRE
            else if (set_rate_3a == 'Fynbo') then
               which_rate = use_rate_3a_JR
            else if (set_rate_3a == 'CF88') then
               which_rate = use_rate_3a_CF88
            else if (set_rate_3a == 'FL87') then
               which_rate = use_rate_3a_FL87
            else
               write(*,*) 'invalid string for set_rate_3a ' // trim(set_rate_3a)
               ierr = -1
               return
            end if
            call set_which_rate_3a(s% which_rates, which_rate)
         end if
         
         if (len_trim(set_rate_1212) > 0) then
            if (set_rate_1212 == 'CF88') then
               which_rate = use_rate_1212_CF88
            else if (set_rate_1212 == 'G05') then
               which_rate = use_rate_1212_G05
            else
               write(*,*) 'invalid string for set_rate_1212 ' // trim(set_rate_1212)
               ierr = -1
               return
            end if
            call set_which_rate_1212(s% which_rates, which_rate)
         end if
         
         if (len_trim(set_rate_1212_to_ne20) > 0) then
            if (set_rate_1212_to_ne20 == 'CF88') then
               which_rate = use_rate_1212_CF88
            else if (set_rate_1212_to_ne20 == 'G05') then
               which_rate = use_rate_1212_G05
            else
               write(*,*) 'invalid string for set_rate_1212_to_ne20 ' // trim(set_rate_1212_to_ne20)
               ierr = -1
               return
            end if
            call set_which_rate_1212(s% which_rates, which_rate)
         end if

      end subroutine set_which_rates
      
      
      subroutine set_rate_factors(id, ierr)
         use net_lib, only: get_net_reaction_table_ptr
         use rates_lib, only: rates_reaction_id
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: j, i, ir
         integer, pointer :: net_reaction_ptr(:) 
         
         include 'formats.inc'
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% rate_factors(:) = 1
         if (num_special_rate_factors <= 0) return
         
         call get_net_reaction_table_ptr(s% net_handle, net_reaction_ptr, ierr)
         if (ierr /= 0) return
         
         do i=1,num_special_rate_factors
            if (len_trim(reaction_for_special_factor(i)) == 0) cycle
            ir = rates_reaction_id(reaction_for_special_factor(i))
            j = 0
            if (ir > 0) j = net_reaction_ptr(ir)
            if (j <= 0) cycle
            s% rate_factors(j) = special_rate_factor(i)
            write(*,2) 'set special rate factor for ' // &
                  trim(reaction_for_special_factor(i)), j, special_rate_factor(i)
         end do
         
      end subroutine set_rate_factors
      

      subroutine do_star_job_controls_before(id, s, restart, ierr)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr

         integer :: which_atm

         include 'formats.inc'
      
         ierr = 0
         
         s% set_which_rates => set_which_rates ! will be called after net is defined
         s% set_rate_factors => set_rate_factors ! will be called after net is defined
         
         if (set_tau_factor) then
            write(*,1) 'set_tau_factor', set_to_this_tau_factor
            s% tau_factor = set_to_this_tau_factor
         end if

         if (set_Z_all_HELM) then
            write(*,*) 'set_Z_all_HELM'
            write(*,1) 'Z_all_HELM', Z_all_HELM
            call eos_set_Z_all_HELM(s% eos_handle, Z_all_HELM, ierr)
            if (failed('eos_set_Z_all_HELM')) return
         end if

         if (set_logRho_OPAL_SCVH_limits) then
            write(*,*) 'set_logRho_OPAL_SCVH_limits'
            write(*,1) 'logRho1_OPAL_SCVH_limit', logRho1_OPAL_SCVH_limit
            write(*,1) 'logRho2_OPAL_SCVH_limit', logRho2_OPAL_SCVH_limit
            call eos_set_logRhos_OPAL_SCVH( &
                  s% eos_handle, logRho1_OPAL_SCVH_limit, logRho2_OPAL_SCVH_limit, ierr)
            if (failed('eos_set_logRhos_OPAL_SCVH')) return
         end if

         if (set_HELM_SCVH_lgTs) then
            write(*,*) 'set_HELM_SCVH_lgTs'
            write(*,1) 'logT_low_all_HELM', logT_low_all_HELM
            write(*,1) 'logT_low_all_SCVH', logT_low_all_SCVH
            if (logT_low_all_HELM < 2.1d0) then
               write(*,*) 'for current eos tables, min allowed logT_low_all_HELM is 2.1'
               ierr = -1
               return
            end if
            if (logT_low_all_HELM > logT_low_all_SCVH) then
               write(*,*) 'logT_low_all_HELM must be <= logT_low_all_SCVH'
               ierr = -1
               return
            end if
            call eos_set_HELM_SCVH_lgTs(s% eos_handle, logT_low_all_HELM, logT_low_all_SCVH, ierr)
            if (failed('eos_set_HELM_SCVH_lgTs')) return
         end if

         if (set_HELM_OPAL_lgTs) then
            write(*,*) 'set_HELM_OPAL_lgTs'
            write(*,1) 'logT_all_HELM', logT_all_HELM
            write(*,1) 'logT_all_OPAL', logT_all_OPAL
            if (logT_all_HELM > 7.7d0) then
               write(*,*) 'for current eos tables, max allowed logT_all_HELM is 7.7'
               ierr = -1
               return
            end if
            if (logT_all_HELM < logT_all_OPAL) then
               write(*,*) 'logT_all_HELM must be >= logT_all_OPAL'
               ierr = -1
               return
            end if
            call eos_set_HELM_OPAL_lgTs(s% eos_handle, logT_all_HELM, logT_all_OPAL, ierr)
            if (failed('eos_set_HELM_OPAL_lgTs')) return
         end if

         if (set_eos_PC_parameters) then
            write(*,*) 'set_eos_PC_parameters'
            write(*,1) 'mass_fraction_limit_for_PC', mass_fraction_limit_for_PC
            write(*,1) 'logRho1_PC_limit', logRho1_PC_limit
            write(*,1) 'logRho2_PC_limit', logRho2_PC_limit
            write(*,1) 'log_Gamma_all_HELM', log_Gamma_all_HELM
            write(*,1) 'log_Gamma_all_PC', log_Gamma_all_PC
            write(*,1) 'PC_min_Z', PC_min_Z
            call eos_set_PC_parameters(s% eos_handle, &
               mass_fraction_limit_for_PC, logRho1_PC_limit, logRho2_PC_limit, &
               log_Gamma_all_HELM, log_Gamma_all_PC, PC_min_Z, ierr)
            if (failed('star_set_eos_PC_params')) return
         end if
         
         which_atm = atm_option(s% which_atm_option, ierr)
         if (failed('atm_option')) return
         s% tau_base = atm_tau_base(which_atm, ierr)
         if (failed('atm_tau_base')) return
         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed

      end subroutine do_star_job_controls_before
      

      subroutine do_star_job_controls_after(id, s, restart, ierr)
         use const_def
         use omp_lib
         use net_lib, only: net_turn_on_T_limits, net_turn_off_T_limits
         use reaclib_def, only: pretty_print_format
         use rates_def
         use rates_lib

         integer, intent(in) :: id
         type (star_info), pointer :: s
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         
         real(dp) :: log_m, log_lifetime, max_dt, minq, maxq
         integer :: i, j, k

         include 'formats.inc'
         
         if (set_initial_age .and. .not. restart) then
            write(*,*) 'set_initial_age', initial_age ! in years
            call star_set_age(id, initial_age, ierr)
            if (failed('star_set_age')) return
         end if

         if (set_initial_dt .and. .not. restart) then
            write(*,*) 'set_initial_dt', years_for_initial_dt
            !s% dt_next = min(s% dt_next, years_for_initial_dt*secyer)
            s% dt_next = years_for_initial_dt*secyer
         end if

         if (limit_initial_dt .and. .not. restart) then
            write(*,*) 'limit_initial_dt', years_for_initial_dt
            s% dt_next = min(s% dt_next, years_for_initial_dt*secyer)
         end if

         if (set_initial_model_number .and. .not. restart) then
            write(*,*) 'set_initial_model_number', initial_model_number
            s% model_number = initial_model_number
         end if

         if (steps_to_take_before_terminate > 0) then
            s% max_model_number = s% model_number + steps_to_take_before_terminate
            write(*,2) 'steps_to_take_before_terminate', steps_to_take_before_terminate
            write(*,2) 'max_model_number', s% max_model_number
         end if

         if (change_net) then
            write(*,*) 'change to ' // trim(new_net_name)
            call star_change_to_new_net(id, new_net_name, ierr)
            if (failed('star_change_to_new_net')) return
            write(*,*) 'number of species =', s% species
         end if
         
         if (turn_on_T_limits) then
            call set_net_T_limits
            if (failed('turn_on_T_limits')) return
         end if
         
         if (turn_off_T_limits) then
            call net_turn_off_T_limits(s% net_handle, ierr)
            if (failed('net_turn_off_T_limits')) return
         end if

         if (change_lnPgas_flag) then
            write(*,*) 'new_lnPgas_flag =', new_lnPgas_flag
            call star_set_lnPgas_flag(id, new_lnPgas_flag, ierr)
            if (failed('star_set_lnPgas_flag')) return
         end if

         if (change_v_flag) then
            write(*,*) 'new_v_flag =', new_v_flag
            call star_set_v_flag(id, new_v_flag, ierr)
            if (failed('star_set_v_flag')) return
         end if

         if (change_rotation_flag) then
            write(*,*) 'new_rotation_flag =', new_rotation_flag
            call star_set_rotation_flag(id, new_rotation_flag, ierr)
            if (failed('star_set_rotation_flag')) return
         end if

         if (s% rotation_flag .and. set_omega) then
            write(*,*) 'new_omega =', new_omega
            call star_set_uniform_omega(id, new_omega, ierr)
            if (failed('star_set_uniform_omega')) return
         end if

         if (s% rotation_flag .and. set_initial_omega .and. .not. restart) then
            write(*,*) 'new_omega =', new_omega
            call star_set_uniform_omega(id, new_omega, ierr)
            if (failed('star_set_uniform_omega')) return
         end if

         if (s% rotation_flag .and. set_surface_rotation_v) then
            new_omega = new_surface_rotation_v*1d5/s% r(1)
            write(*,1) 'new_surface_rotation_v =', new_surface_rotation_v, new_omega
            call star_set_uniform_omega(id, new_omega, ierr)
            if (failed('star_set_uniform_omega')) return
         end if

         if (s% rotation_flag .and. set_initial_surface_rotation_v .and. .not. restart) then
            new_omega = new_surface_rotation_v*1d5/s% r(1)
            write(*,2) 'new_surface_rotation_v', &
               s% model_number, new_surface_rotation_v, new_omega
            call star_set_uniform_omega(id, new_omega, ierr)
            if (failed('star_set_uniform_omega')) return
         end if

         if (s% rotation_flag .and. set_omega_div_omega_crit) then
            new_omega = new_omega_div_omega_crit*star_surface_omega_crit(id, ierr)
            if (failed('star_surface_omega_crit')) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, new_omega_div_omega_crit, new_omega
            call star_set_uniform_omega(id, new_omega, ierr)
            if (failed('star_set_uniform_omega')) return
         end if

         if (s% rotation_flag .and. set_initial_omega_div_omega_crit .and. .not. restart) then
            new_omega = new_omega_div_omega_crit*star_surface_omega_crit(id, ierr)
            if (failed('star_surface_omega_crit')) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, new_omega_div_omega_crit, new_omega
            call star_set_uniform_omega(id, new_omega, ierr)
            if (failed('star_set_uniform_omega')) return
         end if
         
         if (set_to_xa_for_accretion) then
            write(*,*) 'set_to_xa_for_accretion'
            call change_to_xa_for_accretion(id, set_nzlo, set_nzhi, ierr)
            if (failed('set_to_xa_for_accretion')) return
         end if
         
         if (first_model_for_timing > 0) write(*,2) 'first_model_for_timing', first_model_for_timing
         
         if (set_uniform_initial_composition .and. .not. restart) then
            write(*,*)
            write(*,1) 'set_uniform_initial_composition'
            write(*,1) 'initial_h1', initial_h1
            write(*,1) 'initial_h2', initial_h2
            write(*,1) 'initial_he3', initial_he3
            write(*,1) 'initial_he4', initial_he4
            select case(initial_zfracs)
               case (AG89_zfracs)
                  write(*,1) 'metals AG89'
               case (GN93_zfracs)
                  write(*,1) 'metals GN93'
               case (GS98_zfracs)
                  write(*,1) 'metals GS98'
               case (L03_zfracs)
                  write(*,1) 'metals L03'
               case (AGS04_zfracs)
                  write(*,1) 'metals AGS04'
               case (AGSS09_zfracs)
                  write(*,1) 'metals AGSS09'
               case (L09_zfracs)
                  write(*,1) 'metals L09'
               case default
                  write(*,2) 'unknown value for initial_zfracs', initial_zfracs
            end select
            call star_set_standard_composition( &
               id, initial_h1, initial_h2, initial_he3, initial_he4, initial_zfracs, ierr)
            if (failed('set_uniform_initial_composition')) return
         end if
         
         if (relax_initial_composition .and. .not. restart) then
            call do_relax_initial_composition(ierr)
            if (failed('do_relax_initial_composition')) return
         end if
         
         if (relax_initial_to_xaccrete .and. .not. restart) then
            call star_relax_to_xaccrete(id, num_steps_to_relax_composition, ierr)
            if (failed('star_relax_to_xaccrete')) return
         end if
         
         if (set_kap_base_CNO_Z_fracs) then
            write(*,1) 'set_kap_base_CNO_Z_fracs - no longer supported'
            ierr = -1
            if (failed('kap_set_base_CNO_Z_fractions')) return
         end if

         if (set_uniform_xa_from_file) then
            call star_uniform_xa_from_file(id, file_for_uniform_xa, ierr)
            if (failed('star_uniform_xa_from_file')) return
         end if

         if (mix_initial_envelope_down_to_T > 0d0 .and. .not. restart) then
            call uniform_mix_envelope_down_to_T(id, mix_initial_envelope_down_to_T, ierr)
            if (failed('uniform_mix_envelope_down_to_T')) return
         end if

         if (mix_envelope_down_to_T > 0d0) then
            call uniform_mix_envelope_down_to_T(id, mix_envelope_down_to_T, ierr)
            if (failed('uniform_mix_envelope_down_to_T')) return
         end if

         if (mix_initial_envelope_down_to_T > 0d0) then
            call uniform_mix_envelope_down_to_T(id, mix_initial_envelope_down_to_T, ierr)
            if (failed('uniform_mix_envelope_down_to_T')) return
         end if

         if (set_uniform_initial_xa_from_file .and. .not. restart) then
            call star_uniform_xa_from_file(id, file_for_uniform_xa, ierr)
            if (failed('star_uniform_xa_from_file')) return
         end if

         if (change_Y) then
            write(*,1) 'new_Y', new_Y
            call star_set_y(id, new_Y, ierr)
            if (failed('change_Y')) return
            write(*, 1) 'new y', get_current_abundance(id, ihe4, ierr)
            if (failed('get_current_abundance')) return
         end if

         if (change_initial_Y .and. .not. restart) then
            write(*,1) 'new_Y', new_Y
            call star_set_y(id, new_Y, ierr)
            if (failed('change_initial_Y')) return
            write(*, 1) 'new y', get_current_abundance(id, ihe4, ierr)
            if (failed('get_current_abundance')) return
         end if

         if (change_Z) then
            write(*,1) 'new_Z', new_Z
            call star_set_z(id, new_Z, ierr)
            if (failed('star_set_z')) return
            write(*, 1) 'new z', get_current_z(id, ierr)
            if (failed('get_current_z')) return
         end if

         if (change_initial_Z .and. .not. restart) then
            write(*,1) 'new_Z', new_Z
            call star_set_z(id, new_Z, ierr)
            if (failed('star_set_z')) return
            write(*, 1) 'new z', get_current_z(id, ierr)
            if (failed('get_current_z')) return
         end if
         
         if (set_abundance) then
            nzhi = set_nzhi
            nzlo = set_nzlo
            if (nzhi == 0) nzhi = s% nz
            write(*, *) 'set_abundance of ', trim(chem_name), new_frac, nzlo, nzhi 
            chem_id = get_nuclide_index(chem_name)
            if (chem_id <= 0) then
               write(*,*) 'failed to find ' // trim(chem_name)
               write(*,*) 'check valid chem_isos% names in chem/public/chem_def.f'
            end if
            call set_abundance_in_section(id, chem_id, new_frac, nzlo, nzhi, ierr)
            if (failed('set_abundance_in_section')) return
         end if
         
         if (do_abundance_smoothing) then
            if (abundance_smoothing_nzlo == 0) abundance_smoothing_nzlo = s% nz
            write(*, *) 'do_abundance_smoothing ', &
               abundance_smoothing_cnt, abundance_smoothing_nzlo, abundance_smoothing_nzhi 
            call smooth_abundances_in_section(id, &
               abundance_smoothing_cnt, abundance_smoothing_nzlo, abundance_smoothing_nzhi, ierr)
            if (failed('smooth_abundances_in_section')) return
         end if
         
         if (replace_element) then
            write(*, *) 'replace_element ', trim(chem_name1), ' by ', trim(chem_name2)
            chem_id1 = get_nuclide_index(chem_name1)
            chem_id2 = get_nuclide_index(chem_name2)
            if (chem_id1 <= 0) then
               write(*,*) 'failed to find ' // trim(chem_name1)
               write(*,*) 'check valid chem_isos% names in chem/public/chem_def.f'
            end if
            if (chem_id2 <= 0) then
               write(*,*) 'failed to find ' // trim(chem_name2)
               write(*,*) 'check valid chem_isos% names in chem/public/chem_def.f'
            end if
            if (nzhi <= 0) nzhi = s% nz
            if (nzlo <= 0) nzlo = 1
            write(*, *) 'in section', nzlo, nzhi
            call replace_element_in_section(id, chem_id1, chem_id2, nzlo, nzhi, ierr)
            if (failed('replace_element_in_section')) return
         end if

         if (set_irradiation) then
            write(*,2) 'set_irradiation'
            s% irradiation_flux = set_to_this_irrad_flux
            s% column_depth_for_irradiation = irrad_col_depth
         end if
         
         if (do_special_test) then
            write(*, *) 'do_special_test'
            call star_special_test(id, ierr)
            if (failed('star_special_test')) return
         end if
         
         if (remove_inner_fraction_q > 0 .and. remove_inner_fraction_q < 1) then
            write(*, 1) 'remove_inner_fraction_q', remove_inner_fraction_q
            remove_center_at_cell_k = s% nz
            do k = 1, s% nz
               if (s% q(k) <= remove_inner_fraction_q) then
                  remove_center_at_cell_k = k
                  exit
               end if
            end do
         end if
         
         if (remove_center_at_cell_k > 0 .and. remove_center_at_cell_k <= s% nz) then
            write(*, 2) 'remove_center_at_cell_k', remove_center_at_cell_k
            call star_remove_center(id, remove_center_at_cell_k, ierr)
            if (failed('star_remove_center')) return
         end if


         ! do "set" before "relax"


         if (relax_Y) then
            write(*,1) 'relax_Y', new_Y
            call star_relax_Y(id, new_Y, s% relax_dY, ierr)
            if (failed('star_relax_Y')) return
            write(*, 1) 'new y', get_current_y(id, ierr)
            if (failed('get_current_y')) return
         end if

         if (relax_Z) then
            minq = 0
            maxq = 1
            call star_relax_Z(id, new_Z, s% relax_dlnZ, minq, maxq, ierr)
            if (failed('star_relax_Z')) return
            write(*, 1) 'new z', get_current_z(id, ierr)
            if (failed('get_current_z')) return
         end if

         if (relax_mass) then
            write(*, 1) 'relax_mass', new_mass
            call star_relax_mass(id, new_mass, lg_max_abs_mdot, ierr)
            if (failed('star_relax_mass')) return
         end if

         if (relax_dxdt_nuc_factor) then
            write(*, 1) 'relax_dxdt_nuc_factor', new_dxdt_nuc_factor
            call star_relax_dxdt_nuc_factor(id, new_dxdt_nuc_factor, dxdt_nuc_factor_multiplier, ierr)
            if (failed('star_relax_dxdt_nuc_factor')) return
         end if

         if (relax_eps_nuc_factor) then
            write(*, 1) 'relax_eps_nuc_factor', new_eps_nuc_factor
            call star_relax_eps_nuc_factor(id, new_eps_nuc_factor, eps_nuc_factor_multiplier, ierr)
            if (failed('star_relax_eps_nuc_factor')) return
         end if

         if (relax_surf_bc_offset_factor) then
            write(*, 1) 'relax_surf_bc_offset_factor', new_bc_offset_factor
            call star_relax_bc_offset(id, new_bc_offset_factor, dlnbc_offset_factor, ierr)
            if (failed('star_relax_bc_offset')) return
         end if

         if (relax_opacity_max) then
            write(*, 1) 'relax_opacity_max', new_opacity_max
            call star_relax_opacity_max(id, new_opacity_max, opacity_max_multiplier, ierr)
            if (failed('star_relax_opacity_max')) return
         end if

         if (relax_max_surf_dq) then
            write(*, 1) 'relax_max_surf_dq', new_max_surf_dq
            call star_relax_max_surf_dq(id, new_max_surf_dq, max_surf_dq_multiplier, ierr)
            if (failed('star_relax_max_surf_dq')) return
         end if

         if (relax_initial_mass .and. .not. restart) then
            write(*, 1) 'relax_initial_mass to new_mass =', new_mass
            call star_relax_mass(id, new_mass, lg_max_abs_mdot, ierr)
            if (failed('relax_initial_mass')) return
         end if

         if (relax_mass_scale) then
            write(*, 1) 'relax_mass_scale', new_mass
            call star_relax_mass_scale( &
               id, new_mass, dlgm_per_step, change_mass_years_for_dt, ierr)
            if (failed('star_relax_mass_scale')) return
         end if

         if (relax_core) then
            write(*, 1) 'relax_core', new_core_mass
            call star_relax_core( &
               id, new_core_mass, dlg_core_mass_per_step, &
               relax_core_years_for_dt, core_avg_rho, core_avg_eps, ierr)
            if (failed('star_relax_core')) return
         end if
         
         if (relax_M_center) then
            write(*, 1) 'relax_M_center', new_mass
            call star_relax_M_center(id, new_mass, dlgm_per_step, relax_M_center_dt, ierr)
            if (failed('star_relax_M_center')) return
         end if
         
         if (relax_R_center) then
            write(*, 1) 'relax_R_center', new_R_center
            call star_relax_R_center(id, new_R_center, dlgR_per_step, relax_R_center_dt, ierr)
            if (failed('star_relax_R_center')) return
         end if
         
         if (relax_L_center) then
            write(*, 1) 'relax_L_center', new_L_center
            call star_relax_L_center(id, new_L_center, dlgL_per_step, relax_L_center_dt, ierr)
            if (failed('star_relax_L_center')) return
         end if
         
         if (relax_tau_factor) then
            write(*,1) 'relax_tau_factor', relax_to_this_tau_factor
            call star_relax_tau_factor(id, relax_to_this_tau_factor, dlogtau_factor, ierr)
            if (failed('star_relax_tau_factor')) return
         end if

         if (relax_irradiation) then
            write(*,2) 'relax_irradiation -- min steps', relax_irradiation_min_steps
            write(*,1) 'relax_irradiation -- max yrs dt', relax_irradiation_max_yrs_dt
            call star_relax_irradiation(id, &
               relax_irradiation_min_steps, relax_to_this_irrad_flux, irrad_col_depth, &
               relax_irradiation_max_yrs_dt, ierr)
            if (failed('star_relax_irradiation')) return
         end if

         if (relax_mass_change) then
            write(*,2) 'relax_mass_change -- min steps', relax_mass_change_min_steps
            write(*,1) 'relax_mass_change -- max yrs dt', relax_mass_change_max_yrs_dt
            write(*,1) 'relax_mass_change -- initial_mass_change', relax_mass_change_init_mdot
            write(*,1) 'relax_mass_change -- final_mass_change', relax_mass_change_final_mdot
            call star_relax_mass_change(id, &
               relax_mass_change_min_steps, relax_mass_change_init_mdot, relax_mass_change_final_mdot, &
               relax_mass_change_max_yrs_dt, ierr)
            if (failed('star_relax_mass_change')) return
         end if

         if (s% rotation_flag .and. relax_omega) then
            write(*,*) 'new_omega =', new_omega
            call star_relax_uniform_omega( &
               id, new_omega, num_steps_to_relax_rotation, ierr)
            if (failed('star_relax_uniform_omega')) return
         end if

         if (s% rotation_flag .and. relax_initial_omega .and. .not. restart) then
            write(*,*) 'new_omega =', new_omega
            call star_relax_uniform_omega( &
               id, new_omega, num_steps_to_relax_rotation, ierr)
            if (failed('star_relax_uniform_omega')) return
         end if

         if (s% rotation_flag .and. relax_surface_rotation_v) then
            new_omega = new_surface_rotation_v*1d5/s% r(1)
            write(*,1) 'new_surface_rotation_v =', new_surface_rotation_v, new_omega
            call star_relax_uniform_omega( &
               id, new_omega, num_steps_to_relax_rotation, ierr)
            if (failed('star_relax_uniform_omega')) return
         end if

         if (s% rotation_flag .and. relax_initial_surface_rotation_v .and. .not. restart) then
            new_omega = new_surface_rotation_v*1d5/s% r(1)
            write(*,2) 'new_surface_rotation_v', &
               s% model_number, new_surface_rotation_v, new_omega
            write(*,1) 'new_omega', new_omega
            write(*,*) 'call star_relax_uniform_omega'
            call star_relax_uniform_omega( &
               id, new_omega, num_steps_to_relax_rotation, ierr)
            if (failed('star_relax_uniform_omega')) return
         end if

         if (s% rotation_flag .and. relax_omega_div_omega_crit) then
            new_omega = new_omega_div_omega_crit*star_surface_omega_crit(id, ierr)
            if (failed('star_surface_omega_crit')) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, new_omega_div_omega_crit, new_omega
            call star_relax_uniform_omega( &
               id, new_omega, num_steps_to_relax_rotation, ierr)
            if (failed('star_relax_uniform_omega')) return
         end if

         if (s% rotation_flag .and. relax_initial_omega_div_omega_crit .and. .not. restart) then
            new_omega = new_omega_div_omega_crit*star_surface_omega_crit(id, ierr)
            if (failed('star_surface_omega_crit')) return
            write(*,2) 'new_omega_div_omega_crit', &
               s% model_number, new_omega_div_omega_crit, new_omega
            call star_relax_uniform_omega( &
               id, new_omega, num_steps_to_relax_rotation, ierr)
            if (failed('star_relax_uniform_omega')) return
         end if

         if (set_max_dt_to_frac_lifetime) then
            log_m = log10(s% star_mass) ! in Msun units
            log_lifetime = 9.921 - 3.6648*log_m + 1.9697*log_m**2 - 0.9369*log_m**3
            ! Iben & Laughlin (1989) as quoted in H&K (eqn 2.3)
            max_dt = max_frac_of_lifetime_per_step*secyer*10**(log_lifetime)
            if (max_dt < s% max_timestep) then
               s% max_timestep = max_dt
               write(*, *) 'set_max_dt_to_frac_lifetime: lg(maxdt/secyer)', &
                  log10(s% max_timestep/secyer)
            end if
         end if
         
         if (len_trim(history_columns_file) > 0) &
            write(*,*) 'read ' // trim(history_columns_file)
         call star_set_history_columns(id, history_columns_file, ierr)
         if (failed('star_set_history_columns')) return
         
         if (len_trim(profile_columns_file) > 0) &
            write(*,*) 'read ' // trim(profile_columns_file)
         call star_set_profile_columns(id, profile_columns_file, ierr)
         if (failed('star_set_profile_columns')) return
         
         ! print out info about selected non-standard parameter settings
         
         write(*,*) 'net name ' // trim(s% net_name)
         
         if (len_trim(s% extra_terminal_output_file) > 0) &
            write(*,*) 'extra_terminal_output_file: ' // trim(s% extra_terminal_output_file)
         
         if (s% do_element_diffusion) &
            write(*,*) 'do_element_diffusion', s% do_element_diffusion
         
         if (s% lnPgas_flag) &
            write(*,*) 'lnPgas_flag =', s% lnPgas_flag
         
         if (s% v_flag) &
            write(*,*) 'v_flag =', s% v_flag
         
         if (s% rotation_flag) &
            write(*,*) 'rotation_flag =', s% rotation_flag
         
         if (s% mix_factor /= 1d0) &
            write(*,*) 'mix_factor =', s% mix_factor

         if (s% hydro_numerical_jacobian) &
            write(*,*) 'hydro_numerical_jacobian', s% hydro_numerical_jacobian
            
         if (abs(s% tau_base - 2d0/3d0) > 1d-4) &
            write(*,1) 'tau_base', s% tau_base
            
         if (abs(s% tau_factor - 1) > 1d-4) &
            write(*,1) 'tau_factor', s% tau_factor
            
         if (s% eps_grav_factor /= 1) &
            write(*,1) 'eps_grav_factor', s% eps_grav_factor
            
         if (s% dxdt_nuc_factor /= 1) &
            write(*,1) 'dxdt_nuc_factor', s% dxdt_nuc_factor
            
         if (s% which_atm_option /= 'simple_photosphere') &
            write(*,1) 'which_atm_option: ' // trim(s% which_atm_option)
           
         if (s% M_center /= 0) then
            write(*,1) 'xmstar/mstar', s% xmstar/s% mstar
            write(*,1) 'xmstar (g)', s% xmstar
            write(*,1) 'M_center (g)', s% M_center
            write(*,1) 'xmstar/Msun', s% xmstar/Msun
            write(*,1) 'M_center/Msun', s% M_center/Msun
         end if
            
         if (s% R_center /= 0) then
            write(*,1) 'R_center (cm)', s% R_center
            write(*,1) 'R_center/Rsun', s% R_center/Rsun
            write(*,1) 'core density', s% M_center/(4*pi/3*s% R_center**3)
         end if
                     
         if (s% opacity_max > 0) &
            write(*,*) 'opacity_max', s% opacity_max
            
         if (s% operator_coupling_choice /= 0) &
            write(*,2) 'operator_coupling_choice', s% operator_coupling_choice
         
         if (show_net_reactions_info) then
            write(*,'(a)') ' net reactions '
            call show_net_reactions_and_info(s% net_handle, 6, ierr)
            if (failed('show_net_reactions_and_info')) return
         end if
         
         if (s% species > s% hydro_decsol_switch &
               .and. s% operator_coupling_choice == 0) then
            write(*,3) 'use large_mtx_decsol ' // trim(s% large_mtx_decsol), &
               s% species, s% hydro_decsol_switch
         else
            write(*,3) 'use small_mtx_decsol ' // trim(s% small_mtx_decsol), &
               s% species, s% hydro_decsol_switch
         end if
         
         if (show_net_species_info) then
            write(*,'(a)') ' species'
            do j=1,s% species
               write(*,'(i6,3x,a)') j, chem_isos% name(s% chem_id(j))
            end do
            write(*,*)
         end if
         
         if (show_eqns_and_vars_names) then
            do i=1,s% nvar
               write(*,*) i, s% nameofvar(i), s% nameofequ(i)
            end do
            write(*,*)
         end if         
         
         write(*,'(a)') ' kappa_file_prefix ' // trim(kappa_file_prefix)
         if (len_trim(other_kappa_file_prefix) > 0) &
            write(*,'(a)') ' other_kappa_file_prefix ' // trim(other_kappa_file_prefix)
         write(*,'(a)') ' kappa_lowT_prefix ' // trim(kappa_lowT_prefix)
         
         write(*,'(a)') '   eos_file_prefix ' // trim(eos_file_prefix)
         if (len_trim(other_eos_file_prefix) > 0) &
            write(*,'(a)') ' other_eos_file_prefix ' // trim(other_eos_file_prefix)

         write(*,2) 'OMP_NUM_THREADS', omp_get_max_threads()

         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed


         subroutine do_relax_initial_composition(ierr)
            use utils_lib
            integer, intent(out) :: ierr
            real(dp), pointer :: xq(:), xa(:,:)
            integer :: num_pts, num_species, i, iounit
            include 'formats.inc'
            
            write(*,*)
            write(*,1) 'relax_initial_composition'

            iounit = alloc_iounit(ierr)
            if (ierr /= 0) return
            open(unit=iounit, file=trim(relax_composition_filename), &
                  status='old', action='read', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'open failed', ierr, iounit
               write(*, '(a)') 'failed to open ' // trim(relax_composition_filename)
               call free_iounit(iounit)
               return
            end if
            read(iounit, *, iostat=ierr) num_pts, num_species
            if (ierr /= 0) then
               close(iounit)
               call free_iounit(iounit)
               write(*, '(a)') &
                  'failed while trying to read 1st line of ' // trim(relax_composition_filename)
               return
            end if
            allocate(xq(num_pts), xa(num_species,num_pts))
            do i = 1, num_pts
               read(iounit,*,iostat=ierr) xq(i), xa(1:num_species,i)
               if (ierr /= 0) then
                  close(iounit)
                  call free_iounit(iounit)
                  write(*, '(a)') &
                     'failed while trying to read ' // trim(relax_composition_filename)
                  write(*,*) 'line', i+1
                  deallocate(xq,xa)
                  return
               end if
            end do
            close(iounit)
            call free_iounit(iounit)
            call star_relax_composition( &
               id, num_steps_to_relax_composition, num_pts, num_species, xa, xq, ierr)
            deallocate(xq,xa)
         end subroutine do_relax_initial_composition
         
         
         subroutine set_net_T_limits
            use net_def
            type (Net_General_Info), pointer :: g
            call get_net_ptr(s% net_handle, g, ierr)
            if (ierr /= 0) return
     
            g% T_lo_2n = T_lo_2n
            g% T_hi_2n = T_hi_2n
            g% lnT_lo_2n = log(g% T_lo_2n)
            g% lnT_hi_2n = log(g% T_hi_2n)
            g% min_factor_2n = min_factor_2n
            g% min_ln_factor_2n = log(g% min_factor_2n)

            g% T_lo_prot = T_lo_prot
            g% T_hi_prot = T_hi_prot
            g% lnT_lo_prot = log(g% T_lo_prot)
            g% lnT_hi_prot = log(g% T_hi_prot)
            g% min_factor_prot = min_factor_prot
            g% min_ln_factor_prot = log(g% min_factor_prot)
     
            g% T_lo_combo_a_capture = T_lo_combo_a_capture
            g% T_hi_combo_a_capture = T_hi_combo_a_capture
            g% lnT_lo_combo_a_capture = log(g% T_lo_combo_a_capture)
            g% lnT_hi_combo_a_capture = log(g% T_hi_combo_a_capture)
            g% min_factor_combo_a_capture = min_factor_combo_a_capture
            g% min_ln_factor_combo_a_capture = log(g% min_factor_combo_a_capture)
     
            g% T_lo_a_cap_high_mass = T_lo_a_cap_high_mass
            g% T_hi_a_cap_high_mass = T_hi_a_cap_high_mass
            g% lnT_lo_a_cap_high_mass = log(g% T_lo_a_cap_high_mass)
            g% lnT_hi_a_cap_high_mass = log(g% T_hi_a_cap_high_mass)
            g% min_factor_a_cap_high_mass = min_factor_a_cap_high_mass
            g% min_ln_factor_a_cap_high_mass = log(g% min_factor_a_cap_high_mass)
     
            g% T_lo_a_cap_intermediate = T_lo_a_cap_intermediate
            g% T_hi_a_cap_intermediate = T_hi_a_cap_intermediate
            g% lnT_lo_a_cap_intermediate = log(g% T_lo_a_cap_intermediate)
            g% lnT_hi_a_cap_intermediate = log(g% T_hi_a_cap_intermediate)
            g% min_factor_a_cap_intermediate = min_factor_a_cap_intermediate
            g% min_ln_factor_a_cap_intermediate = log(g% min_factor_a_cap_intermediate)
         
            g% max_factor_2n = max_factor_2n
            
            write(*,*) 'set_net_T_limits'
         
         end subroutine set_net_T_limits
         

      end subroutine do_star_job_controls_after
      

      subroutine extend_net(s, ierr)
         use net_def
         use chem_def
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         real(dp), parameter :: tiny = 1d-10, small = 1d-2
         
         real(dp) :: cntr_h, cntr_he
         
         include 'formats.inc'
         
         ierr = 0
         
         if (s% net_name == adv_net) return

         if (s% net_name == co_net) then
            if (s% log_center_temperature > 8.8d0 .or. s% log_center_density > 9d0) then
               call change_net(adv_net)
               if (len_trim(profile_columns_file) > 0) &
                  write(*,*) 'read ' // trim(profile_columns_file)
               call star_set_profile_columns(s% id, profile_columns_file, ierr)
            end if
            return
         end if
         
         if (s% net_name == h_he_net) then
            cntr_h = current_abundance_at_point(s% id, ih1, s% nz, ierr)
            if (ierr /= 0) return
            if (cntr_h > tiny) return
            cntr_he = current_abundance_at_point(s% id, ihe4, s% nz, ierr)
            if (ierr /= 0) return
            if (cntr_he > small) return
            if (s% log_center_temperature > 8.3d0 .or. s% log_center_density > 8.5d0) then
               call change_net(co_net)
               if (len_trim(profile_columns_file) > 0) &
                  write(*,*) 'read ' // trim(profile_columns_file)
               call star_set_profile_columns(s% id, profile_columns_file, ierr)
            end if
         end if

         
         contains
                  
         
         subroutine change_net(net_name)
            use const_def
            character (len=*), intent(in) :: net_name
            integer :: j
            
            call star_change_to_new_net(s% id, net_name, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in star_change_to_new_net ' // trim(net_name)
               stop 'change_net'
               return
            end if
            
            if (net_name /= s% net_name) then
               write(*,*) '   new net_name ', trim(net_name)
               write(*,*) 'old s% net_name ', trim(s% net_name)
               write(*,*) 'failed to change'
               stop 'change_net'
            end if

            write(*,'(a)') ' new net = ' // trim(s% net_name)
            !do j=1,s% species
            !   write(*,fmt='(a,x)',advance='no') trim(chem_isos% name(s% chem_id(j)))
            !end do
            !write(*,*)
            s% dt_next = s% dt_next/5
            write(*,*) 'reduce timestep', log10(s% dt_next/secyer)
            write(*,*)
         end subroutine change_net
         
         
      end subroutine extend_net         


      subroutine before_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine before_evolve         
      
                  
      subroutine write_colors_info(id, s, ierr)
         use colors_lib
         use colors_def
         use chem_def, only: zsol
         use num_lib, only: safe_log10
         use utils_lib, only: alloc_iounit, free_iounit
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         integer :: io, i, j
         character (len=strlen) :: fname
			real(dp)  :: log_Teff ! log10 of surface temp
			real(dp)  :: log_L ! log10 of luminosity in solar units
			real(dp)  :: mass ! mass in solar units
			real(dp)  :: Fe_H ! [Fe/H]
			! output
			real(dp) :: results(n_colors)
			real(dp) :: log_g
         
         ierr = 0
         
         io = alloc_iounit(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in alloc_iounit'
            return
         end if
         
         fname = 'colors.log'
         !if (s% doing_first_model_of_run) then
         if (.false.) then
            open(unit=io, file=trim(fname), action='write', status='replace', iostat=ierr)
            ! write column numbers
            j = 1
            write(io,fmt='(i10)',advance='no') j
            j = j+1
            do i=1,4+n_colors
               write(io,fmt='(i25)',advance='no') j
               j = j+1
            end do
            write(io,fmt='(i25)') j
            ! write column labels
            write(io,fmt='(a10)',advance='no') 'model'
            write(io,fmt='(a25)',advance='no') 'log_Teff'
            write(io,fmt='(a25)',advance='no') 'log_L'
            write(io,fmt='(a25)',advance='no') 'mass'
            write(io,fmt='(a25)',advance='no') 'Fe_H'
            do i=1,n_colors
               write(io,fmt='(a25)',advance='no') trim(colors_name(i))
            end do
            write(io,fmt='(a25)') 'log_g'
         else
            open(unit=io, file=trim(fname), action='write', position='append', iostat=ierr)
         end if
         if (ierr /= 0) then
            write(*,*) 'failed to open colors.log'
            call free_iounit(io)
            return
         end if
         
         log_Teff = log10(s% Teff)
         log_L = s% log_surface_luminosity
         mass = s% star_mass
         Fe_H = safe_log10(get_current_z_at_point(id, 1, ierr) / zsol)
         if (ierr /= 0) then
            write(*,*) 'failed in get_current_z_at_point'
            call cleanup
            return
         end if
         
         call colors_get(log_Teff, log_L, mass, Fe_H, results, log_g, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in colors_get'
            call cleanup
            return
         end if
         
         1 format(1x,f24.12)
         write(io,fmt='(i10)',advance='no') s% model_number
         write(io,fmt=1,advance='no') log_Teff
         write(io,fmt=1,advance='no') log_L
         write(io,fmt=1,advance='no') mass
         write(io,fmt=1,advance='no') Fe_H
         do i=1,n_colors
            write(io,fmt=1,advance='no') results(i)
         end do
         write(io,1) log_g
         
         call cleanup
         
         contains
         
         subroutine cleanup
            close(io)
            call free_iounit(io)
         end subroutine cleanup
      
      end subroutine write_colors_info
      
      
      subroutine do_read_star_job(filename, ierr)
         character(*), intent(in) :: filename
         integer, intent(out) :: ierr  
         
         ierr = 0   
         
         include "star_job.defaults"
         
         ! here are a few "private" controls that aren't intended for general users (yet/ever)
   
         nzlo = -1
         nzhi = -1  ! if nzhi = 0, then uses nz, the current number of zones.
   
         ! set mass fraction of chem_id to new_frac uniformly in cells nzlo to nzhi.
         set_abundance = .false.
         chem_id = -1 ! a chem_id such as ihe4.  see chem_def.
         new_frac = -1
   
         ! replace chem1 by chem2 in cells nzlo to nzhi.
         replace_element = .false.
         chem_id1 = -1
         chem_id2 = -1
   
         do_special_test = .false.
         
         ierr = 0
         call read_inlist(filename, 1, ierr)
         
         if (ierr /= 0) then
            write(*,*) 'ierr from read_inlist ' // trim(filename)
            return
         end if
         
         if (save_star_job_namelist) &
            call write_controls(star_job_namelist_name, ierr)

      end subroutine do_read_star_job


      subroutine write_controls(filename_in, ierr)
         use utils_lib
         character(*), intent(in) :: filename_in
         integer, intent(out) :: ierr
         character (len=256) :: message, filename
         integer :: unit
         ierr = 0
         filename = trim(filename_in)
         if (len_trim(filename) == 0) filename = 'star_job_namelist.out'
         unit=alloc_iounit(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to alloc iounit in write_controls'
            return
         end if
         ! NOTE: when open namelist file, must include delim='APOSTROPHE'
         open(unit=unit, file=trim(filename), action='write', delim='APOSTROPHE', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open ', trim(filename)
         else
            !write(*,*) 'call write(unit, nml=star_job)'
            write(unit, nml=star_job)
            !write(*,*) 'done write(unit, nml=star_job)'
            close(unit)
         end if
         call free_iounit(unit)
         write(*,*) 'saved initial &star_job inlist values: ' // trim(filename)

      end subroutine write_controls
      
      
      recursive subroutine read_inlist(filename, level, ierr)
         use utils_lib
         character(*), intent(in) :: filename
         integer, intent(in) :: level  
         integer, intent(out) :: ierr  
         
         logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
         character (len=256) :: message, extra1, extra2, extra3, extra4, extra5
         integer :: unit
         
         if (level >= 10) then
            write(*,*) 'ERROR: too many levels of nested extra star_job inlist files'
            ierr = -1
            return
         end if
         
         ierr = 0
         unit=alloc_iounit(ierr)
         if (ierr /= 0) return
         
         open(unit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=star_job, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(unit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=star_job)
               close(unit)
            end if  
         end if
         call free_iounit(unit)
         if (ierr /= 0) return
         
         ! recursive calls to read other inlists
         
         read_extra1 = read_extra_star_job_inlist1
         read_extra_star_job_inlist1 = .false.
         extra1 = extra_star_job_inlist1_name
         extra_star_job_inlist1_name = 'undefined'
         
         read_extra2 = read_extra_star_job_inlist2
         read_extra_star_job_inlist2 = .false.
         extra2 = extra_star_job_inlist2_name
         extra_star_job_inlist2_name = 'undefined'
         
         read_extra3 = read_extra_star_job_inlist3
         read_extra_star_job_inlist3 = .false.
         extra3 = extra_star_job_inlist3_name
         extra_star_job_inlist3_name = 'undefined'
         
         read_extra4 = read_extra_star_job_inlist4
         read_extra_star_job_inlist4 = .false.
         extra4 = extra_star_job_inlist4_name
         extra_star_job_inlist4_name = 'undefined'
         
         read_extra5 = read_extra_star_job_inlist5
         read_extra_star_job_inlist5 = .false.
         extra5 = extra_star_job_inlist5_name
         extra_star_job_inlist5_name = 'undefined'
         
         if (read_extra1) then
            !write(*,*) 'read extra star_job inlist1 from ' // trim(extra1)
            call read_inlist(extra1, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra2) then
            !write(*,*) 'read extra star_job inlist2 from ' // trim(extra2)
            call read_inlist(extra2, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra3) then
            !write(*,*) 'read extra star_job inlist3 from ' // trim(extra3)
            call read_inlist(extra3, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra4) then
            !write(*,*) 'read extra star_job inlist4 from ' // trim(extra4)
            call read_inlist(extra4, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra5) then
            write(*,*) 'read extra star_job inlist5 from ' // trim(extra5)
            call read_inlist(extra5, level+1, ierr)
            if (ierr /= 0) return
         end if
         
      end subroutine read_inlist
      
      
      subroutine read_masses(filename, masses, nmasses, ierr)
         character (len=*), intent(in) :: filename
         real(dp), pointer, intent(out) :: masses(:)
         integer, intent(out) :: nmasses, ierr
         call read_items(filename, masses, nmasses, 'masses', ierr)
      end subroutine read_masses
      
      
      subroutine read_items(filename, items, nitems, name, ierr)
         use utils_lib
         use utils_def
         character (len=*), intent(in) :: filename, name
         real(dp), pointer, intent(out) :: items(:)
         integer, intent(out) :: nitems, ierr
         
         integer :: iounit, n, i, t, capacity
         character (len=256) :: buffer, string
         
         nitems = 0
         if (.not. associated(items)) then
            capacity = 10
            allocate(items(capacity))
         else
            capacity = size(items,dim=1)
         end if
         
         ierr = 0
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return

         open(unit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            call free_iounit(iounit)
            write(*,*) 'failed to open file ' // trim(filename)
            return
         end if
         
         n = 0
         i = 0
         
         do
            t = token(iounit, n, i, buffer, string)
            select case(t)
               case(name_token)
                  if (string == name) then
                     call do_read_items(ierr)
                     if (ierr /= 0) then
                        call free_iounit(iounit)
                        return
                     end if
                     exit ! for now, nothing else to be read
                  end if
                  call error; return
               case(eof_token)
                  exit
               case default
                  call error; return
            end select
            
         end do
         
         close(iounit)
         call free_iounit(iounit)
         
         contains
         
         
         subroutine error
            ierr = -1
            write(*,*) 'error in reading file' // trim(filename)
            close(iounit)
            call free_iounit(iounit)
         end subroutine error
         
         
         subroutine do_read_items(ierr)
            integer, intent(out) :: ierr
            real(dp) :: mass
            ierr = 0
            t = token(iounit, n, i, buffer, string)
            if (t /= left_paren_token) then
               call error; return
            end if
         mass_loop: do
               t = token(iounit, n, i, buffer, string)
               if (t /= name_token) then
                  call error; return
               end if
               read(string,fmt=*,iostat=ierr) mass
               if (ierr /= 0) then
                  call error; return
               end if
               nitems = nitems+1
               if (nitems > capacity) then
                  capacity = capacity + 10
                  call realloc_double(items,capacity,ierr)
                  if (ierr /= 0) then
                     call error; return
                  end if
               end if
               items(nitems) = mass
               t = token(iounit, n, i, buffer, string)
               if (t == right_paren_token) exit mass_loop
               if (t /= comma_token) then
                  call error; return
               end if
            end do mass_loop
         end subroutine do_read_items
         
      
      end subroutine read_items
      
      
      subroutine do_report_mass_not_fe56(s)
         use const_def
         type (star_info), pointer :: s
         integer :: k, fe56
         real(dp) :: sumdq
         include 'formats.inc'
         fe56 = s% net_iso(ife56)
         if (fe56 == 0) return
         sumdq = 0
         do k = 1, s% nz
            sumdq = sumdq + s% dq(k)*(1-s% xa(fe56,k))
         end do
         write(*,1) 'R', s% r(1)
         write(*,1) 'g', s% cgrav(1)*s% mstar/(s% r(1)**2)
         write(*,1) 'mass non fe56', s% xmstar*sumdq, sumdq
         write(*,1) 'M_center (Msun)', s% M_center/Msun
         write(*,1) 'xmstar (g)', s% xmstar
         do k=1,s% nz
            if (fe56 == maxloc(s% xa(:,k),dim=1)) then
               write(*,2) 'mass exterior to fe56 (g)', k, (1d0 - s% q(k))*s% xmstar
               write(*,2) 'mass coord top of fe56 (g)', k, s% q(k)*s% xmstar
               return
            end if
         end do
      end subroutine do_report_mass_not_fe56
      
      
      subroutine do_report_cell_for_xm(s)
         use const_def
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: sumdq, dq
         include 'formats.inc'
         dq = report_cell_for_xm/s% xmstar
         if (dq > 1) then
            write(*,2) 'report_cell_for_xm > xmstar', s% nz
            return
         end if
         sumdq = 0
         do k = 1, s% nz
            sumdq = sumdq + s% dq(k)
            if (sumdq >= dq) then
               write(*,*)
               write(*,2) 'total mass in cells from 1 to k', k, sumdq*s% xmstar
               write(*,2) 'logT(k)', k, s% lnT(k)/ln10
               write(*,2) 'logRho(k)', k, s% lnd(k)/ln10
               write(*,2) 'entropy(k)', k, exp(s% lnS(k))*amu/kerg
               write(*,2) 'xmstar*q(k)', k, s% xmstar*s% q(k)
               write(*,2) 'q(k)', k, s% q(k)
               write(*,*)
               return
            end if
         end do
         write(*,2) 'total mass in cells from 1 to nz', s% nz, s% xmstar
      end subroutine do_report_cell_for_xm
      
      
      subroutine do_star_init(ierr)
         integer, intent(out) :: ierr
         call star_init( &
            mesa_dir, chem_isotopes_filename, &
            kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix, &
            kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy, &
            kappa_type2_logT_lower_bdy, other_kappa_file_prefix, &
            eos_file_prefix, other_eos_file_prefix, eosDT_Z1_suffix, eosPT_Z1_suffix, &
            net_reaction_filename, rate_tables_dir, rate_cache_suffix, &
            ionization_file_prefix, ionization_Z1_suffix, &
            eosDT_cache_dir, eosPT_cache_dir, ionization_cache_dir, &
            kap_cache_dir, rates_cache_dir, weaklib_cache_dir, &
            ierr)
      end subroutine do_star_init
      
      
      subroutine init_and_alloc(id, s, ierr)
         integer, intent(out) :: id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         call do_star_init(ierr)
         if (failed('star_init')) return
         
         id = alloc_star(ierr)
         if (failed('alloc_star')) return
         
         call star_ptr(id, s, ierr)
         if (failed('star_ptr')) return
         
         contains
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed

      end subroutine init_and_alloc       


      subroutine run1_star( &
            do_init_and_alloc, do_free_star, okay_to_restart, &
            id, restart, &
            extras_controls, &
            extras_startup, &
            extras_check_model, &
            how_many_extra_history_columns, &
            data_for_extra_history_columns, &
            how_many_extra_profile_columns, &
            data_for_extra_profile_columns, &
            extras_finish_step, &
            extras_after_evolve, &
            ierr, &
            inlist_fname_arg)
         use net_def
         use chem_def
         use const_def
         use num_lib, only: safe_log10
         use se_support, only: se_startup, se_finish_step, se_after_evolve
         use utils_lib, only: number_iounits_allocated
         
         logical, intent(in) :: do_init_and_alloc, do_free_star, okay_to_restart
         integer, intent(inout) :: id ! input if not do_init_and_alloc
         logical, intent(inout) :: restart ! input if not do_init_and_alloc
         character (len=32) :: inlist_fname_arg
         integer, intent(out) :: ierr
         optional inlist_fname_arg
         
         interface

            subroutine extras_controls(s, ierr)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(out) :: ierr
            end subroutine extras_controls      
     
            integer function extras_startup(s, id, restart, ierr)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id
               logical, intent(in) :: restart
               integer, intent(out) :: ierr
            end function extras_startup
      
            integer function extras_check_model(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
            end function extras_check_model

            integer function how_many_extra_history_columns(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
            end function how_many_extra_history_columns
            
            subroutine data_for_extra_history_columns(s, id, id_extra, n, names, vals, ierr)
               use const_def, only: dp
               use star_def, only: star_info, maxlen_history_column_name
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra, n
               character (len=maxlen_history_column_name) :: names(n)
               real(dp) :: vals(n)
               integer, intent(out) :: ierr
            end subroutine data_for_extra_history_columns
      
            integer function how_many_extra_profile_columns(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
            end function how_many_extra_profile_columns      
      
            subroutine data_for_extra_profile_columns(s, id, id_extra, n, nz, names, vals, ierr)
               use const_def, only: dp
               use star_def, only: star_info, maxlen_profile_column_name
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra, n, nz
               character (len=maxlen_profile_column_name) :: names(n)
               real(dp) :: vals(nz,n)
               integer, intent(out) :: ierr
            end subroutine data_for_extra_profile_columns
      
            integer function extras_finish_step(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
               integer :: ierr
            end function extras_finish_step     
      
            subroutine extras_after_evolve(s, id, id_extra, ierr)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
               integer, intent(out) :: ierr
            end subroutine extras_after_evolve

         end interface
         
         integer :: id_extra, result, result_reason, model_number, k, source, nsteps, j, ci, nz
         integer :: time0, time1, clock_rate, time0_extra, time1_extra, prev_num_iounits_in_use
         logical :: first_try, just_did_backup, continue_evolve_loop
         type (star_info), pointer :: s
         real(dp) :: gamma1_integral, integral_norm
         character (len=32) :: inlist_fname
      
         integer :: num_log_columns
         integer :: num_profile_columns
         character (len=maxlen_history_column_name), pointer :: log_names(:) ! num_log_columns
         double precision, pointer :: log_vals(:)
         logical, pointer :: log_is_int(:) ! true if the values in the log column are integers
         character (len=maxlen_profile_column_name), pointer :: profile_names(:) ! num_profiles_columns
         double precision, pointer :: profile_vals(:,:) ! (nz,num_profile_columns)
         logical, pointer :: profile_is_int(:) ! true if the values in the profile column are integers
   
         1 format(a35, 99(1pe26.16))
         2 format(a35, i7, 1pe26.16)
         3 format(a15, 2x, f15.6)
         4 format(a15, 2x, e15.6)

         11 format(a35, f20.10)

         ierr = 0
         result_reason = result_reason_normal
         prev_num_iounits_in_use = -1
         
         if (present(inlist_fname_arg)) then
            inlist_fname = inlist_fname_arg
         else
            inlist_fname = 'inlist'
         end if
         
         if (do_init_and_alloc) then  
            call init_and_alloc(id, s, ierr)
            if (failed('init_and_alloc')) return
         else
            call star_ptr(id, s, ierr)
            if (failed('star_ptr')) return
         end if
         
         call star_setup(id, inlist_fname, ierr)
         if (failed('star_setup')) return

         if (okay_to_restart) then
            restart = doing_a_restart()
         else
            restart = .false.
         end if
         
         if (show_log_description_at_start .and. .not. restart) then
            write(*,*)
            call show_log_description(id, ierr)
            if (failed('show_log_description')) return
         end if

         call extras_controls(s, ierr)
         if (ierr /= 0) return
         
         call do_star_job_controls_before(id, s, restart, ierr)
         if (ierr /= 0) return

         call do_load1_star(id, s, restart, 'restart_photo', ierr)
         if (failed('do_load1_star')) return
         
         call do_star_job_controls_after(id, s, restart, ierr)
         if (failed('do_star_job_controls_after')) return

         write(*,*)
         write(*,*)
         
         if (.not. restart) then
            call before_evolve(id, ierr)
            if (failed('before_evolve')) return
            call start_new_run_for_pgstar(s, ierr)
            if (failed('start_new_run_for_pgstar')) return
         else
            call show_terminal_header(id, ierr)
            if (failed('show_terminal_header')) return
            call restart_run_for_pgstar(s, ierr)
            if (failed('restart_run_for_pgstar')) return
         end if
         
         id_extra = extras_startup(s, id, restart, ierr)
         if (failed('extras_startup')) return
         
         ! Get names and values for history
         num_log_columns = num_standard_history_columns(s) + how_many_extra_history_columns(s, id, id_extra)
         allocate(                            &
              log_names(num_log_columns),     &
              log_vals(num_log_columns),      &
              log_is_int(num_log_columns),    &
              stat=ierr)
         call get_data_for_history_columns(s, id_extra, &
              how_many_extra_history_columns, data_for_extra_history_columns, &
              num_log_columns,s% nz, log_names, log_vals, log_is_int, ierr)
         ! Get names and values from profiles
         num_profile_columns = num_standard_profile_columns(s) + how_many_extra_profile_columns(s, id, id_extra)
         allocate(                                     &
              profile_names(num_profile_columns),      &
              profile_vals(s% nz,num_profile_columns), &
              profile_is_int(num_profile_columns),     &
              stat=ierr)
         call get_data_for_profile_columns(s, id_extra, &
              how_many_extra_profile_columns, data_for_extra_profile_columns, &
              num_profile_columns,s% nz, profile_names, profile_vals, profile_is_int, ierr)

         call se_startup(s, restart, use_se_output,                             &
              se_codev, se_modname, se_prefix, se_num_mod_output,               &
              log_names, log_vals, log_is_int, num_log_columns,                 &
              profile_names, profile_vals, profile_is_int, num_profile_columns, &
              ierr)

         deallocate(log_names,log_vals,log_is_int)
         deallocate(profile_names,profile_vals,profile_is_int)

         if (failed('se_startup')) return

         if (profile_starting_model) then
            write(*, '(a, i12)') 'save profile for model number ', s% model_number
            !call star_set_vars(id, 1, s% nz, ierr)
            if (failed('star_set_vars')) return
            call save_profile(id, id_extra, &
               how_many_extra_profile_columns, data_for_extra_profile_columns, &
               3, ierr)
            if (failed('save_profile')) return
         end if

         continue_evolve_loop = .true.
         s% doing_timing = .false.
         check_before_step_timing = 0
         check_step_loop_timing = 0
         check_after_step_timing = 0
         
         evolve_loop: do while(continue_evolve_loop) ! evolve one step per loop
            
            if (first_model_for_timing >= 0 .and. &
                  s% model_number >= first_model_for_timing .and. &
                  .not. s% doing_timing) then
               s% doing_timing = .true.
               write(*,*) 'start timing'
               write(*,*)
               call system_clock(time0,clock_rate)
               step_loop_timing = 0
               after_step_timing = 0
               before_step_timing = 0
            end if
            
            if (s% doing_timing) then
               call system_clock(time0_extra,clock_rate)
               check_time_start = eval_total_times(s% id, ierr)
            end if
            
            if (auto_extend_net) then
               call extend_net(s, ierr)
               if (failed('extend_net')) return
            end if
            
            if (s% center_ye <= center_ye_limit_for_v_flag .and. .not. s% v_flag) then
               write(*,1) 'have reached center ye limit', s% center_ye, center_ye_limit_for_v_flag
               write(*,1) 'set v_flag true'
               call star_set_v_flag(id, .true., ierr)
               if (failed('star_set_v_flag')) return
               if (ierr /= 0) return
            end if
            
            if (s% log_center_temperature > 9d0 .and. .not. s% v_flag) then 
               ! thanks go to Roni Waldman for this
               gamma1_integral = 0
               integral_norm = 0
               do k=1,s% nz
                  integral_norm = integral_norm + s% P(k)*s% dm(k)/s% rho(k)
                  gamma1_integral = gamma1_integral + &
                     (s% gamma1(k)-4.d0/3.d0)*s% P(k)*s% dm(k)/s% rho(k)
               end do
               gamma1_integral = gamma1_integral/max(1d-99,integral_norm)
               if (gamma1_integral <= gamma1_integral_for_v_flag) then
                  write(*,1) 'have reached gamma1 integral limit', gamma1_integral
                  write(*,1) 'set v_flag true'
                  call star_set_v_flag(id, .true., ierr)
                  if (failed('star_set_v_flag')) return
                  if (ierr /= 0) return
               end if
            end if
         
            if (report_mass_not_fe56) call do_report_mass_not_fe56(s)
            if (report_cell_for_xm > 0) call do_report_cell_for_xm(s)
         
            first_try = .true.
            just_did_backup = .false.
            
            model_number = get_model_number(id, ierr)
            if (failed('get_model_number')) return
            
            if (s% doing_timing) then
            
               call system_clock(time1_extra,clock_rate)
               before_step_timing = before_step_timing + dble(time1_extra - time0_extra) / clock_rate
               
               check_time_end = eval_total_times(s% id, ierr)
               check_before_step_timing = check_before_step_timing + (check_time_end - check_time_start)

               time0_extra = time1_extra
               check_time_start = check_time_end

            end if

            step_loop: do ! may need to repeat this loop
            
               if (stop_is_requested()) then
                  continue_evolve_loop = .false.
                  result = terminate
                  exit
               end if
            
               result = star_evolve_step(id, first_try, just_did_backup)
               if (result == keep_going) result = star_check_model(id)
               if (result == keep_going) result = extras_check_model(s, id, id_extra)
               if (result == keep_going) result = star_pick_next_timestep(id)            
               if (result == keep_going) exit step_loop
               
               model_number = get_model_number(id, ierr)
               if (failed('get_model_number')) return
                              
               result_reason = get_result_reason(id, ierr)
               if (failed('get_result_reason')) return
               if (result == retry .and. report_retries) then
                  write(*,'(i6,3x,a,/)') model_number, &
                     'retry reason ' // trim(result_reason_str(result_reason))
               else if (result == backup .and. report_backups) then
                  write(*,'(i6,3x,a,/)') model_number, &
                     'backup reason ' // trim(result_reason_str(result_reason))
               end if
               
               if (result == redo) then
                  result = star_prepare_to_redo(id)
               end if
               if (result == retry) then
                  result = star_prepare_to_retry(id)
               end if
               if (result == backup) then
                  result = star_do1_backup(id)
                  just_did_backup = .true.
               end if
               if (result == terminate) then
                  continue_evolve_loop = .false.
                  exit step_loop
               end if
               first_try = .false.
               
            end do step_loop
            
            if (s% doing_timing) then
            
               call system_clock(time1_extra,clock_rate)
               step_loop_timing = step_loop_timing + dble(time1_extra - time0_extra) / clock_rate
               
               check_time_end = eval_total_times(s% id, ierr)
               check_step_loop_timing = check_step_loop_timing + (check_time_end - check_time_start)

               time0_extra = time1_extra
               check_time_start = check_time_end
               
            end if
         
            if (set_tau_factor_after_core_He_burn > 0 .and. &
                  abs(s% tau_factor - set_to_this_tau_factor) > &
                     1d-6*max(s% tau_factor, set_to_this_tau_factor)) then
               if (check_for_after_He_burn(s, set_tau_factor_after_core_He_burn)) then
                  s% tau_factor = set_to_this_tau_factor
                  write(*,1) 'set_tau_factor_after_core_He_burn', s% tau_factor
               end if
            end if
         
            if (set_tau_factor_after_core_C_burn > 0 .and. &
                  abs(s% tau_factor - set_to_this_tau_factor) > &
                     1d-6*max(s% tau_factor, set_to_this_tau_factor)) then
               if (check_for_after_C_burn(s, set_tau_factor_after_core_C_burn)) then
                  s% tau_factor = set_to_this_tau_factor
                  write(*,1) 'set_tau_factor_after_core_C_burn', s% tau_factor
               end if
            end if
         
            if (relax_tau_factor_after_core_He_burn > 0 .and. &
                  abs(s% tau_factor - relax_to_this_tau_factor) > &
                     1d-6*max(s% tau_factor, relax_to_this_tau_factor)) then
               if (check_for_after_He_burn(s, relax_tau_factor_after_core_He_burn)) &
                  call relax_tau_factor
            end if
         
            if (relax_tau_factor_after_core_C_burn > 0 .and. &
                  abs(s% tau_factor - relax_to_this_tau_factor) > &
                     1d-6*max(s% tau_factor, relax_to_this_tau_factor)) then
               write(*,*) 'call check_for_after_C_burn'
               if (check_for_after_C_burn(s, relax_tau_factor_after_core_C_burn)) &
                  call relax_tau_factor
            end if
            
            if (s% L_nuc_burn_total/s% L_phot >= s% Lnuc_div_L_zams_limit &
                  .and. .not. s% rotation_flag) then      
                         
               if (set_near_zams_surface_rotation_v_steps > 0) then
                  new_rotation_flag = .true.
                  call star_set_rotation_flag(id, new_rotation_flag, ierr)
                  if (failed('star_set_rotation_flag')) return
                  set_surf_rotation_v_step_limit = &
                     s% model_number + set_near_zams_surface_rotation_v_steps - 1
                  write(*,2) 'near zams: set_surf_rotation_v_step_limit', set_surf_rotation_v_step_limit

               else if (set_near_zams_omega_steps > 0) then
                  new_rotation_flag = .true.
                  call star_set_rotation_flag(id, new_rotation_flag, ierr)
                  if (failed('star_set_rotation_flag')) return
                  set_omega_step_limit = &
                     s% model_number + set_near_zams_omega_steps - 1
                  write(*,2) 'near zams: set_omega_step_limit', set_omega_step_limit

               else if (set_near_zams_omega_div_omega_crit_steps > 0) then
                  new_rotation_flag = .true.
                  call star_set_rotation_flag(id, new_rotation_flag, ierr)
                  if (failed('star_set_rotation_flag')) return
                  set_omega_div_omega_crit_step_limit = &
                     s% model_number + set_near_zams_omega_div_omega_crit_steps - 1
                  write(*,2) 'near zams: set_omega_div_omega_crit_step_limit', &
                     set_omega_div_omega_crit_step_limit

               else if (near_zams_relax_omega) then
                  new_rotation_flag = .true.
                  call star_set_rotation_flag(id, new_rotation_flag, ierr)
                  if (failed('star_set_rotation_flag')) return
                  write(*,*) 'new_omega =', new_omega
                  call star_relax_uniform_omega( &
                     id, new_omega, num_steps_to_relax_rotation, ierr)
                  if (failed('star_relax_uniform_omega')) return

               else if (near_zams_relax_initial_surface_rotation_v) then
                  new_rotation_flag = .true.
                  call star_set_rotation_flag(id, new_rotation_flag, ierr)
                  if (failed('star_set_rotation_flag')) return
                  new_omega = new_surface_rotation_v*1d5/s% r(1)
                  write(*,1) 'new_surface_rotation_v', new_surface_rotation_v
                  write(*,1) 's% r(1)/Rsun', s% r(1)/Rsun
                  write(*,1) 'new_omega', new_omega
                  write(*,2) 'new_surface_rotation_v', &
                     s% model_number, new_surface_rotation_v
                  write(*,*) 'near_zams_relax_initial_surface_rotation_v: call star_relax_uniform_omega'
                  call star_relax_uniform_omega( &
                     id, new_omega, num_steps_to_relax_rotation, ierr)
                  if (failed('star_relax_uniform_omega')) return

               else if (near_zams_relax_omega_div_omega_crit) then
                  new_rotation_flag = .true.
                  call star_set_rotation_flag(id, new_rotation_flag, ierr)
                  if (failed('star_set_rotation_flag')) return
                  new_omega = new_omega_div_omega_crit*star_surface_omega_crit(id, ierr)
                  if (failed('star_surface_omega_crit')) return
                  write(*,2) 'new_omega_div_omega_crit', &
                     s% model_number, new_omega_div_omega_crit
                  call star_relax_uniform_omega( &
                     id, new_omega, num_steps_to_relax_rotation, ierr)
                  if (failed('star_relax_uniform_omega')) return

               end if     
                            
            end if
            
            if (s% rotation_flag) then      
                  
               if (s% model_number <= set_surf_rotation_v_step_limit) then
                  new_omega = new_surface_rotation_v*1d5/(s% photosphere_r*Rsun)
                  write(*,2) 'surface_rotation_v', s% model_number, new_surface_rotation_v
                  write(*,2) 'omega', s% model_number, new_omega
                  call star_set_uniform_omega(id, new_omega, ierr)
                  if (failed('star_set_uniform_omega')) return
                  
               else if (s% model_number <= set_omega_step_limit) then
                  write(*,2) 'omega', s% model_number, new_omega
                  if (failed('star_surface_omega_crit')) return
                  call star_set_uniform_omega(id, new_omega, ierr)
                  if (failed('star_set_uniform_omega')) return
                  
               else if (s% model_number <= set_omega_div_omega_crit_step_limit) then
                  new_omega = new_omega_div_omega_crit*star_surface_omega_crit(id, ierr)
                  write(*,2) 'omega_div_omega_crit', s% model_number, new_omega_div_omega_crit
                  write(*,2) 'omega', s% model_number, new_omega
                  if (failed('star_surface_omega_crit')) return
                  call star_set_uniform_omega(id, new_omega, ierr)
                  if (failed('star_set_uniform_omega')) return
               end if      
                        
            end if 
                        
            if (result == keep_going) then
               ! if you have data that needs to be saved and restored for restarts, 
               ! save it in s% extra_iwork and s% extra_work
               ! before calling star_finish_step
               if (pgstar_flag) call read_pgstar_controls(s, ierr) 
                  ! do this before call extras_finish_step
               if (failed('read_pgstar_controls')) return               
               result = extras_finish_step(s, id, id_extra)
               if (result == terminate) exit evolve_loop
               if (result /= keep_going) then
                  write(*,*) 'extras_finish_step must return either keep_going or terminate'
                  return
               end if
               
               !call test_get_history_data
               !call test_get_profile_data
               
               ! Get names and values for history
               num_log_columns = num_standard_history_columns(s) + how_many_extra_history_columns(s, id, id_extra)
               allocate(                            &
                    log_names(num_log_columns),     &
                    log_vals(num_log_columns),      &
                    log_is_int(num_log_columns),    &
                    stat=ierr)
               call get_data_for_history_columns(s, id_extra, &
                    how_many_extra_history_columns, data_for_extra_history_columns, &
                    num_log_columns,s% nz, log_names, log_vals, log_is_int, ierr)
               ! Get names and values from profiles
               num_profile_columns = num_standard_profile_columns(s) + how_many_extra_profile_columns(s, id, id_extra)
               allocate(                                     &
                    profile_names(num_profile_columns),      &
                    profile_vals(s% nz,num_profile_columns), &
                    profile_is_int(num_profile_columns),     &
                    stat=ierr)
               call get_data_for_profile_columns(s, id_extra, &
                    how_many_extra_profile_columns, data_for_extra_profile_columns, &
                    num_profile_columns,s% nz, profile_names, profile_vals, profile_is_int, ierr)
               
               result = se_finish_step(s, use_se_output, use_se_names,                &
                    se_codev, se_modname, se_prefix, se_num_mod_output,               &
                    how_many_extra_history_columns, data_for_extra_history_columns,   &
                    how_many_extra_profile_columns, data_for_extra_profile_columns,   &
                    log_names, log_vals, log_is_int, num_log_columns,                 &
                    profile_names, profile_vals, profile_is_int, num_profile_columns)

               deallocate(log_names,log_vals,log_is_int)
               deallocate(profile_names,profile_vals,profile_is_int)

               if (result == terminate) exit evolve_loop
               if (result /= keep_going) then
                  write(*,*) 'se_finish_step must return either keep_going or terminate'
                  return
               end if
               result = star_finish_step(id, id_extra, .false., &
                  how_many_extra_profile_columns, data_for_extra_profile_columns, &
                  how_many_extra_history_columns, data_for_extra_history_columns, ierr)
               if (failed('star_finish_step')) return
               if (result /= keep_going) exit evolve_loop               
               if (pgstar_flag) call update_pgstar_plots(s,ierr)
               if (failed('update_pgstar_plots')) return
            else if (result == terminate) then
               if (result_reason == result_reason_normal) then
                  !write(*, '(a, i12)') 'terminate: save profile for model number ', s% model_number
                  call save_profile(id, id_extra, &
                     how_many_extra_profile_columns, data_for_extra_profile_columns, &
                     3, ierr)
                  s% need_to_save_profiles_now = .false.
                  s% need_to_update_history_now = .true.
                  result = star_finish_step(id, id_extra, save_photo_when_terminate, &
                     how_many_extra_profile_columns, data_for_extra_profile_columns, &
                     how_many_extra_history_columns, data_for_extra_history_columns, ierr)
                  if (failed('star_finish_step')) return
                  if (save_model_when_terminate) save_model_number = s% model_number                  
                  if (save_pulsation_info_when_terminate) &
                     save_pulsation_info_for_model_number = s% model_number
                  call do_saves( &
                     id, id_extra, s, &
                     how_many_extra_history_columns, &
                     data_for_extra_history_columns, &
                     how_many_extra_profile_columns, &
                     data_for_extra_profile_columns)
               end if
               exit evolve_loop
            end if
            
            call do_saves( &
               id, id_extra, s, &
               how_many_extra_history_columns, &
               data_for_extra_history_columns, &
               how_many_extra_profile_columns, &
               data_for_extra_profile_columns)

            if (s% doing_timing) then
               call system_clock(time1_extra,clock_rate)
               after_step_timing = after_step_timing + dble(time1_extra - time0_extra) / clock_rate
               check_time_end = eval_total_times(s% id, ierr)
               check_after_step_timing = check_after_step_timing + (check_time_end - check_time_start)
            end if
            
         end do evolve_loop

         if (s% doing_timing) call show_times(id,s,time0)
         
         result_reason = get_result_reason(id, ierr)
         if (result_reason /= result_reason_normal) then
            write(*, *) 
            write(*, '(a)') 'terminated evolution: ' // trim(result_reason_str(result_reason))
            write(*, *)
         end if
         
         if (s% termination_code > 0 .and. s% termination_code <= num_termination_codes) then
            write(*, '(a)') 'termination code: ' // trim(termination_code_str(s% termination_code))
         end if
         
         if (pause_before_terminate) then
            write(*,'(a)') 'pause_before_terminate: hit RETURN to continue'
            read(*,*)
         end if

         call extras_after_evolve(s, id, id_extra, ierr)
         if (failed('after_evolve_extras')) return

         call se_after_evolve(s, ierr)
         if (failed('se_after_evolve')) return

         if (pgstar_flag) call update_pgstar_plots(s,ierr)
         if (failed('update_pgstar_plots')) return

         call star_after_evolve(id, ierr)
         if (failed('star_after_evolve')) return

         call write_terminal_summary(id, ierr)
         if (failed('write_terminal_summary')) return
         
         if (do_free_star) then
            call free_star(id, ierr)
            if (failed('free_star')) return
         end if
         
         
         contains

            
         subroutine relax_tau_factor
            real(dp) :: next
            include 'formats.inc'
            write(*,*) 'relax_to_this_tau_factor < s% tau_factor', relax_to_this_tau_factor < s% tau_factor
            write(*,1) 'relax_to_this_tau_factor', relax_to_this_tau_factor
            write(*,1) 's% tau_factor', s% tau_factor
            if (relax_to_this_tau_factor < s% tau_factor) then
               next = 10**(safe_log10(s% tau_factor) - dlogtau_factor)
               if (next < relax_to_this_tau_factor) next = relax_to_this_tau_factor
            else
               next = 10**(safe_log10(s% tau_factor) + dlogtau_factor)
               if (next > relax_to_this_tau_factor) next = relax_to_this_tau_factor
            end if
            s% tau_factor = next
            write(*,1) 'relax_tau_factor', next, relax_to_this_tau_factor
         end subroutine relax_tau_factor
      
      
         subroutine check_num_units(n)
            integer, intent(in) :: n
            integer :: current_num_iounits_in_use
            integer :: i
            include 'formats.inc'
            current_num_iounits_in_use = number_iounits_allocated()
            if (prev_num_iounits_in_use >= 3 .and. &
               prev_num_iounits_in_use < current_num_iounits_in_use) then
               do i=1,10
                  write(*,*)
               end do
               write(*,2) 'n', n
               write(*,2) 'prev_num_iounits_in_use', prev_num_iounits_in_use
               write(*,2) 'current_num_iounits_in_use', current_num_iounits_in_use
               stop 'check_num_units' 
            
            end if
            prev_num_iounits_in_use = current_num_iounits_in_use      
         end subroutine check_num_units
         
         
         logical function stop_is_requested()
            use utils_lib
            integer :: iounit, ierr
            stop_is_requested = .false.
            if (len_trim(stop_if_this_file_exists) == 0) return
            ierr = 0
            iounit = alloc_iounit(ierr); if (ierr /= 0) return
            open(unit=iounit, file=trim(stop_if_this_file_exists), &
               status='old', action='read', iostat=ierr)
            call free_iounit(iounit)
            if (ierr /= 0) return
            close(iounit)
            write(*,*) 'stopping because found file ' // trim(stop_if_this_file_exists)
            stop_is_requested = .true.
         end function stop_is_requested


         subroutine test_get_history_data
            character (len=maxlen_profile_column_name), pointer :: col_names(:) ! (num_columns)
            real(dp), pointer :: col_vals(:) ! (num_columns)
            logical, pointer :: is_int(:) ! (num_columns) true iff the values in the column are integers
            integer :: num_history_columns, ierr, i
            include 'formats.inc'
            num_history_columns = &
               num_standard_history_columns(s) + how_many_extra_history_columns(s, id, id_extra)
            write(*,2) 'num_history_columns', num_history_columns
            ierr = 0
            allocate(col_names(num_history_columns), col_vals(num_history_columns), &
               is_int(num_history_columns), stat=ierr)
            if (ierr /= 0) then
               write(*,*) 'allocate failed for test_get_history_data'
               stop 1
            end if
            call get_data_for_history_columns(s, id_extra, &
               how_many_extra_history_columns, data_for_extra_history_columns, &
               num_history_columns, s% nz, col_names, col_vals, is_int, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in get_data_for_history_columns'
               stop 1
            end if
            do i = 1, num_history_columns
               if (is_int(i)) then
                  write(*,3) 'INT ' // trim(col_names(i)), i, int(col_vals(i))
               else
                  write(*,2) trim(col_names(i)), i, col_vals(i)
               end if
            end do
            deallocate(col_names, col_vals, is_int)
            write(*,*) 'done test_get_history_data'
            !stop 
         end subroutine test_get_history_data


         subroutine test_get_profile_data
            character (len=maxlen_profile_column_name), pointer :: col_names(:) ! (num_columns)
            real(dp), pointer :: col_vals(:,:) ! (nz,num_columns)
            logical, pointer :: is_int(:) ! (num_columns) true iff the values in the column are integers
            integer :: num_profile_columns, ierr, i
            include 'formats.inc'
            num_profile_columns = &
               num_standard_profile_columns(s) + how_many_extra_profile_columns(s, id, id_extra)
            write(*,2) 'num_profile_columns', num_profile_columns
            ierr = 0
            allocate(col_names(num_profile_columns), is_int(num_profile_columns), &
               col_vals(s% nz,num_profile_columns), stat=ierr)
            if (ierr /= 0) then
               write(*,*) 'allocate failed for test_get_profile_data'
               stop 1
            end if
            call get_data_for_profile_columns(s, id_extra, &
               how_many_extra_profile_columns, data_for_extra_profile_columns, &
               num_profile_columns, s% nz, col_names, col_vals, is_int, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in get_data_for_profile_columns'
               stop 1
            end if
            do i = 1, num_profile_columns
               if (is_int(i)) then
                  write(*,3) 'INT ' // trim(col_names(i)), i, int(col_vals(1,i))
               else
                  write(*,2) trim(col_names(i)), i, col_vals(1,i)
               end if
            end do
            deallocate(col_names, col_vals, is_int)
            write(*,*) 'done test_get_profile_data'
            !stop 
         end subroutine test_get_profile_data
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed


      end subroutine run1_star     
      
      
      subroutine show_times(id, s, time0)
         use utils_lib, only: utils_OMP_GET_MAX_THREADS
         
         integer, intent(in) :: id, time0
         type (star_info), pointer :: s
         integer :: time1, clock_rate, ierr, pass, omp_num_threads
         real(dp) :: elapsed_time, total, misc
         include 'formats.inc'
         ierr = 0
         omp_num_threads = utils_OMP_GET_MAX_THREADS()
         elapsed_time = before_step_timing + step_loop_timing + after_step_timing
         s% time_total = check_before_step_timing + check_step_loop_timing + check_after_step_timing

         if (.true.) then
            write(*,*)
            write(*,*) 'here are the items'
            total = 0
            misc = 0
            call do1_time('run_star', elapsed_time - s% time_total)
            
            call do1_time('time_do_mesh', s% time_do_mesh)
            call do1_time('time_do_mesh_plan', s% time_do_mesh_plan)
            call do1_time('time_do_mesh_adjust', s% time_do_mesh_adjust)
            call do1_time('time_newton_xscale', s% time_newton_xscale)
            call do1_time('time_newton_eval_eqn', s% time_newton_eval_eqn)
            call do1_time('time_newton_size_equ', s% time_newton_size_equ)
            call do1_time('time_newton_size_B', s% time_newton_size_B)
            call do1_time('time_newton_enter_setmatrix', s% time_newton_enter_setmatrix)
            call do1_time('time_do_adjust_mass', s% time_do_adjust_mass)
            call do1_time('time_do_report', s% time_do_report)
            call do1_time('time_next_timestep', s% time_next_timestep)
            call do1_time('time_write_profile', s% time_write_profile)
            call do1_time('time_write_log', s% time_write_log)
            call do1_time('time_write_photo', s% time_write_photo)
            call do1_time('time_pgstar', s% time_pgstar)
            call do1_time('time_set_basic_vars', s% time_set_basic_vars)
            call do1_time('time_set_rotation_vars', s% time_set_rotation_vars)
            call do1_time('time_set_mlt_vars', s% time_set_mlt_vars)
            call do1_time('time_set_eqn_vars', s% time_set_eqn_vars)
            call do1_time('time_eval_eqns', s% time_eval_eqns)
            call do1_time('time_set_newton_vars', s% time_set_newton_vars)
            call do1_time('time_newton_mtx', s% time_newton_mtx)
            call do1_time('time_newton_self', s% time_newton_self)
            call do1_time('time_newton_test', s% time_newton_test)
            call do1_time('time_solve_burn_non_net', s% time_solve_burn_non_net)
            call do1_time('time_solve_burn_in_net', s% time_solve_burn_in_net)
            call do1_time('time_solve_mix', s% time_solve_mix)
            call do1_time('time_op_split_control', s% time_op_split_control)
            call do1_time('time_check_model', s% time_check_model)
            call do1_time('time_prep_new_step', s% time_prep_new_step)
            call do1_time('time_prep_new_try', s% time_prep_new_try)
            call do1_time('time_prep_for_retry', s% time_prep_for_retry)
            call do1_time('time_do_winds', s% time_do_winds)
            call do1_time('time_save_for_d_dt', s% time_save_for_d_dt)
            call do1_time('time_diffusion', s% time_diffusion)
            call do1_time('time_evolve_set_vars', s% time_evolve_set_vars)
            call do1_time('time_save_start', s% time_save_start)
            call do1_time('time_check_newly_non_conv', s% time_check_newly_non_conv)
            call do1_time('time_struct_burn_mix', s% time_struct_burn_mix)
            
            call do1_time('misc', misc)
            call do1_time('sum', total+misc)
            call do1_time('time_total', s% time_total)

            write(*,*)
         end if
         
         write(*,*)
         write(*,2) 'nz', s% nz
         write(*,2) 'nvar', s% nvar
         write(*,2) 'species', s% species
         write(*,2) 'reactions', s% num_reactions
         write(*,*)
         if (.false.) then
            write(*,1) 'check_before_step_timing', check_before_step_timing
            write(*,1) 'check_step_loop_timing', check_step_loop_timing
            write(*,1) 'check_after_step_timing', check_after_step_timing
            write(*,*)
            write(*,1) 'before_step_timing', before_step_timing
            write(*,1) 'step_loop_timing', step_loop_timing
            write(*,1) 'after_step_timing', after_step_timing
            write(*,*)
         end if

         if (s% time_op_split_control == 0) then ! normal case
            do pass = 1, 2
               if (pass == 1) then
                  write(*,'(a8)',advance='no') 'threads'
               else
                  write(*,'(i5,3x)',advance='no') omp_num_threads
               end if
               call show1(pass, 'total', elapsed_time, total)
               total = 0
               call show1(pass, 'net', s% time_net, total)
               call show1(pass, 'eos', s% time_eos, total)
               call show1(pass, 'neu+kap', s% time_neu + s% time_kap, total)
               call show1(pass, 'mlt', s% time_set_mlt_vars, total)
               call show1(pass, 'eqns', s% time_eval_eqns, total)
               call show1(pass, 'mesh', s% time_do_mesh + s% time_do_mesh_plan + s% time_do_mesh_adjust, total)
               call show1(pass, 'matrix', s% time_newton_mtx, total)
               call show1(pass, 'newton', s% time_newton_self, total)
               call show1(pass, 'profile', s% time_write_profile, total)
               call show1(pass, 'history', s% time_write_log, total)
               call show1(pass, 'photo', s% time_write_photo, total)
               call show1(pass, 'other', elapsed_time - total, total)
               write(*,*)
            end do
            write(*,'(99a8)') 'threads', 'msh', 'mshpln', 'mshadj'
            write(*,'(i5,3x,99f8.4)') omp_num_threads, &
               s% time_do_mesh, s% time_do_mesh_plan, s% time_do_mesh_adjust
            write(*,*)
            write(*,'(a,i8)') 'total_num_jacobians', s% total_num_jacobians
            write(*,*)
               
               
         else
            write(*,'(99a10)') &
               'threads', 'total', 'net', 'burn_net', 'burn_slv', 'mix', 'eos', 'kap', 'eqns', &
               'matrix', 'callbck', 'newton', 'mesh', 'profile', 'history', 'photo', 'other'
            write(*,'(i5,3x,99f10.4)') omp_num_threads, elapsed_time, s% time_net, &
               s% time_solve_burn_in_net, s% time_solve_burn_non_net, s% time_solve_mix, &
               s% time_eos, s% time_kap, s% time_eval_eqns, &
               s% time_newton_mtx, &
               s% time_newton_xscale + s% time_newton_eval_eqn + s% time_newton_size_equ + &
                  s% time_newton_size_B + s% time_newton_enter_setmatrix, &
               s% time_newton_self, &
               s% time_do_mesh + s% time_do_mesh_plan + s% time_do_mesh_adjust, &
               s% time_write_profile, s% time_write_log, s% time_write_photo, elapsed_time - s% time_total
         end if
         write(*,*)
         
         
         contains
         
         
         subroutine show1(pass, name, value, total)
            integer, intent(in) :: pass
            character (len=*), intent(in) :: name
            real(dp), intent(in) :: value
            real(dp), intent(inout) :: total
            if (pass == 1) then
               write(*,'(a8)',advance='no') trim(name)
            else
               write(*,'(f8.4)',advance='no') value
               total = total + value
            end if
         end subroutine show1
     
         subroutine do1_time(str, t)
            character (len=*), intent(in) :: str
            real(dp), intent(in) :: t
            if (.false. .and. t < 1d-3*elapsed_time) then
               misc = misc + t
               return
            end if
            write(*,'(a30,f16.4,f18.6)') str, t!, t/elapsed_time
            total = total + t
         end subroutine do1_time
         

      end subroutine show_times

      
      
      subroutine do_saves( &
            id, id_extra, s, &
            how_many_extra_history_columns, &
            data_for_extra_history_columns, &
            how_many_extra_profile_columns, &
            data_for_extra_profile_columns)
         integer, intent(in) :: id, id_extra
         type (star_info), pointer :: s
         interface

            integer function how_many_extra_history_columns(s, id, id_extra)
               use const_def, only: dp
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
            end function how_many_extra_history_columns
            
            subroutine data_for_extra_history_columns(s, id, id_extra, n, names, vals, ierr)
               use const_def, only: dp
               use star_def, only: star_info, maxlen_history_column_name
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra, n
               character (len=maxlen_history_column_name) :: names(n)
               real(dp) :: vals(n)
               integer, intent(out) :: ierr
            end subroutine data_for_extra_history_columns
      
            integer function how_many_extra_profile_columns(s, id, id_extra)
               use star_def, only: star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra
            end function how_many_extra_profile_columns      
      
            subroutine data_for_extra_profile_columns(s, id, id_extra, n, nz, names, vals, ierr)
               use const_def, only: dp
               use star_def, only: star_info, maxlen_profile_column_name
               type (star_info), pointer :: s
               integer, intent(in) :: id, id_extra, n, nz
               character (len=maxlen_profile_column_name) :: names(n)
               real(dp) :: vals(nz,n)
               integer, intent(out) :: ierr
            end subroutine data_for_extra_profile_columns

         end interface

         integer :: ierr
         ierr = 0
      
         if (s% model_number == save_model_number) then
            call star_write_model(id, save_model_filename, save_prev_along_with_current, ierr)
            if (failed('star_write_model')) return
            write(*, *) 'saved to ' // trim(save_model_filename)
         end if
         
         if (s% model_number == save_pulsation_info_for_model_number) then
            call star_write_pulsation_info(id, &
               s% add_center_point_to_pulse_info, &
               s% keep_surface_point_for_pulse_info, &
               s% add_atmosphere_to_pulse_info, &
               s% pulse_info_format, &
               save_pulsation_info_filename, ierr)
            if (failed('star_write_pulsation_info')) return
            write(*, *) 'pulsation data saved to ' // trim(save_pulsation_info_filename)
         end if
         
         if (s% model_number == save_short_format_for_model_number) then
            call star_write_short_format(id, save_short_format_filename, ierr)
            if (failed('star_write_short_format')) return
            write(*, *) 'short_format saved to ' // trim(save_short_format_filename)
         end if
         
         if (s% model_number == profile_model_number) then
            write(*, '(a, i7)') 'save profile for model number', s% model_number
            call save_profile(id, id_extra, &
               how_many_extra_profile_columns, data_for_extra_profile_columns, &
               3, ierr)
            if (failed('save_profile')) return
         end if
         
         if (internals_num >= 0) then
            write(*, '(a, i7)') 'write internals for model number', s% model_number
            call std_write_internals(id, internals_num)
            stop 'finished std_write_internals'
         end if
         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed
         
      end subroutine do_saves


      end module run_star_support
      
      
      
      
