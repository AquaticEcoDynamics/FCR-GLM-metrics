ptf %

&aed_models
   models =
   'aed_sedflux',
   'aed_tracer',
   'aed_noncohesive',
   'aed_oxygen',
   'aed_carbon',
   'aed_silica',
   'aed_nitrogen',
   'aed_phosphorus',
   'aed_organic_matter',
   'aed_phytoplankton',
   'aed_totals'
/
&aed_sedflux
   sedflux_model = 'Constant2D'
/
&aed_sed_const2d
   n_zones = 4
   active_zones = 1, 2, 3, 4
   Fsed_oxy = %        fsed_oxy1%, %        fsed_oxy2%, %        fsed_oxy3%, %        fsed_oxy4%
   !Fsed_oxy = -29.29916397869219, -21.66869075617110, -16.39155735828627, -10.90665072320945
   !Fsed_oxy = -32.76435121, -5.203754493, -0.181529858
   Fsed_ch4 = 15, 15, 15, 232.6331575
   Fsed_amm =  2.848790019, 1.0, 1.0, 1.0 !2.848790019, 1.7, 0.5, 1.59739687  !1,0.1,1
   !Fsed_amm =  2.848790019, 2.848790019, 2.848790019, 1.59739687!1,0.1,1
   Fsed_nit = -0.01,-0.001,-0.001,-0.0001
   !Fsed_nit = -0.13, -0.13, -0.13, -0.001
   !Fsed_frp = 0.009,0.009, 0.04, 0.25
   Fsed_frp = 0.0048, 0.0048, 0.001, 0.25
/
!&aed_sed_const2d
!   n_zones = 3
!   active_zones = 1, 2, 3
!   Fsed_oxy = -28, -7, -1
!   !Fsed_oxy = -32.76435121, -5.203754493, -0.181529858
!   Fsed_ch4 = 15, 15, 232.6331575
!   Fsed_amm = 2.848790019, 2.848790019, 1.59739687!1,0.1,1
!   Fsed_nit = -0.13, -0.13, -0.001
!   Fsed_frp = 0.0048, 0.001, 0.25
!/
&aed_tracer
   retention_time = .true.
   num_tracers = 1
/
&aed_noncohesive
   num_ss = 1
   ss_initial = 1, 1
   Ke_ss = 0.006, 0.063
   settling = 1
   w_ss = -0.001, -0.001
   d_ss = 2e-06, 1e-05
   rho_ss = 1500, 1800
   resuspension = 0
   simSedimentMass = .true.
   fs = 0.4, 0.4
   sed_porosity = 0.6
/
&aed_oxygen
   oxy_initial = 225
   Fsed_oxy = -10
   Ksed_oxy = %         ksed_oxy%
   theta_sed_oxy = %        thsed_oxy%
   Fsed_oxy_variable = 'SDF_Fsed_oxy'
   oxy_min = 0
   oxy_max = 500
/
&aed_carbon
   dic_initial = 91
   Fsed_dic = 3!0.131348518!0.001
   Ksed_dic = 98.87064869!100
   theta_sed_dic = 1.197284209!1.199999
   pH_initial = 6.2
   atm_co2 = 0.00041
   co2_model = 1
   alk_mode = 0
   ionic = 0.1
   co2_piston_model = 1
   ch4_initial = 5
   Rch4ox = 0.2
   Kch4ox = 0.2
   vTch4ox = 1.2
   Ksed_ch4 = 3.437
   theta_sed_ch4 = 1.2
   methane_reactant_variable = 'OXY_oxy'
   atm_ch4 = 1.76e-06
   ch4_piston_model = 1
   Fsed_ch4_variable = 'SDF_Fsed_ch4'
/
&aed_silica
   rsi_initial = 100
   Fsed_rsi = 1.001876
   Ksed_rsi = 1.002457
   theta_sed_rsi = 1.037413
   silica_reactant_variable = 'OXY_oxy'
/
&aed_nitrogen
   amm_initial = 2.6
   nit_initial = 0.1
   n2o_initial = 0.1
   Rnitrif = 0.01358676
   Knitrif = 62.02209
   theta_nitrif = 1.08
   nitrif_reactant_variable = 'OXY_oxy'
   nitrif_ph_variable = ''
   simNitrfpH = .false.
   Rnh4o2 = 1
   Rno2o2 = 1
   simN2O = 0
   Rn2o = 0.05
   Kpart_ammox = 1
   Kin_deamm = 1
   atm_n2o = 3.2e-07
   n2o_piston_model = 4
   Rnh4no2 = 0.001
   kanammox = 0.001
   Kanmx_nit = 1.320344
   Kanmx_amm = 0.8666655
   Rdenit = 9.968717
   Kdenit = 29.86566
   theta_denit = 1.062862
   Rdnra = 0.01123021
   Kdnra_oxy = 0.360534
   Ksed_amm = 41.25!30!74.02427825!30
   Ksed_nit = 73.26015
   Fsed_n2o = 0
   Ksed_n2o = 100
   theta_sed_amm = 1.068994!1.08
   theta_sed_nit = 1.068994
   Fsed_amm = 1
   Fsed_nit = -0.05
   Fsed_amm_variable = 'SDF_Fsed_amm'
   Fsed_nit_variable = 'SDF_Fsed_nit'
/
&aed_phosphorus
   frp_initial = 0.05
   Ksed_frp = 6.907107
   theta_sed_frp = 1.06609
   phosphorus_reactant_variable = 'OXY_oxy'
   Fsed_frp_variable = 'SDF_Fsed_frp'
   simPO4Adsorption = .false.
   ads_use_external_tss = .false.
   po4sorption_target_variable = 'NCS_ss1'
   PO4AdsorptionModel = 1
   Kpo4p = 0.1
   ads_use_pH = .false.
   Kadsratio = 1
   Qmax = 1
   w_po4ads = -9999
/
&aed_organic_matter
   poc_initial = 15
   doc_initial = 15
   pon_initial = 2
   don_initial = 1.1
   pop_initial = 0.1
   dop_initial = 0.01
   docr_initial = 150
   donr_initial = 9
   dopr_initial = 0.15
   cpom_initial = 0
   Rdom_minerl = 0.013
   Rpoc_hydrol = 0.001
   Rpon_hydrol = 0.001
   Rpop_hydrol = 1e-04
   theta_hydrol = 1.07
   theta_minerl = 1.07
   Kpom_hydrol = 33.66593
   Kdom_minerl = 22.36079
   simDenitrification = 1
   dom_miner_oxy_reactant_var = 'OXY_oxy'
   doc_miner_product_variable = 'CAR_dic'
   don_miner_product_variable = 'NIT_amm'
   dop_miner_product_variable = 'PHS_frp'
   dom_miner_nit_reactant_var = 'NIT_nit'
   f_an = 0.2063885
   K_nit = 10
   simRPools = .true.
   Rdomr_minerl = 0.000102192
   Rcpom_bdown = 0.05350772
   X_cpom_n = 0.005
   X_cpom_p = 0.001
   KeDOM = 0.005
   KePOM = 0.00096
   KeDOMR = 0.1
   KeCPOM = 0.00096
   simphotolysis = .false.
   photo_c = 0.75
   settling = 1
   w_pom = -0.01
   d_pom = 1e-05
   rho_pom = 1200
   w_cpom = -0.01
   d_cpom = 1e-05
   rho_cpom = 1400
   resuspension = 0
   resus_link = ''
   sedimentOMfrac = 2e-04
   Xsc = 0.5
   Xsn = 0.05
   Xsp = 0.005
   Fsed_doc = 0.1!1.408053
   Fsed_don = 0
   Fsed_dop = 0
   Ksed_dom = 93.12891
   theta_sed_dom = 1.057064
   diag_level = 10
/
&aed_phytoplankton
   num_phytos = 3
   the_phytos = 1, 2, 3
   settling = 1, 1, 1
   do_mpb = 0
   p_excretion_target_variable = 'OGM_dop'
   n_excretion_target_variable = 'OGM_don'
   c_excretion_target_variable = 'OGM_doc'
   si_excretion_target_variable = ''
   p_mortality_target_variable = 'OGM_pop'
   n_mortality_target_variable = 'OGM_pon'
   c_mortality_target_variable = 'OGM_poc'
   si_mortality_target_variable = ''
   p1_uptake_target_variable = 'PHS_frp'
   n1_uptake_target_variable = 'NIT_amm'
   n2_uptake_target_variable = 'NIT_nit'
   si_uptake_target_variable = 'SIL_rsi'
   do_uptake_target_variable = 'OXY_oxy'
   c_uptake_target_variable = 'CAR_dic'
   !dbase = 'aed/aed2_phyto_pars_28April2022_constant.nml'
   !dbase = 'aed/aed2_phyto_pars_20April2022_constant.nml'
   dbase = 'aed/aed2_phyto_pars_2May2022_RQT.nml'
   diag_level = 10
   min_rho = 900
   max_rho = 1200
/
&aed_totals
   outputLight = .true.
   TN_vars = 'NIT_nit','NIT_amm','OGM_don','OGM_donr','OGM_pon','PHY_in'
   TN_varscale = 1, 1, 1, 1, 1
   TP_vars = 'PHS_frp','OGM_dop','OGM_dopr','OGM_pop','PHY_ip'
   TP_varscale = 1, 1, 1, 1
   TOC_vars = 'OGM_doc','OGM_docr','OGM_poc','OGM_cpom','PHY_tphy'
   TOC_varscale = 1, 1, 1, 1
/
&aed_sediment
   sediment_model = 'DYNAMIC'
   mpb_link_variable = ''
   mag_link_variable = ''
   root_link_variable = ''
/
&aed_sed_candi
   spinup_days = 90
   spinup_dt = 0.25
   driver_dt = 900
   substep = 8
   n_zones = 1
   active_zones = 0
   zone_types = 1
   dbase = 'aed/aed_candi_params.csv'
   vars_files = 'aed/aed_sdg_vars.csv'
   geochem_file = 'aed/aed_geochem_pars.dat'
   swibc_mode = 0
   deepbc_mode = 1
   swibc_file = 'aed/aed_sediment_swibc.dat'
   deepbc_file = 'aed/aed_sediment_deepbc.dat'
   swibc_filevars = ''
   deepbc_filevars = ''
   flux_scale = 1
   SolidInitialUnit = 'mmolLsolid'
   OMInitMethodL = 'LI_I'
   OM_topL = 1
   OM_minL = 0.9
   OM_cfL = 0.6
   InitMinDepthL = 99
   OMInitMethodR = 'LI_I'
   OM_topR = 1
   OM_minR = 0.9
   OM_cfR = 0.6
   InitMinDepthR = 99
   POMVR = 0.3
   diag_level = 10
   output_profiles = NA
   morevariables = 'Rgpp','Rrsp','FO2'
   output_diag_vars = 'oxy','amm','docl','docr'
   n_ddpths = 1
   output_diag_depths = 1
/
&aed_zooplankton
   num_zoops = 1
   the_zoops = 1
   dn_target_variable = 'OGM_don'
   pn_target_variable = 'OGM_pon'
   dp_target_variable = 'OGM_dop'
   pp_target_variable = 'OGM_pop'
   dc_target_variable = 'OGM_doc'
   pc_target_variable = 'OGM_poc'
   dbase = 'aed/aed_zoop_pars.nml'
/
&aed_bivalve
   num_biv = 1
   the_biv = 1
   X_c = 1
   dbase = 'aed/aed_bivalve_pars.nml'
   n_zones = 1
   active_zones = 1
   initFromDensity = .false.
   simBivTracer = .true.
   bt_renewal = 10
   simStaticBiomass = .false.
   simBivFeedback = .true.
   pn_target_variable = 'OGM_pon'
   dp_target_variable = 'PHS_frp'
   pp_target_variable = 'OGM_pop'
   dc_target_variable = 'OGM_doc'
   pc_target_variable = 'OGM_poc'
   do_uptake_variable = 'OXY_oxy'
   ss_uptake_variable = 'NCS_ss1'
   simFixedEnv = .false.
   fixed_temp = 20
   fixed_oxy = 400
   fixed_food = 20
   diag_level = 10
/
