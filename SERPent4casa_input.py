##############################################################################
# SERPent4casa Input File                                                    #
#                                                                            #
##############################################################################
# Original SERPent version 13/09/2011 - Luke Peck                            #
#                                                                            #
# Last update 11/2017 - Danielle Fenech                                      #
#                                                                            #
##############################################################################
# This file contains all the information the user must input before running  #
# SERPent4casa.                                                              #
#                                                                            #
# To run SERPent, type in terminal: casa -c SERPent4casa.py                  #
##############################################################################


#############################################
#### Input and output data and directory ####
#### information                         ####
#############################################


Name = '3C2862014eMERLIN_L_split.ms'                                  #Filename only

path2data = '/import/cobrasscratch1/dmf/python/Anita_CASA_SERPent/'   # Folder path which contains the ms 
                                                                      # file listed


#### Output files directory Information ####

path2folder = '/import/cobrasscratch1/dmf/python/Anita_CASA_SERPent/filestore/'  
                                                                      # Directory path to where the output 
                                                                      # text files (FG text files and
                                                                      # intermediate pickle files) will be 
                                                                      # written to.




#### Parallelisation information ####
NCPU = 16           # Define here the number of CPUs you want to use.
                    # Parallelization is currently implemented on a baseline basis
                    # per source therefore the maximum number of CPUs ultilized 
                    # will be the number of baselines in the data. e.g. for e-MERLIN 
                    # a maximum of 21 CPUs will be used by SERPent (21 baselines) for 
                    # one source.
                    # If you wish to use all available CPUs, set NCPU = 0




#########################################################
#### Phase Cal. and Lovell_dropout settings          ####
#### e-MERLIN specific flagging to deal with Lovell  ####
#### stationary scans when the other telescopes slew ####
#### to observe the phase cal.                       ####
#### If this is not required set phasecal = 'no'     ####
#########################################################

do_Lovell_scan = 'yes'            # If you wish to check for Lovell stationary scans on MERLIN
                                  # or e-MERLIN data use 'yes', else 'no'. If do_Lovell_scan = 'no'
                                  # is used the phascal list will be ignored.

phasecal = ['2007+404']     # If any of the sources (multi files) are phase calibrators
                                  # list them here for checking wth the Lovell stationary scan 
                                  # passage. N.B. for MERLIN/e-MERLIN data only.
                                  # This information is used for the Lovell Stationary Scans passage of
                                  # SERPent if the source is the phasecal designated here, the telescope
                                  # is e-MERLIN and the baseline contains the Lovell telescope.
do_lovell_cross = 'yes'            # If do_lovell_cross = 'yes' the Lovell dropout check will be performed
                                  # on all stokes parameters. If do_lovell_cross = 'no' the Lovell dropout
                                  # check will only be performed on parallel hands (not cross-hands). This
                                  # should be sufficient to catch genuine dropouts.



#########################################################
#### Zero level dropout settings                     ####
#### Executes the Zero Level code to flag telescope  ####
#### dropouts and other system failures where the    ####
#### visibility amplitudes drop off.                 ####
#### If this is not required set zero_level = 'no'   ####
#########################################################


                          # system failures where the visibilities drop to around 0 Jy in the same 
                          # scans where good data is present.
zero_level = 'yes'        # If you want to execute this code set this variable to 'yes'
                          # else set it to 'no'
                          # This setting executes the Zero Level code to flag telescope dropouts and 
                          # other system failures where the visibilities drop to around 0 Jy in the 
                          # same scans where good data is present.

zero_flag_allif = 'yes'   # The zero_level_dropout code operates on each IF individually. If you 
                          # would like to apply the flags from each IF to all IFs, set this variable
                          # to 'yes'. If you would prefer to only apply the flags to the IF in which 
                          # they're found set this variable to 'no'. N.B. 'yes' is the default setting.



#######################################################
#### Settings for merging polarisation information ####
#### SERPent is executed on each polarisation      ####
#### independently. The settings below determine   ####
#### how the flagging information from each        ####
#### polarization is combined.                     ####
#######################################################


select_stokes = 'yes'              # If set to 'yes' only the polarizations listed in stokes below
                                   # will be processed else if set to 'no' all polarisations present
                                   # in the data will be processed.

stokes = ['RR','LL']               # Which polarizations should SERPent include. N.B. This will be 
                                   # ignored if select_stokes == 'no'.

flag_all_stokes = 'yes'            # If select_stokes ='yes' and a subset of polarizations is chosen
                                   # setting flag_all_stokes = 'yes' will apply the flags to ALL 
                                   # polarizations in the data. e.g. if the data contains RR,LL,LR,RL
                                   # and only RR,LL are selected in stokes for processing then setting 
                                   # flag_all_stokes = 'yes' will write flag entries for all four 
                                   # polarizations.
                                   # N.B. if select_stokes = 'yes' and flag_stokes = 'yes' then 
                                   # coadd_polarization_flags will be forced to 'yes'
                                   # i.e. all flags from all polarizations processed will be 
                                   # combined and applied to all polarizations present in the data.
                                   # If select_stokes = 'yes' and flag_all_stokes = 'no' then flag 
                                   # entries will only be written for the polarizations specified 
                                   # in the stokes setting above.

coadd_polarization_flags = 'no'    # Combine flags for all polarizations = 'yes', 'no' flags all
                                   # polarizations separately. If certain polarisations are 
                                   # selected with select_stokes = 'yes', setting 
                                   # coadd_polarizations_flags = 'yes' will combine the flags 
                                   # write entries for only those polarisations selected. To write
                                   # flag entries for ALL polarisation in the data set 
                                   # flag_all_stokes = 'yes'

coadd_zero_flags = 'no'            # If you choose not to combine flags for all polarizations 
                                   # coadd_zero_flags = 'no' will keep the zero dropout flags 
                                   # separate for each stokes, whereas coadd_zero_flags = 'yes'
                                   # will combine them. N.B. This setting will be ignored if 
                                   # zero_level = 'no'

coadd_lovell_flags = 'yes'         # If you choose not to combine flags for all polarizations
                                   # coadd_lovell_flags = 'no' will keep the lovell dropout flags 
                                   # separate for each stokes, whereas coadd_lovell_flags = 'yes'
                                   # will combine them, this should be the default. N.B. This
                                   # setting will be ignored if do_Lovell_scan = 'no'


###############################################
#### Source selection                      ####
#### Choose which sources to flag          ####
#### N.B. If flagging a true single-source ####
#### (i.e. SPLIT) file the source MUST be  ####
#### specified in the flag_list and        ####
#### flag_source = 'choose' used.          ####
###############################################



do_sumthreshold = 'no'             # Choose whether the full SumThreshold flagging sequence runs. 
                                   # Default is 'yes' and will almsot always need to be, but if 
                                   # set to 'no' it will not run. This is useful for running 
                                   # the zero_level_dropouts or the lovell_dropouts without 
                                   # performing full flagging.





flag_sources = 'choose'            # Variable to make source selection available. If 'all' flag all
                                   # sources in the data. If 'choose', the sources must be specified 
                                   # in the 'flag_list' variable below.
                                   # NB: If flagging a single-source file, you MUST set 
                                   # flag_sources='choose' and enter the source name in the flag_list.
flag_list = ['1331+305']
                                   # A list of comma separated sources to flag, will be ignored unless
                                   # flag_sources='choose' above.
                                   # NB: If flagging a single-source file, you MUST set 
                                   # flag_sources='choose' and enter the source name in the flag_list.



###############################################
#### SPW selection                         ####
#### Choose which spws to flag (equivalent ####
#### to IFs)                               ####
###############################################


spwlist = [0,1]         # If spwlist = 'all', all available spws for each source/baseline
                        # will be flagged. Otherwise to specify include a python list of
                        # required SPWs e.g. spwlist = [0,1] will incude the first two SPWs
                        # for flagging.




#####################################################
#### Basic flagging options                      ####
#### Choosing generic flagging parameters to be  ####
#### used on all selected baselines.             ####
#### Separate parameters can be specified for    ####
#### each source if required.                    ####
#### To use the default settings for all sources ####
#### set flagging_options = 'default'            ####
#####################################################


flagging_options = 'choose'         # variable to define whether to use the flagging options
                                    # below or the default options in the SERPent.py file
                                    # 'default' ignores whatever variables are set in this file
                                    # Options are: 'choose', 'source' or 'default'
                                    #
                                    # If 'choose' is specified, one set of options should
                                    # be listed below, which will be used for all sources
                                    #
                                    # If 'source' is specified, a value per source for each 
                                    # parameter should be specified below using a comma-
                                    # separated list.
                                    # If flag_sources = 'all' and flagging_options = 'source'
                                    # the values should be specified in the source order given
                                    # in the AIPS data source table.
                                    # If flag_sources = 'choose' and flagging_options = 'source'
                                    # the values should be specified in the order given in the 
                                    # flag_list option.

#### Flagging parameters  ####


aggressiveness_first_run = [6]      # How aggressive the first run of the flagger is
                                    # A lower number is more aggressive
                                    # N.B. default is 4
max_subset_first_run = [2]          # Maximum subset for first run of the flagger
                                    # This should be a binary number (1,2,4,8,16...)
                                    # N.B. default is 2
aggressiveness_second_run = [6]     # How aggressive the second run of the flagger is
                                    # N.B. default is 4
max_subset_second_run = [2]         # Maximum subset for second run of the flagger
                                    # N.B. default is 2
rho = [1.5]                         # Difference in coarseness between each threshold
                                    # level. 1.5 should be specified here.
                                    # N.B. default is 1.5
kickout_sigma_level = [3.5]         # The kickout clause is tested during before each subset
                                    # run through and takes the form:  
                                    # median + kickout_sigma_level * MAD
                                    # a lower coefficient = deeper flagging/ more aggressive.
                                    # N.B. default is 3.0

                                    # When specifying a value per source use a comma-
                                    # separated list e.g.
                                    # aggressiveness_first_run = [4,8,2]

#### Baseline selection for generic flagging options ####

which_baselines = 'all'             # Variable to define whether to flag all baselines
                                    # or a select few.
                                    # Options are: 'all' or 'choose'.

                                    # If you chose which_baselines = 'choose', add the baselines to 
                                    # the below list:
                                    # Note you'll obviously have to know the antenna numbers for 
                                    # any given baseline...
                                    # If 'all' is chosen then the baseline list below is ignored.

baselines = ['1-5']
                                    # chosen baselines for flagging.
                                    # Format: ['5-7', '7-8'], i.e. for baselines 5-7 and 7-8 
                                    # (strings!) order of baselines does not matter.


#######################################################
#### Advanced flagging options                     ####
#### Choosing baseline specific flagging options   ####
#### These will overide the generic parameter      ####
#### settings above including the flagging_options ####
#### and which_baselines settngs.                  ####
#### If these are not required set                 ####
#### dobline_specific = 'no'                       ####
#######################################################


dobline_specific = 'no'             # Choices are 'yes' or 'no'. 'Yes' will overide the flagging_options 
                                    # and which_baseline options above. At least one list of baselines 
                                    # and respective parameters must be listed if this setting is used. 

                                    # The following lists baseline arrays and corresponding flagging 
                                    # parameters. Multiple baselines can be listed to use the same 
                                    # parameter settings e.g. setting 'baselines_0 = ['2-7','2-9']' 
                                    # 'parameters_0 = [25,16,25,128,1.5,3.0]' will use these parameters 
                                    # for both baselines. Multiple instances can be used. So 
                                    # 'baselines_1 = ['1-5','5-9','6-7']' could be added with a  
                                    # corresponding 'parameters_1' list.
                                    
baselines_0 = ['1-5']               # List of baselines corresponding to each list of parameters
baselines_1 = ['2-5']

parameters_0 = [6,4,6,4,1.5,3.0]   # List of parameters corresponding to each list of baselines
parameters_1 = [6,4,6,4,1.5,3.0]    # parameters MUST be listed in the following order:
                                        # [aggressiveness_first_run, max_subset_first_run, 
                                        # aggressiveness_second_run, max_subset_second_run, rho, 
                                        # kickout_sigma_level]


###########################################################
#### Very advanced flagging options                    ####
#### Choosing baseline and source specific options     ####
#### These will overide the generic flagging           ####
#### parameter settings including the flagging_options ####
#### and which_baselines settngs above. This will also ####
#### overide the advanced baseline-specfic settings    ####
#### above.                                            ####
#### If these are not required set                     ####
#### dosource_specific = 'no'                          ####
###########################################################

dosource_specific = 'no'                 # Choices are 'yes' or 'no'. 'Yes' will overide the flagging_options 
                                          # and the dobline_specific options above. At least one list of source 
                                          # and baseline parameters must be listed if this setting is used. 


source_0_baselines_0 = ['2-5','2-6','2-9','5-6','5-7','5-8','5-9','6-7','6-8','6-9']            
                                             # List of source and baseline combinations corresponding to each 
                                             # list of parameters. The source order will be either the same 
source_0_baselines_1 = ['2-7','7-9','8-9']   # order as the AIPS source table (if flag_sources = 'all' is used) 
                                             # or it will be the order specified in the flag_list (if 
                                             # flag_sources = 'choose' is selected).
                                             # EXAMPLE: if flag_list = ['1407+284','1331+305'] then 
                                             # source_0_baseline_0 should specify baselines for '1407+284'
                                             # with corresponding parameters in source_0_parameters_0
                                             # Further baselines and parameter settings for source '1407+284'
                                             # can be specified in source_0_baselines_n and source_0_parameters_n
                                             # Baseline selections and corresponding parameters for '1331+305' 
                                             # would be specified using source_1_baselines_n and 
                                             # source_1_parameters_n.
                                             # EXAMPLE: if flag_list = ['1407+284','1331+305'] then 
                                             # source_0_baseline_0 should specify baselines for '1407+284'
                                             # with corresponding parameters in source_0_parameters_0
                                             # Further baselines and parameter settings for source '1407+284'
                                             # can be specified in source_0_baselines_n and source_0_parameters_n
                                             # Baseline selections and corresponding parameters for '1331+305' 
                                             # would be specified using source_1_baselines_n and 
                                             # source_1_parameters_n.
source_0_baselines_2 = ['2-8','7-8'] 
source_1_baselines_0 = ['2-5','2-6','2-9','5-6','5-7','5-8','5-9','6-7','6-8','6-9']                                         
source_1_baselines_1 = ['2-7','2-8','7-8','7-9','8-9']                                                  
source_2_baselines_0 = ['2-5','2-6','5-6','5-7','5-8','5-9','6-7','6-8','6-9']                                          
source_2_baselines_1 = ['2-7','2-8','2-9','7-8','7-9','8-9']               
source_3_baselines_0 = ['2-5','2-6','2-9','5-6','5-7','5-8','5-9','6-7','6-8','6-9'] 
source_3_baselines_1 = ['2-7','7-9','8-9']            
source_3_baselines_2 = ['2-8','7-8'] 
source_4_baselines_0 = ['2-5','2-6','2-9','5-6','5-7','5-8','5-9','6-7','6-8','6-9'] 
source_4_baselines_1 = ['2-7','7-9','8-9']            
source_4_baselines_2 = ['2-8','7-8']



source_0_parameters_0 = [5,2,5,2,3.0]        # List of parameters corresponding to each source specific
                                             # list of baselines
                                             # parameters MUST be listed in the following order:
source_0_parameters_1 = [4,2,4,2,2.8]        # [aggressiveness_first_run, max_subset_first_run, 
                                             # aggressiveness_second_run, max_subset_second_run, rho, 
                                             # kickout_sigma_level]
source_0_parameters_2 = [3,2,3,2,2.7] 
source_1_parameters_0 = [10,4,12,16,3.0]                                                   
source_1_parameters_1 = [8,4,10,8,2.8]                                                    
source_2_parameters_0 = [10,4,12,16,3.0]                                                  
source_2_parameters_1 = [8,4,10,8,2.8]                                                    
source_3_parameters_0 = [5,2,5,2,3.0]                                                   
source_3_parameters_1 = [4,2,4,2,2.8]                                                    
source_3_parameters_2 = [3,2,3,2,2.7]
source_4_parameters_0 = [5,2,5,2,3.0]                                                   
source_4_parameters_1 = [4,2,4,2,2.8]                                                   
source_4_parameters_2 = [3,2,3,2,2.7]




####################################################
#### Extra flagging options                     ####
#### These can be used to tidy up the final     ####
#### flag information before writing to a flag  ####
#### table. This will flag any unflagged data   ####
#### found between data that is already flagged ####
#### (akin to the AIPS REFLG task)              ####
#### If these are not required set both         ####
#### options to 0.                              ####
####################################################


flag_coinc_chans = 0                # This will include flag_coinc_chans number of channels between 
                                    # flagged channels when writing the flags to the table akin to REFLG
                                    # i.e. if flag_coinc_chans = 3 and chans 2-3 and 7-11 are flagged, 
                                    # channels 4,5 & 6 will be as well. Default is 0.

flag_coinc_times = 0                # This will include flag_coinc_times number of times between 
                                    # flagged times when writing the flags to the table akin to REFLG
                                    # i.e. if flag_coinc_times = 3 and times 2-3 and 7-11 are flagged, 
                                    # times 4,5 & 6 will be as well. Default is 0.
                                    
                                    




###################################
#### Outputs from SERPent:     ####
###################################


# SERPent will update the flag status of the directory directly in the ms file 



