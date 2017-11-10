#######################################################################
# SERPent4casa.py
# A version of the SERPent RFI mitigation tool for use with the
# CASA data reduction software. 
# Author Danielle Fenech
# Version Nov. 2017
#######################################################################



import LocalProxy, gc
from xmlrpclib import ServerProxy
import copy, optparse, os, sys
import re, string
import os.path
#import numarray.ma
#import numarray.ieeespecial
#import numarray.nd_image
import inspect
import warnings
warnings.defaultaction = "always"
import numpy
import numpy as np
numpy.set_printoptions(threshold='nan')
import platform
import cPickle as pickle
import math, time, datetime
from time import gmtime, strftime, localtime
ti = time.time()    # To time the script
import multiprocessing as multip
from getsize import total_size

numpy.seterr(invalid='ignore')
import warnings
warnings.filterwarnings("ignore", 'Mean of empty slice.')
warnings.filterwarnings("ignore", 'invalid value encountered in double_scalar')
warnings.filterwarnings("ignore", 'Warning: converting a masked element to nan.')


###############################################################################
# Function definitions 
################################################################################


#from serpent4casa_funcs import *
#from SERPent_funcs_IFs_alt import *

def extract_process(send_q, rec_q, cpu):

    for value in iter(send_q.get, 'STOP'):

        start, stop = value[:]
        ti1=time.time()
        data = np.array(tb.getcol('DATA',start,stop))
        ti2 = time.time()
        rec_q.put(([ti2-ti1]))
        

def flagger_process(send_q, rec_q, cpu):

    # print 'Got here'

    for value in iter(send_q.get, 'STOP'):


        ## retrieve inputs from queue ##
        #print value

        # bline,nif,source,nchan,ntimescans,polnames,experiment,name,phasecal,lovell_bline,telescope,srcname,zero_level,do_lovell_cross,coadd_polarization_flags,pickle_list,uvname,uvklass,uvdisk,uvseq, singlesource, source_num, orderedpols, drownums=value[:]

        # bline,nif,source,nchan,pol,ntimescans,polnames,path2folder,experiment,name,phasecal,lovell_bline,telescope,srcname,zero_level,do_lovell_cross,coadd_polarization_flags,pickle_list,flagging_options,parameters,uvname,uvklass,uvdisk,uvseq=value[:]
        #bline,did,source,nchan,polnames,msname,orderedpols,bpolnames,jobcount,totaljobs,path2folder = value[:]

        bline,did,source,nchan,polnames,msname,orderedpols,bpolnames,jobcount,totaljobs,path2folder,phasecal,lovell_bline,telescope,srcname,zero_level,do_lovell_cross,coadd_polarization_flags,flagging_options,parameters,source_num,oldflags_dir,antdict,do_sumthreshold = value[:]



        mycore = cpu
        JOBS = 'JOB %i of %i: ' % (jobcount, totaljobs)
        ap = time.time()
        memory_array = 0
        # print JOBS + " CPU", mycore, " Appending visibilities..."
        print JOBS + " Reading visibilities into memory..."
        # print singlesource

        #bline,nif,nchan,ntimescans=value[:]

        print 'Got values, now loading data...'
        stime = time.time()
        #ms.open(msname)
        #tb.open(msname)
        mycore = cpu
        npol = len(polnames)
        appending_time = 0


        ## Create empty arrays ready to read in the data ##

        visfile = path2folder + "data__"+ str(source_list[srce]) + '__ID' + str(did) + '__' + str(bline)+'.npy'
        # fname = path2folder+old"flags__" + str(source) +"__"+ str(bline)+"--oldflags.npy"
        data_load = np.load(visfile)
        flagsfile = oldflags_dir + "oldflags__"+ str(source_list[srce]) + '__ID' + str(did) + '__' + str(bline)+'.npy'
        old_flags = np.load(flagsfile)
        # print old_flags.shape
        timesinfile = path2folder + "times__"+ str(source_list[srce]) + '__ID' + str(did) + '__' + str(bline)+'.npy'
        times_array = np.load(timesinfile)

        


        '''
        # times_array = numpy.zeros([ntimescans], dtype='|S12')
        # print ntimescans[bline], bline
        # print len(times_array)
        npol = len(polnames)
        appending_time = 0
        ampdatadict = {}
        freqdatadict = {}
        flagdatadict = {}
        tim1 = time.time()
        ampdata =  np.zeros([npol, ntimescans, nchan], dtype=np.complex128)
        for i in xrange(len(drows)):
            flagchunk = np.array(tb.getcol('FLAG',drows[i],1))
            datachunk = np.array(tb.getcol('DATA',drows[i],1))
            ampdata[:,i,:] = datachunk[:,:,0]
        tim2 = time.time()
        # print 'All stokes as one took ', tim2-tim1
        '''
        '''
        tim3 = time.time()
        for i in xrange(len(polnames)):
            ampdatadict[i] = np.zeros([ntimescans, nchan], dtype=np.complex128)
            freqdatadict[i] = np.zeros([ntimescans, nchan], dtype=np.complex128)
            flagdatadict[i] = np.zeros([ntimescans, nchan], dtype=np.complex128)

        for i in xrange(len(drows)):
            flagchunk = np.array(tb.getcol('FLAG',drows[i],1))
            datachunk = np.array(tb.getcol('DATA',drows[i],1))
            for j in xrange(len(polnames)):
                ampdatadict[j][i,:] = datachunk[j,:,0]
            # print datachunk.shape
            # print ampdata.shape
            # print datachunk[:,:,0].shape
            # ampdata[:,i,:] = datachunk[:,:,0]
        tim4 = time.time()
        print 'All stokes as one took ', tim2-tim1, 'but sep. stokes took', tim4-tim3
        #print ampdata.shape
        # print stop
        '''
        
        if 'RR' in polnames or 'XX' in polnames:
            # amp_time_RR = numpy.zeros([ntimescans, nchan])
            # amp_freq_RR = numpy.zeros([ntimescans, nchan])
            if 'RR' in polnames:
                pol = 'RR'
            if 'XX' in polnames:
                pol = 'XX'
            # flag_RR = numpy.zeros([ntimescans, nchan])
            # print polnames[pol]
            stoke = polnames[pol]
            amp_time_RR = copy.deepcopy(np.swapaxes(data_load[stoke],0,1))
            amp_freq_RR = copy.deepcopy(np.swapaxes(data_load[stoke],0,1))
            flag_RR = copy.deepcopy(np.swapaxes(old_flags[stoke],0,1).astype(int))
            # print 'SHAPES', flag_RR.shape, amp_time_RR.shape
            # print flag_RR[2000],amp_time_RR[2000]
            appending_count_RR = 0
            # print amp_time_RR.shape
        if 'LL' in polnames or 'YY' in polnames:
            # amp_time_LL = numpy.zeros([ntimescans, nchan])
            # amp_freq_LL = numpy.zeros([ntimescans, nchan])
            # flag_LL = numpy.zeros([ntimescans, nchan])
            if 'LL' in polnames:
                pol = 'LL'
            if 'YY' in polnames:
                pol = 'YY'
            # print polnames[pol]
            stoke = polnames[pol]
            amp_time_LL = copy.deepcopy(np.swapaxes(data_load[stoke],0,1))
            amp_freq_LL = copy.deepcopy(np.swapaxes(data_load[stoke],0,1))
            flag_LL = copy.deepcopy(np.swapaxes(old_flags[stoke],0,1).astype(int))
            # flag_LL = copy.deepcopy(old_flags[stoke])
            # print flag_LL.shape, amp_time_LL.shape
            appending_count_LL = 0
        if 'RL' in polnames or 'XY' in polnames:
            # amp_time_RL = numpy.zeros([ntimescans, nchan])
            # amp_freq_RL = numpy.zeros([ntimescans, nchan])
            if 'RL' in polnames:
                pol = 'RL'
            if 'XY' in polnames:
                pol = 'XY'
            # flag_RR = numpy.zeros([ntimescans, nchan])
            # print polnames[pol]
            stoke = polnames[pol]
            amp_time_RL = copy.deepcopy(np.swapaxes(data_load[stoke],0,1))
            amp_freq_RL = copy.deepcopy(np.swapaxes(data_load[stoke],0,1))
            flag_RL = copy.deepcopy(np.swapaxes(old_flags[stoke],0,1).astype(int))
            # flag_RL = copy.deepcopy(old_flags[stoke])
            # print flag_RL.shape, amp_time_RL.shape
            appending_count_RL = 0
        if 'LR' in polnames or 'YX' in polnames:
            # amp_time_LR = numpy.zeros([ntimescans, nchan])
            # amp_freq_LR = numpy.zeros([ntimescans, nchan])
            if 'LR' in polnames:
                pol = 'LR'
            if 'YX' in polnames:
                pol = 'YX'
            # flag_RR = numpy.zeros([ntimescans, nchan])
            # print polnames[pol]
            stoke = polnames[pol]
            amp_time_LR = copy.deepcopy(np.swapaxes(data_load[stoke],0,1))
            amp_freq_LR = copy.deepcopy(np.swapaxes(data_load[stoke],0,1))
            flag_LR = copy.deepcopy(np.swapaxes(old_flags[stoke],0,1).astype(int))
            # flag_LR = copy.deepcopy(old_flags[stoke])
            # print flag_LR.shape, amp_time_LR.shape
            appending_count_LR = 0
        del data_load
        os.remove(visfile)
        os.remove(timesinfile)
        
        



        ###  new code to apply previous flags from old_flags array

        if 'RR' in polnames or 'XX' in polnames:
            amp_time_RR[flag_RR!=0] = np.nan
            amp_freq_RR[flag_RR!=0] = np.nan
            flag_RR[flag_RR!=0] = -1
        if 'LL' in polnames or 'YY' in polnames:
            amp_time_LL[flag_LL!=0] = np.nan
            amp_freq_LL[flag_LL!=0] = np.nan
            flag_LL[flag_LL!=0] = -1
        if 'RL' in polnames or 'XY' in polnames:
            amp_time_RL[flag_RL!=0] = np.nan
            amp_freq_RL[flag_RL!=0] = np.nan
            flag_RL[flag_RL!=0] = -1
        if 'LR' in polnames or 'YX' in polnames:
            amp_time_LR[flag_LR!=0] = np.nan
            amp_freq_LR[flag_LR!=0] = np.nan
            flag_LR[flag_LR!=0] = -1





        # print JOBS+"CPU", mycore, "Append time (hh:mm:ss):", time2hms(time.time()-ap)
        print JOBS+ "Append time (hh:mm:ss):", time2hms(time.time()-ap)
        if 'RR' in polnames or 'XX' in polnames:
            memory_array += 3*amp_time_RR.nbytes
        if 'LL' in polnames or 'YY' in polnames:
            memory_array += 3*amp_time_LL.nbytes
        if 'RL' in polnames or 'XY' in polnames:
            memory_array += 3*amp_time_RL.nbytes
        if 'LR' in polnames or 'YX' in polnames:
            memory_array += 3*amp_time_LR.nbytes
        # print "Total memory size of arrays: %s" % array_size(memory_array)
        # print "AMPTIMELL",amp_time_LL
        # print "times", times_array

        

        # Test to see if entire IF is Nan's from the correlator for each stoke
        # print "\n"+JOBS+" CPU", mycore, 'Checking for NaN values in array...'
        print "\n"+JOBS+'Checking for NaN values in array...'
        percentage = 1.0    # 0 to 1.0
        if 'RR' in polnames or 'XX' in polnames:
            nan_count = 0
            nan_array_RR = 0
            nan_count = np.isnan(amp_time_RR).sum()
            if nan_count == (percentage*amp_time_RR.shape[0]*amp_time_RR.shape[1]):
                nan_array_RR = 1
        if 'LL' in polnames or 'YY' in polnames:
            nan_count = 0
            nan_array_LL = 0
            nan_count = np.isnan(amp_time_LL).sum()
            if nan_count == (percentage*amp_time_LL.shape[0]*amp_time_LL.shape[1]):
                nan_array_LL = 1
        if 'RL' in polnames or 'XY' in polnames:
            nan_count = 0
            nan_array_RL = 0
            nan_count = np.isnan(amp_time_RL).sum()
            if nan_count == (percentage*amp_time_RL.shape[0]*amp_time_RL.shape[1]):
                nan_array_RL = 1
        if 'LR' in polnames or 'YX' in polnames:
            nan_count = 0
            nan_array_LR = 0
            nan_count = np.isnan(amp_time_LR).sum()
            if nan_count == (percentage*amp_time_LR.shape[0]*amp_time_LR.shape[1]):
                nan_array_LR = 1


        


        ############################################################################
        # IN SCAN ZERO LEVEL DATA CODE. WHEN TELESCOPES DROPOUT FOR UNKNOWN REASONS#
        ############################################################################
        
        if zero_level == 'yes':
            #### print "\n CPU ", mycore, "Running Zero Level Dropout Passage on unprocessed baseline...", bline
            con = np.zeros([len(times_array), 2], dtype='|S12')
            for t in xrange(len(times_array)):
                con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]
            

            if 'LL' in polnames or 'YY' in polnames:
                amp_time_LL, amp_freq_LL, flag_LL, zeros_LL = inscan_zero_level_func(amp_time_LL, amp_freq_LL, times_array, flag_LL,JOBS)
                # break
                if 'LL' in polnames:
                    pol = 'LL'
                elif 'YY' in polnames:
                    pol = 'YY'
                if len(zeros_LL) > 0:
                    plname = str(name) + "__" + str(source) + "__"+str(bline) + "__IF" + str(did+1)+ "__" + pol + "__inscan_zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_LL, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    del zeros_LL
            if 'RR' in polnames or 'XX' in polnames:
                amp_time_RR, amp_freq_RR, flag_RR, zeros_RR = inscan_zero_level_func(amp_time_RR, amp_freq_RR, times_array, flag_RR,JOBS)
                # break
                if 'RR' in polnames:
                    pol = 'RR'
                elif 'XX' in polnames:
                    pol = 'XX'
                if len(zeros_RR) > 0:
                    plname = str(name) + "__" + str(source) +"__"+ str(bline) + "__IF" + str(did+1)+ "__" + pol + "__inscan_zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_RR, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    del zeros_RR
            if 'RL' in polnames or 'XY' in polnames:
                amp_time_RL, amp_freq_RL, flag_RL, zeros_RL = inscan_zero_level_func(amp_time_RL, amp_freq_RL, times_array, flag_RL,JOBS)
                # break
                if 'RL' in polnames:
                    pol = 'RL'
                elif 'XY' in polnames:
                    pol = 'XY'
                if len(zeros_RL) > 0:
                    plname = str(name) + "__" + str(source) + "__"+str(bline) + "__IF" + str(did+1)+ "__" + pol + "__inscan_zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_RL, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    del zeros_RL
            if 'LR' in polnames or 'YX' in polnames:
                amp_time_LR, amp_freq_LR, flag_LR, zeros_LR = inscan_zero_level_func(amp_time_LR, amp_freq_LR, times_array, flag_LR,JOBS)
                # break
                if 'LR' in polnames:
                    pol = 'LR'
                elif 'YX' in polnames:
                    pol = 'YX'
                if len(zeros_LR) > 0:
                    plname = str(name) + "__"+ str(source) + "__"+ str(bline) + "__IF" + str(did+1)+ "__" + pol + "__inscan_zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_LR, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    del zeros_LR
            gc.collect()


        ##################################################################
        ########## End of in scan zero level dropouts run ################
        ##################################################################





        ########################################################
        #### LOVELL STATIONARY SCAN CODE. ONLY FOR E-MERLIN ####
        ########################################################
        
        
        ## Testing to see whether telescope is e-MERLIN and finding 
        ## Lovell antenna number:
        lovell_num = 0
        covbase = [] 
        '''
        if telescope in ('e-MERLIN', 'E-MERLIN', 'eMERLIN', 'EMERLIN'):
            #### print " Array is %s" % telescope
            for ant in xrange(len(uvdata.antennas)):
                if uvdata.antennas[ant].upper() in ('LO', 'LOV', 'LOVE', 'LOVEL', 'LOVELL'):
                    lovell_num = ant+1
                    #### print "Lovell number:", lovell_num
                    break
        # print lovell_num, bline, lovell_bline, phasecal, srcname
        '''
        ant1 = int(bline.split('-')[0])
        ant2 = int(bline.split('-')[1])
        # print ant1,ant2,antdict
        ## Add pol check, only run on parallel hands unless specify do_lovell_cross

    
        if do_lovell_cross == 'no':
            # if str(lovell_num) in bline and bline not in lovell_bline and srcname in phasecal:
            if ((antdict[ant1] == 'Lo') or (antdict[ant2] == 'Lo')) and (srcname in phasecal):
                # print "\n"+JOBS+" CPU ", mycore, "Running Lovell Stationary Scan Passage on unprocessed Lovell baseline...", bline
                print "\n"+JOBS+"Running Lovell Stationary Scan Passage on unprocessed Lovell baseline...", bline
                if 'LL' in polnames or 'YY' in polnames:
                    amp_time_LL,amp_freq_LL,flag_LL,dropouts_LL,con = lovell_dropout_med(amp_time_LL,amp_freq_LL,times_array,flag_LL,nchan,JOBS)
                    print JOBS+'Running Lovell check on LL or YY' 
                    checkedpol = 1
                    if 'LL' in polnames:
                        pol = 'LL'
                    elif 'YY' in polnames:
                        pol = 'YY'
                    if len(dropouts_LL) > 0:
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(did+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_LL, open(path2folder + str(plname) + ".dropouts", "wb"))
                        del dropouts_LL 
                        del con
                if 'RR' in polnames or 'XX' in polnames:
                    amp_time_RR,amp_freq_RR,flag_RR,dropouts_RR,con = lovell_dropout_med(amp_time_RR,amp_freq_RR,times_array,flag_RR,nchan,JOBS)
                    print JOBS+'Running Lovell check on RR or XX' 
                    checkedpol = 2
                    if 'RR' in polnames:
                        pol = 'RR'
                    elif 'YY' in polnames:
                        pol = 'XX'
                    if len(dropouts_RR) > 0:
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(did+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_RR, open(path2folder + str(plname) + ".dropouts", "wb"))
                        del dropouts_RR
                        del con
                gc.collect()

                        ## below should run on all pols ##


        elif do_lovell_cross == 'yes':
            # if str(lovell_num) in bline and bline not in lovell_bline and srcname in phasecal:
            if ((antdict[ant1] == 'Lo') or (antdict[ant2] == 'Lo')) and srcname in phasecal:
                print "\n"+JOBS+"Running Lovell Stationary Scan Passage on unprocessed Lovell baseline...", bline
                # print "\n"+JOBS+" CPU ", mycore, "Running Lovell Stationary Scan Passage on unprocessed Lovell baseline...", bline
                if 'LL' in polnames or 'YY' in polnames:
                    amp_time_LL,amp_freq_LL,flag_LL,dropouts_LL,con = lovell_dropout_med(amp_time_LL,amp_freq_LL,times_array,flag_LL,nchan,JOBS)
                    print JOBS+'Running Lovell check on LL or YY' 
                    checkedpol = 1
                    if 'LL' in polnames:
                        pol = 'LL'
                    elif 'YY' in polnames:
                        pol = 'YY'
                    if len(dropouts_LL) > 0:
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(did+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_LL, open(path2folder + str(plname) + ".dropouts", "wb"))
                        del dropouts_LL
                        del con
                if 'RR' in polnames or 'XX' in polnames:
                    amp_time_RR,amp_freq_RR,flag_RR,dropouts_RR,con = lovell_dropout_med(amp_time_RR,amp_freq_RR,times_array,flag_RR,nchan,JOBS)
                    print JOBS+'Running Lovell check on RR or XX' 
                    checkedpol = 2
                    if 'RR' in polnames:
                        pol = 'RR'
                    elif 'YY' in polnames:
                        pol = 'XX'
                    if len(dropouts_RR) > 0:
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(did+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_RR, open(path2folder + str(plname) + ".dropouts", "wb"))
                        del dropouts_RR
                        del con
                if 'RL' in polnames or 'XY' in polnames:
                    amp_time_RL,amp_freq_RL,flag_LL,dropouts_RL,con = lovell_dropout_med(amp_time_RL,amp_freq_RL,times_array,flag_RL,nchan,JOBS)
                    print JOBS+'Running Lovell check on LL or YY' 
                    checkedpol = 1
                    if 'RL' in polnames:
                        pol = 'RL'
                    elif 'XY' in polnames:
                        pol = 'XY'
                    if len(dropouts_RL) > 0:
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(did+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_RL, open(path2folder + str(plname) + ".dropouts", "wb"))
                        del dropouts_RL
                        del con
                if 'LR' in polnames or 'YX' in polnames:
                    amp_time_LR,amp_freq_LR,flag_LR,dropouts_LR,con = lovell_dropout_med(amp_time_LR,amp_freq_LR,times_array,flag_LR,nchan,JOBS)
                    print JOBS+'Running Lovell check on RR or XX' 
                    checkedpol = 2
                    if 'RR' in polnames:
                        pol = 'RR'
                    elif 'YY' in polnames:
                        pol = 'XX'
                    if len(dropouts_LR) > 0:
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(did+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_LR, open(path2folder + str(plname) + ".dropouts", "wb"))
                        del dropouts_LR
                        del con
                gc.collect()



        #######################################################
        ################# End of Lovell dropout ###############
        #######################################################




        ####################################################################
        # ZERO LEVEL DATA CODE. WHEN TELESCOPES DROPOUT FOR UNKNOWN REASONS#
        ####################################################################
        
        if zero_level == 'yes':
            print "\n CPU ", mycore, "Running Zero Level Dropout Passage on unprocessed baseline...", bline, did
            con = np.zeros([len(times_array), 2], dtype='|S12')
            for t in xrange(len(times_array)):
                con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]
            

            if 'LL' in polnames or 'YY' in polnames:
                amp_time_LL, amp_freq_LL, flag_LL, zeros_LL = zero_level_func(amp_time_LL, amp_freq_LL, flag_LL,JOBS)
                # break
                if 'LL' in polnames:
                    pol = 'LL'
                elif 'YY' in polnames:
                    pol = 'YY'
                if len(zeros_LL) > 0:
                    plname = str(name) + "__" + str(source) + "__" + str(bline) + "__IF" + str(did+1)+ "__" + pol + "__zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_LL, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    del zeros_LL
            if 'RR' in polnames or 'XX' in polnames:
                amp_time_RR, amp_freq_RR, flag_RR, zeros_RR = zero_level_func(amp_time_RR, amp_freq_RR, flag_RR,JOBS)
                # break
                if 'RR' in polnames:
                    pol = 'RR'
                elif 'XX' in polnames:
                    pol = 'XX'
                if len(zeros_RR) > 0:
                    plname = str(name) + "__" + str(source) + "__" + str(bline) + "__IF" + str(did+1)+ "__" + pol + "__zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_RR, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    del zeros_RR
            if 'RL' in polnames or 'XY' in polnames:
                amp_time_RL, amp_freq_RL, flag_RL, zeros_RL = zero_level_func(amp_time_RL, amp_freq_RL, flag_RL,JOBS)
                # break
                if 'RL' in polnames:
                    pol = 'RL'
                elif 'XY' in polnames:
                    pol = 'XY'
                if len(zeros_RL) > 0:
                    plname = str(name) + "__" + str(source) + "__" + str(bline) + "__IF" + str(did+1)+ "__" + pol + "__zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_RL, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    del zeros_RL
            if 'LR' in polnames or 'YX' in polnames:
                amp_time_LR, amp_freq_LR, flag_LR, zeros_LR = zero_level_func(amp_time_LR, amp_freq_LR, flag_LR,JOBS)
                # break
                if 'LR' in polnames:
                    pol = 'LR'
                elif 'YX' in polnames:
                    pol = 'YX'
                if len(zeros_LR) > 0:
                    plname = str(name) + "__" + str(source) + "__" + str(bline) + "__IF" + str(did+1)+ "__" + pol + "__zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    # pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    # pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_LR, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    del zeros_LR
            gc.collect()




        ##################################################################
        ############## End of zero level dropouts run ####################
        ##################################################################
        



        ##############################################################
        ############ Sum Threshold flagging sequence: ################
        ### Execute the Flagger for each polarisation and baseline ##

        if do_sumthreshold == 'yes':

            if 'RR' in polnames or 'XX' in polnames:
                if 'RR' in polnames:
                    pol = 'RR'
                if 'XX' in polnames:
                    pol = 'XX'
                # print " \n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
                print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
                if amp_time_RR.shape[0] > 0 and amp_freq_RR.shape[0] > 0:
                    # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                    print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and SPW: %i" % (pol, bline, did+1)
                    flag_RR = flagger_new(amp_time_RR, amp_freq_RR, flag_RR,flagging_options,parameters,JOBS)
                    '''
                    if write_olds:
                    flag_RR[flag_RR==-9] = 9
                    '''
                    # only need one direction test as the other will also be equal to 0.
                del amp_time_RR
                del amp_freq_RR
                if nan_array_RR == 0:
                    ##### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                    pname = str(name) +"__"+ str(source) + "__IF" + str(did+1) + "__" + str(bline) + "__" + pol
                    pickle.dump(flag_RR, open(path2folder + str(pname)+".p", "wb"))
                    pickle_list[pname] = [source, did+1, bline, times_array]
                    pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                    # np.save(path2folder+str(pname)+".p",flag_RR)
                    # pickle_list[pname] = [source, nif+1, bline, times_array]
                    # np.save(path2folder+str(pname)+".info",times_array)
                    # print len(flag)
                    del flag_RR

            if 'LL' in polnames or 'YY' in polnames:
                if 'LL' in polnames:
                    pol = 'LL'
                if 'YY' in polnames:
                    pol = 'YY'
                print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
                # print "\n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
                if amp_time_LL.shape[0] > 0 and amp_freq_LL.shape[0] > 0:
                    # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                    print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and SPW: %i" % (pol, bline, did+1)
                    flag_LL = flagger_new(amp_time_LL, amp_freq_LL, flag_LL,flagging_options,parameters,JOBS)
                    '''
                    if write_olds:
                    flag_LL[flag_LL==-9] = 9
                    '''
                del amp_time_LL
                del amp_freq_LL
                if nan_array_LL == 0:
                    #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                    pname = str(name)  +"__"+ str(source) + "__IF" + str(did+1) + "__" + str(bline) + "__" + pol
                    pickle.dump(flag_LL, open(path2folder + str(pname)+".p", "wb"))
                    pickle_list[pname] = [source, did+1, bline, times_array]
                    pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                    # np.save(path2folder+str(pname)+".p",flag_LL)
                    # pickle_list[pname] = [source, nif+1, bline, times_array]
                    # np.save(path2folder+str(pname)+".info",times_array)
                    # print len(flag)
                    # print len(flag)
                    del flag_LL

            if 'RL' in polnames or 'XY' in polnames:
                if 'RL' in polnames:
                    pol = 'RL'
                if 'XY' in polnames:
                    pol = 'XY'
                print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
                # print "\n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
                if amp_time_RL.shape[0] > 0 and amp_freq_RL.shape[0] > 0:
                    # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                    print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and SPW: %i" % (pol, bline, did+1)
                    flag_RL = flagger_new(amp_time_RL, amp_freq_RL, flag_RL,flagging_options,parameters,JOBS)
                    '''
                    if write_olds:
                    flag_RL[flag_RL==-9] = 9
                    '''
                del amp_time_RL
                del amp_freq_RL
                if nan_array_RL == 0:
                    #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                    pname = str(name)  +"__"+ str(source) + "__IF" + str(did+1) + "__" + str(bline) + "__" + pol
                    pickle.dump(flag_RL, open(path2folder+str(pname)+".p", "wb"))
                    pickle_list[pname] = [source, did+1, bline, times_array]
                    pickle.dump(times_array, open(path2folder+str(pname)+".info", "wb"))
                    # print len(flag)
                    del flag_RL

            if 'LR' in polnames or 'YX' in polnames:
                if 'LR' in polnames:
                    pol = 'LR'
                if 'YX' in polnames:
                    pol = 'YX'
                # print "\n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
                print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
                if amp_time_LR.shape[0] > 0 and amp_freq_LR.shape[0] > 0:
                    # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                    print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and SPW: %i" % (pol, bline, did+1)
                    flag_LR = flagger_new(amp_time_LR,amp_freq_LR,flag_LR,flagging_options,parameters,JOBS)
                    '''
                    if write_olds:
                    flag_LR[flag_LR==-9] = 9
                    '''
                del amp_time_LR
                del amp_freq_LR
                if nan_array_LR == 0:
                    #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                    pname = str(name) +"__"+ str(source)  + "__IF" + str(did+1) + "__" + str(bline) + "__" + pol
                    pickle.dump(flag_LR, open(path2folder+str(pname)+".p", "wb"))
                    pickle_list[pname] = [source, did+1, bline, times_array]
                    pickle.dump(times_array, open(path2folder+str(pname)+".info", "wb"))
                    # print len(flag)
                    del flag_LR


            # num_flagged_vis = (flag != 0).sum()
            #### print "Total number of flags: %i for IF %i and stokes %s" % (num_flagged_vis, nif+1,pol)
            # tsum = 0
            # tsum += ntimescans[bline]
            # total_tscans = 0
            # total_tscans += tsum
            # print "Total number of amplitudes: %i" % (tsum*nif*nchan*npol)
            # if num_flagged_vis > 0:
            #  print 'tsum is', tsum
            #print "Percentage of amplitudes flagged for IF %i, baseline %s: %f%%, stokes %s" % (nif+1, bline, (100.0*num_flagged_vis/(tsum*nchan)),pol)
                
            finish_str = "\n"+JOBS+" Finished flagging Source :" + source + ", SPW:" +str(did)+", Baseline:"+bline
            print finish_str

            # return finish_str
        
            rec_q.put(finish_str, cpu)
            gc.collect()



        else:
            finish_str = "\n"+JOBS+" Not running SumThreshold flagger on :" + source + ", SPW:" +str(did)+", Baseline:"+bline
            rec_q.put(finish_str, cpu)
            gc.collect()


        '''
        if 'RR' in polnames or 'XX' in polnames:
            if 'RR' in polnames:
                pol = 'RR'
            if 'XX' in polnames:
                pol = 'XX'
            # print " \n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
            print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
            if amp_time_RR.shape[0] > 0 and amp_freq_RR.shape[0] > 0:
                # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and SPW: %i" % (pol, bline, did+1)
                flag_RR = flagger_new(amp_time_RR, amp_freq_RR, flag_RR,flagging_options,parameters,JOBS)

                # only need one direction test as the other will also be equal to 0.
            del amp_time_RR
            del amp_freq_RR
            if nan_array_RR == 0:
                #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                pname = str(name) +"__"+ str(source) + "__IF" + str(did+1) + "__" + str(bline) + "__" + pol
                pickle.dump(flag_RR, open(path2folder + str(pname)+".p", "wb"))
                pickle_list[pname] = [source, did+1, bline, times_array]
                pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                #np.save(path2folder+str(pname)+".p",flag_RR)
                #pickle_list[pname] = [source, nif+1, bline, times_array]
                #np.save(path2folder+str(pname)+".info",times_array)
                # print len(flag)
                del flag_RR

        if 'LL' in polnames or 'YY' in polnames:
            if 'LL' in polnames:
                pol = 'LL'
            if 'YY' in polnames:
                pol = 'YY'
            print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
            # print "\n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
            if amp_time_LL.shape[0] > 0 and amp_freq_LL.shape[0] > 0:
                # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and SPW: %i" % (pol, bline, did+1)
                flag_LL = flagger_new(amp_time_LL, amp_freq_LL, flag_LL,flagging_options,parameters,JOBS)

            del amp_time_LL
            del amp_freq_LL
            if nan_array_LL == 0:
                #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                pname = str(name)  +"__"+ str(source) + "__IF" + str(did+1) + "__" + str(bline) + "__" + pol
                pickle.dump(flag_LL, open(path2folder + str(pname)+".p", "wb"))
                pickle_list[pname] = [source, did+1, bline, times_array]
                pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                #np.save(path2folder+str(pname)+".p",flag_LL)
                #pickle_list[pname] = [source, nif+1, bline, times_array]
                #np.save(path2folder+str(pname)+".info",times_array)
                # print len(flag)
                # print len(flag)
                del flag_LL

        if 'RL' in polnames or 'XY' in polnames:
            if 'RL' in polnames:
                pol = 'RL'
            if 'XY' in polnames:
                pol = 'XY'
            print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
            # print "\n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
            if amp_time_RL.shape[0] > 0 and amp_freq_RL.shape[0] > 0:
                # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and SPW: %i" % (pol, bline, did+1)
                flag_RL = flagger_new(amp_time_RL, amp_freq_RL, flag_RL,flagging_options,parameters,JOBS)

            del amp_time_RL
            del amp_freq_RL
            if nan_array_RL == 0:
                #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                pname = str(name)  +"__"+ str(source) + "__IF" + str(did+1) + "__" + str(bline) + "__" + pol
                pickle.dump(flag_RL, open(path2folder+str(pname)+".p", "wb"))
                pickle_list[pname] = [source, did+1, bline, times_array]
                pickle.dump(times_array, open(path2folder+str(pname)+".info", "wb"))
                # print len(flag)
                del flag_RL

        if 'LR' in polnames or 'YX' in polnames:
            if 'LR' in polnames:
                pol = 'LR'
            if 'YX' in polnames:
                pol = 'YX'
            # print "\n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
            print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
            if amp_time_LR.shape[0] > 0 and amp_freq_LR.shape[0] > 0:
                # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and SPW: %i" % (pol, bline, did+1)
                flag_LR = flagger_new(amp_time_LR,amp_freq_LR,flag_LR,flagging_options,parameters,JOBS)

            del amp_time_LR
            del amp_freq_LR
            if nan_array_LR == 0:
                #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                pname = str(name) +"__"+ str(source)  + "__IF" + str(did+1) + "__" + str(bline) + "__" + pol
                pickle.dump(flag_LR, open(path2folder+str(pname)+".p", "wb"))
                pickle_list[pname] = [source, did+1, bline, times_array]
                pickle.dump(times_array, open(path2folder+str(pname)+".info", "wb"))
                # print len(flag)
                del flag_LR


        # num_flagged_vis = (flag != 0).sum()
        #### print "Total number of flags: %i for IF %i and stokes %s" % (num_flagged_vis, nif+1,pol)
        # tsum = 0
        # tsum += ntimescans[bline]
        # total_tscans = 0
        # total_tscans += tsum
        # print "Total number of amplitudes: %i" % (tsum*nif*nchan*npol)
        # if num_flagged_vis > 0:
           #  print 'tsum is', tsum
            #print "Percentage of amplitudes flagged for IF %i, baseline %s: %f%%, stokes %s" % (nif+1, bline, (100.0*num_flagged_vis/(tsum*nchan)),pol)
                
        finish_str = "\n"+JOBS+" Finished flagging Source :" + source + ", SPW:" +str(did)+", Baseline:"+bline
        print finish_str

        # return finish_str
        
        rec_q.put(finish_str, cpu)
        gc.collect()
        '''








#############################################################################
######### End of main parallel flagger_process definition ###################
#############################################################################




#############################################################################
################## Start of New SumThreshold flagger definition #################

def flagger_new(amp_time, amp_freq, flag, flagging_options, parameters,JOBS):

    '''
    SumThreshold flagging process. Takes in desired amplitude and flag data arrays 
    then uses the flagging_options and parameters to flag the data.
    Outputs updated amplitude and flag arrays.
    '''

    global_total = 0
    global_time = 0
    global_freq = 0



    if flagging_options == 'default':
        aggressiveness_first_run = 4
        max_subset_first_run = 2
        aggressiveness_second_run = 4
        max_subset_second_run = 2
        rho = 1.5
        kickout_sigma_level = 3.0
    else:
        aggressiveness_first_run = parameters[0]
        max_subset_first_run = parameters[1]
        aggressiveness_second_run = parameters[2]
        max_subset_second_run = parameters[3]
        rho = parameters[4]
        kickout_sigma_level = parameters[5]
        # print parameters[5]
         
    # Initial run with subset size = 1 to remove v. strong RFI
    # Note only one direction needed window = 1

    # Need to consider any values the correlator has flagged. 
    # Written them in SERPent as 'nan'


    # Updated so that read_in_flags will turn old_flag data to nan's but maintain 
    # their positions in the flag array so they can still be written out if required
    
    # flag[np.isnan(amp_time)] = np.where(flag[np.isnan(amp_time)]==-9, flag[np.isnan(amp_time)],-1)
    flag[np.isnan(amp_time)] = np.where(np.logical_or(flag[np.isnan(amp_time)]==-6,flag[np.isnan(amp_time)]==-7), flag[np.isnan(amp_time)],-1)
    # flag[np.isnan(amp_time)] = -1.0

    # flag = np.where(amp_time!="nan",flag,-1.0)
    #print np.where(amp_time[0]!="nan")

    median = np.median(amp_time[~np.isnan(amp_time)])
    final_mad = mad(amp_time)
    chi_1 = median + (20*final_mad)

    i = 1
    if flagging_options == 'default':
        if freq_band(frequency) in ('L Band' or 'S Band'):
            kickout = median + 3.0*final_mad
            # print "Kickout aggressive for low frequency..."
        else:
            kickout = median + 3.0*final_mad
            # print "Kickout less aggressive for high frequency..."
            ### User defined kickout sigma level:
    else:
        # print kickout_sigma_level
        kickout = median + kickout_sigma_level*final_mad



    while i < amp_time.shape[0]+1 and i <= 1:
        timing = time.time()
        limit = chi_1/(rho**(math.log(i,2)))
        ## The below condition is to avoid flagging good data if the 
        ## thresholding gets too close to the mean
        if limit <= kickout:
            break
        # print amp_time, '\n \n \n \n'
        amp_time = np.where(flag!=1,amp_time,limit)
        # print amp_time
        amp_freq = np.where(flag!=1,amp_freq,limit)
        for row in xrange(1,amp_time.shape[0]-(i-2)):
            if i == 1:
                flag[row-1,:][np.where(amp_time[row-1,:]>limit)] = np.where(flag[row-1,:][np.where(amp_time[row-1,:]>limit)]<0,flag[row-1,:][np.where(amp_time[row-1,:]>limit)],1)
        i *= 2

    
    count_i = 0
    # count_i = np.count_nonzero(flag)
    count_i = (flag > 0).sum()
    old_count = (flag != 0).sum()
    # print "Flagged with initial run:", count_i, 'old_count', old_count

    ## Run the full Threshold algorithm for the first time.
    ## Flagged vis will be at lowest threshold for the following statistics:

    hm = time.time()
    median = np.median(amp_time[~np.isnan(amp_time)])
    final_mad = mad(amp_time)
    chi_1 = median + (aggressiveness_first_run*final_mad)

    ## Reset the amplitude arrays for the new statistics.
    ## Flagged visibilities get set to the first threshold level because 
    ## the one within the while loop is for setting limits within that 
    ## sumthreshold run
    
    amp_time = np.where(flag<=0,amp_time,chi_1)
    amp_freq = np.where(flag<=0,amp_freq,chi_1)

    ###amp_time = np.where(flag==0.0,amp_time,chi_1)
    ###amp_freq = np.where(flag==0.0,amp_freq,chi_1)

    '''
    testflag2out = open('/local/python/RFI/SERPent/memory_tests_0316/testflag2out_file.txt','ab')
    print >> testflag2out, JOBS , 'flag', flag
    testflag2out.close()
    
    testamp2out = open('/local/python/RFI/SERPent/memory_tests_0316/testamp2out_file.txt','ab')
    print >> testamp2out, JOBS , 'amp', amp_time
    testamp2out.close()
    
    print 'type', amp_time.dtype
    '''
    i = 1
    # print "\n FLAGGING IN TIME..."
    while i < amp_time.shape[0]+1 and i < max_subset_first_run+1:
        timing = time.time()
        limit = chi_1/(rho**(math.log(i,2)))
        ## The below condition is to avoid flagging good data if the 
        ## thresholding gets too close to the mean
        if limit <= kickout:
            break
        amp_time = np.where(flag!=2,amp_time,limit)
        amp_time = np.where(flag!=-2,amp_time,limit)
        for row in xrange(1, amp_time.shape[0]-(i-2)):
            if i == 1:
                flag[row-1,:][np.where(amp_time[row-1,:]>limit)] = np.where(flag[row-1,:][np.where(amp_time[row-1,:]>limit)]<0,flag[row-1,:][np.where(amp_time[row-1,:]>limit)],2)
                flag[row-1,:][np.where(amp_time[row-1,:]>limit)] = np.where(flag[row-1,:][np.where(amp_time[row-1,:]>limit)]>=0,flag[row-1,:][np.where(amp_time[row-1,:]>limit)],-2)
            else:
                nsum = amp_time[row-1:row-1+i,:].sum(0)
                flag[row-1:row-1+i,np.where(nsum/i>limit)] = np.where(flag[row-1:row-1+i,np.where(nsum/i>limit)]<0,flag[row-1:row-1+i,np.where(nsum/i>limit)],2)
                flag[row-1:row-1+i,np.where(nsum/i>limit)] = np.where(flag[row-1:row-1+i,np.where(nsum/i>limit)]>=0,flag[row-1:row-1+i,np.where(nsum/i>limit)],-2)
        i *= 2

    '''
    testflag2aftout = open('/local/python/RFI/SERPent/memory_tests_0316/testflag2aftout_file.txt','ab')
    print >> testflag2aftout, JOBS , 'flag', flag
    testflag2aftout.close()
    '''

    count_t = (flag == 2).sum()
    print "Flagged time:", count_t

    i = 1
    # print "\n FLAGGING IN FREQUENCY..."
    while i < amp_freq.shape[1]+1 and i < max_subset_first_run+1:
        timing = time.time()
        limit = chi_1/(rho**(math.log(i,2)))
        ## The below condition is to avoid flagging good data if the 
        ## thresholding gets too close to the mean
        if limit <= kickout:
            break
        amp_freq = np.where(flag!=3,amp_freq,limit)
        amp_freq = np.where(flag!=-3,amp_freq,limit)
        for col in xrange(1, amp_freq.shape[1]-(i-2)):
            if i == 1:
                flag[:,col-1][np.where(amp_freq[:,col-1]>limit)] = np.where(flag[:,col-1][np.where(amp_freq[:,col-1]>limit)]<0,flag[:,col-1][np.where(amp_freq[:,col-1]>limit)],3)
                flag[:,col-1][np.where(amp_freq[:,col-1]>limit)] = np.where(flag[:,col-1][np.where(amp_freq[:,col-1]>limit)]>=0,flag[:,col-1][np.where(amp_freq[:,col-1]>limit)],-3)
            else:  
                nsum = amp_freq[:,col-1:col-1+i].sum(1)
                flag[np.where(nsum/i>limit),col-1:col-1+i] = np.where(flag[np.where(nsum/i>limit),col-1:col-1+i]<0,flag[np.where(nsum/i>limit),col-1:col-1+i],3)
                flag[np.where(nsum/i>limit),col-1:col-1+i] = np.where(flag[np.where(nsum/i>limit),col-1:col-1+i]>=0,flag[np.where(nsum/i>limit),col-1:col-1+i],-3)
        i *= 2


    count_f = (flag == 3).sum()
    print "Flagged freq:", count_f

    
    ## second run of flagger removing already found strong RFI flags 
    ## and calculating thresholds with new stats
    ## if (really high max value a certain magnitude above average/median 
    ## then run flagger again removing the strong RFIs from the arrays for 
    ## the sumThreshold...)

    median1 = np.median(amp_time)
    median = np.median(amp_time[~np.isnan(amp_time)])
    final_mad = mad(amp_time)
    maxvalue = np.amax(amp_time)
    chi_1 = median + (aggressiveness_second_run*final_mad)


    #print median1, median
    
    amp_time = np.where(flag<=0,amp_time,chi_1)
    amp_freq = np.where(flag<=0,amp_freq,chi_1)


    ###amp_time = np.where(flag==0.0,amp_time,chi_1)
    ###amp_freq = np.where(flag==0.0,amp_freq,chi_1)

    i = 1
    # print "\n FLAGGING IN TIME..."
    if flagging_options == 'default':
        if freq_band(frequency) in ('L Band' or 'S Band'):
            kickout = median + 3.0*final_mad
            # print "Kickout aggressive for low frequency..."
        else:
            kickout = median + 3.0*final_mad
            # print "Kickout less aggressive for high frequency..."
    ## User defined kickout sigma level:
    else:
        kickout = median + kickout_sigma_level*final_mad
    while i < amp_time.shape[0]+1 and i < max_subset_second_run+1:
        if maxvalue > kickout:
            if count_t > 0 or count_f > 0:
                pass
            else:
                print JOBS+" Nothing flagged. Skipping run..."
                break
        else:
            print JOBS+" Max value < kickout. Skipping run..."
            break
        timing = time.time()
        limit = chi_1/(rho**(math.log(i,2)))
 
        ## The below condition is to avoid flagging good data if the 
        ## thresholding gets too close to the mean
        if limit <= kickout:
            #print 'breaking...'
            break
        amp_time = np.where(flag!=4,amp_time,limit)
        amp_time = np.where(flag!=-4,amp_time,limit)
        
        for row in xrange(1, amp_time.shape[0]-(i-2)):
            if i == 1:
                flag[row-1,:][np.where(amp_time[row-1,:]>limit)] = np.where(flag[row-1,:][np.where(amp_time[row-1,:]>limit)]<0,flag[row-1,:][np.where(amp_time[row-1,:]>limit)],4)
                flag[row-1,:][np.where(amp_time[row-1,:]>limit)] = np.where(flag[row-1,:][np.where(amp_time[row-1,:]>limit)]>=0,flag[row-1,:][np.where(amp_time[row-1,:]>limit)],-4)
            else:
                nsum = amp_time[row-1:row-1+i,:].sum(0)
                flag[row-1:row-1+i,np.where(nsum/i>limit)] = np.where(flag[row-1:row-1+i,np.where(nsum/i>limit)]<0,flag[row-1:row-1+i,np.where(nsum/i>limit)],4)
                flag[row-1:row-1+i,np.where(nsum/i>limit)] = np.where(flag[row-1:row-1+i,np.where(nsum/i>limit)]>=0,flag[row-1:row-1+i,np.where(nsum/i>limit)],-4)
        i *= 2


    count_t = (flag == 4).sum()
    #print "Flagged time:", count_t

    i = 1
    # print "\n FLAGGING IN FREQUENCY..."
    if flagging_options == 'default':
        if freq_band(frequency) in ('L Band' or 'S Band'):
            kickout = median + 3.0*final_mad
            # print "Kickout aggressive for low frequency..."
        else:
            kickout = median + 3.0*final_mad
            # print "Kickout less aggressive for high frequency..."
    ## User defined kickout sigma level:
    else:
        kickout = median + kickout_sigma_level*final_mad
    while i < amp_freq.shape[1]+1 and i < max_subset_second_run+1:
        if maxvalue > kickout:
            if count_t > 0 or count_f > 0:
                pass
            else:
                print JOBS+" Nothing flagged. Skipping run..."
                break
        else:
            print JOBS+" Max value < kickout. Skipping run..."
            break
        timing = time.time()
        limit = chi_1/(rho**(math.log(i,2)))
        ## The below condition is to avoid flagging good data if the 
        ## thresholding gets too close to the mean
        if limit <= kickout:
            break
        amp_freq = np.where(flag!=5,amp_freq,limit)
        amp_freq = np.where(flag!=-5,amp_freq,limit)
        for col in xrange(1, amp_freq.shape[1]-(i-2)):
            if i == 1:
                flag[:,col-1][np.where(amp_freq[:,col-1]>limit)] = np.where(flag[:,col-1][np.where(amp_freq[:,col-1]>limit)]<0,flag[:,col-1][np.where(amp_freq[:,col-1]>limit)],5)
                flag[:,col-1][np.where(amp_freq[:,col-1]>limit)] = np.where(flag[:,col-1][np.where(amp_freq[:,col-1]>limit)]>=0,flag[:,col-1][np.where(amp_freq[:,col-1]>limit)],-5)
            else:  
                nsum = amp_freq[:,col-1:col-1+i].sum(1)
                flag[np.where(nsum/i>limit),col-1:col-1+i] = np.where(flag[np.where(nsum/i>limit),col-1:col-1+i]<0,flag[np.where(nsum/i>limit),col-1:col-1+i],5)
                flag[np.where(nsum/i>limit),col-1:col-1+i] = np.where(flag[np.where(nsum/i>limit),col-1:col-1+i]>=0,flag[np.where(nsum/i>limit),col-1:col-1+i],-5)
        i *= 2


    count_f = (flag == 5).sum()
    # print "Flagged freq:", count_f
    

    # flagged_total = np.count_nonzero(flag)
    flagged_total = (flag > 0.0).sum()
    #### print "Total number of flagged amplitudes: %i" % flagged_total
    # global global_total
    global_total += flagged_total
    # global global_time
    global_time += count_t
    # global global_freq
    global_freq += count_f
    # print "CPU", mycore, "time (secs):", time.time() - hm
    #### print "time (secs):", time.time() - hm

    return (flag)


################ End of New SumThreshold flagger definition ##################
##########################################################################




##########################################################################
############ Start of in scan zero-level drop-out definition #############
##########################################################################


def inscan_zero_level_func(amp_time,amp_freq,times_array,flag_arr,JOBS):
    
    '''
    Inscan_zero_level dropout function. Takes as inputs the amplitude and 
    associated time and flag arrays. Performs a running mad stastical analysis 
    looking for data consistently less than 3*mad within a single scan.
    outputs updated amplitude and flag arrays.

    '''

    zeros = []

 
    con = np.zeros([len(times_array), 2], dtype='|S12')
    # print times_array
    test = len(times_array)
    # print test
        
    for t in xrange(len(times_array)):
        con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]
            
    int_time = np.zeros([len(times_array)-1])
    for t in xrange(1, len(times_array)):
        int_time[t-1] = float(times_array[t]) - float(times_array[t-1])
                        
    step = float(np.median(int_time))
    del int_time
    # print JOBS+" Correlator/ Averaged integration time is:", time2dhms(step), 'seconds'
        
    dropouts = []

    scantimes = []
    start = 0
    for j in xrange(1,len(times_array)):
        if float(con[j][0]) - float(con[j-1][0]) > 5*step:
            end = j
            scantimes.append([start,end])
            start = end
        elif j == len(times_array)-1:
            end = j+1
            # print 'yes',end
            scantimes.append([start,end])
    
    scantimes = np.array(scantimes)
    # print scantimes

    for i in xrange(len(scantimes)):
        start = scantimes[i][0]
        end = scantimes[i][1]

    
        newamp_time = amp_time[start:end]
        mdat = np.ma.masked_array(newamp_time,np.isnan(newamp_time))
        # print 'MDAT is',len(mdat),mdat.shape, mdat
        single = np.ma.median(mdat, axis=1)*10**6

        #single = single

        zero_avg = np.ma.average(single[~np.isnan(single)])
        zero_std = np.ma.std(single[~np.isnan(single)])
        zero_median = np.ma.median(single[~np.isnan(single)])
        zero_mad = mad_1d(single[~np.isnan(single)])
        # print zero_avg,zero_std,zero_median,zero_mad

        # print single
        # print 'checking zeros',zeros
        for v in xrange(len(single)):
            # print v, v+start, single[v],zero_avg,zero_std,zero_median,zero_mad
            if single[v] < (0 + 3*zero_std) and single[v] < (zero_avg - 2*zero_std):
                # print 'time1:', times_array[v+start],v; alt_times(times_array[v+start])
                if v+start not in zeros:
                    # print v+start
                    zeros.append(v+start)
            if single[v] < (zero_median - 3*zero_mad):
                # print 'time2:', times_array[v+start],v; alt_times(times_array[v+start])
                if v+start not in zeros:
                    zeros.append(v+start)
            if single[v] < (zero_median - 3*zero_mad) and single[v] < (0 + 2*zero_mad):
                # print 'time3:', times_array[v+start],v; alt_times(times_array[v+start])
                if v not in zeros:
                    zeros.append(v+start)
            
        ## To get the Cambridge zero levels which affect 75% of observation?
        #if single[v] < (0 + 5*zero_std):
        # if single[v] < (zero_median - zero_mad):
        # if v not in zeros:
        # zeros.append(v)
        # print "Dropouts flagged:", pol, 'IF'+str(j+1), len(zeros)
        # print sorted(zeros), len(zeros)

    #for i in xrange(len(single)):
     #   if i in zeros:
      #      amp_time[i,:] = zero_median
       #     amp_freq[i,:] = zero_median
            # print  bline,str(j+1),pol+':',sorted(zeros), len(zeros)

    zeros = sorted(zeros)
    
    amp_time[zeros] = 'nan'
    amp_freq[zeros] = 'nan'
    flag_arr[zeros] = -6.
    #print UHOH
    # print len(amp_time)
    # print len(flag_arr)
    # return(amp_time, amp_freq, times, flag_arr, zeros)
    # print amp_time,amp_freq,flag_arr
    del con
    return(amp_time, amp_freq, flag_arr, zeros)



########## End of in-scan zero-level drop-out definition #############
######################################################################





##########################################################################
################ Start of zero-level drop-out definition #################
##########################################################################


def zero_level_func(amp_time,amp_freq,flag_arr,JOBS):

    '''
    Zero_level dropout function. Takes as inputs the amplitude and 
    associated time and flag arrays. Performs a running mad stastical analysis 
    looking for data consistently less than 3*mad over all scans.
    Outputs updated amplitude and flag arrays.

    '''

    
    zeros = []

    mdat = np.ma.masked_array(amp_time,np.isnan(amp_time))
    # print 'MDAT is',len(mdat),mdat.shape, mdat
    single = np.ma.median(mdat, axis=1)*10**6
    del mdat

    # con = np.zeros([len(times_array), 2], dtype='|S12')
    # for t in xrange(len(times_array)):
        # con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]

    #single = single

    zero_avg = np.ma.average(single[~np.isnan(single)])
    zero_std = np.ma.std(single[~np.isnan(single)])
    zero_median = np.ma.median(single[~np.isnan(single)])
    zero_mad = mad_1d(single[~np.isnan(single)])
    # print 'zero_avg,zero_std,zero_median,zero_mad', zero_avg,zero_std,zero_median,zero_mad,pol
    # print (0 + 3*zero_std), (zero_avg - 2*zero_std),pol
    # print (zero_median - 2*zero_mad),pol
    # print (zero_median - 2*zero_mad), (0 + 2*zero_mad),pol
    # print (zero_median - zero_mad)

    # print single
    # print 'checking zeros',zeros
    for v in xrange(len(single)):
        if single[v] < (0 + 3*zero_std) and single[v] < (zero_avg - 2*zero_std):
            #print 'time1:', times_array[v],v; alt_times(times_array[v])
            if v not in zeros:
                zeros.append(v)
        if single[v] < (zero_median - 3*zero_mad):
            #print 'time2:', times_array[v],v; alt_times(times_array[v])
            if v not in zeros:
                zeros.append(v)
        if single[v] < (zero_median - 3*zero_mad) and single[v] < (0 + 2*zero_mad):
            #print 'time3:', times_array[v],v; alt_times(times_array[v])
            if v not in zeros:
                zeros.append(v)
            
        ## To get the Cambridge zero levels which affect 75% of observation?
        #if single[v] < (0 + 5*zero_std):
        # if single[v] < (zero_median - zero_mad):
        # if v not in zeros:
        # zeros.append(v)

    # print "Dropouts flagged:", pol, 'IF'+str(j+1), len(zeros)
    # print sorted(zeros), len(zeros)

    for i in xrange(len(single)):
        if i in zeros:
            amp_time[i,:] = zero_median
            amp_freq[i,:] = zero_median
            # print  bline,str(j+1),pol+':',sorted(zeros), len(zeros)

    zeros = sorted(zeros)
    # print amp_time
    # print zeros
    # times_array = np.delete(times_array, zeros)
    # amp_time = np.delete(amp_time, zeros, axis=0)
    # amp_freq = np.delete(amp_freq, zeros, axis=0)
    # flag_arr = np.delete(flag_arr, zeros, axis=0)

    amp_time[zeros] = 'nan'
    amp_freq[zeros] = 'nan'
    flag_arr[zeros] = -6.


    #print len(amp_time)
    #print len(flag_arr)
    # return(amp_time, amp_freq, times, flag_arr, zeros)
    del single
    return(amp_time, amp_freq, flag_arr, zeros)



############## End of zero-level drop-out definition #################
######################################################################





#####################################################################
############# Start of mean Lovell drop-out definition ##############
#####################################################################


def lovell_dropout(amp_time,amp_freq,times_array,flag_arr,nchan,JOBS):

        '''
        Lovell dropout function. Takes as inputs the amplitude and 
        associated time and flag arrays. Performs a running mad stastical analysis 
        looking for data consistently low over all scans. Using the mean as reference.
        This is designed to pick up instances where the eMERLIN Lovell telescope skips 
        phase calibration observations becasue of slow-slewing times.
        Outputs updated amplitude and flag arrays.
        
        '''


        ## Select a single channel to perform test on. Choosing channel in the 
        ## middle of the IF for best quality amps...
        # test_chan = int(nchan/2)
        
        ## choose central 50% chans to average ##
        chan_start = int(nchan/4)
        chan_end = chan_start + int(nchan/2)
        print JOBS+'Checking channels:', chan_start, chan_end
        
        con = np.zeros([len(times_array), 2], dtype='|S12')
          
        for t in xrange(len(times_array)):
            # try:
                # print float(times_array[t])
            # except ValueError,e:
                # print "error",e,"on line",t,'value',times_array[t]
            con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]
            
        int_time = np.zeros([len(times_array)-1])
        for t in xrange(1, len(times_array)):
            int_time[t-1] = float(times_array[t]) - float(times_array[t-1])
                        
        # print 'Lovell polarisation is', pol, j+1 
        step = float(np.median(int_time))
        # print JOBS+"Correlator/ Averaged integration time is:", time2dhms(step), 'seconds'
        del int_time

        checkedpol = 0
        dropouts = []
        
        onechans = np.ma.masked_array(amp_time[:,chan_start:chan_end],np.isnan(amp_time[:,chan_start:chan_end]))
        one = np.ma.median(onechans, axis=1)*10**6
        del onechans

        # print 'continue'
        # print 'checking pol', checkedpol
        
        levels = {'high': [0,0], 'low': [10**10,0]}
        
        ## Firstly cycle through baseline and find mean, std of scans and 
        ## set the lowest and highest levels.
        
        # print "ONE:",mycore, one
        count_drop = 0
        begin_step = 0
        end_step = 0
        for i in xrange(1,len(one)):
            # print i, one[i], con[i][0], con[i][1], con[i-1][0]
            # Some test level of 5 * integration time...
            if float(con[i][0]) - float(con[i-1][0]) > 5*step:
                # print "Scan has changed level"
                # print "loop1"
                end_step = i
                newo = one[begin_step:end_step]
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                if avg > levels['high'][0]:
                    levels['high'][:] = [avg, s]
                if avg < levels['low'][0]:
                    levels['low'][:] = [avg, s]
                begin_step = end_step
                count_drop += 1
                # print 'loop1 levels:', levels
            if i == len(one)-1:
                # print "Last scan"
                # print "loop2"
                end_step = i
                newo = one[begin_step:end_step+1]
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                if avg > levels['high'][0]:
                    levels['high'][:] = [avg, s]
                if avg < levels['low'][0]:
                    levels['low'][:] = [avg, s]
                begin_step = end_step
                # print 'loop2 levels:', levels
                # print "avg:", avg,s


        # now flag all scan levels within lowest avg + lowest s
        count_drop = 0
        begin_step = 0
        end_step = 0
        # print len(one)
        for i in xrange(1,len(one)):
            # print i, one[i], con[i][0], con[i][1], one[0],one[1]
            if float(con[i][0]) - float(con[i-1][0]) > 5*step:
                # print "loop3"
                end_step = i
                newo = one[begin_step:end_step]
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                newmed = np.ma.median(newo[~np.isnan(newo)])
                newmad = mad_1d(newo[~np.isnan(newo)])
                print JOBS+"Average: %f, std: %f" % (avg, s), newmed,newmad
                # print levels
                if levels['low'][0] + levels['low'][1] < levels['high'][0]-(2*levels['high'][1]):
                    if avg < levels['low'][0] + (6*levels['low'][1]) and avg > levels['low'][0] - (6*levels['low'][1]):
                        one[begin_step:end_step] = 0.0
                        # print 'flagging'
                        count_drop += 1
                begin_step = end_step
            if i == len(one)-1:
                # print "loop4"
                end_step = i
                newo = one[begin_step:end_step+1]
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                # print "Average: %f, std: %f" % (avg, s)
                # print levels
                #### print "Average: %f, std: %f" % (avg, s)
                if levels['low'][0] + levels['low'][1] < levels['high'][0]-(2*levels['high'][1]):
                    if avg < levels['low'][0] + (6*levels['low'][1]) and avg > levels['low'][0] - (6*levels['low'][1]):
                        one[begin_step:end_step+1] = 0.0
                        # print 'flagging'
                        count_drop += 1
                begin_step = end_step
                # print end_step, one
                
                # Now collect indexes for amp == 0.0
        
        # print count_drop
        if count_drop > 0:
            dropouts = []
            
            # print 'len one is', len(one)
            #dropoutsf = np.array(np.where(one==0.0))
            #dropouts = dropoutsf[0]
            #print dropoutsf

            for i in xrange(len(one)):
                #print i,one[i]
                # print i+1, one[i], con[i][0], con[i][1]
                if one[i] == 0.0:
                    dropouts.append(i)
                #if i == len(one)-1:
                 #   if one[i] == 0.0:
                  #      dropouts.append(i+1)
        # print dropouts

        del newo
        del one
        dropouts = sorted(dropouts)
        # print dropouts
        # print len(dropouts)
        # old_times = copy.copy(times_array)

        # times_array = np.delete(times_array, dropouts)
        # amp_time = np.delete(amp_time, dropouts, axis=0)
        # amp_freq = np.delete(amp_freq, dropouts, axis=0)
        # flag_arr = np.delete(flag_arr, dropouts, axis=0)
        # print len(amp_time), len(amp_freq), len(flag_arr),len(one)

        # old_times = np.delete(old_times, dropouts)
        amp_time[dropouts] = 'nan'
        amp_freq[dropouts] = 'nan'
        flag_arr[dropouts] = -7.
 
                

        #### print "Before Lovell:", test, "after Lovell:", len(old_times)
        return(amp_time,amp_freq,flag_arr,dropouts,con)
                                



############## End of mean Lovell droput definition ################
####################################################################





#####################################################################
########### Start of median Lovell drop-out definition ##############
#####################################################################


def lovell_dropout_med(amp_time,amp_freq,times_array,flag_arr,nchan,JOBS):


        '''
        Lovell dropout median function. Takes as inputs the amplitude and 
        associated time and flag arrays. Performs a running sad stastical analysis 
        looking for data consistently low over all scans. Using the median as reference.
        This is designed to pick up instances where the eMERLIN Lovell telescope skips 
        phase calibration observations becasue of slow-slewing times.
        Outputs updated amplitude and flag arrays.
        
        '''



        ## Select a single channel to perform test on. Choosing channel in the 
        ## middle of the IF for best quality amps...
        # test_chan = int(nchan/2)
        
        ## choose central 50% chans to average ##
        chan_start = int(nchan/4)
        chan_end = chan_start + int(nchan/2)
        print JOBS+'Checking channels:', chan_start, chan_end
        
        con = np.zeros([len(times_array), 2], dtype='|S12')
         
        for t in xrange(len(times_array)):
            # try:
                # print float(times_array[t])
            # except ValueError,e:
                # print "error",e,"on line",t,'value',times_array[t]
            con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]
            
        int_time = np.zeros([len(times_array)-1])
        for t in xrange(1, len(times_array)):
            int_time[t-1] = float(times_array[t]) - float(times_array[t-1])
                        
        # print 'Lovell polarisation is', pol, j+1 
        step = float(np.median(int_time))
        # print JOBS+"Correlator/ Averaged integration time is:", time2dhms(step), 'seconds'
        del int_time

        checkedpol = 0
        dropouts = []
        
        onechans = np.ma.masked_array(amp_time[:,chan_start:chan_end],np.isnan(amp_time[:,chan_start:chan_end]))
        one = np.ma.median(onechans, axis=1)*10**6
        del onechans

        # print 'continue'
        # print 'checking pol', checkedpol
        
        levels = {'high': [0,0], 'low': [10**10,0]}
        
        ## Firstly cycle through baseline and find mean, std of scans and 
        ## set the lowest and highest levels.
        
        # print "ONE:",mycore, one
        count_drop = 0
        begin_step = 0
        end_step = 0
        for i in xrange(1,len(one)):
            # print i, one[i], con[i][0], con[i][1], con[i-1][0]
            # Some test level of 5 * integration time...
            if float(con[i][0]) - float(con[i-1][0]) > 5*step:
                # print "Scan has changed level"
                # print "loop1"
                end_step = i
                newo = one[begin_step:end_step]
                avg = np.ma.median(newo[~np.isnan(newo)])
                s = mad_1d(newo[~np.isnan(newo)])
                if avg > levels['high'][0]:
                    levels['high'][:] = [avg, s]
                if avg < levels['low'][0]:
                    levels['low'][:] = [avg, s]
                begin_step = end_step
                count_drop += 1
                # print 'loop1 levels:', levels
            if i == len(one)-1:
                # print "Last scan"
                # print "loop2"
                end_step = i
                newo = one[begin_step:end_step+1]
                avg = np.ma.median(newo[~np.isnan(newo)])
                s = mad_1d(newo[~np.isnan(newo)])
                if avg > levels['high'][0]:
                    levels['high'][:] = [avg, s]
                if avg < levels['low'][0]:
                    levels['low'][:] = [avg, s]
                begin_step = end_step
                # print 'loop2 levels:', levels
                # print "avg:", avg,s


        # now flag all scan levels within lowest avg + lowest s
        count_drop = 0
        begin_step = 0
        end_step = 0
        # print len(one)
        for i in xrange(1,len(one)):
            # print i, one[i], con[i][0], con[i][1], one[0],one[1]
            if float(con[i][0]) - float(con[i-1][0]) > 5*step:
                # print "loop3"
                end_step = i
                newo = one[begin_step:end_step]
                med = np.ma.median(newo[~np.isnan(newo)])
                mad = mad_1d(newo[~np.isnan(newo)])
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                # print "Average: %f, std: %f" % (avg, s)
                # print levels
                if levels['low'][0] + levels['low'][1] < levels['high'][0]-(2*levels['high'][1]):
                    # print levels['low'][0], levels['high'][0]-levels['high'][1]
                    # print levels['low'][0] + levels['low'][1], levels['low'][0] - levels['low'][1]
                    if avg < levels['low'][0] + (9*levels['low'][1]) and avg > levels['low'][0] - (9*levels['low'][1]):
                        # print avg, levels['low'][0] + levels['low'][1] , levels['low'][0] - levels['low'][1]
                        one[begin_step:end_step] = 0.0
                        count_drop += 1
                begin_step = end_step
            if i == len(one)-1:
                # print "loop4"
                end_step = i
                newo = one[begin_step:end_step+1]
                med = np.ma.median(newo[~np.isnan(newo)])
                mad = mad_1d(newo[~np.isnan(newo)])
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                # print "Average: %f, std: %f" % (avg, s)
                # print levels
                #### print "Average: %f, std: %f" % (avg, s)
                if levels['low'][0] + levels['low'][1] < levels['high'][0]-(2*levels['high'][1]):
                    if avg < levels['low'][0] + (9*levels['low'][1]) and avg > levels['low'][0] - (9*levels['low'][1]):
                        one[begin_step:end_step+1] = 0.0
                        count_drop += 1
                begin_step = end_step
                # print end_step, one
                
                # Now collect indexes for amp == 0.0

        # print count_drop
        if count_drop > 0:
            dropouts = []
            
            # print 'len one is', len(one)
            for i in xrange(len(one)):
                # print i,one[i]
                # print i+1, one[i], con[i][0], con[i][1]
                if one[i] == 0.0:
                    dropouts.append(i)
                # if i == len(one)-1:
                #   if one[i] == 0.0:
                #    dropouts.append(i+1)
            # print dropouts

        del newo
        del one
        dropouts = sorted(dropouts)
        # old_times = copy.copy(times_array)

        # times_array = np.delete(times_array, dropouts)
        # amp_time = np.delete(amp_time, dropouts, axis=0)
        # amp_freq = np.delete(amp_freq, dropouts, axis=0)
        # flag_arr = np.delete(flag_arr, dropouts, axis=0)

        # old_times = np.delete(old_times, dropouts)
        amp_time[dropouts] = 'nan'
        amp_freq[dropouts] = 'nan'
        flag_arr[dropouts] = -7.
        # del onechans
                

        #### print "Before Lovell:", test, "after Lovell:", len(old_times)
        return(amp_time,amp_freq,flag_arr,dropouts,con)
                                



############ End of median Lovell droput definition ################
####################################################################




###############################################################################
# Global Functions
################################################################################



## Definition to find computer memory stats (works only for csh and bash 
## at the moment)...
def computer_memory(path2folder):
    '''
    Will find the available computer memory information where possible
    '''
    import os
    memfile = open(path2folder + 'mem_stats.txt', 'wr')
    memfile.flush()
    file = path2folder + 'mem_stats.txt'
    os.system('free -k >' + file)    # need to test this!
    #os.system('free -k > mem_stats.txt')    # need to put in file path
    memfile = open(path2folder + 'mem_stats.txt', 'r')
    stats = []
    for line in memfile:
        string = ''
        for char in xrange(len(line)):
            if line[char].isdigit():
                string += line[char]
                if line[char+1] == ' ' or line[char+1] == '\n':
                    stats.append(int(string))
                    string = ''
    global mem_stats
    mem_stats = {}
    mem_stats['mem_total'] = stats[0]
    mem_stats['mem_used'] = stats[1]
    mem_stats['mem_free'] = stats[2]
    mem_stats['mem_shared'] = stats[3]
    mem_stats['mem_buffers'] = stats[4]
    mem_stats['mem_cached'] = stats[5]
    mem_stats['buf/cach_used'] = stats[6]
    mem_stats['buf/cach_free'] = stats[7]
    mem_stats['swap_total'] = stats[8]
    mem_stats['swap_used'] = stats[9]
    mem_stats['swap_free'] = stats[10]
    memfile.close()
    return mem_stats




## Frequency band calculator. Only works for freq between 1.2GHz - 26.5GHz 
## but easy to add on more...
def freq_band(freq):
    '''
    Will determine eMERLIN frequency band for information.
    '''
    if freq > 1.2E+9 and freq < 2.0E+9:
        return "L Band"
    elif freq > 2.0E+9 and freq < 4.0E+9:
        return "S Band"
    elif freq > 4.0E+9 and freq < 8.0E+9:
        return "C Band"
    elif freq > 8.0E+9 and freq < 12.0E+9:
        return "X Band"
    elif freq > 12.0E+9 and freq < 18.0E+9:
        return "Ku Band"
    elif freq > 18.0E+9 and freq < 26.5E+9:
        return "K Band"
    else:
        return "Out of range for calculator"


## Frequency prefix calculator
def freq_prefix(num):
    '''
    Will determine an appropriate formatting 
    for the data frequency.
    '''
    if num/(10**12) > 1:
        return "%.5f THz" % (num / 10.0**12)
    elif num/(10**9) > 1:
        return "%.5f GHz" % (num / 10.0**9)
    elif num/(10**6) > 1:
        return "%.5f MHz" % (num / 10.0**6)
    elif num/(10**3) > 1:
        return "%.5f KHz" % (num / 10.0**3)





## Determines the array size
def array_size(data_mem):
    '''
    Determines the memory required to store the specified
    numpy array
    '''
    ## Returns the size of the array with the appropriate suffix
    if len(str(data_mem)) <= 6:
        return "%.3f KB" % (data_mem / 10.0**3)
    elif len(str(data_mem)) <= 9:
        return "%.3f MB" % (data_mem / 10.0**6)
    elif len(str(data_mem)) <= 12:
        return "%.3f GB" % (data_mem / 10.0**9)
    elif len(str(data_mem)) <= 15:
        return "%.3f TB" % (data_mem / 10.0**12)





## Function to convert time to d/hh:mm:ss.s format which is what is required for FG tables
def time2dhms(seconds):
    '''
    Function to convert time to d/hh:mm:ss.s format
    '''
    conv = seconds * 86400
    d = int(conv/ 86400)
    h=int(conv % 86400)/ 3600
    m=int(conv % 3600)/60
    s=conv-(d*86400)-(h*3600)-(m*60)
    d=`d`
    h=`h`
    m=`m`
    s="%3.1f" % s
    dhms= d.zfill(1) + "/" + h.zfill(2) + ":" + m.zfill(2) + ":" + s.zfill(4)
    return dhms



## Function to convert seconds to hh:mm:ss.ss format, returns a string
def time2hms(seconds):
    '''
    Function to convert seconds to hh:mm:ss.ss format, returns a string
    '''
    h=int(seconds/3600)
    m=int(seconds % 3600)/60
    s=seconds-(h*3600)-(m*60)
    h=`h`
    m=`m`
    s="%4.2f" % s
    hms=h.zfill(2)+":"+m.zfill(2)+":"+s.zfill(4)
    return hms







def mad(array):
    '''
    Calculates the median absolute deviation of the input array
    '''
    ab_dev = np.zeros_like(array)
    coef = 1.4286
    # median_i = float(np.median(array, axis=None))
    median_i = float(np.median(array[~np.isnan(array)], axis=None))
    ab_dev[:,:] = abs(array[:,:] - median_i)
    # median_j = float(np.median(ab_dev, axis=None))
    median_j = float(np.median(ab_dev[~np.isnan(array)], axis=None))
    final_mad = coef*median_j
    # print "The Median Absolute Deviation for this baseline is %f" % final_mad
    return final_mad






def mad_1d(array):
    '''
    Calculates the median absolute deviation of the input array in one dimension
    '''
    ab_dev = np.zeros_like(array)
    coef = 1.4286
    median_i = float(np.ma.median(array[~np.isnan(array)]))
    ab_dev[:] = abs(array[:] - median_i)
    median_j = float(np.ma.median(ab_dev[~np.isnan(ab_dev)]))
    final_mad = coef*median_j
    return final_mad




def alt_times(time):
    '''
    Function to print time information in human-readable format.
    '''

    time=float(time)
    rhrs = (time*24.0)
    hrs = math.floor(rhrs)
    rmins = (rhrs-hrs)*60.0
    mins = math.floor(rmins)
    rsecs = (rmins-mins)*60.0
    secs = math.floor(rsecs)
    # print hrs,mins,rsecs
    # print '%.2f:%.2f:%.2f' %(hrs,mins,rsecs)




def coinc_chan_flags_alt(flagnew,coincs):

    '''
    Function to add flags in input flag array for instances where 'coincs' indices are encompassed
    by True flag entries.
    Outputs an updated flag array.
    '''
    import copy
    ncoincs = coincs+1
    # print 'Flagging %i coincident channels' % coincs
    # new_flag = copy.copy(flagnew)
    
    for col in xrange(flagnew.shape[1]-1):
        if flagnew[np.where((flagnew[:,col] > 0) & (flagnew[:,col+1] <= 0))].size > 0:
            short_arr = flagnew[np.where((flagnew[:,col] > 0) & (flagnew[:,col+1] <= 0)),col+1:col+ncoincs+1]
            if np.any(np.where(short_arr>0)):
                test_arr = np.where(short_arr>0,flagnew[np.where((flagnew[:,col] > 0) & (flagnew[:,col+1] <= 0)),col+1:col+ncoincs+1],1)
                farr = np.argmax(np.any(short_arr>0,axis=0),axis=1)
                for k in xrange(len(farr)):
                    short_arr[0,k,0:farr[k]] = 1
                flagnew[np.where((flagnew[:,col] > 0) & (flagnew[:,col+1] <= 0)),col+1:col+ncoincs+1]=short_arr
    return(flagnew)



def coinc_flags(flagnew,time_coincs,chan_coincs):
    
    '''
    Function to add flags in input flag array for instances where 'coincs' indices are encompassed
    by True flag entries. Performs the process iteratively over each dimension, updating across channels
    first.
    Outputs an updated flag array.
    '''

    import copy
    # print 'Flagging %i coincident channels and %i coincident times' % (chan_coincs,  time_coincs)
    new_flag = copy.copy(flagnew)
    if time_coincs > 0:
        new_flag = coinc_chan_flags_alt(new_flag,chan_coincs)
    if chan_coincs > 0:
        new_flag = new_flag.transpose()
        new_flag = coinc_chan_flags_alt(new_flag,time_coincs)
        new_flag = new_flag.transpose()
    return(new_flag)

def coinc_flags_noncumu(flagnew,time_coincs,chan_coincs):
    
    '''
    Function to add flags in input flag array for instances where 'coincs' indices are encompassed
    by True flag entries. Performs the process iteratively over each dimension but in a non-cumulative
    fashion, updating across channels first.
    Outputs an updated flag array.
    '''


    import copy
    # print 'Flagging %i coincident channels and %i coincident times' % (chan_coincs,  time_coincs)
    new_flag = copy.copy(flagnew)
    if time_coincs > 0:
        new_flag = coinc_chan_flags_alt(new_flag,chan_coincs)
    if chan_coincs > 0:
        new_flag2 = copy.copy(flagnew.transpose())
        new_flag2 = coinc_chan_flags_alt(new_flag2,time_coincs)
        new_flag2 = new_flag2.transpose()
        new_flag[np.where(new_flag2==1)]=1
    return(new_flag)





# Code to flatten lists into single dimension. Code from #
# http://rightfootin.blogspot.co.uk/2006/09/more-on-python-flatten.html #
def flatten(l, ltypes=(list, tuple)):
    '''
    Function to flatten lists into single dimension. Code from 
    http://rightfootin.blogspot.co.uk/2006/09/more-on-python-flatten.html
    '''

    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)




###############################################################################
# End of Global Functions
################################################################################






'''
etime = time.time()
# rec_q.put(([tim2-tim1,tim4-tim3]))
rec_q.put(([etime-stime,etime-stime]))
# etime = time.time()
print 'Total time for job ', etime-stime
'''








###############################################################################
###############################################################################
## Main SERPent Script 
###############################################################################
###############################################################################


version_name = "Aesculapian"
version_number = '1.1'
version_date = "10/11/17"

#print sys.argv
if len(sys.argv) > 3:
    if sys.argv[3] == '--v' or sys.argv[3] == '--version':
        print 'SERPent4casa version is %s: %s dated %s' % (version_name, version_number, version_date)
        sys.exit(0)
    else:
        print 'ERROR: unrecognised option'
        sys.exit(0)


print '\n Started running SERPent4casa version %s: %s' % (version_name, version_number), 'on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),'\n'


## Execute the serpent4casa_input.py file with all the observation information 
## and flagging parameters...
try:
    execfile("SERPent4casa_input.py")
    
except:
    print "\n Could not find SERPent4casa_input.py!"
    print " Make sure input file is in the same folder as this script (SERPent.py)"
    print " Aborting!"
    sys.exit(0)


### Need a declaration of all input variables here to check they exist



## Test to see whether all variables have been set:

## Test to see whether all variables have been set:
try:
    Name
    path2data
    NCPU
    path2folder
    spwlist
    path2folder = os.path.join(path2folder, '')
    oldflags_dir = path2folder
    if os.path.exists(path2folder) == False:
        print " Folder for outputs does not exist! Please check inputs and try again."
        print " Aborting SERPent!\n"
        sys.exit(0)
    do_Lovell_scan
    phasecal
    if do_Lovell_scan == 'no':
        phasecal = []
    do_lovell_cross
    zero_level
    which_baselines
    flagging_options
    flag_sources
    if flag_sources == 'choose':
        flag_list
    do_sumthreshold
    if flagging_options == 'choose' or flagging_options == 'source':
        aggressiveness_first_run
        max_subset_first_run
        aggressiveness_second_run
        max_subset_second_run
        rho
        kickout_sigma_level
    if which_baselines == 'choose':
        baselines
    dobline_specific
    dosource_specific
    select_stokes
    if select_stokes == 'yes':
        flag_all_stokes 
        stokes
    coadd_polarization_flags
    if coadd_polarization_flags == 'no':
        coadd_zero_flags
        coadd_lovell_flags
    flag_coinc_chans
    flag_coinc_times
    zero_flag_allif
    # keep_oldflags
except NameError:
    print " Please check you\'ve set all the variables in the input file."
    print " ABORTING SERPent!\n"
    sys.exit(0)





## Information about the observations:

name=Name
msname=path2data+Name
'''
tempdir = path2folder+'filestore/'
try:
    os.makedirs(tempdir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise
'''
#name = msname
print "\n RUNNING SCRIPT ON DATA: %s" % name

ms.open(msname,nomodify=False)
metadata=ms.metadata()
# print metadata
nvis = int(metadata.nrows())


source_lookup=metadata.namesforfields()
# print source_lookup
source_list={}
# source_lookup = [x for x in source_lookup if x]

for i in xrange(len(source_lookup)):
    srce = source_lookup[i]
    if not srce == '':
        source_list[i] = srce

source_lookup = [x for x in source_lookup if x]


# print source_lookup
# print source_list
# print source_list.keys()
source_temp = copy.copy(source_lookup)


if flag_sources == 'choose':
    # print flag_list
    for srce in source_temp:
        if srce in flag_list:
            # print 'yes',srce
            continue
        else:
            # print 'no',srce
            source_lookup.remove(srce)
            for k in source_list.keys():
                if source_list[k] == srce:
                    # print k, source_list[k]
                    del source_list[k]
del(source_temp)
#print flag_list
#print source_lookup
#print source_list
#print stop

print " Flagging sources: %s" % ', '.join(str(x) for x in source_lookup)
#print stop



## Name of Telescope:
# NEED TO CHANGE
# telescope = uvdata.header['telescop']
# print " Name of Telescope: %s" % telescope

telescope = metadata.observatorynames()

# sets up local variable
#NEED TO CHANGE
# nvis = len(uvdata)    ## prints the number of visibilities
# print " Total number of visibilities (nvis): %i" % nvis

## Find the Frequency of the observation
#NEED TO CHANGE
# frequency = uvdata.header['crval'][uvdata.header['ctype'].index('FREQ')]
# print " Frequency of the first spw: %s" % freq_prefix(frequency)
# print " Frequency band of the Observation: %s" % freq_band(frequency)

'''
## Figure out number of Stokes and polarizations
npol = len(uvdata.stokes)
print " Number of Polarizations (npol): %i" % npol
polnames = {}
for x in xrange(npol):
    polnames[str(uvdata.stokes[x])] = x
## Below is just to print out the stokes in a nice way :)
string = ""
for x in xrange(npol):
    if x != npol-1:
        string += (str(uvdata.stokes[x])+ ", ")
    else:
        string += str(uvdata.stokes[x])
print " Polarizations: %s" % string

# Get the stokes in AIPS array order #
orderedpols = get_ordered_pols(uvdata)
'''


## Figure out polarisations ##

# CASA metadata.corrtypesforpol(0) i.e. correlation types for polarisation setup 0
# CASA has inbuilt reference list for corrtype numbers and stokes parms.
# Taken from here: 
# http://casa.nrao.edu/active/docs/doxygen/html/classcasa_1_1Stokes.html#ae3cb0ef26262eb3fdfbef8273c455e0c
# Dictionary created from above list

casacorrdir = {0:'', 1:'I', 2:'Q', 3:'U', 4:'V', 5:'RR', 6:'RL', 7:'LR', 8:'LL', 9:'XX', 10:'XY', 11:'YX', 12:'YY', 13:'RX', 14:'RY', 15:'LX', 16:'LY', 17:'XR', 18:'XL', 19:'YR', 20:'YL', 21:'PP', 22:'PQ', 23:'QP', 24:'QQ', 25:'RCircular', 26:'LCircular', 27:'Linear', 28:'Ptotal', 29:'Plinear', 30:'PFtotal', 31:'PFlinear',32:'Pangle'}

# Assuming only one polarisation setup for now #
# Also assuming data ordered in polarisation in the same order as given in corrtypes #
corrtypes = metadata.corrtypesforpol(0)
polnames = {}
bpolnames = {}
orderedpols = {}
for i in xrange(len(corrtypes)):
    polnames[casacorrdir[corrtypes[i]]] = i
    bpolnames[i] = casacorrdir[corrtypes[i]]
    orderedpols[casacorrdir[corrtypes[i]]] = corrtypes[i]

# print polnames
# print bpolnames
# print orderedpols

## If selected stokes requested, remove extra from polnames, bpolnames and orderedpols
if select_stokes == 'yes':
    if not stokes:
        print ' No stokes specified, processing all'
    else:
        for pol in polnames.keys():
            if pol not in stokes:
                # print pol
                # print 'no'
                del polnames[pol]
                del orderedpols[pol]
        for key in list(bpolnames.keys()):
            if bpolnames[key] not in stokes:
                del bpolnames[key]


# print polnames
# print bpolnames
# print orderedpols
print " Flagging pols: %s" % ', '.join(str(x) for x in polnames.keys())


## Figure out number of IFs
nif = metadata.nspw()
if spwlist == 'all':
    descids = metadata.datadescids()
else:
    descids = []
    for spw in spwlist:
        niflist = metadata.datadescids(spw=spw)
        #print niflist, 'niflist[i]'
        for polset in niflist:
            descids.append(polset) 
# print descids, 'descids'


chanlist=[]
for i in xrange(nif):
    chanlist.append(metadata.nchan(i))
# print chanlist
print " Flagging spws: %s" % ', '.join(str(x) for x in descids)
print " Number of spws with channels: %i" % nif, '--', ', '.join(str(x) for x in chanlist)




### Load in flagging parameters ###


## If basic flagging options used, load in parameters ##

source_base_parms = []
if flagging_options == 'choose':
    if len(aggressiveness_first_run) > 1 or len(max_subset_first_run) > 1 or len(aggressiveness_second_run) > 1 or len(max_subset_second_run) > 1 or len(rho) > 1 or len(kickout_sigma_level) > 1:
        print "\n ERROR: flagging_options = 'choose' set but multiple instead of single flagging \n ERROR: parameters given, please check inputs. \n ERROR: If separate options required for multiple sources, \n ERROR: set flagging_options = 'source'\n"
        sys.exit(0)
    else:
        source_parms = [aggressiveness_first_run[0], max_subset_first_run[0], aggressiveness_second_run[0], max_subset_second_run[0], rho[0], kickout_sigma_level[0]]
        for i in xrange(len(source_lookup)):
            source_base_parms.append(source_parms)
        
elif flagging_options == 'source':
    noptions = len(source_lookup)
    if len(aggressiveness_first_run) != noptions or len(max_subset_first_run) != noptions or len(aggressiveness_second_run) != noptions or len(max_subset_second_run) != noptions or len(rho) != noptions or len(kickout_sigma_level) != noptions:
        print "\n ERROR: flagging_options = 'source' set but number of parameters given does not \n ERROR: match the number of sources specified for flagging. \n ERROR: If flag_sources = 'all', the number of parameters should match \n ERROR: the number of sources in the data. \n ERROR: If flag_sources = 'choose', the number of parameters should match \n ERROR: the number of sources given in flag_list. \n"
        sys.exit(0)
    else:
        source_base_parms = []
        for i in xrange(noptions):
            newparm = [aggressiveness_first_run[i],max_subset_first_run[i],aggressiveness_second_run[i],max_subset_second_run[i],rho[i],kickout_sigma_level[i]]
            source_base_parms.append(newparm)


## If using source and/or baseline specific options, load in the parameters ###

if dobline_specific == 'no':
    src_baselines = [[] for i in range(len(source_list))]
    for srcething in source_list:
        sdi = source_lookup.index(source_list[srcething])
        source_num = srcething
        src_baselines[sdi] = baselines
        #    for i in xrange(len(source_lookup)):
        #       src_baselines.append(baselines)
    # print src_baselines

elif dobline_specific == 'yes' and dosource_specific == 'no':
    source_base_parms = [[] for i in range(len(source_list))]
    src_baselines = [[] for i in range(len(source_list))]
    for srcething in source_list:
        sdi = source_lookup.index(source_list[srcething])
        lbaselines = []
        full_parms_dict = {}
        for i in xrange(1000):
            try:
                lbaselines.extend(globals()['baselines_'+str(i)]) 
                for base in globals()['baselines_'+str(i)]:
                    full_parms_dict[base] = globals()['parameters_'+str(i)]
            except: continue
        if full_parms_dict:
            source_base_parms[sdi]=full_parms_dict
            src_baselines[sdi]=lbaselines
    # print src_baselines, source_base_parms



elif dosource_specific == 'yes':
    source_base_parms = [[] for i in range(len(source_list))]
    src_baselines = [[] for i in range(len(source_list))]
    for srcething in source_list:
        sdi = source_lookup.index(source_list[srcething])
        base_parms_dict = {}
        lbaselines=[]
        # print sdi
        for i in xrange(1000):
            try:
                lbaselines.extend(globals()['source_'+str(sdi)+'_baselines_'+str(i)])
                for base in globals()['source_'+str(sdi)+'_baselines_'+str(i)]:
                    base_parms_dict[base] = globals()['source_'+str(sdi)+'_parameters_'+str(i)] 
                # print base_parms_dict
            except:
                continue
        if base_parms_dict:
            source_base_parms[sdi]=base_parms_dict
            src_baselines[sdi]=lbaselines







#### Backing up old flag state before do any flagging #####

print '\nBacking up old flag state before doing any flagging'
ctim = time.localtime()
flagfilename = 'pre_SERPent_'+str(ctim.tm_year)+'-'+str(ctim.tm_mon)+'-'+str(ctim.tm_mday)+'-'+str(ctim.tm_hour)+'.'+str(ctim.tm_min)+'.'+str(ctim.tm_sec)+'.flag'
print '\nSaving flags with flagmanager as '+ flagfilename
flagmanager(vis=msname,mode='save',versionname=flagfilename)







## Figure out number of baselines
# NEED TO CHANGE
# print " Number of Possible Baselines (nbase): %i" % nacbase



### REMOVE AUTO-CORRELATIONS AND LO-MK2 ###
# may/may not need code here, but may be better further in #




#### SECTION WHERE GET ROW NUMBERS FOR EACH SOURCE/BASELINE ####

## Dictionary containing number of time scans for each baseline.
if dobline_specific == 'no' and dosource_specific == 'no' and which_baselines == 'all':
    
    print "\n ALL baselines have been selected for flagging..." 
    
    ## List of baselines. Note that this uses the actual visibility data to find the 
    ## antenna numbers and doesn't assume they start at 1 and increase linearly by 
    ## the same amount. This has been written because the VLA and early commissioning 
    ## e-MERLIN data have a unique (i.e. stupid) numbering system.
    ## Note that this is only in the 'all' baseline choice as it then needs to find 
    ## all the baselines but doesn't if they've been stated before...

    # print 'flag_list is:', flag_list


    ## New 11/2015 code, assume have removed any sources from source_list and source_lookup
    ## that are not in the requested source list from the inputs file.

    sbline = [[] for i in range(len(source_list))]
    nbases = [0 for i in range(len(source_list))]
    baseline_dict = [{} for i in range(len(source_list))]
    dtimescans = [[] for i in range(len(source_list))]
    drownums = [[] for i in range(len(source_list))]

    ti=time.time()
    tb.open(msname)
    ants1 = tb.getcol('ANTENNA1',0,nvis)
    ants2 = tb.getcol('ANTENNA2',0,nvis)
    sources = tb.getcol('FIELD_ID',0,nvis)
    spws = tb.getcol('DATA_DESC_ID',0,nvis)
    # times = tb.getcol('TIME',0,nvis)
    # tflags = tb.getcol('FLAG',0,nvis)
    time1=time.time()
    
    antdict = {}
    ti1 = time.time()
    # print ti1-ti
    for i in xrange(nvis):
        srce = sources[i]
        did = spws[i]
        if srce in source_list:
            sdi = source_lookup.index(source_list[srce])
            bline =  "%i-%i" % (ants1[i], ants2[i])
            #            bline = "%i-%i" % (vis.baseline[0], vis.baseline[1])
            #print srce, sdi, source_list[srce], bline
            if bline not in sbline[sdi]:
                if ants1[i] not in antdict:
                    antdict[ants1[i]] = metadata.antennanames(ants1[i])[0]
                if ants2[i] not in antdict:
                    antdict[ants2[i]] = metadata.antennanames(ants2[i])[0]
                sbline[sdi].append(bline)
                nbases[sdi] += 1
                baseline_dict[sdi][bline] = 1
                '''
                dtimescans[sdi].append([])
                drownums[sdi].append([])
                for j in xrange(len(descids)):
                    dtimescans[sdi][sbline[sdi].index(bline)].append([])
                    drownums[sdi][sbline[sdi].index(bline)].append([])
                if did in descids:
                    dtimescans[sdi][sbline[sdi].index(bline)][descids.index(did)].append(times[i])
                    #                dtimescans[sdi][sbline[sdi].index(bline)].append(vis.time)
                    drownums[sdi][sbline[sdi].index(bline)][descids.index(did)].append(i)
                '''
            else:
                baseline_dict[sdi][bline] += 1
                '''
                if did in descids:
                    dtimescans[sdi][sbline[sdi].index(bline)][descids.index(did)].append(times[i])
                    drownums[sdi][sbline[sdi].index(bline)][descids.index(did)].append(i)
                '''




    # print time.time()-ti1
    #print sys.getsizeof(dtimescans[0]), sys.getsizeof(drownums[0])
    #print total_size(dtimescans, verbose=False)
    # print times.nbytes
    # print len(times)
    

## Testing for specified baselines.
if which_baselines == 'choose' or dobline_specific == 'yes' or dosource_specific == 'yes':
    print "\n SPECIFIC baselines selected for flagging..."
    
    sbline = [[] for i in range(len(source_list))]
    nbases = [0 for i in range(len(source_list))]
    baseline_dict = [{} for i in range(len(source_list))]
    dtimescans = [[] for i in range(len(source_list))]
    drownums = [[] for i in range(len(source_list))]

    tb.open(msname)
    ants1 = np.array(tb.getcol('ANTENNA1',0,nvis))
    ants2 = tb.getcol('ANTENNA2',0,nvis)
    sources = tb.getcol('FIELD_ID',0,nvis)
    spws = tb.getcol('DATA_DESC_ID',0,nvis)
    times = tb.getcol('TIME',0,nvis)
    # may need to convert times to seconds.... nope they should be fine
    antdict = {}

    for i in xrange(nvis):
        srce = sources[i]
        did = spws[i]
        if srce in source_list:
            sdi = source_lookup.index(source_list[srce])
            bline  = "%i-%i" % (ants1[i],ants2[i])
            if (bline not in sbline[sdi]) and (bline in src_baselines[sdi]):
                if ants1[i] not in antdict:
                    antdict[ants1[i]] = metadata.antennanames(ants1[i])[0]
                if ants2[i] not in antdict:
                    antdict[ants2[i]] = metadata.antennanames(ants2[i])[0]
                sbline[sdi].append(bline)
                nbases[sdi] += 1
                baseline_dict[sdi][bline] = 1
                '''
                dtimescans[sdi].append([])
                drownums[sdi].append([])
                # dtimescans[sdi][sbline[sdi].index(bline)].append(times[i])
                # drownums[sdi][sbline[sdi].index(bline)].append(i)
                for j in xrange(len(descids)):
                    dtimescans[sdi][sbline[sdi].index(bline)].append([])
                    drownums[sdi][sbline[sdi].index(bline)].append([])
                if did in descids:
                    dtimescans[sdi][sbline[sdi].index(bline)][descids.index(did)].append(times[i])
                    drownums[sdi][sbline[sdi].index(bline)][descids.index(did)].append(i)
                '''
            elif bline in src_baselines[sdi]:
                baseline_dict[sdi][bline] += 1
                '''
                if did in descids:
                    dtimescans[sdi][sbline[sdi].index(bline)][descids.index(did)].append(times[i])
                    drownums[sdi][sbline[sdi].index(bline)][descids.index(did)].append(i)
                # dtimescans[sdi][sbline[sdi].index(bline)].append(times[i])
                # drownums[sdi][sbline[sdi].index(bline)].append(i)
                '''
 

#print drownums
#print drownums.shape

# Check data exists for all sources #
nodata_list = []
for srce in source_list:
    sdi = source_lookup.index(source_list[srce])
    # print srce, source_list[srce]
    if not baseline_dict[sdi]:
        print "\n WARNING: No data found for source %s, it will not be flagged!\n" % source_list[srce]
        nodata_list.append(source_list[srce])
# print baseline_dict




        


if not np.all(nbases):
    print 'WARNING: Different number of baselines per source!'
    #print np.all(nbases)

else:

    #### Hopefully new bit of code #####
    print ' Writing data and flags to temporary files...\n'

    srce_arr_posn = [[] for i in range(len(source_list))]
    for srce in source_list:
        sdi = source_lookup.index(source_list[srce])
        # print sdi
        nblines = len(sbline[sdi])
        bline_posn =  {}
        for bline in sbline[sdi]:
            ant1 = int(bline.split('-')[0])
            ant2 = int(bline.split('-')[1])
            did_arr = [[] for i in range(len(descids))]
            for did in descids:
                ifp = descids.index(did)
                data = []
                # datalist = []
                # print srce,sdi,bline,did,source_list,ant1,ant2
                # print drownums[sdi][sbline[sdi].index(bline)][descids.index(did)][0]
                # print baseline_dict[sdi]
                # print len(sources),len(spws),len(ants1),len(ants2)
                # print np.where((sources==srce) & (spws==did) & (ants1[:]==ant1))
                select_arr = np.where((sources==srce) & (spws==did) & (ants1[:]==ant1) & (ants2[:] == ant2))[0]
                testdiffs = np.diff(select_arr)
                # print select_arr
                # print select_arr.shape
                # print testdiffs.shape, testdiffs.max()
                if (testdiffs>nblines).any():
                    # print 'Multiple scans detected'
                    list_of_starts = np.where(testdiffs>nblines)[0]
                    # print list_of_starts
                    # print select_arr[0], (testdiffs[0]+1-select_arr[0]), len()
                    # data.append(tb.getcol('DATA',select_arr[0],(testdiffs[0]+1-select_arr[0])/nblines,nblines))
                    array_posn = []
                    for i in xrange(len(list_of_starts)+1):
                        # data = []
                        if i == 0:
                            start=select_arr[0]
                            end = select_arr[list_of_starts[i]]
                            patht = 1
                            # print start,((end-start)/nblines)+1,nblines
                            if (antdict[ant1] == 'Lo') and (antdict[ant2] == 'Mk2'):
                                array_posn.append([0,((end-start)/nblines),i])
                                prev_end = ((end-start)/nblines)+1
                                lovelltotal = select_arr.shape[0]
                            elif (antdict[ant1] == 'Mk2') and (antdict[ant2] == 'Lo'):
                                array_posn.append([0,((end-start)/nblines),i])
                                prev_end = ((end-start)/nblines)+1
                                lovelltotal = select_arr.shape[0]
                            else:
                                data = tb.getcol('DATA',start,((end-start)/nblines)+1,nblines)
                                old_flags = tb.getcol('FLAG',start,((end-start)/nblines)+1,nblines)
                                times_array = tb.getcol('TIME',start,((end-start)/nblines)+1,nblines)
                                # print end-start,data.shape[2]
                                array_posn.append([0,data.shape[2]-1,i])
                                # print old_flags.shape
                                # print 'CHECK1'
                                # print data.shape
                        elif i==len(list_of_starts):
                            start = select_arr[list_of_starts[i-1]+1]
                            end=select_arr[-1]
                            patht = 2
                            if (antdict[ant1] == 'Lo') and (antdict[ant2] == 'Mk2'):
                                array_posn.append([prev_end,((end-start)/nblines),i])
                                prev_end = ((end-start)/nblines)+1
                                lovelltotal = select_arr.shape[0]
                            elif (antdict[ant1] == 'Mk2') and (antdict[ant2] == 'Lo'):
                                array_posn.append([prev_end,((end-start)/nblines),i])
                                prev_end = ((end-start)/nblines)+1
                                lovelltotal = select_arr.shape[0]
                            else:
                                pre_arr_end = data.shape[2]
                                data = np.append(data,tb.getcol('DATA',start,((end-start)/nblines)+1,nblines),axis=2)
                                old_flags = np.append(old_flags,tb.getcol('FLAG',start,((end-start)/nblines)+1,nblines),axis=2)
                                times_array = np.append(times_array,tb.getcol('TIME',start,((end-start)/nblines)+1,nblines),axis=2)
                                array_posn.append([pre_arr_end,data.shape[2]-1,i])
                                # print end-start,data.shape[2]
                                # print data.shape
                        else:
                            start=select_arr[list_of_starts[i-1]+1]
                            end = select_arr[list_of_starts[i]]
                            patht=3
                            if (antdict[ant1] == 'Lo') and (antdict[ant2] == 'Mk2'):
                                array_posn.append([0,((end-start)/nblines),i])
                                prev_end = ((end-start)/nblines)+1
                                lovelltotal = select_arr.shape[0]
                            elif (antdict[ant1] == 'Mk2') and (antdict[ant2] == 'Lo'):
                                array_posn.append([0,((end-start)/nblines),i])
                                prev_end = ((end-start)/nblines)+1
                                lovelltotal = select_arr.shape[0]
                            else:
                                pre_arr_end = data.shape[2]
                                data = np.append(data,tb.getcol('DATA',start,((end-start)/nblines)+1,nblines),axis=2)
                                old_flags = np.append(old_flags,tb.getcol('FLAG',start,((end-start)/nblines)+1,nblines),axis=2)
                                times_array = np.append(times_array,tb.getcol('TIME',start,((end-start)/nblines)+1,nblines),axis=2)
                                array_posn.append([pre_arr_end,data.shape[2]-1,i])
                                # print end-start,data.shape[2]
                                # print data.shape
                        # print patht,bline,did,'if',start,((end-start)/nblines)+1,len(select_arr[start:end]),end-start,end
                        # arrdata = tb.getcol('DATA',start,((end-start)/nblines)+1,nblines)
                        # data.extend(tb.getcol('DATA',start,((end-start)/nblines)+1,nblines))
                        # print arrdata.shape
                        # data.extend(arrdata)
                        data=np.array(data)
                        old_flags=np.array(old_flags)
                        # print old_flags
                        # print end,data[0,0,-1]
                        # print data.shape
                        did_arr[ifp] = array_posn
                        # print len(array_posn),ifp,len(did_arr)
                else:
                    data = tb.getcol('DATA',select_arr[0],len(select_arr)+1,nblines)
                    old_flags = tb.getcol('FLAG',select_arr[0],len(select_arr)+1,nblines)
                    # print old_flags
                    times_array = tb.getcol('TIME',select_arr[0],len(select_arr)+1,nblines)
                    # print data.shape
                    # print 'else',select_arr[0],select_arr[1],(select_arr[-1]-select_arr[0])/nblines
                #data=np.array(data)
                #testarr = np.array(data)
                #print len(data),len(data[0][0][0]),select_arr.shape,testarr.shape
                # print select_arr.shape,testarr.shape
                # print tflags.nbytes,tflags.shape
                # print times.nbytes,times.shape,times[-1]
                # print old_flags.shape
                visoutfile = open(path2folder + "data__"+ str(source_list[srce]) + '__ID' + str(did) + '__' + str(bline)+'.npy',"wb")
                np.save(visoutfile,np.absolute(data))
                visoutfile.close()
                flagoutfile = open(oldflags_dir + "oldflags__"+ str(source_list[srce]) + '__ID' + str(did) + '__' + str(bline)+'.npy',"wb")
                np.save(flagoutfile,old_flags)
                flagoutfile.close()
                timesoutfile = open(path2folder + "times__"+ str(source_list[srce]) + '__ID' + str(did) + '__' + str(bline)+'.npy',"wb")
                np.save(timesoutfile,times_array)
                timesoutfile.close()

            bline_posn[bline] = did_arr
        srce_arr_posn[sdi] = bline_posn

        # print atop
          
##################################



# print len(srce_arr_posn)



'''
## Checks for Linux-based OS then runs memory definition
try:
    mem_stats = computer_memory(path2folder)
except:
    print "\n Sorry, computer_memory() definition does not work on your system!"
    
    
if 'mem_total' in mem_stats:
    print "System Memory Information:"
    print "Total Memory  :    %s" % array_size(mem_stats['mem_total']*1000)
    print "Used Memory   :    %s" % array_size(mem_stats['mem_used']*1000)
    print "Free Memory   :    %s" % array_size(mem_stats['mem_free']*1000)
    print "Total Swap    :    %s" % array_size(mem_stats['swap_total']*1000)
    print "Used Swap     :    %s" % array_size(mem_stats['swap_used']*1000)
    print "Free Swap     :    %s" % array_size(mem_stats['swap_free']*1000)


## The predicted array sizes, these predict the arrays when flagging all baselines 
## simultaneously
pred_array = 0
# print 'source_list is:', source_list
for srce in source_list:
    sdi = source_lookup.index(source_list[srce])
    # print sum(baseline_dict[sdi].itervalues())
    pred_array += sum(baseline_dict[sdi].itervalues())*nif*nchan*npol*3*8
#pred_array = nsource*nif*nchan*ntime*npol*3*8
print "Total predicted memory usage: %s" % array_size(int(pred_array))
'''





###############################################################################
#
#   Flagging algorithm split into baselines and IFs (for each baseline work 
#   through the IFs individually then move on to the next baseline)
#   This now uses the process/queue approach to multiprocessing from the 
#   multiprocessing python module.
#
###############################################################################

# Begin the parallelization block:

total_tscans = 0
pickle_list = {}
lovell_bline = []
zeros_bline = []

# nbline = 15
nsource=len(source_list)
#nbline = len(baselines)    ## need to change
total_jobs = nblines*nif*nsource



### Create array of jobs to do ####
jobcount = 1
joblist = []
for srcething in source_list:
    sdi = source_lookup.index(source_list[srcething])
    source_num = srcething
    source = source_list[srcething]
    srcname = source_list[srcething]
    for i in xrange(len(sbline[sdi])):
        # start = start_rownums[i]
        bline = sbline[sdi][i]
        ant1 = int(bline.split('-')[0])
        ant2 = int(bline.split('-')[1])
        # ntimescans = baseline_dict[sdi][bline]
        # rownums = drownums[sdi][i]
        if (antdict[ant1] == 'Lo') and (antdict[ant2] == 'Mk2'):
            print 'Bulk flagging Lo-Mk2'
            break
        elif (antdict[ant1] == 'Mk2') and (antdict[ant2] == 'Lo'):
            print 'Bulk flagging Mk2-Lo'
            break

        elif dobline_specific == 'yes' or dosource_specific == 'yes':
            parameters = source_base_parms[sdi][bline]
        else:
            parameters = source_base_parms[sdi]

        for j in xrange(len(descids)):
            nchan = chanlist[j]
            did = descids[j]
            # oldtimes = times[start::nblines]
            # times_array = np.array(oldtimes[np.where((sources[start::nblines]==srce) & (spws[start::nblines]==did))[0]])
            #print bline,j,source,nchan,ntimescans,polnames,path2folder,experiment,name,phasecal,lovell_bline,telescope,srcname,zero_level,do_lovell_cross,coadd_polarization_flags,pickle_list,flagging_options,parameters,uvdata.name,uvdata.klass,uvdata.disk,uvdata.seq,singlesource, source_num
            # joblist.append([bline,j,source,nchan,ntimescans,polnames,experiment,name,phasecal,lovell_bline,telescope,srcname,zero_level,do_lovell_cross,coadd_polarization_flags,pickle_list,ms,singlesource,source_num,orderedpols, drownums])
            

            joblist.append([bline,descids[j],source,nchan,polnames,msname,orderedpols,bpolnames,jobcount,total_jobs,path2folder,phasecal,lovell_bline,telescope,srcname,zero_level,do_lovell_cross,coadd_polarization_flags,flagging_options,parameters,source_num, oldflags_dir,antdict,do_sumthreshold])
            # print bline,descids[j],source, nchan,polnames,orderedpols
            #joblist.append([bline,j,nchan,ntimescans])
            jobcount += 1



########################################################################
### Parallelisation block for full flagger process, calls function 
### which runs lovell dropouts, zero-level dropouts and flagger on 
### each source, baseline and IF separately
########################################################################


print '\n Beginning flagging processes...\n'

send_q = multip.Queue()
rec_q = multip.Queue()
ncpus = multip.cpu_count()

for i in xrange(len(joblist)):
    print 'Sending to queue job', i+1, 'of', len(joblist)
    # print "CPU", c, "will flag baselines and IFs: \n"
    send_q.put(joblist[i])


for i in xrange(ncpus):
    proc = multip.Process(target=flagger_process, args=(send_q,rec_q,i))
    proc.start()
    print 'Starting process on CPU:', i






results = []
for i in xrange(len(joblist)):
    results.append(rec_q.get())

for i in xrange(ncpus):
    send_q.put('STOP')


results=np.array(results)
# print results.shape

#np.reshape(results,(2,len(results)))
#print results.shape

'''
pname = str(name) +"__"+ str(source)  + "__IF" + str(nif+1) + "__" + str(bline) + "__" + pol
pickle.dump(flag_LR, open(path2folder+str(pname)+".p", "wb"))
pickle_list[pname] = [source, nif+1, bline, times_array]
pickle.dump(times_array, open(path2folder+str(pname)+".info", "wb"))
'''

# coadd_polarization_flags = 'yes'
# zero_flag_allif = True


# for pol in polnames:
    # print pol


# print stop


npol = len(corrtypes)
# print npol

try:
    time_limit = flag_coinc_times
except NameError:
    time_limit = 0

try:
    count_limit = flag_coinc_chans
except NameError:
    count_limit = 0
   




print "\n Loading new flags to data....\n"
preloadtime = time.time()
for srce in source_list:
    sdi = source_lookup.index(source_list[srce])
    nblines = len(sbline[sdi])
    for bline in sbline[sdi]:
        ant1 = int(bline.split('-')[0])
        ant2 = int(bline.split('-')[1])
        lovellextras = [[] for i in range(npol)]
        zeroextras = [[] for i in range(npol)]
        if ((do_Lovell_scan == 'yes') & ((antdict[ant1] == 'Lo') or (antdict[ant2] == 'Lo'))) or (zero_flag_allif == 'yes'):
            print 'Loading dropouts across spws'
            for did in descids:
                for pol in polnames:
                    stoke = polnames[pol]
                    plname = str(name) + "__" + str(source_list[srce]) + "__" + str(bline) + "__IF" + str(did+1)+ "__" + pol + "__zeros_dummy"
                    if os.path.exists(path2folder+str(plname)+'.zeros') == True:
                        # print plname
                        zeroextras[stoke] = zeroextras[stoke] + pickle.load(open(path2folder + str(plname) + ".zeros", "rb"))
                        os.remove(path2folder+str(plname)+'.zeros')
                    plname = str(name) + "__" + str(source_list[srce]) + "__" + str(bline) + "__IF" + str(did+1)+ "__" + pol + "__inscan_zeros_dummy"
                    if os.path.exists(path2folder+str(plname)+'.zeros') == True:
                        # print plname
                        zeroextras[stoke] = zeroextras[stoke] + pickle.load(open(path2folder + str(plname) + ".zeros", "rb"))
                        os.remove(path2folder+str(plname)+'.zeros')
                    plname = str(name) + "__" + str(source_list[srce]) + "__" + str(bline) + "__IF" + str(did+1)+ "__" + pol + "__lovell_dummy"
                    if os.path.exists(path2folder+str(plname)+'.dropouts') == True:
                        # print plname
                        lovellextras[stoke] = lovellextras[stoke] + pickle.load(open(path2folder + str(plname) + ".dropouts", "rb"))
                        os.remove(path2folder+str(plname)+'.dropouts')
        for did in descids:
            nchan = chanlist[descids.index(did)]
            # Reload old_flags temp file as contains all stokes including those not flagged
            if (antdict[ant1] == 'Lo') and (antdict[ant2] == 'Mk2'):
                nchan = chanlist[descids.index(did)]
                new_flags = np.ones((npol,nchans,lovelltotal),int)
            elif (antdict[ant1] == 'Mk2') and (antdict[ant2] == 'Lo'):
                new_flags = np.ones((npol,nchans,lovelltotal),int)
                nchan = chanlist[descids.index(did)]
            else:
                flagsfile = oldflags_dir + "oldflags__"+ str(source_list[srce]) + '__ID' + str(did) + '__' + str(bline)+'.npy'
                new_flags = np.load(flagsfile).astype(int)
                os.remove(flagsfile)
            # print '#orig'
            # print new_flags[:,:,0]
            new_flags = np.where(new_flags!=1,new_flags,-15)
            # print '#15s'
            # print new_flags[:,:,0]
            # allpols = [[] for i in range(len(polnames))]
            # newpolflag = [[] for i in range(len(polnames))]
            for pol in polnames:
                stoke = polnames[pol]
                # print name, source_list[source], j+1,bline,pol
                pname = str(name)  +"__"+ str(source_list[srce]) + '__IF' + str(did+1) + '__' + str(bline) +  '__' + pol
                # print "Reading file:", pname
                # print path2folder+str(pname)+'.p'
                if os.path.exists(path2folder+str(pname)+'.p') == True:
                    print "Reading file:", pname
                    polflag = np.swapaxes(pickle.load(open(path2folder + str(pname) +".p", "rb")),0,1)
                    # new_flags[stoke] = polflag
                    if ((do_Lovell_scan == 'yes') & ((antdict[ant1] == 'Lo') or (antdict[ant2] == 'Lo'))):
                        polflag[:,lovellextras[stoke]] = -7
                    if (zero_flag_allif == 'yes'):
                        polflag[:,zeroextras[stoke]] = -6
                    # print '#polflag'
                    polflag = np.where(polflag!=-6,polflag,6)
                    polflag = np.where(polflag!=-7,polflag,7)
                    # print polflag.shape, new_flags[stoke].shape
                    # print polflag[:,0]
                    new_flags[stoke] = np.select([new_flags[stoke]<0,polflag>0],[new_flags[stoke],polflag])
                    # print '#select'
                    # print np.select([new_flags[stoke]<0,polflag>0],[new_flags[stoke],polflag])[:,0]
                    # print polflag.shape
                    if flag_all_stokes == 'yes':
                        # print new_flags[:,:,0]
                        # print np.logical_and(new_flags[:]<=0,polflag<=0)[:,:,0]
                        new_flags[:] = np.where(np.logical_and(new_flags[:]!=-1,polflag<=0),new_flags[:],-1)
                        new_flags[:] = np.where(np.logical_and(new_flags[:]!=6,polflag!=6),new_flags[:],-6) #zeros
                        new_flags[:] = np.where(np.logical_and(new_flags[:]!=7,polflag!=7),new_flags[:],-7) #Lovell
                        # new_flags[:][np.logical_or(new_flags[:]>0,polflag>0)]=1
                        # print new_flags[:,:,0]
                    else:
                        if coadd_polarization_flags == 'yes':
                            new_flags[bpolnames.keys()] = np.where(np.logical_and(new_flags[bpolnames.keys()]!=-1,polflag<=0),new_flags[bpolnames.keys()],-1) #flagger
                            new_flags[bpolnames.keys()] = np.where(np.logical_and(new_flags[bpolnames.keys()]!=6,polflag!=6),new_flags[bpolnames.keys()],-6) #zeros
                            new_flags[bpolnames.keys()] = np.where(np.logical_and(new_flags[bpolnames.keys()]!=7,polflag!=7),new_flags[bpolnames.keys()],-7) #Lovell
                        else:
                            if zero_level == 'yes' and coadd_zero_flags == 'yes':
                                new_flags[bpolnames.keys()] = np.where(np.logical_and(new_flags[bpolnames.keys()]!=6,polflag!=6),new_flags[bpolnames.keys()],-6) #zeros
                                # new_flags[bpolnames.keys()]= np.where(np.logical_or(new_flags[bpolnames.keys()]==-6,polflag==-6),new_flags[bpolnames.keys()],1) #zeros
                            if do_lovell_scan == 'yes' and coadd_lovell_flags == 'yes':
                                new_flags[bpolnames.keys()] = np.where(np.logical_and(new_flags[bpolnames.keys()]!=7,polflag!=7),new_flags[bpolnames.keys()],-7) #Lovell
                                # new_flags[bpolnames.keys()]= np.where(np.logical_or(new_flags[bpolnames.keys()]==-7,polflag==-7),new_flags[bpolnames.keys()],1) #Lovell
                    os.remove(path2folder+str(pname)+'.p')
                    os.remove(path2folder+str(pname)+'.info')
            new_flags = np.where(new_flags!=-15,new_flags,1)
            new_flags = np.where(new_flags!=-1,new_flags,1)
            new_flags = np.where(new_flags!=-6,new_flags,1)
            new_flags = np.where(new_flags!=-7,new_flags,1)
            if count_limit > 0 or time_limit > 0:
                print '\nExtra flagging requested'
                print 'Flagging every %i coincident channels and every %i coincident times\n' % (count_limit, time_limit)
                for pol in polnames:
                    stoke = polnames[pol]
                    new_flags[stoke] = coinc_flags_noncumu(new_flags[stoke],time_limit, count_limit)
            # print '#final'
            # print new_flags[:,:,0]
            new_flags=np.where(new_flags==0,new_flags,1).astype(bool)
            select_arr = np.where((sources==srce) & (spws==did) & (ants1[:]==ant1) & (ants2[:] == ant2))[0]
            testdiffs = np.diff(select_arr)
            if (testdiffs>nblines).any():
                # print 'Multiple scans detected'
                list_of_starts = np.where(testdiffs>nblines)[0]
                for i in xrange(len(list_of_starts)+1):
                    if i == 0:
                        start=select_arr[0]
                        end = select_arr[list_of_starts[i]]
                        patht = 1
                        arr_start = srce_arr_posn[sdi][bline][did][i][0]
                        arr_end = srce_arr_posn[sdi][bline][did][i][1]
                        # old_flags = tb.getcol('FLAG',start,((end-start)/nblines)+1,nblines)
                        # print start,end,nblines,new_flags.shape,arr_start,arr_end
                        # print start,((end-start)/nblines)+1,nblines,new_flags[:,:,arr_start:arr_end+1].shape
                        tb.putcol('FLAG',new_flags[:,:,arr_start:arr_end+1],start,((end-start)/nblines)+1,nblines)
                    elif i==len(list_of_starts):
                        start = select_arr[list_of_starts[i-1]+1]
                        end=select_arr[-1]
                        patht = 2
                        arr_start = srce_arr_posn[sdi][bline][did][i][0]
                        arr_end = srce_arr_posn[sdi][bline][did][i][1]
                        # old_flags = np.append(new_flags,tb.getcol('FLAG',start,((end-start)/nblines)+1,nblines),axis=2)
                        tb.putcol('FLAG',new_flags[:,:,arr_start:arr_end+1],start,((end-start)/nblines)+1,nblines)
                    else:
                        start=select_arr[list_of_starts[i-1]+1]
                        end = select_arr[list_of_starts[i]]
                        patht=3
                        arr_start = srce_arr_posn[sdi][bline][did][i][0]
                        arr_end = srce_arr_posn[sdi][bline][did][i][1]
                        # old_flags = np.append(old_flags,tb.getcol('FLAG',start,((end-start)/nblines)+1,nblines),axis=2)
                        tb.putcol('FLAG',new_flags[:,:,arr_start:arr_end+1],start,((end-start)/nblines)+1,nblines)
            else:
                old_flags = tb.getcol('FLAG',select_arr[0],len(select_arr)+1,nblines)
                tb.putcol('FLAG',new_flags,start,select_arr[0],len(select_arr)+1,nblines)



ms.close()
print ' Time taken to load flags %s' % time2hms(time.time()-preloadtime)


# print "\nTime taken to append visibilities (hh:mm:ss):", time2hms(time.time()-ti)
# Calculate time taken and write log file
print '\n---------------------------------------------\n'
print 'Finished on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),''
print "Time taken (hh:mm:ss):", time2hms(time.time()-ti),"\n"
if nodata_list:
    print 'Source(s) %s not flagged, no data found' % ', '.join(nodata_list)
    
print '\n---------------------------------------------\n'


