# =============================================================================
#              meclas.py
# Package for running various laser utilities in MEC
# Apologies on behalf of: Eric Cunningham (and others)
#
# To load: use import meclas or use IPython's %run magic function
#
# Class list and brief description:
#    LPL -- routines for LPL pulse shaping (with some aux functions), data acquisition, etc.
#    efc -- extra function class, holds many useful utilities and shortcuts
#    LOSC -- LeCroy oscilloscope trace read-out, plotting, saving, etc.
#    EMeters -- LPL and SPL energy meters
#    MBC -- LPL bias controller utilities
#    YFE -- LPL YLF Front End seed laser utilities
#    PFN -- LPL Pulse Forming Network cap bank charging utilities
#    HWP -- LPL Half Wave Plate motor utilities
#    **Stage -- Newport and SmarAct stage utilities
#    **Timing -- ns and fs timing utilities
#    CAM -- functions for GigE camera acquisition, configuration, etc.
#    TTL_shutter -- Beckhoff utilities for controlling/tracking Thorlabs TTL shutters
#    DG645 -- functions for DG645 operation, parameter backup and restoration, etc.
#    **SPL -- routines for SPL alignment, etc. 
#    UNIBLITZ -- UNIBLITZ shutter utilities for setting SPL trigger modes, etc.
#    **Spectrometer -- functions for Qmini and Ocean Optics USB4000 spectrometers
#    **VISAR -- routines for VISAR timing, streak camera configuration, laser control, etc.
#    CtrlSys -- routines for checking responsivity of PVs, hosts, hutch computers, etc.
#    **SCALLOPS -- routines for LPL pulse shaping simulations
#    **LabEnv -- functions for interfacing with lab environment monitors
#    **RIS -- functions for monitoring RIS-related systems and PVs
#    **PDU -- functions for operating power distribution units
#    GLOBAL -- home for global constants, PV definitions, etc.
#    
#    ** = FUTURE DEVELOPMENT
# =============================================================================


### load packages
import socket
import time
import math
import numpy as np
import struct
#from scipy import signal
#from scipy import stats
import matplotlib.pyplot as plt
import pickle
from datetime import date, datetime
from binascii import hexlify
from binascii import unhexlify
import csv
import os.path
import sys
from ophyd.signal import EpicsSignal
import elog
import pcdsdaq.ext_scripts 
import glob
#import pandas as pd
import stat
import getpass
import multiprocessing
#import threading
import termios, tty
import importlib.util
import select
import regex as re



class LPL:
    """
    Stores functions related to LPL pulse shaping, data acquisition, etc.
    """
    def LinearWave(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height):
        """
        Generates linearly-interpolated waveform of length 140 between two input points
        Primarily intended for use with the Highland AWG, so max value capped at 65535
        
        LinearWave(1,10000,120,28000) returns a waveform of length 140 with linear
            ramp from (pixel 1, height 10000) to (pixel 120, height 28000)
        """
        itt=0
        NewString=''
        if Edge1Height>65535:
            print('Edge1 height exceeds max value of 65535')
            h1=65535
        elif Edge1Height<0:
            print('Edge1 height must be positive')
            h1=0
        else:
            h1=int(Edge1Height)
        if Edge2Height>65535:
            print('Edge2 height exceeds max value of 65535')
            h2=65535
        elif Edge2Height<0:
            print('Edge2 height must be positive')
            h2=0
        else:
            h2=int(Edge2Height)
        #
        if Edge1PixNo>Edge2PixNo:
            print('Edge1 must come before Edge2')
            Dummy=int(Edge1PixNo)
            Edge1PixNo=int(Edge2PixNo)
            Edge2PixNo=Dummy
        if Edge1PixNo<1:
            print('Edge1 pixel number must be >=1')
            p1=0
        elif Edge1PixNo>140:
            print('Edge1 pixel number must be <=140')
            p1=139
        else:
            p1=int(Edge1PixNo)-1
        if Edge2PixNo<1:
            print('Edge2 pixel number must be >=1')
            p2=0
        elif Edge2PixNo>140:
            print('Edge2 pixel number must be <=140')
            p2=139
        else:
            p2=int(Edge2PixNo)-1
        #
        if p1==p2:
            print('Warning: pulse width specified as single pixel.')
            return 140*[0]
        #
        while itt<140:
            if itt<p1:
                NewString+='0000'
            elif p1<=itt<=p2:
                NewString+=HAWG._Hex2Byte(int(h2+((itt-p2)*(h2-h1)/float(p2-p1))))
            else:
                NewString+='0000'
            itt+=1
        return NewString


    def LinearWave2(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height,offsetQ,arraylenQ):
        """
        Generates linearly-interpolated waveform between two points (and 0 outside those points)
        with requested linear offset and array length
        
        Useful in part for specifying YFE waveforms at 10Hz, generating "goal" waveforms, etc.
        
        LinearWave2(500,.1,1025,.8,0,5002) returns an array of length 5002 with a linear ramp
            from (pixel 500, height 0.1) to (pixel 1025, height 0.8) with no vertical offset
        """
        itt=0
        h1=Edge1Height-offsetQ
        h2=Edge2Height-offsetQ
        p1=int(Edge1PixNo)-1
        p2=int(Edge2PixNo)-1
        NewList=[]
            #
        while itt<arraylenQ:
            if itt<p1:
                NewList.append(offsetQ-offsetQ)
            elif p1<=itt<=p2:
                nextval=h2+((itt-p2)*(h2-h1)/float(p2-p1))#h1*((h2/float(h1))**((itt-p1)/float(p2-p1)))
                NewList.append(nextval)
            else:# itt>p2:
                NewList.append(offsetQ-offsetQ)
            itt+=1
        return np.array(NewList)+offsetQ

    def ParabolicWave2(Edge1PixNo,Edge1Height,MidPixNo,MidHeight,Edge2PixNo,Edge2Height,offsetQ,arraylenQ):
        """
        Generates parabolically-interpolated waveform using three points (and 0 outside those points)
        with requested linear offset and array length
        
        Only rarely used for specifying YFE waveforms at 10Hz, generating "goal" waveforms, etc.
        
        ParabolicWave2(500,.1,800,.15,1025,.8,0,5002) returns an array of length 5002 with a
            parabola fit to (pixel 500, height 0.1), (pixel 800, height 0.15), and 
            (pixel 1025, height 0.8) with no vertical offset
        """
        itt=0
        h1=Edge1Height-offsetQ
        h2=Edge2Height-offsetQ
        h3=MidHeight-offsetQ
        p1=int(Edge1PixNo)-1
        p2=int(Edge2PixNo)-1
        p3=int(MidPixNo)-1
        NewList=[]

        while itt<arraylenQ:
            if itt<p1:
                NewList.append(offsetQ-offsetQ)
            elif p1<=itt<=p2:
                nextval=(h1*(itt-p2)*(itt-p3)/float((p2-p1)*(p3-p1)))+(h2*(itt-p1)*(itt-p3)/float((p2-p1)*(p2-p3)))+(h3*(itt-p1)*(itt-p2)/float((p3-p1)*(p3-p2)))
                NewList.append(nextval)
            else:# itt>p2:
                NewList.append(offsetQ-offsetQ)
            itt+=1
        return np.array(NewList)+offsetQ

    def ExponentialWave2(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height,offsetQ,arraylenQ):
        """
        Generates exponentially-interpolated waveform using two points (and 0 outside those points)
        with requested linear offset and array length
        
        Most-used function for specifying YFE waveforms at 10Hz (exponential seed waveforms
            become ~linear after laser amplification)

        ExponentialWave2(500,.1,1025,.8,0,5002) returns an array of length 5002 with an exponential ramp
            from (pixel 500, height 0.1) to (pixel 1025, height 0.8) with no vertical offset
        """
        itt=0
        h1=Edge1Height-offsetQ
        h2=Edge2Height-offsetQ
        p1=int(Edge1PixNo)-1
        p2=int(Edge2PixNo)-1
        NewList=[]
            #
        while itt<arraylenQ:
            if itt<Edge1PixNo:
                NewList.append(offsetQ-offsetQ)
            elif p1<=itt<=p2:
                nextval=h1*((h2/float(h1))**((itt-p1)/float(p2-p1)))
                NewList.append(nextval)
            else:# itt>p2:
                NewList.append(offsetQ-offsetQ)
            itt+=1
        return np.array(NewList)+offsetQ

    def LogWave2(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height,LogBase,offsetQ,arraylenQ):
        """
        Generates logrithmically-interpolated waveform using two points (and 0 outside those points)
        with requested log base, linear offset, and array length
        
        Essentially never used but left in the code just in case it finds use later
        """
        itt=0
        h1=Edge1Height-offsetQ
        h2=Edge2Height-offsetQ
        p1=int(Edge1PixNo)-1
        p2=int(Edge2PixNo)-1
        NewList=[]
        #
        while itt<140:
            if itt<p1:
                NewList.append(offsetQ-offsetQ)
            elif p1<=itt<=p2:
                NewList.append(((h1-h2)/math.log(((p1+.0001)/float(p2)),LogBase))*math.log((LogBase**(((h2*math.log(p1+.0001,LogBase))-(h1*math.log(p2+.0001,LogBase)))/float(h1-h2)))*(itt+.0001),LogBase))
            else:# itt>p2:
                NewList.append(offsetQ-offsetQ)
            itt+=1
        return np.array(NewList)+offsetQ

    def ComboWave(WList): #accept list or csv of 140 pts
        """
        Combines list of 140-pt arrays into single waveform scaled to have a maximum value of 1
        """
        PreNormL=[]
        for DesiredOutputPulseShapeQ in WList:
            if len(DesiredOutputPulseShapeQ) == 140*4:#will accept pre-formatted Hex2Byte text
                PreNormL.append(np.array([int(DesiredOutputPulseShapeQ[4*ii:4*ii+4],16) for ii in range(len(DesiredOutputPulseShapeQ)//4)]))
            elif len(DesiredOutputPulseShapeQ)==140:#will accept a straight list
                PreNormL.append(np.array(DesiredOutputPulseShapeQ))
            elif DesiredOutputPulseShapeQ.endswith(('.txt','.csv','.dat')):#will accept file
                with open(DesiredOutputPulseShapeQ,'r') as filehead:
                    RawListQ=filehead.read()
                    if '\r\n' in RawListQ:
                        ListedValues=RawListQ.split('\r\n')
                    elif '\n' in RawListQ:
                        ListedValues=RawListQ.split('\n')
                    elif ',' in RawListQ:
                        ListedValues=RawListQ.split(',')
                    else:
                        print('Unrecognized format on input file.')
                        return
                if len(ListedValues) != 140:
                    print('File must have 140 entries; entry count: '+str(len(ListedValues)))
                    return
                PreNormL.append(np.array(ListedValues))
        CPreNormL=np.sum(PreNormL,0)
        return [1.*entry/float(max(CPreNormL)) for entry in CPreNormL]

    def TraceFormatting(PDTrace, PDFETMap, MaxPDValue, AvgRange=25, FWHM=4):
        """
        Takes PD trace and mapping function to generate windowed, averaged, scaled list of 140 pts
        """
        MeasuredOutputPulseShape=[]
        WeightList=[math.exp(-4*math.log(2)*((ii+1-round(AvgRange/2.))/FWHM)**2) for ii in range(AvgRange)]
        WSum = sum(WeightList)
        MX , B = PDFETMap
        for FETNo in range(140):
            Loc = round(MX*FETNo + B)
            WSample=sum([PDTrace[int(Loc+(ii+1-round(AvgRange/2.)))]*WeightList[ii] for ii in range(AvgRange)])/WSum
            MeasuredOutputPulseShape+=[WSample/MaxPDValue]
        return MeasuredOutputPulseShape 

    def UpdatingShapingAlgorithm(DesiredOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape, StepQ):
        """
        Accepts pre-formatted input, measurement, and goal waveforms to calculate next-iteration input waveform using specified step size
        """
        G, M, I = DesiredOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape
        NewInputPulseShape=np.clip([abs((StepQ*(G[ii]-M[ii]))+I[ii])*math.ceil(G[ii]) for ii in range(len(G))],0,1)#math.ceil(G) is a mask that disallows values outside the goal
        return NewInputPulseShape

    def FixEdges(WavF,DurationListQ,StartStopListQ,CorrFactorFront=1.02,CorrFactorBack=1.02,PtNumFront=1,PtNumBack=2):
        """
        Applies fixed relationship between points close to waveform discontinuities (e.g. at beginning, end, step, etc.)
        specified by Duration and Start/Stop lists
        """
        FirstPix=0#was 50
        DurListQ=np.cumsum([0]+DurationListQ)
        fWavF=WavF[:]
        for ii in range(PtNumFront):
            fWavF[FirstPix+(PtNumFront-ii-1)]=fWavF[FirstPix+(PtNumFront-ii)]*CorrFactorFront
        DisconCount=0
        ContCount=0
        if len(StartStopListQ)>1:
            for ii in range(len(StartStopListQ)-1):
                if StartStopListQ[ii][-1] != StartStopListQ[ii+1][0]:
                    DisconCount+=1
                    fWavF[FirstPix+int(4*DurListQ[ii+1])-1-ContCount]=fWavF[FirstPix+int(4*DurListQ[ii+1])-2-ContCount]
                    fWavF[FirstPix+int(4*DurListQ[ii+1])-ContCount]=fWavF[FirstPix+int(4*DurListQ[ii+1])+1-ContCount]
                else:
                    ContCount+=1
        #fWavF[FirstPix]=fWavF[FirstPix+1]*1.1
        #try fixing last THREE pixels to help back edge
        for ii in range(PtNumBack):
            fWavF[FirstPix+int(4*DurListQ[-1])-(PtNumBack-1-ii)-ContCount]=fWavF[FirstPix+int(4*DurListQ[-1])-(PtNumBack-ii)-ContCount]*CorrFactorBack
        for ii in range(len(fWavF)):
            if np.mean(fWavF) < 1:
                if fWavF[ii] > 1:
                    fWavF[ii]=1.0
            else:
                if fWavF[ii]>28000:
                    fWavF[ii]=28000
        return fWavF
    
    def smooth_wvfm(wvfm_in):
        """
        Performs rudimentary smoothing of waveform by looking at neighboring pixels
        """
        wvfm_out=wvfm_in[:]
        for ii in range(len(wvfm_in)-2):
            wvfm_out[ii+1]=.25*wvfm_in[ii]+.5*wvfm_in[ii+1]+.25*wvfm_in[ii+2]
        return wvfm_out

    @classmethod
    def PulseGoal(cls,DurationListQ,StartStopListQ):#140 pt list with max at 1
        """
        Generates 140-pt list according to provided Duration and Start/Stop lists
        """
        BeginPix=1#was 51
        #PulseGoal([3,5],[[15,25],[75,100]],['l','l'])
        #PulseGoal([3,5],[[15,25],[75,100]],['l','p'])
        DurListQ=np.cumsum([0]+DurationListQ)
        SSListQ=StartStopListQ[:]
        SegTotQ=len(DurListQ)-1
        if SegTotQ!=len(SSListQ):
            print('Error')
            return
        DelListQ=[]
        for ii in range(-1+len(SSListQ)):
            if len(SSListQ[ii])!=2:
                print('Error')
                return
            if SSListQ[ii][1]==SSListQ[ii+1][0]:
                DelListQ.append(BeginPix+(DurListQ[ii+1]*4)-1)
        SegmentsQ=[]
        for ii in range(SegTotQ):
            #SegmentsQ.append(LinearWave(51+(DurListQ[ii]*4),int(20000.*SSListQ[ii][0]/100.),51+(DurListQ[ii+1]*4),int(20000.*SSListQ[ii][1]/100.)))
            SegmentsQ.append(cls.LinearWave(int(BeginPix+(DurListQ[ii]*4)),int(20000.*SSListQ[ii][0]/100.),int(BeginPix+(DurListQ[ii+1]*4)-1),int(20000.*SSListQ[ii][1]/100.)))
        return np.append(np.delete(np.array(cls.ComboWave(SegmentsQ)),np.array(DelListQ).astype(int)),[0]*len(DelListQ))
    
    @classmethod
    def PulseGoal10Hz(cls):
        """One day this will make it easier to use psefc10Hz and psrefrwvfm by making the specification of 'pwt=' less cryptic"""
        pass
    
    @classmethod
    def SmartPulseGoal10Hz(cls):
        """One day this will help determine the 10Hz PulseGoal using input from the SCALLOPS class"""
        pass

    @classmethod
    def PulseMax(cls,DurationListQ,StartStopListQ,zzJQ):#fixed normalization part
        """Gets amplitude setting for segmented, arbitrary PulseGoal using Duration and Start/Stop lists and targeted energy"""
        return (1.*StartStopListQ[-1][-1]/100.)*(50.*zzJQ/(5.*500.*np.sum(cls.PulseGoal(DurationListQ,StartStopListQ))))

    @classmethod
    def wIter2(cls,sQ,wQ,DurationListQ,StartStopListQ,zzJQ,mapnowQ,stepQQ):
        """Calculates next suggested AWG input given 1) a previous full-energy waveform (+ mapping) and its corresponding AWG input,
        2) the Duration and Start/Stop lists to specify the goal, and 3) the requested step size of the correction"""
        DurListQ=np.cumsum([0]+DurationListQ)
        w1,w2=0,int(DurListQ[-1]*4)+5 # 50-5, 50+int(DurListQ[-1]*4)+5
        PGQ=cls.PulseGoal(DurationListQ,StartStopListQ)
        if np.abs(len(sQ)-10000)<10:
            PMQcorr=1
        elif np.abs(len(sQ)-1000)<10:
            PMQcorr=10
        else:
            print('Warning: unanticipated pulse shape array length of '+str(len(sQ))+', PulseMax scaling may be off...')
            PMQcorr=1
        PMQ=cls.PulseMax(DurationListQ,StartStopListQ,zzJQ)*PMQcorr
        wnew2=cls.FixEdges(cls.UpdatingShapingAlgorithm(PGQ,cls.TraceFormatting(sQ,mapnowQ,PMQ,AvgRange=1,FWHM=1), wQ,stepQQ),DurationListQ,StartStopListQ)
        #epll([wQ[w1:w2],wnew2[w1:w2],np.array(TraceFormatting2(sQ,mapnowQ,PMQ))[w1:w2]*.6,np.array(PGQ)[w1:w2]*.6])
        efc.epll([0.*np.array(wnew2[w1:w2]),np.array(wnew2[w1:w2])-np.array(wQ[w1:w2]),np.array(cls.TraceFormatting(sQ,mapnowQ,PMQ,AvgRange=1,FWHM=1))[w1:w2]*.6,np.array(PGQ)[w1:w2]*.6])
        return wnew2


    def weichall():
        """
        Generates weighted waveforms for YFE, 1in1w, 4x2in1w, and 4x2in2w outputs using scopes and energy meters 
        """
        try:
            LAchall=LOSC('a').rchall();LBchall=LOSC('b').rchall();L2chall=LOSC('2').rchall();
            allwvfm=[*LAchall[:2],*LBchall,*L2chall];
            allenergy=[*EMeters.EG1wYFE1in(),*EMeters.EG1w2in()[0],*EMeters.EG()[0][0]]
            allweich=[]
            for ii in range(10):
                templistQ=allwvfm[ii]
                bkgrdbuffer=int(0.038*len(templistQ))
                bkgrdQ=np.mean(templistQ[:bkgrdbuffer])
                ensampQ=allenergy[ii]
                weightQ=ensampQ/np.sum(np.array(templistQ)-bkgrdQ)
                allweich.append(np.array(weightQ*(np.array(templistQ)-bkgrdQ)))
            return allweich
        except:
            print('Weighted waveform generation failed!')
            return False

    def weichToPowerVsTime(weiarr):
        """Converts energy-weighted scope channels into a tuple of (time_values, instantaneous_power_in_watts)"""
        return (np.linspace(-5,45,len(weiarr)), np.array(weiarr)/(50e-9/len(weiarr)))
    
    def Deconv():
        """One day this will account for the instrument response of the PD, amp, scope, etc. in determining a detected waveform"""
        pass
        
    def pshostcheck():
        """This warns users if they're trying to do sg that requires the use specific hosts (e.g. to reach the ICS subnet)"""
        try:
            hostname=socket.gethostname()
            if hostname not in GLOBAL.OKHOSTS:
                print('Host must be one of the following:')
                for eahost in GLOBAL.OKHOSTS:
                    print(eahost)
                print('Current host: '+hostname)
                raise Exception
        except Exception:
            print('EXIT')
            os._exit(1)
        try:
            curruser=getpass.getuser()
            if curruser not in GLOBAL.OKUSERS:
                print('Warning: you are logged in as '+curruser+'. Beware of permissions issues... You may even unwittingly cause some yourself!')
                print('Suggested users:')
                for eauser in GLOBAL.OKUSERS:
                    print(eauser)
        except:
            print('Failed: could not ID current user!')
        return

    def LMap2():
        """"Pixel mapping from LeCroy2 horizontal axis (10002px) to Highland (140px)"""
        return GLOBAL.LMap2
    #            if len(stoday[-1]) == 10002:
    #            mapnow=LMap()
    #        elif len(stoday[-1]) == 1002:
    #            mapnow=[5,100]

    def LMapAB():
        """"Pixel mapping from LeCroyA and LeCroyB horizontal axis (1002px) to Highland (140px)"""
        return GLOBAL.LMapAB

    def DateString():
        """Shorthand way of getting the current date in format YYYYMMDD"""
        qdate=date.today()
        return qdate.strftime('%Y%m%d')

    def psfilepath():
        """Sets the file path for everything related to pulse shaping"""
        return GLOBAL.PSFILEPATH

    def get_curr_exp(timeout=15, hutch_name='mec'): 
        """Returns the name of the current experiment running in MEC"""
        script=pcdsdaq.ext_scripts.SCRIPTS.format(hutch_name,'get_curr_exp') 
        exp=pcdsdaq.ext_scripts.cache_script(script,timeout=timeout) 
        curr_exp=exp.lower().strip('\n')
        try:
            GLOBAL.CurrExp.put(curr_exp)
        except:
            GLOBAL.CurrExp.put(hutch_name+'xx####')
        finally:
            print('Failed to write current experiment to notepad PV!')
        return curr_exp
    
    def get_curr_run(timeout=15, hutch_name='mec'): 
        """Returns the current run number in MEC"""
        try:
            curr_run=pcdsdaq.ext_scripts.get_run_number(hutch=hutch_name,timeout=10)
        except:
            print('Failed to retrieve run number, setting to 9000!')
            curr_run=9000
        try:
            GLOBAL.CurrRun.put(int(curr_run))
        except:
            print('Failed to write run number to notepad PV!')
        return curr_run

    @classmethod
    def psheaders(cls):
        """Checks to see if shotlog files exist for the day, returns the date string and current run number"""
        cls.pshostcheck()
        DateStr=cls.DateString()
        for head in ['w','y','s1in1w','s42in1w','s42in2w','s']:
            if not os.path.exists(cls.psfilepath()+head+DateStr+'.p'):
                if head == 'w':
                    print('No laser file found -- probably the first shot of the day.')
                try:
                    efc.pickledump2([],cls.psfilepath()+head+DateStr+'.p')
                except:
                    print('Could not create file {}!'.format(cls.psfilepath()+head+DateStr+'.p'))
        curr_run=cls.get_curr_run()
        return (DateStr, curr_run)

    @classmethod
    def psacqx(cls, save_flag=True, RunNumQQ=False):
        """Acquisition sequence after shooting the LPL, primary component of pspreshot()"""
        (DateStr, curr_run) = cls.psheaders()
        psfpQ=cls.psfilepath()
        RunNumStr=str(curr_run).zfill(4)
        RunName='run'+str(RunNumStr)+'_'
        headlist=['AB','EF','GH','IJ']
        #get the current experiment name
        try:
            ExpName=cls.get_curr_exp()
        except:
            ExpName=GLOBAL.CurrExp.get()
        finally:
            print('Could not retrieve current experiment name!')
            ExpName='temp'
        #check for folders for saving scope data; if they don't exist, it makes them
        fpQ='/reg/neh/operator/mecopr/experiments/'+ExpName+'/lecroy/'
        if not os.path.exists(fpQ[-7]):
            print('File path '+fpQ[-7]+' does not exist! Trying to create it...')
            try:
                os.makedirs(fpQ[-7]);print('Folder created successfully!');
                os.chmod(fpQ[-7],stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO);
            except:
                print('Failed to create '+fpQ[-7]+'!')

        if not os.path.exists(fpQ):
            print('File path '+fpQ+' does not exist! Trying to create it...')
            try:
                os.makedirs(fpQ);print('Folder created successfully!');
                os.chmod(fpQ,stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO);
            except:
                print('Failed to create '+fpQ+'!')
        
        #check which heads are enabled
        for ii in range(4):
            if PFN.HeadENB()[ii]:
                RunName=RunName+headlist[ii]
                          
        weiarr=cls.weichall()#fix to give np.zeros(5002) when broken

        total_print=''
        #check so don't overwrite if someone forgets to change run number
        wtoday=pickle.load(open(psfpQ+'w'+DateStr+'.p','rb'))
        if os.path.exists(fpQ+RunName+'energies.txt'):
            RunName=str(RunName+'DUPLICATE')
            print(str('This run number already exists; run name '+RunName+' used'))
            total_print+=str('This run number already exists; run name '+RunName+' used')
            total_print+='\n'
        else:
            print(str('Run name: '+RunName+', shot number: '+str(len(wtoday)+1)))
            total_print+=str('Run name: '+RunName+', shot number: '+str(len(wtoday)+1))
            total_print+='\n'

        print(datetime.now().strftime('%A, %d. %B %Y %I:%M%p'))
        total_print+=str(datetime.now().strftime('%A, %d. %B %Y %I:%M%p'))
        total_print+='\n'

        if not save_flag:
            print('This is a test run. Use save_flag=True if you want to save data.')

        if save_flag:
            for ii in range(4):
                np.savetxt(str(fpQ+RunName+'ch'+str(ii)+'.txt'),cls.weichToPowerVsTime(PFN.HeadENB()[ii]*weiarr[-4+ii]))

        WeightedSum=np.sum(np.expand_dims(PFN.HeadENB(),axis=1)*weiarr[-4:],axis=0)

        if save_flag:
            np.savetxt(str(fpQ+RunName+'chsum.txt'),cls.weichToPowerVsTime(WeightedSum))

        PulseEnergies=[GLOBAL.EAB2w.get(),GLOBAL.EEF2w.get(),GLOBAL.EGH2w.get(),GLOBAL.EIJ2w.get()]
        EnMess='***'
        for ii in range(len(PulseEnergies)):
            EnMess+=' {}: '.format(headlist[ii])
            EnMess+=str(PulseEnergies[ii])
            EnMess+=' J ***'
        EnMess+=str(' total: '+str(np.sum(PulseEnergies))+' J ***')
        #print(EnMess)
        total_print+=EnMess
        total_print+='\n'#start over...
        total_print='Run: '+str(curr_run)+'\n'
        
        wppvlist=[GLOBAL.HWPAB, GLOBAL.HWPEF, GLOBAL.HWPGH, GLOBAL.HWPIJ];
        headstr='';wpstr='';
        for ii in range(4):
            if PFN.HeadENB()[ii]:
                headstr+=headlist[ii]
            wpstr=wpstr+headlist[ii]+': '+str(round(wppvlist[ii].get(),3))+', '
        wpen=np.mean([PFN.HeadENB()[ii]*np.cos((np.pi/180)*2*wppvlist[ii].get())**2 for ii in range(4)])
        total_print+=str('The following heads are enabled: '+headstr+'\n')
        total_print+=str('The waveplate settings are: '+wpstr[:-2]+'\n')
        total_print+=str('This is ~'+str(round(100*wpen,3))+'% of max energy.'+'\n')

        total_print+=EMeters.EGall(return_txt=True);

        if save_flag:
            np.savetxt(str(fpQ+RunName+'energies.txt'),PulseEnergies)
            efc.eplxysav(*cls.weichToPowerVsTime(WeightedSum),str(fpQ+RunName+'_'+str(int(round(np.sum(PulseEnergies))))+'J'),
                     abs_path=True,xlb='Time (ns)',ylb='Power (W)')
         

        fileheads=['w','y','s1in1w','s42in1w','s42in2w','s'];
        prepickle=[HAWG().ReadPulseHeights(),weiarr[0],weiarr[1],weiarr[2:6],weiarr[6:10],WeightedSum]
        Psns=GLOBAL.PSNS.get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        SSs=[GLOBAL.SSS.get()[2*ii:2*ii+2] for ii in range(len(Psns))]#pickle.load(open(psfpQ+'SSs.p','rb'))
        GLOBAL.WVFMHAWG.put(prepickle[0])
        GLOBAL.WVFMYFE.put(LPL.TraceFormatting(prepickle[1],cls.LMapAB(),PMQ=LPL.PulseMax(Psns,SSs,np.sum(prepickle[1]))*10))
        GLOBAL.WVFM1IN1w.put(LPL.TraceFormatting(prepickle[2],cls.LMapAB(),PMQ=LPL.PulseMax(Psns,SSs,np.sum(prepickle[2]))*10))
        GLOBAL.WVFM2IN1w.put(LPL.TraceFormatting(np.sum(prepickle[3],axis=0),cls.LMapAB(),PMQ=LPL.PulseMax(Psns,SSs,np.sum(prepickle[3]))*10))
        GLOBAL.WVFM2IN2w.put(LPL.TraceFormatting(np.sum(prepickle[5],axis=0),cls.LMap2(),PMQ=LPL.PulseMax(Psns,SSs,np.sum(prepickle[5]))*1))
        for ii in range(fileheads):
            templist=pickle.load(open(psfpQ+fileheads[ii]+DateStr+'.p','rb'))
            templist.append(prepickle[ii])
            efc.pickledump2(templist,psfpQ+fileheads[ii]+DateStr+'.p')

        if save_flag:
            mecel = elog.ELog({'experiment':ExpName},user='mecopr',pw=pickle.load(open(psfpQ+'elogauth.p','rb'))) 
            try:
                mecel.post(total_print,attachments=[str(fpQ+RunName+'_'+str(int(round(np.sum(PulseEnergies))))+'J.png')],run=curr_run,tags=['laser'])
                print('Auto-saved to eLog with run '+str(curr_run)) 
            except:
                try:
                    mecel.post(total_print,attachments=[str(fpQ+RunName+'_'+str(int(round(np.sum(PulseEnergies))))+'J.png')],tags=['laser']) 
                    print('Auto-saved to eLog')
                except:
                    print('Failed to auto-save to eLog')

    @classmethod
    def psefc(cls,JreqQ=0,AQQ=0.0):
        """Provides new update and plot based on single-shot full-energy output pulse"""
        (DateStr, curr_run) = cls.psheaders()

        psfpQ=cls.psfilepath()
        Psns=GLOBAL.PSNS.get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        SSs=[GLOBAL.SSS.get()[2*ii:2*ii+2] for ii in range(len(Psns))]#pickle.load(open(psfpQ+'SSs.p','rb'))
        Jscale=1
        PulseEnergies=[GLOBAL.EAB2w.get(),GLOBAL.EEF2w.get(),GLOBAL.EGH2w.get(),GLOBAL.EIJ2w.get()]
        wtoday=pickle.load(open(psfpQ+'w'+DateStr+'.p','rb'))
        stoday=pickle.load(open(psfpQ+'s'+DateStr+'.p','rb'))
        if int(JreqQ) < 5:
            Jreq=np.sum(PulseEnergies)*Jscale #65
        else:
            Jreq=JreqQ
        if Jreq>42:
            BumpQ=2#3 or 4
        else:
            BumpQ=1.5
        if np.abs(len(stoday[-1])-10000)<10:
            bkgrdbuffer=380
        elif np.abs(len(stoday[-1])-1000)<10:
            bkgrdbuffer=38
        else:
            print('Warning: unanticipated pulse shape array length of '+str(len(stoday[-1]))+', bkgrd subtraction may be off...')
            bkgrdbuffer=1
        try:#this is to check for bad scope traces; rewrite this later
            if np.sum(abs(stoday[-1])[:bkgrdbuffer]) > 1:
                BumpQ=0
                print('To whom it may concern: \n The next shot will not include an update to the pulse shaper because the saved scope trace seems abnormal. \n Did you switch to 10Hz and clear the scope traces before you saved them? or maybe did you disable or enable some amplifier arms before you saved the scope traces? or maybe you accidentally told the script to take more than 1 shot? \n If you answered yes to any of these questions, please don\'t do that again. (If you did something else out of the ordinary that could be linked to this anomaly, please let someone know.) \n\n Sincerely, \nThe Laser :)')
        except:
            pass
        AQ=AQQ*BumpQ#.03,.035,.05 seems too big most of the time
        if len(wtoday) > 0:
            if len(stoday[-1]) == 10002:
                mapnow=cls.LMap2()#[50,1000]
            elif len(stoday[-1]) == 1002:
                mapnow=cls.LMapAB()#[5,100]
            else:
                print('Unanticipated pulse shot array length: '+str(len(stoday[-1])));
                print('Aborting...');return
            wupd=cls.wIter2(stoday[-1],np.array(wtoday[-1])/28000.,Psns,SSs,Jreq,mapnow,AQ)
        else:
            print('No shots yet today; beginning with pre-loaded shape')
            try:
                wupd=HAWG().ReadPulseHeights()
            except:
                print('Error! HAWG')
        return wupd
        ##EXECUTE THIS FILE FIRST TO DETERMINE THE UPDATED WAVEFORM

    @classmethod
    def psefc10Hz(cls,pwt=0,numIterQ=50,AQQ=0.03,displayPlot=True,reloopPrompt=True):
        """Performs 10Hz feedback according to user-supplied YFE waveform target"""
        GLOBAL.EVRLPLSSEC.put(43)
        #psfpQ=cls.psfilepath()   
        Psns=GLOBAL.PSNS.get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        SSs=[GLOBAL.SSS.get()[2*ii:2*ii+2] for ii in range(len(Psns))]#pickle.load(open(psfpQ+'SSs.p','rb'))
        pwtF=np.array(cls.TraceFormatting(pwt,[25,500],1,AvgRange=1,FWHM=1))
        GLOBAL.WVFMYFEGOAL.put(pwtF)

        try:
            SLA=LOSC('A');SLA.Open();SH=HAWG();SH.Open();#replaced w/LeCroyA
            ops00=SLA._rch(1);time.sleep(0.1);

            meanerr=[]
            meanerr.append(np.sum(np.abs(pwtF[:26]-cls.TraceFormatting(ops00,cls.LMapAB(),1,AvgRange=1,FWHM=1)[:26])/pwtF[:26])/len(pwtF[:26]))
            ops00F=cls.TraceFormatting(ops00,cls.LMapAB(),1,AvgRange=1,FWHM=1)
            print('Start loop')
            if displayPlot:
                plt.ion()
                fig,axs=plt.subplots(1,2,figsize=(10,5),gridspec_kw={'hspace':0.4,'wspace':0.2});
                #xdat=[[startposlist[ii]+alphlist[ii]*(-.1+.02*(jj)) for jj in range(11)] for ii in range(4)]
                xdat=np.linspace(.25,35,140)
                #ydat=[[0]*11 for ii in range(4)]
                ax1,=axs[0].plot(xdat,ops00F); axs[0].set_xlabel('curr vs goal'); plt.pause(0.01); axs[0].plot(xdat,pwtF); plt.pause(0.01);
                ax2,=axs[1].plot(list(range(len(meanerr))),meanerr); axs[1].set_xlabel('mean error'); plt.pause(0.01);
                axss=[ax1,ax2]
                #xdat=[]
            
            LoopIsDone=False
            while not LoopIsDone:
                for ii in range(numIterQ):
                    if (ii+1)%50 == 0:
                        print(str('Iter:'+str(ii+1)))
                    ops0=SLA._rch(1);time.sleep(0.025);####added.215 when 200mV/div instead of 100mV/div
                    if all(ops0 == ops00):
                        print('No scope update detected... no feedback applied!')
                    else:
                        rph=SH._ReadPulseHeights();time.sleep(0.025);
                        #pwtF=np.array(TraceFormatting2(pwt,[25,500],1))
                        ##ops0F=TraceFormatting(ops0,[25,500],1)
                        ops0F=cls.TraceFormatting(ops0,cls.LMapAB(),1,AvgRange=1,FWHM=1)
                        #epll([pwtF,ops0F])
                        meanerr.append(np.sum(np.abs(pwtF[:26]-ops0F[:26])/pwtF[:26])/len(pwtF[:26]));
                        if displayPlot:
                            axss[0].set_data(xdat,ops0F); axs[0].relim(); axs[0].autoscale_view(True,True,True);
                            axss[1].set_data(list(range(len(meanerr))),meanerr); axs[1].relim(); axs[1].autoscale_view(True,True,True);
                            fig.canvas.draw_idle(); plt.pause(0.01);
                        usa0=cls.UpdatingShapingAlgorithm(pwtF,ops0F,np.array(rph)/28000.,AQQ)#.075#.25
                        usa0FE=cls.FixEdges(usa0,Psns,SSs)
                        #usa0FE=FixEdges(usa0,[3,4.25],[[.98*100/8.0,100/8.0],[98,100]])
                        #epll([rph,usa0FE*28000.])
                        SH._WritePulseHeights(usa0FE*28000.);time.sleep(0.05);
                        ops00=ops0[:]
    #        if displayPlot:
    #            epll([pwtF,ops0F])
    #            epl(meanerr)
                ######check and aim
                if reloopPrompt:
                    print('Would you like to try another 50 iterations? [enter y/n]',end='',flush=True)
                    checkprompt=input();
                    if checkprompt.lower().startswith('y'):
                        numIterQ=50
                        if len(checkprompt)>1:
                            try:
                                numIterQ=int(checkprompt[1:]) if int(checkprompt[1:]) < 250 else 250
                            except:
                                print('Could not parse loop # instruction -- using 50 instead.')
                        print('Performing {} more iterations!'.format(str(numIterQ)))
                    else:
                        LoopIsDone=True
                else:
                    LoopIsDone=True
            SH.Close();SLA.Close();
        except:
            print('Failed')
            SH.Close();SLA.Close();
        if displayPlot:
            plt.ioff()

    @classmethod
    def psupd(cls,newwavQ):
        """Shortcut to updating the Highland based according to provided waveform"""
        cls.pshostcheck()
        wupdt=newwavQ[:]
        if max(wupdt) < 1.5:
            wupdt=28000.0*np.array(wupdt)
        try:
            HAWG().WritePulseHeights(wupdt);
        except:
            print('Error, check HAWG!')

    @classmethod
    def psloadwvfm(cls,RecipeStrQ,WvGoal10HzHint=False):
        """Loads a new waveform according to previously-saved recipe"""
        cls.pshostcheck()
        print('Loading timestamp: '+datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
        try:
            [Psns,SSs,YFE02mmCurr,YFE06mmCurr,YFE10mmCurr,NewWvfm,WvGoal10Hz] = pickle.load(open(cls.psfilepath()+'recipes/load'+RecipeStrQ+'.p','rb'))
        except:
            print('Recipe file '+cls.psfilepath()+'recipes/load'+RecipeStrQ+'.p\' not found.')
            return
        efc.pickledump2(Psns,cls.psfilepath()+'Psns.p')
        GLOBAL.PSNS.put(Psns)
        efc.pickledump2(SSs,cls.psfilepath()+'SSs.p')
        GLOBAL.SSS.put([elem for sublist in SSs for elem in sublist])
        YFE.Set(2,YFE02mmCurr);time.sleep(.15);
        YFE.Set(6,YFE06mmCurr);time.sleep(.15);
        YFE.Set(10,YFE10mmCurr);time.sleep(.15);
        if WvGoal10HzHint:
            if len(WvGoal10Hz) > 5:
                print('Hint for 10Hz waveform: '+WvGoal10Hz)
            else:
                print('No hint for 10Hz waveform available. MSG: \''+WvGoal10Hz+'\'')
        pseparams=np.array(re.findall('ExponentialWave2\((\d+),(\d+\.\d+|\.\d+),(\d+),(\d+\.\d+|\.\d+),0,5002\)',WvGoal10Hz),dtype=np.float32);
        pslparams=np.array(re.findall('LinearWave2\((\d+),(\d+\.\d+|\.\d+),(\d+),(\d+\.\d+|\.\d+),0,5002\)',WvGoal10Hz),dtype=np.float32);  
        yfegoal=cls.LinearWave2(500,0,1025,0,0,5002);
        if (len(pslparams) > 0) or (len(pseparams) > 0):
            for ii in range(len(pslparams)): 
                yfegoal+=cls.LinearWave2(pslparams[ii][0],pslparams[ii][1],pslparams[ii][2],pslparams[ii][3],0,5002)
            for ii in range(len(pseparams)): 
                yfegoal+=cls.ExponentialWave2(pseparams[ii][0],pseparams[ii][1],pseparams[ii][2],pseparams[ii][3],0,5002)
        else:
            print('No wave extracted: '+WvGoal10Hz)
        yfegoalF=np.array(cls.TraceFormatting(yfegoal,[25,500],1,AvgRange=1,FWHM=1))
        GLOBAL.WVFMYFEGOAL.put(yfegoalF)
        GLOBAL.WVFM2IN2wGOAL.put(cls.PulseGoal(Psns,SSs))
        try:
            cls.psupd(NewWvfm)
            print('New waveform loaded! ')
        except:
            print('Failed to load new waveform.')

    @classmethod
    def pssavewvfm(cls,PsnsQ=0,SSsQ=0,YFEgetQ=0,TargetwlistDateQ='curr',TargetwindexQ=0,RecipeStrQ=0,WvGoal10HzQ='none'):
        """Saves a new pulse shape recipe using user-provided parameters"""
        cls.pshostcheck()
        print('Saving timestamp: '+datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
        if TargetwlistDateQ == 'curr':
            print('Using current Highland waveform...')
            try:
                NewWvfmQ=HAWG().ReadPulseHeights();
            except:
                print('Failed! HAWG');
                return False
            try:
                wlastarr=pickle.load(open(cls.psfilepath()+'w'+cls.DateString()+'.p','rb'))
                if wlastarr[-1] == NewWvfmQ:
                    print('Pulse looks equivalent to the most recent pulse, w'+cls.DateString()+'['+str(len(wlastarr)-1)+']')
                    WvGoal10HzQ=WvGoal10HzQ+';; w'+cls.DateString()+'['+str(len(wlastarr)-1)+']'
                else:
                    WvGoal10HzQ=WvGoal10HzQ+';; sometime after most recent w'+cls.DateString()+'['+str(len(wlastarr)-1)+']'
            except:
                print('Failed to load most recent amplified shot.')
        else:
            wavehistQ=pickle.load(open(cls.psfilepath()+'w'+TargetwlistDateQ+'.p','rb'))
            NewWvfmQ=wavehistQ[TargetwindexQ][:]
            WvGoal10HzQ=WvGoal10HzQ+', w'+TargetwlistDateQ+'['+str(TargetwindexQ)+']'

        if PsnsQ == 0:
            #PsnsQ = pickle.load(open(cls.psfilepath()+'Psns.p','rb'))
            PsnsQ=GLOBAL.PSNS.get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        if SSsQ == 0:
            #SSsQ = pickle.load(open(cls.psfilepath()+'SSs.p','rb'))
            SSsQ=[GLOBAL.SSS.get()[2*ii:2*ii+2] for ii in range(len(PsnsQ))]#pickle.load(open(psfpQ+'SSs.p','rb'))
        if YFEgetQ == 0:
            YFEgetQ=YFE.Get(display=False)
        [YFE02mmCurrQ,YFE02mmCurrQ,YFE02mmCurrQ,YFE02mmCurrQ,YFE06mmCurrQ,YFE10mmCurrQ]=YFEgetQ
        if RecipeStrQ == 0:
            #learn formatting for pulse from Psns and SSs
            pass
        oldfilefound=True

        try:
            dummyQ = pickle.load(open(cls.psfilepath()+'recipes/load'+RecipeStrQ+'.p','rb'))
            print('Old recipe found with same name: '+cls.psfilepath()+'recipes/load'+RecipeStrQ+'.p')
        except:
            oldfilefound=False
        iiQ=0
        while oldfilefound:
            iiQ=iiQ+1
            try:
                pickle.load(open(cls.psfilepath()+'recipes/load'+RecipeStrQ+'_'+str(iiQ).zfill(2)+'.p','rb'));
            except:
                oldfilefound=False
                dummyQ[-1]='**replaced on '+cls.DateString()+'** '+dummyQ[-1]
                efc.pickledump2(dummyQ,cls.psfilepath()+'recipes/load'+RecipeStrQ+'_'+str(iiQ).zfill(2)+'.p')
                print('Saved old recipe as '+cls.psfilepath()+'recipes/load'+RecipeStrQ+'_'+str(iiQ).zfill(2)+'.p')
        try:
            efc.pickledump2([PsnsQ,SSsQ,YFE02mmCurrQ,YFE06mmCurrQ,YFE10mmCurrQ,NewWvfmQ,WvGoal10HzQ],cls.psfilepath()+'recipes/load'+RecipeStrQ+'.p')
            print('Saved new recipe as '+cls.psfilepath()+'recipes/load'+RecipeStrQ+'.p')
        except:
            print('Failed to save new recipe.')

    @classmethod
    def psviewwvfm(cls,RecipeStrQ='none',TargetwlistDateQ='curr',TargetwindexQ=0,WvGoal10HzHint=False):
        """Displays waveform and parameters that are part of a previously-saved recipe"""
        cls.pshostcheck()
        if RecipeStrQ == 'none':
            if TargetwlistDateQ == 'curr':
                foundlastshot=False
                iidQ=0
                while not foundlastshot:
                    try:
                        wlastarr=pickle.load(open(cls.psfilepath()+'w'+str(int(cls.DateString())-iidQ)+'.p','rb'))
                        slastarr=pickle.load(open(cls.psfilepath()+'s'+str(int(cls.DateString())-iidQ)+'.p','rb'))
                        wlast=wlastarr[-1][:]
                        slast=slastarr[-1][:]
                        print('Retrieving most recent shot: w'+str(int(cls.DateString())-iidQ)+'['+str(len(wlastarr)-1)+']')
                        foundlastshot=True
                    except:
                        iidQ=iidQ+1
            else:
                try:
                    wlastarr=pickle.load(open(cls.psfilepath()+'w'+TargetwlistDateQ+'.p','rb'))
                    slastarr=pickle.load(open(cls.psfilepath()+'s'+TargetwlistDateQ+'.p','rb'))
                    wlast=wlastarr[TargetwindexQ]
                    slast=slastarr[TargetwindexQ]
                except:
                    print('Failed to load at given date and index: '+TargetwlistDateQ+', '+str(TargetwindexQ))
            efc.epl(wlast)
            efc.epl(slast)
        else:
            try:
                [Psns,SSs,YFE02mmCurr,YFE06mmCurr,YFE10mmCurr,NewWvfm,WvGoal10Hz] = pickle.load(open(cls.psfilepath()+'recipes/load'+RecipeStrQ+'.p','rb'))
            except:
                print('Recipe file '+cls.psfilepath()+'recipes/load'+RecipeStrQ+'.p\' not found.')
                return
            print('Retrieved recipe: load'+RecipeStrQ)
            print('Psns: '+str(Psns)+', SSs: '+str(SSs))
            print('YFEcurr:: 2mm: '+'{:5.1f}'.format(YFE02mmCurr)+', 6mm: '+'{:5.1f}'.format(YFE06mmCurr)+', 10mm: '+'{:5.1f}'.format(YFE10mmCurr))
            print('Extra info: '+WvGoal10Hz)
            efc.epl(NewWvfm)
            try:
                tempstr=WvGoal10Hz[-18:]
                wstrinx=tempstr.index('w')
                wstrinx0=wstrinx-len(tempstr)
                loaddate=tempstr[wstrinx0+1:wstrinx0+9]
                loadindx=int(tempstr[wstrinx0+10:-1])
                saskarr=pickle.load(open(cls.psfilepath()+'s'+loaddate+'.p','rb'))
                efc.epl(saskarr[loadindx])
            except:
                print('Failed to load 2w waveform for display.')
            return NewWvfm

    @classmethod
    def psrefrwvfm(cls,RecipeStrQ,numStepsQ=50,stepSizeQ=0.25,displayPlot=True,reloopPrompt=True):
        """Refreshes 10Hz YFE waveform according to target shape given by previously-saved recipe"""
        cls.pshostcheck()
        print('Loading timestamp: '+datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
        #load and extract the pulse target from the desired recipe
        try:
            [Psns,SSs,YFE02mmCurr,YFE06mmCurr,YFE10mmCurr,NewWvfm,WvGoal10Hz] = pickle.load(open(cls.psfilepath()+'recipes/load'+RecipeStrQ+'.p','rb'))
        except:
            print('Recipe file '+cls.psfilepath()+'recipes/load'+RecipeStrQ+'.p\' not found.')
            return
        print('Hint text: '+WvGoal10Hz)
        pseparams=np.array(re.findall('ExponentialWave2\((\d+),(\d+\.\d+|\.\d+),(\d+),(\d+\.\d+|\.\d+),0,5002\)',WvGoal10Hz),dtype=np.float32);
        pslparams=np.array(re.findall('LinearWave2\((\d+),(\d+\.\d+|\.\d+),(\d+),(\d+\.\d+|\.\d+),0,5002\)',WvGoal10Hz),dtype=np.float32);  
        yfegoal=cls.LinearWave2(500,0,1025,0,0,5002);
        if (len(pslparams) > 0) or (len(pseparams) > 0):
            for ii in range(len(pslparams)): 
                yfegoal+=cls.LinearWave2(pslparams[ii][0],pslparams[ii][1],pslparams[ii][2],pslparams[ii][3],0,5002)
            for ii in range(len(pseparams)): 
                yfegoal+=cls.ExponentialWave2(pseparams[ii][0],pseparams[ii][1],pseparams[ii][2],pseparams[ii][3],0,5002)
        else:
            print('No wave extracted: '+WvGoal10Hz)
            return
        #close the shutters
        #print('NOT closing the shutters... hope you\'re protecting your sample!')
        print('Closing all shutters...')
        TTL_shutter.Toggle('closeall',display=False)#close all the shutters
        time.sleep(4)
        try:#used to make sure shutters re-open even in case of error or KeyboardInterrupt
            #check if laser is on
            if np.sum(YFE.Get(display=False)) < 400:
                if np.sum(YFE.Get(display=False)) < 20:
                    print('WARNING: YFE seems to be off... Attempting to turn on YFE...')
                    YFE.On();
                else:
                    print('WARNING: eDrive currents seem low...')
                    #prechk = False
            #turn off bias dither; AUTO is 0; MAN is 1
            if GLOBAL.MBCmode.get() != 1:
                print('Turning off MBC bias dither...',end='',flush=True);
                currbias=GLOBAL.MBCbias.get();
                GLOBAL.MBCmode.put(1);#set bias mode to MAN
                efc.dotsleep(3);
                GLOBAL.MBCbias.put(currbias+10);
                print('Testing bias responsivity...',end='',flush=True);
                for ii in range(2):
                    time.sleep(1);print('..',end='',flush=True);
                if np.abs(GLOBAL.MBCbias.get() - (currbias+10)) < 3:
                    GLOBAL.MBCbias.put(currbias);
                    efc.dotsleep(2);
                else:
                    print('*')
                    print('WARNING: MBC not responding!!')
            #else:
            #    print('MBC is not safe! Resetting the MBC...')#problem is here...
            #    resetMBC();
            #set and enable 10Hz output
            if GLOBAL.EVRLPLSSEC.get() != 43:
                GLOBAL.EVRLPLSSEN.put(0)
                GLOBAL.EVRLPLSSEC.put(43)
            if GLOBAL.EVRLPLSSEN.get() != 1:
                GLOBAL.EVRLPLSSEN.put(1)
            #run the update code
            print('Refreshing the YFE wavefront...')
            cls.psefc10Hz(pwt=yfegoal,numIterQ=numStepsQ,AQQ=stepSizeQ,displayPlot=displayPlot,reloopPrompt=reloopPrompt)
            #reset to single shot on pulse picker
            GLOBAL.EVRLPLSSEN.put(0);time.sleep(0.5);
            GLOBAL.EVRLPLSSEC.put(182);time.sleep(0.5);
            GLOBAL.EVRLPLSSEN.put(1)
            #re-open shutters
            print('Opening all shutters...')
            TTL_shutter.Toggle('openall',display=False);#open all the shutters
            MBC.resetMBC();
            YFE.SetAll(True,displayQ=False);
        except:#used to make sure shutters re-open even in case of error or KeyboardInterrupt
            #reset to single shot on pulse picker
            GLOBAL.EVRLPLSSEN.put(0);time.sleep(0.5);
            GLOBAL.EVRLPLSSEC.put(182);time.sleep(0.5);
            GLOBAL.EVRLPLSSEN.put(1)
            #re-open shutters
            #print('Opening all shutters...')
            #toggle_TTL_shutter('openall',display=False);#open all the shutters
            MBC.Reset();
            YFE.SetAll(True,displayQ=False);

    @classmethod
    def psrecipes(cls):
        """Prints and returns a list of all previously-saved pulse recipes"""
        allrec=glob.glob(cls.psfilepath()+'recipes/*.p');
        oldrec=glob.glob(cls.psfilepath()+'recipes/*_*.p')
        currec=[ext[60:-2] for ext in allrec if ext not in oldrec];currec.sort();
        return currec

    @classmethod
    def pspreshot(cls):
        """Prepares state of the laser for taking a single full-energy shot"""
        cls.pshostcheck()
        if not YFE.OnCheck(display=False):
            print('WARNING: YFE seems to be off... Attempting to turn on YFE...')
            YFE.On();
        if np.sum(YFE.Get(display=False)) < 550:
            if np.sum(YFE.Get(display=False)) < 20:
                print('WARNING: YFE seems to be turned down. Attempting to turn up YFE...')
                YFE.SetAll(True,displayQ=True)
        if MBC.IsSafe():
            print('Turning off MBC bias dither...',end='',flush=True);
            currbias=GLOBAL.MBCbias.get();
            GLOBAL.MBCmode.put(1);#set bias mode to MAN
            efc.dotsleep(3);
            GLOBAL.MBCbias.put(currbias+10);
            print('Testing bias responsivity...',end='',flush=True);
            for ii in range(2):
                time.sleep(1);print('..',end='',flush=True);
            if np.abs(GLOBAL.MBCbias.get() - (currbias+10)) < 3:
                GLOBAL.MBCbias.put(currbias);
                efc.dotsleep(2);
            else:
                print('*')
                print('WARNING: MBC not responding!!')
                prechk = False
        else:
            print('MBC is not safe! Resetting the MBC...')
            MBC.Reset();##up above, check emission AND check currents???
            YFE.SetAll(True,displayQ=False)
            cls.pspreshot()
            prechk = False
        if GLOBAL.EVRLPLLAMPEC.get() != 182:
            GLOBAL.EVRLPLLAMPEN.put(0)
            GLOBAL.EVRLPLLAMPEC.put(182)
        if GLOBAL.EVRLPLSSEC.get() != 182:
            GLOBAL.EVRLPLSSEN.put(0)
            GLOBAL.EVRLPLSSEC.put(182)
        if GLOBAL.EVRLPLLAMPEN.get() != 1:
            GLOBAL.EVRLPLLAMPEN.put(1)
        if GLOBAL.EVRLPLSSEN.get() != 1:
            GLOBAL.EVRLPLSSEN.put(1)
        
        wppvlist=[GLOBAL.HWPAB, GLOBAL.HWPEF, GLOBAL.HWPGH, GLOBAL.HWPIJ];
        headstr='';wpstr='';headlist=['AB','EF','GH','IJ'];
        for ii in range(4):
            if PFN.HeadENB()[ii]:
                headstr+=headlist[ii]
            wpstr=wpstr+headlist[ii]+': '+str(round(wppvlist[ii].get(),3))+', '
        wpen=np.mean([PFN.HeadENB()[ii]*np.cos((np.pi/180)*2*wppvlist[ii].get())**2 for ii in range(4)])
        print('The following heads are enabled: '+headstr)
        print('The waveplate settings are: '+wpstr[:-2])
        print('This is ~'+str(round(100*wpen,3))+'% of max energy.')
        Psns=GLOBAL.PSNS.get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        SSs=[GLOBAL.SSS.get()[2*ii:2*ii+2] for ii in range(len(Psns))]#pickle.load(open(psfpQ+'SSs.p','rb'))
        print('Current pulse target is: '+str(Psns)+' ns, '+str(SSs)+' % of max power.')
        not_charging='';not_enabled='';yes_enabled='';headchanlist=['CD','A','B','E','F','G','H','I','J'];
        temppv1=[GLOBAL.PFNCDEN,GLOBAL.PFNAEN,GLOBAL.PFNBEN,GLOBAL.PFNEEN,GLOBAL.PFNFEN,GLOBAL.PFNGEN,GLOBAL.PFNHEN,GLOBAL.PFNIEN,GLOBAL.PFNJEN]
        temppv2=[GLOBAL.PFNCDCS,GLOBAL.PFNACS,GLOBAL.PFNBCS,GLOBAL.PFNECS,GLOBAL.PFNFCS,GLOBAL.PFNGCS,GLOBAL.PFNHCS,GLOBAL.PFNICS,GLOBAL.PFNJCS]
        for ii in range(9):
            if temppv1[ii].get() == 0:
                not_enabled+=headchanlist[ii]
            if temppv1[ii].get() == 1:
                yes_enabled+=headchanlist[ii]
            if (temppv1[ii].get() == 1) and (temppv2[ii].get() == 0):
                not_charging+=headchanlist[ii]
        if len(not_charging)>0:
            print('** WARNING: The following heads are enabled but NOT charging: '+not_charging) 

        TTL_shutter.Toggle('open'+yes_enabled+'wwxx',display=False);time.sleep(1);#make sure all shutters are open...
        TTL_shutter.Toggle('close'+not_enabled,display=False)#close shutters that aren't enabled
        return prechk
        #waveform pre-check? verify shutters are open?

    @classmethod
    def pspostshot(cls,save_flag=True,RunNumQQ=9000):
        """Executes post-shot routine for saving data and returning laser to appropriate state"""
        GLOBAL.EVRLPLLAMPEN.put(0)
        GLOBAL.EVRLPLSSEN.put(0)
        cls.pshostcheck()
        cls.psacqx(save_flag=save_flag,RunNumQQ=RunNumQQ)#took out _noLecroyA
        #psefc();
        TTL_shutter.Toggle('openall',display=False);#make sure all shutters are open again...
        print('Resetting bias tracking...')
        MBC.Reset();
        YFE.SetAll(True);
        EMeters.E_synth_refresh();

    @classmethod
    def SHG_opt(cls,armsQ='ABEFGHIJ'):#check for trace height;#All shutters must start in the open state... 
        """Optimizes the tuning of the doubler angles to maximize the conversion efficiency of the arms of the LPL"""
        print('Running this routine requires ALL TTL shutters to begin in the open state! The YFE must be on with the bias dither initially enabled!')
        if np.sum(TTL_shutter.Status(display=False)[-1]) > 0:
            print('Warning! The shutters don\'t all appear to be open! ',end='',flush=True);TTL_shutter.Status(display=True);
        else:
            print('(Shutters seem OK...)')
        if not YFE.OnCheck(display=False):
            print('Warning! The YFE doesn\'t appear to be on! ',end='',flush=True);YFE.OnCheck(display=True);
        else:
            print('(YFE emission seems OK...)')
        if np.sum(YFE.Get(display=False)) < 550:
            print('Warning! The YFE doesn\'t appear to be turned up! ');YFE.Get(display=True);
        else:
            print('(YFE current seems OK...)')
        if MBC.ModeCheck() != 0:
            print('(Warning! The MBC doesn\'t appear to be in AUTO mode!)')
        else:
            print('(MBC mode seems OK...)')
        print('Are you sure you are ready to proceed? [enter y/n]',end='',flush=True)
        checkprompt=input();
        if checkprompt.lower() != 'y':
            print('Try again later then!');
            return
        else:
            print('OK, I hope you know what you\'re doing!')
        HWP.On('all',set_T=1)
        armlist=['AB','EF','GH','IJ']
        #YFEoff();YFEon();
        GLOBAL.MBCmode.put(1)#set MAN mode on MBC
        if np.sum(YFE.Get(display=False)) < 100:
            print('Check YFE before optimizing!')
        optwvfm=pickle.load(open(cls.psfilepath()+'opttrace.p','rb'));
        try:
            oldwvfm=HAWG().ReadPulseHeights();
            HAWG().WritePulseHeights(optwvfm);
        except:
            print('Failed! HAWG')
        SHGpvlist=[GLOBAL.SHGABmot, GLOBAL.SHGEFmot, GLOBAL.SHGGHmot, GLOBAL.SHGIJmot];
        print('Closing all shutters...')
        TTL_shutter.Toggle('closeall',display=False);#close all the shutters
        time.sleep(4)
        GLOBAL.EVRLPLSSEC.put(43);GLOBAL.EVRLPLSSEN.put(1);#enable these...
        try:
            tempchk1=LOSC('a').rch(1);time.sleep(.15);tempchk2=LOSC('a').rch(1);
            if np.sum(np.abs(tempchk1-tempchk2))<1e-6:
                print('Warning: scope trace doesn\'t appear to be updating, please check scope! Abort? [enter y/n]')
                checkprompt=input();
                if checkprompt.lower() != 'y':
                    print('Try again later then!');
                    HAWG().WritePulseHeights(oldwvfm);
                    return
                else:
                    print('OK, I hope you know what you\'re doing!')
        except:
            print('Scope error, check scope status! Aborting...')
            HAWG().WritePulseHeights(oldwvfm);
            return
        startposlist=[SHGrbv.get() for SHGrbv in SHGpvlist];
        newposlist=startposlist[:]
        alphlist=[1,0.5,0.5,1];
        for ii in range(4):
            if armlist[ii] in armsQ:#only prep the stage if it's going to be used
                SHGpvlist[ii].put(startposlist[ii]+alphlist[ii]*(-.1+.01*0))
        currentshutter=0;#trying to re-open a shutter in case of failure...
        #set up all the plotting stuff
        plt.ion()
        fig,axs=plt.subplots(2,2,gridspec_kw={'hspace':0.4,'wspace':0.3})
        xdat=[[startposlist[ii]+alphlist[ii]*(-.1+.02*(jj)) for jj in range(11)] for ii in range(4)]
        ydat=[[0]*11 for ii in range(4)]
        ax1,=axs[0,0].plot(xdat[0],ydat[0]); axs[0,0].set_xlabel('AB'); plt.pause(0.01);
        ax2,=axs[0,1].plot(xdat[1],ydat[1]); axs[0,1].set_xlabel('EF'); plt.pause(0.01);
        ax3,=axs[1,0].plot(xdat[2],ydat[2]); axs[1,0].set_xlabel('GH'); plt.pause(0.01);
        ax4,=axs[1,1].plot(xdat[3],ydat[3]); axs[1,1].set_xlabel('IJ'); plt.pause(0.01);
        axss=[ax1,ax2,ax3,ax4]
        try:
            SLA=LOSC('A');SLA.Open();#changed to LecroyA since repair
            for ii in range(4):
                if armlist[ii] in armsQ:
                    print('Begin optimizing '+armlist[ii]+'... ',end='',flush=True);
                    shgarmdatax,shgarmdatay=[],[]
                    TTL_shutter.Toggle('open'+armlist[ii],display=False);currentshutter=ii;time.sleep(4);print('Shutter opened!');#open one shutter
                    for jj in range(11):
                        print('.',end='',flush=True)
                        SHGpvlist[ii].put(startposlist[ii]+alphlist[ii]*(-.1+.02*(jj)));time.sleep(2.5);#step to new position
                        curr_x=SHGpvlist[ii].get();curr_y=np.max(SLA._rch(3));time.sleep(.15);#in testing, max is more stable than sum
                        if curr_y > 0.005:#threshold so don't skew fit with noise; max is ~~10x this
                            shgarmdatax.append(curr_x);shgarmdatay.append(curr_y);#save x and y
                        print('.',end='',flush=True)
                        ydat[ii][jj]=curr_y
                        axss[ii].set_data(xdat[ii],ydat[ii])
                        axs[ii//2,ii%2].set_ylim((min(ydat[ii]),max(ydat[ii])))
                        plt.pause(0.01)
                        #axs[ii//2,ii%2].autoscale(True,True,True)
                        fig.canvas.draw_idle()
                        plt.pause(0.01)
                    print('*')
                    qfit=np.polyfit(shgarmdatax,shgarmdatay,2);newpos=qfit[1]/(-2*qfit[0]);#find fit and new max
                    if np.abs(startposlist[ii]-newpos)<.15:
                        SHGpvlist[ii].put(newpos);newposlist[ii]=newpos;
                        print('SHG position on arm '+armlist[ii]+' changed from '+str(round(startposlist[ii],4))+' to '+str(round(newpos,4)))
                    else:
                        print('Failed! New SHG position on arm '+armlist[ii]+' seems too far off... '+str(round(newpos,4))+' from '+str(round(startposlist[ii],4))+'... Restoring...')
                        SHGpvlist[ii].put(startposlist[ii])
                    TTL_shutter.Toggle('close'+armlist[ii],display=False);currentshutter=0;#close that shutter;
                    #xpq=np.arange(startposlist[ii]+alphlist[ii]*(-.1+.02*(-1)),startposlist[ii]+alphlist[ii]*(-.1+.02*(11)),.0001);
                    qfitp=np.poly1d(qfit);
                    axs[ii//2,ii%2].plot(xdat[ii],qfitp(xdat[ii]))
                    axs[ii//2,ii%2].relim();plt.pause(0.01);
                    axs[ii//2,ii%2].autoscale(True,True,True)
                    fig.canvas.draw_idle()
                    plt.pause(0.01)                
                    #epllxy([[shgarmdatax,shgarmdatay],[xpq,qfitp(xpq)]],xlb=armlist[ii])
                else:
                    print('Skipping '+armlist[ii]+'...')
                    pass
            SLA.Close();time.sleep(.15);#changed to LeCroyA
        except:
            print('Failed! Restoring original values and attempting to re-open most-recent shutter... you should verify!')
            SLA.Close();time.sleep(.15);#changed to LeCroyA
            if currentshutter > 0:
                TTL_shutter.Toggle('open'+armlist[currentshutter],display=False);
            for ii in range(4):
                SHGpvlist[ii].put(startposlist[ii]);newposlist[ii]=startposlist[ii];
        time.sleep(2);#need time so that last shutter trigger ends before trying to open IJ
        try:
            HAWG().WritePulseHeights(oldwvfm);
        except:
            print('Error! Check waveform!')
        GLOBAL.EVLPLSSEN.put(0);#disable PC before re-opening shutters
        datestamp=int(datetime.now().strftime('%Y%m%d%H%M%S'))
        SHGlog=pickle.load(open(cls.psfilepath()+'SHG_opt_log.p','rb'))
        SHGlog.append([datestamp,[newposlist[ii] for ii in range(4)]])
        efc.pickledump2(SHGlog,cls.psfilepath()+'SHG_opt_log.p')
        TTL_shutter.Toggle('openall',display=False);#open all the shutters
        MBC.Reset();YFE.SetAll(True);#reset bias...
        plt.ioff()
        
        
        

        
class efc:
    """Extra Function Class for convenient shorthand ways of doing common things like plotting data, printing in color, and other frequent actions"""
    def epl(listq):
        """Shorthand plotting function for a single list of y-values"""
        df1=plt.figure()
        plt.plot(listq);
        df1.show()
        return

    def eplxy(listxq,listyq):
        """Shorthand plotting function for a single list of x-values and y-values"""
        df1=plt.figure()
        plt.plot(listxq,listyq);
        df1.show()
        return

    def eplxyloglog(listxq,listyq):
        """Shorthand plotting function for a single list of x-values and y-values with logrithmic axes in both directions"""
        df1=plt.figure()
        plt.loglog(listxq,listyq);
        df1.show()
        return

    def eplsav(listq,FileNameQ,blockdisplay=True):
        """Shorthand function for saving a plot of a single list of y-values"""
        df1=plt.figure()
        plt.plot(listq);
        df1.savefig(str(FileNameQ+'.png'))
        if blockdisplay:
            plt.close(df1)
        return
        
    def eplxysav(listxq,listyq,FileNameQ,abs_path=False):
        """Shorthand function for saving a plot of a single list of x-values and y-values"""
        df1=plt.figure()
        plt.plot(listxq,listyq);
        if abs_path:
            figfilename=FileNameQ;
        else:
            figfilename=str(LPL.psfilepath()+FileNameQ+'.png')
        df1.savefig(figfilename)        
        plt.close(df1)
        return
        
    def eplcomp(listq,goalq,Map,tMax):
        """Shorthand function for comparing a pulse shape to its targeted pulse shape; not used much currently"""
        formtra=[]
        formtra.append(LPL.TraceFormatting(listq,Map,tMax))
        formtra.append(goalq)
        efc.epll(formtra)
        return


    def eplcsv(CSVname):
        """Shorthand plotting function for a single list of y-values loaded from a CSV file; hasn't been used for ages"""
        with open(LPL.psfilepath()+'data/'+CSVname+'.csv','r') as filehead:
            RawListQ=filehead.read()
            ListedValues=RawListQ.split('\n')
        efc.epl(ListedValues[:-1])
        return 

    def epllcsv(CSVHeadname):
        """Shorthand plotting function for a nested list of four sets of y-values loaded from a CSV file; hasn't been used for ages"""
        ListofListedValues=[]
        for ii in range(1,5):
            with open(LPL.psfilepath()+'data/'+CSVHeadname+'_ch'+str(ii)+'.csv','r') as filehead:
                RawListQ=filehead.read()
                ListedValues=RawListQ.split('\n')
            ListofListedValues.append(ListedValues[:-1])
        efc.epll(ListofListedValues)
        return 

    def rcsv(CSVname):
        """Shorthand function for reading a list of y-values from a CSV file; hasn't been used for ages"""
        with open(LPL.psfilepath()+'data/'+CSVname+'.csv','r') as filehead:
            RawListQ=filehead.read()
            if '\n' in RawListQ:
                ListedValues=RawListQ.split('\n')
            elif '\r\n' in RawListQ:
                ListedValues=RawListQ.split('\r\n')
            elif ',' in RawListQ:
                ListedValues=RawListQ.split(',')
            else:
                print('Unrecognized format on input file.')
        return ListedValues

    def epll(llist):
        """Shorthand plotting function for a list of multiple lists of y-values"""
        df1=plt.figure()
        for ii in range(len(llist)):
            plt.plot(llist[ii]);
        df1.show()
        
    def epllxy(llistxyq,xlb='none',ylb='none'):
        """Shorthand plotting function for a list of multiple lists of x-values and y-values"""
        df1=plt.figure()
        for ii in range(len(llistxyq)):
            plt.plot(llistxyq[ii][0],llistxyq[ii][1]);
        if xlb != 'none':
            plt.xlabel(xlb)
        if ylb != 'none':
            plt.ylabel(ylb)
        df1.show()
        return

    def epllxyloglog(llistxyq):
        """Shorthand plotting function for a list of multiple lists of x-values and y-values with logrithmic axes in both directions"""
        df1=plt.figure()
        for ii in range(len(llistxyq)):
            plt.loglog(llistxyq[ii][0],llistxyq[ii][1]);
        df1.show()
        return

    def epllt(listq,Map):
        formtra=[]
        for ii in range(len(listq)):
            formtra.append(LPL.TraceFormatting(listq[ii],Map[ii],1))
        efc.epll(formtra)
        return
        
    def epllcomp(listq,goalq,Map,tMax):
        """Shorthand function for comparing a list of pulse shapes to a targeted pulse shape; not used much currently"""
        formtra=[]
        formtra.append(LPL.TraceFormatting(listq[-1],Map,tMax))
        formtra.append(goalq)
        efc.epll(formtra)
        return
        
    def eplfft(errlistQ,time_stepQ):
        """Shorthand function for loglog-plotting the normalized power spectrum of the Fourier transform of a waveform, given the temporal spacing of the waveform"""
        #time_step1=0.1
        freqs1=np.fft.fftfreq(np.array(errlistQ).size, time_stepQ)
        idx1=np.argsort(freqs1)
        fftd1=np.fft.fft(errlistQ)
        ps1=np.abs(fftd1/max(fftd1))**2
        efc.eplxyloglog(freqs1[idx1],ps1[idx1])
        return [freqs1[idx1],ps1[idx1]]
        
    def reloadchk():
        """Shorthand sanity check for the current version of the code"""
        print('Last stamped: 20220201')
        
    def reloadpkg(pkgname):
        """Used to be a poor attempt at reloading packages while under development; much easier to use %run IPython commands instead"""
        #spec = importlib.util.spec_from_file_location("mecps4", "/reg/neh/operator/mecopr/mecpython/pulseshaping/mecps4.py")
        #mecps4 = importlib.util.module_from_spec(spec)
        #spec.loader.exec_module(mecps4);
        importlib.reload(pkgname)

    @staticmethod
    def tc():
        """A list of color codes used for color printing to terminal, used by cprint()"""
        colors=['ENDC','BLINK','K','R','G','Y','B','M','C','W','BK','BR','BG','BY','BB','BM','BC','BW','BRK','BRR','BRG','BRY','BRB','BRM','BRC','BRW','BBRK','BBRR','BBRG','BBRY','BBRB','BBRM','BBRC','BBRW']#B- for background, BR+ for 
        colorcodes=['\033['+str(ii)+'m' for ii in [0,5,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,90,91,92,93,94,95,96,97,100,101,102,103,104,105,106,107]];
        return dict(zip(colors,colorcodes))

    @classmethod
    def cprint(cls,strQ,paramsQ):#delim with commas
        """Prints to terminal using provided parameters for color, etc."""
        prargs=''
        if len(paramsQ) == 0:
            paramsQ='ENDC'
        for eaarg in paramsQ.split(','):
            prargs+=cls.tc()[eaarg.upper()]
        print(f"{prargs}"+strQ+f"{cls.tc()['ENDC']}")
        return

    @classmethod
    def cstr(cls,strQ,paramsQ):
        """Prepares a string (i.e. for later printing to terminal) using provided parameters for color, etc."""
        prargs=''
        if len(paramsQ) == 0:
            paramsQ='ENDC'
        for eaarg in paramsQ.split(','):
            prargs+=cls.tc()[eaarg.upper()]
        return f"{prargs}"+strQ+f"{cls.tc()['ENDC']}"

    def keybd():
        """"Prepares some keyboard input interpretation parameters needed for interpreting certain keystrokes returned by getch()"""
        return dict(zip(['key_Enter','key_Esc','key_Up','key_Dn','key_Rt','key_Lt'],[13,27,'\033[A','\033[B','\033[C','\033[D']))

    def getch():
        """Similar to input() but takes only a single character and doesn't require hitting Enter"""
        fdInput = sys.stdin.fileno()
        termAttr = termios.tcgetattr(0)
        tty.setraw(fdInput)
        ch = sys.stdin.buffer.raw.read(4).decode(sys.stdin.encoding)
        if len(ch) == 1:
            if ord(ch) < 32 or ord(ch) > 126:
                ch = ord(ch)
        elif ord(ch[0]) == 27:
            ch = '\033' + ch[1:]
        termios.tcsetattr(fdInput, termios.TCSADRAIN, termAttr)
        return ch     

    def getch_with_TO(TOsec):
        print("You have {} seconds to answer!".format(str(TOsec)))
        fdInput = sys.stdin.fileno()
        termAttr = termios.tcgetattr(0)
        tty.setraw(fdInput)
        i, o, e = select.select( [sys.stdin], [], [], TOsec)
        if (i):
            ch = sys.stdin.buffer.raw.read(4).decode(sys.stdin.encoding)
            if len(ch) == 1:
                if ord(ch) < 32 or ord(ch) > 126:
                    ch = ord(ch)
            elif ord(ch[0]) == 27:
                ch = '\033' + ch[1:]
            termios.tcsetattr(fdInput, termios.TCSADRAIN, termAttr)
            return ch     
        else:
            return False

    def input_with_timeout(TOsec):
        """Similar to input() but includes a user timeout so the window for input doesn't stay open forever and block the terminal"""
        print("You have {} second{} to answer!".format(str(TOsec), '' if TOsec==1 else 's'))
        i, o, e = select.select( [sys.stdin], [], [], TOsec)
        if (i):
            return sys.stdin.readline().strip()
        else:
            return False
    
    def pickledump2(objQ,fullFileNameQ):
        """Shortcut for generating pickle files and setting file access permissions as liberally as possible"""
        pickle.dump(objQ,open(fullFileNameQ,'wb'));
        os.chmod(fullFileNameQ,stat.S_IRUSR|stat.S_IWUSR|stat.S_IRGRP|stat.S_IWGRP|stat.S_IROTH|stat.S_IWOTH);#
        #os.chmod(fullFileNameQ,stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO);#
        return
    
    def dotsleep(tSEC):
        """Similar to time.sleep() but also prints a . character every second for the entire duration, finishing with a * character"""
        for ii in range(tSEC):
            print('.',end='',flush=True);time.sleep(1);
        print('*')
        return

    def ssleep():
        """A convenient shortcut for pausing 150ms for certain types of socket functions to wrap up"""
        time.sleep(0.15);
        return

    @classmethod                  
    def pldemo(cls):
        """An old demo meant to demonstrate how Python could be used to tune a laser based on camera feedback and EPICS actuators"""
        plt.ion() 
        fig,axs=plt.subplots(1,1) 
        plt.show()
        xpos=85;ypos=65;xdel=20;ydel=20;
        Z=[[np.exp(-((ii-xpos)/10)**2-((jj-ypos)/7.5)**2) for ii in range(100)] for jj in range(100)]
        Zref=[[1 if np.exp(-((ii-50)/2)**2-((jj-50)/2)**2) > 0.5 else 0 for ii in range(100)] for jj in range(100)]
        ax1=axs.imshow(Z,origin='lower')
        axs.imshow(Zref,alpha=0.1,origin='lower')
        #ax1,=axs[0].plot(xdat,ydat)
        #ax2,=axs[1].plot(xdat,ydat) 
        #ax4,=axs[1,1].plot(xdat,ydat) 
        axss=[ax1]#[ax1,ax2,ax3,ax4]
        cont=True
        while cont:
            axss[0].set_data(Z)
            fig.canvas.draw_idle()
            plt.pause(0.025)
            qq=cls.getch()
            if qq==cls.keybd()['key_Dn']:
                ypos-=ydel
            elif qq==cls.keybd()['key_Up']:
                ypos+=ydel
            elif qq==cls.keybd()['key_Rt']:
                xpos+=xdel
            elif qq==cls.keybd()['key_Lt']:
                xpos-=xdel
            elif qq=='w':
                ydel=ydel*2
            elif qq=='s':
                ydel=ydel/2
            elif qq=='a':
                xdel=xdel/2
            elif qq=='d':
                xdel=xdel*2
            elif qq==cls.keybd()['key_Esc']:
                cont=False
            else:
                pass
            Z=[[np.exp(-((ii-xpos)/10)**2-((jj-ypos)/7.5)**2) for ii in range(100)] for jj in range(100)]
            print('['+str(xpos)+','+str(ypos)+']')
        plt.ioff()

    def EZeLog():
        """One day this will take care of all the hassle of posting to the LCLS eLog from Python so people can do so very easily"""
        pass
    
    def rPV(yourPV):
        try:
            temppv=EpicsSignal(yourPV);
            tempval=temppv.get()
            return tempval
        except:
            print('Failed: {}!'.format(yourPV))
            return False
        
    def wPV(yourPV, yourVal):
        try:
            temppv=EpicsSignal(yourPV);
            temppv.put(yourVal)
        except:
            print('Failed:{} and {}!'.format(yourPV, yourVal))
        return 
    
# =============================================================================
#     def Thread_Function(thrEvent, interfunc, args, timeout=50, loopmax=10):
#         #definitely the wrong function, just forgot where the right one is
#         loopiter=0
#         t0=time.time()
#         while thrEvent.is_set():
#             interfunc(thrEvent, *args)
#             loopiter+=1
#             if (loopiter >= loopmax) or (t0 + timeout <= time.time()):
#                 thrEvent.clear()
# 
#     def Thread_Function3(thrEvent, interfunc, args, timeout=50, loopmax=10):
#         loopiter=0
#         t0=time.time()
#         while thrEvent.is_set():
#             interfunc(thrEvent, *args)
#             loopiter+=1
#             if (loopiter >= loopmax) or (t0 + timeout <= time.time()):
#                 thrEvent.clear()
#                 
#     def interfunc1(thrEvent, nameQ, name2Q):
#         uin=input()
#         print('{} and {}, Here it is big: {}'.format(nameQ, name2Q,uin.upper()))
#         if uin.lower() == 'q':
#             thrEvent.clear()
#         time.sleep(1)
#         
#     @classmethod
#     def threaddemo(cls):
#         running = threading.Event()
#         running.set()
#     
#         thread = threading.Thread(target=cls.Thread_Function, args=(running,'t1'))
#         thread2 = threading.Thread(target=cls.Thread_Function3, args=(running, cls.interfunc1, ('sam','dave'), 3, 200))
#     
#         thread.start()
#         time.sleep(0.5)
#         thread2.start()
#         time.sleep(0.25)
#     
#         ppp=0def CAMview(*CAMargs,ImageNo=2,LIVE=False,MAXLOOPS=10):
#     if CAMargs == ('all',):
#         CAMargs = ('Regen', 'Trap', 'StrInA', 'StrInB', 'MPA1In', 'MPA1Out', 'MPA2In', 'MPA2Out', 'MPA2Xtal', 'CompIn', 'CompOutNF', 'CompOutFF')
#     if len(CAMargs) == 1:
#         quickCAM(*CAMargs,LIVE=LIVE,MAXLOOPS=MAXLOOPS)
#         return
#     plt.ion()
#     subply=len(CAMargs)//2 + len(CAMargs)%2;subplx=2;
#     fig,axs=plt.subplots(subply,subplx,figsize=(5,2*subply));
#     axss=[];tres1L=[];tPVheadL=[]
#     for ii in range(len(CAMargs)):
#         tidx=(ii//2,ii%2) if len(CAMargs) > 2 else (ii%2)
#         axs[tidx].axes.xaxis.set_ticklabels([]);
#         axs[tidx].axes.yaxis.set_ticklabels([]);
#         axs[tidx].tick_params(direction='in'); 
#         axs[tidx].set_ylabel(CAMargs[ii]); 
#         try:
#             tPVhead=CAMname(CAMargs[ii])
#             tPVheadL.append(tPVhead);
#             tres1=rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
#             tres2=rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArraySize0_RBV')
#             tres3=rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
#             tresL=sorted([tres1,tres2,tres3],reverse=True) 
#             twf=rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArrayData')
#             if len(twf) != tresL[0]*tresL[1]:
#                 twf = list(twf) + (tresL[0]*tresL[1] - len(twf))*[0]
#                 tres1=tres3
#             tempax=axs[tidx].imshow(np.array_split(np.array(twf),tres1)); 
#             tres1L.append(tres1);
#             axss.append(tempax)
#         except:
#             print('Error occured when plotting {}!'.format(CAMargs[ii]))
#     if (len(CAMargs) > 2) and (len(CAMargs)%2 > 0):
#         iit=len(CAMargs)
#         tidx=(iit//2,iit%2) 
#         axs[tidx].axis('off');
#     fig.tight_layout()
#     plt.show();
#     waittime=.01;
#     plt.pause(waittime);
#     time.sleep(waittime)
#     loopcount=0
#     if LIVE:
#         while loopcount<MAXLOOPS:
#             for ii in range(len(CAMargs)):
#                 tidx=(ii//2,ii%2) if len(CAMargs) > 2 else (ii%2)
#                 try:
#                     twf=rPV(tPVheadL[ii]+':IMAGE'+str(ImageNo)+':ArrayData')
#                     if tPVheadL[ii] == 'MEC:GIGE:29':
#                         tres1=rPV(tPVheadL[ii]+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
#                         tres2=rPV(tPVheadL[ii]+':IMAGE'+str(ImageNo)+':ArraySize0_RBV')
#                         tres3=rPV(tPVheadL[ii]+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
#                         tresL=sorted([tres1,tres2,tres3],reverse=True) 
#                         if len(twf) != tresL[0]*tresL[1]:
#                             twf = list(twf) + (tresL[0]*tresL[1] - len(twf))*[0]
#                             tres1=tres3
#                         tres1L[ii]=tres1
#                     axss[ii].set_data(np.array_split(np.array(twf),tres1L[ii])); 
#                     fig.canvas.draw_idle()
#                     plt.pause(waittime)
#                     time.sleep(waittime)
#                 except:
#                     print('Error occured when plotting {}!'.format(CAMargs[ii]))
#             loopcount+=1
#         while running.is_set():
#             print('Waiting...')
#             time.sleep(1)
#             ppp+=1
#             if ppp > 20:
#                 running.clear()
#             
#         print('Wait until Thread is terminating')
#         thread.join()
#         thread2.join()
#         print("EXIT __main__")
# =============================================================================

    def reprintdemo():
        for x in range(10):
            print('{:>5}'.format(x*10**int(5*np.random.rand())), end='\r');time.sleep(1);
            print()





class HAWG:
    """Class containing all the necessary functions for running the Highland Arbitrary Waveform Generator for LPL pulse shaping"""
    def __init__(self):
        """Initializes the object; only one should be instantiated at a time"""
        self._HighlandSocket=None
        self._HIGHLAND_SLAVE_ADDRESS=0 #arbitrary, I think    

    def Open(self):
        """Takes care of opening the socket to the Highland; if called explicitly like this, it MUST be followed by a Close() statement or else you'll block the socket and need to power cycle the unit"""
        if not self._HighlandSocket:
            try:
                self._HighlandSocket=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                self._HighlandSocket.settimeout(1.0)
                self._HighlandSocket.connect(('highland-mec-01', 2000))
                #Highland's IP address can be changed using the Lantronix DeviceInstaller
            except:
                print('HIGHLAND NOT CONNECTED')
        else:
            print('Socket may already be open!')

    def Close(self):
        """Takes care of closing the socket to the Highland"""
        try:
            self._HighlandSocket.close()
            self._HighlandSocket=None
        except:
            print('Unable to close socket -- it may already be closed')
        return
        
    def Reset(self):
        """Power-cycles the Highland; this is necessary if the socket was broken"""
        print('Powering off Highland AWG, waiting 3sec...',end='',flush=True);
        GLOBAL.LPLHAWGpwr.put(2);
        efc.dotsleep(3);
        print('Rebooting Highland AWG, waiting 10sec...',end='',flush=True)
        GLOBAL.LPLHAWGpwr.put(1);
        efc.dotsleep(10);
        
    @staticmethod
    def _Hex1Byte(num):
        """Internal function for formatting a number as a single byte; used for communicating with the Highland"""
        return '{0:02x}'.format(int(num)%(0xff+1))
        
    @staticmethod
    def _Hex2Byte(num):
        """Internal function for formatting a number as two bytes; used for communicating with the Highland"""
        return '{0:04x}'.format(int(num)%(0xffff+1))

    @staticmethod
    def _ByteSum(Datastr):
        """Internal byte-wise addition used for the Highland checksum
        accepts string in 'xx' format, e.g. '1e....'
        input must be in format returned by hexlify(data)"""
        bytesum=0
        for byte in range(len(Datastr)//2):
            bytesum+=int(Datastr[2*byte:2*byte+2],16)
        return bytesum #returns an integer

    @classmethod
    def _PollConstructor(cls,COMMAND_CODE,POLL_LENGTH,SLAVE_ADDRESS,DATA):
        """Generic constructor for all polls to the internal Highland processor
        enter COMMAND_CODE, POLL_LENGTH, SLAVE_ADDRESS as integer values
        enter DATA as a string, e.g. 'ffff0000' or empty string '' for no data"""
        ProtoCommand=''
        ProtoCommand+='0B'#MSYN: all commands begin with this byte (0B)
        ProtoCommand+=cls._Hex2Byte(POLL_LENGTH) #BC:BC: number of bytes in message
        ProtoCommand+=cls._Hex2Byte(SLAVE_ADDRESS) #RA:RA: slave address
        ProtoCommand+=cls._Hex1Byte(COMMAND_CODE) #CMD: command code
        ProtoCommand+=DATA #<data> must already be formatted properly 'xxxx' or ''
        BYTE_SUM=cls._ByteSum(ProtoCommand) #compute the sum
        ProtoCommand+=cls._Hex2Byte(BYTE_SUM) #CKS:CKS: 16-bit sum of all preceding bytes
        ProtoCommand+='17' #ETB: end of message byte, 17 hex
        Command=unhexlify(ProtoCommand)
        return Command

    @classmethod
    def _ReplyInterpreter(cls,REPLY_LENGTH,SLAVE_ADDRESS,REPLY_STRING):
        """Interpreter for replies from the internal Highland processor
        input REPLY_STRING already formatted using hexlify"""
        HError=''
        if int(REPLY_STRING[0:2],16)==int('1e',16): HError+='0'
        else: HError+='1' #wrong start-of-message byte, 1E hex
        if int(REPLY_STRING[2:6],16)==REPLY_LENGTH: HError+='0'
        else: HError+='1' #wrong reply length; should never happen, as we recv(expected #)
        if int(REPLY_STRING[6:10],16)==SLAVE_ADDRESS: HError+='0'
        else: HError+='1' #slave address not echoed
        HStatus=REPLY_STRING[10:12] #will return status as string, interpret later
        HData=REPLY_STRING[12:-6] #cuts off SSYN,BC:BC,RA:RA,STS and CKS:CKS,ETB bytes
        # leaves only the data string; leaves empty string '' for no data
        if cls._ByteSum(REPLY_STRING[:-6])==int(REPLY_STRING[-6:-2],16): HError+='0'
        else: HError+='1' #checksum error
        if int(REPLY_STRING[-2:],16)==int('17',16): HError+='0'
        else: HError+='1' #wrong end-of-message byte, 17 hex
        return HStatus, HData, HError

    @classmethod
    def _SendPollRecvReply(cls,MySocketQ,COMMAND_CODE,POLL_LENGTH,REPLY_LENGTH,SLAVE_ADDRESS,DATA):
        """Generic utility for sending a poll to the internal Highland processor and receiving its reply"""
        MyPollQ=cls._PollConstructor(COMMAND_CODE,POLL_LENGTH,SLAVE_ADDRESS,DATA)
        MySocketQ.send(MyPollQ)
        MyRawReplyQ=MySocketQ.recv(REPLY_LENGTH)
        HStatusQ, HDataQ, HErrorQ = cls._ReplyInterpreter(REPLY_LENGTH,SLAVE_ADDRESS,hexlify(MyRawReplyQ))
        return HStatusQ, HDataQ, HErrorQ

    @staticmethod
    def _StatusInterpreter(HError, HStatus, Quiet=True):
        """Interpreter for the part of the Highland processor response indicating internal status"""
        if HError[0]=='1': print('WARNING: Wrong start-of-message byte received')
        if HError[1]=='1': print('WARNING: Reply length discrepancy')
        if HError[2]=='1': print('WARNING: Slave address not echoed')
        if HError[3]=='1': print('WARNING: Checksum error')
        if HError[4]=='1': print('WARNING: Wrong end-of-message byte received')
        if not Quiet:
            if int(HStatus,16)==0: print('STATUS: NORMAL')
            else:
                print('STATUS: ERROR FLAG(S) RECEIVED')
                if ((int(HStatus,16))&(2**(8-1)))!=0: print('-trigger/bias timing error')
                if ((int(HStatus,16))&(2**(8-3)))!=0: print('-backup RAM data/calibrations lost flag')
                if ((int(HStatus,16))&(2**(8-4)))!=0: print('-powerfail/restart flag')
                if ((int(HStatus,16))&(2**(8-7)))!=0: print('-trigger/bias timing error')

    @classmethod
    def _BaseFunction(cls,MySocketQ, SlaveAddress, CommandCode, PollLength, ReplyLength, MyData='',**kwargs):
        """Generic function used to construct the different individual commands to the Highland;
        each individual command simply supplies its own unique code and data specific to that command's protocol"""
        try:
            HStatusQ, HDataQ, HErrorQ = cls._SendPollRecvReply(
            MySocketQ,CommandCode,PollLength,ReplyLength,SlaveAddress,MyData)
            cls._StatusInterpreter(HErrorQ, HStatusQ)
            return HDataQ
        except:
            print('Failed!')
            return False

    def _FunctionWrapper(self,FuncQ,kwargs={}):
        """A function wrapper that allows one to call each Highland command directly without having to worry about opening/closing sockets;
        if issuing one-off commands that don't require high-frequency consistent execution, this is sufficient"""
        try:
            self.Open();time.sleep(0.15);
            HDataQList=FuncQ(**kwargs);time.sleep(0.15);
            self.Close();time.sleep(0.15);
        except:
            self.Close();time.sleep(0.15);
        return HDataQList

    def _ReadStatus(self):
        """
        Bare function for the ReadStatus command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        POLL command 0 (with no <data> field) invokes a 'status data' reply from a slave. 
        The reply message contains a data field having the following subfields... 
        
        PROGRAM ID 8-byte ASCII firmware ID/revision field. 
         
        UPTIME 4-byte uptime, as 32-bit value, in seconds; 
        cleared at powerup time. 
         
        ENABLE 1-byte field identifying subsystems which 
        are enabled. See 'WRITE ENABLE 
        COMMAND' below for bit assignments. 
        """
        print('**READ STATUS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=0,PollLength=9, ReplyLength=22,MyData='')
        print('PROGRAM ID: '+unhexlify(HDataQ[:16]).decode())
        print('UPTIME: '+str(int(HDataQ[16:24],16))+' seconds')
        print('ENABLE:')
        if ((int(HDataQ[-2:],16))&(2**(0)))!=0:
            print('-CPU self-trigger, 960 Hz (test mode)')
        if ((int(HDataQ[-2:],16))&(2**(1)))!=0:
            print('-self-trigger, 20 kHz')
        if ((int(HDataQ[-2:],16))&(2**(2)))!=0:
            print('-external triggers')
        if ((int(HDataQ[-2:],16))&(2**(3)))!=0:
            print('-the BIAS generators')
        print('****')
        return HDataQ
    def ReadStatus(self):
        """
        Wrapped function for the ReadStatus command; socket is automatically opened/closed
        ...
        From the T400B manual:
        POLL command 0 (with no <data> field) invokes a 'status data' reply from a slave. 
        The reply message contains a data field having the following subfields... 
        
        PROGRAM ID 8-byte ASCII firmware ID/revision field. 
         
        UPTIME 4-byte uptime, as 32-bit value, in seconds; 
        cleared at powerup time. 
         
        ENABLE 1-byte field identifying subsystems which 
        are enabled. See 'WRITE ENABLE 
        COMMAND' below for bit assignments. 
        """
        HDataQ=self._FunctionWrapper(self._ReadStatus);
        return HDataQ

    def _ClearStatus(self):
        """
        Bare function for the ClearStatus command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        Execution of this command clears the slave STATUS byte. If any error conditions
        persist, status error bits may reappear immediately. 
        """
        print('**CLEAR STATUS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=1,PollLength=9, ReplyLength=9,MyData='')
        print('****')
        return HDataQ
    def ClearStatus(self):
        """
        Wrapped function for the ClearStatus command; socket is automatically opened/closed
        ...
        From the T400B manual:
        Execution of this command clears the slave STATUS byte. If any error conditions
        persist, status error bits may reappear immediately. 
        """        
        HDataQ=self._FunctionWrapper(self._ClearStatus);
        return HDataQ

    def _ReadPulseHeights(self,ShowPlot=False):
        """
        Bare function for the ReadPulseHeights command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        Poll command 2 (no <data> field) reads the 140 programmable waveform impulse 
        heights, each an unsigned 16-bit value (0 is zero height, 65535 is max). The slave 
        thus returns a 280-byte <data> field. The first two <data> bytes are the programmed 
        pulse height for the first 250 ps impulse segment, sent MS byte first.   
        """
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=2,PollLength=9, ReplyLength=289,MyData='')
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(int(len(HDataQ)//4))]
        if ShowPlot:
            efc.epl(HDataQList)
        return HDataQList
    def ReadPulseHeights(self,ShowPlot=False):
        """
        Wrapped function for the ReadPulseHeights command; socket is automatically opened/closed
        ...
        From the T400B manual:
        Poll command 2 (no <data> field) reads the 140 programmable waveform impulse 
        heights, each an unsigned 16-bit value (0 is zero height, 65535 is max). The slave 
        thus returns a 280-byte <data> field. The first two <data> bytes are the programmed 
        pulse height for the first 250 ps impulse segment, sent MS byte first.   
        """      
        HDataQList=self._FunctionWrapper(self._ReadPulseHeights,{'ShowPlot':ShowPlot});
        return HDataQList
        
    def _WritePulseHeights(self, FileNameOrStringOrList=140*[0]):
        """
        Bare function for the WritePulseHeights command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        Poll command 3 writes the 140 programmable pulse heights, each an unsigned  
        16-bit value. The poll <data> field is 280 bytes in size, with the first two bytes being 
        the programmed pulse height for the first 250 ps impulse segment, sent MS byte 
        first. The slave responds with the standard 9-byte ACK (no data) message. 
        """
        MyDataQ=''
        if len(FileNameOrStringOrList) == 140*4:#will accept pre-formatted Hex2Byte text
            MyDataQ=FileNameOrStringOrList
        elif len(FileNameOrStringOrList)==140:#will accept a straight list
            for value in range(len(FileNameOrStringOrList)):
                MyDataQ+=self._Hex2Byte(int(FileNameOrStringOrList[value]))
        elif FileNameOrStringOrList.endswith(('.txt','.csv','.dat')):
            with open(FileNameOrStringOrList,'r') as filehead:
                RawListQ=filehead.read()
                if '\r\n' in RawListQ:
                    ListedValues=RawListQ.split('\r\n')
                elif '\n' in RawListQ:
                    ListedValues=RawListQ.split('\n')
                elif ',' in RawListQ:
                    ListedValues=RawListQ.split(',')
                else:
                    print('Unrecognized format on input file.')
                    return False
            if len(ListedValues) != 140:
                print('File must have 140 entries; entry count: '+str(len(ListedValues)))
                return False
            for value in range(len(ListedValues)):
                MyDataQ+=self._Hex2Byte(int(ListedValues[value]))
        else:
            print('Bad file entry count: '+str(len(FileNameOrStringOrList)))
            return False
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=3,PollLength=289, ReplyLength=9,MyData=MyDataQ)
        return HDataQ
    def WritePulseHeights(self,FileNameOrStringOrList=140*[0]):
        """
        Wrapped function for the WritePulseHeights command; socket is automatically opened/closed
        ...
        From the T400B manual:
        Poll command 3 writes the 140 programmable pulse heights, each an unsigned  
        16-bit value. The poll <data> field is 280 bytes in size, with the first two bytes being 
        the programmed pulse height for the first 250 ps impulse segment, sent MS byte 
        first. The slave responds with the standard 9-byte ACK (no data) message. 
        """  
        HDataQ=self._FunctionWrapper(self._WritePulseHeights,{'FileNameOrStringOrList':FileNameOrStringOrList});
        return HDataQ
        
    def _ReadFiducialImpulseSettings(self):
        """
        Bare function for the ReadFiducialImpulseSettings command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        Poll command 4 (no <data> field) reads two unsigned integers: the  16-bit value 
        which determines the height of the auxiliary fiducial impulse generator and a second 
        16-bit value which determines its delay. Each value is returned in a two-byte reply 
        data field. Values range from 0 to +65535. The LSB of the time setting is 1 ps. If the 
        first value (impulse amplitude) is zero, the impulse circuitry will be disabled. 
        
        The fiducial impulse is of fixed width (nom 100 ps) and is summed into the main 
        140-point modulator waveform. 
        """
        print('**READ FIDUCIAL IMPULSE SETTINGS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=4,PollLength=9, ReplyLength=13,MyData='')
        print('Fiducial pulse height (max 65535): '+str(int(HDataQ[:4],16)))
        print('Fiducial pulse delay: '+str(int(HDataQ[4:8],16)))
        print('****')
        return HDataQ
    def ReadFiducialImpulseSettings(self):
        """
        Wrapped function for the ReadFiducialImpulseSettings command; socket is automatically opened/closed
        ...
        From the T400B manual:
        Poll command 4 (no <data> field) reads two unsigned integers: the  16-bit value 
        which determines the height of the auxiliary fiducial impulse generator and a second 
        16-bit value which determines its delay. Each value is returned in a two-byte reply 
        data field. Values range from 0 to +65535. The LSB of the time setting is 1 ps. If the 
        first value (impulse amplitude) is zero, the impulse circuitry will be disabled. 
        
        The fiducial impulse is of fixed width (nom 100 ps) and is summed into the main 
        140-point modulator waveform. 
        """  
        HDataQ=self._FunctionWrapper(self._ReadFiducialImpulseSettings);
        return HDataQ        
        
    def _WriteFiducialImpulseSettings(self,AmpReq=0,TimeReq=0):
        """
        Bare function for the WriteFiducialImpulseSettings command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        Poll command 5 writes the fiducial impulse settings, described above. The 4-byte 
        <data> field sets the impulse height and delay. 
        """
        print('**WRITE FIDUCIAL IMPULSE SETTINGS**')
        MyDataQ=''
        MyDataQ+=self._Hex2Byte(int(AmpReq))
        MyDataQ+=self._Hex2Byte(int(TimeReq))
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=5,PollLength=13, ReplyLength=9,MyData=MyDataQ)
        print('****')
        return HDataQ
    def WriteFiducialImpulseSettings(self,AmpReq=0,TimeReq=0):
        """
        Wrapped function for the WriteFiducialImpulseSettings command; socket is automatically opened/closed
        ...
        From the T400B manual:
        Poll command 5 writes the fiducial impulse settings, described above. The 4-byte 
        <data> field sets the impulse height and delay. 
        """  
        HDataQ=self._FunctionWrapper(self._WriteFiducialImpulseSettings,{'AmpReq':AmpReq,'TimeReq':TimeReq});
        return HDataQ
        
    def _WriteEnableByte(self,EnableTotal=4):
        """
        Bare function for the WriteEnableByte command; socket must be explicitly opened/closed
        to determine EnableTotal input, start from 0 and:
        +1 for Enable CPU self-trigger, 960 Hz (test mode)
        +2 for Enable self-trigger, 20 kHz
        +4 for Enable external triggers
        +8 for Enable the BIAS generators
        ...
        From the T400B manual:
        The poll message <data> field contains the single ENABLE byte. Bits in this byte 
        enable/disable subsystems. Bits are... 
        
         BIT FUNCTION WHEN SET  
         0 Enable CPU self-trigger, 960 Hz (test mode) 
         1 Enable self-trigger, 20 KHz 
         2 Enable external triggers 
         3 Enable the BIAS generators 
        """
        MyDataQ=''
        MyDataQ+=self._Hex1Byte(EnableTotal)
        print('**WRITE ENABLE BYTE**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=10,PollLength=10, ReplyLength=9,MyData=MyDataQ)
        print('****')
        return HDataQ
    def WriteEnableByte(self,EnableTotal=4):
        """
        Wrapped function for the WriteEnableByte command; socket is automatically opened/closed
        to determine EnableTotal input, start from 0 and:
        +1 for Enable CPU self-trigger, 960 Hz (test mode)
        +2 for Enable self-trigger, 20 kHz
        +4 for Enable external triggers
        +8 for Enable the BIAS generators
        ...
        From the T400B manual:
        The poll message <data> field contains the single ENABLE byte. Bits in this byte 
        enable/disable subsystems. Bits are... 
        
         BIT FUNCTION WHEN SET  
         0 Enable CPU self-trigger, 960 Hz (test mode) 
         1 Enable self-trigger, 20 KHz 
         2 Enable external triggers 
         3 Enable the BIAS generators 
        """
        HDataQ=self._FunctionWrapper(self._WriteEnableByte,{'EnableTotal':EnableTotal});
        return HDataQ
        
    def _ReadT0Delay(self):
        """
        Bare function for the ReadT0Delay command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        The <data> field returns the current 16-bit T0 delay value, with decimal value 0 
        corresponding to minimum delay, with LSB resolution of 1 ps nominal. The 
        maximum legal value is 50,000, corresponding to 50 ns nominal delay. This is the 
        delay applied to the T400 trigger, and shifts all other timed events.  
        
        The earliest possible square pulse or fiducial pulse edge will occur 20 ns after the 
        T0 delay setting, and the first arb glitch will be centered at T0 + 25 ns. 
        """
        print('**READ T0 DELAY**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=13,PollLength=9, ReplyLength=11,MyData='')
        print('T0 delay (max 50000 (50ns)): '+str(int(HDataQ,16)))
        print('****')
        return int(HDataQ,16)
    def ReadT0Delay(self):
        """
        Wrapped function for the ReadT0Delay command; socket is automatically opened/closed
        ...
        From the T400B manual:
        The <data> field returns the current 16-bit T0 delay value, with decimal value 0 
        corresponding to minimum delay, with LSB resolution of 1 ps nominal. The 
        maximum legal value is 50,000, corresponding to 50 ns nominal delay. This is the 
        delay applied to the T400 trigger, and shifts all other timed events.  
        
        The earliest possible square pulse or fiducial pulse edge will occur 20 ns after the 
        T0 delay setting, and the first arb glitch will be centered at T0 + 25 ns. 
        """
        HDataQ=self._FunctionWrapper(self._ReadT0Delay);
        return HDataQ
        
    def _ReadWaveAmplitudeCalibrations(self):
        """
        Bare function for the ReadWaveAmplitudeCalibrations command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        This command reads the waveform impulse calibration table, a list of 140 integers 
        corresponding to the 140 impulses that form the arbitrary waveform. Each integer 
        scales the amplitude of its corresponding impulse. The nominal value of each integer 
        is about 2800 decimal. The T400 is factory calibrated, so these values should not 
        need to be altered unless a segment board is replaced. 
        """
        print('**READ WAVE AMPLITUDE CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=15,PollLength=9, ReplyLength=289,MyData='')
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
        print('Wave amplitude calibration (nom 2800): '+str(HDataQList))
        print('****')
        return HDataQList
    def ReadWaveAmplitudeCalibrations(self):
        """
        Wrapped function for the ReadWaveAmplitudeCalibrations command; socket is automatically opened/closed
        ...
        From the T400B manual:
        This command reads the waveform impulse calibration table, a list of 140 integers 
        corresponding to the 140 impulses that form the arbitrary waveform. Each integer 
        scales the amplitude of its corresponding impulse. The nominal value of each integer 
        is about 2800 decimal. The T400 is factory calibrated, so these values should not 
        need to be altered unless a segment board is replaced. 
        """
        HDataQList=self._FunctionWrapper(self._ReadWaveAmplitudeCalibrations);
        return HDataQList
        
    def _WriteWaveAmplitudeCalibrations(self, StringOrList):
        """
        Bare function for the WriteWaveAmplitudeCalibrations command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        (nothing listed; see description for ReadWaveAmplitudeCalibrations)
        """
        MyDataQ=''
        if len(StringOrList) == 140*4:#will accept pre-formatted Hex2Byte text
            MyDataQ=StringOrList
        elif len(StringOrList)==140:#will accept a straight list
            for value in range(len(StringOrList)):
                MyDataQ+=self._Hex2Byte(int(StringOrList[value]))
        else:
            print('Bad file entry count: '+str(len(StringOrList)))
            return 
        print('**WRITE WAVE AMPLITUDE CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=16,PollLength=289, ReplyLength=9,MyData=MyDataQ)
        print('****')
        return HDataQ
    def WriteWaveAmplitudeCalibrations(self,StringOrList):
        """
        Wrapped function for the WriteWaveAmplitudeCalibrations command; socket is automatically opened/closed
        ...
        From the T400B manual:
        (nothing listed; see description for ReadWaveAmplitudeCalibrations)
        """
        HDataQ=self._FunctionWrapper(self._WriteWaveAmplitudeCalibrations,{'StringOrList':StringOrList});
        return HDataQ
        
    def _ReadWaveTimeCalibrations(self):
        """
        Bare function for the ReadWaveTimeCalibrations command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        This command reads the waveform glitch timing table, a list of 140 integers 
        corresponding to the time of each of the 140 glitches that form the arbitrary 
        waveform. Each of the ten segment boards generates 14 glitches which must be 
        spaced 250 ps apart. The table is thus organized as ten sets of 14 integers, with 
        each set having approximately the same ascending pattern, from about 1000 to 
        2755 in 14 steps of about 135, with LSB weight of about 1.8 ps. 
        
        In order to ensure the best possible arbitrary waveform matching, users should 
        consider performing occasional online recalibration of the glitch spacings; drifts of as 
        little as 20 ps can result in significant deviations of synthesized waveforms.
        """
        print('**READ WAVE TIME CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=17,PollLength=9, ReplyLength=289,MyData='')
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
        print('Wave time calibrations (max 65535): '+str(HDataQList))
        print('****')
        return HDataQList
    def ReadWaveTimeCalibrations(self):
        """
        Wrapped function for the ReadWaveTimeCalibrations command; socket is automatically opened/closed
        ...
        From the T400B manual:
        This command reads the waveform glitch timing table, a list of 140 integers 
        corresponding to the time of each of the 140 glitches that form the arbitrary 
        waveform. Each of the ten segment boards generates 14 glitches which must be 
        spaced 250 ps apart. The table is thus organized as ten sets of 14 integers, with 
        each set having approximately the same ascending pattern, from about 1000 to 
        2755 in 14 steps of about 135, with LSB weight of about 1.8 ps. 
        
        In order to ensure the best possible arbitrary waveform matching, users should 
        consider performing occasional online recalibration of the glitch spacings; drifts of as 
        little as 20 ps can result in significant deviations of synthesized waveforms.
        """
        HDataQList=self._FunctionWrapper(self._ReadWaveTimeCalibrations);
        return HDataQList  
        
    def _WriteWaveTimeCalibrations(self, StringOrList):
        """
        Bare function for the WriteWaveTimeCalibrations command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        (nothing listed; see description for ReadWaveTimeCalibrations)
        """
        MyDataQ=''
        if len(StringOrList) == 140*4:#will accept pre-formatted Hex2Byte text
            MyDataQ=StringOrList
        elif len(StringOrList)==140:#will accept a straight list
            for value in range(len(StringOrList)):
                MyDataQ+=self._Hex2Byte(int(StringOrList[value]))
        else:
            print('Bad file entry count: '+str(len(StringOrList)))
            return 
        print('**WRITE WAVE TIME CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=18,PollLength=289, ReplyLength=9,MyData=MyDataQ)
        print('****')
        return HDataQ
    def WriteWaveTimeCalibrations(self,StringOrList):
        """
        Wrapped function for the WriteWaveTimeCalibrations command; socket is automatically opened/closed
        ...
        From the T400B manual:
        (nothing listed; see description for ReadWaveTimeCalibrations)
        """  
        HDataQ=self._FunctionWrapper(self._WriteWaveTimeCalibrations,{'StringOrList':StringOrList});
        return HDataQ

    def _ReadMiscellaneousCalibrations(self):
        """
        Bare function for the ReadMiscellaneousCalibrations command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        This command reads a list of 36 unsigned 16-bit integers, referred to as MC0 
        through MC35. The first 20 are paired (coarse, fine) time settings for the positions in 
        time of the ten waveform segment boards, which must be spaced 3.5 ns apart. The 
        remaining cal factors are described below. 
        
        MC0-19 paired coarse:fine board delays 
        
        MC20 T0 delay time slope 
        MC21 T0 delay base time offset 
        
        MC22 Square pulse amplitude calibrator 
        MC23 Square pulse amplitude baseline 
        MC24 Square pulse delay time slope cal 
        MC25 Square pulse base time offset 
        
        MC26 Fiducial impulse amplitude calibrator 
        MC27 Fiducial impulse amplitude baseline 
        MC28 Fiducial impulse time slope 
        MC29 Fiducial impulse base time offset 
        
        MC30-35 Spares 
        """
        print('**READ MISCELLANEOUS CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=19,PollLength=9, ReplyLength=81,MyData='')
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
        print('Miscellaneous calibrations: '+str(HDataQList))
        print('****')
        return HDataQList
    def ReadMiscellaneousCalibrations(self):
        """
        Wrapped function for the ReadMiscellaneousCalibrations command; socket is automatically opened/closed
        ...
        From the T400B manual:
        This command reads a list of 36 unsigned 16-bit integers, referred to as MC0 
        through MC35. The first 20 are paired (coarse, fine) time settings for the positions in 
        time of the ten waveform segment boards, which must be spaced 3.5 ns apart. The 
        remaining cal factors are described below. 
        
        MC0-19 paired coarse:fine board delays 
        
        MC20 T0 delay time slope 
        MC21 T0 delay base time offset 
        
        MC22 Square pulse amplitude calibrator 
        MC23 Square pulse amplitude baseline 
        MC24 Square pulse delay time slope cal 
        MC25 Square pulse base time offset 
        
        MC26 Fiducial impulse amplitude calibrator 
        MC27 Fiducial impulse amplitude baseline 
        MC28 Fiducial impulse time slope 
        MC29 Fiducial impulse base time offset 
        
        MC30-35 Spares 
        """  
        HDataQList=self._FunctionWrapper(self._ReadMiscellaneousCalibrations);
        return HDataQList  
            
    def _WriteMiscellaneousCalibrations(self,StringOrList):
        """
        Bare function for the WriteMiscellaneousCalibrations command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        (nothing listed; see description for ReadMiscellaneousCalibrations)
        """
        MyDataQ=''
        if len(StringOrList) == 36*4:#will accept pre-formatted Hex2Byte text
            MyDataQ=StringOrList
        elif len(StringOrList)==36:#will accept a straight list
            for value in range(len(StringOrList)):
                MyDataQ+=self._Hex2Byte(int(StringOrList[value]))
        else:
            print('Bad file entry count: '+str(len(StringOrList)))
            return
        print('**WRITE MISCELLANEOUS CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=20,PollLength=81, ReplyLength=9,MyData=MyDataQ)
        print('****')
        return HDataQ
    def WriteMiscellaneousCalibrations(self,StringOrList):
        """
        Wrapped function for the WriteMiscellaneousCalibrations command; socket is automatically opened/closed
        ...
        From the T400B manual:
        (nothing listed; see description for ReadMiscellaneousCalibrations)
        """  
        HDataQ=self._FunctionWrapper(self._WriteMiscellaneousCalibrations,{'StringOrList':StringOrList});
        return HDataQ
    
    def _ReadWalkTable(self):
        """
        Bare function for the ReadWalkTable command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        This command reads the time walk compensation table. This table contains 
        32 unsigned integers which compensate for a small interaction between arb glitch 
        heights and glitch centroid timing. 
        """
        print('**READ WALK TABLE**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=25,PollLength=9, ReplyLength=73,MyData='')
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
        print('Walk Table: '+str(HDataQList))
        print('****')
        return HDataQList    
    def ReadWalkTable(self):
        """
        Wrapped function for the ReadWalkTable command; socket is automatically opened/closed
        ...
        From the T400B manual:
        This command reads the time walk compensation table. This table contains 
        32 unsigned integers which compensate for a small interaction between arb glitch 
        heights and glitch centroid timing. 
        """  
        HDataQList=self._FunctionWrapper(self._ReadWalkTable);
        return HDataQList  
        
    def _WriteWalkTable(self,StringOrList):
        """
        Bare function for the WriteWalkTable command; socket must be explicitly opened/closed
        ...
        From the T400B manual:
        (nothing listed; see description for ReadWalkTable)
        """
        MyDataQ=''
        if len(StringOrList) == 32*4:#will accept pre-formatted Hex2Byte text
            MyDataQ=StringOrList
        elif len(StringOrList)==32:#will accept a straight list
            for value in range(len(StringOrList)):
                MyDataQ+=self._Hex2Byte(int(StringOrList[value]))
        else:
            print('Bad file entry count: '+str(len(StringOrList)))
            return
        print('**WRITE WALK TABLE**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=26,PollLength=73, ReplyLength=9,MyData=MyDataQ)
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
        print('****')
        return HDataQList
    def WriteWalkTable(self,StringOrList):
        """
        Wrapped function for the WriteWalkTable command; socket is automatically opened/closed
        ...
        From the T400B manual:
        (nothing listed; see description for ReadWalkTable)
        """
        HDataQ=self._FunctionWrapper(self._WriteWalkTable,{'StringOrList':StringOrList});
        return HDataQ
    
    def IndFETWave(self,ListOfPixels,WriteValue):#only used in FET survey?? put inside HAWG??
        """Primarily used to write patterns of individual pixels on the Highland for the sake of calibration"""
        itt=0
        NewString=''
        while itt<140:
            if (itt+1) in ListOfPixels:
                NewString+=self._Hex2Byte(WriteValue)
            else:
                NewString+='0000'
            itt+=1
        return NewString
    
    def FidOn(self):
        """Shortcut function for turning on the Highland's fiducial impulse used for proper timing of the oscilloscope delay"""
        try:
            self.WriteFiducialImpulseSettings(20000,45000);
        except:
            print('Error!')
            return False
        return

    def FidOff(self):
        """Shortcut function for turning off the Highland's fiducial impulse"""
        try:
            self.WriteFiducialImpulseSettings(0,0);
        except:
            print('Error!')
            return False
        return


# #FIX THIS!!!
#     def ___FETsurvey(HSQ,LS1Q):
#         #HSQ=HOpen()
#         #time.sleep(.15)
#         #LS1Q=LOpen()
#         #time.sleep(.15)
#         qdatalist=[]
#         for ii in range(140):
#             WritePulseHeights(HSQ,0,IndFETWave([ii+1],28000))#could improve by doing a few spread-out points at a time
#             time.sleep(6)
#             qpixnodata=readchan(1,LS1Q)['DATA']
#             qpixnomax=max(qpixnodata)
#             qpixnomaxindex=np.mean([i for i,j in enumerate(qpixnodata) if j == qpixnomax])##could improve by changing to abbreviated centroid around peak, avoiding tail-end bump
#             qdatalist.append([qpixnomaxindex,qpixnomax])
#         #time.sleep(.15)
#         #HClose(HSQ)
#         #time.sleep(.15)
#         #LClose(LS1Q)
#         #time.sleep(.15)
#         return qdatalist
# =============================================================================
# =============================================================================
# #FIX THIS!!!
#     def ___FastFETSurvey(HSQ,LS1Q):
#         #HSQ=HOpen()
#         #time.sleep(.15)
#         #LS1Q=LOpen()
#         #time.sleep(.15)
#         qdatalist=[]
#         for ii in range(10):
#             WritePulseHeights(HSQ,0,IndFETWave([jj*10+ii+1 for jj in range(14)],28000))#could improve by doing a few spread-out points at a time
#             time.sleep(6)
#             qpixnodata=readchan(1,LS1Q)['DATA']
#             qpeakcoords=signal.find_peaks_cwt(np.clip(qpixnodata,max(qpixnodata)/5,max(qpixnodata)),np.arange(180,200))#threshold->clip
#             #this is 475ps/(2.5ps/pix)=190pix expected for S1 scope at max sampling; S2 scope needs different
#             if len(qpeakcoords) != 15:
#                 print('\nWrong number of peaks detected!\n',len(qpeakcoords))
#                 return
#             qdatalist.append([qpeakcoords[1:],qpixnodata[qpeakcoords[1:]]])
#         #time.sleep(.15)
#         #HClose(HSQ)
#         #time.sleep(.15)
#         #LClose(LS1Q)
#         #time.sleep(.15)
#         qdatalist2=np.array(qdatalist).transpose(2,0,1).reshape(140,2)
#         qTimeErrInHpix=np.array([qdatalist2[ii,0]-np.mean(qdatalist2[:14,0])+100*(6.5-ii) for ii in range(140)])*2.5/1.8
#         qTimeErrInHpixBoardAvg=np.array([np.mean(qTimeErrInHpix[14*ii:14*ii+14]) for ii in range(10)])
#         epl(qTimeErrInHpix)
#         epl(qTimeErrInHpixBoardAvg)
#         return np.array([qTimeErrInHpix,qTimeErrInHpixBoardAvg])
# =============================================================================
# =============================================================================
# #FIX THIS!!!
#     def ___VeryFastFETSurvey(HSQ,LS1Q):#returns np.array([qTimeErrInHpix,qTimeErrInHpixBoardAvg])
#         #HSQ=HOpen()
#         #time.sleep(.15)
#         #LS1Q=LOpen()
#         #time.sleep(.15)
#         qdatalist=[]
#         qrfis=ReadFiducialImpulseSettings(HSQ,0)
#         WriteFiducialImpulseSettings(HSQ,0,4000,45000)#try full-time 15000,45000
#         time.sleep(.15)
#         for ii in range(4):
#             WritePulseHeights(HSQ,0,IndFETWave([jj*4+ii+1 for jj in range(35)],15000))#could improve by doing a few spread-out points at a time
#             time.sleep(6)
#             qpixnodata=readchan(1,LS1Q)['DATA']
#             #qpeakcoords=signal.find_peaks_cwt(np.clip(qpixnodata,max(qpixnodata)/4,max(qpixnodata)),np.arange(90,110))#threshold->clip
#             qpeakcoordspre=signal.find_peaks_cwt(np.clip(qpixnodata,max(qpixnodata)/4,max(qpixnodata))-max(qpixnodata)/4,np.arange(90,110));
#             qpeakcoords=[pt for pt in qpeakcoordspre if qpixnodata[pt]>max(qpixnodata)/4]
#             #this is 475ps/(2.5ps/pix)=190pix expected for S1 scope at max sampling; S2 scope needs different
#             if len(qpeakcoords) != 36:
#                 print('\nWrong number of peaks detected!\n',len(qpeakcoords))
#                 #time.sleep(.15)
#                 #HClose(HSQ)
#                 #time.sleep(.15)
#                 #LClose(LS1Q)
#                 #time.sleep(.15)
#                 return
#             #qdatalist.append([qpeakcoords[1:],qpixnodata[qpeakcoords[1:]]])
#             #try post-pulse trigger instead
#             qdatalist.append([qpeakcoords[:-1],qpixnodata[qpeakcoords[:-1]]])
#         #time.sleep(.15)
#         #HClose(HSQ)
#         #time.sleep(.15)
#         #LClose(LS1Q)
#         #time.sleep(.15)
#         qdatalist2=np.array(qdatalist).transpose(2,0,1).reshape(140,2)
#         qTimeErrInHpix=np.array([qdatalist2[ii,0]-np.mean(qdatalist2[:14,0])+100*(6.5-ii) for ii in range(140)])*2.5/1.8
#         qTimeErrInHpixBoardAvg=np.array([np.mean(qTimeErrInHpix[14*ii:14*ii+14]) for ii in range(10)])
#         epl(qTimeErrInHpix)
#         epl(qTimeErrInHpixBoardAvg)
#         return np.array([qTimeErrInHpix,qTimeErrInHpixBoardAvg])
# =============================================================================
# =============================================================================
# #FIX THIS!!!
#     def ___MiscCalCorrection(HSQ,MiscCalOldQ,TimeErrInHpixBoardAvgQ):
#         MiscCalNewQ=MiscCalOldQ[:]
#         for ii in range(10):
#             MiscCalNewQ[1+2*ii]=MiscCalNewQ[1+2*ii]-int(round(TimeErrInHpixBoardAvgQ[ii]))
#         WriteMiscellaneousCalibrations(HSQ,0,MiscCalNewQ)
#         time.sleep(.15)
#         return MiscCalNewQ
# =============================================================================
# =============================================================================
# #FIX THIS!!!        
#     def ___WaveTimeCalCorrection(HSQ,WaveTimeCalOldQ,TimeErrInHpixQ):
#         WaveTimeCalNewQ=WaveTimeCalOldQ[:]
#         for ii in range(140):
#             WaveTimeCalNewQ[ii]=WaveTimeCalNewQ[ii]-int(round(TimeErrInHpixQ[ii]))
#         WriteWaveTimeCalibrations(HSQ,0,WaveTimeCalNewQ)
#         time.sleep(.15)
#         return WaveTimeCalNewQ
# =============================================================================
# =============================================================================
# #FIX THIS!!!        
#     def ___ScanAndShift(HSQ,LS1Q):
#         #ideally: auto handling of sockets
#         #need error handling/better stability of code first
#         #HSQ=HOpen()
#         #time.sleep(.15)
#         #LS1Q=LOpen()
#         #time.sleep(.15)
#         #way to check scope settings? currently did this with the following:
#         #5 mV/div, 10sweeps, 5ns/div, -33ns delay, 13ns deskew, -13.8mV offset
#         ##### use on 20210209: on 5n/div on LeCroy1, trigger off another channel, set fiducial at 9th div
#         #ideally: read YFE settings, turn down, turn back up after done
#         #set fiducial and everything like that
#         scanresults=[]
#         #we don't care about historical values; we just want things fixed
#         #so don't read in/pass in parameters; just get them straight from HSQ
#         PulseHeightQ=ReadPulseHeights(HSQ,0)
#         MiscCalQ=ReadMiscellaneousCalibrations(HSQ,0)
#         time.sleep(.15)
#         WaveTimeCalQ=ReadWaveTimeCalibrations(HSQ,0)
#         time.sleep(.15)
#         scanresults.append(VeryFastFETSurvey(HSQ,LS1Q))
#         #test if need correction
#         if any(abs(elem)>2.5 for elem in scanresults[-1][1]):
#             print('Adjusting MiscCal\n')
#             MiscCalQ=MiscCalCorrection(HSQ,MiscCalQ,scanresults[-1][1])
#             time.sleep(.15)
#             scanresults.append(VeryFastFETSurvey(HSQ,LS1Q))
#         if any(abs(elem)>5.5 for elem in scanresults[-1][0]):
#             #this is factor of 2 away from "bad" ("20ps"=11.1 Hpix of error)
#             print('Adjusting WaveTimeCal\n')
#             WaveTimeCalQ=WaveTimeCalCorrection(HSQ,WaveTimeCalQ,scanresults[-1][0])
#             time.sleep(.15)
#             scanresults.append(VeryFastFETSurvey(HSQ,LS1Q))
#         if any(abs(elem)>2.5 for elem in scanresults[-1][1]) or any(abs(elem)>5.5 for elem in scanresults[-1][0]):
#             print('Consider running a second iteration')
#             #ideally: re-cal with for loops and iterate until corrected
#         #time.sleep(.15)
#         #HClose(HSQ)
#         #time.sleep(.15)
#         #LClose(LS1Q)
#         #time.sleep(.15)
#         WritePulseHeights(HSQ,0,PulseHeightQ)
#         return
# =============================================================================
# =============================================================================
# #FIX THIS!!!        
#     def ___FETsurveyfull():
#         HSQ=HOpen()
#         time.sleep(.15)
#         LS1Q=LXOpen('1')
#         time.sleep(.15)
#         qdatalist=[]
#         for ii in range(140):
#             qptdatalist=[]
#             for jj in range(5):
#                 WritePulseHeights(HSQ,0,IndFETWave([ii+1],int((jj+1)*65535/5)))#could improve by doing a few spread-out points at a time
#                 time.sleep(24)
#                 qpixnodata=readchan(1,LS1Q)['DATA'][2400:]
#                 qpixnomax=max(qpixnodata)
#                 qpixnomaxindex=np.mean([i for i,j in enumerate(qpixnodata) if j == qpixnomax])##could improve by changing to abbreviated centroid around peak, avoiding tail-end bump
#                 qptdatalist.append([2400+qpixnomaxindex,qpixnomax])
#             qdatalist.append(qptdatalist)
#         efc.pickledump2(qdatalist,psfilepath()+'fullFETsurvey20181106.p')
#         time.sleep(.15)
#         HClose(HSQ)
#         time.sleep(.15)
#         LXClose(LS1Q)
#         time.sleep(.15)
#         return qdatalist
# =============================================================================
# =============================================================================
# #FIX THIS!!!
#     def HParamReset():##need to fix pickle... !!!
#         HSQ=HOpen()
#         time.sleep(.15)
#         [qrphQ,qrwacQ,qrmcQ,qrwtcQ,qrwtQ]=pickle.load(open(psfilepath()+'HighlandParameterSnapshot20181116.p','rb'))#1108 original
#         WriteFiducialImpulseSettings(HSQ,0,0,0) #turn off fiducial; was at (HSQ,0,65535,0)
#         time.sleep(.15)
#         WritePulseHeights(HSQ,0,qrphQ)
#         #time.sleep(.15)
#         #WriteWaveAmplitudeCalibrations(HSQ,0,qrwacQ)
#         #time.sleep(.15)
#         #WriteMiscellaneousCalibrations(HSQ,0,qrmcQ)
#         #time.sleep(.15)
#         #WriteWaveTimeCalibrations(HSQ,0,qrwtcQ)
#         #time.sleep(.15)
#         #WriteWalkTable(HSQ,0,qrwtQ)
#         #time.sleep(.15)
#         time.sleep(.15)
#         HClose(HSQ)
#         time.sleep(.15) 
#         return
# =============================================================================
# =============================================================================
# #FIX THIS!!!
#     def findfid(TraceInQ): 
#         #make it fast by narrowing area where we know peak should be
#         TQ=TraceInQ[:]; maxTQ=max(TQ); minTQ=min(TQ);
#         TQP=signal.find_peaks_cwt(np.clip(TQ,(maxTQ-.8*(maxTQ-minTQ)),maxTQ),np.arange(5,15));
#         if (TQP[-1]-TQP[-2]) > 1000:
#             return TQP[-1]
#         else:
#             print('check your answer...')
#             return TQP[-1]
# =============================================================================
    
    def HParamSnapshot():
        """Once revived again, function will serve to save and preserve the various settings of the Highland"""
        pass
    


    

class LOSC:
    """Class containing all the necessary functions for running the LeCroy oscilloscopes
    Because we have several such scopes, instantiation of a certain device is required"""
    def __init__(self, LStrQ):
        """Initialize a LeCroy oscilloscope for use; possible choices are:
            LStrQ='A' #for the YFE, 1in1w, and SHG_opt diodes
            LStrQ='B' #for the four 2in1w diodes
            LStrQ='1' #for the instruments scientists' oscilloscope
            LStrQ='2' #for the four 2in2w diodes"""
        if   str(LStrQ).lower() == 'a':
            self._hostIP = '172.21.46.100'#'scope-ics-meclas-lecroy-a'
            self._name = 'LeCroyA'
        elif str(LStrQ).lower() == 'b':
            self._hostIP = '172.21.46.120'#'scope-ics-meclas-lecroy-b'#
            self._name = 'LeCroyB'
        elif str(LStrQ).lower() == '1':
            self._hostIP = '172.21.46.60'#'scope-ics-mectc1-1'
            self._name = 'LeCroy1'
        elif str(LStrQ).lower() == '2':
            self._hostIP = '172.21.46.128'#'scope-ics-meclas-lecroy-02'
            self._name = 'LeCroy2'
        else:
            print('Invalid scope name! Choose 1, 2, A, or B!!')
            return False
        self._LSock = None
        self._port = 1861

    def Open(self):
        """Takes care of opening the socket to the specified LeCroy; if called explicitly like this, it MUST be followed by a Close() statement or else you'll block the socket and need to locally disable/enable its networking card or power cycle the unit"""
        if not self._LSock:
            try:
                self._LSock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                self._LSock.settimeout(1.0)
                self._LSock.connect((self._hostIP,self._port))
            except:
                print(self._name +' NOT CONNECTED!')
        else:
            print('Socket may already be open!')

    def Close(self):
        """Takes care of closing the socket to the specified LeCroy"""
        try:
            self._LSock.close()
            self._LSock=None
        except:
            print('Unable to close socket -- it may already be closed')
        return

    def _send_and_reply(self,msg,SendOnly=False):
        """Generic utility for sending a poll to the specified LeCroy's internal processor and receiving its reply"""
        try:
            x = bytearray()
            msg=bytearray(msg,'utf8') ##NEW FIX ATTEMPT
            x.append(0x81)   # Data with EOI terminator
            x.append(1)      # Header v1
            x.append(0)      # Sequence Number
            x.append(0)      # Spare
            l = len(msg) + 1
            x.append((l >> 24) & 0xff)  # MSB!
            x.append((l >> 16) & 0xff)
            x.append((l >> 8) & 0xff)
            x.append((l >> 0) & 0xff)
            x.extend(msg)
            ##x.append('\n')
            x.extend(bytearray('\n','utf8'))#WAS APPEND
            self._LSock.sendall(x)
            if not SendOnly:
                data = ""
                done = False
                while not done:
                    hdr = self._LSock.recv(8) ##a bytes object
                    hdr = hdr.decode('latin1')##try sg...
                    done = (ord(hdr[0]) & 1) == 1  ##had an ord(hdr[0])
                    l = struct.unpack(">i", bytes(hdr[4:8],encoding='latin1'))[0]##ADDED bytes(...)
                    while (l != 0):
                        d = self._LSock.recv(l)##########################################
                        d = d.decode('latin1')##try sg....
                        data = data + d#.decode('utf-8')
                        l -= len(d)
                return data
        except:
            print('Send and reply failed!')

    @staticmethod
    def _LFields():
        """List of fields needed for interpreting read-out data of LeCroy fields"""
        fields = [
        [0,"DESCRIPTOR_NAME","string"],[16,"TEMPLATE_NAME","string"],[32,"COMM_TYPE","enum", {0:"byte",1: "word",}], 
        [34,"COMM_ORDER","enum", {0:"HIFIRST",1:"LOFIRST",}],[36,"WAVE_DESCRIPTOR","long"],[40,"USER_TEXT","long"], 
        [44,"RES_DESC1","long"],[48,"TRIGTIME_ARRAY","long"],[52,"RIS_TIME_ARRAY","long"],[56,"RES_ARRAY1","long"], 
        [60,"WAVE_ARRAY_1","long"],[64,"WAVE_ARRAY_2","long"],[68,"RES_ARRAY2","long"],[72,"RES_ARRAY3","long"], 
        [76,"INSTRUMENT_NAME","string"],[92,"INSTRUMENT_NUMBER","long"],[96,"TRACE_LABEL","string"], 
        [112,"RESERVED1","word"],[114,"RESERVED2","word"],[116,"WAVE_ARRAY_COUNT","long"],[120,"PNTS_PER_SCREEN","long"], 
        [124,"FIRST_VALID_PNT","long"],[128,"LAST_VALID_PNT","long"],[132,"FIRST_POINT","long"],
        [136,"SPARSING_FACTOR","long"],[140,"SEGMENT_INDEX","long"],[144,"SUBARRAY_COUNT","long"], 
        [148,"SWEEPS_PER_ACQ","long"],[152,"POINTS_PER_PAIR","word"],[154,"PAIR_OFFSET","word"], 
        [156,"VERTICAL_GAIN","float"],[160,"VERTICAL_OFFSET","float"],[164,"MAX_VALUE","float"], 
        [168,"MIN_VALUE","float"],[172,"NOMINAL_BITS","word"],[174,"NOM_SUBARRAY_COUNT","word"], 
        [176,"HORIZ_INTERVAL","float"],[180,"HORIZ_OFFSET","double"],[188,"PIXEL_OFFSET","double"], 
        [196,"VERTUNIT","unit_definition"],[244,"HORUNIT","unit_definition"],[292,"HORIZ_UNCERTAINTY","float"], 
        [296,"TRIGGER_TIME","time_stamp"],[312,"ACQ_DURATION","float"],
        [316,"RECORD_TYPE","enum",{0:"single_sweep",1:"interleaved",2:"histogram",3:"graph",4:"filter_coefficient", 
            5:"complex",6:"extrema",7:"sequence_obsolete",8:"centered_RIS",9:"peak_detect",}], 
        [318,"PROCESSING_DONE","enum",{0:"no_processing",1:"fir_filter",2:"interpolated",3:"sparsed", 
            4:"autoscaled",5:"no_result",6:"rolling",7:"cumulative",}], 
        [320,"RESERVED5","word"],[322,"RIS_SWEEPS","word"], 
        [324,"TIMEBASE","enum",{0:"1_ps/div",1:"2_ps/div",2:"5_ps/div",3:"10_ps/div",4:"20_ps/div",5:"50_ps/div", 
            6:"100_ps/div",7:"200_ps/div",8:"500_ps/div",9:"1_ns/div",10:"2_ns/div",11:"5_ns/div",12:"10_ns/div", 
            13:"20_ns/div",14:"50_ns/div",15:"100_ns/div",16:"200_ns/div",17:"500_ns/div",18:"1_us/div",19:"2_us/div", 
            20:"5_us/div",21:"10_us/div",22:"20_us/div",23:"50_us/div",24:"100_us/div",25:"200_us/div",26:"500_us/div", 
            27:"1_ms/div",28:"2_ms/div",29:"5_ms/div",30:"10_ms/div",31:"20_ms/div",32:"50_ms/div",33:"100_ms/div", 
            34:"200_ms/div",35:"500_ms/div",36:"1_s/div",37:"2_s/div",38:"5_s/div",39:"10_s/div",40:"20_s/div", 
            41:"50_s/div",42:"100_s/div",43:"200_s/div",44:"500_s/div",45:"1_ks/div",46:"2_ks/div",47:"5_ks/div", 
            100: "EXTERNAL",}], 
        [326,"VERT_COUPLING","enum",{0:"DC_50_Ohms",1:"ground",2:"DC_1MOhm",3:"ground",4:"AC_1MOhm",}], 
        [328,"PROBE_ATT","float"], 
        [332,"FIXED_VERT_GAIN","enum",{0:"1_uV/div",1:"2_uV/div",2:"5_uV/div",3:"10_uV/div",4:"20_uV/div", 
            5:"50_uV/div",6:"100_uV/div",7:"200_uV/div",8:"500_uV/div",9:"1_mV/div",10:"2_mV/div",11:"5_mV/div", 
            12:"10_mV/div",13:"20_mV/div",14:"50_mV/div",15:"100_mV/div",16:"200_mV/div",17:"500_mV/div", 
            18:"1_V/div",19:"2_V/div",20:"5_V/div",21:"10_V/div",22:"20_V/div",23:"50_V/div",24:"100_V/div", 
            25:"200_V/div",26:"500_V/div",27:"1_kV/div",}], 
        [334,"BANDWIDTH_LIMIT","enum",{0:"off",1:"on",}],[336,"VERTICAL_VERNIER","float"],[340,"ACQ_VERT_OFFSET","float"], 
        [344,"WAVE_SOURCE","enum",{0:"CHANNEL_1",1:"CHANNEL_2",2:"CHANNEL_3",3:"CHANNEL_4",9:"UNKNOWN",}],]
        return fields
        
    @classmethod
    def _parsewf(cls, data, verbose=False):
        """Internal function needed for parsing read-out data of LeCroy fields"""
        fields = cls._LFields()
        x = data.find(",#9")
        l = int(data[x+3:x+12])##
        data = data[x+12:x+12+l]
        d = {}
        for f in fields:
            if f[2] == "string" or f[2] == "unit_definition" or f[2] == "text":
                d[f[1]] = data[f[0]:f[0]+16].rstrip('\0')
                if (verbose): print("%30s    %s" % (f[1], d[f[1]]))
            elif f[2] == "enum":
                d[f[1]] = f[3][struct.unpack("<h", bytes(data[f[0]:f[0]+2],encoding='latin1'))[0]]##bytes(...,encoding='latin1')
                if (verbose): print("%30s    %s" % (f[1], d[f[1]]))
            elif f[2] == "word":
                d[f[1]] = struct.unpack("<h", bytes(data[f[0]:f[0]+2],encoding='latin1'))[0]##bytes...
                if (verbose): print("%30s    %s" % (f[1], d[f[1]]))
            elif f[2] == "long":
                d[f[1]] = struct.unpack("<i", bytes(data[f[0]:f[0]+4],encoding='latin1'))[0]##bytes...
                if (verbose): print("%30s    %i" % (f[1], d[f[1]]))
            elif f[2] == "float":
                d[f[1]] = struct.unpack("<f", bytes(data[f[0]:f[0]+4],encoding='latin1'))[0]##bytes...
                if (verbose): print("%30s    %g" % (f[1], d[f[1]]))
            elif f[2] == "double":
                d[f[1]] = struct.unpack("<d", bytes(data[f[0]:f[0]+8],encoding='latin1'))[0]##bytes...
                if (verbose): print("%30s    %g" % (f[1], d[f[1]]))
            elif f[2] == "time_stamp":
                d[f[1]] = "{}:{}:{} {}/{}/{}".format(data[f[0]+9],
                                                       data[f[0]+8],
                                                       struct.unpack("<d", bytes(data[f[0]:f[0]+8],encoding='latin1'))[0],##bytes...
                                                       data[f[0]+11],
                                                       data[f[0]+10],
                                                       struct.unpack("<h", bytes(data[f[0]+12:f[0]+14],encoding='latin1'))[0])
                if (verbose): print("%30s    %s" % (f[1], d[f[1]]))
            else:
                if (verbose): print("***** %24s    %s" % (f[1], f[2]))
        if struct.unpack("<h", bytes(data[32:34],encoding='latin1'))[0] == 0:##bytes...
            d['RAW'] = np.frombuffer(bytes(data[346:],encoding='latin1'), dtype=np.int8)
        else:
            d['RAW'] = np.frombuffer(bytes(data[346:],encoding='latin1'), dtype=np.int16)
        d['DATA'] = d['VERTICAL_GAIN'] * d['RAW'] - d['VERTICAL_OFFSET']
        return d
        
    def _FunctionWrapper(self,FuncQ,kwargs={}):
        """A function wrapper that allows one to call each LeCroy command directly without having to worry about opening/closing sockets;
        if issuing one-off commands that don't require high-frequency consistent execution, this is sufficient"""
        try:
            self.Open();time.sleep(0.15);
            LData=FuncQ(**kwargs);time.sleep(0.15);
            self.Close();time.sleep(0.15);
        except:
            self.Close();time.sleep(0.15);
            LData=False
        return LData

    def _waitrch(self,ChannelNo,verbose=False):
        """
        Bare function for reading voltage data from the specified channel after reading the internal state change register;
        socket must be explicitly opened/closed
        """
        while True:
            ready = (int(self._send_and_reply("INR?").split()[1]) & 1) == 1
            if ready:
                rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL".format(str(ChannelNo)))
                fullaq = self._parsewf(rawdataq, verbose)
                return fullaq
    def waitrch(self,ChannelNo,verbose=False):
        """
        Wrapped function for reading voltage data from the specified channel after reading the internal state change register;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._waitrch,{'ChannelNo':ChannelNo,'verbose':verbose});
        return LData
        
    def _waitrchxy(self,ChannelNo,verbose=False):
        """
        Bare function for reading voltage data and corresponding time data from the specified channel after reading the internal state change register;
        socket must be explicitly opened/closed
        """
        while True:
            ready = (int(self._send_and_reply("INR?").split()[1]) & 1) == 1
            if ready:
                rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL".format(str(ChannelNo)))
                fullaq = self._parsewf(rawdataq, verbose)
                yvals=fullaq['DATA'];xvals=[fullaq['HORIZ_OFFSET'] + ii*fullaq['HORIZ_INTERVAL'] for ii in range(len(fullaq['DATA']))];
                return [xvals,yvals]
    def waitrchxy(self,ChannelNo,verbose=False):
        """
        Wrapped function for reading voltage data and corresponding time data from the specified channel after reading the internal state change register;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._waitrchxy,{'ChannelNo':ChannelNo,'verbose':verbose});
        return LData
                
    def _rch(self,OChan):
        """
        Bare function for reading voltage data from the specified channel;
        socket must be explicitly opened/closed
        """
        rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL".format(str(OChan)))
        return self._parsewf(rawdataq, False)['DATA']
    def rch(self,OChan):
        """
        Wrapped function for reading voltage data from the specified channel;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._rch,{'OChan':OChan});
        return LData

    def _rchxy(self,OChan):
        """
        Bare function for reading voltage data and corresponding time data from the specified channel;
        socket must be explicitly opened/closed
        """
        rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL".format(str(OChan)))
        fullaq=self._parsewf(rawdataq, False);
        yvals=fullaq['DATA'];xvals=[fullaq['HORIZ_OFFSET'] + ii*fullaq['HORIZ_INTERVAL'] for ii in range(len(fullaq['DATA']))];
        return [xvals,yvals]
    def rchxy(self,OChan):
        """
        Wrapped function for reading voltage data and corresponding time data from the specified channel;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._rchxy,{'OChan':OChan});
        return LData
        
    def _rchall(self):
        """
        Bare function for reading voltage data from all channels of specified scope;
        socket must be explicitly opened/closed
        """
        rchans=[]
        for OChan in range(1,5):
            rchans.append(self._rch(OChan))
        return rchans
    def rchall(self):
        """
        Wrapped function for reading voltage data from all channels of specified scope;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._rchall);
        return LData
                
    def _rchallxy(self):
        """
        Bare function for reading voltage data and corresponding time data from all channels of specified scope;
        socket must be explicitly opened/closed
        """
        rchans=[]
        for OChan in range(1,5):
            rchans.append(self._rchxy(OChan))
        return rchans
    def rchallxy(self):
        """
        Wrapped function for reading voltage data and corresponding time data from all channels of specified scope;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._rchallxy);
        return LData

    def _sch(self,OChan,FileName):
        """
        Bare function for saving voltage data from the specified channel;
        socket must be explicitly opened/closed
        """
        rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL".format(OChan))
        parseddataq=self._parsewf(rawdataq, False)
        with open(LPL.psfilepath()+'data/'+str(FileName)+'.csv','w',newline='') as f:
            writer=csv.writer(f, delimiter='\n')
            writer.writerow(parseddataq['DATA'])
        with open(LPL.psfilepath()+'data/'+str(FileName)+'-h.csv','w',newline='') as f:
            writer=csv.DictWriter(f, parseddataq.keys())
            writer.writeheader()
            writer.writerow(parseddataq)
        return 
    def sch(self,OChan,FileName):
        """
        Wrapped function for saving voltage data from the specified channel;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._sch,{'OChan':OChan,'FileName':FileName});
        return LData
                
    def _schall(self,FileName):
        """
        Bare function for saving voltage data from all channels of specified scope;
        socket must be explicitly opened/closed
        """
        for OChan in range(1,5):
            self._sch(OChan,FileName+'_ch'+str(OChan))
        return
    def schall(self,FileName):
        """
        Wrapped function for saving voltage data from all channels of specified scope;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._schall,{'FileName':FileName});
        return LData

    def _pch(self,OChan):
        """
        Bare function for plotting voltage data from the specified channel;
        socket must be explicitly opened/closed
        """
        pdata=self._rch(OChan)
        efc.epl(pdata)
        return 
    def pch(self,OChan):
        """
        Wrapped function for plotting voltage data from the specified channel;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._pch,{'OChan':OChan});
        return LData
        
    def _pchxy(self,OChan):
        """
        Bare function for plotting voltage data and corresponding time data from the specified channel;
        socket must be explicitly opened/closed
        """
        pdata=self._rchxy(OChan)
        efc.eplxy(*pdata)
        return 
    def pchxy(self,OChan):
        """
        Wrapped function for plotting voltage data and corresponding time data from the specified channel;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._pchxy,{'OChan':OChan});
        return LData
        
    def _pchall(self):
        """
        Bare function for plotting voltage data from all channels of specified scope;
        socket must be explicitly opened/closed
        """
        pdata=self._rchall()
        efc.epll(pdata)
        return
    def pchall(self):
        """
        Wrapped function for plotting voltage data from all channels of specified scope;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._pchall);
        return LData
        
    def _pchallxy(self):
        """
        Bare function for plotting voltage data and corresponding time data from all channels of specified scope;
        socket must be explicitly opened/closed
        """
        pdata=self._rchallxy()
        efc.epllxy(pdata)
        return
    def pchallxy(self):
        """
        Wrapped function for plotting voltage data and corresponding time data from all channels of specified scope;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._pchallxy);
        return LData
        
    def _sumch(self,OChan):
        """
        Bare function for summing all voltage values from the specified channel;
        socket must be explicitly opened/closed
        """
        templistQ=self._rch(OChan)
        return np.sum(templistQ)
    def sumch(self,OChan):
        """
        Wrapped function for summing all voltage values from the specified channel;
        socket is automatically opened/closed
        """
        LData=self._FunctionWrapper(self._sumch,{'OChan':OChan});
        return LData
        
    def save_scope_to_eLog(self,chan_to_eLog=2):
        """Saves specified scope channel (data + voltage vs time plot) to the eLog"""
        ExpName=LPL.get_curr_exp()
        fpQ=str('/reg/neh/operator/mecopr/experiments/'+ExpName+'/lecroy_xray/')
        chan_to_eLog = int(chan_to_eLog)
        if chan_to_eLog not in [1,2,3,4]:
            print('Channel must be 1, 2, 3, or 4! Using channel 2 as default.')
            chan_to_eLog=2;
        if not os.path.exists(fpQ[-13]):
            print('File path '+fpQ[-13]+' does not exist! Trying to create it...')
            try:
                os.makedirs(fpQ[-13]);print('Folder created successfully!');
                os.chmod(fpQ[-13],stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO);
            except:
                print('Failed to create '+fpQ[-13]+'!')
        if not os.path.exists(fpQ):
            print('File path '+fpQ+' does not exist! Trying to create it...')
            try:
                os.makedirs(fpQ);print('Folder created successfully!');
                os.chmod(fpQ,stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO);
            except:
                print('Failed to create '+fpQ+'!')
        try:
            mecel = elog.ELog({'experiment':ExpName},user='mecopr',pw=pickle.load(open(LPL.psfilepath()+'elogauth.p','rb')))
            #print('Connected to eLog of current MEC experiment, which is: '+ExpName)
        except:
            print('Failed to connect to eLog!')
        try:
            chdata=self.rchallxy();
            timestamp=datetime.now().strftime('%Y%m%d_%H%M%S')
        except:
            print('Failed to read out data!')
        for ii in range(4):
            np.savetxt(fpQ+'lecroy1_ch'+str(ii+1)+'_'+timestamp+'.dat',tuple(chdata[ii]))
        efc.eplxysav(chdata[chan_to_eLog-1][0],chdata[chan_to_eLog-1][1],fpQ+'lecroy1_ch'+str(chan_to_eLog)+'_'+timestamp+'.png',abs_path=True)
        fullmsg=str('Scope trace data for all 4 channels saved to '+fpQ+' with time stamp '+timestamp+'. Attached are the data and plot files for channel '+str(chan_to_eLog)+'.')
        #eplxy(chdata[chan_to_eLog-1][0],chdata[chan_to_eLog-1][1])
        try:
            mecel.post(fullmsg,attachments=[fpQ+'lecroy1_ch'+str(chan_to_eLog)+'_'+timestamp+'.dat', fpQ+'lecroy1_ch'+str(chan_to_eLog)+'_'+timestamp+'.png'], tags=['scope_trace'])
            print('Auto-saved to eLog.') 
        except:
            print('Failed to auto-save to eLog!')
        
    def _ctrl(self,msg,SendOnly=True):
        """
        Once completed, this function will be used for sending specialized commands to the specified LeCroy; this can be found in the manual
        other old scribbles:
        SetScopeParameters(kwargs) --> scan for different inputs like those listed below; 
           could also specify use case and let it set or RCL appropriately
        (msg='TIME_DIV?',SendOnly=False)
        (msg='*RCL 3',SendOnly=False)
        (msg='TDIV 100E-9',SendOnly=True)
        (msg='C1:VDIV 500E-3',SendOnly=True)
        (msg=r\"""vbs 'app.acquisition.triggermode = "single" ' \""",SendOnly=True)
        (msg=r\"""vbs 'app.acquisition.c1.deskew = 0 ' \""",SendOnly=True)
        """
        self._send_and_reply(msg,SendOnly=SendOnly)





class EMeters:
    """Class containing readout functions for energy meters on all MEC laser systems"""
    def LPLInChamber(printDisplay=False):
        tvalW=GLOBAL.EGLPLWest.get();
        tvalE=GLOBAL.EGLPLEast.get();
        if printDisplay:
            print('Power meter: WEST: ' + str(tvalW) + ', EAST: ' + str(tvalE) + ', TOTAL: ' + str(tvalW+tvalE))
        return [tvalW, tvalE]

                  
    def EG():
        eab=GLOBAL.EGLPL2in2wAB.get()
        eef=GLOBAL.EGLPL2in2wEF.get()
        egh=GLOBAL.EGLPL2in2wGH.get()
        eij=GLOBAL.EGLPL2in2wIJ.get()
        EAB=round(GLOBAL.Ecoeff2in2wAB.get()*eab,4)
        EEF=round(GLOBAL.Ecoeff2in2wEF.get()*eef,4)
        EGH=round(GLOBAL.Ecoeff2in2wGH.get()*egh,4)
        EIJ=round(GLOBAL.Ecoeff2in2wIJ.get()*eij,4)
        guessarray=[[EAB,EEF,EGH,EIJ],round(EAB+EEF,4),round(EGH+EIJ,4),round(EAB+EEF+EGH+EIJ,4)]
        tempglobarr=[GLOBAL.EAB2w,GLOBAL.EEF2w,GLOBAL.EGH2w,GLOBAL.EIJ2w]
        for ii in range(4):
            try:
                tempglobarr[ii].put(guessarray[0][ii]*PFN.HeadENB()[ii])
            except:
                print('Failed to update notepad PV!')
        try:
            EABEF=GLOBAL.EGLPLWest.get()
            EGHIJ=GLOBAL.EGLPLEast.get()
        except:
            EABEF=-1;EGHIJ=-1
        realarray=[EABEF,EGHIJ,EABEF+EGHIJ]
        #print(realarray)
        return [guessarray,realarray]
        
    def EG1wYFE1in():
        eyfe=GLOBAL.EGLPLYFE.get()
        e1in=GLOBAL.EGLPL1in1w.get()
        EYFE,E1IN=GLOBAL.EcoeffYFE.get()*eyfe,GLOBAL.Ecoeff1in1wCD.get()*e1in#was 0.3578
        guessarray=[round(EYFE,4),round(E1IN,4)]
        tempglobarr=[GLOBAL.EYFE,GLOBAL.ECD1w]
        for ii in range(2):
            try:
                tempglobarr[ii].put(guessarray[ii])
            except:
                print('Failed to update notepad PV!')
        return guessarray

    def EG1w2in():
        eab=GLOBAL.EGLPL2in1wAB.get()
        eef=GLOBAL.EGLPL2in1wEF.get()
        egh=GLOBAL.EGLPL2in1wGH.get()
        eij=GLOBAL.EGLPL2in1wIJ.get()
        EAB=round(GLOBAL.Ecoeff2in1wAB.get()*eab,4)
        EEF=round(GLOBAL.Ecoeff2in1wEF.get()*eef,4)
        EGH=round(GLOBAL.Ecoeff2in1wGH.get()*egh,4)
        EIJ=round(GLOBAL.Ecoeff2in1wIJ.get()*eij,4)
        guessarray=[[EAB,EEF,EGH,EIJ],round(EAB+EEF,4),round(EGH+EIJ,4),round(EAB+EEF+EGH+EIJ,4)]
        tempglobarr=[GLOBAL.EAB1w,GLOBAL.EEF1w,GLOBAL.EGH1w,GLOBAL.EIJ1w]
        for ii in range(4):
            try:
                tempglobarr[ii].put(guessarray[0][ii]*PFN.HeadENB()[ii])
            except:
                print('Failed to update notepad PV!')
        return guessarray

                  
    def EGall(return_txt=False,chamber_meter_in=False):
        [en1wYFE, en1w1in] = EMeters.EG1wYFE1in()
        [en1wAB, en1wEF, en1wGH, en1wIJ] = EMeters.EG1w2in()[0]
        [en2wAB, en2wEF, en2wGH, en2wIJ] = EMeters.EG()[0][0]
        [enWEST, enEAST]=EMeters.EG()[1][:2]
        [cAB,cEF,cGH,cIJ]=PFN.HeadENB()
        ###
        wppvlist=[GLOBAL.HWPAB, GLOBAL.HWPEF, GLOBAL.HWPGH, GLOBAL.HWPIJ];
        headstr='';wpstr='';headlist=['AB','EF','GH','IJ'];
        for ii in range(4):
            wpstr=headstr+wpstr+headlist[ii]+': '+str(round(wppvlist[ii].get(),3))+', '
        wpenlist=[1e-12 + PFN.HeadENB()[ii]*np.cos((np.pi/180)*2*wppvlist[ii].get())**2 for ii in range(4)]
        ###
        strlist=[]
        strlist.append('Date: '+ datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
        print(strlist[-1]);strlist.append('\n');
        strlist.append('Time since last shot: '+GLOBAL.PFNSS.get())
        print(strlist[-1]);strlist.append('\n');
        Psns=GLOBAL.PSNS.get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        SSs=[GLOBAL.SSS.get()[2*ii:2*ii+2] for ii in range(len(Psns))]#pickle.load(open(psfpQ+'SSs.p','rb'))
        strlist.append('Current pulse target is: '+str(Psns)+' ns, '+str(SSs)+' % of max power.')
        print(strlist[-1]);strlist.append('\n');
        strlist.append('YFE: '+'{:5.1f}'.format(1000*round(en1wYFE,4))+'mJ, 1"@1w: '+'{:5.2f}'.format(en1w1in)+'J')
        print(strlist[-1]);strlist.append('\n');
        strlist.append('2"@1w: AB: '+'{:5.2f}'.format(round(cAB*en1wAB,4))+'J, EF: '+'{:5.2f}'.format(round(cEF*en1wEF,4))+'J, GH: '+'{:5.2f}'.format(round(cGH*en1wGH,4))+'J, IJ: '+'{:5.2f}'.format(round(cIJ*en1wIJ,4))+'J')
        print(strlist[-1]);strlist.append('\n');
        strlist.append('2"@2w: AB: '+'{:5.2f}'.format(round(cAB*en2wAB,4))+'J, EF: '+'{:5.2f}'.format(round(cEF*en2wEF,4))+'J, GH: '+'{:5.2f}'.format(round(cGH*en2wGH,4))+'J, IJ: '+'{:5.2f}'.format(round(cIJ*en2wIJ,4))+'J')
        print(strlist[-1]);strlist.append('\n');
        if np.sum(wpenlist) < 4:
            #strlist.append('2"@2w/HWP: AB: '+'{:5.2f}'.format(round(cAB*en2wAB/wpenlist[0],4))+'J, EF: '+'{:5.2f}'.format(round(cEF*en2wEF/wpenlist[1],4))+'J, GH: '+'{:5.2f}'.format(round(cGH*en2wGH/wpenlist[2],4))+'J, IJ: '+'{:5.2f}'.format(round(cIJ*en2wIJ/wpenlist[3],4))+'J')
            #print(strlist[-1]);strlist.append('\n');
            #strlist.append('Conv%/HWP: AB: '+'{:5.2f}'.format(round(100*cAB*en2wAB/en1wAB/wpenlist[0],4))+'%, EF: '+'{:5.2f}'.format(round(100*cEF*en2wEF/en1wEF/wpenlist[1],4))+'%, GH: '+'{:5.2f}'.format(round(100*cGH*en2wGH/en1wGH/wpenlist[2],4))+'%, IJ: '+'{:5.2f}'.format(round(100*cIJ*en2wIJ/en1wIJ/wpenlist[3],4))+'%')
            #print(strlist[-1]);strlist.append('\n');
            pass
        else:
            strlist.append('Conv%: AB: '+'{:5.2f}'.format(round(100*cAB*en2wAB/en1wAB,4))+'%, EF: '+'{:5.2f}'.format(round(100*cEF*en2wEF/en1wEF,4))+'%, GH: '+'{:5.2f}'.format(round(100*cGH*en2wGH/en1wGH,4))+'%, IJ: '+'{:5.2f}'.format(round(100*cIJ*en2wIJ/en1wIJ,4))+'%')
            print(strlist[-1]);strlist.append('\n');
        if chamber_meter_in:
            strlist.append('Measured energy: WEST: '+'{:5.2f}'.format(round(enWEST,4))+'J, EAST: '+'{:5.2f}'.format(round(enEAST,4))+'J')
            print(strlist[-1]);strlist.append('\n');
        strlist.append('Inferred energy: WEST: '+'{:5.2f}'.format(round(cAB*en2wAB+cEF*en2wEF,4))+'J, EAST: '+'{:5.2f}'.format(round(cGH*en2wGH+cIJ*en2wIJ,4))+'J, TOTAL: '+'{:5.2f}'.format(round(cAB*en2wAB+cEF*en2wEF+cGH*en2wGH+cIJ*en2wIJ,4)))
        print(strlist[-1]);strlist.append('\n');
        if return_txt:
            tot_msg='';
            for seg in strlist:
                tot_msg+=seg
            return tot_msg
        return [en1wYFE, en1w1in, cAB*en1wAB, cEF*en1wEF, cGH*en1wGH, cIJ*en1wIJ, cAB*en2wAB, cEF*en2wEF, cGH*en2wGH, cIJ*en2wIJ]

    def SPLEG():
        #later add in GLOBALS.EREGEN .ETOPAS .EMPA1 .EMPA2
        reen=GLOBAL.EGSPLRE.get() 
        toen=GLOBAL.EGSPLTO.get() 
        m1en=GLOBAL.EGSPLM1.get() 
        m2en=GLOBAL.EGSPLM2.get() 
        regenen=GLOBAL.EcoeffRE1*reen + GLOBAL.EcoeffRE0
        topasen=GLOBAL.EcoeffTO1*toen + GLOBAL.EcoeffTO0
        mpa1en=GLOBAL.EcoeffM11*m1en + GLOBAL.EcoeffM10 #1.76e5, -2.22e0
        mpa2en=GLOBAL.EcoeffM21*m2en + GLOBAL.EcoeffM20
        return np.round([regenen,topasen,mpa1en,mpa2en],2)

                  
    def gentec_refresh():
        pvhead='MEC:LAS:GENTEC:';pvtails=['DESCRIPTION','SET_WAVLEN','SET_SCALE','SET_TRIGMODE','SET_TRIGLVL','SET_ATTENUATOR'];#'SET_SCALE' was 'SCALE'
        pvids=['02:CH1:','02:CH2:','01:CH1:','01:CH2:','03:CH1:','03:CH2:','04:CH1:','04:CH2:'];
        abirvals=['AB IRsamp',1053,24,1,2,0]#24 was '1' with 'SCALE'
        efirvals=['EF IRsamp',1053,24,1,2,0]
        ghirvals=['GH IRsamp',1053,24,1,2,0]
        ijirvals=['IJ IRsamp',1053,24,1,2,0]

        ab2wvals=['AB 2wsamp',527,23,1,2,0]#23 was '300m' with 'SCALE'
        ef2wvals=['EF 2wsamp',527,23,1,2,0]
        gh2wvals=['GH 2wsamp',527,23,1,2,0]
        ij2wvals=['IJ 2wsamp',527,23,1,2,0]

        pvwest=['MEC:GENTEC:01:CH2:DESCRIPTION','MEC:GENTEC:01:CH2:SET_WAVLEN','MEC:GENTEC:01:CH2:SET_SCALE','MEC:GENTEC:01:CH2:SET_TRIGMODE','MEC:GENTEC:01:CH2:SET_TRIGLVL','MEC:GENTEC:01:CH2:SET_ATTENUATOR']
        westvals=['West ABEF',527,28,1,2,1]#28 was '100' with 'SCALE'

        pveast=['MEC:GENTEC:01:CH1:DESCRIPTION','MEC:GENTEC:01:CH1:SET_WAVLEN','MEC:GENTEC:01:CH1:SET_SCALE','MEC:GENTEC:01:CH1:SET_TRIGMODE','MEC:GENTEC:01:CH1:SET_TRIGLVL','MEC:GENTEC:01:CH1:SET_ATTENUATOR']
        eastvals=['East GHIJ, short pulse',527,28,1,2,1]#28 was '100' with 'SCALE'

        pvgroups = [[pvhead+pvid+pvtail for pvtail in pvtails] for pvid in pvids];pvgroups.append(pvwest);pvgroups.append(pveast);
        valgroups = [abirvals,efirvals,ghirvals,ijirvals,ab2wvals,ef2wvals,gh2wvals,ij2wvals,westvals,eastvals]
        for pvgroup,valgroup in zip(pvgroups,valgroups):
            for pvx,valx in zip(pvgroup,valgroup):
                temppv=EpicsSignal(pvx)
                temppv.put(valx)       

    def E_coeff_refresh():
        GLOBAL.notepadPVreset()
        #pvlist=['MEC:LAS:FLOAT:'+str(ii) for ii in range(31,41)];
        #inddesclist=['YFE','CD1w','AB1w','EF1w','GH1w','IJ1w','AB2w','EF2w','GH2w','IJ2w']
        #desclist=['E_coeff_'+inddesc for inddesc in inddesclist]
        #valulist=[.3578,0.5971,224.0,177.5,307.4*0.849,113.2,111.0*1.17,187.9*0.860,182.1*0.897,123.5*1.25]
        #for jj in range(len(pvlist)):
        #    temppv1=EpicsSignal(str(pvlist[jj]+'.DESC'));temppv2=EpicsSignal(pvlist[jj]);
        #    temppv1.put(desclist[jj]);temppv2.put(valulist[jj]);

                  
    def E_synth_refresh():
        pvlist=['MEC:LAS:FLOAT:'+str(ii).zfill(2) for ii in range(1,11)];
        inddesclist=['YFE','CD1w','AB1w','EF1w','GH1w','IJ1w','AB2w','EF2w','GH2w','IJ2w']
        desclist=['E_synth_'+inddesc for inddesc in inddesclist]
        #add in again for SPL meters once coefficient is stable
        #pvlist2=['MEC:LAS:FLOAT:'+str(ii).zfill(2) for ii in range(21,25)];
        #inddesclist2=['REGEN','TOPAS','MPA1','MPA2']
        #desclist2=['E_synth_'+inddesc2 for inddesc2 in inddesclist2]

        eyfe=GLOBAL.EGLPLYFE.get();
        e1in=GLOBAL.EGLPL1in1w.get();
        eab1w=GLOBAL.EGLPL2in1wAB.get();
        eef1w=GLOBAL.EGLPL2in1wEF.get();
        egh1w=GLOBAL.EGLPL2in1wGH.get();
        eij1w=GLOBAL.EGLPL2in1wIJ.get();
        eab2w=GLOBAL.EGLPL2in2wAB.get();
        eef2w=GLOBAL.EGLPL2in2wEF.get();
        egh2w=GLOBAL.EGLPL2in2wGH.get();
        eij2w=GLOBAL.EGLPL2in2wIJ.get();
        energyarr=np.array([eyfe,e1in,eab1w,eef1w,egh1w,eij1w,eab2w,eef2w,egh2w,eij2w])

        coefflist=[]
        for ii in range(31,41):
            temppv=EpicsSignal('MEC:LAS:FLOAT:'+str(ii));
            coefflist.append(temppv.get());

        valulist=energyarr*np.array(coefflist)
        for jj in range(len(pvlist)):
            temppv1=EpicsSignal(str(pvlist[jj]+'.DESC'));temppv2=EpicsSignal(pvlist[jj]);
            temppv1.put(desclist[jj]);temppv2.put(valulist[jj]);
#        for jj in range(len(pvlist2)):
#            temppv1=EpicsSignal(str(pvlist2[jj]+'.DESC'));temppv2=EpicsSignal(pvlist2[jj]);
#            temppv1.put(desclist2[jj]);temppv2.put(valulist2[jj]);
        return




                  
class MBC:
    def MBCmodecheck():
        return GLOBAL.MBCmode.get()

    def isMBCsafe():#re-write checks as individual functions
        status = True
        print('Checking MBC status...')
        if GLOBAL.MBCpwr.get() != 1:
            print('MBC is not on!!')
            status*=False
        if GLOBAL.MBCmode.get() != 0:
            print('MBC is not in AUTO mode!!')
            status*=False
        if GLOBAL.MBCsetpt.get() != 1:
            print('MBC is not in MIN mode!!')
            status*=False
        if not -7000<GLOBAL.MBCbias.get()<7000:
            print('MBC is out of range!')
            status*=False
        if GLOBAL.MBCfault.get() != 0:
            print('MBC fault detected!')
            status*=False
        if status:
            biaschk=[]
            print('Checking MBC bias level...',end='',flush=True) 
            for ii in range(3):
                biaschk.append(GLOBAL.MBCbias.get())
                time.sleep(1);print('..',end='',flush=True);time.sleep(1);print('..',end='',flush=True);
            print('*')
            if np.max(np.abs(np.diff(biaschk))) > 5:
                print('MBC bias level unstable!')
                return False
            else:
                return True
        else:
            return False

    def resetMBC():
        YFE.SetAll(False);
        #add KeyboardInterrupt?
        print('Begin resetting the MBC...')
        if GLOBAL.MBCpwr.get() != 1:
            print('Powering on MBC, starting scan...',end='',flush=True)
            GLOBAL.MBCpwr.put(1);time.sleep(1);print('.',end='',flush=True);GLOBAL.MBCmode.put(0);
            efc.dotsleep(8);
        if GLOBAL.MBCfault.get() != 0:
            print('Attempting to reset MBC fault...',end='',flush=True)
            GLOBAL.MBCfault.put(1)
            time.sleep(2);print('*');
        if GLOBAL.MBCmode.get() != 0:
            print('Setting MBC to AUTO mode, starting scan...',end='',flush=True)
            GLOBAL.MBCmode.put(0)
            efc.dotsleep(8);
        if GLOBAL.MBCsetpt.get() != 1:
            print('Setting MBC to MIN mode,starting scan...',end='',flush=True)
            GLOBAL.MBCsetpt.put(1)
            time.sleep(2);print('*');
        inibias=GLOBAL.MBCbias.get()
        if not -7000<inibias<7000:
            print('MBC is out of range! Aborting and power-cycling...');
            GLOBAL.MBCbias.put((np.round(time.time()*1000)%2)*9000*np.sign(inibias));time.sleep(1);
            GLOBAL.MBCpwr.put(2);time.sleep(2);
            MBC.resetMBC();return
        biaschk=[]
        print('Checking the initial MBC bias level...',end='',flush=True)
        for ii in range(3):
            biaschk.append(GLOBAL.MBCbias.get())
            time.sleep(1);print('..',end='',flush=True);time.sleep(1);print('..',end='',flush=True);
        print('*')
        waitloop=True;loopcnt=0;
        biaschklog=[]
        biaschklog.append(np.abs(np.diff(biaschk)))
        while waitloop:
            newchk=np.abs(np.diff(biaschk))
            biaschklog.append(newchk)
            if np.sum(np.abs(np.diff(biaschk))) > 3:
                print('MBC bias level unstable... '+str(biaschk),end='',flush=True)
                biaschk=[]
                for ii in range(3):
                    biaschk.append(GLOBAL.MBCbias.get())
                    time.sleep(1);print('..',end='',flush=True);time.sleep(1);print('..',end='',flush=True);
                print('')
                loopcnt+=1
                if (loopcnt >= 15) and (biaschklog[-1] > biaschklog[-1]):
                    print('MBC bias level stability fail. Aborting and power-cycling...')
                    GLOBAL.MBCbias.put((np.round(time.time()*1000)%2)*9000*np.sign(biaschk[-1]));time.sleep(1);
                    GLOBAL.MBCpwr.put(2);time.sleep(2);
                    MBC.resetMBC();return
            else:
                print('MBC bias level stabilized... '+str(biaschk))
                waitloop = False
        return





class YFE:
    """Class for organizing functions associated with the YLF Front End (YFE) laser system"""
    def OnCheck(display=True):
        YFEadd='MEC:LPL:LCO:0'
        YFEamp=['2','3','5','6','1','4']
        YFEsuf=[':SensedCurrent',':ActiveCurrent',':PowerSupply',':Temperature',':Emission_RBV',':Emission',':FaultState.RVAL',':ClearFault']
        statuslist=[]
        for ii in range(len(YFEamp)):
            temprbvemispv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[4])
            statuslist.append(temprbvemispv.get())#check emission
        if np.sum(statuslist)==0:
            YFEonbool=False
        if 0<np.sum(statuslist)<6:
            print('Warning: YFE seems to be partially on/off.')
            YFEonbool=False
        if np.sum(statuslist)==6:
            YFEonbool=True
        if display:
            print('Current status: '+str(statuslist))
        return YFEonbool

    @classmethod
    def On(cls):
        if YFE.OnCheck(dsplay=False):
            print('YFE emission already enabled.');return
        if GLOBAL.LPLPCpwr.get() != 1:
            GLOBAL.LPLPCpwr.put(1)
        if GLOBAL.LPLVACpwr.get() != 1:
            GLOBAL.LPLVACpwr.put(1)
        if GLOBAL.LPLPS1pwr.get() != 1:
            GLOBAL.LPLPS1pwr.put(1)
        YFEadd='MEC:LPL:LCO:0'
        YFEamp=['2','3','5','6','1','4']
        YFEsuf=[':SensedCurrent',':ActiveCurrent',':PowerSupply',':Temperature',':Emission_RBV',':Emission',':FaultState.RVAL',':ClearFault']
        print('Turning on YFE! Checking for faults...')
        faultpvlist=[EpicsSignal(YFEadd+amplabel+YFEsuf[6], write_pv=YFEadd+amplabel+YFEsuf[7]) for amplabel in YFEamp]
        emisspvlist=[EpicsSignal(YFEadd+amplabel+YFEsuf[4]) for amplabel in YFEamp]
        faultstatlist=[faultpv.get() for faultpv in faultpvlist]
        if any(faultstatlist):
            print('YFE fault state detected! Trying to reset...')
            for faultpv in faultpvlist:
                faultpv.put(1)
            time.sleep(3)
            faultstatlist=[faultpv.get() for faultpv in faultpvlist]
            if any(faultstatlist):
                print('YFE fault state still detected, turn-on failed.')
                return False
            else:
                print('Fault cleared!')
        if not MBC.IsSafe():
            cls.SetAll(False,displayQ=False)
            print('MBC not configured properly!')
            MBC.Reset()
            cls.SetAll(True,displayQ=False);
        else: #later add check to avoid over-energizing by reading power meter
            for ii in range(len(YFEamp)):
                tempsetcurrpv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[1])
                tempsetemispv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[5])
                tempsetcurrpv.put(0)#set current to 0
                tempsetemispv.put(1)#turn on emission
            print('Initializing eDrives...',end='',flush=True)
            efc.dotsleep(10);
            emissstatlist=[emisspv.get() for emisspv in emisspvlist]
            if all(emissstatlist):
                print('Ramping up currents...')
                cls.SetAll(True);
                print('YFE LASER ON')
            else:
                print('Turn on sequence failed. Check emission!')
                cls.Off();
                return False
        if GLOBAL.LPLPCpwr.get() != 1:
            print('Failed to turn on Pockels cell!')
        if GLOBAL.LPLVACpwr.get() != 1:
            print('Failed to turn on scroll pump!')
        if GLOBAL.LPLPS1pwr.get() != 1:
            print('Failed to turn on YFE PS1!')
        return True

    def Off():
        GLOBAL.LPLPCpwr.put(2)
        GLOBAL.LPLVACpwr.put(2)
        GLOBAL.LPLPS1pwr.put(2)
        YFEadd='MEC:LPL:LCO:0'
        YFEamp=['2','3','5','6','1','4']
        YFEsuf=[':SensedCurrent',':ActiveCurrent',':PowerSupply',':Temperature',':Emission_RBV',':Emission',':FaultState.RVAL']
        for ii in range(len(YFEamp)):
            tempsetcurrpv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[1])
            tempsetemispv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[5])
            tempsetcurrpv.put(0)#set current to 0
            tempsetemispv.put(0)#turn on emission
        GLOBAL.PFNmode.put(0)#makes sure glass chillers turn off when YFEoff... people have been leaving them on
        GLOBAL.MBCpwr.put(2)
        print('YFE LASER OFF')
        return
        
    def Get(display=True):
        YFEadd='MEC:LPL:LCO:0'
        YFEamp=['2','3','5','6','1','4']
        YFEsuf=[':SensedCurrent',':ActiveCurrent',':PowerSupply',':Temperature',':Emission_RBV',':Emission',':FaultState.RVAL']
        currreqQ=[]
        curractQ=[]
        for ii in range(len(YFEamp)):
            tempsetcurrpv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[1])
            tempactcurrpv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[0])
            currreqQ.append(round(tempsetcurrpv.get(),4))
            curractQ.append(round(tempactcurrpv.get(),4))
        if display:
            print('Requested currents: '+str(currreqQ))
            print('Actual currents:    '+str(curractQ))
        return currreqQ

    def Set(mmQ,currQ,display=True):
        YFEadd='MEC:LPL:LCO:0'
        YFEamp=['2','3','5','6','1','4']
        YFEsuf=[':SensedCurrent',':ActiveCurrent',':PowerSupply',':Temperature',':Emission_RBV',':Emission',':FaultState.RVAL']
        changelist=[]
        if mmQ==2:
            changelist=range(4)
            tempoldcurr02pv=EpicsSignal(YFEadd+YFEamp[changelist[0]]+YFEsuf[1])
            oldcurrQ=str(tempoldcurr02pv.get())
            if currQ>88:
                print('Too high!')
                return
        elif mmQ==6:
            changelist=[4]
            tempoldcurr06pv=EpicsSignal(YFEadd+YFEamp[changelist[0]]+YFEsuf[1])
            oldcurrQ=str(tempoldcurr06pv.get())
            if currQ>135:
                print('Too high!')
                return
        elif mmQ==10:
            changelist=[5]
            tempoldcurr10pv=EpicsSignal(YFEadd+YFEamp[changelist[0]]+YFEsuf[1])
            oldcurrQ=str(tempoldcurr10pv.get())
            if currQ>140:
                print('Too high!')
                return
        else:
            print('No such head!')
            return
        currmeapv=EpicsSignal(YFEadd+YFEamp[changelist[0]]+YFEsuf[0])
        currmeaQ=currmeapv.get()
        for amphead in changelist:
            tempnewcurrpv=EpicsSignal(YFEadd+YFEamp[amphead]+YFEsuf[1])
            if currQ-currmeaQ>20:
                nostep=int((currQ-currmeaQ)/20.0)
                for ii in range(nostep):
                    tempnewcurrpv.put(currQ-int(((1.0*nostep-ii)/(nostep+1.0))*(currQ-currmeaQ)))
                    time.sleep(1)
            tempnewcurrpv.put(currQ)
        if display:
            print(str(mmQ)+' mm changed from ' + oldcurrQ + ' to ' + str(currQ))
        return

    @classmethod
    def SetAll(cls,IOBool,displayQ=False):
        YFEadd='MEC:LPL:LCO:0'
        YFEamp=['2','3','5','6','1','4']
        YFEsuf=[':SensedCurrent',':ActiveCurrent',':PowerSupply',':Temperature',':Emission_RBV',':Emission',':FaultState.RVAL']
        YFElvl=[85,85,85,85,130,124];
        temppvlist=[]
        for ii in range(6):
            temppvlist.append(EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[1]))
        if IOBool:
            print('Ramping currents...',end='',flush=True)
            for jj in range(5):
                for eachpv,pvlvl in zip(temppvlist,YFElvl):
                    eachpv.put(pvlvl*(jj+1)/(5.0));
                time.sleep(1);print('..',end='',flush=True);
            print('*')
        else:
            print('Zeroing currents...',end='',flush=True);
            for eachpv in temppvlist:
                eachpv.put(0)
            time.sleep(1.5);print('*');
        cls.Get(display=displayQ);
        return

    def YFEtrace():#fix this
        try:
            LOSC('a').pch(1)
        except:
            print('Failed to display trace')
            return False





class PFN:
    def HeadENB():
        cAB=GLOBAL.PFNAEN.get()*GLOBAL.PFNBEN.get()
        cEF=GLOBAL.PFNEEN.get()*GLOBAL.PFNFEN.get()
        cGH=GLOBAL.PFNGEN.get()*GLOBAL.PFNHEN.get()
        cIJ=GLOBAL.PFNIEN.get()*GLOBAL.PFNJEN.get()
        return np.array([cAB,cEF,cGH,cIJ])
        
# =============================================================================
#     def PFNon(ArmStrQ):
#         if ArmStrQ.lower() == 'all':
#             ArmStrQ = 'ABCDEFGHIJ'
#         GLOBAL.PFNmode.put(0)
#         time.sleep(2)
#         if ('A' in ArmStrQ) or ('a' in ArmStrQ):
#             GLOBAL.PFNAEN.put(1)
#         if ('B' in ArmStrQ) or ('b' in ArmStrQ):
#             GLOBAL.PFNBEN.put(1)
#         if ('CD' in ArmStrQ) or ('cd' in ArmStrQ):
#             GLOBAL.PFNCDEN.put(1)
#         if ('E' in ArmStrQ) or ('e' in ArmStrQ):
#             GLOBAL.PFNEEN.put(1)
#         if ('F' in ArmStrQ) or ('f' in ArmStrQ):
#             GLOBAL.PFNFEN.put(1)
#         if ('G' in ArmStrQ) or ('g' in ArmStrQ):
#             GLOBAL.PFNGEN.put(1)
#         if ('H' in ArmStrQ) or ('h' in ArmStrQ):
#             GLOBAL.PFNHEN.put(1)
#         if ('I' in ArmStrQ) or ('i' in ArmStrQ):
#             GLOBAL.PFNIEN.put(1)
#         if ('J' in ArmStrQ) or ('j' in ArmStrQ):
#             GLOBAL.PFNJEN.put(1)
#         time.sleep(2);GLOBAL.PFNmode.put(1);time.sleep(3.5);GLOBAL.PFNmode.put(2);
#         return
# 
#     def PFNoff(ArmStrQ):
#         if ArmStrQ.lower() == 'all':
#             ArmStrQ = 'ABCDEFGHIJ'
#         GLOBAL.PFNmode.put(0)
#         time.sleep(2)
#         if ('A' in ArmStrQ) or ('a' in ArmStrQ):
#             GLOBAL.PFNAEN.put(0)
#         if ('B' in ArmStrQ) or ('b' in ArmStrQ):
#             GLOBAL.PFNBEN.put(0)
#         if ('CD' in ArmStrQ) or ('cd' in ArmStrQ):
#             GLOBAL.PFNCDEN.put(0)
#         if ('E' in ArmStrQ) or ('e' in ArmStrQ):
#             GLOBAL.PFNEEN.put(0)
#         if ('F' in ArmStrQ) or ('f' in ArmStrQ):
#             GLOBAL.PFNFEN.put(0)
#         if ('G' in ArmStrQ) or ('g' in ArmStrQ):
#             GLOBAL.PFNGEN.put(0)
#         if ('H' in ArmStrQ) or ('h' in ArmStrQ):
#             GLOBAL.PFNHEN.put(0)
#         if ('I' in ArmStrQ) or ('i' in ArmStrQ):
#             GLOBAL.PFNIEN.put(0)
#         if ('J' in ArmStrQ) or ('j' in ArmStrQ):
#             GLOBAL.PFNJEN.put(0)
#         time.sleep(2);GLOBAL.PFNmode.put(1);time.sleep(3.5);GLOBAL.PFNmode.put(2);
#         return
# =============================================================================

    def PFNonly(ArmStrQ):#delete and just use ARMonly?
        AllStrQ='ABEFGHIJ'
        AllStrq='abefghij'
        if ArmStrQ.lower() == 'all':
            ArmStrQ = AllStrQ
        GLOBAL.PFNmode.put(0)
        time.sleep(2)
        for ii in range(len(AllStrQ)):
            if (AllStrQ[ii] in ArmStrQ) or (AllStrq[ii] in ArmStrQ):
                temppv=EpicsSignal(str('MEC:PFN:CH'+str(ii+1)+':ENABLE'))
                temppv.put(1)
            else:
                temppv=EpicsSignal(str('MEC:PFN:CH'+str(ii+1)+':ENABLE'))
                temppv.put(0)
        time.sleep(2);GLOBAL.PFNmode.put(1);time.sleep(3.5);GLOBAL.PFNmode.put(2);
        return

    @classmethod
    def ARMonly(cls,ArmStrQ,set_T=1):
        HWP.On(ArmStrQ,set_T=set_T)
        cls.PFNonly(ArmStrQ)
        return
        
    
    
    
        
class HWP:
    def On(ArmStrQ,set_T=1):#fullON=1;fullOFF=0
        armlist=['AB','EF','GH','IJ'];
        if (set_T<0) or (set_T>1):
            set_T=1;print('Error: set_T must be between 0 and 1! Using set_T=1 instead...')
        if ArmStrQ.lower() == 'all':
            ArmStrQ = 'ABEFGHIJ'
        HWPpvlist=[GLOBAL.HWPAB,GLOBAL.HWPEF,GLOBAL.HWPGH,GLOBAL.HWPIJ];
        motor_angle=np.arccos(np.sqrt(set_T))*(180/np.pi)*0.5 
        print('Moving waveplates '+ArmStrQ+' to '+str(round(motor_angle,4))+'deg...');
        ABchk=('A' in ArmStrQ) or ('a' in ArmStrQ) or ('B' in ArmStrQ) or ('b' in ArmStrQ);
        EFchk=('E' in ArmStrQ) or ('e' in ArmStrQ) or ('F' in ArmStrQ) or ('f' in ArmStrQ);
        GHchk=('G' in ArmStrQ) or ('g' in ArmStrQ) or ('H' in ArmStrQ) or ('h' in ArmStrQ);
        IJchk=('I' in ArmStrQ) or ('i' in ArmStrQ) or ('J' in ArmStrQ) or ('j' in ArmStrQ);
        chklist=[ABchk,EFchk,GHchk,IJchk];
        HWPoldrbvlist=[GLOBAL.HWPAB.get(),GLOBAL.HWPEF.get(),GLOBAL.HWPGH.get(),GLOBAL.HWPIJ.get()];
        for ii in range(4):
            if chklist[ii]:
                HWPpvlist[ii].put(motor_angle)
        time.sleep(2)
        for ii in range(4):
            chklist[ii] = chklist[ii] and ((np.abs(HWPoldrbvlist[ii]-HWPpvlist[ii].get()) < 0.01) and (np.abs(motor_angle-HWPpvlist[ii].get()) > 0.01));
        retryii=0;
        while any(chklist):
            retryii+=1;
            for ii in range(4):
                if chklist[ii]:
                    HWPpvlist[ii].put(motor_angle);print('HWP issue detected on '+armlist[ii]+'. Re-trying...');
                time.sleep(2)
                chklist[ii] = chklist[ii] and ((np.abs(HWPoldrbvlist[ii]-HWPpvlist[ii].get()) < 0.01) and (np.abs(motor_angle-HWPpvlist[ii].get()) > 0.01)) and (retryii<3);
        if any(chklist):
            for ii in range(4):
                if chklist[ii]:
                    print('Re-try on '+armlist[ii]+' failed!')
        return
        
    def ClearStart():
        resetPVs=[GLOBAL.HWPABclr, GLOBAL.HWPEFclr]#no equivalent for GH and IJ yet...
        try:
            for eapv in resetPVs:
                eapv.put(1)
        except:
            print('Failed!')
        
    def HWP_opt(armsQ='ABEFGHIJ'):#check for trace height;#All shutters must start in the open state... 
        print('Running this routine requires ALL TTL shutters to begin in the open state! The YFE must be on with the bias dither initially enabled!')
        if np.sum(TTL_shutter.Status(display=False)[-1]) > 0:
            print('Warning! The shutters don\'t all appear to be open! ',end='',flush=True);TTL_shutter.Status(display=True);
        else:
            print('(Shutters seem OK...)')
        if not YFE.OnCheck(display=False):
            print('Warning! The YFE doesn\'t appear to be on! ',end='',flush=True);YFE.OnCheck(display=True);
        else:
            print('(YFE emission seems OK...)')
        if np.sum(YFE.Get(display=False)) < 550:
            print('Warning! The YFE doesn\'t appear to be turned up! ');YFE.Get(display=True);
        else:
            print('(YFE current seems OK...)')
        if MBC.MBCmodecheck() != 0:
            print('(Warning! The MBC doesn\'t appear to be in AUTO mode!')
        else:
            print('(MBC mode seems OK...)')
        print('Are you sure you are ready to proceed? [enter y/n]',end='',flush=True)
        checkprompt=input();
        if checkprompt.lower() != 'y':
            print('Try again later then!');
            return
        else:
            print('OK, I hope you know what you\'re doing!')
        HWP.On('all',set_T=1)
        armlist=['AB','EF','GH','IJ']
        #YFEoff();YFEon();
        GLOBAL.MBCmode.put(1)#set MAN mode on MBC
        if np.sum(YFE.Get(display=False)) < 100:
            print('Check YFE before optimizing!')
        optwvfm=pickle.load(open(LPL.psfilepath()+'opttrace.p','rb'));
        try:
            oldwvfm=HAWG().ReadPulseHeights();
            HAWG().WritePulseHeights(optwvfm);
        except:
            print('Failed to read old waveform/load new waveform!')
            return
        HWPpvlist=[GLOBAL.HWPAB, GLOBAL.HWPEF, GLOBAL.HWPGH, GLOBAL.HWPIJ];
        print('Closing all shutters...')
        TTL_shutter.Toggle('closeall',display=False);#close all the shutters
        time.sleep(4)
        GLOBAL.EVRLPLSSEC.put(43);GLOBAL.EVRLPLSSEN.put(1);#enable these...
        try:
            tempchk1=LOSC('a').rch(1);time.sleep(.15);tempchk2=LOSC('a').rch(1);
            if np.sum(np.abs(tempchk1-tempchk2))<1e-6:
                print('Warning: scope trace doesn\'t appear to be updating, please check scope! Abort? [enter y/n]')
                checkprompt=input();
                if checkprompt.lower() != 'y':
                    print('Try again later then!');
                    HAWG().WritePulseHeights(oldwvfm);
                    return
                else:
                    print('OK, I hope you know what you\'re doing!')
        except:
            print('Scope error, check scope status! Aborting...')
            HAWG().WritePulseHeights(oldwvfm);
            return
        startposlist=[HWPrbv.get() for HWPrbv in HWPpvlist];
        newposlist=startposlist[:]
        for ii in range(4):
            if armlist[ii] in armsQ:#only prep the stage if it's going to be used
                HWPpvlist[ii].put(startposlist[ii]+(-20.0+4.0*0))
        currentshutter=0;#trying to re-open a shutter in case of failure...
        #set up all the plotting stuff
        plt.ion()
        fig,axs=plt.subplots(2,2,gridspec_kw={'hspace':0.4,'wspace':0.3})
        xdat=[[startposlist[ii]+(-20.0+4.0*(jj)) for jj in range(11)] for ii in range(4)]
        ydat=[[0]*11 for ii in range(4)]
        ax1,=axs[0,0].plot(xdat[0],ydat[0]); axs[0,0].set_xlabel('AB'); plt.pause(0.01);
        ax2,=axs[0,1].plot(xdat[1],ydat[1]); axs[0,1].set_xlabel('EF'); plt.pause(0.01);
        ax3,=axs[1,0].plot(xdat[2],ydat[2]); axs[1,0].set_xlabel('GH'); plt.pause(0.01);
        ax4,=axs[1,1].plot(xdat[3],ydat[3]); axs[1,1].set_xlabel('IJ'); plt.pause(0.01);
        axss=[ax1,ax2,ax3,ax4]
        stepQ=1.0;rangeQ=20.0;
        try:
            SLA=LOSC('A');SLA.Open();#changed to LecroyA since repair
            for ii in range(4):
                if armlist[ii] in armsQ:
                    print('Begin optimizing '+armlist[ii]+'... ',end='',flush=True);
                    hwparmdatax,hwparmdatay=[],[]
                    TTL_shutter.Toggle('open'+armlist[ii],display=False);currentshutter=ii;time.sleep(4);print('Shutter opened!');#open one shutter
                    for jj in range(int(1+(2*rangeQ/stepQ))):
                        print('.',end='',flush=True)
                        HWPpvlist[ii].put(startposlist[ii]+(-rangeQ+stepQ*(jj)));time.sleep(4);#step to new position#was 2.5
                        curr_x=HWPpvlist[ii].get();curr_y=np.max(SLA._rch(3));time.sleep(.15);#in testing, max is more stable than sum
                        if curr_y > 0.005:#threshold so don't skew fit with noise; max is ~~10x this
                            hwparmdatax.append(curr_x);hwparmdatay.append(curr_y);#save x and y
                        print('.',end='',flush=True)
                        ydat[ii][jj]=curr_y
                        axss[ii].set_data(xdat[ii],ydat[ii])
                        axs[ii//2,ii%2].set_ylim((min(ydat[ii]),max(ydat[ii])))
                        plt.pause(0.01)
                        #axs[ii//2,ii%2].autoscale(True,True,True)
                        fig.canvas.draw_idle()
                        plt.pause(0.01)
                    print('*')
                    qfit=np.polyfit(hwparmdatax,hwparmdatay,2);newpos=qfit[1]/(-2*qfit[0]);#find fit and new max
                    if np.abs(startposlist[ii]-newpos)<0.85*rangeQ:
                        HWPpvlist[ii].put(newpos);newposlist[ii]=newpos;
                        print('HWP position on arm '+armlist[ii]+' changed from '+str(round(startposlist[ii],4))+' to '+str(round(newpos,4)))
                    else:
                        print('Failed! New HWP position on arm '+armlist[ii]+' seems too far off... '+str(round(newpos,4))+' from '+str(round(startposlist[ii],4))+'... Restoring...')
                        HWPpvlist[ii].put(startposlist[ii])
                    TTL_shutter.Toggle('close'+armlist[ii],display=False);currentshutter=0;#close that shutter;
                    #xpq=np.arange(startposlist[ii]+(-rangeQ+stepQ*(-1)),startposlist[ii]+(-rangeQ+stepQ*int(1+(2*rangeQ/stepQ))),.1);
                    qfitp=np.poly1d(qfit);
                    axs[ii//2,ii%2].plot(xdat[ii],qfitp(xdat[ii]))
                    axs[ii//2,ii%2].relim();plt.pause(0.01);
                    axs[ii//2,ii%2].autoscale(True,True,True)
                    fig.canvas.draw_idle()
                    plt.pause(0.01)    
                    #efc.epllxy([[hwparmdatax,hwparmdatay],[xpq,qfitp(xpq)]],xlb=armlist[ii])
                else:
                    print('Skipping '+armlist[ii]+'...')
                    pass
            SLA.Close();time.sleep(.15);#changed to LeCroyA
        except:
            print('Failed! Restoring original values and attempting to re-open most-recent shutter... you should verify!')
            SLA.Close();time.sleep(.15);#changed to LeCroyA
            if currentshutter > 0:
                TTL_shutter.Toggle('open'+armlist[currentshutter],display=False);
            for ii in range(4):
                HWPpvlist[ii].put(startposlist[ii]);newposlist[ii]=startposlist[ii];
        time.sleep(2);#need time so that last shutter trigger ends before trying to open IJ
        try:
            HAWG().WritePulseHeights(oldwvfm);
        except:
            print('Error! Check waveform!')
        GLOBAL.EVRLPLSSEN.put(0);#disable PC before re-opening shutters
        datestamp=int(datetime.now().strftime('%Y%m%d%H%M%S'))
        HWPlog=pickle.load(open(LPL.psfilepath()+'HWP_opt_log.p','rb'))
        HWPlog.append([datestamp,[newposlist[ii] for ii in range(4)]])
        efc.pickledump2(HWPlog,LPL.psfilepath()+'HWP_opt_log.p')
        TTL_shutter.Toggle('openall',display=False);#open all the shutters
        MBC.resetMBC();YFE.SetAll(True);#reset bias...
        plt.ioff()
        motnamelist=[GLOBAL.HWPABoff, GLOBAL.HWPEFoff, GLOBAL.HWPGHoff, GLOBAL.HWPIJoff]#.OFF
        #add adjustment to offset automatically....
        for temppv in motnamelist:
            tempval=temppv.get();
            temppv.put(tempval-newposlist[ii]);





class Stage:
    def NewportInitRefAll():
        ipvlist=[GLOBAL.XPS1IALL, GLOBAL.XPS2IALL, GLOBAL.XPS3IALL, GLOBAL.XPS4IALL]
        rpvlist=[GLOBAL.XPS1RALL, GLOBAL.XPS2RALL, GLOBAL.XPS3RALL, GLOBAL.XPS4RALL]        
        try:
            for eapv in ipvlist:
                eapv.put(1)
            efc.dotsleep(10)
            for eapv in rpvlist:
                eapv.put(1)
        except:
            print('Failed!') 
                  
    def SmarAct():
        pass
                  
    def Snapshot():
        pass
                  
    def Restore():
        pass
                  



                  
class Timing:
    def fstiming():
        pass
                  
    def Nstiming():
        pass                  
                  
    def EVR():
        pass
                  
    def Vitara():
        pass
                  




class CAM:
    @classmethod
    def quickCAM(cls,CAMreq,ImageNo=2,xlb='none',ylb='none',LIVE=False,MAXLOOPS=25):#like MEC:GIGE:31
        PVhead=cls.CAMname(CAMreq)
        try:#if True:#try:
            tres1=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
            tres2=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
            twf=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArrayData')
            if len(twf) != tres1*tres2:
                twf=list(twf)+(tres1*tres2-len(twf))*[0]
            fig,axs=plt.subplots(1,1)
            ax1=axs.imshow(np.array_split(np.array(twf),tres1));
            axs.axes.xaxis.set_ticklabels([]);
            axs.axes.yaxis.set_ticklabels([]);
            axs.tick_params(direction='in'); 
            axs.set_ylabel(CAMreq); 
        except:#if True:#except:
            print('Failed!')
            return
        fig.tight_layout()
        fig.show();
        waittime=.01;
        plt.pause(waittime);
        time.sleep(waittime)
        loopcount=0
        if LIVE:
            while loopcount<MAXLOOPS:
                try:
                    cls.CAMconfig(PVhead,LIVE=True)
                    twf=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArrayData')
                    if PVhead == 'MEC:GIGE:29':
                        tres1=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
                        tres2=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize0_RBV')
                        tres3=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
                        tresL=sorted([tres1,tres2,tres3],reverse=True) 
                        if len(twf) != tresL[0]*tresL[1]:
                            twf = list(twf) + (tresL[0]*tresL[1] - len(twf))*[0]
                            tres1=tres3
                    ax1.set_data(np.array_split(np.array(twf),tres1));
                    fig.canvas.draw_idle()
                    plt.pause(waittime)
                    time.sleep(waittime)
                except:
                    print('Error occured when plotting {}!'.format(CAMreq))
                loopcount+=1    
        cls.CAMconfig(PVhead,LIVE=False)
        return
    
    
    def CAMname(GIGEnam='none',returnAll=False):
        NickNameList=[['Legend','Regen'], ['StrInA','TopasA'], ['StrInB','TopasB'], ['MPA1In','MPA1A'], ['MPA1Out','MPA1B'], ['MPA2In','MPA2A'], ['MPA2Xtal','MPA2F'], ['MPA2Out','MPA2B'], ['CompIn', 'CompA'], ['CompOutFF', 'CompBFF'], ['CompOutNF', 'CompBNF'], ['Trap', 'Mousetrap']]
        SPLNameList=['MEC_SPL_3', 'MEC_SPL_1', 'MEC_SPL_7', 'MEC_SPL_2', 'MEC_SPL_4', 'MEC_SPL_5', 'GigE17_TimeTool_Diag', 'MEC_SPL_6', 'MEC_SPL_9', 'MEC_SPL_10', 'MEC_SPL_11', 'MEC_SPL_8']
        PVNameList=['MEC:GIGE:24', 'MEC:GIGE:22', 'MEC:GIGE:28', 'MEC:GIGE:23', 'MEC:GIGE:25', 'MEC:GIGE:26', 'MEC:GIGE:17', 'MEC:GIGE:27', 'MEC:GIGE:29', 'MEC:GIGE:30', 'MEC:GIGE:31', 'MEC:GIGE:16']
        SmarActList=['SM0','SM1','SM2','SM3','SM4','SM5','SM6','SM7','SM8','XPS Mid','XPS In','None']
        pos=np.where(np.array([list(map(lambda name: name.casefold(), NickName)) for NickName in  np.array(NickNameList)])==GIGEnam.casefold())
        if len(pos[0]) > 0:
            if returnAll==True:
                return [NickNameList[pos[0][0]][0], SPLNameList[pos[0][0]], PVNameList[pos[0][0]]] 
            else:
                return PVNameList[pos[0][0]]
        pos2=np.where(np.array(list(map(lambda name: name.casefold(), np.array(SPLNameList))))==GIGEnam.casefold())
        if len(pos2[0]) > 0:
            if returnAll==True:
                return [NickNameList[pos2[0][0]][0], SPLNameList[pos2[0][0]], PVNameList[pos2[0][0]]] 
            else:
                return PVNameList[pos2[0][0]]
        pos3=np.where(np.array(list(map(lambda name: name.casefold(), np.array(PVNameList))))==GIGEnam.casefold())
        if len(pos3[0]) > 0:
            if returnAll==True:
                return [NickNameList[pos3[0][0]][0], SPLNameList[pos3[0][0]], PVNameList[pos3[0][0]]] 
            else:
                return PVNameList[pos3[0][0]]
        else:
            print('{:<30} {:<21} {:<11}  {:<10}'.format('Camera NickName', 'SPL CAM Name', 'PV Name', 'Motor'))
            for ii in range(len(NickNameList)):
                print('{:<10} or   {:<12}   {:<21} {:<11}  {:<10}'.format(str(NickNameList[ii][0]), str(NickNameList[ii][1]), str(SPLNameList[ii]), str(PVNameList[ii]), str(SmarActList[ii])))
            return False

    @classmethod
    def CAMview(cls, *CAMargs,ImageNo=2,LIVE=False,MAXLOOPS=10):
        if CAMargs == ('all',):
            CAMargs = ('Regen', 'Trap', 'StrInA', 'StrInB', 'MPA1In', 'MPA1Out', 'MPA2In', 'MPA2Out', 'MPA2Xtal', 'CompIn', 'CompOutNF', 'CompOutFF')
        if len(CAMargs) == 1:
            cls.quickCAM(*CAMargs,LIVE=LIVE,MAXLOOPS=MAXLOOPS)
            return
        plt.ion()
        subply=len(CAMargs)//2 + len(CAMargs)%2;subplx=2;
        fig,axs=plt.subplots(subply,subplx,figsize=(5,2*subply));
        axss=[];tres1L=[];tPVheadL=[]
        for ii in range(len(CAMargs)):
            tidx=(ii//2,ii%2) if len(CAMargs) > 2 else (ii%2)
            axs[tidx].axes.xaxis.set_ticklabels([]);
            axs[tidx].axes.yaxis.set_ticklabels([]);
            axs[tidx].tick_params(direction='in'); 
            axs[tidx].set_ylabel(CAMargs[ii]); 
            try:
                tPVhead=cls.CAMname(CAMargs[ii])
                tPVheadL.append(tPVhead);
                tres1=efc.rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
                tres2=efc.rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArraySize0_RBV')
                tres3=efc.rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
                tresL=sorted([tres1,tres2,tres3],reverse=True) 
                twf=efc.rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArrayData')
                if len(twf) != tresL[0]*tresL[1]:
                    twf = list(twf) + (tresL[0]*tresL[1] - len(twf))*[0]
                    tres1=tres3
                tempax=axs[tidx].imshow(np.array_split(np.array(twf),tres1)); 
                tres1L.append(tres1);
                axss.append(tempax)
            except:
                print('Error occured when plotting {}!'.format(CAMargs[ii]))
        if (len(CAMargs) > 2) and (len(CAMargs)%2 > 0):
            iit=len(CAMargs)
            tidx=(iit//2,iit%2) 
            axs[tidx].axis('off');
        fig.tight_layout()
        plt.show();
        waittime=.01;
        plt.pause(waittime);
        time.sleep(waittime)
        loopcount=0
        if LIVE:
            while loopcount<MAXLOOPS:
                for ii in range(len(CAMargs)):
                    tidx=(ii//2,ii%2) if len(CAMargs) > 2 else (ii%2)
                    try:
                        twf=efc.rPV(tPVheadL[ii]+':IMAGE'+str(ImageNo)+':ArrayData')
                        if tPVheadL[ii] == 'MEC:GIGE:29':
                            tres1=efc.rPV(tPVheadL[ii]+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
                            tres2=efc.rPV(tPVheadL[ii]+':IMAGE'+str(ImageNo)+':ArraySize0_RBV')
                            tres3=efc.rPV(tPVheadL[ii]+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
                            tresL=sorted([tres1,tres2,tres3],reverse=True) 
                            if len(twf) != tresL[0]*tresL[1]:
                                twf = list(twf) + (tresL[0]*tresL[1] - len(twf))*[0]
                                tres1=tres3
                            tres1L[ii]=tres1
                        axss[ii].set_data(np.array_split(np.array(twf),tres1L[ii])); 
                        fig.canvas.draw_idle()
                        plt.pause(waittime)
                        time.sleep(waittime)
                    except:
                        print('Error occured when plotting {}!'.format(CAMargs[ii]))
                loopcount+=1
                
# =============================================================================
#     def _tiptilt(cls,camPVhead):#allow toggle btw xyz, getch btw step/etc, plot view of beam on CRIX:GIGE:06
#         plt.ion() 
#         fig,axs=plt.subplots(1,1)
#         if camPVhead == 'sg':
#             camPV='head'; motorPVx='head'; motorPVy='head';
#         else:
#             print('No head found!');#return
#         camPV=EpicsSignal(camPVhead+':IMAGE2:ArrayData');
#         camPVres=cls.rPV(camPVhead+':IMAGE2:ArraySize1_RBV');
#         pdat=np.array(camPV.get())
#         ax1=axs.imshow(np.array_split(pdat,camPVres),origin='upper',cmap='prism');plt.pause(0.01);#gist_ncar
#         plt.show();
#         waittime=.01;
#         plt.pause(waittime);
#         time.sleep(waittime)
#         #axs.imshow(Zref,alpha=0.1,origin='lower')
#         while True:
#             pdatn=np.array(camPV.get())
#             if True:#not np.lib.arraysetops.setdiff1d(pdatn[0], pdat[0]):
#                 ax1.set_data(np.array_split(pdat,camPVres))
#                 fig.canvas.draw_idle()
#                 plt.pause(waittime)
#                 pdat=pdatn
#                 time.sleep(waittime);
#             else:
#                 print('Same frame!')
# =============================================================================
    def setdynxhair(PVhead):
        pass

    @classmethod
    def CAMconfig(cls, CAMreq, RefXhairXY=[-1,-1], InstructionStr='', MisalignmentTolerance=-1, ImageNo=2, RefXhairXYsize=[-1,-1],LIVE=False):
        #default RefXhairXY=[-1,-1] indicates not to change anything, it was set up previously
        #default InstructionStr='' indicates not to change anything, it was set up previously
        #default MisalignmentTolerance=-1 indicates not to change anything, it was set up previously
        #default RefXhairXYsize=[-1,-1] indicates not to change anything, it was set up previously
        PVhead=cls.CAMname(CAMreq)
        efc.wPV('{}:Acquire'.format(PVhead),0)#hopefully prevents IOC crashes, per TylerJ
        ArraySize=[efc.rPV('{}:IMAGE{}:ArraySize0_RBV'.format(PVhead,ImageNo)),efc.rPV('{}:IMAGE{}:ArraySize1_RBV'.format(PVhead,ImageNo))]
        NameList=['Ref X-hair','DynCentroid','CAMname','Instructions','TimeStamp','Stats','Counts',
                  'AlignmentOK'+str(str(MisalignmentTolerance).zfill(3) if MisalignmentTolerance > 0 else efc.rPV('{}:IMAGE{}:Cross4:Name'.format(PVhead,ImageNo))[-3:])]
        ShapeList=[0,2,3,3,3,3,3,3]#0=Cross,1=Rect,2=Text,3=Ellipse;;;DOCUMENTATION WRONG!! 2=Ellipse, 3=Text!!!!
        #each overlay needs the following information:
        #[SizeXLink.DOL, SizeX, PositionXLink.DOL, PositionX, CenterXLink.DOL, CenterX, 
        # SizeYLink.DOL, SizeY, PositionYLink.DOL, PositionY, CenterYLink.DOL, CenterY, 
        # TimeStampFormat]
        #Size affects Position/Center, so set it first
        #specify for each overlay in ArgsList
        #in specifying SizeX or SizeY: 9999 indicates it should be set according to string length
        OverlayName='{}:IMAGE{}:{}'.format(PVhead,ImageNo,'XXXX')
        ArgsList=[
        #Ref X-hair
        [OverlayName+':SizeXDummy', RefXhairXYsize[0], '', 0, OverlayName+':CenterXDummy', RefXhairXY[0], 
         OverlayName+':SizeYDummy', RefXhairXYsize[1], '', 0, OverlayName+':CenterYDummy', RefXhairXY[1]],
        #DynCentroid
        [PVhead+':Stats2:SigmaX_RBV CP', 0, '', 0, PVhead+':Stats2:CentroidX_RBV CP', 0,
         PVhead+':Stats2:SigmaY_RBV CP', 0, '', 0, PVhead+':Stats2:CentroidY_RBV CP', 0],
        #CAMname
        [OverlayName+':SizeXDummy', 9999, OverlayName+':PositionXDummy', 10, '', 0, 
         OverlayName+':SizeYDummy', 9999, OverlayName+':PositionYDummy', 10, '', 0,
         '{:<10.9} {:<11.10} {:<11.11}'.format(*cls.CAMname(PVhead,returnAll=True))],
        #Instructions
        [OverlayName+':SizeXDummy', 9999, OverlayName+':PositionXDummy', 10, '', 0, 
         OverlayName+':SizeYDummy', 9999, OverlayName+':PositionYDummy', -20+ArraySize[1], '', 0, 
         InstructionStr if len(InstructionStr)>1 else efc.rPV('{}:IMAGE{}:Box4:DisplayText'.format(PVhead,ImageNo))],
        #TimeStamp
        [OverlayName+':SizeXDummy', 9999, OverlayName+':XPosDummy', -205+ArraySize[0], '', 0, 
         OverlayName+':SizeYDummy', 9999, OverlayName+':YPosDummy', -20+ArraySize[1], '', 0, 
         '%Y-%m-%d %H:%M:%S.%03f'],
        #Stats
        [OverlayName+':SizeXDummy', 9999, OverlayName+':XPosDummy', -300+ArraySize[0], '', 0, 
         OverlayName+':SizeYDummy', 9999, OverlayName+':YPosDummy', 10, '', 0, 
         'Stats: X:{:<3.3} Y:{:<3.3} sX:{:<3.3} sY:{:<3.3}'.format(str(efc.rPV(PVhead+':Stats2:CentroidX_RBV')),
                                                                   str(efc.rPV(PVhead+':Stats2:CentroidY_RBV')), 
                                                                   str(efc.rPV(PVhead+':Stats2:SigmaX_RBV')), 
                                                                   str(efc.rPV(PVhead+':Stats2:SigmaY_RBV')))],
        #Counts
        [OverlayName+':SizeXDummy', 9999, OverlayName+':XPosDummy', -282+ArraySize[0], '', 0, 
         OverlayName+':SizeYDummy', 9999, OverlayName+':YPosDummy', 25, '', 0, 
         'Max:{:<5} Counts:{:<12}'.format(efc.rPV(PVhead+':Stats2:MaxValue_RBV'), int(efc.rPV(PVhead+':Stats2:Net_RBV')))],
        #Alignment OK
        [OverlayName+':SizeXDummy', 9999, '', 0, OverlayName+':XPosDummy', efc.rPV('{}:IMAGE{}:Box1:CenterX_RBV'.format(PVhead,ImageNo)), 
         OverlayName+':SizeYDummy', 9999, '', 0, OverlayName+':YPosDummy', efc.rPV('{}:IMAGE{}:Box1:CenterY_RBV'.format(PVhead,ImageNo))+efc.rPV('{}:IMAGE{}:Box1:SizeY_RBV'.format(PVhead,ImageNo))+40, 
         'AlignChk']
        ]
    
        for ii in range(8):#eight overlays possible
            if (LIVE==True and ii<5):
                pass
            else:
                try:
                    #turn on stats...
                    efc.wPV(PVhead+':Stats{}:EnableCallbacks'.format(ImageNo), 1)
                    efc.wPV(PVhead+':Stats{}:ComputeCentroid'.format(ImageNo), 1)
                    efc.wPV(PVhead+':Stats{}:ComputeStatistics'.format(ImageNo), 1)
                    #set up everything...
                    OverlayNameField='{}{}'.format('Box' if ii<4 else 'Cross',(ii%4)+1)
                    OverlayName='{}:IMAGE{}:{}'.format(PVhead,ImageNo,OverlayNameField)
                    efc.wPV(OverlayName+':Name', NameList[ii])#simple description field
                    efc.wPV(OverlayName+':Use', 1)#0=off,1=on
                    efc.wPV(OverlayName+':Shape', ShapeList[ii])#0=Cross,1=Rect,2=Text,3=Ellipse;;;DOCUMENTATION WRONG!! 2=Ellipse, 3=Text!!!!
                    efc.wPV(OverlayName+':DrawMode', 1)#0=set,1=XOR
                    efc.wPV(OverlayName+':Green', 2000 if (ii < 5 or LIVE == True) else 0)#set mono color; don't let last ones appear yet
                    if ShapeList[ii] == 0:
                        Width=3#1
                    elif ShapeList[ii] == 1:
                        Width=3
                    elif ShapeList[ii] == 2:
                        Width=3#1
                    elif ShapeList[ii] == 3:
                        Width=3
                    else:
                        pass
                        
                    efc.wPV(OverlayName+':WidthX', Width)#width param for shape lines
                    efc.wPV(OverlayName+':WidthY', Width)#width param for shape lines
        
                    #all params needed for size and positioning of overlays 
                    #try to filter out writing zeroes, as this seems to recalculate stuff inadvertantly??
                    #also: display text repeats if size is larger than needed for that string
                    #   -->use -1 to indicate that the size needs to be calculated from the instruction length
                    #      in order to show only the intended message 
                    efc.wPV(OverlayName+':SizeXLink.DOL', ArgsList[ii][0].replace('XXXX',OverlayNameField))
                    if ArgsList[ii][1] > 0:
                        if ii==7:#special case for Alignment OK
                            ArgsList[ii][12]='Alignment {}'.format('NOT OK' if np.sqrt((efc.rPV(PVhead+':Stats2:CentroidX_RBV')-efc.rPV('{}:IMAGE{}:Box1:CenterX_RBV'.format(PVhead,ImageNo)))**2+(efc.rPV(PVhead+':Stats2:CentroidY_RBV')-efc.rPV('{}:IMAGE{}:Box1:CenterY_RBV'.format(PVhead,ImageNo)))**2) > int(efc.rPV('{}:IMAGE{}:Cross4:Name'.format(PVhead,ImageNo))[-3:]) else 'OK')
                        efc.wPV(OverlayName+':SizeX', ArgsList[ii][1] if ArgsList[ii][1] < 9999 else 9*len(ArgsList[ii][12]))
                    efc.wPV(OverlayName+':PositionXLink.DOL', ArgsList[ii][2].replace('XXXX',OverlayNameField))
                    if ArgsList[ii][3] > 0:
                        efc.wPV(OverlayName+':PositionX', ArgsList[ii][3])
                        
                    efc.wPV(OverlayName+':CenterXLink.DOL', ArgsList[ii][4].replace('XXXX',OverlayNameField))
                    if ArgsList[ii][5] > 0:
                        efc.wPV(OverlayName+':CenterX', ArgsList[ii][5])
                    
                    efc.wPV(OverlayName+':SizeYLink.DOL', ArgsList[ii][6].replace('XXXX',OverlayNameField))
                    if ArgsList[ii][7] > 0:
                        efc.wPV(OverlayName+':SizeY', ArgsList[ii][7] if 9999 > ArgsList[ii][7] >= 0 else 20)
                    efc.wPV(OverlayName+':PositionYLink.DOL', ArgsList[ii][8].replace('XXXX',OverlayNameField))
                    if ArgsList[ii][9] > 0:
                        efc.wPV(OverlayName+':PositionY', ArgsList[ii][9])
                    efc.wPV(OverlayName+':CenterYLink.DOL', ArgsList[ii][10].replace('XXXX',OverlayNameField))
                    if ArgsList[ii][11] > 0:
                        efc.wPV(OverlayName+':CenterY', ArgsList[ii][11])
                    
                    if ShapeList[ii] == 3:#only needed for text boxes; NOT 2!!!
                        #use :DisplayText??
                        if ArgsList[ii][12][0] != '%':
                            if (len(ArgsList[ii][12]) > 1):
                                efc.wPV(OverlayName+':DisplayText', ArgsList[ii][12])
                        else:
                            efc.wPV(OverlayName+':TimeStampFormat', ArgsList[ii][12])
                        efc.wPV(OverlayName+':Font',3)
                except:
                    print('Error setting up {}!'.format(NameList[ii]))

# =============================================================================
#     def GigE_toggle_trigger(ENBINHReq, GigEreq='all', RepRateReq=5):
#         """
#         INHall, ENBall
#         5, 0
#         """
#         #check to see if UNI_EVR output is enabled or disabled before changing trigger timing of UNI_EVR
#         #may also check for 
#         # individually set free run or fixed rate or ext trig in? use CAMconfig for this?
#         if UNI_TRIG_IN_EVR in [44,45,46]:
#             EFdelay=0e-3
#         elif UNI_TRIG_IN_EVR in [177]:
#             EFdelay=100e-3
#         else:
#             print('Unanticipated UNI_TRIG_IN_EVR case!')
#             
#         if command == 'INH':
#             DG8_EF_Polarity = 'POS'
#         elif command == 'ENA':
#             DG_8_EF_Polarity = 'NEG'
#         else:
#             print('Unanticipated command case!')
#         return
# =============================================================================

    @classmethod
    def CAMconfigCurrent(cls):
        cls.CAMconfig(CAMreq='Legend',RefXhairXY=[359,251],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[100,100],LIVE=False,InstructionStr='DO NOT tune to xhair!')
        cls.CAMconfig(CAMreq='StrInA',RefXhairXY=[335,257],MisalignmentTolerance=25,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM0 (NOT SM1!) to put centroid on xhair!')
        cls.CAMconfig(CAMreq='StrInB',RefXhairXY=[334,243],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM2 to put centroid on xhair!')
        cls.CAMconfig(CAMreq='MPA1In',RefXhairXY=[152,147],MisalignmentTolerance=15,ImageNo=2, RefXhairXYsize=[25,25],LIVE=False,InstructionStr='Use SM3 to put centroid on xhair!')
        cls.CAMconfig(CAMreq='MPA1Out',RefXhairXY=[411,249],MisalignmentTolerance=15,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM4 to maximize spot!')
        cls.CAMconfig(CAMreq='MPA2In',RefXhairXY=[366,276],MisalignmentTolerance=15,ImageNo=2, RefXhairXYsize=[50,50],LIVE=False,InstructionStr='Use SM5(V) to put centroid on xhair!')
        cls.CAMconfig(CAMreq='MPA2Out',RefXhairXY=[385,259],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[190,190],LIVE=False,InstructionStr='Use SM7 (and SM6) to put centroid on xhair!')
        cls.CAMconfig(CAMreq='CompIn',RefXhairXY=[276/2,252/2],MisalignmentTolerance=50,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM8 to put centroid on xhair!')
        cls.CAMconfig(CAMreq='CompOutFF',RefXhairXY=[292,300],MisalignmentTolerance=25,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use XPS2 Mid (and Input) Mirror to align to xhair!')
        cls.CAMconfig(CAMreq='CompOutNF',RefXhairXY=[364,277],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use XPS2 Input(/Mid) Mirror to align to xhair!')






                  
class TTL_shutter:
    def Status(display=True):
        pvstatuslist=['MEC:LAS:FLOAT:'+str(ii) for ii in range(14,21)];
        statuslist=[]
        for eapv in pvstatuslist:
            temppv=EpicsSignal(eapv);
            statuslist.append(int(temppv.get()))#0=open,1=closed
        if display:
            shutlist=['AB','EF','GH','IJ','WEST 527','EAST 527','REGEN']
            print('(0=out,1=in) '+', '.join([ea_shut+':'+str(ea_stat) for ea_shut,ea_stat in zip(shutlist,statuslist)]))
        return statuslist

    @classmethod
    def Toggle(cls,ArmStrQ,display=True):
        """x"""
        #could add GLOBALS.TTLAB etc later
        AllStrq='abefghijwwxxzz';#ww is WEST 527, xx is EAST 527, zz is REGEN SHUT
        ArmStrQ=ArmStrQ.lower().replace('all','abefghijwwxx')
        currTTLstate=cls.Status(display=False)[:-1];
        if ArmStrQ.lower()[:4] == 'clos':
            ArmStrQ=''.join(sorted(ArmStrQ[5:].lower()))
            ArmStrQ=''.join([shutt for shutt,oshutt in zip(re.findall('..','abefghijwwxx'),currTTLstate) if oshutt==0 and shutt in ArmStrQ])
        elif ArmStrQ.lower()[:4] == 'open':
            ArmStrQ=''.join(sorted(ArmStrQ[4:].lower()))
            ArmStrQ=''.join([shutt for shutt,oshutt in zip(re.findall('..','abefghijwwxx'),currTTLstate) if oshutt==1 and shutt in ArmStrQ])
        else:
            ArmStrQ=''.join(sorted(ArmStrQ.lower()))
        if display:
            print('Initially: ',end='',flush=True)
        statuslist=cls.Status(display=display)
        for ii in range(len(AllStrq)//2):
            if (AllStrq[2*ii] in ArmStrQ.lower()) or (AllStrq[2*ii+1] in ArmStrQ.lower()):
                temppv1=EpicsSignal('MEC:LAS:TTL:0'+str(ii+1));temppv2=EpicsSignal('MEC:LAS:FLOAT:'+str(ii+14));
                temppv1.put(1);
        temppvdur=EpicsSignal('MEC:EK9K1:BO5:1.HIGH');
        time.sleep(temppvdur.get()*0.5);
        for ii in range(len(AllStrq)//2):
            if (AllStrq[2*ii] in ArmStrQ.lower()) or (AllStrq[2*ii+1] in ArmStrQ.lower()):
                temppv1=EpicsSignal('MEC:LAS:TTL:0'+str(ii+1));temppv2=EpicsSignal('MEC:LAS:FLOAT:'+str(ii+14));
                if temppv1.get():
                    temppv2.put((1+temppv2.get())%2);
                else:
                    print('Warning: no shutter toggle detected; status not updated! '+AllStrq[2*ii:2*ii+2])
        if display:
            print('Finally:   ',end='',flush=True)
        statuslist=cls.Status(display=display)
        return statuslist

    def Refresh():
        pvlist=['MEC:LAS:FLOAT:'+str(ii) for ii in range(14,21)];
        inddesclist=['AB','EF','GH','IJ','WEST (ABEF)','EAST (GHIJ)','Regen']
        desclist=[inddesc+' shutter state' for inddesc in inddesclist]
        valulist=[0,0,0,0,0,0,0];
        print('This will reset the shutter counter! Are you sure you want to continue? [y/n]',end='',flush=True)
        checkprompt=input();
        if checkprompt.lower() != 'y':
            print('Try again later then!');
            return
        else:
            print('OK, I hope you know what you\'re doing!')
        print('Resetting shutter PV statuses... Please insure all shutters are actually open!!')
        for jj in range(len(pvlist)):
            temppv1=EpicsSignal(str(pvlist[jj]+'.DESC'));temppv2=EpicsSignal(pvlist[jj]);
            temppv1.put(desclist[jj]);temppv2.put(valulist[jj]);





class DG645:
# =============================================================================
#     def Refresh():
#         pvlist=[['MEC:LAS:DDG:0'+str(numii)+':'+chii+'DelaySI.DESC' for chii in ['a','c','e','g']] for numii in [1,2,6,8]];
#         desclist=[['A:PS LATE','C:INH PS EARLY','E:unused','G: unused'],['A:PS EARLY','C:unused','E:EvoHE1','G:EvoHE2'],['A:GaiaQSW','C:GaiaLamp','E:INH UNI','G:GigE TRIG IN'],['A:BIG UNI','C:small UNI','E:INH GigE','G:unused']]
#         for eachboxii in range(len(pvlist)):
#             for eachentryii in range(len(eachboxii)):
#                 temppv=EpicsSignal(pvlist[eachboxii][eachentryii])
#                 temppv.put(desclist[eachboxii][eachentryii])
# 
#     def Snapshot():
#         pass
#     
#     def DGrecordall():
#         fullstr=datetime.now().strftime('%Y%m%d_%H%M%S')+'\n'#date.today.strftime('%Y%m%d')
#         DGnam=['SLICER','MPA1','MASTER','UP/DOWN','STREAK','MPA2','USR','UNIBLITZ']
#         for DGno in [1,2,3,4,5,6,8]:
#             fullstr+=DGnam[DGno-1]+'\n'
#             for lett in ['a','b','c','d','e','f','g','h']:
#                 head='MEC:LAS:DDG:{:0>2}:{}DelaySI'.format(DGno,lett)
#                 tail='.DESC'
#                 tpv1=EpicsSignal(head);tpv2=EpicsSignal(head+tail);
#                 try:
#                     fullstr+='{}: {}'.format(tpv2.get(),tpv1.get())+'\n'
#                 except:
#                     fullstr+='ERROR'+'\n'
#         with open(str(LPL.psfilepath()+'DG/snapshot'+LPL.DateString()+'.txt'),'a') as out:
#             out.write(fullstr)
#     
#     def DGsaveref(Laser):
#         if Laser.lower() == 'spl':
#             DGnam=['SLICER','MPA1','MPA2','UNIBLITZ']
#             DGno=[1,2,6,8]
#         elif Laser.lower() == 'lpl':
#             DGnam=['MASTER','UP/DOWN','STREAK']
#             DGno=[3,4,5]
#         else:
#             print('Choose LPL or SPL!'); return
#         #MEC:LAS:DDG:03:gReferenceMO
#         #re.findall('(.+) \+ (.+)','A + 1.005e-6')
#         return (DGnam, DGno)
# =============================================================================
    pass





class SPL:#pointing, alignment, and other macros + automation routines
    pass





class UNIBLITZ:
# =============================================================================
#     def UNIBLITZconfig(modeReq):
#         """
#         Configures UNIBLITZ shutter state and triggers for different MEC SPL modes:
#         
#         modeReq=      6mm state  6mm TRIG?  65mm state  65mm TRIG?
#         ----------------------------------------------------------
#         'alignment'-->   OPEN       INH        OPEN        INH
#         'blocked'  -->  CLOSED      INH       CLOSED       INH
#         'MPA1SS'   -->  CLOSED    ENABLED      OPEN        INH
#         'MPA1Burst'-->  CLOSED**    INH**      OPEN        INH
#         'MPA2SS'   -->  CLOSED    ENABLED     CLOSED     ENABLED
#         'MPA2Burst'-->  CLOSED    ENABLED      OPEN        INH    
#         """
#         if   modeReq.casefold() == 'alignment':
#             UNI_toggle_trigger('INHall')
#             UNI_toggle_shutter('openall')
#         elif modeReq.casefold() == 'blocked':
#             UNI_toggle_trigger('INHall')
#             UNI_toggle_shutter('closeall')
#         elif modeReq.casefold() == 'mpa1ss':
#             UNI_toggle_trigger('INH65mm')
#             UNI_toggle_shutter('open65mm')
#             UNI_toggle_trigger('ENA06mm')
#             UNI_toggle_shutter('close06mm')
#         elif modeReq.casefold() == 'mpa1burst':
#             UNI_toggle_trigger('INHall')
#             UNI_toggle_shutter('open65mm')
#             UNI_toggle_shutter('close06mm')
#         elif modeReq.casefold() == 'mpa2ss':
#             UNI_toggle_trigger('ENAall')
#             UNI_toggle_shutter('closeall')
#         elif modeReq.casefold() == 'mpa2burst':
#             UNI_toggle_trigger('INH65mm')
#             UNI_toggle_shutter('open65mm')
#             UNI_toggle_trigger('ENA06mm')
#             UNI_toggle_shutter('close06mm')
#         else:
#             print('Please choose a valid configuration mode!')
#             print(UNIBLITZconfig.__doc__)
#         
#     def UNI_toggle_shutter(ENAINHReq):
#         """
#         [command=open|close][shutdia=06|65|all]mm
#         """
#         #get a status readback capability...
#         if shutdia=='all':
#             templist=['06','65']
#         else:
#             if shutdia in ['06','65']:
#                 templist=[shutdia]
#             else:
#                 print("Only accepts '06' or '65' or 'all' for shutter specification!")
#                 return False
#         if command=='open':
#             polarity='POS'#get actual val
#         elif command=='close':
#             polarity='NEG'#get actual val
#         else:
#             print("Only accepts 'open' or 'close' for shutter command!")
#             return False
#         for eashutt in templist:
#             wPV('PV:CH{}:Polarity'.format(chAB if eashutt=='65' else chCD), polarity)#get actual chan val and PV
#         return
#         
#     
#     def UNI_toggle_trigger(ENAINHReq, RepRateReq=0):
#         """
#         INHall, INH65mm == ENA06mm, ENBall
#         RepRateReq
#         all need to address UNI_TRIG_IN EVR channel, which needs to be reestablished... (stolen!)
#         """
#         #get a status readback capability...
#         if 'INH':
#             pass
#         return
# =============================================================================
    pass





class Spectrometer:#Qminis
    pass





class VISAR:#(laser, streak cameras, etc.)
    pass





class CtrlSys:
    def plzchkpv(inpvnam):
        try:
            qp=EpicsSignal(inpvnam);
            currval=qp.get()
            msgout = str(currval)+' vs oldval'
        except TimeoutError as err:
            msgout = 'Timeout!';err;
        return msgout

    def plzchksrv(inpvnam,printall=False):
        pvnam=re.findall('^([^\.]+)',inpvnam);
        if len(pvnam)>0:
            iocnam=re.findall('/reg/d/iocData/(.+)/iocInfo',os.popen('grep_pv '+pvnam[0]).read());
        else:
            pvnam=['Fail!']
            iocnam=[]
        if len(iocnam)>0:
            hostnam=re.findall('host: \'([^,]+)\', ', os.popen('grep_ioc '+iocnam[0]).read());
        else:
            iocnam=['Fail!']
            hostnam=[]
        if len(hostnam)>0:
            netstat=os.system('ping -c 1 -w2 '+hostnam[0]+' > /dev/null 2>&1')
            locstat=re.findall('Location: (.+)\n', os.popen('netconfig search '+hostnam[0]).read());
        else:
            hostnam=['Fail!']
            netstat=-1
            locstat=['Fail!']

        #msgout='PV: '+'{:<34}'.format(pvnam[0])+'IOC: '+'{:<26}'.format(iocnam[0])+'Host: '+'{:<16}'.format(hostnam[0])+'Host rack:'+'{:<19}'.format(locstat[0])+' Host Ping?: '+str('true ' if netstat==0 else 'false')
        msgout=['{:<34}'.format(pvnam[0]),'{:<26}'.format(iocnam[0]),'{:<19}'.format(hostnam[0]),'{:<25}'.format(locstat[0]),efc.cstr('{:<6}'.format(str('true ' if netstat==0 else 'false')),str('blink,r' if netstat!=0 else ''))]
        return msgout

    @classmethod
    def pv_checker(cls,las='lpl'):
        if las.lower() == 'lpl':
            qqpl=np.genfromtxt(LPL.psfilepath()+'_ps_pvlist_lpl.txt',delimiter='\n',dtype=str);pvmsg=[];
        elif las.lower() == 'spl':
            qqpl=np.genfromtxt(LPL.psfilepath()+'_ps_pvlist_spl.txt',delimiter='\n',dtype=str);pvmsg=[];
        else:
            print('No such laser!');return
        for eapv in qqpl:
            pvmsg.append(cls.plzchkpv(eapv))
        with multiprocessing.Pool() as pool:
            srvmsg=pool.map(cls.plzchksrv,qqpl)
        print(''.join(['{:<34}'.format('PV name'),'{:<26}'.format('IOC name'),'{:<19}'.format('Host name'),'{:<25}'.format('Host location'),'{:<6}'.format('Ping?'),'{:<15}'.format('PV value?')]))
        for ii in range(len(pvmsg)):
            print(''.join(srvmsg[ii])+pvmsg[ii])
        return #datq

    def plzchkcmp(cmpnam):
        try:
            netstat=os.system('ping -c 1 -w2 '+cmpnam[0]+' > /dev/null 2>&1')
            msgout=['{:<15}'.format(cmpnam[2])+'{:<28}'.format(cmpnam[0])+efc.cstr('{:<6}'.format(str('true ' if netstat==0 else 'false')),str('blink,r' if netstat!=0 else ''))]
        except:
            msgout=['{:<15}'.format(cmpnam[1])+'{:<28}'.format(cmpnam[0])+efc.cstr('{:<6}'.format('fail!'),'blink,r')]
        return msgout

# =============================================================================
#     def plzchkcmp2(cmpnam):
#         netstat=os.system('ping -c 1 -w2 '+cmpnam[0]+' > /dev/null 2>&1')
#         msgout=['{:<15}'.format(cmpnam[1])+'{:<28}'.format(cmpnam[0])+efc.cstr('{:<6}'.format(str('true ' if netstat==0 else 'false')),str('blink,r' if netstat!=0 else ''))]
#         return msgout
# =============================================================================
# =============================================================================
#     @staticmethod
#     def MECcompylist000():
#         qqip=['172.21.46.147','172.21.46.148','172.21.46.146','172.21.46.60','172.21.46.128','172.21.46.100', '172.21.46.120','172.21.46.159', '172.21.46.197','172.21.46.142','172.21.46.70','172.21.46.78', '172.21.46.71','172.21.46.88','172.21.46.198','172.21.46.213','172.21.46.215','172.21.46.136', '172.21.46.218','172.21.46.219','172.21.46.182','172.21.46.144'];
#         qqn=['evo1','evo2','gaia','lecroy1','lecroy2','lecroya','lecroyb','PIMikroMove','spider','spectrometer', 'tundra','topas','visar1','visar2','vitara','rga','emp','phasicslaptop','phasics1','phasics2','dacage','legend']
#         nmlist=['mec-las-laptop06','mec-las-laptop07','mec-las-laptop05','scope-ics-mectc1-1','scope-ics-meclas-lecroy01','scope-ics-meclas-lecroy-a','scope-ics-meclas-lecroy-b','mec-las-laptop09','mec-las-laptop11','mec-las-laptop01','win-ics-mec-tundra','mec-las-laptop12','win-ics-mec-visar1','win-ics-mec-visar2','mec-las-vitara','mec-rga-laptop','scope-ics-mec-tektronix','mec-phasics-laptop01','win-ics-mec-phasics01','win-ics-mec-phasics02','mec-visar-cage','mec-las-laptop03']
#         return list(zip(nmlist,qqip,qqn))
# =============================================================================
# =============================================================================
#     @classmethod
#     def cmp2_checker(cls):
#         cmpmsg=[];
#         #qqpl=[eacmp[0] for eacmp in MECcompylist()]
#         qqpl=cls.MECcompylist2()
#         with multiprocessing.Pool() as pool:
#             cmpmsg=pool.map(cls.plzchkcmp2,qqpl)
#         print('{:<15}'.format('Computer name')+'{:<28}'.format('IP shorthand')+'{:<6}'.format('Ping?'))
#         for ii in range(len(cmpmsg)):
#             print(''.join(cmpmsg[ii]))
#         return
# =============================================================================
#also consider adding ping, netconfig search, grep_pv, grep_ioc, serverStat, imgr, etc.
    @staticmethod
    def MECcompylist():
        qqn=['mec-laser','mec-monitor','mec-daq','mec-control','mec-console',
    'vitara','legend','evo1','evo2','gaia','topas',
    'PIMikroMove','spider','spectrometer','tundra',
    'lecroy1','lecroy2','lecroya','lecroyb',
    'dacage','visar1','visar2','rga','emp','phasicslaptop','phasics1','phasics2']
        nmlist=['mec-laser','mec-monitor','mec-daq','mec-control','mec-console',
    'mec-las-vitara','mec-las-laptop03','mec-las-laptop06','mec-las-laptop07','mec-las-laptop05','mec-las-laptop12',
    'mec-las-laptop09','mec-las-laptop11','mec-las-laptop01','win-ics-mec-tundra',
    'scope-ics-mectc1-1','scope-ics-meclas-lecroy01','scope-ics-meclas-lecroy-a','scope-ics-meclas-lecroy-b',
    'mec-visar-cage','win-ics-mec-visar1','win-ics-mec-visar2','mec-rga-laptop','scope-ics-mec-tektronix','mec-phasics-laptop01','win-ics-mec-phasics01','win-ics-mec-phasics02']
        return list(zip(nmlist,qqn))

    @classmethod
    def cmp_checker(cls):
        cmpmsg=[];
        #qqpl=[eacmp[0] for eacmp in MECcompylist()]
        qqpl=cls.MECcompylist()
        with multiprocessing.Pool() as pool:
            cmpmsg=pool.map(cls.plzchkcmp,qqpl)
        print('{:<15}'.format('Computer name')+'{:<28}'.format('IP shorthand')+'{:<6}'.format('Ping?'))
        for ii in range(len(cmpmsg)):
            print(''.join(cmpmsg[ii]))
        return




                  
class SCALLOPS:#port everything over from Bethany and Patin
    pass





class LabEnv:#(air/rack/etc. temp/humidity/etc.) 
    pass





class RIS:#etc.
    pass




                  
class PDU:
    pass





class GLOBAL:
    EYFE=EpicsSignal('MEC:LAS:FLOAT:01');
    ECD1w=EpicsSignal('MEC:LAS:FLOAT:02');
    EAB1w=EpicsSignal('MEC:LAS:FLOAT:03');
    EEF1w=EpicsSignal('MEC:LAS:FLOAT:04');
    EGH1w=EpicsSignal('MEC:LAS:FLOAT:05');
    EIJ1w=EpicsSignal('MEC:LAS:FLOAT:06');
    EAB2w=EpicsSignal('MEC:LAS:FLOAT:07');
    EEF2w=EpicsSignal('MEC:LAS:FLOAT:08');
    EGH2w=EpicsSignal('MEC:LAS:FLOAT:09');
    EIJ2w=EpicsSignal('MEC:LAS:FLOAT:10');
    CurrExp=EpicsSignal('MEC:LAS:FLOAT:11.DESC')
    CurrRun=EpicsSignal('MEC:LAS:FLOAT:11')
    #=EpicsSignal('MEC:LAS:FLOAT:12');
    #=EpicsSignal('MEC:LAS:FLOAT:13');
    TTLAB=EpicsSignal('MEC:LAS:FLOAT:14');
    TTLEF=EpicsSignal('MEC:LAS:FLOAT:15');
    TTLGH=EpicsSignal('MEC:LAS:FLOAT:16');
    TTLIJ=EpicsSignal('MEC:LAS:FLOAT:17');
    TTLWW=EpicsSignal('MEC:LAS:FLOAT:18');
    TTLXX=EpicsSignal('MEC:LAS:FLOAT:19');
    TTLREGEN=EpicsSignal('MEC:LAS:FLOAT:20');
    #
    EREGEN=EpicsSignal('MEC:LAS:FLOAT:21')
    ETOPAS=EpicsSignal('MEC:LAS:FLOAT:22')
    EMPA1=EpicsSignal('MEC:LAS:FLOAT:23')
    EMPA2=EpicsSignal('MEC:LAS:FLOAT:24')
    #=EpicsSignal('MEC:LAS:FLOAT:25')
    #=EpicsSignal('MEC:LAS:FLOAT:26')
    #=EpicsSignal('MEC:LAS:FLOAT:27')
    #=EpicsSignal('MEC:LAS:FLOAT:28')
    #=EpicsSignal('MEC:LAS:FLOAT:29')
    #=EpicsSignal('MEC:LAS:FLOAT:30')
    EcoeffYFE=EpicsSignal('MEC:LAS:FLOAT:31')
    Ecoeff1in1wCD=EpicsSignal('MEC:LAS:FLOAT:32')
    Ecoeff2in1wAB=EpicsSignal('MEC:LAS:FLOAT:33')
    Ecoeff2in1wEF=EpicsSignal('MEC:LAS:FLOAT:34')
    Ecoeff2in1wGH=EpicsSignal('MEC:LAS:FLOAT:35')
    Ecoeff2in1wIJ=EpicsSignal('MEC:LAS:FLOAT:36')
    Ecoeff2in2wAB=EpicsSignal('MEC:LAS:FLOAT:37')
    Ecoeff2in2wEF=EpicsSignal('MEC:LAS:FLOAT:38')
    Ecoeff2in2wGH=EpicsSignal('MEC:LAS:FLOAT:39')
    Ecoeff2in2wIJ=EpicsSignal('MEC:LAS:FLOAT:40')
    #(w/y/s1in1w/s42in1w/s42in2w/s + DateStr/today, RunNum, RunFilePath, PulseEnergies; notepadPVs: HAWG; YFE; 1w,2w,etc.; recipe  better than pickle?)
    #MEC:LAS:ARRAY:01, Desc: loaded pulse segment lengths, len: 10
    #MEC:LAS:ARRAY:02, Desc: loaded pulse segment endpoints, len: 1
    #MEC:LAS:ARRAY:03, Desc: Highland grafana, len: 1
    #MEC:LAS:ARRAY:04, Desc: YFE grafana, len: 1
    #MEC:LAS:ARRAY:05, Desc: YFEgoal grafana, len: 1
    #MEC:LAS:ARRAY:06, Desc: 1in1w grafana, len: 1
    #MEC:LAS:ARRAY:07, Desc: 2in1w grafana, len: 1
    #MEC:LAS:ARRAY:08, Desc: 2in2w grafana, len: 1
    #MEC:LAS:ARRAY:09, Desc: 2in2wgoal grafana, len: 1
    #MEC:LAS:ARRAY:10, Desc: Spare grafana, len: 1
    #MEC:LAS:ARRAY:11, Desc: Grafana, len: 1
    #MEC:LAS:ARRAY:12, Desc: Grafana, len: 1
    #MEC:LAS:ARRAY:13, Desc: Grafana, len: 1
    #MEC:LAS:ARRAY:14, Desc: Grafana, len: 1

    
    EcoeffRE1 = 1.64e5
    EcoeffRE0 = 1.03156061e-01 
    EcoeffTO1 = 3.48e7
    EcoeffTO1 = - 1.63e1 
    EcoeffM11 = 1.81e5
    EcoeffM10 = - 0.301
    EcoeffM21 = 1.05e5
    EcoeffM20 = - 1.39e-1 
    
    
    PSNS=EpicsSignal('MEC:LAS:ARRAY:01')
    SSS=EpicsSignal('MEC:LAS:ARRAY:02')
    WVFMHAWG=EpicsSignal('MEC:LAS:ARRAY:03')
    WVFMYFE=EpicsSignal('MEC:LAS:ARRAY:04')
    WVFMYFEGOAL=EpicsSignal('MEC:LAS:ARRAY:05')
    WVFM1IN1w=EpicsSignal('MEC:LAS:ARRAY:06')
    WVFM2IN1w=EpicsSignal('MEC:LAS:ARRAY:07')
    WVFM2IN2w=EpicsSignal('MEC:LAS:ARRAY:08')
    WVFM2IN2wGOAL=EpicsSignal('MEC:LAS:ARRAY:09')
    
    
    
    HWPAB=EpicsSignal(read_pv='MEC:NS1:MMS:02.RBV',write_pv='MEC:NS1:MMS:02.VAL');
    HWPEF=EpicsSignal(read_pv='MEC:NS1:MMS:01.RBV',write_pv='MEC:NS1:MMS:01.VAL');
    HWPGH=EpicsSignal(read_pv='MEC:LAS:MMN:30.RBV',write_pv='MEC:LAS:MMN:30.VAL');
    HWPIJ=EpicsSignal(read_pv='MEC:LAS:MMN:29.RBV',write_pv='MEC:LAS:MMN:29.VAL');

    HWPABoff=EpicsSignal('MEC:NS1:MMS:02.OFF');#dial offset
    HWPEFoff=EpicsSignal('MEC:NS1:MMS:01.OFF');#dial offset
    HWPGHoff=EpicsSignal('MEC:LAS:MMN:30.OFF');#dial offset
    HWPIJoff=EpicsSignal('MEC:LAS:MMN:29.OFF');#dial offset
    
    HWPABclr=EpicsSignal('MEC:NS1:MMS:02:SEQ_SELN');#clear start for mforce chassis
    HWPEFclr=EpicsSignal('MEC:NS1:MMS:01:SEQ_SELN');#clear start for mforce chassis
    
    
    EVRLPLLAMPEC=EpicsSignal(read_pv='EVR:MEC:USR01:TRIG6:EC_RBV',write_pv='EVR:MEC:USR01:TRIG6:TEC')#LPL lamp event code; needs 182
    EVRLPLLAMPEN=EpicsSignal('EVR:MEC:USR01:TRIG6:TCTL') #LPL lamp enable;
    EVRLPLSSEC=EpicsSignal(read_pv='EVR:MEC:USR01:TRIG7:EC_RBV', write_pv='EVR:MEC:USR01:TRIG7:TEC')#LPL slicer event code; needs 182 or 43 typically
    EVRLPLSSEN=EpicsSignal('EVR:MEC:USR01:TRIG7:TCTL') #LPL slicer enable;
    
    EGSPLRE=EpicsSignal('MEC:LAS:GENTEC:07:CH1:MEAS')
    EGSPLTO=EpicsSignal('MEC:LAS:GENTEC:07:CH2:MEAS')
    EGSPLM1=EpicsSignal('MEC:LAS:GENTEC:06:CH1:MEAS')
    EGSPLM2=EpicsSignal('MEC:LAS:GENTEC:06:CH2:MEAS')
    
    EGLPLYFE=EpicsSignal('MEC:LAS:LEM:03:A:CUR_DISP')
    EGLPL1in1w=EpicsSignal('MEC:LAS:LEM:03:B:CUR_DISP')
    
    EGLPL2in1wAB=EpicsSignal('MEC:LAS:GENTEC:02:CH1:MEAS')
    EGLPL2in1wEF=EpicsSignal('MEC:LAS:GENTEC:02:CH2:MEAS')
    EGLPL2in1wGH=EpicsSignal('MEC:LAS:GENTEC:01:CH1:MEAS')
    EGLPL2in1wIJ=EpicsSignal('MEC:LAS:GENTEC:01:CH2:MEAS')

    EGLPL2in2wAB=EpicsSignal('MEC:LAS:GENTEC:03:CH1:MEAS')
    EGLPL2in2wEF=EpicsSignal('MEC:LAS:GENTEC:03:CH2:MEAS')
    EGLPL2in2wGH=EpicsSignal('MEC:LAS:GENTEC:04:CH1:MEAS')
    EGLPL2in2wIJ=EpicsSignal('MEC:LAS:GENTEC:04:CH2:MEAS')
    
    EGLPLWest=EpicsSignal('MEC:GENTEC:01:CH2:MEAS');
    EGLPLEast=EpicsSignal('MEC:GENTEC:01:CH1:MEAS');

    MBCpwr=EpicsSignal('MEC:64B:PWR:2:Outlet:8:SetControlAction')#WAS 'MEC:S60:PWR:01:Outlet:7:SetControlAction'#read AND write:1=ON,2=OFF
    MBCmode=EpicsSignal('MEC:LPL:MBC:01:RunningMode_RBV',write_pv='MEC:LPL:MBC:01:RunningMode')#AUTO=0,MAN=1
    MBCsetpt=EpicsSignal('MEC:LPL:MBC:01:AutoCalibration.VAL',write_pv='MEC:LPL:MBC:01:AutoCalibration') #QUAD=0,MIN=1,MAX=2
    MBCbias=EpicsSignal('MEC:LPL:MBC:01:BiasValue_RBV',write_pv='MEC:LPL:MBC:01:BiasValue')
    MBCfault=EpicsSignal('MEC:LPL:MBC:01:ErrorStatus',write_pv='MEC:LPL:MBC:01:ClearErrors')
    
    SHGABmot=EpicsSignal(read_pv='MEC:LAS:MMN:22.RBV',write_pv='MEC:LAS:MMN:22.VAL')
    SHGEFmot=EpicsSignal(read_pv='MEC:LAS:MMN:24.RBV',write_pv='MEC:LAS:MMN:24.VAL')
    SHGGHmot=EpicsSignal(read_pv='MEC:LAS:MMN:17.RBV',write_pv='MEC:LAS:MMN:17.VAL')
    SHGIJmot=EpicsSignal(read_pv='MEC:LAS:MMN:18.RBV',write_pv='MEC:LAS:MMN:18.VAL')
    
    
    PFNSS=EpicsSignal('MEC:PFN:SINCESHOT')
    PFNmode=EpicsSignal('MEC:PFN:MODE')
    PFNCDEN=EpicsSignal(read_pv='MEC:PFN:CH0:ENABLE_RBV',write_pv='MEC:PFN:CH0:ENABLE');
    PFNCDCS=EpicsSignal('MEC:PFN:CH0:CHARGE_STATE');
    PFNAEN=EpicsSignal(read_pv='MEC:PFN:CH1:ENABLE_RBV',write_pv='MEC:PFN:CH1:ENABLE');
    PFNACS=EpicsSignal('MEC:PFN:CH1:CHARGE_STATE');
    PFNBEN=EpicsSignal(read_pv='MEC:PFN:CH2:ENABLE_RBV',write_pv='MEC:PFN:CH2:ENABLE');
    PFNBCS=EpicsSignal('MEC:PFN:CH2:CHARGE_STATE');
    PFNEEN=EpicsSignal(read_pv='MEC:PFN:CH3:ENABLE_RBV',write_pv='MEC:PFN:CH3:ENABLE');
    PFNECS=EpicsSignal('MEC:PFN:CH3:CHARGE_STATE');
    PFNFEN=EpicsSignal(read_pv='MEC:PFN:CH4:ENABLE_RBV',write_pv='MEC:PFN:CH4:ENABLE');
    PFNFCS=EpicsSignal('MEC:PFN:CH4:CHARGE_STATE');
    PFNGEN=EpicsSignal(read_pv='MEC:PFN:CH5:ENABLE_RBV',write_pv='MEC:PFN:CH5:ENABLE');
    PFNGCS=EpicsSignal('MEC:PFN:CH5:CHARGE_STATE');
    PFNHEN=EpicsSignal(read_pv='MEC:PFN:CH6:ENABLE_RBV',write_pv='MEC:PFN:CH6:ENABLE');
    PFNHCS=EpicsSignal('MEC:PFN:CH6:CHARGE_STATE');
    PFNIEN=EpicsSignal(read_pv='MEC:PFN:CH7:ENABLE_RBV',write_pv='MEC:PFN:CH7:ENABLE');
    PFNICS=EpicsSignal('MEC:PFN:CH7:CHARGE_STATE');
    PFNJEN=EpicsSignal(read_pv='MEC:PFN:CH8:ENABLE_RBV',write_pv='MEC:PFN:CH8:ENABLE');
    PFNJCS=EpicsSignal('MEC:PFN:CH8:CHARGE_STATE');

    LPLPCpwr=EpicsSignal('MEC:S60:PWR:01:Outlet:6:SetControlAction')
    LPLVACpwr=EpicsSignal('MEC:S60:PWR:01:Outlet:7:SetControlAction')
    LPLPS1pwr=EpicsSignal('MEC:S60:PWR:01:Outlet:1:SetControlAction')
    LPLHAWGpwr=EpicsSignal('MEC:64B:PWR:2:Outlet:1:SetControlAction')

    XPS3pwr=EpicsSignal('MEC:64A:PWR:2:Outlet:5:SetControlAction')
    XPS4pwr=EpicsSignal('MEC:64B:PWR:1:Outlet:1:SetControlAction')
    
    XPS1IALL=EpicsSignal('MEC:LAS:MMN_0108.IALL')
    XPS1RALL=EpicsSignal('MEC:LAS:MMN_0108.RALL')
    XPS2IALL=EpicsSignal('MEC:LAS:MMN_0916.IALL')
    XPS2RALL=EpicsSignal('MEC:LAS:MMN_0916.RALL')
    XPS3IALL=EpicsSignal('MEC:LAS:MMN_1724.IALL')
    XPS3RALL=EpicsSignal('MEC:LAS:MMN_1724.RALL')
    XPS4IALL=EpicsSignal('MEC:LAS:MMN_2532.IALL')
    XPS4RALL=EpicsSignal('MEC:LAS:MMN_2532.RALL')

    OKHOSTS=['mec-monitor', 'mec-daq', 'mec-laser']
    OKUSERS=['mecopr']

    LMapAB=[5,100]
    LMap2=[50,1000]
    PSFILEPATH='/reg/neh/operator/mecopr/mecpython/pulseshaping/'
    
    @classmethod
    def notepadPVreset(cls):
        efc.wPV('MEC:LAS:FLOAT:01.DESC', 'E_synth_YFE');
        efc.wPV('MEC:LAS:FLOAT:02.DESC', 'E_synth_CD1w');
        efc.wPV('MEC:LAS:FLOAT:03.DESC', 'E_synth_AB1w');
        efc.wPV('MEC:LAS:FLOAT:04.DESC', 'E_synth_EF1w');
        efc.wPV('MEC:LAS:FLOAT:05.DESC', 'E_synth_GH1w');
        efc.wPV('MEC:LAS:FLOAT:06.DESC', 'E_synth_IJ1w');
        efc.wPV('MEC:LAS:FLOAT:07.DESC', 'E_synth_AB2w');
        efc.wPV('MEC:LAS:FLOAT:08.DESC', 'E_synth_EF2w');
        efc.wPV('MEC:LAS:FLOAT:09.DESC', 'E_synth_GH2w');
        efc.wPV('MEC:LAS:FLOAT:10.DESC', 'E_synth_IJ2w');
        #11: DESC is CurrExp
        efc.wPV('MEC:LAS:FLOAT:12.DESC', 'reserved');
        efc.wPV('MEC:LAS:FLOAT:13.DESC', 'reserved');
        efc.wPV('MEC:LAS:FLOAT:14.DESC', 'AB shutter state');
        efc.wPV('MEC:LAS:FLOAT:15.DESC', 'EF shutter state');
        efc.wPV('MEC:LAS:FLOAT:16.DESC', 'GH shutter state');
        efc.wPV('MEC:LAS:FLOAT:17.DESC', 'IJ shutter state');
        efc.wPV('MEC:LAS:FLOAT:18.DESC', 'WEST (ABEF) shutter state');
        efc.wPV('MEC:LAS:FLOAT:19.DESC', 'EAST (GHIJ)shutter state');
        efc.wPV('MEC:LAS:FLOAT:20.DESC', 'Regen shutter state');
        #
        efc.wPV('MEC:LAS:FLOAT:21.DESC', 'E_synth_regen');
        efc.wPV('MEC:LAS:FLOAT:22.DESC', 'E_synth_TOPAS');
        efc.wPV('MEC:LAS:FLOAT:23.DESC', 'E_synth_MPA1');
        efc.wPV('MEC:LAS:FLOAT:24.DESC', 'E_synth_MPA2');
        efc.wPV('MEC:LAS:FLOAT:25.DESC', 'reserved');
        efc.wPV('MEC:LAS:FLOAT:26.DESC', 'reserved');
        efc.wPV('MEC:LAS:FLOAT:27.DESC', 'reserved');
        efc.wPV('MEC:LAS:FLOAT:28.DESC', 'reserved');
        efc.wPV('MEC:LAS:FLOAT:29.DESC', 'reserved');
        efc.wPV('MEC:LAS:FLOAT:30.DESC', 'reserved');
        efc.wPV('MEC:LAS:FLOAT:31.DESC', 'E_coeff_YFE');
        efc.wPV('MEC:LAS:FLOAT:32.DESC', 'E_coeff_CD1w');
        efc.wPV('MEC:LAS:FLOAT:33.DESC', 'E_coeff_AB1w');
        efc.wPV('MEC:LAS:FLOAT:34.DESC', 'E_coeff_EF1w');
        efc.wPV('MEC:LAS:FLOAT:35.DESC', 'E_coeff_GH1w');
        efc.wPV('MEC:LAS:FLOAT:36.DESC', 'E_coeff_IJ1w');
        efc.wPV('MEC:LAS:FLOAT:37.DESC', 'E_coeff_AB2w');
        efc.wPV('MEC:LAS:FLOAT:38.DESC', 'E_coeff_EF2w');
        efc.wPV('MEC:LAS:FLOAT:39.DESC', 'E_coeff_GH2w');
        efc.wPV('MEC:LAS:FLOAT:40.DESC', 'E_coeff_IJ2w');
        #last updated 20220128
        efc.wPV('MEC:LAS:FLOAT:31', 0.3285);#was .3578
        efc.wPV('MEC:LAS:FLOAT:32', 0.5871);
        efc.wPV('MEC:LAS:FLOAT:33', 224.0);
        efc.wPV('MEC:LAS:FLOAT:34', 177.5);
        efc.wPV('MEC:LAS:FLOAT:35', 260.9826);
        efc.wPV('MEC:LAS:FLOAT:36', 113.2);
        efc.wPV('MEC:LAS:FLOAT:37', 134.0135);
        efc.wPV('MEC:LAS:FLOAT:38', 165.2398);
        efc.wPV('MEC:LAS:FLOAT:39', 194.1412);
        efc.wPV('MEC:LAS:FLOAT:40', 156.9307);
        
        efc.wPV('MEC:LAS:ARRAY:01', 'Psns pulse segment lengths:10')
        efc.wPV('MEC:LAS:ARRAY:02', 'SSs pulse segment endpoint pairs:20')
        efc.wPV('MEC:LAS:ARRAY:03', 'Highland Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:04', 'YFE Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:05', 'YFEgoal Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:06', '1in1w Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:07', '2in1w Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:08', '2in2w Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:09', '2in2wgoal Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:10', 'Spare Grafana')
        efc.wPV('MEC:LAS:ARRAY:11', 'Spare Grafana')
        efc.wPV('MEC:LAS:ARRAY:12', 'Spare Grafana')
        efc.wPV('MEC:LAS:ARRAY:13', 'Spare Grafana')
        efc.wPV('MEC:LAS:ARRAY:14', 'Spare Grafana')
        





