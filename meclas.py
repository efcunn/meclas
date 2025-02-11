"""
             meclas.py
Package for running various laser utilities in MEC
Apologies on behalf of: Eric Cunningham (and others)

To load: use import meclas or use IPython's %run magic function

Class list and brief description:
     LPL -- routines for LPL pulse shaping (with some aux functions),
            data acquisition, etc.
     efc -- extra function class, holds many useful utilities and shortcuts
     ep -- shortcut plotting functions for the lazy
     HAWG -- Highland AWG control and readout
     LOSC -- LeCroy oscilloscope trace read-out, plotting, saving, etc.
     EMeters -- LPL and SPL energy meters
     MBC -- LPL bias controller utilities
     YFE -- LPL YLF Front End seed laser utilities
     PFN -- LPL Pulse Forming Network cap bank charging utilities
     HWP -- LPL Half Wave Plate motor utilities
   **Stage -- Newport and SmarAct stage utilities
   **Timing -- ns and fs timing utilities
     CAM -- functions for GigE camera acquisition, configuration, etc.
     TTL_shutter -- Beckhoff utilities for controlling/tracking Thorlabs TTL
                    shutters
   **DG645 -- functions for DG645 operation, parameter backup and restoration,
              etc.
   **SPL -- routines for SPL alignment, etc.
   **UNIBLITZ -- UNIBLITZ shutter utilities for setting SPL trigger modes, etc.
   **Spectrometer -- functions for Qmini and Ocean Optics USB4000 spectrometers
   **VISAR -- routines for VISAR timing, streak camera configuration,
              laser control, etc.
     CtrlSys -- routines for checking responsivity of PVs, hosts,
                hutch computers, etc.
   **SCALLOPS -- routines for LPL pulse shaping simulations
   **LabEnv -- functions for interfacing with lab environment monitors
   **RIS -- functions for monitoring RIS-related systems and PVs
   **PDU -- functions for operating power distribution units
     GLOBAL -- home for global constants, PV definitions, etc.

   ** = FUTURE DEVELOPMENT
"""

# load packages
import socket
import time
import math
import numpy as np
import struct
# from scipy import signal
# from scipy import stats
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
# import pandas as pd
import stat
import getpass
import multiprocessing
# import threading
import termios
import tty
import select
import re
# import regex as re



class LPL:
    """
    Stores functions related to LPL pulse shaping, data acquisition, etc.
    Functions include:
        :_LinearWave #function for linear waveforms appropriate for loading to
            the Highland AWG
        :_LinearWave2 #function for linear waveforms of variable length,
            appropriate for specifying targeted pulse shapes
        :_ParabolicWave2 #function for parabolic waveforms of variable length,
            appropriate for specifying targeted pulse shapes
        :_ExponentialWave2 #function for exponential waveforms of variable
            length, appropriate for specifying targeted pulse shapes
        :_EW2 #shorthand version of _ExponentialWave2 tailored for facilitating
            LPL shaping at 10Hz
        :_EW2stringhint #same as above but produces strings for the benefit of
            saving wave hints in recipes
        :_ComboWave #combines waveforms for pulse shaping goals
        :_TraceFormatting #projects/downsamples waveforms from one horizontal
            base to another
        :_UpdatingShapingAlgorithm #calculates new Highland input based on old
            input, its corresponding output, and a target
        :_FixEdges #allows tweaking of behavior of points near waveform
            discontinuities (edges and steps)
        :_SmoothWvfm #locally smooths an input waveform as an option as part of
            pulse shaping
        :_PulseGoal #helps define a targeted waveform for the full energy
            output
        :_PulseMax #helps set the desired amplitude of _PulseGoal based on
            shape and desired output energy
        :_Psns_get #retrieves the list of pulse segment durations of the
            current targeted output pulse
        :_Psns_set #sets the list of pulse segment durations of a new targeted
            output pulse
        :_SSs_get #retrieves the list of pulse segment start/stop heights of
            the current targeted output pulse
        :_SSs_set #sets the list of pulse segment start/stop heights of a new
            targeted output pulse
        :_YSSs_get #retrieves the list of YFE exponential pulse segment
            start/stop heights of the current targeted 10Hz output pulse
        :_YSSs_set #sets the list of YFE expnential pulse segment start/stop
            heights of a new targeted 10Hz output pulse
        :_wIter2 #wrapper function for using _UpdatingShapingAlgorithm in
            context
        :_weichall #generates weighted waveforms for YFE, 1in1w, 4x2in1w, and
            4x2in2w outputs using scopes and energy meters
        :_weichToPowerVsTime #converts energy-weighted waveforms into a tuple
            of instantaneous power vs time
        :_PGToPowerVsTime #converts scaled _PulseGoal waveforms into a tuple of
            instantaneous power vs time
        :_pshostcheck #checks the host of the current computer to verify it is
            a machine with all the proper network access for pulse shaping
        :_DateString #shortcut for generating a string of today's date
        :get_curr_exp #retrieves the current experiment name in MEC
        :get_curr_run #retrieves the current run number in MEC
        :get_curr_shape #retrieves the last loaded pulse shape in MEC
        :_psheaders #prepares the groundwork for pulse shaping exercises
        :_psacqx #acquires data after taking a shot
        :_psefc #plots most-recent shot compared to its goal and returns a
            suggestion for a new input waveform
        :psefc10Hz #performs pulse shaping at 10Hz to converge towards an input
            goal for the YFE output
        :_psupd #shortcut for updating the Highland waveform
        :psloadwvfm #loads a previously-saved pulse shaping recipe
        :pssavewvfm #saves a new pulse shaping recipe
        :psviewwvfm #displays a previously-saved pulse shaping recipe
        :psrefrwvfm #refreshes a previously-saved pulse shaping recipe to
            account for system drift
        :psrecipes #lists all previously-saved pulse shaping recipes
        :psmenu #allows loading or viewing of previously-saved pulse shaping
            recipes from a list
        :pspreshot #executes routine that prepares the state of the laser and
            its diagnostics for a full-energy shot
        :pspostshot #executes routine that records all data after a full-energy
            shot and returns the laser to a stand-by state
        :On #turns the long-pulse laser system on
        :SHG_opt #executes the optimization routine for the SHG crystal angles
            of all four arms of the LPL

    Potential future work includes:
       - put all the stuff Galtier wants to do inside lpl
           - LPL.On() instead of YFE.On(), e.g.
           - same for seeing shapes, changing energy, etc.
           - help avoid TAB-complete problems because of object instantiation
               - like LOSC('a'), for example
       - change names to be lowercases so Galtier can TAB-complete easier
       - add underscores to methods not generally intended for non-expert usage
       - consider consolidating some functions into one
           - e.g. _LinearWave, _ExponentialWave, etc.
       - add Deconvolution function -- help shape, seems that 100ps length also
         affected
           - account for the instrument response of the PD, amp, scope, etc. in
             determining a detected waveform
       - save PFN voltages on-shot too?
       - Scope vertical resolution problem -- casting somewhere? better to
           motorize characterized ND filters...
       - Slow feedback to 10Hz pulse goal based on full shot rbv? Esp to help
           front edge problems?
           - Need this combined with SCALLOPS? Need a fit function to data for
               the full rbv?
       - YSSs: YFE equivalent of SSs
           - Steal a notepad pv array
       - _EW2(Psns, YSSs)
       - Have a version of pulsegoal but 10hz using _ew2 to interpolate instead
           of lw, still have dellist
       - SmartPulseGoal10Hz using output of SCALLOPS
       """

    def _LinearWave(Edge1PixNo, Edge1Height, Edge2PixNo, Edge2Height):
        """
        Generates hex string-formatted linearly-interpolated waveform of
            length 140 between two input points
        Primarily intended for use with the Highland AWG, so max value capped
            at 65535

        _LinearWave(1,10000,120,28000) returns a string of length 140
            (560 chars) with linear ramp from (pixel 1, height 10000) to
            (pixel 120, height 28000)
        """
        itt = 0
        NewString = ''
        if Edge1Height > 65535:
            print('Edge1 height exceeds max value of 65535')
            h1 = 65535
        elif Edge1Height < 0:
            print('Edge1 height must be positive')
            h1 = 0
        else:
            h1 = int(Edge1Height)
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

    def _LinearWave2(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height,offsetQ=0,arraylenQ=5002):
        """
        Generates linearly-interpolated waveform between two points (and 0 outside those points)
        with requested linear offset and array length
        
        Useful in part for specifying YFE waveforms at 10Hz, generating "goal" waveforms, etc.
        
        _LinearWave2(500,.1,1025,.8,0,5002) returns an array of length 5002 with a linear ramp
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

    def _ParabolicWave2(Edge1PixNo,Edge1Height,MidPixNo,MidHeight,Edge2PixNo,Edge2Height,offsetQ=0,arraylenQ=5002):
        """
        Generates parabolically-interpolated waveform using three points (and 0 outside those points)
        with requested linear offset and array length
        
        Only rarely used for specifying YFE waveforms at 10Hz, generating "goal" waveforms, etc.
        
        _ParabolicWave2(500,.1,800,.15,1025,.8,0,5002) returns an array of length 5002 with a
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

    def _ExponentialWave2(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height,offsetQ=0,arraylenQ=5002):
        """
        Generates exponentially-interpolated waveform using two points (and 0 outside those points)
        with requested linear offset and array length
        
        Most-used function for specifying YFE waveforms at 10Hz (exponential seed waveforms
            become ~linear after laser amplification)

        _ExponentialWave2(500,.1,1025,.8,0,5002) returns an array of length 5002 with an exponential ramp
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


    @classmethod
    def _EW2(cls,Psns='curr',SSs='curr',YSSs='curr',offsetQ=0,YFEbkgrdY=-.004):
        """
        Shorthand for generating most common _ExponentialWave2 output (or even combination of several _ExponentialWave2 outputs)
        
        Preferred use is to use LPL._Psns_set(), LPL._SSs_set(), and LPL._YSSs_set() for _EW2 to find and use
        
        When using Psns='curr' | SSs='curr' | YSSs='curr': Psns, SSs, and YSSs all are loaded with the corresponding ***_get() commands
        
        Alternatively they can be specified explicitly, i.e. Psns=[10.25], SSs=[[98,100]], YSSs=[[.02,.114]]
        
        Compared to _ExponentialWave2, _EW2:
            :can combine several summed _ExponentialWave2 outputs into one convenient function
                :Example: _EW2(Psns=[5,5.25],SSs=[[49,50],[98,100]],YSSs=[[.01,.03],[.055,.16]]) is equivalent to
                 _ExponentialWave2(500,.01,1000,.03,0,5002)+_ExponentialWave2(1000,.055,1525,.16,0,5002)
            :infers Edge1PixNo and Edge2PixNo of each _ExponentialWave2 segment from Psns and SSs
            :specifies Edge1Height and Edge2Height of each _ExponentialWave2 segment using YSSs
            :still permits offset specification using offsetQ, if desired
            :automatically uses the standard arraylenQ of 5002
        """
        if Psns=='curr':
            Psns=cls._Psns_get()
        if SSs=='curr':
            SSs=cls._SSs_get()
        if YSSs=='curr':
            YSSs=cls._YSSs_get()
        if len(Psns) != len(YSSs):
            print('Lengths of Psns ({}) and YSSs ({}) do not match! Exiting...'.format(len(Psns),len(YSSs)))
            return False
        outwvfm=cls._LinearWave2(500,0,1025,0,0,5002);#500, 1000 1000, 1525
        YPsns=np.cumsum([0]+[Psns[ii]-0.25*(1 if SSs[ii][1] == SSs[ii+1][0] else 0) for ii in range(len(Psns)-1)] + [Psns[-1]])
        YPsnsPairs=[[int(500+YPsns[ii]*100), int(500+YPsns[ii+1]*100)] for ii in range(len(YPsns)-1)]
        try:
            for ii in range(len(YSSs)):
                outwvfm+=cls._ExponentialWave2(Edge1PixNo=YPsnsPairs[ii][0],Edge1Height=YSSs[ii][0],
                                          Edge2PixNo=YPsnsPairs[ii][1],Edge2Height=YSSs[ii][1],offsetQ=offsetQ,arraylenQ=5002)
            return outwvfm
        except:
            print('Failed to generate waveform!')
            return False

    @classmethod        
    def _EW2stringhint(cls,Psns='curr',SSs='curr',YSSs='curr',YFEbkgrdY=-.004):
        """
        Shorthand for generating string hint based most common _ExponentialWave2 output (or even combination of several _ExponentialWave2 outputs)
        
        Preferred use is to use LPL._Psns_set(), LPL._SSs_set(), and LPL._YSSs_set() for _EW2 to find and use
        
        When using Psns='curr' | SSs='curr' | YSSs='curr': Psns, SSs, and YSSs all are loaded with the corresponding ***_get() commands
        
        Alternatively they can be specified explicitly, i.e. Psns=[10.25], SSs=[[98,100]], YSSs=[[.02,.114]]
        
        Compared to _ExponentialWave2, _EW2:
            :can combine several summed _ExponentialWave2 outputs into one convenient function
                :Example: _EW2(Psns=[5,5.25],SSs=[[49,50],[98,100]],YSSs=[[.01,.03],[.055,.16]]) is equivalent to
                 _ExponentialWave2(500,.01,1000,.03,0,5002)+_ExponentialWave2(1000,.055,1525,.16,0,5002)
            :infers Edge1PixNo and Edge2PixNo of each _ExponentialWave2 segment from Psns and SSs
            :specifies Edge1Height and Edge2Height of each _ExponentialWave2 segment using YSSs
            :still permits offset specification using YFEbkgrdY, if desired
            :automatically uses the standard arraylenQ of 5002 and EW offset of 0
        """
        if Psns=='curr':
            Psns=cls._Psns_get()
        if SSs=='curr':
            SSs=cls._SSs_get()
        if YSSs=='curr':
            YSSs=cls._YSSs_get()
        if len(Psns) != len(YSSs):
            print('Lengths of Psns ({}) and YSSs ({}) do not match! Exiting...'.format(len(Psns),len(YSSs)))
            return False
        outwvfm=''
        YPsns=np.cumsum([0]+[Psns[ii]-0.25*(1 if SSs[ii][1] == SSs[ii+1][0] else 0) for ii in range(len(Psns)-1)] + [Psns[-1]])
        YPsnsPairs=[[int(500+YPsns[ii]*100), int(500+YPsns[ii+1]*100)] for ii in range(len(YPsns)-1)]
        try:
            for ii in range(len(YSSs)):
                outwvfm+='_ExponentialWave2({},{},{},{},0,5002)+'.format(YPsnsPairs[ii][0],YSSs[ii][0],YPsnsPairs[ii][1],YSSs[ii][1])
            outwvfm=outwvfm[:-1]
            outwvfm+=';; YFEbkgrdY={}'.format(YFEbkgrdY)
            return outwvfm
        except:
            print('Failed to generate string!')
            return False

    def _ComboWave(WList): #accept list or csv of 140 pts
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

    def _TraceFormatting(PDTrace, PDFETMap, MaxPDValue, AvgRange=25, FWHM=4):
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

    def _UpdatingShapingAlgorithm(DesiredOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape, StepQ):
        """
        Accepts pre-formatted input, measurement, and goal waveforms to calculate next-iteration input waveform using specified step size
        """
        G, M, I = DesiredOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape
        NewInputPulseShape=np.clip([abs((StepQ*(G[ii]-M[ii]))+I[ii])*np.abs(np.sign(G[ii])) for ii in range(len(G))],0,1)#np.abs(np.sign(G[ii])) is a mask that disallows values outside the goal
        return NewInputPulseShape

    def _FixEdges(WavF,DurationListQ,StartStopListQ,PtNumFront=3,PtNumBack=2,CorrFactorFront=.97,CorrFactorBack=1.):
        """
        Applies fixed relationship between points close to waveform discontinuities (e.g. at beginning, end, step, etc.)
        specified by Duration and Start/Stop lists
        
        WavF is the input waveform 
        DurationListQ is the pulse segment durations (e.g. Psns, like [10.25] or [5,5.25])
        StartStopListQ is the pulse segment stop/stop list (e.g. SSs, like [[98,100]] or [[24,25],[98,100]])
        PtNumFront is the number of points at the front edge of the pulse to be fixed
        PtNumBack is the number of points at the back edge of the pulse to be fixed
        CorrFactorFront is the fixed multiplicative factor applied from point to point moving towards the front edge of the pulse
           :for the example of PtNumFront=2, CorrFactorFront=0.97 sets P_2=0.97*P_3 then P1=0.97*P_2
        CorrFactorBack is the fixed multiplicative factor applied from point to point moving towards the back edge of the pulse
           :for the example of PtNumBack=3, CorrFactorBack=1.01 sets P_{K-2}=1.01*P_{K-3} then P_{K-1}=1.01*P_{K-2} 
            then P_{K}=1.01*P_{K-1}, where P_{K} is the pulse's last point, P_{K-1} is the pulse's second-to-last point, etc.
        
        New waveform with fixed edges is returned
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
                if fWavF[ii]>GLOBAL.PDMAX_TEST:
                    fWavF[ii]=GLOBAL.PDMAX_TEST
        return fWavF
    
    def _SmoothWvfm(wvfm_in):
        """
        Performs rudimentary smoothing of waveform by looking at neighboring pixels
        
        Accepts single waveform as input and return smoothed output waveform
        """
        wvfm_out=wvfm_in[:]
        for ii in range(len(wvfm_in)-2):
            wvfm_out[ii+1]=.25*wvfm_in[ii]+.5*wvfm_in[ii+1]+.25*wvfm_in[ii+2]
        return wvfm_out

    @classmethod
    def _PulseGoal(cls,DurationListQ,StartStopListQ):#140 pt list with max at 1
        """
        Generates 140-pt list according to provided Duration and Start/Stop lists
        """
        BeginPix=1
        DurListQ=np.cumsum([0]+list(DurationListQ))
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
            SegmentsQ.append(cls._LinearWave(int(BeginPix+(DurListQ[ii]*4)),int(20000.*SSListQ[ii][0]/100.),int(BeginPix+(DurListQ[ii+1]*4)-1),int(20000.*SSListQ[ii][1]/100.)))
        return np.append(np.delete(np.array(cls._ComboWave(SegmentsQ)),np.array(DelListQ).astype(int)),[0]*len(DelListQ))

    @classmethod
    def _PulseMax(cls,DurationListQ,StartStopListQ,zzJQ):
        """Gets amplitude setting for segmented, arbitrary _PulseGoal using Duration and Start/Stop lists and targeted energy"""
        return (1.0*StartStopListQ[-1][-1]/100.)*(50.0*zzJQ/(5.0*500.0*np.sum(cls._PulseGoal(DurationListQ,StartStopListQ))))
    
    def _Psns_get():
        """
        Return the current value of Psns, which is an array of pulse duration segments in nanoseconds
        Example: after loading a 15ns flat-top pulse, LPL._Psns_get() will return [15.25]
        Example: after loading a 5ns-5ns step pulse, LPL._Psns_get() will return [5,5.25]
        """
        return list(GLOBAL.PSNS.get())
        
    def _Psns_set(NewPsnsArr):
        """
        Sets the current value of Psns, which is an array of pulse duration segments in nanoseconds
        Each individual value in the array must be a multiple of 0.25 (spacing of Highland AWG pixels is 250ps or 0.25ns)
        Psns must have one segment duration for every pair of segment heights in the segment heights parameter SSs
        Total summed length of array values should be less than 35ns (total shaping window of Highland AWG)
        If the back end of a particular segment is continuous in height with the front end of the following segment
            (e.g. for smooth ramp pulses), you must add 0.25 to the duration of the first segment
            
        This function can be used in preparation for making a new pulse recipe
        Example: use LPL._Psns_set([7.25]) to prepare to make a new 7ns pulse
        Example: use LPL._Psns_set([5,8.25]) to prepare to make a new 5ns-8ns step pulse
        """
        if type(NewPsnsArr) is not list:
            NewPsnsArr=[NewPsnsArr]
        for val in NewPsnsArr:
            if val%.25 != 0:
                print('Psns values need to be multiples of 0.25ns! Update failed!')
                return False
        try:
            GLOBAL.PSNS.set(NewPsnsArr)
            return
        except:
            print('Psns update failed!')
            return False
        
    @classmethod    
    def _SSs_get(cls):
        """
        Return the current value of SSs, which is an array of pairs of pulse segments heights as percentages of maximum (100%)
        Example: after loading a 15ns flat-top pulse (Psns=[15.25]), LPL._SSs_get() may return something like [[98,100]]
        Example: after loading a 3ns-7ns pulse with a 4x ratio (Psns=[3,7.25]), LPL._SSs_get() may return something like [[24,25],[98,100]]
        """
        return [list(GLOBAL.SSS.get()[2*ii:2*ii+2]) for ii in range(len(cls._Psns_get()))]
        
    def _SSs_set(NewSSsArr):
        """
        Sets the current value of SSs, which is an array of pairs of pulse segments heights as percentages of maximum (100%)
        SSs must have one pair of heights for every duration segment in the segment durations parameter Psns
        This function can be used in preparation for making a new pulse recipe
        Example: use LPL._SSs_set([[98,100]]) to prepare to make a new flat-top
        Example: use LPL._Psns_set([[32,33],[75,100]]) to prepare to make a new step pulse with a flat first step and a
                 15%-gradient second step and a 3x ratio between the maxima of the two steps
        """
        if len(np.array(NewSSsArr).shape) == 1:
            try:
                if len(NewSSsArr)%2 == 0:
                    GLOBAL.SSS.set(NewSSsArr)
                    return
                else:
                    raise Exception
            except:
                print('Failed to write new SSs values! SSs length: {}'.format(len(NewSSsArr)))
                return False
        elif len(np.array(NewSSsArr).shape) == 2:
            try:
                GLOBAL.SSS.set([val for SSpair in NewSSsArr for val in SSpair])
                return
            except:
                print('Parsing failure! Failed to write new SSs values!')
                return False
        else:
            print('Unexpected shape! Failed to write new SSs values!')
            return False
        
    @classmethod    
    def _YSSs_get(cls):
        """
        Return the current value of YSSs, which is an array of pairs of pulse segments heights as YFE diode amplitude
        These parameters are used as shorthand in _EW2()
        """
        return [list(GLOBAL.YSSS.get()[2*ii:2*ii+2]) for ii in range(len(cls._Psns_get()))]
        
    def _YSSs_set(NewYSSsArr):
        """
        Sets the current value of YSSs, which is an array of pairs of pulse segments heights as YFE diode amplitude
        YSSs must have one pair of heights for every duration segment in the segment durations parameter Psns
        This function can be used in specifying the YFE goal segment heights
        """
        if len(np.array(NewYSSsArr).shape) == 1:
            try:
                if len(NewYSSsArr)%2 == 0:
                    GLOBAL.YSSS.set(NewYSSsArr)
                    return
                else:
                    raise Exception
            except:
                print('Failed to write new YSSs values! YSSs length: {}'.format(len(NewYSSsArr)))
                return False
        elif len(np.array(NewYSSsArr).shape) == 2:
            try:
                GLOBAL.YSSS.set([val for YSSpair in NewYSSsArr for val in YSSpair])
                return
            except:
                print('Parsing failure! Failed to write new YSSs values!')
                return False
        else:
            print('Unexpected shape! Failed to write new YSSs values!')
            return False

    @classmethod
    def _wIter2(cls,sQ,wQ,DurationListQ,StartStopListQ,zzJQ,mapnowQ,stepQQ,Hdisplay=False):
        """Calculates next suggested AWG input given 1) a previous full-energy waveform (+ mapping) and its corresponding AWG input,
        2) the Duration and Start/Stop lists to specify the goal, and 3) the requested step size of the correction"""
        avgfwhm=90;avgrange=11;#250;
        DurListQ=np.cumsum([0]+DurationListQ)
        w1,w2=0,int(DurListQ[-1]*4)+5 # 50-5, 50+int(DurListQ[-1]*4)+5
        PGQ=cls._PulseGoal(DurationListQ,StartStopListQ)
        if np.abs(len(sQ)-10000)<10:
            PMQcorr=1
        elif np.abs(len(sQ)-1000)<10:
            PMQcorr=10
        else:
            print('Warning: unanticipated pulse shape array length of '+str(len(sQ))+', _PulseMax scaling may be off...')
            PMQcorr=1
        PMQ=cls._PulseMax(DurationListQ,StartStopListQ,zzJQ)*PMQcorr
        wnew2=cls._FixEdges(cls._UpdatingShapingAlgorithm(PGQ,cls._TraceFormatting(sQ,mapnowQ,PMQ,AvgRange=avgrange,FWHM=avgfwhm), wQ,stepQQ),DurationListQ,StartStopListQ)
        if Hdisplay == True:
            ep.ll([0.*np.array(wnew2[w1:w2]),np.array(wnew2[w1:w2])-np.array(wQ[w1:w2]),np.array(cls._TraceFormatting(sQ,mapnowQ,PMQ,AvgRange=avgrange,FWHM=avgfwhm))[w1:w2]*.6,np.array(PGQ)[w1:w2]*.6])
        ep.llxy([cls._weichToPowerVsTime(sQ),cls._PGToPowerVsTime(Psns=DurationListQ, SSs=StartStopListQ, zzJQ=zzJQ)],
                   xlb='Time (ns)',ylb='Power (W)',
                   xlim=[-1,1+np.sum(DurationListQ)-0.25*np.sum([1 if StartStopListQ[ii][1] == StartStopListQ[ii][0] else 0 for ii in range(len(StartStopListQ)-1)])])
        return wnew2

    def _weichall():
        """
        Generates weighted waveforms for YFE, 1in1w, 4x2in1w, and 4x2in2w outputs using scopes and energy meters 
        Returns an array of energy-weighted waveforms in the order listed above (with the 2in heads in order AB, EF, GH, IJ)
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

    def _weichToPowerVsTime(weiarr):
        """
        Converts energy-weighted scope channels into a tuple of (time_values, instantaneous_power_in_watts)
        (Assumes a fixed time window of 50ns)
        """
        return (np.linspace(-5,45,len(weiarr)), np.array(weiarr)/(50e-9/len(weiarr)))
    
    @classmethod
    def _PGToPowerVsTime(cls, Psns, SSs, zzJQ):
        """
        Converts energy-weighted _PulseGoal into a tuple of (time_values, instantaneous_power_in_watts) 
        (time window chosen to match the same 50ns window of diagnostic fast photodiodes + oscilloscopes)
        Accepts pulse-specifying input parameters Psns (pulse duration segments in ns), SSs (pulse segment start/stop heights in %),
            and overal pulse energy zzJQ
        """
        PGlistwvfm=[0]*20 + list(cls._PulseGoal(DurationListQ=Psns,StartStopListQ=SSs)) + [0]*40
        PGwvfm = np.array(PGlistwvfm) * (zzJQ)/(np.sum(PGlistwvfm)) / (50e-9 / len(PGlistwvfm))
        return (np.linspace(-5,45,len(PGwvfm)),PGwvfm)
        
    def _pshostcheck():
        """
        This warns users if they're trying to do sg that requires the use specific hosts (e.g. to reach the ICS subnet)
        List of approved hosts can be found in GLOBAL.OKHOSTS
        """
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

    def _DateString():
        """Shorthand way of getting the current date in format YYYYMMDD"""
        qdate=date.today()
        return qdate.strftime('%Y%m%d')

    def get_curr_exp(timeout=15, hutch_name='mec'): 
        """Returns the name of the current experiment running in MEC and adds it to notepad PV GLOBAL.CurrExp"""
        script=pcdsdaq.ext_scripts.SCRIPTS.format(hutch_name,'get_curr_exp') 
        exp=pcdsdaq.ext_scripts.cache_script(script,timeout=timeout) 
        curr_exp=exp.lower().strip('\n')
        try:
            GLOBAL.CurrExp.put(curr_exp)
        except:
            try:
                GLOBAL.CurrExp.put(hutch_name+'xx####')
            except:
                print('Failed to write current experiment to notepad PV!')
        return curr_exp
    
    def get_curr_run(timeout=15, hutch_name='mec'): 
        """Returns the current run number in MEC and adds it to notepad PV GLOBAL.CurrRun"""
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
    
    def get_curr_shape(display=True): 
        """Returns tuple of the name of the last loaded pulse and the last time it was loaded or refreshed (as a str in '%H%M%S.%Y%m%d' format)"""
        try:
            curr_shape=GLOBAL.CurrShape.get()
            shape_time=str(GLOBAL.CurrShapeLoadTime.get())#%H%M%S.%Y%m%d
            if display:
                print('Last loaded pulse shape: {}'.format(curr_shape))
                print('Last loaded or refreshed: {}:{}:{} on {}{}-{}-{}'.format(*re.findall('..', shape_time.split('.')[0]),*re.findall('..', shape_time.split('.')[1])))
            return (curr_shape, shape_time)
        except:
            print('Failed to retrieve current shape!')
            return False

    @classmethod
    def _psheaders(cls):
        """Checks to see if shotlog files exist for the day, returns the date string and current run number"""
        cls._pshostcheck()
        DateStr=cls._DateString()
        for head in ['w','y','s1in1w','s42in1w','s42in2w','s']:
            if not os.path.exists(GLOBAL.PSFILEPATH+head+DateStr+'.p'):
                if head == 'w':
                    print('No laser file found -- probably the first shot of the day.')
                try:
                    efc.pickledump2([],GLOBAL.PSFILEPATH+head+DateStr+'.p')
                except:
                    print('Could not create file {}!'.format(GLOBAL.PSFILEPATH+head+DateStr+'.p'))
        curr_run=cls.get_curr_run()
        return (DateStr, curr_run)

    @classmethod
    def _psacqx(cls, save_flag=True, display=False):#RunNumQQ=False, 
        """
        Acquisition sequence after shooting the LPL, primary component of psposthot()
        (pspostshot() is the preferred function to use before taking a shot)
        Actions taken include preparing save folders, gathering shot data, preparing mini shot report, saving data to file, etc.
        
        save_flag=True means that the shot data will be saved to the eLog of the current experiment (in addition to internal laser records)
        save_flag=False means that the shot data will be saved to internal laser records only, NOT to user eLog
        display=True means that the acquired shot's energy-weighted, combined 2in2w waveform will be plotted as power vs. time
        display=False means that no waveform plot will be generated upon execution of the function
        """
        (DateStr, curr_run) = cls._psheaders()
        psfpQ=GLOBAL.PSFILEPATH
        RunNumStr=str(curr_run).zfill(4)
        RunName='run'+str(RunNumStr)+'_'
        headlist=['AB','EF','GH','IJ']
        #get the current experiment name
        try:
            ExpName=cls.get_curr_exp()
        except:
            ExpName=GLOBAL.CurrExp.get()
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
                          
        weiarr=cls._weichall()#fix to give np.zeros(5002) when broken

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
                np.savetxt(str(fpQ+RunName+'ch'+str(ii)+'.txt'),cls._weichToPowerVsTime(PFN.HeadENB()[ii]*weiarr[-4+ii]))

        WeightedSum=np.sum(np.expand_dims(PFN.HeadENB(),axis=1)*weiarr[-4:],axis=0)

        if save_flag:
            np.savetxt(str(fpQ+RunName+'chsum.txt'),cls._weichToPowerVsTime(WeightedSum))

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
        
        old_energies = pickle.load(open(GLOBAL.PSFILEPATH+'preshot_energies.p','rb'))
        new_energies = EMeters.EGall(return_energy_only=True);
        energy_warning = False
        for ii in range(len(old_energies)):
            if new_energies[ii] == old_energies[ii]:
                if new_energies[ii] != 0:
                    energy_warning = True
        if energy_warning:
            total_print+=efc.cstr('Caution: at least one non-zero pulse energy did not update on the previous shot!','BRY')

        if save_flag:
            np.savetxt(str(fpQ+RunName+'energies.txt'),PulseEnergies)
            ep.lxysav(*cls._weichToPowerVsTime(WeightedSum),str(fpQ+RunName+'_'+str(int(round(np.sum(PulseEnergies))))+'J'),
                     abs_path=True,xlb='Time (ns)',ylb='Power (W)')

        fileheads=['w','y','s1in1w','s42in1w','s42in2w','s'];
        prepickle=[HAWG().ReadPulseHeights(),weiarr[0],weiarr[1],weiarr[2:6],weiarr[6:10],WeightedSum]
        Psns=cls._Psns_get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        SSs=cls._SSs_get()#pickle.load(open(psfpQ+'SSs.p','rb'))
        
        if display:
            ep.llxy([cls._weichToPowerVsTime(WeightedSum),cls._PGToPowerVsTime(Psns=Psns, SSs=SSs, zzJQ=np.sum(PulseEnergies))],
                       xlb='Time (ns)',ylb='Power (W)',
                       xlim=[-1,1+np.sum(Psns)-0.25*np.sum([1 if SSs[ii][1] == SSs[ii][0] else 0 for ii in range(len(SSs)-1)])])
            
        GLOBAL.WVFMHAWG.put(prepickle[0])
        GLOBAL.WVFMYFE.put(LPL._TraceFormatting(prepickle[1],GLOBAL.LMapAB,LPL._PulseMax(Psns,SSs,np.sum(prepickle[1]))*10,AvgRange=1,FWHM=1))
        GLOBAL.WVFM1IN1w.put(LPL._TraceFormatting(prepickle[2],GLOBAL.LMapAB,LPL._PulseMax(Psns,SSs,np.sum(prepickle[2]))*10,AvgRange=1,FWHM=1))
        GLOBAL.WVFM2IN1w.put(LPL._TraceFormatting(np.sum(prepickle[3],axis=0),GLOBAL.LMapAB,LPL._PulseMax(Psns,SSs,np.sum(prepickle[3]))*10,AvgRange=1,FWHM=1))
        GLOBAL.WVFM2IN2w.put(LPL._TraceFormatting(prepickle[5],GLOBAL.LMap2,LPL._PulseMax(Psns,SSs,np.sum(prepickle[5]))*1,AvgRange=1,FWHM=1))
        for ii in range(len(fileheads)):
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
    def _psefc(cls,JreqQ=0,AQQ=0.0):
        """
        Looks at last combined 2in2w waveform and calculates new suggested update to Highland AWG to step towards _PulseGoal
        Also generates plot based on single-shot full-energy output pulse
        (function not used much for shaping anymore since we predominantly shape at 10Hz now using the function psefc10Hz)
        JreqQ is the energy scaling for the goal waveform (found using GLOBAL.PSNS and GLOBAL.SSS), used for calculating feedback and for plotting
        AQQ is the step size used in the feedback calculation; leave set to 0 unless you really know what you are doing (number should be ~<0.1)!
        """
        (DateStr, curr_run) = cls._psheaders()

        psfpQ=GLOBAL.PSFILEPATH
        Psns=cls._Psns_get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        SSs=cls._SSs_get()#pickle.load(open(psfpQ+'SSs.p','rb'))
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
                mapnow=GLOBAL.LMap2#[50,1000]
            elif len(stoday[-1]) == 1002:
                mapnow=GLOBAL.LMapAB#[5,100]
            else:
                print('Unanticipated pulse shot array length: '+str(len(stoday[-1])));
                print('Aborting...');return
            wupd=cls._wIter2(stoday[-1],np.array(wtoday[-1])/(GLOBAL.PDMAX_TEST*1.0),Psns,SSs,Jreq,mapnow,AQ)
        else:
            print('No shots yet today; beginning with pre-loaded shape')
            try:
                wupd=HAWG().ReadPulseHeights()
            except:
                print('Error! HAWG')
        return wupd
        ##EXECUTE THIS FILE FIRST TO DETERMINE THE UPDATED WAVEFORM

    @classmethod
    def psefc10Hz(cls,pwt='curr',numIterQ=50,AQQ=0.1,displayPlot=True,reloopPrompt=True,YFEbkgrdY=-.004,PtNumFront=3,PtNumBack=2,CorrFactorFront=.97,CorrFactorBack=1.0,avgfwhm=9,avgrange=1):
        """
        Looks at last 10Hz YFE waveform, calculates new suggested update to Highland AWG to step towards pulse waveform target pwt
        Also generates side-by-side plots of current YFE waveform & pwt (left) and of residual difference between the two (right)
        Use the left plot for monitoring performance, use the right plot to watch for convergence/etc.
        (This function is the underlying workhorse of the psrefrwvfm function)
        pwt is the pulse waveform target, most often _ExponentialWave2(...) or the like
        If pwt is zero then it is assumed that pwt is provided by _EW2(_Psns_get(), _SSs_get(), _YSSs_get(), offsetQ=0)
        numIterQ is the number of iterations used in the convergence loop
        AQQ is the step size used in each iteration of the convergence loop; value should be kept <<1 (usu. 0.03-0.2)
        displayPlot allows one to suppress the plots described above (=False) or allow them to be displayed (=True)
        reloopPrompt allows one the option of adding 50 more iterations after numIterQ has been reached (=True) or not (=False)
        YFEbkgrdY applies an offset of the YFE diode trace vs pwt; use if noise floor of diode trace does not line up with pwt background
        CorrFactorFront,CorrFactorBack,PtNumFront,PtNumBack are all useful if one desires to use LPL._FixEdges while converging on pwt
        """
        #avgfwhm=9;avgrange=1;#25;
        if GLOBAL.EVRLPLSSEC.get() != 43:
            GLOBAL.EVRLPLSSEN.put(0);time.sleep(0.5);
            GLOBAL.EVRLPLSSEC.put(43);time.sleep(0.5);
        if GLOBAL.EVRLPLSSEN.get() != 1:
            print('Pulse slicer not enabled! Enable now? [y/n]')
            checkprompt=efc.getch_with_TO(TOsec=10,display=False)
            if checkprompt in ('y','Y'):
                GLOBAL.EVRLPLSSEN.put(1);time.sleep(0.5);
            else:
                print('Please try again later then!')
                return False
        if GLOBAL.MBCmode.get() == 0:
            print('Bias dither enabled! Disable now? [y/n]')
            checkprompt=efc.getch_with_TO(TOsec=10,display=False)
            if checkprompt in ('y','Y'):
                GLOBAL.MBCmode.put(1);time.sleep(0.5);
            else:
                print('No? Enjoy your energy and shape fluctuations then, I guess... :/')
        Psns=cls._Psns_get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        SSs=cls._SSs_get()#pickle.load(open(psfpQ+'SSs.p','rb'))
        if pwt == 'curr': #if pwt=0 then assume pwt specified by _EW2
            pwt = cls._EW2(Psns=cls._Psns_get(), SSs=cls._SSs_get(), YSSs=cls._YSSs_get(), offsetQ=0)
        pwtF=np.array(cls._TraceFormatting(pwt,GLOBAL.pwttfmap,1,AvgRange=1,FWHM=1))
        GLOBAL.WVFMYFEGOAL.put(pwtF)

        try:
            SLA=LOSC('A');efc.ssleep();SLA._Open();efc.ssleep();SH=HAWG();efc.ssleep();SH._Open();efc.ssleep();#replaced w/LeCroyA
            tempwv=SH._ReadPulseHeights();efc.ssleep();
            SH._WritePulseHeights(140*[0]);efc.ssleep();
            bkgrdnum=20;
            iiter=0;
            while (np.sum(SLA._rch(1)) > 1) and iiter < 20:
                print('Warning: looks like the background did not clear!!')
                iiter+=1
            time.sleep(.1);
            YFEbkgrdY=0*np.array(SLA._rch(1));efc.ssleep();
            print('Acquiring background...')
            for ii in range(bkgrdnum):
                YFEbkgrdY += (np.array(SLA._rch(1))/bkgrdnum); time.sleep(0.1);
            SH._WritePulseHeights(tempwv);time.sleep(1);
                    
            ops00=SLA._rch(1)-YFEbkgrdY;time.sleep(0.1);
            print('..')
            meanerr=[]
            meanerr.append(np.sum(np.abs(pwtF[:26]-cls._TraceFormatting(ops00,GLOBAL.LMapAB,1,AvgRange=avgrange,FWHM=avgfwhm)[:26])/(pwtF[:26]+1e-3))/len(pwtF[:26]))
            ops00F=cls._TraceFormatting(ops00,GLOBAL.LMapAB,1,AvgRange=avgrange,FWHM=avgfwhm)
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
                    ops0=SLA._rch(1)-YFEbkgrdY;time.sleep(0.025);####added.215 when 200mV/div instead of 100mV/div
                    if all(ops0 == ops00):
                        print('No scope update detected... no feedback applied!')
                    else:
                        rph=SH._ReadPulseHeights();time.sleep(0.025);
                        #pwtF=np.array(TraceFormatting2(pwt,[25,500],1))
                        ##ops0F=TraceFormatting(ops0,[25,500],1)
                        ops0F=cls._TraceFormatting(ops0,GLOBAL.LMapAB,1,AvgRange=avgrange,FWHM=avgfwhm)
                        #epll([pwtF,ops0F])
                        meanerr.append(np.sum(np.abs(pwtF[:26]-ops0F[:26])/(pwtF[:26]+1e-3))/len(pwtF[:26]));
                        if displayPlot:
                            axss[0].set_data(xdat,ops0F); axs[0].relim(); axs[0].autoscale_view(True,True,True);
                            axss[1].set_data(list(range(len(meanerr))),meanerr); axs[1].relim(); axs[1].autoscale_view(True,True,True);
                            fig.canvas.draw_idle(); plt.pause(0.01);
                        usa0=cls._UpdatingShapingAlgorithm(pwtF,ops0F,np.array(rph)/(GLOBAL.PDMAX_TEST*1.0),AQQ)#.075#.25
                        usa0FE=cls._FixEdges(usa0,Psns,SSs,CorrFactorFront=CorrFactorFront,CorrFactorBack=CorrFactorBack,PtNumFront=PtNumFront,PtNumBack=PtNumBack)
                        #usa0FE=_FixEdges(usa0,[3,4.25],[[.98*100/8.0,100/8.0],[98,100]])
                        #epll([rph,usa0FE*(GLOBAL.PDMAX_TEST * 1.0)])
                        SH._WritePulseHeights(usa0FE*(GLOBAL.PDMAX_TEST *1.0));time.sleep(0.05);
                        ops00=ops0[:]
    #        if displayPlot:
    #            epll([pwtF,ops0F])
    #            epl(meanerr)
                ######check and aim
                if reloopPrompt:
                    print('Would you like to try another 50 iterations? [enter y/n]',end='',flush=True)
                    checkprompt=efc.getch_with_TO(TOsec=30,display=False);
                    if checkprompt in ('y','Y'):
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
            SH._Close();efc.ssleep();SLA._Close();efc.ssleep();
        except:
            print('Failed')
            SH._Close();efc.ssleep();SLA._Close();efc.ssleep();
        GLOBAL.CurrShapeLoadTime.put(float(datetime.now().strftime('%H%M%S.%Y%m%d')))
        if displayPlot:
            plt.ioff()

    @classmethod
    def _psupd(cls,newwavQ):
        """
        Shortcut to updating the Highland AWG using the provided input waveform newwavQ
        """
        cls._pshostcheck()
        wupdt=newwavQ[:]
        if max(wupdt) < 1.5:
            wupdt=(GLOBAL.PDMAX_TEST * 1.0)*np.array(wupdt)
        try:
            HAWG().WritePulseHeights(wupdt);
        except:
            print('Error, check HAWG!')

    @classmethod
    def psloadwvfm(cls,RecipeStrQ,WvGoal10HzHint=False):
        """
        Loads a new waveform according to previously-saved recipe RecipeStrQ
        Options for RecipeStrQ can be found using the LPL.psrecipes() command
        (LPL.psmenu() allows you to choose a recipe to load without needing to type in the full name of it)
        
        WvGoal10HzHint=True means that you will see the YFE waveform hint needed to make the pulse
            and upon which the execution of LPL.psrefrwvfm() would be based; usually based on _ExponentialWave2
            
        WvGoal10HzHint=False means that the YFE waveform hint will not be printed to terminal
        """
        cls._pshostcheck()
        print('Loading timestamp: '+datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
        try:
            [Psns,SSs,YFE02mmCurr,YFE06mmCurr,YFE10mmCurr,NewWvfm,WvGoal10Hz] = pickle.load(open(GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'.p','rb'))
        except:
            print('Recipe file '+GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'.p\' not found.')
            return
        efc.pickledump2(Psns,GLOBAL.PSFILEPATH+'Psns.p')
        cls._Psns_set(Psns)
        efc.pickledump2(SSs,GLOBAL.PSFILEPATH+'SSs.p')
        cls._SSs_set([elem for sublist in SSs for elem in sublist])
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
        yfegoal=cls._LinearWave2(500,0,1025,0,0,5002);
        if (len(pslparams) > 0) or (len(pseparams) > 0):
            for ii in range(len(pslparams)): 
                yfegoal+=cls._LinearWave2(pslparams[ii][0],pslparams[ii][1],pslparams[ii][2],pslparams[ii][3],0,5002)
            for ii in range(len(pseparams)): 
                yfegoal+=cls._ExponentialWave2(pseparams[ii][0],pseparams[ii][1],pseparams[ii][2],pseparams[ii][3],0,5002)
        else:
            print('No wave extracted: '+WvGoal10Hz)
        cls._YSSs_set([[pseparams[ii][1],pseparams[ii][3]] for ii in range(len(pseparams))])
        yfegoalF=np.array(cls._TraceFormatting(yfegoal,GLOBAL.pwttfmap,1,AvgRange=1,FWHM=1))
        GLOBAL.WVFMYFEGOAL.put(yfegoalF)
        GLOBAL.WVFM2IN2wGOAL.put(cls._PulseGoal(Psns,SSs))
        try:
            cls._psupd(NewWvfm)
            GLOBAL.CurrShape.put(RecipeStrQ)
            GLOBAL.CurrShapeLoadTime.put(float(datetime.now().strftime('%H%M%S.%Y%m%d')))
            print('New waveform loaded! ')
        except:
            print('Failed to load new waveform.')

    @classmethod
    def pssavewvfm(cls,RecipeStrQ=0,PsnsQ=0,SSsQ=0,YFEgetQ=0,TargetwlistDateQ='curr',TargetwindexQ=0,WvGoal10HzHint='none',YFEbkgrdY=-.004):
        """
        Saves a new pulse shape recipe using user-provided parameters
        if setting RecipeStrQ='0' or = of pattern (^E\d+J) (e.g., 'E40J'): function will guess the apropriate name for the recipe based on PsnsQ and SSsQ
            otherwise: set RecipeStrQ equal to the name of the recipe you want
            Example: leaving RecipeStrQ blank or using RecipeStrQ=0 results in a name of '10ns00grad' if PsnsQ=[10.25] and SSsQ=[[98,100]]
            Example: using RecipeStrQ='E20J' results in a name of '04ns06ns020to10stepE20J' if PsnsQ=[4,6.25] and SSsQ=[[19,20],[98,100]]
            Example: using RecipeStrQ='ThanosPulse1337' results in a name of 'ThanosPulse1337'
        
        if setting PsnsQ=0 and SSsQ=0: function will load PsnsQ and SSsQ from GLOBAL.PSNS and GLOBAL.SSS, respectively
            otherwise: set PsnsQ and SSsQ equal to the values you would like included in your recipe
            Example: leaving PsnsQ and SSsQ blank or setting them equal to 0 results in 
                     PsnsQ=[20.25] and SSsQ=[[85,100]] if GLOBAL.PSNS is [20.25] and GLOBAL.SSS is [[85,100]]
            Example: setting PsnsQ=[2.25,2.25,2.25,2.25,2.25] and SSsQ=[[10,15],[15,25],[25,40],[40,80],[80,100]] sets the
                     values of PsnsQ and SSsQ explicitly without regard for GLOBAL.PSNS and GLOBAL.SSS
        
        if setting YFEgetQ=0: function will read out and save current YFE eDrive currents as part of the recipe
            otherwise: these can be saved explicitly, but there is almost never a good reason to do this anymore
        
        if setting TargetwlistDateQ='curr': current Highland AWG waveform shot will be used to make the recipe
            otherwise: set TargetwlistDateQ equal to the desired date and TargetwindexQ equal to the desired shot number
            Example: TargetwlistDateQ='20220211' and TargetwindexQ=3 would make a recipe using the fourth shot on 2022Feb11
            
        IMPORTANT: please set WvGoal10HzHint equal to the target waveform for the YFE output
            Failure to do so will prevent psrefrwvfm() from working later, as it will have no provided target waveform 
            Example: use WvGoal10HzHint='_ExponentialWave2(500,.1,1025,.8,0,5002)' if _ExponentialWave2(500,.1,1025,.8,0,5002) was
                     the waveform used in psefc10Hz while you were creating the waveform for your recipe
        """
        cls._pshostcheck()
        if not isinstance(WvGoal10HzHint, str):
            print('Warning: WvGoal10HzHint must be a string. Please put quotes around your shaping hint and try again!')
            return False
        if PsnsQ == 0:
            PsnsQ=cls._Psns_get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        if SSsQ == 0:
            SSsQ=cls._SSs_get()#pickle.load(open(psfpQ+'SSs.p','rb'))
            
        if WvGoal10HzHint=='none':#try to figure out WvGoal10HzHint from YSSs
            YSSsQ=cls._YSSs_get()
            WvGoal10HzHint=cls._EW2stringhint(Psns=PsnsQ,SSs=SSsQ,YSSs=YSSsQ,YFEbkgrdY=YFEbkgrdY)

        print('Saving timestamp: '+datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
        if TargetwlistDateQ == 'curr':
            print('Using current Highland waveform...')
            try:
                NewWvfmQ=HAWG().ReadPulseHeights();
            except:
                print('Failed! HAWG');
                return False
            try:
                wlastarr=pickle.load(open(GLOBAL.PSFILEPATH+'w'+cls._DateString()+'.p','rb'))
                if wlastarr[-1] == NewWvfmQ:
                    print('Pulse looks equivalent to the most recent pulse, w'+cls._DateString()+'['+str(len(wlastarr)-1)+']')
                    WvGoal10HzHint=WvGoal10HzHint+';; w'+cls._DateString()+'['+str(len(wlastarr)-1)+']'
                else:
                    WvGoal10HzHint=WvGoal10HzHint+';; sometime after most recent w'+cls._DateString()+'['+str(len(wlastarr)-1)+']'
            except:
                print('Failed to load most recent amplified shot.')
        else:
            wavehistQ=pickle.load(open(GLOBAL.PSFILEPATH+'w'+str(TargetwlistDateQ)+'.p','rb'))
            NewWvfmQ=wavehistQ[TargetwindexQ][:]
            WvGoal10HzHint=WvGoal10HzHint+', w'+str(TargetwlistDateQ)+'['+str(TargetwindexQ)+']'


        if YFEgetQ == 0:
            YFEgetQ=YFE.Get(display=False)
        [YFE02mmCurrQ,YFE02mmCurrQ,YFE02mmCurrQ,YFE02mmCurrQ,YFE06mmCurrQ,YFE10mmCurrQ]=YFEgetQ
        
        if (RecipeStrQ == 0) or (len(re.findall('(^E\d+J)',RecipeStrQ))>0):#learn formatting for pulse from Psns and SSs
            if isinstance(RecipeStrQ, str):
                if (RecipeStrQ[0] in ('e','E')):#allows to tag on an energy description, e.g. E20J for a 20J version of a pulse
                    RecipeStrQend=RecipeStrQ
                else:
                    RecipeStrQend=''
            else:
                RecipeStrQend=''
            RecipeStrQ=''
            if len(PsnsQ) == 1:#if only one segment, it'll be ##ns##grad [+any energy description]
                RecipeStrQ+='{:02}ns'.format(round(PsnsQ[0]))
                gradval=round(SSsQ[0][1]-SSsQ[0][0])
                RecipeStrQ+='{:02}grad'.format(0 if gradval < 3 else 99 if gradval>99 else gradval)
            elif len(PsnsQ) == 2:#if two segments, it'll be ##ns##ns###to100step [+any grad info] [+any energy description]
                RecipeStrQ+='{:02}ns{:02}ns'.format(round(PsnsQ[0]),round(PsnsQ[0]))
                RecipeStrQ+='{:03}to100step'.format(round(SSsQ[0][1]))
                if ((SSsQ[0][1] - SSsQ[0][0]) > 2) or ((SSsQ[1][1] - SSsQ[1][0]) > 2):
                    gradval=round(SSsQ[0][1]-SSsQ[0][0])
                    RecipeStrQ+='{:02}grad'.format(0 if gradval < 3 else 99 if gradval>99 else gradval)
                    gradval2=round(SSsQ[1][1]-SSsQ[1][0])
                    RecipeStrQ+='{:02}grad'.format(0 if gradval2 < 3 else 99 if gradval2>99 else gradval2)
            else:#otherwise it'll be ##ns###to100ramp [+any energy description]
                RecipeStrQ+='{:02}ns'.format(round(np.sum(PsnsQ)-0.25*np.sum([1 if SSsQ[ii][1] == SSsQ[ii][0] else 0 for ii in range(len(SSsQ)-1)])))
                RecipeStrQ+='{:03}to100ramp'.format(round(SSsQ[0][0]))
            RecipeStrQ+=RecipeStrQend
        oldfilefound=True

        try:
            dummyQ = pickle.load(open(GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'.p','rb'))
            print('Old recipe found with same name: '+GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'.p')
        except:
            oldfilefound=False
        iiQ=0
        while oldfilefound:
            iiQ=iiQ+1
            try:
                pickle.load(open(GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'_'+str(iiQ).zfill(2)+'.p','rb'));
            except:
                oldfilefound=False
                dummyQ[-1]='**replaced on '+cls._DateString()+'** '+dummyQ[-1]
                efc.pickledump2(dummyQ,GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'_'+str(iiQ).zfill(2)+'.p')
                print('Saved old recipe as '+GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'_'+str(iiQ).zfill(2)+'.p')
        try:
            efc.pickledump2([PsnsQ,SSsQ,YFE02mmCurrQ,YFE06mmCurrQ,YFE10mmCurrQ,NewWvfmQ,WvGoal10HzHint],GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'.p')
            print('Saved new recipe as '+GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'.p')
        except:
            print('Failed to save new recipe.')

    @classmethod
    def psviewwvfm(cls,RecipeStrQ='none',TargetwlistDateQ='curr',TargetwindexQ=0,WvGoal10HzHint=False):
        """
        Displays waveform and parameters that are part of a previously-saved recipe
        (LPL.psmenu() allows you to choose a recipe to view without needing to type in the full name of it)
        
        Set RecipeStrQ equal to the recipe you want to view, e.g. RecipeStrQ='10ns00grad'
            :if no RecipeStrQ is provided, displayed pulse will be chosen according to TargetwlistDateQ and TargetwindexQ
        
        if using TargetwlistDateQ='curr' and TargetwindexQ=0 with RecipeStrQ='none' or blank: most recent pulse will be shown
            :otherwise choose the date and shot number of the pulse desired
            :Example: TargetwlistDateQ='20220214' and TargetwindexQ=6 will display the seventh saved pulse from 2022Feb14
            
        if using WvGoal10HzHint=True: any saved waveform hints will be printed to terminal; setting WvGoal10HzHint=False results in no printout
            :this is useful primarily when hoping to make adjustments to a recipe or to create a completely new recipe
            :i.e. use the hint to make modifications to the returned _ExponentialWave2(...) parameters suitable for your new needs
        """
        cls._pshostcheck()
        if RecipeStrQ == 'none':
            if TargetwlistDateQ == 'curr':
                foundlastshot=False
                iidQ=0
                while not foundlastshot:
                    try:
                        wlastarr=pickle.load(open(GLOBAL.PSFILEPATH+'w'+str(int(cls._DateString())-iidQ)+'.p','rb'))
                        slastarr=pickle.load(open(GLOBAL.PSFILEPATH+'s'+str(int(cls._DateString())-iidQ)+'.p','rb'))
                        wlast=wlastarr[-1][:]
                        slast=slastarr[-1][:]
                        print('Retrieving most recent shot: w'+str(int(cls._DateString())-iidQ)+'['+str(len(wlastarr)-1)+']')
                        foundlastshot=True
                    except:
                        iidQ=iidQ+1
            else:
                try:
                    wlastarr=pickle.load(open(GLOBAL.PSFILEPATH+'w'+str(TargetwlistDateQ)+'.p','rb'))
                    slastarr=pickle.load(open(GLOBAL.PSFILEPATH+'s'+str(TargetwlistDateQ)+'.p','rb'))
                    wlast=wlastarr[int(TargetwindexQ)]
                    slast=slastarr[int(TargetwindexQ)]
                except:
                    print('Failed to load at given date and index: '+TargetwlistDateQ+', '+str(TargetwindexQ))
            Psns=cls._Psns_get()
            SSs=cls._SSs_get()
            ep.l(wlast)
            ep.llxy([cls._weichToPowerVsTime(slast), cls._PGToPowerVsTime(Psns=Psns, SSs=SSs, zzJQ=np.sum(slast))],
                       xlb='Time (ns)',ylb='Power (W)',
                       xlim=[-1,1+np.sum(Psns)-0.25*np.sum([1 if SSs[ii][1] == SSs[ii][0] else 0 for ii in range(len(SSs)-1)])])
            return 
        else:
            try:
                [Psns,SSs,YFE02mmCurr,YFE06mmCurr,YFE10mmCurr,NewWvfm,WvGoal10Hz] = pickle.load(open(GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'.p','rb'))
            except:
                print('Recipe file '+GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'.p\' not found.')
                return False
            print('Retrieved recipe: load'+RecipeStrQ)
            print('Psns: '+str(list(Psns))+', SSs: '+str([list(SSs_sub) for SSs_sub in SSs]))
            print('YFEcurr:: 2mm: '+'{:5.1f}'.format(YFE02mmCurr)+', 6mm: '+'{:5.1f}'.format(YFE06mmCurr)+', 10mm: '+'{:5.1f}'.format(YFE10mmCurr))
            print('Extra info: '+WvGoal10Hz)
            ep.l(NewWvfm)
            try:
                tempstr=WvGoal10Hz[-18:]
                wstrinx=tempstr.index('w')
                wstrinx0=wstrinx-len(tempstr)
                loaddate=tempstr[wstrinx0+1:wstrinx0+9]
                loadindx=int(tempstr[wstrinx0+10:-1])
                saskarr=pickle.load(open(GLOBAL.PSFILEPATH+'s'+loaddate+'.p','rb'))
                ep.llxy([cls._weichToPowerVsTime(saskarr[loadindx]),cls._PGToPowerVsTime(Psns=Psns, SSs=SSs, zzJQ=np.sum(saskarr[loadindx]))],
                           xlb='Time (ns)',ylb='Power (W)', 
                           xlim=[-1,1+np.sum(Psns)-0.25*np.sum([1 if SSs[ii][1] == SSs[ii][0] else 0 for ii in range(len(SSs)-1)])])
                print('Pulse energy was ~{}J'.format(np.sum(saskarr[loadindx])))
            except:
                print('Failed to load 2w waveform for display.')
            return NewWvfm

    @classmethod
    def psrefrwvfm(cls,RecipeStrQ='latest',numStepsQ=50,stepSizeQ=0.1,YFEbkgrdY=-.004,displayPlot=True,reloopPrompt=True):
        """
        Refreshes 10Hz YFE waveform according to target shape given by previously-saved recipe
        
        If you do not specify a recipe in RecipeStrQ, it will automatically begin refreshing the shape loaded most recently
            :otherwise: use RecipeStrQ='10ns00grad' to begin refreshing from the 
             current waveform towards the YFE goal waveform for the 10ns00grad pulse
        
        Use numStepsQ to specify the number of 10Hz iterations to take 
        Use stepSizeQ to specify the size of the corrective step per iteration; should be <<1 (usu. 0.01-0.3)
        Use displayPlot=True to display the waveform as it converges to the goal (as with psefc10Hz)
            :using displayPlot=False will not show the plot
        Use reloopPrompt=True to give the operator an opportunity to add more iterations after the first numStepsQ have elapsed
            :using reloopPrompt=False will cause the function to terminate as soon as numStepsQ have elapsed
        """
        cls._pshostcheck()
        print('Loading timestamp: '+datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
        if not YFE.OnCheck(display=False):
            print('YFE does not appear to be on! Check YFE status first!')
            return False
        #load and extract the pulse target from the desired recipe
        if RecipeStrQ=='latest':#if refreshing most recent wvfm, just use _EW2 and RBV of Psns, SSs, and YSSs
            yfegoal=cls._EW2(Psns=cls._Psns_get(), SSs=cls._SSs_get(), YSSs=cls._YSSs_get())
        else:#if not refreshing most recent wvfm
            try:
                [Psns,SSs,YFE02mmCurr,YFE06mmCurr,YFE10mmCurr,NewWvfm,WvGoal10HzHint] = pickle.load(open(GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'.p','rb'))
            except:
                print('Recipe file '+GLOBAL.PSFILEPATH+'recipes/load'+RecipeStrQ+'.p\' not found.')
                return
            print('Hint text: '+WvGoal10HzHint)
            pseparams=np.array(re.findall('ExponentialWave2\((\d+),(\d+\.\d+|\.\d+),(\d+),(\d+\.\d+|\.\d+),0,5002\)',WvGoal10HzHint),dtype=np.float32);
            pslparams=np.array(re.findall('LinearWave2\((\d+),(\d+\.\d+|\.\d+),(\d+),(\d+\.\d+|\.\d+),0,5002\)',WvGoal10HzHint),dtype=np.float32);
            try:
                YFEbkgrdY=float(re.findall('YFEbkgrdY\s?=\s?([-+]?\d+\.?\d*|[-+]?\.\d+)',WvGoal10HzHint)[0]);
            except:
                YFEbkgrdY=YFEbkgrdY
            yfegoal=cls._LinearWave2(500,0,1025,0,0,5002);
            if (len(pslparams) > 0) or (len(pseparams) > 0):
                for ii in range(len(pslparams)): 
                    yfegoal+=cls._LinearWave2(pslparams[ii][0],pslparams[ii][1],pslparams[ii][2],pslparams[ii][3],0,5002)
                for ii in range(len(pseparams)): 
                    yfegoal+=cls._ExponentialWave2(pseparams[ii][0],pseparams[ii][1],pseparams[ii][2],pseparams[ii][3],0,5002)
            else:
                print('No wave extracted: '+WvGoal10HzHint)
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
                    YFE.On(CtrlChk=False);
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
            print('Set and enable 10Hz output')
            if GLOBAL.EVRLPLSSEC.get() != 43:
                GLOBAL.EVRLPLSSEN.put(0);time.sleep(0.75);
                GLOBAL.EVRLPLSSEC.put(43);time.sleep(0.75);
            if GLOBAL.EVRLPLSSEN.get() != 1:
                GLOBAL.EVRLPLSSEN.put(1);time.sleep(0.75);
            #run the update code
            print('Refreshing the YFE wavefront...')
            cls.psefc10Hz(pwt=yfegoal,numIterQ=numStepsQ,AQQ=stepSizeQ,displayPlot=displayPlot,reloopPrompt=reloopPrompt,YFEbkgrdY=YFEbkgrdY)
            #reset to single shot on pulse picker
            GLOBAL.EVRLPLSSEC.put(182);time.sleep(0.75);
            GLOBAL.EVRLPLSSEN.put(1);time.sleep(0.75);
            #re-open shutters
            print('Opening all shutters...')
            TTL_shutter.Toggle('openall',display=False);#open all the shutters
            MBC.Reset();
            YFE.SetAll(True,displayQ=False);
        except:#used to make sure shutters re-open even in case of error or KeyboardInterrupt
            #reset to single shot on pulse picker
            GLOBAL.EVRLPLSSEC.put(182);time.sleep(0.75);
            GLOBAL.EVRLPLSSEN.put(1);time.sleep(0.75);
            #re-open shutters
            #print('Opening all shutters...')
            #toggle_TTL_shutter('openall',display=False);#open all the shutters
            MBC.Reset();
            YFE.SetAll(True,displayQ=False);

    @classmethod
    def psrecipes(cls):
        """
        Prints and returns a list of all previously-saved pulse recipes
        LPL.psmenu() shows this list but also allows you to load or view a recipe
        """
        allrec=glob.glob(GLOBAL.PSFILEPATH+'recipes/*.p');
        oldrec=glob.glob(GLOBAL.PSFILEPATH+'recipes/*_*.p')
        currec=[ext[60:-2] for ext in allrec if ext not in oldrec];currec.sort();
        return currec
    
    @classmethod
    def psmenu(cls):
        """
        Allows you to select a recipe to load or view
        """
        recipelist=cls.psrecipes()
        print('Pulse recipes:')
        for pair in list(zip(list(range(len(recipelist))),recipelist)):
            if pair[0]:
                print(pair)
        print('{}{}{}'.format('Would you like to load or view a recipe? Enter L7 to load recipe 7.\n',
                        'Enter V16 to view recipe 16 first before having the option to load it.\n',
                        'Enter Q to quit.'))
        checkprompt=efc.input_with_TO(TOsec=30,display=False)
        if checkprompt:
            if checkprompt[0] not in ('v','V','l','L'):
                print('Exiting...')
                return
            try:
                recipenum=int(checkprompt[1:])
                if recipenum > len(recipelist):
                    raise Exception
            except:
                print('I\'m sorry, that\'s not a valid recipe number. Try again!')
                return False
            if checkprompt.lower().startswith('l'):
                print('Loading {} pulse!'.format(recipelist[recipenum]))
                cls.psloadwvfm(recipelist[recipenum])
                return
            elif checkprompt.lower().startswith('v'):
                print('Displaying {} pulse!'.format(recipelist[recipenum]))
                cls.psviewwvfm(recipelist[recipenum],WvGoal10HzHint=True);
                #add load option back in after fixing plot blocking...
                #print('Would you like to load it? [y/n]')
                #checkprompt2=efc.getch_with_TO(20,display=False)
                #if checkprompt2 not in ('y','Y'):
                #    print('Hope you can find someething else you like! :D')
                #    return
                #else:
                #    print('Loading {} pulse!'.format(recipelist[recipenum]))
                #    cls.psloadwvfm(recipelist[recipenum])
                #    return
                return
            else:
                print('Exiting menu...')
                return
        else:
            print('Call me back when you\'re ready to order! :D')#handles initial prompt timeout
    

    @classmethod
    def pspreshot(cls,MBC_bypass=False):
        """
        Prepares and checks state of the laser for taking a single full-energy shot
        Because of checks/safeguards, this is the preferred function to use before taking a shot
        """
        cls._pshostcheck()
        if not YFE.OnCheck(display=False):
            print('WARNING: YFE seems to be off... Attempting to turn on YFE...')
            YFE.On(CtrlChk=False);
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
                prechk = True
            else:
                print('*')
                print('WARNING: MBC not responding!!')
                prechk = False
        else:
            if not MBC_bypass:
                print('MBC is not safe! Resetting the MBC...')
                MBC.Reset();##up above, check emission AND check currents???
                YFE.SetAll(True,displayQ=False)
                returnval=cls.pspreshot()
                return returnval
            else:
                prechk = True
        if GLOBAL.CurrShape.get().lower() == 'yfeedge':
            print('{}{}'.format(efc.cstr('WARNING: ','BRW,BBRR'),
                                efc.cstr('last loaded pulse was YFEedge! Are you SURE you want to continue? [y/n]','BRW,BBRR')))
            resp=efc.getch()
            if resp.lower()[0] == 'y':
                print('I '+efc.cstr('REALLY','BRR')+' hope you know what you\'re doing...')
            else:
                print('Good choice. Please try again when you are ready.')
                return False
        #if GLOBAL.EVRLPLLAMPEC.get() != 182:
            #GLOBAL.EVRLPLLAMPEN.put(0)
        GLOBAL.EVRLPLLAMPEC.put(182)
        #if GLOBAL.EVRLPLSSEC.get() != 182:
            #GLOBAL.EVRLPLSSEN.put(0)
        GLOBAL.EVRLPLSSEC.put(182)
        #if GLOBAL.EVRLPLLAMPEN.get() != 1:
        GLOBAL.EVRLPLLAMPEN.put(1)
        #if GLOBAL.EVRLPLSSEN.get() != 1:
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
        Psns=cls._Psns_get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        SSs=cls._SSs_get()#pickle.load(open(psfpQ+'SSs.p','rb'))
        print('Current pulse target is: '+str(list(Psns))+' ns, '+str([list(SSs_sub) for SSs_sub in SSs])+' % of max power.')
        print('The most recently loaded pulse shape is: '+GLOBAL.CurrShape.get()+' (load/refresh stamp: {:015.8f})'.format(GLOBAL.CurrShapeLoadTime.get()))
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

        TTL_shutter.Toggle('open'+yes_enabled+'wwxx',display=False);time.sleep(1.5);#make sure all shutters are open...
        TTL_shutter.Toggle('close'+not_enabled,display=False);time.sleep(1.5);#close shutters that aren't enabled
        print('Current shutter status:')
        TTL_shutter.Status(display=True);
        efc.pickledump2(EMeters.EGall(return_energy_only=True), GLOBAL.PSFILEPATH+'preshot_energies.p')#used to check if energy meters update or not
        return prechk
        #waveform pre-check? verify shutters are open?

    @classmethod
    def pspostshot(cls,save_flag=True,display=False):
        """
        Executes post-shot routine for saving data and returning laser to appropriate state
        This is the preferred function to use after taking a shot (_psacqx() is inside)
        
        save_flag=True means that the shot data will be saved to the eLog of the current experiment (in addition to internal laser records)
        save_flag=False means that the shot data will be saved to internal laser records only, NOT to user eLog
        display=True means that the acquired shot's energy-weighted, combined 2in2w waveform will be plotted as power vs. time
        display=False means that no waveform plot will be generated upon execution of the function
        """
        GLOBAL.EVRLPLLAMPEN.put(0)
        GLOBAL.EVRLPLSSEN.put(0)
        cls._pshostcheck()
        cls._psacqx(save_flag=save_flag,display=display)#took out _noLecroyA
        #cls._psefc();
        TTL_shutter.Toggle('openall',display=False);#make sure all shutters are open again...
        print('Resetting bias tracking...')
        MBC.Reset();
        YFE.SetAll(True);
        EMeters.E_synth_refresh();

    @classmethod
    def SHG_opt(cls,armsQ='ABEFGHIJ'):#check for trace height;#All shutters must start in the open state... 
        """
        Optimizes the tuning of the doubler angles to maximize the conversion efficiency of the arms of the LPL
        (SHG_opt is useful especially if ambient conditions change; sometimes needs to be run fairly regularly)
        Prompts included to help insure system is in appropriate state for optimization
        Entire function should take only a couple of minutes; live display allows you to monitor progress
        If new value is too close to edge of window, it is recommended to optimize that arm again
        (If the doublers are way far off (no detection), tune the Newport motor back until it is closer)
        
        If leaving armsQ blank: all arms are optimized in the order of AB, EF, GH, IJ
            :EXAMPLE: SHG_opt() is equivalent to SHG_opt('ABEFGHIJ') and optimizes all arms
        Instead, one may choose specific arms to optimize
            :EXAMPLE: SHG_opt(armsQ='ABIJ') or SHG_opt('ABIJ') optimizes only arms AB and IJ
            :EXAMPLE: SHG_opt(armsQ='EF') or SHG_opt('EF') optimizes only arm EF
        """
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
        checkprompt=efc.getch_with_TO(TOsec=10,display=False);
        if checkprompt not in ('y','Y'):
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
        optwvfm=pickle.load(open(GLOBAL.PSFILEPATH+'opttrace.p','rb'));
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
                checkprompt=efc.getch_with_TO(TOsec=10,display=False);
                if checkprompt not in ('n','N'):
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
            SLA=LOSC('A');SLA._Open();#changed to LecroyA since repair
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
            SLA._Close();time.sleep(.15);#changed to LeCroyA
        except:
            print('Failed! Restoring original values and attempting to re-open most-recent shutter... you should verify!')
            SLA._Close();time.sleep(.15);#changed to LeCroyA
            if currentshutter > 0:
                TTL_shutter.Toggle('open'+armlist[currentshutter],display=False);
            for ii in range(4):
                SHGpvlist[ii].put(startposlist[ii]);newposlist[ii]=startposlist[ii];
        time.sleep(2);#need time so that last shutter trigger ends before trying to open IJ
        try:
            HAWG().WritePulseHeights(oldwvfm);
        except:
            print('Error! Check waveform!')
        GLOBAL.EVRLPLSSEN.put(0);#disable PC before re-opening shutters
        datestamp=int(datetime.now().strftime('%Y%m%d%H%M%S'))
        SHGlog=pickle.load(open(GLOBAL.PSFILEPATH+'SHG_opt_log.p','rb'))
        SHGlog.append([datestamp,[newposlist[ii] for ii in range(4)]])
        efc.pickledump2(SHGlog,GLOBAL.PSFILEPATH+'SHG_opt_log.p')
        TTL_shutter.Toggle('openall',display=False);#open all the shutters
        MBC.Reset();YFE.SetAll(True);#reset bias...
        plt.ioff()
    
    def On():
        """Turns on LPL front-end laser system; shortcut for YFE.On() command"""
        result = YFE.On()
        return result

    def Off():
        """Turns off LPL front-end laser system; shortcut for YFE.Off() command"""
        result = YFE.Off()
        return result
        
        
class ep:
    """
    Easy Plotting class for convenient shorthand ways of plotting data; nothing special, just lazy!
    Typical usage via ep.[command]
    Possible commands include:
        :l #quick plot of single array of y values
        :lxy #quick plot of single array of y values with x values also specified
        :lxyloglog #quick xy plot with log x- and y-axes
        :lsav #save plot of single array of y values
        :lxysav #save plot of single array of y values with x values also specified
        :lcomp #quick plot meant to compare diode waveform to pulse goal
        :lcsv #quick plot of single array from csv file
        :llcsv #quick plot of list of arrays from csv file
        :rcsv #read in array from csv file
        :ll #quick plot of list of arrays of y values
        :llxy #quick plot of list of arrays of y values with x values also specified
        :llxyloglog #quick xy plot of several lists with log x- and y-axes
        :llt #plot list of arrays of values according to a time mapping 
        :llcomp #plot list of diode waveforms along with target waveform
        :lfft #quick plotting of FFT of provided waveform
        :ShotVsGoal #quick plotting of LPL shot vs goal waveform
    potential future work
        - consider consolidating some functions into one
            - e.g. getch with and without TimeOut?
    """
    def l(listq,xlb='none',ylb='none',xlim='none',ylim='none'):
        """
        Shorthand plotting function for a single list of y-values
        Example: ep.l([1,2,4]) generates a plot of [1,2,4]
        Optional: Set the labels for the x-axis and y-axis with xlb and ylb, respectively
        Optional: Set the limits for the x-axis and y-axis with xlim and ylim, respectively
        """
        df1=plt.figure()
        plt.plot(listq);
        if xlb != 'none':
            plt.xlabel(xlb)
        if ylb != 'none':
            plt.ylabel(ylb)
        if xlim != 'none':
            try:
                plt.xlim(xlim)
            except:
                pass
        if ylim != 'none':
            try:
                plt.ylim(ylim)
            except:
                pass
        df1.show()
        return

    def lxy(listxq,listyq,xlb='none',ylb='none',xlim='none',ylim='none'):
        """
        Shorthand plotting function for a single list of x-values and y-values
        Be sure to include the list of (equal-length) x-values and y-values separately
        Example: ep.lxy([10,20,30],[1,5,18]) generates a plot of points (10,1), (20,5), and (30,18)
        Optional: Set the labels for the x-axis and y-axis with xlb and ylb, respectively
        Optional: Set the limits for the x-axis and y-axis with xlim and ylim, respectively
        """
        df1=plt.figure()
        plt.plot(listxq,listyq);
        if xlb != 'none':
            plt.xlabel(xlb)
        if ylb != 'none':
            plt.ylabel(ylb)
        if xlim != 'none':
            try:
                plt.xlim(xlim)
            except:
                pass
        if ylim != 'none':
            try:
                plt.ylim(ylim)
            except:
                pass
        df1.show()
        return

    def lxyloglog(listxq,listyq,xlb='none',ylb='none',xlim='none',ylim='none'):
        """
        Shorthand plotting function for a single list of x-values and y-values with logrithmic axes in both directions
        Be sure to include the list of (equal-length) x-values and y-values separately
        Example: ep.lxyloglog([10,100,1000],[.3,.05,.001]) generates a log-log plot of points (10,.3), (100,.05), and (1000,.001)
        Optional: Set the labels for the x-axis and y-axis with xlb and ylb, respectively
        Optional: Set the limits for the x-axis and y-axis with xlim and ylim, respectively
        """
        df1=plt.figure()
        plt.loglog(listxq,listyq);
        if xlb != 'none':
            plt.xlabel(xlb)
        if ylb != 'none':
            plt.ylabel(ylb)
        if xlim != 'none':
            try:
                plt.xlim(xlim)
            except:
                pass
        if ylim != 'none':
            try:
                plt.ylim(ylim)
            except:
                pass
        df1.show()
        return

    def lsav(listq,FileNameQ,blockdisplay=True):
        """
        Shorthand function for saving a plot of a single list of y-values as a .png file
        Example: ep.lsav([1,2,4,8],'/reg/neh/operator/mecopr/my_plot') saves the plot of [1,2,4,8]
                 to the file name and path '/reg/neh/operator/mecopr/my_plot'
                 (there is no need to include the file extension -- it is added inside the code)
        Optional: blockdisplay=True closes the figure after creating and saving; 
                  blockdisplay=False means the figure will be saved AND displayed on the screen
        """
        df1=plt.figure()
        plt.plot(listq);
        df1.savefig(str(FileNameQ+'.png'))
        if blockdisplay:
            plt.close(df1)
        return
        
    def lxysav(listxq,listyq,FileNameQ,abs_path=False,xlb='none',ylb='none',xlim='none',ylim='none'):
        """
        Shorthand function for saving a plot of a single list of x-values and y-values
        Be sure to include the list of (equal-length) x-values and y-values separately
        Example: ep.lxysav([10,20,30],[1,5,18],'/reg/neh/operator/mecopr/my_plot',abs_path=True)
                 saves the plot of points (10,1), (20,5), and (30,18) to the file name and
                 path '/reg/neh/operator/mecopr/my_plot'
                 (there is no need to include the file extension -- it is added inside the code)
        Optional: Set abs_path=False to save the plot to your current working directory;
                  Set abs_path=True to specify the absolute path where you will save your plot
        Optional: Set the labels for the x-axis and y-axis with xlb and ylb, respectively
        Optional: Set the limits for the x-axis and y-axis with xlim and ylim, respectively
        """
        df1=plt.figure()
        plt.plot(listxq,listyq);
        if abs_path:
            figfilename=FileNameQ;
        else:
            figfilename=str(GLOBAL.PSFILEPATH+FileNameQ+'.png')
        if xlb != 'none':
            plt.xlabel(xlb)
        if ylb != 'none':
            plt.ylabel(ylb)
        if xlim != 'none':
            try:
                plt.xlim(xlim)
            except:
                pass
        if ylim != 'none':
            try:
                plt.ylim(ylim)
            except:
                pass
        df1.savefig(figfilename)        
        plt.close(df1)
        return
        
    def lcomp(listq,goalq,Map,tMax):
        """
        Shorthand function for comparing a pulse shape to its targeted pulse shape; not used much currently
        Map and tMax are the parameters used for formatting the data in listq according to the desired new x- and y-range
        goalq is the waveform to which the formatted listq will be compared
        """
        formtra=[]
        formtra.append(LPL._TraceFormatting(listq,Map,tMax, AvgRange=1, FWHM=1))
        formtra.append(goalq)
        ep.ll(formtra)
        return

    @classmethod
    def lcsv(cls,CSVname):
        """
        Shorthand plotting function for a single list of y-values loaded from a CSV file; hasn't been used for ages
        File should be located within the data folder along GLOBAL.PSFILEPATH, i.e. '/reg/neh/operator/mecopr/mecpython/pulseshaping/data'
        """
        with open(GLOBAL.PSFILEPATH+'data/'+CSVname+'.csv','r') as filehead:
            RawListQ=filehead.read()
            ListedValues=RawListQ.split('\n')
        cls.l(ListedValues[:-1])
        return 

    @classmethod
    def llcsv(cls,CSVHeadname):
        """
        Shorthand plotting function for a nested list of four sets of y-values loaded from a CSV file; hasn't been used for ages
        File should be located within the data folder along GLOBAL.PSFILEPATH, i.e. '/reg/neh/operator/mecopr/mecpython/pulseshaping/data'
        Channels formatted to match output of LeCroy schall function, but this also isn't used much anymore at all
        Only the header of the name needs to be specified, i.e. not counting _ch1, etc.
        """
        ListofListedValues=[]
        for ii in range(1,5):
            with open(GLOBAL.PSFILEPATH+'data/'+CSVHeadname+'_ch'+str(ii)+'.csv','r') as filehead:
                RawListQ=filehead.read()
                ListedValues=RawListQ.split('\n')
            ListofListedValues.append(ListedValues[:-1])
        cls.ll(ListofListedValues)
        return 

    def rcsv(CSVname):
        """
        Shorthand function for reading a list of y-values from a CSV file; hasn't been used for ages
        File should be located within the data folder along GLOBAL.PSFILEPATH, i.e. '/reg/neh/operator/mecopr/mecpython/pulseshaping/data'
        Returns the values read out from the CSV file as an array
        """
        with open(GLOBAL.PSFILEPATH+'data/'+CSVname+'.csv','r') as filehead:
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

    def ll(llist,xlb='none',ylb='none',xlim='none',ylim='none'):
        """
        Shorthand plotting function for a list of multiple lists of y-values
        Example: ep.ll([[1,2,4],[2,7,5],[8,4,2]]) plots three traces over the top of each other:
                 [1,2,4] in one trace, [2,7,5] in one trace, and [8,4,2] in one trace
        """
        df1=plt.figure()
        for ii in range(len(llist)):
            plt.plot(llist[ii]);
        if xlb != 'none':
            plt.xlabel(xlb)
        if ylb != 'none':
            plt.ylabel(ylb)
        if xlim != 'none':
            try:
                plt.xlim(xlim)
            except:
                pass
        if ylim != 'none':
            try:
                plt.ylim(ylim)
            except:
                pass
        df1.show()
        
    def llxy(llistxyq,xlb='none',ylb='none',xlim='none',ylim='none'):
        """
        Shorthand plotting function for a list of multiple lists of x-values and y-values
        Example: ep.llxy([[[1,2,3],[1,2,4]],[[1,2,3],[2,7,5]],[[.9,1.7,2.8],[8,4,2]]]) plots three traces over the top of each other:
                 points (1,1), (2,2), and (3,4) in one trace, points (1,2), (2,7), and (3,5) in one trace, and
                 points (0.9,8), (1.7,4), and (2.8,2) in one trace
        Optional: Set the labels for the x-axis and y-axis with xlb and ylb, respectively
        Optional: Set the limits for the x-axis and y-axis with xlim and ylim, respectively
        """
        df1=plt.figure()
        for ii in range(len(llistxyq)):
            plt.plot(llistxyq[ii][0],llistxyq[ii][1]);
        if xlb != 'none':
            plt.xlabel(xlb)
        if ylb != 'none':
            plt.ylabel(ylb)
        if xlim != 'none':
            try:
                plt.xlim(xlim)
            except:
                pass
        if ylim != 'none':
            try:
                plt.ylim(ylim)
            except:
                pass
        df1.show()
        return

    def llxyloglog(llistxyq,xlb='none',ylb='none',xlim='none',ylim='none'):
        """
        Shorthand plotting function for a list of multiple lists of x-values and y-values with logrithmic axes in both directions
        Example: ep.llxyloglog([[[10,100,1000],[1,20,4]],[[10,100,1000],[.2,70,5]],[[9,170,2800],[80,4,.02]]]) plots three traces over the top of each other:
                 points (10,1), (100,20), and (1000,4) in one trace, points (10,0.2), (100,70), and (1000,5) in one trace, and
                 points (9,80), (170,4), and (2800,0.02) in one trace
        Optional: Set the labels for the x-axis and y-axis with xlb and ylb, respectively
        Optional: Set the limits for the x-axis and y-axis with xlim and ylim, respectively
        """
        df1=plt.figure()
        for ii in range(len(llistxyq)):
            plt.loglog(llistxyq[ii][0],llistxyq[ii][1]);
        if xlb != 'none':
            plt.xlabel(xlb)
        if ylb != 'none':
            plt.ylabel(ylb)
        if xlim != 'none':
            try:
                plt.xlim(xlim)
            except:
                pass
        if ylim != 'none':
            try:
                plt.ylim(ylim)
            except:
                pass
        df1.show()
        return

    @classmethod
    def llt(cls,listq,Map):
        """
        Shorthand function for plotting a list of lists projected according to a map; not used much currently
        """
        formtra=[]
        for ii in range(len(listq)):
            formtra.append(LPL._TraceFormatting(listq[ii],Map[ii],1, AvgRange=1, FWHM=1))
        cls.ll(formtra)
        return

    @classmethod        
    def llcomp(cls,listq,goalq,Map,tMax):
        """Shorthand function for comparing a list of pulse shapes to a targeted pulse shape; not used much currently"""
        formtra=[]
        formtra.append(LPL._TraceFormatting(listq[-1],Map,tMax, AvgRange=1, FWHM=1))
        formtra.append(goalq)
        cls.ll(formtra)
        return
        
    @classmethod
    def lfft(cls,errlistQ,time_stepQ):
        """
        Shorthand function for loglog-plotting the normalized power spectrum of the Fourier transform of a waveform,
        given the temporal spacing of the waveform using time_stepQ in seconds
        Example: for a waveform of energy data sampled every 1min, errlistQ would be the data and time_stepQ would be 60
        Example: for a simulated laser field with point spacing of 0.1fs, errlistQ would be the laser field and time_stepQ would be 1e-16
        """
        #time_step1=0.1
        freqs1=np.fft.fftfreq(np.array(errlistQ).size, time_stepQ)
        idx1=np.argsort(freqs1)
        fftd1=np.fft.fft(errlistQ)
        ps1=np.abs(fftd1/max(fftd1))**2
        cls.lxyloglog(freqs1[idx1],ps1[idx1])
        return [freqs1[idx1],ps1[idx1]]
    
    @classmethod
    def ShotVsGoal(cls,shotDateStr,dateShotIndex,Psns,SSs,zzJQ):
        """
        Shorthand function for plotting the power vs. time of a specified shot of the LPL compared to a specified goal
        - specify the shot you would like to plot using shotDateStr (e.g. '20240116') and dateShotIndex (e.g. 20)
          to plot the index=20 (21st) shot from shot date 2024 January 16
        - specify the goal you would like to plot against using Psns (e.g. [2.25, 2.25, 2.25, 2.25, 12.25]),
          SSs (e.g. [[0.0, 4.5], [4.5, 10.0], [10.0, 17.0], [17.0, 25.0], [25.0, 100.0]]), and zzJQ (e.g. 45)
          to plot the segment durations, segment endpoint heights, and pulse energy of a particular 45J ramp pulse
        """
        shotData = LPL._weichToPowerVsTime(efc.pickleload2(GLOBAL.PSFILEPATH + "s"+str(shotDateStr)+".p")[int(dateShotIndex)])
        goalData = LPL._PGToPowerVsTime(Psns=Psns,SSs=SSs,zzJQ=zzJQ,)
        cls.llxy([shotData,goalData], xlim=[-1, 1+0.25*np.count_nonzero(goalData)], xlb="Time (ns)", ylb="Power (W)",)
        return

        
class efc:
    """
    Extra Function Class for convenient shorthand ways of doing common things like printing in color, accepting keyboard input, interfacing with PVs, etc.
    Typical usage via efc.[command]
    Possible commands include:
        :reloadchk #checks the version of the code being used
        :reloadpkg #reloads the meclas package from file, in case of version update
        :cstr #generates string with color text/background options
        :cprint #prints strings with color text/background options
        :getch #prompts user for a single input character
        :getch_with_TO #getch but with TimeOut
        :input_with_TO #input but with TimeOut
        :pickledump2 #shortcut for saving objects to file
        :pickleload2 #shortcut for reloading objects from file
        :dotsleep #shortcut for printing dots while waiting
        :ssleep #shortcut for standard socket waiting time
        :rPV #shortcut for reading values from PVs
        :wPV #shortcut for writing values to PVs
    potential future work
        - consider consolidating some functions into one
            - e.g. getch with and without TimeOut?
        - add EZeLog function
            - take care of all the hassle of posting to the LCLS eLog from Python so people can do so very easily
        - add threading helper functions?
    """
        
    def reloadchk():
        """
        Shorthand sanity check for the current version of the code
        When making code edits, the author typically administratively writes the date and maybe a unqiue and helpful message
        """
        print('Last stamped: 20250210')
        
    def reloadpkg():
        """
        Shorthand way of reloading the meclas package using the %run IPython command
        Example: after making and saving edits to the file while testing, simply run reloadpkg() to load the latest changes
                 into the current session of hutch python (rather than needing to exit and start a new session)
        """
        from IPython import get_ipython
        #spec = importlib.util.spec_from_file_location("mecps4", "/reg/neh/operator/mecopr/mecpython/pulseshaping/mecps4.py")
        #mecps4 = importlib.util.module_from_spec(spec)
        #spec.loader.exec_module(mecps4);
        #importlib.reload(pkgname)
        ipy = get_ipython()
        ipy.magic("run /reg/g/pcds/pyps/apps/hutch-python/mec/mec/macros/meclas.py")
        #/reg/g/pcds/pyps/apps/hutch-python/mec/mec/macros/

    @staticmethod
    def _tc():
        """An internal list of color codes used for color printing to terminal, used by cprint(), based on ANSI escape codes"""
        #future add: RGB 0-255 foreground: "\033[38;2;" + Rch + ";" + Gch + ";" + Bch + "m"
        #future add: RGB 0-255 background: "\033[48;2;" + Rch + ";" + Gch + ";" + Bch + "m"
        colors=['ENDC','BLINK','K','R','G','Y','B','M','C','W','BK','BR','BG','BY','BB','BM','BC','BW','BRK','BRR','BRG','BRY','BRB','BRM','BRC','BRW','BBRK','BBRR','BBRG','BBRY','BBRB','BBRM','BBRC','BBRW']#B- for background, BR+ for 
        colorcodes=['\033['+str(ii)+'m' for ii in [0,5,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,90,91,92,93,94,95,96,97,100,101,102,103,104,105,106,107]];
        return dict(zip(colors,colorcodes))

    @classmethod
    def cprint(cls,strQ,paramsQ):#delim with commas
        """
        Prints to terminal using provided parameters for color, etc.
        Param color choices are 'K','R','G','Y','B','M','C','W'; also available is 'BLINK'
            :Colors above correspond to blacK, Red, Green, Yellow, Blue, Magenta, Cyan, White
            :Add a 'BR' before each color to specify it to the BRight
            :Add a 'B' on the very front to specify the color as Background
            :Example: 'W' is White (text)
            :Example: 'BK' is Background blacK
            :Example: 'BRR' is BRight Red (text)
            :Example: 'BBRB' is Background BRight Blue
        Multiple parameters can be delimited using commas (order does not matter)
            :Example: 'BLINK,BRY,BB' is BLINKing BRight Yellow text on Background Blue 
        Example: cprint('Warning!','BLINK,BRW,BBRR') prints 'Warning!' in blinking BRight White text on Background BRight Red
        Alternatively, color strings can be created using cstr() and passed straight to the standard print() function
        """
        prargs=''
        if len(paramsQ) == 0:
            paramsQ='ENDC'
        for eaarg in paramsQ.split(','):
            prargs+=cls._tc()[eaarg.upper()]
        print(f"{prargs}"+strQ+f"{cls._tc()['ENDC']}")
        return

    @classmethod
    def cstr(cls,strQ,paramsQ):
        """
        Prepares a string (i.e. for later printing to terminal) using provided parameters for color, etc.
        Param color choices are 'K','R','G','Y','B','M','C','W'; also available is 'BLINK'
            :Colors above correspond to blacK, Red, Green, Yellow, Blue, Magenta, Cyan, White
            :Add a 'BR' before each color to specify it to the BRight
            :Add a 'B' on the very front to specify the color as Background
            :Example: 'W' is White (text)
            :Example: 'BK' is Background blacK
            :Example: 'BRR' is BRight Red (text)
            :Example: 'BBRB' is Background BRight Blue
        Multiple parameters can be delimited using commas (order does not matter)
            :Example: 'BLINK,BRY,BB' is BLINKing BRight Yellow text on Background Blue
        Example: cstr('Warning!','BLINK,BBRY,K') creates a string 'Warning!' that, when printed to terminal
                 using print(), appears in blinking blacK text on Background BRight Yellow
        Note: print(cstr('abc','R,BBRB')) is equivalent to cprint('abc','R,BBRB'), so cstr() is mostly useful
              when auto-generating colored strings to be printed later
        Note: if concaenating colored strings with normal strings, the normal strings will not be affected
        Note: concatenating multiple color strings is necessary for creating combined strings with varying coloration
        """
        prargs=''
        if len(paramsQ) == 0:
            paramsQ='ENDC'
        for eaarg in paramsQ.split(','):
            prargs+=cls._tc()[eaarg.upper()]
        return f"{prargs}"+str(strQ)+f"{cls._tc()['ENDC']}"

    @staticmethod
    def _keybd():
        """"Prepares some keyboard input interpretation parameters needed for interpreting certain keystrokes returned by getch()"""
        return dict(zip(['key_Enter','key_Esc','key_Up','key_Dn','key_Rt','key_Lt'],[13,27,'\033[A','\033[B','\033[C','\033[D']))

    def getch():
        """
        Similar to input() but takes only a single character and doesn't require hitting Enter
        Example: my_char = getch() will cause the terminal to wait for keyboard input and then record
                 the first keypress into my_char and then return to terminal or the next line in a function/script
        """
        fdInput = sys.stdin.fileno()
        termAttr = termios.tcgetattr(0);#fdInput);#0) test to fix print problems
        tty.setraw(fdInput)
        ch = sys.stdin.buffer.raw.read(4).decode(sys.stdin.encoding)
        if len(ch) == 1:
            if ord(ch) < 32 or ord(ch) > 126:
                ch = ord(ch)
        elif ord(ch[0]) == 27:
            ch = '\033' + ch[1:]
        termios.tcsetattr(fdInput, termios.TCSADRAIN, termAttr)
        return ch     

    def getch_with_TO(TOsec,display=True):
        """
        Same as getch() except with a useful timeout period so as to not block a terminal waiting for input
        Use TOsec to set the time-out window in seconds
        Use display=True if you want a message printed to terminal that lets the operator know the time-out duration
        Use display=False to avoid printing a message to terminal telling how many seconds one has to enter a character
        Example: my_char = getch_with_TO(5,display=False) will wait for the next keyboard input for five seconds
               : if a keystroke is not recorded within 5 seconds, function returns False
               : if a keystroke *is* recorded within 5 seconds, function returns the first key pressed
        """
        if display:
            print("You have {} seconds to answer! ".format(str(TOsec)),end='',flush=True)
        fdInput = sys.stdin.fileno()
        termAttr = termios.tcgetattr(0);#fdInput);#0) test to fix print problems
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
            print('Timed out!',end='\r\n',flush=False)
            termios.tcsetattr(fdInput, termios.TCSADRAIN, termAttr)
            return False

    def input_with_TO(TOsec, display=True):
        """
        Similar to input() but includes a user timeout so the window for input doesn't stay open forever and block the terminal
        Use TOsec to set the time-out window in seconds
        Use display=True if you want a message printed to terminal that lets the operator know the time-out duration
        Use display=False to avoid printing a message to terminal telling how many seconds one has to enter input
        Example: my_input = input_with_TO(30,display=False) will wait for keyboard input to be entered for thirty seconds
               : if input is not entered within 30 seconds, the function returns False
               : if input *is* enetered within 30 seconds, the function returns the entered text
        """
        if display:
            print("You have {} second{} to answer! ".format(str(TOsec), '' if TOsec==1 else 's'),end='',flush=True)
        i, o, e = select.select( [sys.stdin], [], [], TOsec)
        if (i):
            return sys.stdin.readline().strip()
        else:
            print('Timed out!',end='\r\n',flush=False)
            return False
    
    def pickledump2(objQ,fullFileNameQ):
        """
        Shortcut for generating pickle files and setting file access permissions as liberally as possible
        The first argument objQ is the python object you want to save to file
        The second argument fullFileNameQ is the full file path and file name for your file
            :note that for pickle files it is best to end the file with the '.p' file extension
        Example: pickledump2(thanos_data, '/reg/neh/operator/mecopr/thanos_data.p') saves the contents of thanos_data
                 to the file '/reg/neh/operator/mecopr/thanos_data.p' and sets file access permissions liberally
        """
        pickle.dump(objQ,open(fullFileNameQ,'wb'));
        os.chmod(fullFileNameQ,stat.S_IRUSR|stat.S_IWUSR|stat.S_IRGRP|stat.S_IWGRP|stat.S_IROTH|stat.S_IWOTH);#
        #os.chmod(fullFileNameQ,stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO);#
        return

    def pickleload2(fullFileNameQ):
        """
        Shortcut for loading pickle files without having to type the open() command
        The only argument fullFileNameQ is the full file path and file name for your file
            :note that pickle files tend to have the '.p' file extension
        Example: thanos_data = pickleload2('/reg/neh/operator/mecopr/thanos_data.p') loads the contents of
                 the file '/reg/neh/operator/mecopr/thanos_data.p' into thanos_data
        """
        return pickle.load(open(fullFileNameQ,'rb'));

    
    def dotsleep(tSEC):
        """
        Similar to time.sleep() but also prints a . character every second for the entire duration, finishing with a * character
        Example: dotsleep(5) prints a '.' character to screen once per second for 5 seconds, at which point it prints a '*'
        Note: printing all happens in-line rather than on separate lines
        """
        for ii in range(tSEC):
            print('.',end='',flush=True);time.sleep(1);
        print('*')
        return

    def ssleep():
        """
        A convenient shortcut for pausing 150ms for certain types of socket functions to wrap up
        ssleep() is equivalent to time.sleep(0.15)
        """
        time.sleep(0.15);
        return

# =============================================================================
#     @classmethod                  
#     def pldemo(cls):
#         """An old demo meant to demonstrate how Python could be used to tune a laser based on camera feedback and EPICS actuators"""
#         plt.ion() 
#         fig,axs=plt.subplots(1,1) 
#         plt.show()
#         xpos=85;ypos=65;xdel=20;ydel=20;
#         Z=[[np.exp(-((ii-xpos)/10)**2-((jj-ypos)/7.5)**2) for ii in range(100)] for jj in range(100)]
#         Zref=[[1 if np.exp(-((ii-50)/2)**2-((jj-50)/2)**2) > 0.5 else 0 for ii in range(100)] for jj in range(100)]
#         ax1=axs.imshow(Z,origin='lower')
#         axs.imshow(Zref,alpha=0.1,origin='lower')
#         #ax1,=axs[0].plot(xdat,ydat)
#         #ax2,=axs[1].plot(xdat,ydat) 
#         #ax4,=axs[1,1].plot(xdat,ydat) 
#         axss=[ax1]#[ax1,ax2,ax3,ax4]
#         cont=True
#         while cont:
#             axss[0].set_data(Z)
#             fig.canvas.draw_idle()
#             plt.pause(0.025)
#             qq=cls.getch()
#             if qq==cls._keybd()['key_Dn']:
#                 ypos-=ydel
#             elif qq==cls._keybd()['key_Up']:
#                 ypos+=ydel
#             elif qq==cls._keybd()['key_Rt']:
#                 xpos+=xdel
#             elif qq==cls._keybd()['key_Lt']:
#                 xpos-=xdel
#             elif qq=='w':
#                 ydel=ydel*2
#             elif qq=='s':
#                 ydel=ydel/2
#             elif qq=='a':
#                 xdel=xdel/2
#             elif qq=='d':
#                 xdel=xdel*2
#             elif qq==cls._keybd()['key_Esc']:
#                 cont=False
#             else:
#                 pass
#             Z=[[np.exp(-((ii-xpos)/10)**2-((jj-ypos)/7.5)**2) for ii in range(100)] for jj in range(100)]
#             print('['+str(xpos)+','+str(ypos)+']')
#         plt.ioff()
# =============================================================================
    
    def rPV(yourPV,display=True):
        """
        Convenient short-hand way to read out and return the value from a PV
        Input yourPV needs to be formatted as a string
        If display=True then any readout failures are printed to terminal
        If display=False then there is no printout (used mainly in report generation while checking many PVs)
        """
        try:
            temppv=EpicsSignal(yourPV);
            tempval=temppv.get()
            return tempval
        except:
            if display:
                print('Read-out failed: {}!'.format(yourPV))
            return False
        
    def wPV(yourPV, yourVal, display=True):
        """
        Convenient short-hand way to write a value to a PV
        Input yourPV needs to be formatted as a string
        Input yourVal needs to fit the type (e.g. string vs number, etc.) expected by yourPV
        """
        try:
            temppv=EpicsSignal(yourPV);
            temppv.put(yourVal)
        except:
            if display:
                print('Write failed:{} and {}!'.format(yourPV, yourVal))
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
#         ppp=0def View(*CAMargs,ImageNo=2,LIVE=False,MAXLOOPS=10):
#     if CAMargs == ('all',):
#         CAMargs = ('Regen', 'Trap', 'StrInA', 'StrInB', 'MPA1In', 'MPA1Out', 'MPA2In', 'MPA2Out', 'MPA2Xtal', 'CompIn', 'CompOutNF', 'CompOutFF')
#     if len(CAMargs) == 1:
#         _QuickView(*CAMargs,LIVE=LIVE,MAXLOOPS=MAXLOOPS)
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

# =============================================================================
#     def reprintdemo():
#     #does not work in IPython on RHEL machines
#         for x in range(10):
#             print('{:>5}'.format(x*10**int(5*np.random.rand())), end='\r');time.sleep(1);
#             print()
# =============================================================================





class HAWG:
    """
    Class containing all the necessary functions for running the Highland Arbitrary Waveform Generator for LPL pulse shaping
    
    Unless speed is necessary, it is usually most appropriate to interface with the Highland simply by using HAWG().[command].
    This will take care of all of the socket opening/closing by itself.
    Example: read out the current Highland waveform using HAWG().ReadPulseHeights()
    Example: reset the Highland using HAWG().Reset()
    (Alternatively, save the initialized object via SH=HAWG() and then use SH.ReadPulseHeights(), etc.)
    List of possible commands includes:
        :ReadStatus
        :ClearStatus
        :ReadPulseHeights
        :WritePulseHeights
        :ReadFiducialImpulseSettings
        :WriteFiducialImpulseSettings
        :WriteEnableByte
        :ReadT0Delay
        :ReadWaveAmplitudeCalibrations
        :WriteWaveAmplitudeCalibrations
        :ReadWaveTimeCalibrations
        :WriteWaveTimeCalibrations
        :ReadMiscellaneousCalibrations
        :WriteMiscellaneousCalibrations
        :ReadWalkTable
        :WriteWalkTable
        :FidOn
        :FidOff
    
    Most of these functions are for expert use only, so please be careful with them!
    Potential future work:
        : FET surveys
        : Highland internal parameter saving and restoring
        : find fiducial bump
        : other outstanding functions not yet programmed from the T400B manual
    """
    def __init__(self):
        """Initializes the object; only one should be instantiated at a time"""
        self._HighlandSocket=None
        self._HIGHLAND_SLAVE_ADDRESS=0 #arbitrary, I think    

    def _Open(self):
        """
        Takes care of opening the socket to the Highland; if called explicitly like this, it MUST be followed by a _Close()
            statement or else you'll block the socket and need to power cycle the unit using HAWG().Reset()
            
        Takes care of opening the socket to the Highland; if called explicitly like this, 
        it MUST be followed by a _Close() statement or else you'll block the socket and need to 
        power cycle the unit using HAWG().Reset()
        
        Using this function allows one to leave the socket open, which allows for quicker access to Highland functions.
            :Example: after SH = HAWG() and SH._Open() then one may use functions inside a loop, e.g.
             for loops in range(50):
                 curr_wvfm = SH._ReadPulseHeights()
                 #calculate some change to the waveform based on feedback from diodes, etc.
                 SH._WritePulseHeights(new_wvfm)
             reads and writes 50 Highland waveforms withotu opening and closing the socket in between each loop
             WARNING: DO NOT FORGET to close the socket at the end of your loop (etc.) using SH._Close()
        In general, there is an "underscored" version of most of the functions mentioned in the SH docstring that
            can be used in the way described above (e.g. _ReadPulseHeights, _GetStatus, etc.).
        """
        if not self._HighlandSocket:
            try:
                self._HighlandSocket=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                self._HighlandSocket.settimeout(1.0)
                self._HighlandSocket.connect((GLOBAL.HIGHLAND_IP, 2000))#'highland-mec-01'
                #Highland's IP address can be changed using the Lantronix DeviceInstaller
            except:
                print('HIGHLAND NOT CONNECTED')
        else:
            print('Socket may already be open!')

    def _Close(self):
        """
        Takes care of closing the socket to the Highland if first explicitly opened using _Open()
        Use this after you have taken care of all of your business that you started with _Open()
        Example: SH=HAWG()
                 SH._Open()
                 #a bunch of consecutive Highland calls using underscored commands like _ReadPulseHeights, _WritePulseHeights, etc.
                 SH._Close()
        """
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
        """        
        A function wrapper that allows one to call each Highland command directly without having to worry about opening/closing sockets;
        if issuing one-off commands that don't require high-frequency consistent execution, this is sufficient
        In short, this automatically wraps an _Open() and _Close() statement around an "underscored" function (like _ReadPulseHeights)
           in order to create a function (like ReadPulseHeights) that can be used without explicitly worrying about good socket habits
        """
        try:
            self._Open();time.sleep(0.15);
            HDataQList=FuncQ(**kwargs);time.sleep(0.15);
            self._Close();time.sleep(0.15);
        except:
            self._Close();time.sleep(0.15);
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
            ep.l(HDataQList)
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
    
    def IndFETWave(self,ListOfPixels,WriteValue):
        """
        Primarily used to write patterns of individual pixels on the Highland for the sake of calibration
        Historically only really used when recalibrating the Highland, which we don't do that often
        """
        itt=0
        NewString=''
        while itt<140:
            if (itt+1) in ListOfPixels:
                NewString+=self._Hex2Byte(WriteValue)
            else:
                NewString+='0000'
            itt+=1
        return NewString
    
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
    
# =============================================================================
#     def HParamSnapshot():
#         """Once revived again, function will serve to save and preserve the various settings of the Highland"""
#         pass
# =============================================================================
    
    def FidOn(self,amplitude=6000,delay_ps=45125):
        """
        Shortcut function for turning on the Highland's fiducial impulse used for proper timing of the oscilloscope delay
        If necessary, non-standard values can be passed into the function using
        amplitude (usu. ~20000, depending on energy saturation in YFE; range 0-65000) and
        delay_ps (usu. 45125 to be in correct location, but can be scanned to help e.g. calibrate a streak camera)
        """
        try:
            self.WriteFiducialImpulseSettings(amplitude,delay_ps);
        except:
            print('Error!')
            return False
        return

    def FidOff(self):
        """
        Shortcut function for turning off the Highland's fiducial impulse
        """
        try:
            self.WriteFiducialImpulseSettings(0,0);
        except:
            print('Error!')
            return False
        return

    

class LOSC:
    """
    Class containing all the necessary functions for running the LeCroy oscilloscopes
    Because we have several such scopes, instantiation of a certain device is required
    Possible scope choices are:
    my_scope=LOSC('A') #for the YFE, 1in1w, and SHG_opt diodes
    my_scope=LOSC('B') #for the four 2in1w diodes
    my_scope=LOSC('1') #for the instruments scientists' oscilloscope
    my_scope=LOSC('2') #for the four 2in2w diodes
    
    Unless speed is necessary, it is usually most appropriate to interface with a LeCroy simply by using LOSC('[1/2/A/B]').[command]
    This will take care of all of the socket opening/closing by itself.
    Example: read out all the current waveform amplitudes on scope '2' using wvfm4=LOSC('2').rchall()
    Example: wait for a fresh acquisition before reading out channel 2's voltage vs time on scope '1' using ch2_wvfm=LOSC('1').waitrchxy(2)
    (Alternatively, use the approach above: my_scope = LOSC('A') and then wvfm4=my_scope.rchall() or my_scope.pchxy(4) or whatever)
    
    Possible commands that can be used as described above include:
        :waitrch #wait and read specified channel amplitude
        :waitrchxy #wait and read specified channel amplitude and time
        :rch #immediately read specified channel amplitude
        :rchxy #immediately read specified channel amplitude and time
        :rchall #immediately read amplitude for all channels
        :rchallxy #immediately read amplitude and time for all channels
        :sch #save plot of specified channel amplitude to file
        :schall #save plot of amplitude for all channels to file
        :pch #plot specified channel amplitude
        :pchxy #plot specified channel amplitude vs time
        :pchall #plot amplitude of all channels
        :pchallxy #plot amplitude vs time of all channels
        :sumch #sum apmlitude of specified channel
        :save_scope_to_eLog #save specified channel amplitude to eLog
        :RestoreConfig #restore acquisition settings according to internal memory
    There are also functions available for use with bare sockets; these tend to start with underscores.
    See the docstrings for more guidance.
    
    Potential future work:
        : adding more functionality from the LeCroy manual using the _ctrl function
    """
    def __init__(self, LStrQ):
        """
        Initialize a LeCroy oscilloscope for use; possible choices are:
        LStrQ='A' #for the YFE, 1in1w, and SHG_opt diodes
        LStrQ='B' #for the four 2in1w diodes
        LStrQ='1' #for the instruments scientists' oscilloscope
        LStrQ='2' #for the four 2in2w diodes
        
        ***If hacking the function to use a non-MEC LeCroy scope, use LStrQ='custom:IP_address' (use at your own risk)***

        """
        if   str(LStrQ).lower() == 'a':
            self._hostIP = GLOBAL.LECROY_A_IP #'172.21.46.100'#'scope-ics-meclas-lecroy-a'
            self._name = 'LeCroyA'
        elif str(LStrQ).lower() == 'b':
            self._hostIP = GLOBAL.LECROY_B_IP #'172.21.46.120'#'scope-ics-meclas-lecroy-b'#
            self._name = 'LeCroyB'
        elif str(LStrQ).lower() == '1':
            self._hostIP = GLOBAL.LECROY_1_IP #'172.21.46.60'#'scope-ics-mectc1-1'
            self._name = 'LeCroy1'
        elif str(LStrQ).lower() == '2':
            self._hostIP = GLOBAL.LECROY_2_IP #'172.21.46.128'#'scope-ics-meclas-lecroy-02'
            self._name = 'LeCroy2'
        elif str(LStrQ).lower() == 'l':
            self._hostIP = GLOBAL.LECROY_L_IP #'172.21.160.252'#'scope-ics-meclas-lecroy-02'
            self._name = 'LeCroyL'
        elif str(LStrQ).lower()[:6] == 'custom':
            self._hostIP = str(str(LStrQ).lower()[7:])
            self._name = 'LeCroy'+str(LStrQ)
        else:
            print('Invalid scope name! Choose 1, 2, A, or B!!')
            return False
        self._LSock = None
        self._port = 1861

    def _Open(self):
        """
        Takes care of opening the socket to the specified LeCroy; if called explicitly like this, 
        it MUST be followed by a _Close() statement or else you'll block the socket and need to 
        locally disable/enable its networking card (or power cycle the unit, but this is not preferred!!)
        
        Using this function allows one to leave the socket open, which allows for quicker access to scope functions.
            :Example: after my_scope = LOSC('1') and my_scope._Open() then one may use functions inside a loop, e.g.
             save_traces = [my_scope._waitrchxy(4) for ii in range(25)] to save 25 consecutive fresh traces without
             opening and closing the socket in between each acquistion;
             WARNING: DO NOT FORGET to close the socket at the end of your loop (etc.) using my_scope._Close()
        In general, there is an "underscored" version of most of the functions mentioned in the LOSC docstring that
            can be used in the way described above (e.g. _waitrchxy, _rchall, _pch, etc.). There are also some
            specialty functions like _ctrl that allow for programming of the oscilloscope using VB commands that
            can be found in the LeCroy manual.
        """
        if not self._LSock:
            try:
                self._LSock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                self._LSock.settimeout(1.0)
                self._LSock.connect((self._hostIP,self._port))
            except:
                print(self._name +' NOT CONNECTED!')
        else:
            print('Socket may already be open!')

    def _Close(self):
        """
        Takes care of closing the socket to the specified LeCroy if preceded by an _Open() statement
        Use this after you have taken care of all of your business that you started with _Open()
        Example: my_scope=LOSC('1')
                 my_scope._Open()
                 #a bunch of consecutive scope calls using underscored commands like _rch, _pchall, etc.
                 my_scope._Close()
        """
        try:
            self._LSock.close()
            self._LSock=None
        except:
            print('Unable to close socket -- it may already be closed')
        return

    def _send_and_reply(self,msg,SendOnly=False):
        """
        Generic utility for sending a poll to the specified LeCroy's internal processor and receiving its reply
        Not typically used in external scripts -- mostly just for creating the functions in this class
        """
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
        """
        A function wrapper that allows one to call each LeCroy command directly without having to worry about opening/closing sockets;
        if issuing one-off commands that don't require high-frequency consistent execution, this is sufficient
        In short, this automatically wraps an _Open() and _Close() statement around an "underscored" function (like _rchall)
           in order to create a function (like rchall) that can be used without explicitly worrying about good socket habits
        """
        try:
            self._Open();time.sleep(0.15);
            LData=FuncQ(**kwargs);time.sleep(0.15);
            self._Close();time.sleep(0.15);
        except:
            self._Close();time.sleep(0.15);
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
        with open(GLOBAL.PSFILEPATH+'data/'+str(FileName)+'.csv','w',newline='') as f:
            writer=csv.writer(f, delimiter='\n')
            writer.writerow(parseddataq['DATA'])
        with open(GLOBAL.PSFILEPATH+'data/'+str(FileName)+'-h.csv','w',newline='') as f:
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
        ep.l(pdata)
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
        ep.lxy(*pdata)
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
        ep.ll(pdata)
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
        ep.llxy(pdata)
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
        """
        Saves specified scope channel (data + voltage vs time plot) to the eLog
        Example: LOSC('1').save_scope_to_eLog(chan_to_eLog=2) saves channel 2 of LeCroy1
                 to the current eLog
        """
        ExpName=LPL.get_curr_exp()
        RunNumber=LPL.get_curr_run()
        print(RunNumber)
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
            mecel = elog.ELog({'experiment':ExpName},user='mecopr',pw=pickle.load(open(GLOBAL.PSFILEPATH+'elogauth.p','rb')))
            #print('Connected to eLog of current MEC experiment, which is: '+ExpName)
        except:
            print('Failed to connect to eLog!')
        try:
            chdata=self.rchallxy();
            timestamp=datetime.now().strftime('%Y%m%d_%H%M%S')
        except:
            print('Failed to read out data!')
        for ii in range(4):
            # save file in raw
            #np.savetxt(fpQ+'lecroy1_ch'+str(ii+1)+'_'+timestamp+'_run'+str(RunNumber)+'.dat',tuple(chdata[ii]))
            # save file in column 
            np.savetxt(fpQ+'lecroy1_ch'+str(ii+1)+'_'+timestamp+'_run'+str(RunNumber)+'.txt',np.c_[chdata[ii][0], chdata[ii][1]])
        ep.lxysav(chdata[chan_to_eLog-1][0],chdata[chan_to_eLog-1][1],fpQ+'lecroy1_ch'+str(chan_to_eLog)+'_'+timestamp+'_run'+str(RunNumber)+'.png',abs_path=True)
        fullmsg=str('Scope trace data for all 4 channels saved to '+fpQ+' with time stamp '+timestamp+' from run '+str(RunNumber)+'. Attached are the data and plot files for channel '+str(chan_to_eLog)+'.')
        #eplxy(chdata[chan_to_eLog-1][0],chdata[chan_to_eLog-1][1])
        try:
            mecel.post(fullmsg,attachments=[fpQ+'lecroy1_ch'+str(chan_to_eLog)+'_'+timestamp+'_run'+str(RunNumber)+'.txt', fpQ+'lecroy1_ch'+str(chan_to_eLog)+'_'+timestamp+'_run'+str(RunNumber)+'.png'], tags=['scope_trace'])
            print('Auto-saved to eLog.') 
        except:
            print('Failed to auto-save to eLog!')
        
    def _ctrl(self,msg,SendOnly=True):
        """
        Once completed, this function will be used for sending specialized commands to the specified LeCroy; this can be found in the manual
        Other old scribbles:
        SetScopeParameters(kwargs) --> scan for different inputs like those listed below; 
           could also specify use case and let it set or RCL appropriately
        (msg='TIME_DIV?',SendOnly=False)
        (msg='*RCL 3',SendOnly=False)
        (msg='TDIV 100E-9',SendOnly=True)
        (msg='C1:VDIV 500E-3',SendOnly=True)
        NOTE: ignore the backslash before the triple quotes below if you are reading the bare code... it is there to avoid closing the docstring
        (msg=r\"""vbs 'app.acquisition.triggermode = "single" ' \""",SendOnly=True)
        (msg=r\"""vbs 'app.acquisition.c1.deskew = 0 ' \""",SendOnly=True)
        (msg=r\"""vbs 'app.acquisition.c2.AverageSweeps = 10 ' \""",SendOnly=True)
        (msg=r\"""vbs 'app.acquisition.c2.EnhancedResType = 2 ' \""",SendOnly=True) #0:None, 1:0.5bits, 2: 1.0bits, ..., 6: 3.0bits
        """
        self._send_and_reply(msg,SendOnly=SendOnly)
    
    def RestoreConfig(self):
        """
        Demonstration of one way to use the _ctrl function by means of _FunctionWrapper
        Function recalls scope configuration from scope internal memory file #1
            (note: administratively this is the memory file where we keep the latest config)
        Example: LOSC('A').RestoreConfig() resets LeCroyA's scope setup according to internal memory file #1
        This is especially useful if scope settings were mistakenly changed and need to be restored
            in order to allow for appropriate pulse shaping performance
        """
        LData=self._FunctionWrapper(self._ctrl,{'msg':'*RCL 1','SendOnly':False});
        return LData





class EMeters:
    """
    Class containing readout functions for energy meters on all MEC laser systems
    Typically used via EMeters.[command]
    
    Possible commands include:
        :LPLInChamber #returns in-chamber energy meter read-outs
        :EG #returns 2in2w energy meter read-outs (scaled by calibration factors)
        :EG1wYFE1in #returns 2in2w energy meter read-outs (scaled by calibration factors)
        :EG1w2in #returns 2in2w energy meter read-outs (scaled by calibration factors)
        :EGall #returns 2in2w energy meter read-outs (scaled by calibration factors)
        :SPLEG #returns energy of SPL energy meter diagnostics
        :gentec_refresh #refreshes settings of LPL Gentec meters
        :E_coeff_refresh #refreshes energy coefficients
        :E_synth_refresh #refreshes synthesized energies written to notepad PVs
        
    See the docstrings for more guidance.
    
    Potential future work:
        : improve SPLEG accuracy currently affected by changing ambient conditions
    """
    def LPLInChamber(printDisplay=False):
        """
        Returns energy values for in-chamber LPL Gentec meters
        printDisplay=True means the values will be printed to terminal in addition to being returned in an array
        printDisplay=False means the values will not be printed to terminal while still being returned in an array
        """
        tvalW=GLOBAL.EGLPLWest.get();
        tvalE=GLOBAL.EGLPLEast.get();
        if printDisplay:
            print('Power meter: WEST: ' + str(tvalW) + ', EAST: ' + str(tvalE) + ', TOTAL: ' + str(tvalW+tvalE))
        return [tvalW, tvalE]

                  
    def EG():
        """
        Gathers all the energy meter readouts for the 2in2w outputs and scales them appropriately
        Returns an array containing an array of the guessed values and an array of the in-chamber measured values
        """
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
        """
        Gathers the energy meter readouts for the YFE and 1in outputs and scales them appropriately
        Returns an array containing the guess for each
        """
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
        """
        Gathers all the energy meter readouts for the 2in1w outputs and scales them appropriately
        Returns an array containing the guessed values of each
        """
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

                  
    def EGall(return_txt=False,chamber_meter_in=False, return_energy_only=False):
        """
        Gathers all the energy meter outputs from EG1eYFE1in(), EG1w2in(), and EG()
        Use return_text=True to return the energy report as a string
        Use return_text=False to return the numerical energies in an array instead
        Use chamber_meter_in=True to add a line in the print-out about the in-chamber LPL energy read-outs
        Use chamber_meter_in=False to suppress a line in the print-out about the in-chamber LPL energy read-outs
        Use return_energy_only=True to return ONLY the numerical energies in an array -- no terminal printing, no message
        Use return_energy_only=False to allow for the terminal printing and the message return (if return_text=True)
        """
        [en1wYFE, en1w1in] = EMeters.EG1wYFE1in()
        [en1wAB, en1wEF, en1wGH, en1wIJ] = EMeters.EG1w2in()[0]
        [en2wAB, en2wEF, en2wGH, en2wIJ] = EMeters.EG()[0][0]
        [enWEST, enEAST]=EMeters.EG()[1][:2]
        [cAB,cEF,cGH,cIJ]=PFN.HeadENB()
        if return_energy_only:
            return [en1wYFE, en1w1in, cAB*en1wAB, cEF*en1wEF, cGH*en1wGH, cIJ*en1wIJ, cAB*en2wAB, cEF*en2wEF, cGH*en2wGH, cIJ*en2wIJ]
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
        Psns=LPL._Psns_get()#pickle.load(open(psfpQ+'Psns.p','rb'))
        SSs=LPL._SSs_get()#pickle.load(open(psfpQ+'SSs.p','rb'))
        strlist.append('Current pulse target is: '+str(list(Psns))+' ns, '+str([list(SSs_sub) for SSs_sub in SSs])+' % of max power.')
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
        """ 
        returns the scaled energy meter read-outs of the shot-pulse laser
        energy array output in order of [regen, TOPAS, MPA1, MPA2]
        (currently not widely used due to sampling percentage variability due to changing ambient conditions)
        """
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
        """ 
        Shorthand way of resetting all the desired input parameters of several different Gentec meters (mostly for 2in beams)
        """
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
        """Resets E_coefficients in notepad PVs based on coefficients contained in GLOBAL"""
        GLOBAL.notepadPVreset()
        #pvlist=['MEC:LAS:FLOAT:'+str(ii) for ii in range(31,41)];
        #inddesclist=['YFE','CD1w','AB1w','EF1w','GH1w','IJ1w','AB2w','EF2w','GH2w','IJ2w']
        #desclist=['E_coeff_'+inddesc for inddesc in inddesclist]
        #valulist=[.3578,0.5971,224.0,177.5,307.4*0.849,113.2,111.0*1.17,187.9*0.860,182.1*0.897,123.5*1.25]
        #for jj in range(len(pvlist)):
        #    temppv1=EpicsSignal(str(pvlist[jj]+'.DESC'));temppv2=EpicsSignal(pvlist[jj]);
        #    temppv1.put(desclist[jj]);temppv2.put(valulist[jj]);

                  
    def E_synth_refresh():
        """Updates the synthetic energy guesses (based on energy readouts) stored in notepad PVs"""
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
        #Add in SPL meters later'
        return




                  
class MBC:
    """
    Class for controlling different functions of the MBC bias control box used with the LPL front-end seed laser
    These functions are meant mostly for laser experts.
    Usage can simply proceed via MBC.[command]
    Potential commands include:
        :ModeCheck()
        :IsSafe()
        :Reset()
    Check docstrings of individual functions for more details
    Potential future improvements:
        - dither parameter control (if possible -- would need expansion of IOC capability)
    """
    def ModeCheck():
        """Returns the current mode of the MBC, either 0 (AUTO) or 1 (MANUAL)"""
        return GLOBAL.MBCmode.get()

    def IsSafe():#re-write checks as individual functions
        """
        Verifies that MBC is operating safely -- powered on, in AUTO/MIN mode, not out of voltage range, no fault detected
        Returns True if all systems nominal, returns False if something was amiss
        """
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

    def Reset():
        """
        Resets the bias control box to make sure it is working properly
        """
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
            MBC.Reset();return
        biaschk=[]
        print('Checking the initial MBC bias level...',end='',flush=True)
        for ii in range(3):
            biaschk.append(GLOBAL.MBCbias.get())
            time.sleep(1);print('..',end='',flush=True);time.sleep(1);print('..',end='',flush=True);
        print('*')
        waitloop=True;loopcnt=0;
        biaschklog=[]
        biaschklog.append(np.sum(np.abs(np.diff(biaschk))))
        while waitloop:
            newchk=np.abs(np.diff(biaschk))
            biaschklog.append(np.sum(newchk))
            if np.sum(newchk) > 3:
                print('MBC bias level unstable... '+str(biaschk),end='',flush=True)
                biaschk=[]
                for ii in range(3):
                    biaschk.append(GLOBAL.MBCbias.get())
                    time.sleep(1);print('..',end='',flush=True);time.sleep(1);print('..',end='',flush=True);
                print('')
                loopcnt+=1
                if ((loopcnt >= 10) and (biaschklog[-1] > biaschklog[-2] + 1)) or (loopcnt >= 20):
                    print('MBC bias level stability fail. Aborting and power-cycling...')
                    GLOBAL.MBCbias.put((np.round(time.time()*1000)%2)*9000*np.sign(biaschk[-1]));time.sleep(1);
                    GLOBAL.MBCpwr.put(2);time.sleep(2);
                    MBC.Reset();return
            else:
                print('MBC bias level stabilized... '+str(biaschk))
                waitloop = False
        return





class YFE:
    """
    Class for organizing functions associated with the YLF Front End (YFE) laser system
    Usage can simply proceed via YFE.[command]
    Potential commands include:
        :OnCheck() #checks whether YFE is on or off
        :CommCheck() #checks whether the eDrives are communicating or not
        :On() #initiates turn-on sequence
        :Off() #initiates shut-off sequence
        :Get() #retrieves eDrive current sepoints and RBVs
        :Set(mmQ,currQ) #changes current setpoint corresponding to certain rod diameters
        :SetAll(IOBool) #turns on or off all eDrive currents without turning off emission
        :Trace() #plots oscilloscope trace of YFE output
    Check docstrings of individual functions for more details
    """
    def OnCheck(display=True):
        """
        Checks whether YFE is turned on or off
        Returns False if turned off; returns True if turned on
        Use display=True to print current YFE status to terminal
        Use display=False to avoid printing current YFE status to terminal
        """
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
    
    def CommCheck(display=True):
        """
        Checks whether eDrives are communicating or not, using the readback on
            the temperature/voltage/pulse width as feedback (typical symptom of comm issues)
        Use display=True to print current comm check to terminal
        Use display=False to avoid printing current comm check to terminal 
        """
        YFEadd='MEC:LPL:LCO:0'
        YFEamp=['2','3','5','6','1','4']
        YFEsuf=[':SensedCurrent',':ActiveCurrent',':PowerSupply',':Temperature',':Emission_RBV',':Emission',':FaultState.RVAL',':ClearFault',':PulseWidth_RBV']
        print('Checking for eDrive comm issues...')
        tempdegCpvlist=[]
        voltageVpvlist=[]
        pulsewidthuspvlist=[]
        for ii in range(len(YFEamp)):
            temprbvemispv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[2])
            voltageVpvlist.append(temprbvemispv.get())#check voltage
            temprbvemispv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[3])
            tempdegCpvlist.append(temprbvemispv.get())#check temperature
            temprbvemispv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[8])
            pulsewidthuspvlist.append(temprbvemispv.get())#check pulse width
        commsbool=True#innocent until proven guilty
        if min(tempdegCpvlist) < 15:
            print('Warning! Either a chiller is VERY cold, or at least one of the eDrives is not communicating!')
            commsbool=False
        if min(voltageVpvlist) < 1:
            print('Warning! Voltage of at least one eDrive is low! Check the comms or the TDKLambda.')
            commsbool=False
        if min(pulsewidthuspvlist) < 350:
            print('Warning! Either a pulse width is WAY too short, or at least one of the eDrives is not communicating!')
            commsbool=False
        if display:
            print('Current status:')
            print('temp (degC): ' + str(tempdegCpvlist))
            print('voltage (V): ' + str(voltageVpvlist))
            print('pulse width (us): ' + str(pulsewidthuspvlist)) 
        return commsbool

    @classmethod
    def On(cls,CtrlChk=True):
        """
        Initiates turn-on procedure for YFE laser system, including several preliminary equipment checks; may take a minute or so
        Use CtrlChk=True to check the LPL control system (pertinent hosts, IOCs, PVs, etc.) as part of the turn-on procedure
        Use CtrlChk=False to avoid check the LPL control system as part of the turn-on procedure
        """
        if YFE.OnCheck(display=False):
            print('YFE emission already enabled.');return True
        if YFE.CommCheck(display=False) == False:
           print('ERROR! eDrive communication error suspected... Check the eDrive EDM screen.'); return False
        if GLOBAL.LPLPCpwr.get() != 1 and GLOBAL.LPLPCpwr.get() != 3:
            GLOBAL.LPLPCpwr.put(1)
        if GLOBAL.LPLVACpwr.get() != 1 and GLOBAL.LPLPCpwr.get() != 3:
            GLOBAL.LPLVACpwr.put(1)
        if GLOBAL.LPLPS1pwr.get() != 1 and GLOBAL.LPLPCpwr.get() != 3:
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
            time.sleep(5)#3 sec seems to be too short sometimes for it to clear
            faultstatlist=[faultpv.get() for faultpv in faultpvlist]
            if any(faultstatlist):
                print('YFE fault state still detected, turn-on failed.')
                return False
            else:
                print('Fault cleared!')
        if CtrlChk == True:
            print('Checking LPL controls status (~10sec)!')
            PVsuccess = CtrlSys.pv_checker(pv='lpl')
            if not PVsuccess:
                print('Control error detected! Continue with system turn-on? [y/n]')
                checkprompt = efc.getch_with_TO(TOsec=10,display=False)
                if checkprompt not in ('y','Y'):
                    print('Try again later then!');
                    return False
                else:
                    print('OK, I hope you know what you\'re doing!')
        if not MBC.IsSafe():
            cls.SetAll(False,displayQ=False)
            print('MBC not configured properly!')
            MBC.Reset()
            tempbool=cls.On(CtrlChk=False);
            return tempbool #avoid re-running the bottom text in true or false case
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
        if GLOBAL.LPLPCpwr.get() != 1 and GLOBAL.LPLPCpwr.get() != 3:
            print('Failed to turn on Pockels cell!')
        if GLOBAL.LPLVACpwr.get() != 1 and GLOBAL.LPLPCpwr.get() != 3:
            print('Failed to turn on scroll pump!')
        if GLOBAL.LPLPS1pwr.get() != 1 and GLOBAL.LPLPCpwr.get() != 3:
            print('Failed to turn on YFE PS1!')
        return True

    def Off():
        """
        Initiates turn-off procedure for YFE laser system
        """
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
        """
        Retrieves eDrive current sepoints and RBVs for all six diode-pumped heads within the YFE
        Use display=True to print the requested and actual eDrive currents to terminal
        Use display=False to avoid printing the requested and actual eDrive currents to terminal
        """
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
        """
        Changes current setpoint corresponding to certain rod diameters mmQ
        Choices for rod diameters mmQ are 2, 6, or 10
        Current limits for currQ are 88, 135, and 140, respectively
        Nominal values for currQ are 85, 130, and 124, respectively
        Use display=True to print message to terminal of what values were changed
        Use display=False to avoid writing such a message to terminal
        Example: Set(10,100) changes the 10mm head to 100A
        """
        YFEadd='MEC:LPL:LCO:0'
        YFEamp=['2','3','5','6','1','4']
        YFEsuf=[':SensedCurrent',':ActiveCurrent',':PowerSupply',':Temperature',':Emission_RBV',':Emission',':FaultState.RVAL']
        changelist=[]
        if mmQ==2:
            changelist=range(4)
            tempoldcurr02pv=EpicsSignal(YFEadd+YFEamp[changelist[0]]+YFEsuf[1])
            oldcurrQ=str(tempoldcurr02pv.get())
            if currQ>GLOBAL.YFE02MM_MAX:
                print('Too high!')
                return
        elif mmQ==6:
            changelist=[4]
            tempoldcurr06pv=EpicsSignal(YFEadd+YFEamp[changelist[0]]+YFEsuf[1])
            oldcurrQ=str(tempoldcurr06pv.get())
            if currQ>GLOBAL.YFE06MM_MAX:
                print('Too high!')
                return
        elif mmQ==10:
            changelist=[5]
            tempoldcurr10pv=EpicsSignal(YFEadd+YFEamp[changelist[0]]+YFEsuf[1])
            oldcurrQ=str(tempoldcurr10pv.get())
            if currQ>GLOBAL.YFE10MM_MAX:
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
        """
        Shorthand way of adjusting all eDrive currents to nominal setpoints all together
        Use IOBool=True to turn all eDrives up to their nominal active setpoint
            :(reminder: 2mm --> 85A, 6mm --> 130A, 10mm --> 124A)
        Use IOBool=False to turn all eDrives down to zero without turning off the emission
            :(very convenient when resetting the MBC, for example; takes much less time than turning off and on eDrive emission)
        Use displayQ=True to print out current setpoint and RBV after pushing changes to setpoints
        Use displayQ=False to avoid printing out any such message to terminal
        """
        YFEadd='MEC:LPL:LCO:0'
        YFEamp=['2','3','5','6','1','4']
        YFEsuf=[':SensedCurrent',':ActiveCurrent',':PowerSupply',':Temperature',':Emission_RBV',':Emission',':FaultState.RVAL']
        YFElvl=[GLOBAL.YFE02MM_SETPT,GLOBAL.YFE02MM_SETPT,GLOBAL.YFE02MM_SETPT,GLOBAL.YFE02MM_SETPT,GLOBAL.YFE06MM_SETPT,GLOBAL.YFE10MM_SETPT];
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

    def Trace(goal_display=False):
        """
        Plots the voltages of the diode trace for the YFE output waveform
        """
        try:
            LOSC('a').pch(1)
        except:
            print('Failed to display trace')
            return False





class PFN:
    """
    Class for controlling the PFNs of the MEC Nd:glass LPL
    Typical usage is as follows: PFN.[command]
    Possible commands include:
        :HeadENB #reads back which PFNs are currently enabled
        :EnableOnly #enables PFNs of specified arms and disables PFNs of unmentioned arms
        :ArmOnly #same as EnableOnly but also allows specification of the HWP (throttle) on the enabled arms
    """
    def HeadENB():
        """returns array of length 4, with each entry corresponding to charge state of one pair of heads
        [ABstate, EFstate, GHstate, IJstate]
        0 returned if either head in the pair of heads is disabled
        1 returned if both heads in the pair of heads are enabled
        Example: HeadENB() returns [1,1,1,1] if AB, EF, GH, and IJ PFNs are all enabled
        Example: HeadENB() returns [0,0,1,1] if only GH and IJ (east side) PFNs are enabled
        """
        cAB=GLOBAL.PFNAEN.get()*GLOBAL.PFNBEN.get()
        cEF=GLOBAL.PFNEEN.get()*GLOBAL.PFNFEN.get()
        cGH=GLOBAL.PFNGEN.get()*GLOBAL.PFNHEN.get()
        cIJ=GLOBAL.PFNIEN.get()*GLOBAL.PFNJEN.get()
        return np.array([cAB,cEF,cGH,cIJ])

    def EnableOnly(ArmStrQ):
        """
        enables just the heads listed while disabling all the heads that are not listed
        Example: EnableOnly('all') enables AB, EF, GH, and IJ PFNs
        Example: EnableOnly('ABIJ') enables AB and IJ PFNs and disables EF and GH PFNs
        """
        AllStrQ='ABEFGHIJ'
        if ArmStrQ.lower() == 'all':
            ArmStrQ = AllStrQ
        ArmStrQ = ArmStrQ.upper()
        GLOBAL.PFNmode.put(0)
        time.sleep(2)
        for ii in range(len(AllStrQ)):
            if (AllStrQ[ii] in ArmStrQ):
                temppv=EpicsSignal(str('MEC:PFN:CH'+str(ii+1)+':ENABLE'))
                temppv.put(1)
            else:
                temppv=EpicsSignal(str('MEC:PFN:CH'+str(ii+1)+':ENABLE'))
                temppv.put(0)
        time.sleep(2);GLOBAL.PFNmode.put(1);time.sleep(3.5);GLOBAL.PFNmode.put(2);
        return

    @classmethod
    def ArmOnly(cls,ArmStrQ,set_T=1):
        """
        enables just the heads listed AND changes HWP settings of specified arms to set_T;
        all the heads that are not listed are disabled (although their HWP values are not changed)
        Example: ArmOnly('all') enables all PFNs and sets all HWPs for full transmission
        Example: ArmOnly('EFGH',set_T=0.5) enables EF and GH PFNs, disables AB and IJ PFNs,
                 and sets EF and GH HWPs for 50% transmission (AB and IJ HWPs are not touched)
        """
        HWP.On(ArmStrQ,set_T=set_T)
        cls.EnableOnly(ArmStrQ)
        return
        
    
    
    
        
class HWP:
    """
    Class for controlling the HWPs of the MEC Nd:glass LPL
    Typical usage is as follows: HWP.[command]
    Possible commands include:
        :On #sets the desired transmission level of the specified arms
        :Status #prints current HWP settings and returns an array of anticipated transmission levels
        :ClearStart #attempts to clear motor restart for AB and EF HWPs
        :HWP_opt #attempts to calibrate the HWP settings such that 0deg corresponds to full transmission
    """
    def On(ArmStrQ,set_T=1):
        """
        Adjusts the waveplates of the specified arms in order to provide the specified transmission
        Example: On('all') sets the transmission for AB, EF, GH, and IJ to 100% (default for set_T is 1)
        Example: On('ABIJ',set_T=0.75) sets the transmission for AB and IJ to 75%
        """
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
    
    def Status():
        """
        prints current HWP settings to terminal and returns an array of anticipated transmission levels
        Example: if all HWPs are set to 0deg, returned array has form [1,1,1,1] showing 100% transmission on all arms
        Example: if all AB and IJ HWPs are set to 22.5deg, and EF and GH HWPs are set to 45deg,
                 returned array has form [0.5,0,0,0.5] showing AB and IJ HWPs are set for 50% transmission,
                 and EF and GH are set for 0% transmission
        (Both examples will also include a print-out to terminal showing current angle setting and corresponding transmission)
        """
        armlist=['AB','EF','GH','IJ'];
        HWPpvlist=[GLOBAL.HWPAB,GLOBAL.HWPEF,GLOBAL.HWPGH,GLOBAL.HWPIJ];
        set_Tlist=[round((np.cos(2*motor_angleRBV.get()*(np.pi/180)))**2,4) for motor_angleRBV in HWPpvlist]
        print('Current waveplate settings:');
        for ii in range(4):
            print(armlist[ii]+': '+str(round(HWPpvlist[ii].get(),1))+'deg --> T~'+str(round(100*set_Tlist[ii],1))+'%')
        return set_Tlist
        
    def ClearStart():
        """
        attempts to clear motor restart for AB and EF HWPs (running on MFORCE chassis)
        (no equivalent function for GH and IJ waveplates yet (running on Newport crate))
        """
        resetPVs=[GLOBAL.HWPABclr, GLOBAL.HWPEFclr]#no equivalent for GH and IJ yet...
        try:
            for eapv in resetPVs:
                eapv.put(1)
        except:
            print('Failed!')
        
    def HWP_opt(armsQ='ABEFGHIJ'):#check for trace height;#All shutters must start in the open state... 
        """
        Like SHG_opt() but instead optimizes the HWP angles to set maximum transmission at 0deg
        (HWP_opt is really only useful if someone manually rotates the knob inside the enclosure, so usage is very infrequent)
        Prompts included to help insure system is in appropriate state for optimization
        Entire function takes several minutes (longer than SHG_opt); live display allows you to monitor progress
        If new value is too close to edge of window, it is recommended to optimize that arm again
        (If the angles are way far off (no detection), tune the Newport motor back until it is closer)
        
        If leaving armsQ blank: all arms are optimized in the order of AB, EF, GH, IJ
            :EXAMPLE: HWP_opt() is equivalent to HWP_opt('ABEFGHIJ') and optimizes all arms
        Instead, one may choose specific arms to optimize
            :EXAMPLE: HWP_opt(armsQ='ABIJ') or HWP_opt('ABIJ') optimizes only arms AB and IJ
            :EXAMPLE: HWP_opt(armsQ='EF') or HWP_opt('EF') optimizes only arm EF
        """
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
            print('(Warning! The MBC doesn\'t appear to be in AUTO mode!')
        else:
            print('(MBC mode seems OK...)')
        print('Are you sure you are ready to proceed? [enter y/n]',end='',flush=True)
        checkprompt=efc.getch_with_TO(TOsec=10,display=False);
        if checkprompt not in ('y','Y'):
            print('Try again later then!');
            return
        else:
            print('OK, I hope you know what you\'re doing!')
        HWPpvlist=[GLOBAL.HWPAB, GLOBAL.HWPEF, GLOBAL.HWPGH, GLOBAL.HWPIJ];
        HWP.On('all',set_T=1)
        currHWPrbv=np.sum([np.abs(HWP_rbv.get()) for HWP_rbv in HWPpvlist])
        if  currHWPrbv > 10:
            hwploopiter=0
            while currHWPrbv > 1:
                time.sleep(1)
                currHWPrbv=np.sum([np.abs(HWP_rbv.get()) for HWP_rbv in HWPpvlist])
                hwploopiter+=1
                if hwploopiter > 10:
                    print('HWP settings do not seem to be settling at 0deg! Please try again later!')
                    return False
        armlist=['AB','EF','GH','IJ']
        #YFEoff();YFEon();
        GLOBAL.MBCmode.put(1)#set MAN mode on MBC
        if np.sum(YFE.Get(display=False)) < 100:
            print('Check YFE before optimizing!')
        optwvfm=pickle.load(open(GLOBAL.PSFILEPATH+'opttrace.p','rb'));
        try:
            oldwvfm=HAWG().ReadPulseHeights();
            HAWG().WritePulseHeights(optwvfm);
        except:
            print('Failed to read old waveform/load new waveform!')
            return

        print('Closing all shutters...')
        TTL_shutter.Toggle('closeall',display=False);#close all the shutters
        time.sleep(4)
        GLOBAL.EVRLPLSSEC.put(43);GLOBAL.EVRLPLSSEN.put(1);#enable these...
        try:
            tempchk1=LOSC('a').rch(1);time.sleep(.15);tempchk2=LOSC('a').rch(1);
            if np.sum(np.abs(tempchk1-tempchk2))<1e-6:
                print('Warning: scope trace doesn\'t appear to be updating, please check scope! Abort? [enter y/n]')
                checkprompt=efc.getch_with_TO(TOsec=10,display=False);
                if checkprompt not in ('y','Y'):
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
        stepQ=1.0;rangeQ=20.0;
        for ii in range(4):
            if armlist[ii] in armsQ:#only prep the stage if it's going to be used
                HWPpvlist[ii].put(startposlist[ii]+(-rangeQ+stepQ*0))
        currentshutter=0;#trying to re-open a shutter in case of failure...
        #set up all the plotting stuff
        plt.ion()
        fig,axs=plt.subplots(2,2,gridspec_kw={'hspace':0.4,'wspace':0.3})
        xdat=[[startposlist[ii]+(-rangeQ+stepQ*(jj)) for jj in range(int(1+(2*rangeQ/stepQ)))] for ii in range(4)]
        ydat=[[0]*int(1+(2*rangeQ/stepQ)) for ii in range(4)]
        ax1,=axs[0,0].plot(xdat[0],ydat[0]); axs[0,0].set_xlabel('AB'); plt.pause(0.01);
        ax2,=axs[0,1].plot(xdat[1],ydat[1]); axs[0,1].set_xlabel('EF'); plt.pause(0.01);
        ax3,=axs[1,0].plot(xdat[2],ydat[2]); axs[1,0].set_xlabel('GH'); plt.pause(0.01);
        ax4,=axs[1,1].plot(xdat[3],ydat[3]); axs[1,1].set_xlabel('IJ'); plt.pause(0.01);
        axss=[ax1,ax2,ax3,ax4]
        try:
            SLA=LOSC('A');SLA._Open();#changed to LecroyA since repair
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
                    #ep.llxy([[hwparmdatax,hwparmdatay],[xpq,qfitp(xpq)]],xlb=armlist[ii])
                else:
                    print('Skipping '+armlist[ii]+'...')
                    pass
            SLA._Close();time.sleep(.15);#changed to LeCroyA
        except:
            print('Failed! Restoring original values and attempting to re-open most-recent shutter... you should verify!')
            SLA._Close();time.sleep(.15);#changed to LeCroyA
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
        HWPlog=pickle.load(open(GLOBAL.PSFILEPATH+'HWP_opt_log.p','rb'))
        HWPlog.append([datestamp,[newposlist[ii] for ii in range(4)]])
        efc.pickledump2(HWPlog,GLOBAL.PSFILEPATH+'HWP_opt_log.p')
        TTL_shutter.Toggle('openall',display=False);#open all the shutters
        MBC.Reset();YFE.SetAll(True);#reset bias...
        plt.ioff()
        motnamelist=[GLOBAL.HWPABoff, GLOBAL.HWPEFoff, GLOBAL.HWPGHoff, GLOBAL.HWPIJoff]#.OFF
        #add adjustment to offset automatically....
        for temppv in motnamelist:
            tempval=temppv.get();
            temppv.put(tempval-newposlist[ii]);





class Stage:
    """
    Stores functions related to stages and motors
    Current functions include:
        :NewportInitRefAll #useful routine for recovering from Newport outage -- initializes and references all channels
        :NewportBrowserCtrl #internal function meant for launching the Newport broswer controller
    Potential future work:
        -SmarAct motor utilities
        -motor setting snapshots
        -motor setting restoration
    """
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
    
    def NewportBrowserCtrl(inpvnam):
        """
        internal function meant for launching the Newport broswer controller
        takes a single Newport motor PV as input, opens a firefox session of the crate with the corresponding XPS driver
        Example: Stage.NewportBrowserCtrl('MEC:LAS:MMN:25.RBV') opens the admin controls for mcn-mec-las4 in a firefox browser
        """
        
        pvno=re.findall('^MEC:(LAS|PPL):MMN:([\d]+)',inpvnam);
        if len(pvno)>0:
            try:
                pvnoint=int(pvno[0][1])
                #print('PV number int: {}, attempting to read {}'.format(pvnoint,'MEC:LAS:MMN_{:02}{:02}.CTRL'.format((8*(pvnoint//8))+1,8*((pvnoint//8)+1))))
                ipaddr=efc.rPV('MEC:{}:MMN_{:02}{:02}.CTRL'.format(pvno[0][0],(8*(pvnoint//8))+1,8*((pvnoint//8)+1)))
                #print('IP address shorthand: {}'.format(ipaddr))
                rawipaddr=re.findall('IP: (.+)\n', os.popen('netconfig search '+ipaddr).read());
                #print('Raw IP address shorthand: {}'.format(rawipaddr[0]))
                #print('firefox --new-window http://{}'.format(rawipaddr[0]))
                os.system('firefox --new-window http://{}'.format(rawipaddr[0]))
                
                print('Reminder: username/password are Administrator/Administrator')
            except:
                print('Process failed!')
                return False
        else:
            print('Bad PV, please try again!')


                  



                  
class Timing:
    """
    Potential future class containing timing-related utilities such as:
        - EVR snapshot and refresh functions
        - fstiming utilities
        - Nstiming utilities
        - Vitara utilities      
    """
    pass



class CAM:
    """
    Class mostly containing utilities for GigE cameras for MEC
    Typical usage via CAM.[command]
    List of possible commands includes:
        :Name #helps with names and naming conventions of all cameras in MEC
        :View #quickly view a single camera or multiple cameras in a python window
        :QuickSave #saves the data from a specified GigE camera to a PNG image
        :QuickSaveData #saves the data from a specified GigE camera to a 2D txt file
        :Config #configures plug-ins for a given camera
        :ConfigReset #refreshes all current camera configurations
    Potential future work:
        - save images to file easily 
        - (also in connection with scan parameters)
        - GigE_toggle_trigger for configuring SPL camera triggers
        - setdynxhair for quickly setting a custom temporary crosshair position, etc.
        - combined utility for tuning a beam while watching a live image of the camera
    """
    def Name(GIGEnam='none',returnAll=False):
        """
        There are three types of camera names: the PV Name, the SPL CAM Name, and the NickName
        This function helps display the table of names and also translates NickNames and SPL CAM Names
            into the PV Names, which are needed to run the different camera viewer utilities (etc.)
            
        With no argument given, Name() returns a list of names of all cameras
        Example: Name() prints out a table with all supported GigE cameras -- something like this:
            Camera NickName                CamViewer             PV Name      Motor
            Legend     or   Regen          Regen                 MEC:GIGE:24  SM0
            StrInA     or   TopasA         StrinA                MEC:GIGE:22  SM1
            StrInB     or   TopasB         StrinB                MEC:GIGE:28  SM2
            MPA1In     or   MPA1A          MPA1In                MEC:GIGE:23  SM3
            MPA1Out    or   MPA1B          MPA1Out               MEC:GIGE:25  SM4
            MPA2In     or   MPA2A          MPA2In                MEC:GIGE:26  SM5
            MPA2Xtal   or   MPA2F          GigE9                 MEC:GIGE:09  SM6     
            MPA2Out    or   MPA2B          MPA2Out               MEC:GIGE:27  SM7
            CompInNF   or   CompANF        CompInNF              MEC:GIGE:33  SM7
            CompInFF   or   CompAFF        CompInFF              MEC:GIGE:32  SM8
            CompOutNF  or   CompBNF        CompOutNF             MEC:GIGE:31  XPS?
            CompOutFF  or   CompBFF        CompOutFF             MEC:GIGE:30  XPS?
            Trap       or   Mousetrap      MEC_SPL_8             MEC:GIGE:16  None
            TTIn       or   TTA            TTIn                  MEC:GIGE:34  None
            TTOut      or   TTB            TTOut                 MEC:GIGE:35  None
        
        For every PV name (e.g. 'MEC:GIGE:24'), there is a corresponding SPL CAM Name (e.g. 'MEC_SPL_3')
            and two Camera NickNames (e.g. 'Legend' and 'Regen'); there is also a motor that is
            used to tune the beam to the crosshair shown on the camera's image (e.g. SmarAct SM0)
        
        With a single argument for GIGEnam given, function returns the corresponding PV name
            :EXAMPLE: Name('Legend') and Name('Regen') and Name('MEC_SPL_3') all return 'MEC:GIGE:24'
            
        ***If hacking the function to use a non-MEC laser camera, use GIGEnam='custom:NAME:OF:GIGE:PV:HEAD'***
            
        Option: using returnAll=True will return not only the PV Name but also the first Camera NickName
                and the SPL CAM Name in an array of format [NickName_0, SPL CAM Name, PV Name]
              : using returnAll=False is the typical usage that returns just the PV Name
        """
        NickNameList=[['Legend','Regen'], ['StrInA','TopasA'], ['StrInB','TopasB'], ['MPA1In','MPA1A'], ['MPA1Out','MPA1B'], ['MPA2In','MPA2A'], ['MPA2Xtal','MPA2F'],
                      ['MPA2Out','MPA2B'], ['CompInNF', 'CompANF'], ['CompInFF', 'CompAFF'], ['CompOutNF', 'CompBNF'], ['CompOutFF', 'CompBFF'], ['Trap', 'Mousetrap'],
                      ['TTIn','TTA'], ['TTOut','TTB']]
        SPLNameList=['Regen', 'StrinA', 'StrinB', 'MPA1In', 'MPA1Out', 'MPA2In', 'GigE9',
                     'MPA2Out', 'CompInNF', 'CompInFF', 'CompOutNF', 'CompOutFF',
                     'MEC_SPL_8', 'TTIn', 'TTOut']
        PVNameList=['MEC:GIGE:24', 'MEC:GIGE:22', 'MEC:GIGE:28', 'MEC:GIGE:23', 'MEC:GIGE:25', 'MEC:GIGE:26', 'MEC:GIGE:09',
                    'MEC:GIGE:27', 'MEC:GIGE:33', 'MEC:GIGE:32', 'MEC:GIGE:31', 'MEC:GIGE:30',
                    'MEC:GIGE:16', 'MEC:GIGE:34', "MEC:GIGE:35"]
        SmarActList=['SM0','SM1','SM2','SM3','SM4','SM5','SM6','SM7','SM7','SM8','XPS?',
                     'XPS?','None','None','None']
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
            if GIGEnam[:7]=='custom:':
                #print("Using custom camera name: "+GIGEnam[7:])
                if returnAll==True:
                    return [GIGEnam[7:],'','']
                else:
                    return GIGEnam[7:]
            else:
                print('{:<30} {:<21} {:<11}  {:<10}'.format('Camera NickName', 'CamViewer   ', 'PV Name', 'Motor'))
                for ii in range(len(NickNameList)):
                    print('{:<10} or   {:<12}   {:<21} {:<11}  {:<10}'.format(str(NickNameList[ii][0]), str(NickNameList[ii][1]), str(SPLNameList[ii]), str(PVNameList[ii]), str(SmarActList[ii])))
                return False

    def _NameOLD(GIGEnam='none',returnAll=False):
        """
        There are three types of camera names: the PV Name, the SPL CAM Name, and the NickName
        This function helps display the table of names and also translates NickNames and SPL CAM Names
            into the PV Names, which are needed to run the different camera viewer utilities (etc.)
            
        With no argument given, Name() returns a list of names of all cameras
        Example: Name() prints out a table with all supported GigE cameras -- something like this:
            Camera NickName                SPL CAM Name          PV Name      Motor
            Legend     or   Regen          MEC_SPL_3             MEC:GIGE:24  SM0
            StrInA     or   TopasA         MEC_SPL_1             MEC:GIGE:22  SM1
            StrInB     or   TopasB         MEC_SPL_7             MEC:GIGE:28  SM2
            MPA1In     or   MPA1A          MEC_SPL_2             MEC:GIGE:23  SM3
            MPA1Out    or   MPA1B          MEC_SPL_4             MEC:GIGE:25  SM4
            MPA2In     or   MPA2A          MEC_SPL_5             MEC:GIGE:26  SM5
            MPA2Xtal   or   MPA2F          GigE17_TimeTool_Diag  MEC:GIGE:17  SM6     NOW GigE9 MEC:GIGE:09
            MPA2Out    or   MPA2B          MEC_SPL_6             MEC:GIGE:27  SM7
            CompIn     or   CompA          MEC_SPL_9             MEC:GIGE:29  SM8     NOW NO LONGER EXISTS.
            CompOutFF  or   CompBFF        MEC_SPL_10            MEC:GIGE:30  XPS Mid
            CompOutNF  or   CompBNF        MEC_SPL_11            MEC:GIGE:31  XPS In
            Trap       or   Mousetrap      MEC_SPL_8             MEC:GIGE:16  None    look for MEC_SPL_8
            SPFloat1   or   SPFloater1     MEC_SPL_12            MEC:GIGE:32          NOW CompInFF
            CompInNF   or   CompANF        MEC_SPL_13            MEC:GIGE:33
            TTIn       or   TTA            MEC_SPL_14            MEC:GIGE:34
            TTOut      or   TTB            MEC_SPL_15            MEC:GIGE:35
        
        For every PV name (e.g. 'MEC:GIGE:24'), there is a corresponding SPL CAM Name (e.g. 'MEC_SPL_3')
            and two Camera NickNames (e.g. 'Legend' and 'Regen'); there is also a motor that is
            used to tune the beam to the crosshair shown on the camera's image (e.g. SmarAct SM0)
        
        With a single argument for GIGEnam given, function returns the corresponding PV name
            :EXAMPLE: Name('Legend') and Name('Regen') and Name('MEC_SPL_3') all return 'MEC:GIGE:24'
            
        ***If hacking the function to use a non-MEC laser camera, use GIGEnam='custom:NAME:OF:GIGE:PV:HEAD'***
            
        Option: using returnAll=True will return not only the PV Name but also the first Camera NickName
                and the SPL CAM Name in an array of format [NickName_0, SPL CAM Name, PV Name]
              : using returnAll=False is the typical usage that returns just the PV Name
        """
        NickNameList=[['Legend','Regen'], ['StrInA','TopasA'], ['StrInB','TopasB'], ['MPA1In','MPA1A'], ['MPA1Out','MPA1B'], ['MPA2In','MPA2A'], ['MPA2Xtal','MPA2F'],
                      ['MPA2Out','MPA2B'], ['CompIn', 'CompA'], ['CompOutFF', 'CompBFF'], ['CompOutNF', 'CompBNF'], ['Trap', 'Mousetrap'],
                      ['SPFloat1','SPFloater1'], ['CompInNF','CompANF'], ['TTIn','TTA'], ['TTOut','TTB']]
        SPLNameList=['MEC_SPL_3', 'MEC_SPL_1', 'MEC_SPL_7', 'MEC_SPL_2', 'MEC_SPL_4', 'MEC_SPL_5', 'GigE17_TimeTool_Diag',
                     'MEC_SPL_6', 'MEC_SPL_9', 'MEC_SPL_10', 'MEC_SPL_11', 'MEC_SPL_8',
                     'MEC_SPL_12', 'MEC_SPL_13', 'MEC_SPL_14', 'MEC_SPL_15']
        PVNameList=['MEC:GIGE:24', 'MEC:GIGE:22', 'MEC:GIGE:28', 'MEC:GIGE:23', 'MEC:GIGE:25', 'MEC:GIGE:26', 'MEC:GIGE:17',
                    'MEC:GIGE:27', 'MEC:GIGE:29', 'MEC:GIGE:30', 'MEC:GIGE:31', 'MEC:GIGE:16',
                    'MEC:GIGE:32', 'MEC:GIGE:33', 'MEC:GIGE:34', "MEC:GIGE:35"]
        SmarActList=['SM0','SM1','SM2','SM3','SM4','SM5','SM6','SM7','SM8','XPS Mid','XPS In','None',
                     '(placeholder)', '(placeholder)', '(placeholder)', '(placeholder)']
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
            if GIGEnam[:7]=='custom:':
                #print("Using custom camera name: "+GIGEnam[7:])
                if returnAll==True:
                    return [GIGEnam[7:],'','']
                else:
                    return GIGEnam[7:]
            else:
                print('{:<30} {:<21} {:<11}  {:<10}'.format('Camera NickName', 'SPL CAM Name', 'PV Name', 'Motor'))
                for ii in range(len(NickNameList)):
                    print('{:<10} or   {:<12}   {:<21} {:<11}  {:<10}'.format(str(NickNameList[ii][0]), str(NickNameList[ii][1]), str(SPLNameList[ii]), str(PVNameList[ii]), str(SmarActList[ii])))
                return False


    @classmethod
    def Jitter(cls,CAMreq,AcqNum=100,RepRateHz=-1,BeamSig=1):
        """
        Utility for taking a quick beam pointing stability measurement on a GigE camera.
        Prints several statistics and generates plot for centroid jitter relative to the 
            average position, normalized by the spatial beam dimensions sigma_X and sigma_Y

        Parameters
        ----------
        CAMreq : supply the name of the camera using one of the options in CAM.Name(), 
            e.g. 'Regen' or 'CompOutFF'
        AcqNum : number of acquisitions in the scan; the default is 100
        RepRateHz : anticipated repetition rate of the GigE camera in Hz, e.g. RepRateHz=5 for 5Hz operation;
            the default is -1, which instead looks up the current refresh rate using the ArrayRate_RBV PV
        BeamSig : allows custom multiplicative definition of beam sigma relative to the standard
            deviation, with the default of 1; e.g. BeamSig=2.355 would calculate sigma as the FWHM instead

        Returns
        -------
        [centoid_X, centroid_Y, sigma_X, sigma_Y]

        """
        PVhead=cls.Name(CAMreq)
        efc.wPV('{}:Acquire'.format(PVhead),1)#tries to start the camera first
        if RepRateHz<=0:
            RepRateHz=efc.rPV('{}:ArrayRate_RBV'.format(PVhead))#tries to start the camera first
        cX=[]; cY=[]; sX=[]; sY=[];
        for ii in range(AcqNum):
            if (ii+1)%int(AcqNum//5) == 0:
                print('Shot number: {}, Percent complete: {}'.format(str(ii+1), str((ii+1)/AcqNum)))
            cX.append(efc.rPV('{}:Stats2:CentroidX_RBV'.format(PVhead)))
            cY.append(efc.rPV('{}:Stats2:CentroidY_RBV'.format(PVhead)))
            sX.append(BeamSig*efc.rPV('{}:Stats2:SigmaX_RBV'.format(PVhead)))
            sY.append(BeamSig*efc.rPV('{}:Stats2:SigmaY_RBV'.format(PVhead)))
            time.sleep(1/RepRateHz)
        fracX=(np.array(cX) - np.mean(cX)) / np.array(sX);
        fracY=(np.array(cY) - np.mean(cY)) / np.array(sY);
        fracTOT=np.sqrt(fracX**2 + fracY**2)
        print('Number of shots: {}, beam sigma: {} x sigma'.format(AcqNum,BeamSig))
        print('Average centroid_X (cX) and centroid_Y (cY) in pixels: cXavg={}px and cYavg={}px'.format(np.mean(cX),np.mean(cY)))
        print('Average sigma_x (sX) and sigma_y (sY) in pixels: sXavg={}px and sYavg={}px'.format(np.mean(sX),np.mean(sY)))
        print('StDev of cX relative to cXavg as \% of sX, i.e. fracX = (cX - cXavg)/sXavg): {}\%'.format(np.std(fracX)))
        print('StDev of cY relative to cYavg as \% of sY, i.e. fracY = (cY - cYavg)/sYavg): {}\%'.format(np.std(fracY)))
        print('StDev of total fractional relative drift, i.e. fracTOT=np.sqrt(fracX**2 + fracY**2): {}/%'.format(fracTOT))
        ep.ll([fracX, fracY, fracTOT])
        return [cX, cY, sX, sY]


        
    @classmethod
    def _QuickView(cls,CAMreq,ImageNo=2,LIVE=False,MAXLOOPS=25,endfreeze=False,reConfig=False):#like MEC:GIGE:31
        """
        returns 2D plot of the camera specified in CAMreq
        specify CAMreq using a Camera NickName, SPL Cam Name, or PV Name (i.e. from Name())
        Option: ImageNo should be 2 typically but could be set to 1 potentially
                (recall that GigE IOCs produce multiple image streams)
                (for AD cameras, try ImageNo=0 to trying talking to their ancient decrepit IOCs)
        Option: use LIVE=True to show a live image rather than a still image with LIVE=False
        Option: if LIVE=True then use MAXLOOPS to set the number of loops on the live view before finishing
        Option: reConfig=True means that the CAM.Config function will be re-executed at the conclusion of the loop
              : reConfig=False means that CAM.Config will not execute; tries to be as "read-only" as possible 
              : (advisable to set reConfig=False when hacking for use with non-MEC cameras!!)
        Option: endfreeze=True means that the acquisition is stopped after reaching MAXLOOPS
              : endfreeze=False means that the acquisition will continue even after the live view ends
        Example: use _QuickView('Regen') to get just a quick current look at the regen camera output
        Note: to view more than one camera at a time, use View() instead
        Note: using View() with only one argument is equivalent to using _QuickView() so using View() is preferred
        """
        PVhead=cls.Name(CAMreq)
        if reConfig:
            efc.wPV('{}:Acquire'.format(PVhead),1)#tries to start the camera first
        try:#if True:#try:
            if ImageNo in [1,2]:
                tres1=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
                tres2=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
                twf=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArrayData')
            elif ImageNo == 0:
                if int(CAMreq.split(":")[-1]) in [100,90,285,295,287,297,320,423]:
                    tres1=efc.rPV(PVhead+':ArraySizeY_RBV')
                    tres2=efc.rPV(PVhead+':ArraySizeX_RBV')
                    twf=efc.rPV(PVhead+':Image:ArrayData')
                elif int(CAMreq.split(":")[-1]) in [186,461,469]:
                    tres1=len(efc.rPV(PVhead+':PROJ_V'))
                    tres2=len(efc.rPV(PVhead+':PROJ_H'))
                    twf=efc.rPV(PVhead+':IMAGE')
                else:
                    print('Could not find the camera!: '+CAMreq)
            else:
                print('Invalid ImageNo!: '+str(ImageNo))
                return False
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
                    if reConfig:
                        cls.Config(CAMreq,LIVE=True)
                        efc.wPV('{}:Acquire'.format(PVhead),1)
                    twf=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArrayData')
                    if ImageNo in [1,2]:
                        twf=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArrayData')
                    else:
                        if int(CAMreq.split(":")[-1]) in [100,90,285,295,287,297,320,423]:
                            twf=efc.rPV(PVhead+':Image:ArrayData')
                        elif int(CAMreq.split(":")[-1]) in [186,461,469]:
                            twf=efc.rPV(PVhead+':IMAGE')
                        else:
                            print('Could not find the camera!: '+CAMreq)
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
        if reConfig:
            cls.Config(CAMreq,LIVE=False)
            if endfreeze == False:
                efc.wPV('{}:Acquire'.format(PVhead),1)
        return
    
    @classmethod
    def QuickSave(cls,CAMreq,ImageNo=1,FileNameQ='default',AbsFilePathQ='default'):#like MEC:GIGE:31
        """
        saves 2D plot of the camera specified in CAMreq without trying to change acquisition status
        specify CAMreq using a Camera NickName, SPL Cam Name, or PV Name (i.e. from Name())
        Option: ImageNo should be 2 typically but could be set to 1 potentially
                (recall that GigE IOCs produce multiple image streams)
                (for AD cameras, try ImageNo=0 to trying talking to their ancient decrepit IOCs)
        """
        PVhead=cls.Name(CAMreq)
        timestamp=datetime.now().strftime('%Y%m%d_%H%M%S')
        try:#if True:#try:
            if ImageNo in [1,2]:
                tres1=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
                tres2=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
                twf=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArrayData')
            elif ImageNo == 0:
                if int(CAMreq.split(":")[-1]) in [100,90,285,295,287,297,320,423]:
                    tres1=efc.rPV(PVhead+':ArraySizeY_RBV')
                    tres2=efc.rPV(PVhead+':ArraySizeX_RBV')
                    twf=efc.rPV(PVhead+':Image:ArrayData')
                elif int(CAMreq.split(":")[-1]) in [186,461,469]:
                    tres1=len(efc.rPV(PVhead+':PROJ_V'))
                    tres2=len(efc.rPV(PVhead+':PROJ_H'))
                    twf=efc.rPV(PVhead+':IMAGE')
                else:
                    print('Could not find the camera!: '+CAMreq)
            else:
                print('Invalid ImageNo!: '+str(ImageNo))
                return False
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
        if FileNameQ=='default':
            FileNameQ='GigE_'+cls.Name(CAMreq,returnAll=True)[0]
        if AbsFilePathQ=='default':
            AbsFilePathQ=GLOBAL.PSFILEPATH
        try:
            fig.savefig(AbsFilePathQ+FileNameQ+'_'+timestamp+'.png');
            print('File saved as '+AbsFilePathQ+FileNameQ+'_'+timestamp+'.png')
        except:
            print('Save failed!')
        return
    
    @classmethod
    def QuickSaveData(cls,CAMreq,ImageNo=1,FileNameQ='default',AbsFilePathQ='default'):#like MEC:GIGE:31
        """
        saves 2D array of the camera specified in CAMreq without trying to change acquisition status
        specify CAMreq using a Camera NickName, SPL Cam Name, or PV Name (i.e. from Name())
        Option: ImageNo should be 2 typically but could be set to 1 potentially
                (recall that GigE IOCs produce multiple image streams)
                (for AD cameras, try ImageNo=0 to trying talking to their ancient decrepit IOCs)
        """
        PVhead=cls.Name(CAMreq)
        timestamp=datetime.now().strftime('%Y%m%d_%H%M%S')
        try:#if True:#try:
            if ImageNo in [1,2]:
                tres1=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
                tres2=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
                twf=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArrayData')
            elif ImageNo == 0:
                if int(CAMreq.split(":")[-1]) in [100,90,285,295,287,297,320,423]:
                    tres1=efc.rPV(PVhead+':ArraySizeY_RBV')
                    tres2=efc.rPV(PVhead+':ArraySizeX_RBV')
                    twf=efc.rPV(PVhead+':Image:ArrayData')
                elif int(CAMreq.split(":")[-1]) in [186,461,469]:
                    tres1=len(efc.rPV(PVhead+':PROJ_V'))
                    tres2=len(efc.rPV(PVhead+':PROJ_H'))
                    twf=efc.rPV(PVhead+':IMAGE')
                else:
                    print('Could not find the camera!: '+CAMreq)
            else:
                print('Invalid ImageNo!: '+str(ImageNo))
                return False
            if len(twf) != tres1*tres2:
                twf=list(twf)+(tres1*tres2-len(twf))*[0]
            dat2d=np.array_split(np.array(twf),tres1);
        except:#if True:#except:
            print('Failed!')
            return
        if FileNameQ=='default':
            FileNameQ='GigE_'+cls.Name(CAMreq,returnAll=True)[0]
        if AbsFilePathQ=='default':
            AbsFilePathQ=GLOBAL.PSFILEPATH
        try:
            np.savetxt(AbsFilePathQ+FileNameQ+'_'+timestamp+'.txt', dat2d);
            print('Array file saved as '+AbsFilePathQ+FileNameQ+'_'+timestamp+'.txt')
        except:
            print('Array save failed!')
        return

    @classmethod
    def View(cls, *CAMargs,ImageNo=2,LIVE=False,MAXLOOPS=10,endfreeze=False,reConfig=False):
        """
        returns 2D plot of the cameras specified in *CAMargs
        specify *CAMargs using a Camera NickName, SPL Cam Name, or PV Name, listing all desired cameras
            separated only by commas (see Names() list for options)
        Example: View('Legend','StrInA','StrInB') plots a static view of the 'Legend', 'StrInA', and
                 'StrInB' cameras in subplots within a single frame
        Example: View('all') plots a static view of all cameras in subplots within a single frame
                 ('all' is equivalent to 'Regen', 'Trap', 'StrInA', 'StrInB', 'MPA1In', 'MPA1Out', 'MPA2In',
                  'MPA2Out', 'MPA2Xtal', 'CompIn', 'CompOutNF', and 'CompOutFF')
        Note: if only one camera is specified, then View('Trap') is equivalent to _QuickView('Trap')
        Note: due to bandwidth limitations and other various performance issues, refresh rate may suffer
              according to how many cameras are selected!
        Option: ImageNo should be 2 typically but could be set to 1 potentially
                (recall that GigE IOCs produce multiple image streams)
                (for AD cameras, try ImageNo=0 to trying talking to their ancient decrepit IOCs)
        Option: use LIVE=True to show a live image rather than a still image with LIVE=False
        Option: if LIVE=True then use MAXLOOPS to set the number of loops on the live view before finishing
        Option: reConfig=True means that the CAM.Config function will be re-executed at the conclusion of the loop
              : reConfig=False means that CAM.Config will not execute; tries to be as "read-only" as possible 
              : (advisable to set reConfig=False when hacking for use with non-MEC cameras!!)
        Option: endfreeze=True means that the acquisition is stopped after reaching MAXLOOPS
              : endfreeze=False means that the acquisition will continue even after the live view ends
        """
        if isinstance(CAMargs[0],tuple) or isinstance(CAMargs[0],list):
            CAMargs=tuple(CAMargs[0])#try to catch case of someone accidentally entering input as tuple or list
        if CAMargs == ('all',):
            CAMargs = ('Regen', 'Trap', 'StrInA', 'StrInB', 'MPA1In', 'MPA1Out', 'MPA2In', 'MPA2Out', 'MPA2Xtal', 'CompIn', 'CompOutNF', 'CompOutFF')
        if len(CAMargs) == 1:
            cls._QuickView(*CAMargs,ImageNo=ImageNo,LIVE=LIVE,MAXLOOPS=MAXLOOPS,reConfig=reConfig)
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
                tPVhead=cls.Name(CAMargs[ii])
                tPVheadL.append(tPVhead);
                if reConfig:
                    efc.wPV('{}:Acquire'.format(tPVhead),1)#tries to start each camera               
                if ImageNo in [1,2]:
                    tres1=efc.rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
                    tres2=efc.rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArraySize0_RBV')
                    tres3=efc.rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
                    tresL=sorted([tres1,tres2,tres3],reverse=True) 
                    twf=efc.rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArrayData')
                    #
                    #tres1=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize1_RBV')
                    #tres2=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArraySize2_RBV')
                    #twf=efc.rPV(PVhead+':IMAGE'+str(ImageNo)+':ArrayData')
                elif ImageNo == 0:
                    if int(tPVhead.split(":")[-1]) in [100,90,285,295,287,297,320,423]:
                        tres1=efc.rPV(tPVhead+':ArraySizeY_RBV')
                        tres2=efc.rPV(tPVhead+':ArraySizeX_RBV')
                        twf=efc.rPV(tPVhead+':Image:ArrayData')
                    elif int(tPVhead.split(":")[-1]) in [186,461,469]:
                        tres1=len(efc.rPV(tPVhead+':PROJ_V'))
                        tres2=len(efc.rPV(tPVhead+':PROJ_H'))
                        twf=efc.rPV(tPVhead+':IMAGE')
                    else:
                        print('Could not find the camera!: '+tPVhead)
                    #tres3=tres2[:]
                    tresL=[tres1,tres2];#sorted([tres1,tres2],reverse=True) 
                else:
                    print('Invalid ImageNo!: '+str(ImageNo))
                    return False
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
                        if ImageNo in [1,2]:
                            twf=efc.rPV(tPVhead+':IMAGE'+str(ImageNo)+':ArrayData')
                        else:
                            if int(tPVhead.split(":")[-1]) in [100,90,285,295,287,297,320,423]:
                                twf=efc.rPV(tPVhead+':Image:ArrayData')
                            elif int(tPVhead.split(":")[-1]) in [186,461,469]:
                                twf=efc.rPV(tPVhead+':IMAGE')
                            else:
                                print('Could not find the camera!: '+tPVhead)
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
        if reConfig:
            for tPVhead in CAMargs:
                cls.Config(cls.Name(tPVhead),LIVE=False)
                if endfreeze == False:
                    efc.wPV('{}:Acquire'.format(cls.Name(tPVhead)),1)
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

    @classmethod
    def DisableOverlay(cls, CAMreq, ImageNo=2):
        PVhead=cls.Name(CAMreq)
        efc.wPV('{}:Acquire'.format(PVhead),0)#hopefully prevents IOC crashes, per TylerJ
        for ii in range(8):#eight overlays possible
            try:
                #set up everything...
                OverlayNameField='{}{}'.format('Box' if ii<4 else 'Cross',(ii%4)+1)
                OverlayName='{}:IMAGE{}:{}'.format(PVhead,ImageNo,OverlayNameField)
                efc.wPV(OverlayName+':Use', 1)#0=off,1=on
            except:
                print('Failed to disable overlay!')
        return



    @classmethod
    def Config(cls, CAMreq, RefXhairXY=[-1,-1], InstructionStr='', MisalignmentTolerance=-1, ImageNo=2, RefXhairXYsize=[-1,-1],LIVE=False):
        """ 
        configures stats, overlay plug-ins, etc. of a given camera specified by CAMreq
        specify CAMreq using a Camera NickName, SPL Cam Name, or PV Name (i.e. from Name())
        Option: RefXhairXY sets the desired static reference crosshair position,
                for example indicating the "targeted" location for a laser beam
              : default RefXhairXY=[-1,-1] indicates not to change anything, as it was set up previously
        Option: InstructionStr sets the desired instruction string to be displayed in the corner of the camera image,
                for example giving a helpful tip/reminder to the viewer (e.g. which motor to use for tuning, etc.)
              : default InstructionStr='' indicates not to change anything, as it was set up previously
        Option: MisalignmentTolerance specifies a threshold for the absolute distance (in pixels) between RefXhaiXY
                and the current beam centroid; if difference > or < MisalignmentTolerance then different messages
                will be printed to the output screen
              : default MisalignmentTolerance=-1 indicates not to change anything, as it was set up previously
        Option: RefXhairXYsize sets the desired size of the static reference crosshair position,
                for example indicating the approximate size of the laser beam when centered on the "targeted" location
              : default RefXhairXYsize=[-1,-1] indicates not to change anything, as it was set up previously
        Option: ImageNo should be 2 typically but could be set to 1 potentially
                (recall that GigE IOCs produce multiple image streams)
        Option: LIVE will cause certain overlays to appear or to be hidden, depending on if it is True or False
              : This is because certain chosen overlays require Python calculations that only occur while the loop is running
              : In order to avoid confusion, these Python-calculated overlays are hidden when LIVE=False and visible when LIVE=True
              : Examples of such overlays are the Stats, Counts, and AlignmentOK overlays
        
        The GigE cam IOC typically supports 8 overlays on a single image; this function configures the 8 overlays as follows:
            1: Ref X-hair #static reference crosshair, typically gives a target for the live beam
            2: DynCentroid #dynamic centroid whose position and size follow the position and size of the laser
            3: CAMname #text overlay of the name of the camera, meant to help avoid confusing cameras for each other
            4: Instructions #meant to provide helpful tips to viewers, e.g. on what action to take if the beam is misaligned
            5: TimeStamp #adds the current timestamp, useful for checking how recent the visible image is, etc.
            6: Stats #uses Python to calculate and display centroid X and Y and sigma X and Y on-screen
            7: Counts #uses Python to calculate and display camera's Max and Total counts on-screen
            8: AlignmentOK #displays a conditional message depending on if the misalignment exceeds the MisalignmentTolerance or not
        
        Nominal values for the MEC SPL camera configurations are kept in ConfigReset(); if things get screwed up (because sometimes
        the IOC just has a bad day and screws everything up), you can run ConfigReset() and *hopefully* get back to normal
        
        (not advised to use this on custom cameras unless you are REALLY brave and daring)
        """
        PVhead=cls.Name(CAMreq)
        if LIVE==False:
            efc.wPV('{}:Acquire'.format(PVhead),0)#hopefully prevents IOC crashes, per TylerJ
        ArraySize=[efc.rPV('{}:IMAGE{}:ArraySize0_RBV'.format(PVhead,ImageNo)),efc.rPV('{}:IMAGE{}:ArraySize1_RBV'.format(PVhead,ImageNo))]
        NameList=['Ref X-hair','DynCentroid','CAMname','Instructions','TimeStamp','Stats','Counts',
                  'AlignmentOK'+str(str(MisalignmentTolerance).zfill(3) if MisalignmentTolerance > 0 else efc.rPV('{}:IMAGE{}:Cross4:Name'.format(PVhead,ImageNo))[-3:])]
        ShapeList=[0,2,3,3,3,3,3,3]#0=Cross,1=Rect,2=Text,3=Ellipse;;;DOCUMENTATION WRONG!! 2=Ellipse, 3=Text!!!!
        #each overlay needs the following information:
        #[SizeXLink.DOL, SizeX, PositionXLink.DOL, PositionX, CenterXLink.DOL, CenterX, 
        # SizeYLink.DOL, SizeY, PositionYLink.DOL, PositionY, CenterYLink.DOL, CenterY, 
        # TimeStampFormat]
        #Size affects Position/Center, so set Size first then Position/Center
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
         '{:<10.9} {:<11.10} {:<11.11}'.format(*cls.Name(CAMreq,returnAll=True))],
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
                    efc.wPV('{}:IMAGE{}:NDArrayPort'.format(PVhead,ImageNo),'IMAGE{}:Over'.format(ImageNo))
                    efc.wPV('{}:IMAGE{}:Over:EnableCallbacks'.format(PVhead,ImageNo),1)
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
#         elif command == 'ENA':1c1t5
#             DG_8_EF_Polarity = 'NEG'
#         else:
#             print('Unanticipated command case!')
#         return
# =============================================================================

    @classmethod
    def ConfigReset(cls,CAMreq='all'):#not sure why these change sometimes -- config file somewhere? camera problem? 
        """
        if default CAMreq='all', executes the Config function serially on SPL1-10 with their own default values
        can also accept a single specified camera input to reset
        """

        if CAMreq == 'all':
            ResetList=['Legend','StrInA','StrInB','MPA1In','MPA1Out','MPA2In','MPA2Xtal','MPA2Out','CompInNF','CompInFF','CompOutNF','CompOutFF','TTIn','TTOut']
        else:
            ResetList=[cls.Name(CAMreq,returnAll=True)[0]]

        for eaEntry in ResetList:
            if eaEntry == 'Legend':
                cls.Config(CAMreq='Legend',RefXhairXY=[351,222],MisalignmentTolerance=40,ImageNo=1, RefXhairXYsize=[100,100],LIVE=False,InstructionStr='DO NOT tune to xhair!')
            elif eaEntry == 'StrInA':
                cls.Config(CAMreq='StrInA',RefXhairXY=[279,241],MisalignmentTolerance=25,ImageNo=1, RefXhairXYsize=[50,50],LIVE=False,InstructionStr='Use SM0 (NOT SM1!) to put centroid on xhair!')
            elif eaEntry == 'StrInB':
                cls.Config(CAMreq='StrInB',RefXhairXY=[342,221],MisalignmentTolerance=40,ImageNo=1, RefXhairXYsize=[5,5],LIVE=False,InstructionStr='Use SM2 to put centroid on xhair!')
            elif eaEntry == 'MPA1In':
                cls.Config(CAMreq='MPA1In',RefXhairXY=[139,221],MisalignmentTolerance=15,ImageNo=1, RefXhairXYsize=[80,100],LIVE=False,InstructionStr='Use SM3 to put centroid on xhair!')
            elif eaEntry == 'MPA1Out':
                cls.Config(CAMreq='MPA1Out',RefXhairXY=[172,200],MisalignmentTolerance=15,ImageNo=1, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM4 to maximize spot!')
            elif eaEntry == 'MPA2In':
                cls.Config(CAMreq='MPA2In',RefXhairXY=[301,226],MisalignmentTolerance=15,ImageNo=1, RefXhairXYsize=[220,180],LIVE=False,InstructionStr='Use SM5(V) to put centroid on xhair!')
            elif eaEntry == 'MPA2Xtal':
                cls.Config(CAMreq='MPA2Xtal',RefXhairXY=[419,201],MisalignmentTolerance=40,ImageNo=1, RefXhairXYsize=[90,80],LIVE=False,InstructionStr='Use SM7 (and SM6) to put centroid on xhair!')
            elif eaEntry == 'MPA2Out':
                cls.Config(CAMreq='MPA2Out',RefXhairXY=[377,212],MisalignmentTolerance=40,ImageNo=1, RefXhairXYsize=[320,280],LIVE=False,InstructionStr='Use SM7 to put centroid on xhair!')
            elif eaEntry == 'CompOutNF':
                cls.Config(CAMreq='CompOutNF',RefXhairXY=[230,230],MisalignmentTolerance=40,ImageNo=1, RefXhairXYsize=[150,150],LIVE=False,InstructionStr='Do something!!')
            elif eaEntry == 'TTIn':
                cls.Config(CAMreq='TTIn',RefXhairXY=[377,263],MisalignmentTolerance=50,ImageNo=1, RefXhairXYsize=[20,20],LIVE=False,InstructionStr='Do something!!')
            elif eaEntry == 'TTOut':
                cls.Config(CAMreq='TTOut',RefXhairXY=[377,263],MisalignmentTolerance=50,ImageNo=1, RefXhairXYsize=[20,20],LIVE=False,InstructionStr='Do something!!')
            else:
                print('Camera configuration not found: '+eaEntry)
        return

#            elif eaEntry == 'CompInNF':
#                cls.Config(CAMreq='CompInNF',RefXhairXY=[276/2,252/2],MisalignmentTolerance=50,ImageNo=1, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM8 to put centroid on xhair!')
#            elif eaEntry == 'CompInFF':
#                cls.Config(CAMreq='CompInFF',RefXhairXY=[276/2,252/2],MisalignmentTolerance=50,ImageNo=1, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM8 to put centroid on xhair!')
#            elif eaEntry == 'CompOutFF':
#                cls.Config(CAMreq='CompOutFF',RefXhairXY=[292,300],MisalignmentTolerance=25,ImageNo=1, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use XPS2 Mid (and Input) Mirror to align to xhair!')
#            elif eaEntry == 'Trap':
#                cls.Config(CAMreq='Trap',RefXhairXY=[215,156],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Watch for back reflections!!')


    @classmethod
    def _ConfigResetOLD(cls,CAMreq='all'):#not sure why these change sometimes -- config file somewhere? camera problem? 
        """
        if default CAMreq='all', executes the Config function serially on SPL1-10 with their own default values
        can also accept a single specified camera input to reset
        """

        if CAMreq == 'all':
            ResetList=['Legend','StrInA','StrInB','MPA1In','MPA1Out','MPA2In','MPA2Xtal','MPA2Out','CompIn','CompOutFF','CompOutNF','Trap','SPFloat1','CompInNF','TTIn','TTOut']
        else:
            ResetList=[cls.Name(CAMreq,returnAll=True)[0]]

        for eaEntry in ResetList:
            if eaEntry == 'Legend':
                cls.Config(CAMreq='Legend',RefXhairXY=[359,251],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[100,100],LIVE=False,InstructionStr='DO NOT tune to xhair!')
            elif eaEntry == 'StrInA':
                cls.Config(CAMreq='StrInA',RefXhairXY=[335,257],MisalignmentTolerance=25,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM0 (NOT SM1!) to put centroid on xhair!')
            elif eaEntry == 'StrInB':
                cls.Config(CAMreq='StrInB',RefXhairXY=[334,243],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM2 to put centroid on xhair!')
            elif eaEntry == 'MPA1In':
                cls.Config(CAMreq='MPA1In',RefXhairXY=[152,147],MisalignmentTolerance=15,ImageNo=2, RefXhairXYsize=[25,25],LIVE=False,InstructionStr='Use SM3 to put centroid on xhair!')
            elif eaEntry == 'MPA1Out':
                cls.Config(CAMreq='MPA1Out',RefXhairXY=[411,249],MisalignmentTolerance=15,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM4 to maximize spot!')
            elif eaEntry == 'MPA2In':
                cls.Config(CAMreq='MPA2In',RefXhairXY=[366,276],MisalignmentTolerance=15,ImageNo=2, RefXhairXYsize=[50,50],LIVE=False,InstructionStr='Use SM5(V) to put centroid on xhair!')
            elif eaEntry == 'MPA2Out':
                cls.Config(CAMreq='MPA2Out',RefXhairXY=[385,259],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[190,190],LIVE=False,InstructionStr='Use SM7 (and SM6) to put centroid on xhair!')
            elif eaEntry == 'CompIn':
                cls.Config(CAMreq='CompIn',RefXhairXY=[276/2,252/2],MisalignmentTolerance=50,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use SM8 to put centroid on xhair!')
            elif eaEntry == 'CompOutFF':
                cls.Config(CAMreq='CompOutFF',RefXhairXY=[292,300],MisalignmentTolerance=25,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use XPS2 Mid (and Input) Mirror to align to xhair!')
            elif eaEntry == 'CompOutNF':
                cls.Config(CAMreq='CompOutNF',RefXhairXY=[364,277],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Use XPS2 Input(/Mid) Mirror to align to xhair!')
            elif eaEntry == 'SPFloat1':
                cls.Config(CAMreq='SPFloat1',RefXhairXY=[360,230],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='test')
            elif eaEntry == 'Trap':
                cls.Config(CAMreq='Trap',RefXhairXY=[215,156],MisalignmentTolerance=40,ImageNo=2, RefXhairXYsize=[40,40],LIVE=False,InstructionStr='Watch for back reflections!!')
            elif eaEntry == 'CompInNF':
                cls.Config(CAMreq='CompInNF',RefXhairXY=[304,210],MisalignmentTolerance=50,ImageNo=2, RefXhairXYsize=[180,180],LIVE=False,InstructionStr='Use SM7 to center beam on xhair!')
            else:
                print('Camera configuration not found: '+eaEntry)
        return
                






























                  
class TTL_shutter:
    """
    Class for organizing utilities for the TTL-triggerable beam shutters in the LPL and SPL
    Current shutter locations (and names) are as follows:
        :LPL: 4 x immediately before the 2" heads on 'AB', 'EF', 'GH', 'IJ'
            : 2 x immediately before the periscopes going into the chamber, 'WEST 527' ('WW') and 'EAST 527' ('XX')
        :SPL: 1 x before the compressor grating of the Legend ('REGEN' or 'ZZ')
    Shutters do not possess internal state sensor, so state of each shutter (open vs closed) tracked using notepad PVs
    Typical usage via TTL_shutter.[command]
    Possible commands include:
        :Status #provides the current status (open vs closed) of each shutter
        :Toggle #allows control of shutters to open/close/toggle state
        :Refresh #refresh state tracker if shutters are touched manually or Beckhoff has problems, etc.
    """
    def Status(display=True):
        """
        Provides the current status (open vs closed) of all shutters, returning an array of 0s and 1s
        Array values correspond to shutters in the following order:
            : ['AB','EF','GH','IJ','WEST 527','EAST 527','REGEN']
        Example: if all shutters are closed except for 'AB', 'EF', and 'WEST 527', the output will be
                 [0,0,1,1,0,1,1]
        Option: if display=True then a status report is printed out to terminal in addition to returning the status array
              : if display=False then only the status array is returned; nothing is printed to terminal
        """
        pvstatuslist=['MEC:LAS:FLOAT:'+str(ii) for ii in range(14,21)];
        statuslist=[]
        for eapv in pvstatuslist:
            temppv=EpicsSignal(eapv);
            statuslist.append(int(temppv.get()))#0=open,1=closed
        if display:
            shutlist=['AB','EF','GH','IJ','WEST 527','EAST 527','REGEN']
            print('(0=out,1=in) '+', '.join([ea_shut+':'+efc.cstr(ea_stat,'BRR' if ea_stat else '') for ea_shut,ea_stat in zip(shutlist,statuslist)]))
        return statuslist

    @classmethod
    def Toggle(cls,ArmStrQ,display=True):
        """
        Toggles the current state of all shutters mentioned in ArmStrQ and then returns the shutter status at the end
        Using 'open' or 'close' at the beginning of ArmStrQ simply makes sure that all specified arms end up in the specified state
            rather than just toggling the state, i.e. it will toggle specified arms only if the specified state is not already satisfied
        Shutter options are 'AB', 'EF', 'GH', 'IJ', 'WW' (west 527), 'XX' (east 527), and 'ZZ' (SPL regen)
        Using 'all' is equivalent to using 'ABEFGHIJWWXX' (note that the regen shutter 'ZZ' is left out)
        Example: Toggle('ABEF') will toggle the 'AB' and 'EF' shutters from whatever their current state is,
                 i.e. if 'AB' was open and 'EF' was closed then Toggle('ABEF') would close 'AB' and open 'EF'
        Example: Toggle('openABIJ') will make sure 'AB' and 'IJ' shutters end up in the open position,
                 i.e. if 'AB' was closed and 'IJ' was open then 'AB' would be toggled to open and 'IJ' would remain open
        Example: Toggle('closeall') will make sure all LPL shutters end up in the closed position, 
                 i.e. if 'AB', 'EF', 'GH', and 'IJ' were initially open and 'WW' and 'XX' were initially closed then
                 'AB', 'EF', 'GH', and 'IJ' would be toggled to closed and 'WW' and 'XX' would remain in the closed position
        Example: Toggle('all') toggles the state of all LPL shutters so that the configuration is exactly the opposite of the initial state,
                 i.e. if 'AB', 'EF', 'GH', and 'IJ' were initially open and 'WW' and 'XX' were initially closed then
                 'AB', 'EF', 'GH', and 'IJ' would be toggled to closed and 'WW' and 'XX' would be toggled to opened
        Option: display=True reads back and prints the initial state, makes the toggle, and then reads back and prints the final state
        Option: display=False still returns the shutter status at the end but does not print anything to terminal
        """
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
        """
        Refreshes state tracker if shutters are touched manually or Beckhoff has problems, etc.
        Ensure all shutters are open before executing the function in order to insure that the notepad PVs reset correctly
        """
        pvlist=['MEC:LAS:FLOAT:'+str(ii) for ii in range(14,21)];
        inddesclist=['AB','EF','GH','IJ','WEST (ABEF)','EAST (GHIJ)','Regen']
        desclist=[inddesc+' shutter state' for inddesc in inddesclist]
        valulist=[0,0,0,0,0,0,0];
        print('This will reset the shutter counter! Are you sure you want to continue? [y/n]',end='',flush=True)
        checkprompt=efc.getch_with_TO(TOsec=10,display=False);
        if checkprompt not in ('y','Y'):
            print('Try again later then!');
            return
        else:
            print('OK, I hope you know what you\'re doing!')
        print('Resetting shutter PV statuses... Please insure all shutters are actually open!!')
        for jj in range(len(pvlist)):
            temppv1=EpicsSignal(str(pvlist[jj]+'.DESC'));temppv2=EpicsSignal(pvlist[jj]);
            temppv1.put(desclist[jj]);temppv2.put(valulist[jj]);





class DG645:
    """
    Potential future class containing DG645-related utilities such as:
        - take a snapshot of settings to be restored later with a restore function
        - restore saved data from a snapshot function into the corresponding DG645
        - perform snapshots of multiple DG boxes at the same, etc.
        - change labels of channels to be more descriptive
    """      
# =============================================================================
#     def Refresh():
#         """
#         
#         """
#         pvlist=[['MEC:LAS:DDG:0'+str(numii)+':'+chii+'DelaySI.DESC' for chii in ['a','c','e','g']] for numii in [1,2,6,8]];
#         desclist=[['A:PS LATE','C:INH PS EARLY','E:unused','G: unused'],['A:PS EARLY','C:unused','E:EvoHE1','G:EvoHE2'],['A:GaiaQSW','C:GaiaLamp','E:INH UNI','G:GigE TRIG IN'],['A:BIG UNI','C:small UNI','E:INH GigE','G:unused']]
#         for eachboxii in range(len(pvlist)):
#             for eachentryii in range(len(eachboxii)):
#                 temppv=EpicsSignal(pvlist[eachboxii][eachentryii])
#                 temppv.put(desclist[eachboxii][eachentryii])
# 
#     def Snapshot():
#         """
#         saves data into file that has array structure for easy unpacking
#         """
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
#         with open(str(LPL.psfilepath()+'DG/snapshot'+LPL._DateString()+'.txt'),'a') as out:
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
# 
#     def ChannelDescription(DGbox):
#         """  """
#         pass
    pass





class SPL:#pointing, alignment, and other macros + automation routines
    """
    Potential future class containing SPL-related utilities such as:
        - pointing functions
        - aligment functions
        - other macros + automation routines
    """
    @classmethod
    def alignment_mode(cls):
        efc.wPV('MEC:LAS:EVR:01:TRIG3:TCTL',0)#disable SPL Slicer trigger from New LAS EVR Front Panel 3       
        efc.wPV('MEC:LAS:EVR:01:TRIG4:TCTL',0)#disable UNIBLITZ trigger from New LAS EVR Front Panel 4
        efc.wPV('MEC:LAS:DDG:08:abOutputPolarityBO',1)#open 65mm shutter by changing DG645 output polarity to POS
        efc.wPV('MEC:LAS:DDG:08:cdOutputPolarityBO',1)#open 6mm shutter by changing DG645 output polarity to POS
        efc.wPV('MEC:LAS:DDG:08:triggerInhibitMO',3)#enable inhibit on 65mm shutter and 6mm shutter DG645 outputs
        cls.gige_trig_5Hz()

    @classmethod
    def mpa2ss_mode(cls):
        efc.wPV('MEC:LAS:EVR:01:TRIG3:TCTL',0)#disable SPL Slicer trigger from New LAS EVR Front Panel 3       
        efc.wPV('MEC:LAS:EVR:01:TRIG4:TCTL',1)#enable UNIBLITZ trigger from New LAS EVR Front Panel 4
        efc.wPV('MEC:LAS:EVR:01:TRIG4:TEC',177)#write single-shot UNIBLITZ event code to New LAS EVR Front Panel 4
        efc.wPV('MEC:LAS:DDG:08:cdOutputPolarityBO',0)#close 6mm shutter by changing DG645 output polarity to NEG
        efc.wPV('MEC:LAS:DDG:08:abOutputPolarityBO',0)#close 65mm shutter by changing DG645 output polarity to NEG
        efc.wPV('MEC:LAS:DDG:08:triggerInhibitMO',0)#disable inhibit on 65mm shutter and 6mm shutter DG645 outputs
        cls.gige_trig_ss()

    @staticmethod
    def gige_trig_ss():
        GigElist=['24','22','28','23','25','26','27','30','31','16','32','33','34','35']#skipping :09 MPA2Xtal 
        for camhandle in GigElist:
            try:
                efc.wPV('MEC:GIGE:'+camhandle+':TriggerMode',1)#set camera to external trigger mode
                efc.wPV('MEC:GIGE:'+camhandle+':Acquire',1)#set camera to acquire; worry about timeout??
            except:
                print('Failed to set: MEC:GIGE:'+camhandle)
        efc.wPV('MEC:LAS:DDG:08:efOutputPolarityBO',1)#set POS polarity on Uniblitz DG645 EF which inhibits the GigE trigger server DG535

    @staticmethod
    def gige_trig_5Hz():
        GigElist=['24','22','28','23','25','26','27','30','31','16','32','33','34','35']#skipping :09 MPA2Xtal 
        for camhandle in GigElist:
            try:
                efc.wPV('MEC:GIGE:'+camhandle+':TriggerMode',1)#set camera to external trigger mode
                efc.wPV('MEC:GIGE:'+camhandle+':Acquire',1)#set camera to acquire; worry about timeout??
            except:
                print('Failed to set: MEC:GIGE:'+camhandle)
        efc.wPV('MEC:LAS:DDG:08:efOutputPolarityBO',0)#set NEG polarity on Uniblitz DG645 EF which DOESN'T inhibit the GigE trigger server DG535 -- it will always externally trigger at 5Hz




class UNIBLITZ:
    """
    Potential future class containing UNIBLITZ-related utilities such as:
        - toggle system (DG boxes, etc.) triggering configurations (e.g. alignment vs SS vs burst modes, etc.)
        - toggle shutter state
        - toggle trigger state
    """
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





class Spectrometer:
    """
    Potential future class containing spectrometer-related utilities such as:
        - Qmini configuration
        - Qmini readout
        - floater Ocean Optics/other spectrometers?
    """
    pass





class VISAR:
    """
    Potential future class containing VISAR-related utilities such as:
        - laser triggering/any other EPICS-connected functions (may need to be expanded)
        - streak cameras (timing or configuration with shot or whatever is useful)
    """
    pass





class CtrlSys:
    """
    Class for monitoring (and potentially controlling) laser-related control systems
    Typical usage via CtrlSys.[command]
    Possible commands include:
        :cmp_checker #checks a list of relevant computers to see if they are pinging the network
        :pv_checker #checks a list of relevant PVs for current RBV, to see if their IOCs/host are live, etc.
    Potential future improvements:
        - expand list of PVs checked in pv_checker
        - add functionality for checking current PV RBV vs historical reference or allowed range
        - consider adding functionality of ping, netconfig search, grep_pv, grep_ioc, serverStat, imgr, etc.
        - when pinging motors, check the control box itself too? (see Stage.NewportBrowserCtrl)
    """
    @staticmethod
    def _plzchkpv(inpvnam):
        """
        internal function meant for checking PV read-back values, used as part of pv_checker
        takes a single PV as input, outputs the message to be displayed as part of pv_checker
        """
        try:
            currval=efc.rPV(inpvnam,display=False)
            if currval is False:
                msgout='rPV fail!'
            else:
                currvalstr=str(currval)
                if len(currvalstr) > 10:
                    currvalstr=currvalstr[:10]+'...'
                msgout = currvalstr+' vs oldval'
        except TimeoutError as err:
            msgout = 'Timeout!';err;
        return msgout

    @staticmethod
    def _plzchksrv(inpvnam):
        """
        internal function meant for checking server status, used as part of pv_checker
        takes a single PV as input, outputs the message to be displayed as part of pv_checker
        (internally uses grep_pv, grep_ioc, ping, and netconfig search to gather necessary information)
        """
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

    @staticmethod
    def _plzchkcmp(cmpnam):
        """
        internal function meant for checking computer status as part of cmp_checker
        takes a single computer name as input (e.g. from _MECcompylist()),
            outputs the message to be displayed as part of cmp_checker
        """
        try:
            netstat=os.system('ping -c 1 -w2 '+cmpnam[0]+' > /dev/null 2>&1')
            msgout=['{:<15}'.format(cmpnam[1])+'{:<28}'.format(cmpnam[0])+efc.cstr('{:<6}'.format(str('true ' if netstat==0 else 'false')),str('blink,r' if netstat!=0 else ''))]
        except:
            msgout=['{:<15}'.format(cmpnam[1])+'{:<28}'.format(cmpnam[0])+efc.cstr('{:<6}'.format('fail!'),'blink,r')]
        return msgout

# =============================================================================
#     @staticmethod
#     def MECcompylist000():
#         qqip=['172.21.46.147','172.21.46.148','172.21.46.146','172.21.46.60','172.21.46.128','172.21.46.100', '172.21.46.120','172.21.46.159', '172.21.46.197','172.21.46.142','172.21.46.70','172.21.46.78', '172.21.46.71','172.21.46.88','172.21.46.198','172.21.46.213','172.21.46.215','172.21.46.136', '172.21.46.218','172.21.46.219','172.21.46.182','172.21.46.144'];
#         qqn=['evo1','evo2','gaia','lecroy1','lecroy2','lecroya','lecroyb','PIMikroMove','spider','spectrometer', 'tundra','topas','visar1','visar2','vitara','rga','emp','phasicslaptop','phasics1','phasics2','dacage','legend']
#         nmlist=['mec-las-laptop06','mec-las-laptop07','mec-las-laptop05','scope-ics-mectc1-1','scope-ics-meclas-lecroy01','scope-ics-meclas-lecroy-a','scope-ics-meclas-lecroy-b','mec-las-laptop09','mec-las-laptop11','mec-las-laptop01','win-ics-mec-tundra','mec-las-laptop12','win-ics-mec-visar1','win-ics-mec-visar2','mec-las-vitara','mec-rga-laptop','scope-ics-mec-tektronix','mec-phasics-laptop01','win-ics-mec-phasics01','win-ics-mec-phasics02','mec-visar-cage','mec-las-laptop03']
#         return list(zip(nmlist,qqip,qqn))
# =============================================================================

    @staticmethod
    def _MECcompylist():
        """
        internal list of MEC computers to check as part of cmp_checker()
        """
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
        """
        checks a list of critical MEC computers to see if they are pinging or not, prints the status report to terminal
        (note: just because a machine is not pinging does not mean it is necessarily offline -- it may just not be responding to pings!)
        """
        cmpmsg=[];
        #qqpl=[eacmp[0] for eacmp in MECcompylist()]
        qqpl=cls._MECcompylist()
        with multiprocessing.Pool() as pool:
            cmpmsg=pool.map(cls._plzchkcmp,qqpl)
        print('{:<15}'.format('Computer name')+'{:<28}'.format('IP shorthand')+'{:<6}'.format('Ping?'))
        for ii in range(len(cmpmsg)):
            print(''.join(cmpmsg[ii]))
        return

    @classmethod
    def pv_checker(cls,pv='lpl'):
        """
        Checks a list of PVs critical for running either the MEC SPL or LPL laser systems or any user-supplied PV(s)
        If no argument is provided, it is assumed that one wants to check the PVs for the LPL (hence pv='lpl')
        Instead, if pv='spl' then the PVs for the SPL are checked instead
        Note: list of PVs to be checked for the LPL or SPL laser can be edited in the files
              GLOBAL.PSFILEPATH+'_ps_pvlist_lpl.txt' and GLOBAL.PSFILEPATH+'_ps_pvlist_spl.txt'
        
        The function prints a report of the PV statuses, which includes PV name, IOC name, Host name, Host location, Ping?, and PV value?
        Returns False if any of the PVs/hosts timeout or fail to ping
        Returns True if none of the PVs/hosts timeout or fail to ping
        """

        if isinstance(pv,str):
            if (pv.lower() == 'lpl'):
                qqpl=np.genfromtxt(GLOBAL.PSFILEPATH+'_ps_pvlist_lpl.txt',delimiter='\n',dtype=str);
            elif (pv.lower() == 'spl'):
                qqpl=np.genfromtxt(GLOBAL.PSFILEPATH+'_ps_pvlist_spl.txt',delimiter='\n',dtype=str);
            else:
                qqpl=(pv,);
        elif (isinstance(pv, tuple)) or (not isinstance(pv, list)):
            qqpl=pv;
        else:
            print('Please input a valid PV or list of PVs to check!')
            return False
        
        pvmsg=[];
        for eapv in qqpl:
            pvret=cls._plzchkpv(eapv)
            pvmsg.append(efc.cstr(pvret,'BRR,BLINK' if (pvret in ('rPV fail!','Timeout!')) else ''))
        with multiprocessing.Pool() as pool:
            srvmsg=pool.map(cls._plzchksrv,qqpl)
        print(''.join(['{:<34}'.format('PV name'),'{:<26}'.format('IOC name'),'{:<19}'.format('Host name'),'{:<25}'.format('Host location'),'{:<6}'.format('Ping?'),'{:<15}'.format('PV value?')]))
        for ii in range(len(pvmsg)):
            fullmsg=''.join(srvmsg[ii])+pvmsg[ii]
            print(fullmsg)
        if any(errmsg in fullmsg for errmsg in ('Timeout','Fail', 'false', 'False')):
            return False
        else:
            return True



                  
class SCALLOPS:
    """
    Potential future class containing VISAR-related utilities such as:
        - laser triggering/any other EPICS-connected functions (may need to be expanded)
        - streak cameras (timing or configuration with shot or whatever is useful)
    """
# =============================================================================
#       potential future work
#       - port everything over from Bethany and Patin
#       - improved background subtraction needed for improved transfer function accuracy?
#       - need LPL.Deconv too?
#       - test shot routine with triangle pulse (etc.) to grab day's starting transfer function guess?
#       - fit model that makes everything most self consistent; what gets weighed most heavily?
#           - (note: _EW2 probably not useful in that case)
#       - what will recipe look like when based on SCALOPS? what will recipe contain?
#       - 
#       - 
#       - 
# =============================================================================
    pass





class LabEnv:
    """
    Potential future class containing lab environment-related utilities such as:
        - air/enclosure/rack temperature/humidity/etc.
        - dust/particulate levels
    """
    pass





class RIS:
    """
    Potential future class containing RIS-related utilities such as:
        - easy RBV on state, radiation sensors, etc.
    """
    pass




                  
class PDU:
    """
    Potential future class containing PDU-related utilities such as:
        - quick and unified access/configuration/toggling of lab PDUs
    """
    pass





class GLOBAL:
    """
    Class meant to serve like global variables for all other meclas classes
    Idea is that only the contents of this class need to change if PV names change, calibrations or other
        constants need to be updated, etc.
    Typical usage via GLOBAL.[attribute] where attributes might be PVs, constants, arrays, strings, etc.
    
    Potential future improvements:
        - convert more of the parameters above into GLOBAL attributes for the sake of simplicity, clarity, longevity
        - add more restoration-type functions
    """
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
    CurrShape=EpicsSignal('MEC:LAS:FLOAT:12.DESC');
    CurrShapeLoadTime=EpicsSignal('MEC:LAS:FLOAT:12');
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

    #last updated 20230316
    EcoeffYFEval=0.3285 #was .3578
    Ecoeff1in1wCDval=0.7082#0.7718 # 0.5871
    Ecoeff2in1wABval=108.4#224.0
    Ecoeff2in1wEFval=199.2#177.5
    Ecoeff2in1wGHval=317.9#260.9826
    Ecoeff2in1wIJval=94.6#113.2
    Ecoeff2in2wABval=141.6#138.4#143.0#158.9#134.0135*0.9250
    Ecoeff2in2wEFval=152.2#121.2#174.7#165.2398*0.9978
    Ecoeff2in2wGHval=205.6#121.9#142.6#127.1#172.3#194.1412*1.0653
    Ecoeff2in2wIJval=175.9#142.4#145.3#153.0#147.7#156.9307*0.9198
    
    EcoeffRE1 = 1.64e5
    EcoeffRE0 = 1.03156061e-01 
    EcoeffTO1 = 3.48e7
    EcoeffTO0 = - 1.63e1 
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
    YSSS=EpicsSignal('MEC:LAS:ARRAY:10')

    
    
    
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
    
    
    EVRLPLLAMPEC=EpicsSignal(read_pv='MEC:LAS:EVR:01:TRIG7:EC_RBV',write_pv='MEC:LAS:EVR:01:TRIG7:TEC')#LPL lamp event code; needs 182
    EVRLPLLAMPEN=EpicsSignal('MEC:LAS:EVR:01:TRIG7:TCTL') #LPL lamp enable;
    EVRLPLSSEC=EpicsSignal(read_pv='MEC:LAS:EVR:01:TRIG8:EC_RBV', write_pv='MEC:LAS:EVR:01:TRIG8:TEC')#LPL slicer event code; needs 182 or 43 typically
    EVRLPLSSEN=EpicsSignal('MEC:LAS:EVR:01:TRIG8:TCTL') #LPL slicer enable;
    
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

    LMapAB=[5,100] #Pixel mapping from LeCroyA and LeCroyB horizontal axis (1002px) to Highland (140px)
    LMap2=[50,1000] #Pixel mapping from LeCroy2 horizontal axis (10002px) to Highland (140px)
    pwttfmap=[25,500]
    PSFILEPATH='/reg/neh/operator/mecopr/mecpython/pulseshaping/'
    
    HIGHLAND_IP = 'highland-mec-01'
    LECROY_A_IP = '172.21.46.100' #permanent: '172.21.46.100'#'scope-ics-meclas-lecroy-a'
    LECROY_B_IP = '172.21.46.100'###TEMP SUBSTITUTE### #permanent: '172.21.46.120'#'scope-ics-meclas-lecroy-b'#
    LECROY_1_IP = '172.21.46.60'#'scope-ics-mectc1-1'
    LECROY_2_IP = '172.21.46.128'#'scope-ics-meclas-lecroy-02'
    LECROY_L_IP = '172.21.160.252'#'scope-ics-meclas-lecroy-02'

    PDMAX_TEST = 28000
    
    YFE02MM_MAX = 100
    YFE06MM_MAX = 140
    YFE10MM_MAX = 150
    YFE02MM_SETPT = 90
    YFE06MM_SETPT = 135
    YFE10MM_SETPT = 130
    
    
    
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
        #12: DESC is CurrShape
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

        efc.wPV('MEC:LAS:FLOAT:31', cls.EcoeffYFEval);
        efc.wPV('MEC:LAS:FLOAT:32', cls.Ecoeff1in1wCDval);
        efc.wPV('MEC:LAS:FLOAT:33', cls.Ecoeff2in1wABval);
        efc.wPV('MEC:LAS:FLOAT:34', cls.Ecoeff2in1wEFval);
        efc.wPV('MEC:LAS:FLOAT:35', cls.Ecoeff2in1wGHval);
        efc.wPV('MEC:LAS:FLOAT:36', cls.Ecoeff2in1wIJval);
        efc.wPV('MEC:LAS:FLOAT:37', cls.Ecoeff2in2wABval);
        efc.wPV('MEC:LAS:FLOAT:38', cls.Ecoeff2in2wEFval);
        efc.wPV('MEC:LAS:FLOAT:39', cls.Ecoeff2in2wGHval);
        efc.wPV('MEC:LAS:FLOAT:40', cls.Ecoeff2in2wIJval);
        
        efc.wPV('MEC:LAS:ARRAY:01.DESC', 'Psns pulse segment lengths:10')
        efc.wPV('MEC:LAS:ARRAY:02.DESC', 'SSs pulse segment endpoint pairs:20')
        efc.wPV('MEC:LAS:ARRAY:03.DESC', 'Highland Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:04.DESC', 'YFE Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:05.DESC', 'YFEgoal Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:06.DESC', '1in1w Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:07.DESC', '2in1w Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:08.DESC', '2in2w Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:09.DESC', '2in2wgoal Grafana:140')
        efc.wPV('MEC:LAS:ARRAY:10.DESC', 'YSSs pulse segment endpoint pairs:20(1000)')
        efc.wPV('MEC:LAS:ARRAY:11.DESC', 'Spare Grafana:1000')
        efc.wPV('MEC:LAS:ARRAY:12.DESC', 'Spare Grafana:1000')
        efc.wPV('MEC:LAS:ARRAY:13.DESC', 'Spare Grafana:1000')
        efc.wPV('MEC:LAS:ARRAY:14.DESC', 'Spare Grafana:1000')


class Ben:#routines for LaserNet Ben
    """
    Routines for LaserNet Ben
    """
        
    def save_stuff(Lecroy_scope='A',chan_to_eLog=4):
        """
        Saves LecroyA Ch.4
        Saves camera LPL_100J_1
        """
        ExpName=LPL.get_curr_exp()
        RunNumber=LPL.get_curr_run()
        print(RunNumber)
        path_extension='/lecroy/'
        fpQ=str('/reg/neh/operator/mecopr/experiments/'+ExpName+path_extension)
        chan_to_eLog = int(chan_to_eLog)
        if chan_to_eLog not in [1,2,3,4]:
            print('Channel must be 1, 2, 3, or 4! Using channel 4 as default.')
            chan_to_eLog=4;
        if not os.path.exists(fpQ[-len(path_extension)]):
            print('File path '+fpQ[-len(path_extension)]+' does not exist! Trying to create it...')
            try:
                os.makedirs(fpQ[-len(path_extension)]);print('Folder created successfully!');
                os.chmod(fpQ[-len(path_extension)],stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO);
            except:
                print('Failed to create '+fpQ[-len(path_extension)]+'!')
        if not os.path.exists(fpQ):
            print('File path '+fpQ+' does not exist! Trying to create it...')
            try:
                os.makedirs(fpQ);print('Folder created successfully!');
                os.chmod(fpQ,stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO);
            except:
                print('Failed to create '+fpQ+'!')
        try:
            mecel = elog.ELog({'experiment':ExpName},user='mecopr',pw=pickle.load(open(GLOBAL.PSFILEPATH+'elogauth.p','rb')))
            #print('Connected to eLog of current MEC experiment, which is: '+ExpName)
        except:
            print('Failed to connect to eLog!')
            
        timestamp=datetime.now().strftime('%Y%m%d_%H%M%S')
            
        try:
            chdata=LOSC(Lecroy_scope).rchallxy();
        except:
            print('Failed to read out data!')
        np.savetxt(fpQ+'lecroyA_ch'+str(4)+'_'+timestamp+'.txt',np.c_[chdata[3][0], chdata[3][1]])
        ep.lxysav(chdata[chan_to_eLog-1][0],chdata[chan_to_eLog-1][1],fpQ+'lecroyA_ch'+str(4)+'_'+timestamp+'.png',abs_path=True)
        fullmsg=str('Scope trace data for Ch. 4 saved to '+fpQ+' with time stamp '+timestamp+' from run '+str(RunNumber)+'. Attached are the data and plot files for channel '+str(chan_to_eLog)+'.')
        #eplxy(chdata[chan_to_eLog-1][0],chdata[chan_to_eLog-1][1])
        try:
            CAM.QuickSave("custom:MEC:GIGE:36", ImageNo=1, AbsFilePathQ=fpQ)
            CAM.QuickSaveData("custom:MEC:GIGE:36", ImageNo=1, AbsFilePathQ=fpQ)
        except:
            print('Failed to save camera data!')

        try:
            mecel.post(fullmsg,attachments=[fpQ+'lecroyA_ch'+str(4)+'_'+timestamp+'.txt', fpQ+'lecroyA_ch'+str(4)+'_'+timestamp+'.png'], tags=['scope_trace'])
            print('Auto-saved to eLog.') 
        except:
            print('Failed to auto-save to eLog!')        
        
        
        
        
