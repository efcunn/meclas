# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:25:26 2026

@author: efcunn
"""
import socket
import time
import numpy as np
import struct

class LOSC:
    """
    Class containing all the necessary functions for running the LeCroy oscilloscopes
    Because we have several such scopes, instantiation of a certain device is required
    Possible scope choices are:
    my_scope=LOSC(myIP) #for a custom LeCroy with specified IP address
    
    Unless speed is necessary, it is usually most appropriate to interface with a LeCroy simply by using LOSC(myIP).[Command]
    This will take care of all of the socket opening/closing by itself.
    Example: read out all the current waveform amplitudes on scope myIP1 using wvfm4=LOSC(myIP1).rchall()
    Example: wait for a fresh acquisition before reading out channel 2's voltage vs time on scope myIP2 using ch2_wvfm=LOSC(myIP2).waitrchxy(2)
    (Alternatively, use the approach above: my_scope = LOSC(myIPx) and then wvfm4=my_scope.rchall() or whatever)
    
    Possible commands that can be used as described above include:
        :waitrch #wait and read specified channel amplitude
        :waitrchxy #wait and read specified channel amplitude and time
        :rch #immediately read specified channel amplitude
        :rchxy #immediately read specified channel amplitude and time
        :rchall #immediately read amplitude for all channels
        :rchallxy #immediately read amplitude and time for all channels
        :RestoreConfig #restore acquisition settings according to internal memory
        :WordFormatReset #enforce 16-bit readout of scope waveform, not janky 8-bit :D
    There are also functions available for use with bare sockets; these tend to start with underscores.
    See the docstrings for more guidance.
    
    Potential future work:
        : adding more functionality from the LeCroy manual using the _ctrl function
    """
    def __init__(self, IPAddrStrQ):
        """
        Initialize a LeCroy oscilloscope for use by passing in the scope's IP address

        """
        self._hostIP = IPAddrStrQ
        self._name = 'LeCroy'+str(IPAddrStrQ)
        self._LSock = None
        self._port = 1861

    def _Open(self):
        """
        Takes care of opening the socket to the specified LeCroy; if called explicitly like this, 
        it MUST be followed by a _Close() statement or else you'll block the socket and need to 
        locally disable/enable its networking card (or power cycle the unit, but this is not preferred!!)
        
        Using this function allows one to leave the socket open, which allows for quicker access to scope functions.
            :Example: after my_scope = LOSC(myIP) and my_scope._Open() then one may use functions inside a loop, e.g.
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
        Example: my_scope=LOSC(myIP)
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
                rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL,WORD".format(str(ChannelNo)))
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
                rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL,WORD".format(str(ChannelNo)))
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
        rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL,WORD".format(str(OChan)))
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
        rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL,WORD".format(str(OChan)))
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
            (note: administratively this is the memory file where one may keep the latest config)
        Example: LOSC(myIP).RestoreConfig() resets that scope's setup according to internal memory file #1
        This is especially useful if scope settings were mistakenly changed and need to be restored
            in order to allow for appropriate pulse shaping performance
        """
        LData=self._FunctionWrapper(self._ctrl,{'msg':'*RCL 1','SendOnly':False});
        return LData
    
    def WordFormatReset(self):
        """
        Function for avoiding the low-res readout problem when the scopes give a janky, 8-bit version of the desired waveform; sets the scope readout to 16-bit WORD
        """
        self._FunctionWrapper(self._ctrl,{'msg':"COMM_FORMAT DEF9,WORD,BIN",'SendOnly':True})
        return

