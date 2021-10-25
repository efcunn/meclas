import socket
import time
import math
import numpy as np
import struct
from scipy import signal
from scipy import stats
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
import pcdsdaq.ext_scripts ##find where used -- get_run_number and experiment?
import glob
import pandas as pd
import stat
import getpass
import multiprocessing
import termios, tty
import importlib.util


class efc:
    def dotsleep(tSEC):
        for ii in range(tSEC):
            print('.',end='',flush=True);time.sleep(1);
        print('*')
        return
        
    def IndFETWave(ListOfPixels,WriteValue):#only used in FET survey?? put inside HAWG??
        itt=0
        NewString=''
        while itt<140:
            if (itt+1) in ListOfPixels:
                NewString+=Hex2Byte(WriteValue)
            else:
                NewString+='0000'
            itt+=1
        return NewString

    def LinearWave(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height):
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
            return LinearWave(1,0,140,0)
        #
        while itt<140:
            if itt<p1:
                NewString+='0000'
            elif p1<=itt<=p2:
                NewString+=Hex2Byte(int(h2+((itt-p2)*(h2-h1)/float(p2-p1))))
            else:
                NewString+='0000'
            itt+=1
        return NewString

    def LinearWave2(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height,offsetQ,arraylenQ):
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
                nextval=int((h1*(itt-p2)*(itt-p3)/float((p2-p1)*(p3-p1)))+(h2*(itt-p1)*(itt-p3)/float((p2-p1)*(p2-p3)))+(h3*(itt-p1)*(itt-p2)/float((p3-p1)*(p3-p2))))
                NewList.append(nextval)
            else:# itt>p2:
                NewList.append(offsetQ-offsetQ)
            itt+=1
        return np.array(NewList)+offsetQ

    def ExponentialWave2(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height,offsetQ,arraylenQ):
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
                NewList.append(int(((h1-h2)/math.log(((p1+.0001)/float(p2)),LogBase))*math.log((LogBase**(((h2*math.log(p1+.0001,LogBase))-(h1*math.log(p2+.0001,LogBase)))/float(h1-h2)))*(itt+.0001),LogBase)))
            else:# itt>p2:
                NewList.append(offsetQ-offsetQ)
            itt+=1
        return np.array(NewList)+offsetQ

    def epl(listq):
        df1=plt.figure()
        plt.plot(listq);
        df1.show()
        return

    def eplxy(listxq,listyq):
        df1=plt.figure()
        plt.plot(listxq,listyq);
        df1.show()
        return

    def eplxyloglog(listxq,listyq):
        df1=plt.figure()
        plt.loglog(listxq,listyq);
        df1.show()
        return

    def eplsav(listq,FileNameQ,blockdisplay=True):
        df1=plt.figure()
        plt.plot(listq);
        df1.savefig(str(FileNameQ+'.png'))
        if blockdisplay:
            plt.close(df1)
        return
        
    def eplxysav(listxq,listyq,FileNameQ,abs_path=False):
        df1=plt.figure()
        plt.plot(listxq,listyq);
        if abs_path:
            figfilename=FileNameQ;
        else:
            figfilename=str(psfilepath()+FileNameQ+'.png')
        df1.savefig(figfilename)        
        plt.close(df1)
        return
        
    def eplcomp(listq,goalq,Map,tMax):
        formtra=[]
        formtra.append(TraceFormatting(listq,Map,tMax))
        formtra.append(goalq)
        epll(formtra)
        return


    def eplcsv(CSVname):
        with open(psfilepath()+'data/'+CSVname+'.csv','r') as filehead:
            RawListQ=filehead.read()
            ListedValues=RawListQ.split('\n')
        epl(ListedValues[:-1])
        return 

    def epllcsv(CSVHeadname):
        ListofListedValues=[]
        for ii in range(1,5):
            with open(psfilepath()+'data/'+CSVHeadname+'_ch'+str(ii)+'.csv','r') as filehead:
                RawListQ=filehead.read()
                ListedValues=RawListQ.split('\n')
            ListofListedValues.append(ListedValues[:-1])
        epll(ListofListedValues)
        return 

    def rcsv(CSVname):
        with open(psfilepath()+'data/'+CSVname+'.csv','r') as filehead:
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
        df1=plt.figure()
        for ii in range(len(llist)):
            plt.plot(llist[ii]);
        df1.show()
        
    def epllxy(llistxyq,xlb='none',ylb='none'):
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
        df1=plt.figure()
        for ii in range(len(llistxyq)):
            plt.loglog(llistxyq[ii][0],llistxyq[ii][1]);
        df1.show()
        return

    def epllt(listq,Map):
        formtra=[]
        for ii in range(len(listq)):
            formtra.append(TraceFormatting(listq[ii],Map[ii],1))
        epll(formtra)
        return
        
    def epllcomp(listq,goalq,Map,tMax):
        formtra=[]
        formtra.append(TraceFormatting(listq[-1],Map,tMax))
        formtra.append(goalq)
        epll(formtra)
        return
        
    def sumch(chnoQ,sockQ):
        templistQ=rch(chnoQ,sockQ)
        return np.sum(templistQ)
        
    def reloadchk():
        print('Last stamped: 20211020bbbb')
        
    def reloadpkg(pkgname):
        #spec = importlib.util.spec_from_file_location("mecps4", "/reg/neh/operator/mecopr/mecpython/pulseshaping/mecps4.py")
        #mecps4 = importlib.util.module_from_spec(spec)
        #spec.loader.exec_module(mecps4);
        importlib.reload(pkgname)

    def YFEtrace():#fix this
        try:
            SLA=LXOpen('A');time.sleep(0.15);
            pch(1,SLA);time.sleep(0.15);
            LXClose(SLA);time.sleep(0.15);
        except:
            LXClose(SLA);
            print('Failed to display trace')
            return False






class HAWG:
    def __init__(self):
        self._HighlandSocket=None
        self._HIGHLAND_SLAVE_ADDRESS=0 #arbitrary, I think    

    def Open(self):  
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
        try:
            self._HighlandSocket.close()
            self._HighlandSocket=None
        except:
            print('Unable to close socket -- it may already be closed')
        return
        
    def Reset(self):
        pvAWG=EpicsSignal('MEC:64B:PWR:2:Outlet:1:SetControlAction')#read AND write:1=ON,2=OFF
        print('Powering off Highland AWG, waiting 3sec...',end='',flush=True);
        pvAWG.put(2);
        dotsleep(3);
        print('Rebooting Highland AWG, waiting 10sec...',end='',flush=True)
        pvAWG.put(1);
        dotsleep(10);
        
    @staticmethod
    def _Hex1Byte(num):
        return '{0:02x}'.format(int(num)%(0xff+1))
        
    @staticmethod
    def _Hex2Byte(num):
        return '{0:04x}'.format(int(num)%(0xffff+1))

    @staticmethod
    def _ByteSum(Datastr): #byte-wise addition
        #accepts string in 'xx' format, e.g. '1e....'
        #inmust be in format returned by hexlify(data)
        #for REPLY CKS, input should use
        bytesum=0
        for byte in range(len(Datastr)//2):
            bytesum+=int(Datastr[2*byte:2*byte+2],16)
        return bytesum #returns an integer

    @classmethod
    def _PollConstructor(cls,COMMAND_CODE,POLL_LENGTH,SLAVE_ADDRESS,DATA):
        #enter COMMAND_CODE, POLL_LENGTH, SLAVE_ADDRESS as integer values
        #enter DATA as a string, e.g. 'ffff0000' or empty string '' for no data 
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
        #input REPLY_STRING already formatted using hexlify
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
        MyPollQ=cls._PollConstructor(COMMAND_CODE,POLL_LENGTH,SLAVE_ADDRESS,DATA)
        MySocketQ.send(MyPollQ)
        MyRawReplyQ=MySocketQ.recv(REPLY_LENGTH)
        HStatusQ, HDataQ, HErrorQ = cls._ReplyInterpreter(REPLY_LENGTH,SLAVE_ADDRESS,hexlify(MyRawReplyQ))
        return HStatusQ, HDataQ, HErrorQ

    @staticmethod
    def _StatusInterpreter(HError, HStatus, Quiet=True):
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
        try:
            HStatusQ, HDataQ, HErrorQ = cls._SendPollRecvReply(
            MySocketQ,CommandCode,PollLength,ReplyLength,SlaveAddress,MyData)
            cls._StatusInterpreter(HErrorQ, HStatusQ)
            return HDataQ
        except:
            print('Failed!')
            return False

    def _FunctionWrapper(self,FuncQ,kwargs={}):
        try:
            self.Open();time.sleep(0.15);
            HDataQList=FuncQ(**kwargs);time.sleep(0.15);
            self.Close();time.sleep(0.15);
        except:
            self.Close();time.sleep(0.15);
        return HDataQList

    def _ReadStatus(self):
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
        HDataQ=self._FunctionWrapper(self._ReadStatus);
        return HDataQ

    def _ClearStatus(MySocketQ, SlaveAddress):
        print('**CLEAR STATUS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=1,PollLength=9, ReplyLength=9,MyData='')
        print('****')
        return HDataQ
    def ClearStatus(self):
        HDataQ=self._FunctionWrapper(self._ClearStatus);
        return HDataQ

    def _ReadPulseHeights(self,ShowPlot=False):
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=2,PollLength=9, ReplyLength=289,MyData='')
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(int(len(HDataQ)//4))]
        if ShowPlot:
            efc.epl(HDataQList)
        return HDataQList
    def ReadPulseHeights(self,ShowPlot=False):
        HDataQList=self._FunctionWrapper(self._ReadPulseHeights,{'ShowPlot':ShowPlot});
        return HDataQList
        
    def _WritePulseHeights(self, FileNameOrStringOrList=140*[0]):
        MyDataQ=''
        if len(FileNameOrStringOrList) == 140*4:#will accept pre-formatted Hex2Byte text
            MyDataQ=FileNameOrStringOrList
        elif len(FileNameOrStringOrList)==140:#will accept a straight list
            for value in range(len(FileNameOrStringOrList)):
                MyDataQ+=Hex2Byte(int(FileNameOrStringOrList[value]))
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
                MyDataQ+=Hex2Byte(int(ListedValues[value]))
        else:
            print('Bad file entry count: '+str(len(FileNameOrStringOrList)))
            return False
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=3,PollLength=289, ReplyLength=9,MyData=MyDataQ)
        return HDataQ
    def WritePulseHeights(self,FileNameOrStringOrList=140*[0]):
        HDataQ=self._FunctionWrapper(self._WritePulseHeights,{'FileNameOrStringOrList':FileNameOrStringOrList});
        return HDataQ
        
    def _ReadFiducialImpulseSettings(self):
        print('**READ FIDUCIAL IMPULSE SETTINGS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=4,PollLength=9, ReplyLength=13,MyData='')
        print('Fiducial pulse height (max 65535): '+str(int(HDataQ[:4],16)))
        print('Fiducial pulse delay: '+str(int(HDataQ[4:8],16)))
        print('****')
        return HDataQ
    def ReadFiducialImpulseSettings(self):
        HDataQ=self._FunctionWrapper(self._ReadFiducialImpulseSettings);
        return HDataQ        
        
    def _WriteFiducialImpulseSettings(self,AmpReq=0,TimeReq=0):
        print('**WRITE FIDUCIAL IMPULSE SETTINGS**')
        MyDataQ=''
        MyDataQ+=Hex2Byte(int(AmpReq))
        MyDataQ+=Hex2Byte(int(TimeReq))
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=5,PollLength=13, ReplyLength=9,MyData=MyDataQ)
        print('****')
        return HDataQ
    def WriteFiducialImpulseSettings(self,AmpReq=0,TimeReq=0):
        HDataQ=self._FunctionWrapper(self._WriteFiducialImpulseSettings,{'AmpReq':AmpReq,'TimeReq':TimeReq});
        return HDataQ
        
    def _WriteEnableByte(self,EnableTotal=4):
        #to determine EnableTotal input, start from 0 and:
        #+1 for Enable CPU self-trigger, 960 Hz (test mode)
        #+2 for Enable self-trigger, 20 kHz
        #+4 for Enable external triggers
        #+8 for Enable the BIAS generators
        MyDataQ=''
        MyDataQ+=Hex1Byte(EnableTotal)
        print('**WRITE ENABLE BYTE**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=10,PollLength=10, ReplyLength=9,MyData=MyDataQ)
        print('****')
        return HDataQ
    def WriteEnableByte(self,EnableTotal=4):
        HDataQ=self._FunctionWrapper(self._WriteEnableByte,{'EnableTotal':EnableTotal});
        return HDataQ
        
    def _ReadT0Delay(self):
        print('**READ T0 DELAY**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=13,PollLength=9, ReplyLength=11,MyData='')
        print('T0 delay (max 50000 (50ns)): '+str(int(HDataQ,16)))
        print('****')
        return int(HDataQ,16)
    def ReadT0Delay(self):
        HDataQ=self._FunctionWrapper(self._ReadT0Delay);
        return HDataQ
        
    def _ReadWaveAmplitudeCalibrations(self):
        print('**READ WAVE AMPLITUDE CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=15,PollLength=9, ReplyLength=289,MyData='')
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
        print('Wave amplitude calibration (nom 2800): '+str(HDataQList))
        print('****')
        return HDataQList
    def ReadWaveAmplitudeCalibrations(self):
        HDataQList=self._FunctionWrapper(self._ReadWaveAmplitudeCalibrations);
        return HDataQList
        
    def _WriteWaveAmplitudeCalibrations(self, StringOrList):
        MyDataQ=''
        if len(StringOrList) == 140*4:#will accept pre-formatted Hex2Byte text
            MyDataQ=StringOrList
        elif len(StringOrList)==140:#will accept a straight list
            for value in range(len(StringOrList)):
                MyDataQ+=Hex2Byte(int(StringOrList[value]))
        else:
            print('Bad file entry count: '+str(len(StringOrList)))
            return 
        print('**WRITE WAVE AMPLITUDE CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=16,PollLength=289, ReplyLength=9,MyData=MyDataQ)
        print('****')
        return HDataQ
    def WriteWaveAmplitudeCalibrations(self,StringOrList):
        HDataQ=self._FunctionWrapper(self._WriteWaveAmplitudeCalibrations,{'StringOrList':StringOrList});
        return HDataQ
        
    def _ReadWaveTimeCalibrations(self):
        print('**READ WAVE TIME CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=17,PollLength=9, ReplyLength=289,MyData='')
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
        print('Wave time calibrations (max 65535): '+str(HDataQList))
        print('****')
        return HDataQList
    def ReadWaveTimeCalibrations(self):
        HDataQList=self._FunctionWrapper(self._ReadWaveTimeCalibrations);
        return HDataQList  
        
    def _WriteWaveTimeCalibrations(self, StringOrList):
        MyDataQ=''
        if len(StringOrList) == 140*4:#will accept pre-formatted Hex2Byte text
            MyDataQ=StringOrList
        elif len(StringOrList)==140:#will accept a straight list
            for value in range(len(StringOrList)):
                MyDataQ+=Hex2Byte(int(StringOrList[value]))
        else:
            print('Bad file entry count: '+str(len(StringOrList)))
            return 
        print('**WRITE WAVE TIME CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=18,PollLength=289, ReplyLength=9,MyData=MyDataQ)
        print('****')
        return HDataQ
    def WriteWaveTimeCalibrations(self,StringOrList):
        HDataQ=self._FunctionWrapper(self._WriteWaveTimeCalibrations,{'StringOrList':StringOrList});
        return HDataQ

    def _ReadMiscellaneousCalibrations(self):
        print('**READ MISCELLANEOUS CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=19,PollLength=9, ReplyLength=81,MyData='')
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
        print('Miscellaneous calibrations: '+str(HDataQList))
        print('****')
        return HDataQList
    def ReadMiscellaneousCalibrations(self):
        HDataQList=self._FunctionWrapper(self._ReadMiscellaneousCalibrations);
        return HDataQList  
            
    def _WriteMiscellaneousCalibrations(self,StringOrList):
        MyDataQ=''
        if len(StringOrList) == 36*4:#will accept pre-formatted Hex2Byte text
            MyDataQ=StringOrList
        elif len(StringOrList)==36:#will accept a straight list
            for value in range(len(StringOrList)):
                MyDataQ+=Hex2Byte(int(StringOrList[value]))
        else:
            print('Bad file entry count: '+str(len(StringOrList)))
            return
        print('**WRITE MISCELLANEOUS CALIBRATIONS**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=20,PollLength=81, ReplyLength=9,MyData=MyDataQ)
        print('****')
        return 
    def WriteMiscellaneousCalibrations(self,StringOrList):
        HDataQ=self._FunctionWrapper(self._WriteMiscellaneousCalibrations,{'StringOrList':StringOrList});
        return HDataQ
    
    def _ReadWalkTable(self):
        print('**READ WALK TABLE**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=25,PollLength=9, ReplyLength=73,MyData='')
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
        print('Walk Table: '+str(HDataQList))
        print('****')
        return HDataQList    
    def ReadWalkTable(self):
        HDataQList=self._FunctionWrapper(self._ReadWalkTable);
        return HDataQList  
        
    def _WriteWalkTable(self,StringOrList):
        MyDataQ=''
        if len(StringOrList) == 32*4:#will accept pre-formatted Hex2Byte text
            MyDataQ=StringOrList
        elif len(StringOrList)==32:#will accept a straight list
            for value in range(len(StringOrList)):
                MyDataQ+=Hex2Byte(int(StringOrList[value]))
        else:
            print('Bad file entry count: '+str(len(StringOrList)))
            return
        print('**WRITE WALK TABLE**')
        HDataQ=self._BaseFunction(self._HighlandSocket,self._HIGHLAND_SLAVE_ADDRESS,CommandCode=26,PollLength=73, ReplyLength=9,MyData=MyDataQ)
        HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
        print('****')
        return
    def WriteWalkTable(self,StringOrList):
        HDataQ=self._FunctionWrapper(self._WriteWalkTable,{'StringOrList':StringOrList});
        return HDataQ
        













class LOSC:
    def __init__(self, LStrQ):
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
        try:
            self._LSock.close()
            self._LSock=None
        except:
            print('Unable to close socket -- it may already be closed')
        return

    def _send_and_reply(self,msg,SendOnly=False):
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
        try:
            self.Open();time.sleep(0.15);
            LData=FuncQ(**kwargs);time.sleep(0.15);
            self.Close();time.sleep(0.15);
        except:
            self.Close();time.sleep(0.15);
        return LData

    def _waitrch(self,ChannelNo,verbose=False):
        while True:
            ready = (int(self._send_and_reply("INR?").split()[1]) & 1) == 1
            if ready:
                rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL".format(str(ChannelNo)))
                fullaq = self._parsewf(rawdataq, verbose)
                return fullaq
    def waitrch(self,ChannelNo,verbose=False):
        LData=self._FunctionWrapper(self._waitrch,{'ChannelNo':ChannelNo,'verbose':verbose});
        return LData
        
    def _waitrchxy(self,ChannelNo,verbose=False):
        while True:
            ready = (int(self._send_and_reply("INR?").split()[1]) & 1) == 1
            if ready:
                rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL".format(str(ChannelNo)))
                fullaq = self._parsewf(rawdataq, verbose)
                yvals=fullaq['DATA'];xvals=[fullaq['HORIZ_OFFSET'] + ii*fullaq['HORIZ_INTERVAL'] for ii in range(len(fullaq['DATA']))];
                return [xvals,yvals]
    def waitrchxy(self,ChannelNo,verbose=False):
        LData=self._FunctionWrapper(self._waitrchxy,{'ChannelNo':ChannelNo,'verbose':verbose});
        return LData
                
    def _rch(self,OChan):
        rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL".format(str(OChan)))
        return self._parsewf(rawdataq, False)['DATA']
    def rch(self,OChan):
        LData=self._FunctionWrapper(self._rch,{'OChan':OChan});
        return LData

    def _rchxy(self,OChan):
        rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL".format(str(OChan)))
        fullaq=self._parsewf(rawdataq, False);
        yvals=fullaq['DATA'];xvals=[fullaq['HORIZ_OFFSET'] + ii*fullaq['HORIZ_INTERVAL'] for ii in range(len(fullaq['DATA']))];
        return [xvals,yvals]
    def rchxy(self,OChan):
        LData=self._FunctionWrapper(self._rchxy,{'OChan':OChan});
        return LData
        
    def _rchall(self):
        rchans=[]
        for OChan in range(1,5):
            rchans.append(self._rch(OChan))
        return rchans
    def rchall(self):
        LData=self._FunctionWrapper(self._rchall);
        return LData
                
    def _rchallxy(self):
        rchans=[]
        for OChan in range(1,5):
            rchans.append(self._rchxy(OChan))
        return rchans
    def rchallxy(self):
        LData=self._FunctionWrapper(self._rchallxy);
        return LData

    def _sch(self,OChan,FileName):
        rawdataq = self._send_and_reply("C{}:WAVEFORM? ALL".format(OChan))
        parseddataq=self._parsewf(rawdataq, False)
        with open(psfilepath()+'data/'+str(FileName)+'.csv','w',newline='') as f:
            writer=csv.writer(f, delimiter='\n')
            writer.writerow(parseddataq['DATA'])
        with open(psfilepath()+'data/'+str(FileName)+'-h.csv','w',newline='') as f:
            writer=csv.DictWriter(f, parseddataq.keys())
            writer.writeheader()
            writer.writerow(parseddataq)
        return 
    def sch(self,OChan,FileName):
        LData=self._FunctionWrapper(self._sch,{'OChan':OChan,'FileName':FileName});
        return LData
                
    def _schall(self,FileName):
        for OChan in range(1,5):
            self._sch(OChan,FileName+'_ch'+str(OChan))
        return
    def schall(self,FileName):
        LData=self._FunctionWrapper(self._schall,{'FileName':FileName});
        return LData

    def _pch(self,OChan):
        pdata=self._rch(OChan)
        efc.epl(pdata)
        return 
    def pch(self,OChan):
        LData=self._FunctionWrapper(self._pch,{'OChan':OChan});
        return LData
        
    def _pchxy(self,OChan):
        pdata=self._rchxy(OChan)
        efc.eplxy(*pdata)
        return 
    def pchxy(self,OChan):
        LData=self._FunctionWrapper(self._pchxy,{'OChan':OChan});
        return LData
        
    def _pchall(self):
        pdata=self._rchall()
        efc.epll(pdata)
        return
    def pchall(self):
        LData=self._FunctionWrapper(self._pchall);
        return LData
        
    def _pchallxy(self):
        pdata=self._rchallxy()
        efc.epllxy(pdata)
        return
    def pchallxy(self):
        LData=self._FunctionWrapper(self._pchallxy);
        return LData
        
    def _ctrl(msg):
        #(msg='TIME_DIV?',SendOnly=False)
        #(msg='*RCL 3',SendOnly=False)
        #(msg='TDIV 100E-9',SendOnly=True)
        #(msg='C1:VDIV 500E-3',SendOnly=True)
        #(msg=r"""vbs 'app.acquisition.triggermode = "single" ' """,SendOnly=True)
        #(msg=r"""vbs 'app.acquisition.c1.deskew = 0 ' """,SendOnly=True)
        self._send_and_reply(msg,SendOnly=True)
