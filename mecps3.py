#((once finished, rename as mecps in folder and then transfer to final location))
# -*- coding: utf-8 -*-
###conglomerate of different old files

###shape2.py

import socket
import time
import math
from binascii import hexlify

from binascii import unhexlify
#for the LeCroy
import struct
import numpy as np
#

from scipy import signal
##from scipy.stats import threshold
from scipy import stats
import matplotlib.pyplot as plt
import pickle
import time
from datetime import date, datetime
import math
#for Highland stuff
from binascii import hexlify
from binascii import unhexlify
#for the LeCroy
import struct
import numpy as np
#
import csv
#for checking for file existence
import os.path

import sys
from ophyd.signal import EpicsSignal
import elog
#from pcdsdaq.ext_scripts import get_run_number
from pcdsdaq.ext_scripts import *

import glob
import pandas as pd
import stat
import getpass
import multiprocessing
import termios, tty

#
def Hex1Byte(num):
    return '{0:02x}'.format(int(num)%(0xff+1))
#
def Hex2Byte(num):
    return '{0:04x}'.format(int(num)%(0xffff+1))
#
def ByteSum(datastr): #byte-wise addition
    #accepts string in 'xx' format, e.g. '1e....'
    #input must be in format returned by hexlify(data)
    #for REPLY CKS, input should use
    bytesum=0
    for byte in range(len(datastr)//2):
        bytesum+=int(datastr[2*byte:2*byte+2],16)
    return bytesum #returns an integer
#
def PollConstructor(COMMAND_CODE,POLL_LENGTH,SLAVE_ADDRESS,DATA):
    #enter COMMAND_CODE, POLL_LENGTH, SLAVE_ADDRESS as integer values
    #enter DATA as a string, e.g. 'ffff0000' or empty string '' for no data 
    ProtoCommand=''
    ProtoCommand+='0B'#MSYN: all commands begin with this byte (0B)
    ProtoCommand+=Hex2Byte(POLL_LENGTH) #BC:BC: number of bytes in message
    ProtoCommand+=Hex2Byte(SLAVE_ADDRESS) #RA:RA: slave address
    ProtoCommand+=Hex1Byte(COMMAND_CODE) #CMD: command code
    ProtoCommand+=DATA #<data> must already be formatted properly 'xxxx' or ''
    BYTE_SUM=ByteSum(ProtoCommand) #compute the sum
    ProtoCommand+=Hex2Byte(BYTE_SUM) #CKS:CKS: 16-bit sum of all preceding bytes
    ProtoCommand+='17' #ETB: end of message byte, 17 hex
    Command=unhexlify(ProtoCommand)
    return Command
#
def ReplyInterpreter(REPLY_LENGTH,SLAVE_ADDRESS,REPLY_STRING):
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
    if ByteSum(REPLY_STRING[:-6])==int(REPLY_STRING[-6:-2],16): HError+='0'
    else: HError+='1' #checksum error
    if int(REPLY_STRING[-2:],16)==int('17',16): HError+='0'
    else: HError+='1' #wrong end-of-message byte, 17 hex
    return HStatus, HData, HError
#
def SendPollRecvReply(MySocketQ,COMMAND_CODE,POLL_LENGTH,REPLY_LENGTH,SLAVE_ADDRESS,DATA):
    MyPollQ=PollConstructor(COMMAND_CODE,POLL_LENGTH,SLAVE_ADDRESS,DATA)
    MySocketQ.send(MyPollQ)
    MyRawReplyQ=MySocketQ.recv(REPLY_LENGTH)
    HStatusQ, HDataQ, HErrorQ = ReplyInterpreter(REPLY_LENGTH,SLAVE_ADDRESS,hexlify(MyRawReplyQ))
    return HStatusQ, HDataQ, HErrorQ
#
def StatusInterpreter(HError, HStatus, Quiet=True):
    if HError[0]=='1':
        print('WARNING: Wrong start-of-message byte received')
    if HError[1]=='1':
        print('WARNING: Reply length discrepancy')
    if HError[2]=='1':
        print('WARNING: Slave address not echoed')
    if HError[3]=='1':
        print('WARNING: Checksum error')
    if HError[4]=='1':
        print('WARNING: Wrong end-of-message byte received')
    if not Quiet:
        if int(HStatus,16)==0:
            print('STATUS: NORMAL')
        else:
            print('STATUS: ERROR FLAG(S) RECEIVED')
            if ((int(HStatus,16))&(2**(8-1)))!=0:
                print('-trigger/bias timing error')
            if ((int(HStatus,16))&(2**(8-3)))!=0:
                print('-backup RAM data/calibrations lost flag')
            if ((int(HStatus,16))&(2**(8-4)))!=0:
                print('-powerfail/restart flag')
            if ((int(HStatus,16))&(2**(8-7)))!=0:
                print('-trigger/bias timing error')
    return

##############################################################################################################################################################
#
def IndFETWave(ListOfPixels,WriteValue):
    itt=0
    NewString=''
    while itt<140:
        if (itt+1) in ListOfPixels:
            NewString+=Hex2Byte(WriteValue)
        else:
            NewString+='0000'
        itt+=1
    return NewString
#
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
#
def LinearWave2(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height,offsetQ,arraylenQ):
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
            nextval=h2+((itt-p2)*(h2-h1)/float(p2-p1))#h1*((h2/float(h1))**((itt-p1)/float(p2-p1)))
            NewList.append(nextval)
        else:# itt>p2:
            NewList.append(offsetQ-offsetQ)
        itt+=1
    return np.array(NewList)+offsetQ
#
def StepWave(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height,MidPixNo):
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
        elif p1<=itt<=int(MidPixNo):
            NewString+=Hex2Byte(h1)
        elif int(MidPixNo)+1<=itt<=p2:
            NewString+=Hex2Byte(h2)
        else:# itt>p2:
            NewString+='0000'
        itt+=1
    return NewString
#
def ParabolicWave(Edge1PixNo,Edge1Height,MidPixNo,MidHeight,Edge2PixNo,Edge2Height):
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
    if MidHeight>65535:
        print('Mid height exceeds max value of 65535')
        h3=65535
    elif MidHeight<0:
        print('Mid height must be positive')
        h3=0
    else:
        h3=int(MidHeight)
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
    if p1<MidPixNo<p2:
        p3=int(MidPixNo)
    else:
        print('Middle pixel must come be between two end pixels.')
        return LinearWave(1,0,140,0)
    #
    while itt<140:
        if itt<p1:
            NewString+='0000'
        elif p1<=itt<=p2:
            nextval=int((h1*(itt-p2)*(itt-p3)/float((p2-p1)*(p3-p1)))+(h2*(itt-p1)*(itt-p3)/float((p2-p1)*(p2-p3)))+(h3*(itt-p1)*(itt-p2)/float((p3-p1)*(p3-p2))))
            if nextval>65535:
                nextval=65535
            elif nextval<0:
                nextval=0
            else:
                pass
            NewString+=Hex2Byte(nextval)
        else:# itt>p2:
            NewString+='0000'
        itt+=1
    return NewString
#
def ExponentialWave(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height):
    itt=0
    NewString=''
    if Edge1Height>65535:
        print('Edge1 height exceeds max value of 65535')
        h1=65535
    elif Edge1Height<1:
        print('Edge1 height must be positive and nonzero')
        h1=1
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
            NewString+=Hex2Byte(int(h1*((h2/float(h1))**((itt-p1)/float(p2-p1)))))
        else:# itt>p2:
            NewString+='0000'
        itt+=1
    return NewString
#
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
#
def LogWave(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height,LogBase):
    itt=0
    NewString=''
    if Edge1Height>Edge2Height:
        print('Edge1 must be lower than Edge2')
        Dummy=int(Edge1Height)
        Edge1Height=int(Edge2Height)
        Edge2Height=Dummy
    if Edge1Height==Edge2Height:
        print('Edge1 must be different than Edge2')
        return LinearWave(Edge1PixNo,Edge1Height,Edge2PixNo,Edge2Height)
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
            NewString+=Hex2Byte(int(((h1-h2)/math.log(((p1+.0001)/float(p2)),LogBase))*math.log((LogBase**(((h2*math.log(p1+.0001,LogBase))-(h1*math.log(p2+.0001,LogBase)))/float(h1-h2)))*(itt+.0001),LogBase)))
        else:# itt>p2:
            NewString+='0000'
        itt+=1
    return NewString


##############################################################################################################################################################
#
def ReadStatus(MySocketQ, SlaveAddress):
    CommandCodeQ = 0
    PollLengthQ = 9
    ReplyLengthQ = 22
    MyDataQ = ''
    print('**READ STATUS**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
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
    return
#
def ClearStatus(MySocketQ, SlaveAddress):
    CommandCodeQ = 1
    PollLengthQ = 9
    ReplyLengthQ = 9
    MyDataQ = ''
    print('**CLEAR STATUS**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    print('****')
    return    
#
def ReadPulseHeights(MySocketQ, SlaveAddress):
    CommandCodeQ = 2
    PollLengthQ = 9
    ReplyLengthQ = 289
    MyDataQ = ''
    #print('**READ PULSE HEIGHTS**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(int(len(HDataQ)//4))]
    #print('Raw pulse heights: '+str(HDataQList))
    ##think about plotting?
    #print('****')
    return HDataQList #muted the printing for use in a script
#
def WritePulseHeights(MySocketQ, SlaveAddress, FileNameOrStringOrList):
    CommandCodeQ = 3
    PollLengthQ = 289
    ReplyLengthQ = 9
    #best formatting for data in? .txt or .csv file of values?
    #currently written to accept a single-row .txt/.csv of comma-separated values
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
                return
        if len(ListedValues) != 140:
            print('File must have 140 entries; entry count: '+str(len(ListedValues)))
            return
        for value in range(len(ListedValues)):
            MyDataQ+=Hex2Byte(int(ListedValues[value]))
    else:
        print('Bad file entry count: '+str(len(FileNameOrStringOrList)))
        return 
    #
    #print('**WRITE PULSE HEIGHTS**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    #StatusInterpreter(HErrorQ, HStatusQ)
    #print('****')
    return
#
def ReadFiducialImpulseSettings(MySocketQ, SlaveAddress):
    CommandCodeQ = 4
    PollLengthQ = 9
    ReplyLengthQ = 13
    MyDataQ=''
    print('**READ FIDUCIAL IMPULSE SETTINGS**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    print('Fiducial pulse height (max 65535): '+str(int(HDataQ[:4],16)))
    print('Fiducial pulse delay: '+str(int(HDataQ[4:8],16)))
    print('****')
    return
#    
def WriteFiducialImpulseSettings(MySocketQ, SlaveAddress, AmpReq, TimeReq):
    CommandCodeQ = 5
    PollLengthQ = 13
    ReplyLengthQ = 9
    MyDataQ=''
    print('**WRITE FIDUCIAL IMPULSE SETTINGS**')
    MyDataQ+=Hex2Byte(int(AmpReq))
    MyDataQ+=Hex2Byte(int(TimeReq))
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    print('****')
    return
#
def WriteEnableByte(MySocketQ, SlaveAddress, EnableTotal):
    #to determine EnableTotal input, start from 0 and:
    #+1 for Enable CPU self-trigger, 960 Hz (test mode)
    #+2 for Enable self-trigger, 20 kHz
    #+4 for Enable external triggers
    #+8 for Enable the BIAS generators
    CommandCodeQ = 10
    PollLengthQ = 10
    ReplyLengthQ = 9
    MyDataQ=''
    MyDataQ+=Hex1Byte(EnableTotal)
    print('**WRITE ENABLE BYTE**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    print('****')
    return
##############################
#
def ReadT0Delay(MySocketQ, SlaveAddress):
    CommandCodeQ = 13
    PollLengthQ = 9
    ReplyLengthQ = 11
    MyDataQ=''
    print('**READ T0 DELAY**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    print('T0 delay (max 50000 (50ns)): '+str(int(HDataQ,16)))
    print('****')
    return int(HDataQ,16)
#
def ReadWaveAmplitudeCalibrations(MySocketQ, SlaveAddress):
    CommandCodeQ = 15
    PollLengthQ = 9
    ReplyLengthQ = 289
    MyDataQ=''
    print('**READ WAVE AMPLITUDE CALIBRATIONS**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    #print(HDataQ)
    HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
    print('Wave amplitude calibration (nom 2800): '+str(HDataQList))
    print('****')
    return HDataQList
#    
def WriteWaveAmplitudeCalibrations(MySocketQ, SlaveAddress, StringOrList):
    CommandCodeQ = 16
    PollLengthQ = 289
    ReplyLengthQ = 9
    MyDataQ=''
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
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    #print(HDataQ)
    HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
    print('****')
    return 
#   
def ReadWaveTimeCalibrations(MySocketQ, SlaveAddress):
    CommandCodeQ = 17
    PollLengthQ = 9
    ReplyLengthQ = 289
    MyDataQ=''
    print('**READ WAVE TIME CALIBRATIONS**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    #print(HDataQ)
    HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
    print('Wave time calibrations (max 65535): '+str(HDataQList))
    print('****')
    return HDataQList
#
def WriteWaveTimeCalibrations(MySocketQ, SlaveAddress, StringOrList):
    CommandCodeQ = 18
    PollLengthQ = 289
    ReplyLengthQ = 9
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
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    #print(HDataQ)
    HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
    print('****')
    return 
#    
def ReadMiscellaneousCalibrations(MySocketQ, SlaveAddress):
    CommandCodeQ = 19
    PollLengthQ = 9
    ReplyLengthQ = 81 #not 79!!!
    MyDataQ=''
    print('**READ MISCELLANEOUS CALIBRATIONS**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
    print('Miscellaneous calibrations: '+str(HDataQList))
    print('****')
    return HDataQList
#    
def WriteMiscellaneousCalibrations(MySocketQ, SlaveAddress, StringOrList):
    CommandCodeQ = 20
    PollLengthQ = 81 #not 79!!! 36 * 2 + 9 = 81
    ReplyLengthQ = 9
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
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    print('****')
    return 
#    
def ReadWalkTable(MySocketQ, SlaveAddress):
    CommandCodeQ = 25
    PollLengthQ = 9
    ReplyLengthQ = 73
    MyDataQ=''
    print('**READ WALK TABLE**')
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
    print('Walk Table: '+str(HDataQList))
    print('****')
    return HDataQList
#
def WriteWalkTable(MySocketQ, SlaveAddress, StringOrList):
    CommandCodeQ = 26
    PollLengthQ = 73
    ReplyLengthQ = 9
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
    HStatusQ, HDataQ, HErrorQ = SendPollRecvReply(
        MySocketQ,CommandCodeQ,PollLengthQ,ReplyLengthQ,SlaveAddress,MyDataQ)
    #all return values are strings; HDataQ is already formatted using hexlify
    StatusInterpreter(HErrorQ, HStatusQ)
    HDataQList=[int(HDataQ[ii*4:4+ii*4],16) for ii in range(len(HDataQ)//4)]
    print('****')
    return
#
def HOpen():
    try:
        HIGHLAND_SLAVE_ADDRESS=0 #arbitrary, I think    
        HighlandSocket=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        HighlandSocket.settimeout(1.0)
        HighlandSocket.connect(('highland-mec-01', 2000))##172.21.45.185 allocated on MEC network; 192.168.254.158 originally; pose sometimes as 172.21.46.128 from LeCroy 4GHz oscilloscope
        #controlling computer must be on same subnet as the Highland!
        #Highland's IP address can be changed using the Lantronix DeviceInstaller
        #print('HIGHLAND CONNECTED')
    except:
        print('HIGHLAND NOT CONNECTED')
        return False
    return HighlandSocket
#
def HClose(SocketName):
    SocketName.close()
    #print('HIGHLAND DISCONNECTED')
    return
#
def LXOpen(LStrQ):
    if   str(LStrQ).lower() == 'a':
        host = '172.21.46.100'
    elif str(LStrQ).lower() == 'b':
        host = '172.21.46.120'
    elif str(LStrQ).lower() == '1':
        host = '172.21.46.60'
    elif str(LStrQ).lower() == '2':
        host = '172.21.46.128'
    else:
        print('Invalid scope name! Choose 1, 2, A, or B!!')
        return False
    port = 1861
    try:
        LSock = socket.socket(socket.AF_INET, socket.SOCK_STREAM); LSock.settimeout(1.0);LSock.connect((host, port));
    except:
        print('LECROY '+str(LStrQ)+' NOT CONNECTED!')
        return False
    return LSock
#
def LXClose(SocketName):
    try:
        SocketName.close()
    except:
        print('LECROY SOCKET FAILED!')
    return
#
def send_and_reply(msg,SocketName):
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
    SocketName.sendall(x)
    data = ""
    done = False
    while not done:
        hdr = SocketName.recv(8) ##a bytes object
        hdr = hdr.decode('latin1')##try sg...
        done = (ord(hdr[0]) & 1) == 1  ##had an ord(hdr[0])
        l = struct.unpack(">i", bytes(hdr[4:8],encoding='latin1'))[0]##ADDED bytes(...)
        while (l != 0):
            d = SocketName.recv(l)##########################################
            d = d.decode('latin1')##try sg....
            data = data + d#.decode('utf-8')
            l -= len(d)
    return data
#
def send_only(msg,SocketName):
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
    SocketName.sendall(x)
    return
#
def LFields():
    fields = [
    [0, "DESCRIPTOR_NAME", "string"], 
    [16, "TEMPLATE_NAME", "string"], 
    [32, "COMM_TYPE", "enum", {
        0: "byte", 
        1: "word", 
    }], 
    [34, "COMM_ORDER", "enum", {
        0: "HIFIRST", 
        1: "LOFIRST", 
    }], 
    [36, "WAVE_DESCRIPTOR", "long"], 
    [40, "USER_TEXT", "long"], 
    [44, "RES_DESC1", "long"], 
    [48, "TRIGTIME_ARRAY", "long"], 
    [52, "RIS_TIME_ARRAY", "long"], 
    [56, "RES_ARRAY1", "long"], 
    [60, "WAVE_ARRAY_1", "long"], 
    [64, "WAVE_ARRAY_2", "long"], 
    [68, "RES_ARRAY2", "long"], 
    [72, "RES_ARRAY3", "long"], 
    [76, "INSTRUMENT_NAME", "string"], 
    [92, "INSTRUMENT_NUMBER", "long"], 
    [96, "TRACE_LABEL", "string"], 
    [112, "RESERVED1", "word"], 
    [114, "RESERVED2", "word"], 
    [116, "WAVE_ARRAY_COUNT", "long"], 
    [120, "PNTS_PER_SCREEN", "long"], 
    [124, "FIRST_VALID_PNT", "long"], 
    [128, "LAST_VALID_PNT", "long"], 
    [132, "FIRST_POINT", "long"], 
    [136, "SPARSING_FACTOR", "long"], 
    [140, "SEGMENT_INDEX", "long"], 
    [144, "SUBARRAY_COUNT", "long"], 
    [148, "SWEEPS_PER_ACQ", "long"], 
    [152, "POINTS_PER_PAIR", "word"], 
    [154, "PAIR_OFFSET", "word"], 
    [156, "VERTICAL_GAIN", "float"], 
    [160, "VERTICAL_OFFSET", "float"], 
    [164, "MAX_VALUE", "float"], 
    [168, "MIN_VALUE", "float"], 
    [172, "NOMINAL_BITS", "word"], 
    [174, "NOM_SUBARRAY_COUNT", "word"], 
    [176, "HORIZ_INTERVAL", "float"], 
    [180, "HORIZ_OFFSET", "double"], 
    [188, "PIXEL_OFFSET", "double"], 
    [196, "VERTUNIT", "unit_definition"], 
    [244, "HORUNIT", "unit_definition"], 
    [292, "HORIZ_UNCERTAINTY", "float"], 
    [296, "TRIGGER_TIME", "time_stamp"], 
    [312, "ACQ_DURATION", "float"], 
    [316, "RECORD_TYPE", "enum", {
        0: "single_sweep", 
        1: "interleaved", 
        2: "histogram", 
        3: "graph", 
        4: "filter_coefficient", 
        5: "complex", 
        6: "extrema", 
        7: "sequence_obsolete", 
        8: "centered_RIS", 
        9: "peak_detect", 
    }], 
    [318, "PROCESSING_DONE", "enum", {
        0: "no_processing", 
        1: "fir_filter", 
        2: "interpolated", 
        3: "sparsed", 
        4: "autoscaled", 
        5: "no_result", 
        6: "rolling", 
        7: "cumulative", 
    }], 
    [320, "RESERVED5", "word"], 
    [322, "RIS_SWEEPS", "word"], 
    [324, "TIMEBASE", "enum", {
        0: "1_ps/div", 
        1: "2_ps/div", 
        2: "5_ps/div", 
        3: "10_ps/div", 
        4: "20_ps/div", 
        5: "50_ps/div", 
        6: "100_ps/div", 
        7: "200_ps/div", 
        8: "500_ps/div", 
        9: "1_ns/div", 
        10: "2_ns/div", 
        11: "5_ns/div", 
        12: "10_ns/div", 
        13: "20_ns/div", 
        14: "50_ns/div", 
        15: "100_ns/div", 
        16: "200_ns/div", 
        17: "500_ns/div", 
        18: "1_us/div", 
        19: "2_us/div", 
        20: "5_us/div", 
        21: "10_us/div", 
        22: "20_us/div", 
        23: "50_us/div", 
        24: "100_us/div", 
        25: "200_us/div", 
        26: "500_us/div", 
        27: "1_ms/div", 
        28: "2_ms/div", 
        29: "5_ms/div", 
        30: "10_ms/div", 
        31: "20_ms/div", 
        32: "50_ms/div", 
        33: "100_ms/div", 
        34: "200_ms/div", 
        35: "500_ms/div",
        36: "1_s/div", 
        37: "2_s/div", 
        38: "5_s/div", 
        39: "10_s/div", 
        40: "20_s/div", 
        41: "50_s/div", 
        42: "100_s/div", 
        43: "200_s/div", 
        44: "500_s/div", 
        45: "1_ks/div", 
        46: "2_ks/div", 
        47: "5_ks/div", 
        100: "EXTERNAL", 
    }], 
    [326, "VERT_COUPLING", "enum", {
        0: "DC_50_Ohms", 
        1: "ground", 
        2: "DC_1MOhm", 
        3: "ground", 
        4: "AC_1MOhm", 
    }], 
    [328, "PROBE_ATT", "float"], 
    [332, "FIXED_VERT_GAIN", "enum", {
        0: "1_uV/div", 
        1: "2_uV/div", 
        2: "5_uV/div", 
        3: "10_uV/div", 
        4: "20_uV/div", 
        5: "50_uV/div", 
        6: "100_uV/div", 
        7: "200_uV/div", 
        8: "500_uV/div", 
        9: "1_mV/div", 
        10: "2_mV/div", 
        11: "5_mV/div", 
        12: "10_mV/div", 
        13: "20_mV/div", 
        14: "50_mV/div", 
        15: "100_mV/div", 
        16: "200_mV/div", 
        17: "500_mV/div", 
        18: "1_V/div", 
        19: "2_V/div", 
        20: "5_V/div", 
        21: "10_V/div", 
        22: "20_V/div", 
        23: "50_V/div", 
        24: "100_V/div", 
        25: "200_V/div", 
        26: "500_V/div", 
        27: "1_kV/div", 
    }], 
    [334, "BANDWIDTH_LIMIT", "enum", {
        0: "off", 
        1: "on", 
    }], 
    [336, "VERTICAL_VERNIER", "float"], 
    [340, "ACQ_VERT_OFFSET", "float"], 
    [344, "WAVE_SOURCE", "enum", {
        0: "CHANNEL_1", 
        1: "CHANNEL_2", 
        2: "CHANNEL_3", 
        3: "CHANNEL_4", 
        9: "UNKNOWN", 
    }], ]
    return fields

#
def parsewf(data, verbose=False):
    fields = LFields()
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
        d['RAW'] = np.frombuffer(bytes(data[346:],encoding='latin1'), dtype=np.int8)###buffer should be 5000?? 'utf8'?? print data?
    else:
        d['RAW'] = np.frombuffer(bytes(data[346:],encoding='latin1'), dtype=np.int16) ###buffer should be 5000??
    d['DATA'] = d['VERTICAL_GAIN'] * d['RAW'] - d['VERTICAL_OFFSET']
    return d
#
def readchan(ChannelNo, SocketName, verbose=False):
    while True:
        ready = (int(send_and_reply("INR?", SocketName).split()[1]) & 1) == 1
        if ready:
            data = send_and_reply("C%d:WAVEFORM? ALL" % ChannelNo, SocketName)
            d = parsewf(data, verbose)
            return d

#!!! small block that hasn't been moved yet

#allow input from format of LCWave or ParabolicWave or that kind of thing
def DesiredShapeToList(DesiredOutputPulseShapeQ): #accept list or csv of 140 pts
    if len(DesiredOutputPulseShapeQ) == 140*4:#will accept pre-formatted Hex2Byte text
        PreNormL=[int(DesiredOutputPulseShapeQ[4*ii:4*ii+4],16) for ii in range(len(DesiredOutputPulseShapeQ)//4)]
    elif len(DesiredOutputPulseShapeQ)==140:#will accept a straight list
        PreNormL=DesiredOutputPulseShapeQ
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
                print('Unrecognized format on inrchput file.')
                return
        if len(ListedValues) != 140:
            print('File must have 140 entries; entry count: '+str(len(ListedValues)))
            return
        PreNormL = ListedValues
    DesiredOutputPulseShapePreconv=[entry/float(max(PreNormL)) for entry in PreNormL]
#smoothing/convolution to help ringing at sharp edges go away
    AvgRange=11 #35; set before to 5; must choose odd number 
    FWHM=1.6 #400ps FWHM of pixel supposedly, which is ~16; usually set to 1.6 
    DesiredOutputPulseShape=[]
    WeightList=[math.exp(-4*math.log(2)*((ii+1-round(AvgRange/2.))/FWHM)**2) for ii in range(AvgRange)]
    WSum = sum(WeightList)
    for FETNo in range(140):
        if FETNo<round(AvgRange/2.)-1 or FETNo>140-round(AvgRange/2.):
            DesiredOutputPulseShape+=[DesiredOutputPulseShapePreconv[FETNo]]
        else:
            WSample=sum([DesiredOutputPulseShapePreconv[int(FETNo+(ii+1-round(AvgRange/2.)))]*WeightList[ii] for ii in range(AvgRange)])/WSum
            DesiredOutputPulseShape+=[WSample]
#
    return DesiredOutputPulseShape #list of 140 pts, scaled such that peak is 1
        
def DesiredMaskToList(DesiredOutputPulseShapeQ): #accept list or csv of 140 pts
    if len(DesiredOutputPulseShapeQ) == 140*4:#will accept pre-formatted Hex2Byte text
        PreNormL=[int(DesiredOutputPulseShapeQ[4*ii:4*ii+4],16) for ii in range(len(DesiredOutputPulseShapeQ)//4)]
    elif len(DesiredOutputPulseShapeQ)==140:#will accept a straight list
        PreNormL=DesiredOutputPulseShapeQ
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
        PreNormL = ListedValues
    return [entry/float(max(PreNormL)) for entry in PreNormL] #list of 140 pts, scaled such that peak is 1

def ComboWave(WList): #accept list or csv of 140 pts
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

def FirstInputPulseGuess(DesiredOutputPulseShapeQ, StartingRatio):
    #StartingRatio=.1 #give control over ratio or make it bigger to start with (!!!)
    MaxFETValue=28000
    FirstInputPulseShape=[DesiredOutputPulseShapeQ[ii]*StartingRatio + (((1000*StartingRatio)*(1+math.cos((2*math.pi/14)*1*ii+2)))/MaxFETValue) for ii in range(len(DesiredOutputPulseShapeQ))]#flat + some basic ripple correction
    #can also pre-correct for edge, saturation effects
    #edge may be important to keep pulse from starting in wrap-around region
    #change to 26000 maybe to keep away from roll-over regime?
    return FirstInputPulseShape #list of 140 pts, scaled such that peak is 1

def FirstInputPulseExact(InputPulseQ):
    MaxFETValue=28000
    if len(InputPulseQ) == 140*4:#will accept pre-formatted Hex2Byte text
        PreNormL=[int(InputPulseQ[4*ii:4*ii+4],16) for ii in range(len(InputPulseQ)//4)]
    elif len(InputPulseQ)==140:#will accept a straight list
        PreNormL=InputPulseQ
    elif InputPulseQ.endswith(('.txt','.csv','.dat')):#will accept file
        with open(InputPulseQ,'r') as filehead:
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
        PreNormL = ListedValues
    InputPulseQScaled=[entry/float(max(PreNormL)) for entry in PreNormL]
    return InputPulseQScaled #list of 140 pts, scaled such that peak is 1

def HFormatting(InputPulseShapeQ): #turns a 140-pt list into byte string for Highland WritePulseHeights command
    #scale values by 28000 or whatever "max open" happens to be for FETs
    MaxFETValue=28000
    if np.max(InputPulseShapeQ) <= 1:
        MyDataQ=[int(entry*MaxFETValue) for entry in InputPulseShapeQ]
    else:
        MyDataQ=InputPulseShapeQ[:]
    return MyDataQ #string (w/Hex2Byte)/list/file accepted by WritePulseHeights

def HFormattingPlot(InputPulseShapeQ):
    MaxFETValue=28000.
    ScaledPulse=[(entry/MaxFETValue) for entry in InputPulseShapeQ]
    return ScaledPulse 

def TraceFormatting(PDTrace, PDFETMap, MaxPDValue): #turns PD trace into windowed, averaged list of 140 pts
    #scale values such that maximum energy output would have a value of 1
    #if too hard to know MaxPDValue well, self-scale it down if not converging? (!!!)
    AvgRange=25 #35  ; did 25
    FWHM=4 #400ps FWHM of pixel supposedly, which is ~16; did 4 
    MeasuredOutputPulseShape=[]
    WeightList=[math.exp(-4*math.log(2)*((ii+1-round(AvgRange/2.))/FWHM)**2) for ii in range(AvgRange)]
    WSum = sum(WeightList)
    MX , B = PDFETMap
    for FETNo in range(140):
        Loc = round(MX*FETNo + B)
        WSample=sum([PDTrace[int(Loc+(ii+1-round(AvgRange/2.)))]*WeightList[ii] for ii in range(AvgRange)])/WSum
        MeasuredOutputPulseShape+=[WSample/MaxPDValue]
    return MeasuredOutputPulseShape 
    
def TraceFormatting2(PDTrace, PDFETMap, MaxPDValue): #turns PD trace into windowed, averaged list of 140 pts
    #scale values such that maximum energy output would have a value of 1
    #if too hard to know MaxPDValue well, self-scale it down if not converging? (!!!)
    AvgRange=1#5 #35  ; did 25
    FWHM=1 #400ps FWHM of pixel supposedly, which is ~16; did 4 
    MeasuredOutputPulseShape=[]
    WeightList=[math.exp(-4*math.log(2)*((ii+1-round(AvgRange/2.))/FWHM)**2) for ii in range(AvgRange)]
    WSum = sum(WeightList)
    MX , B = PDFETMap
    for FETNo in range(140):
        Loc = round(MX*FETNo + B)
        WSample=sum([PDTrace[int(Loc+(ii+1-round(AvgRange/2.)))]*WeightList[ii] for ii in range(AvgRange)])/WSum
        MeasuredOutputPulseShape+=[WSample/MaxPDValue]
    return MeasuredOutputPulseShape 

def ErrorSignal(DesiredOutputPulseShape, MeasuredOutputPulseShape):#least-squared difference or something
    #may want to look at percent error in just the region of interest? (!!!)
    #get error better by doing regression or something?
    #this definition is useless when the desired setpoint is 0
    ErrorVals=[abs((DesiredOutputPulseShape[ii]-MeasuredOutputPulseShape[ii])) for ii in range(140)] #/(DesiredOutputPulseShape[ii]+1e-10)
    return sum(ErrorVals)

def UpdatingShapingAlgorithm0(DesiredOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape, DurationListQ, StartStopListQ, StepQ):
    #for arbitary inputs;; duration list clips everything
    G, M, I = DesiredOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape
    PulseMask=np.ceil(PulseGoal(DurationListQ, StartStopListQ))
    NewInputPulseShape=np.clip([abs((StepQ*(G[ii]-M[ii]))+I[ii])*PulseMask[ii] for ii in range(len(G))],0,1)#math.ceil(G) is a mask that disallows values outside the goal
#this helps with problems caused by ripples/noise in the OSC trace
    return NewInputPulseShape #list of 140 pts, maybe written about 10% closer to where we need to go

def UpdatingShapingAlgorithm(DesiredOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape, StepQ):
    #consider edges where there may be undershoot on waveform, so may actually need to turn up to go down? (!!!)
    #StepQ=.1 #step this much of the way; .1 = 10%
    G, M, I = DesiredOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape
    #NewInputPulseShape=[int(MaxFETValue*abs(((1+Q)*G[ii])-(Q*M[ii]))*I[ii]/G[ii]) for ii in range(len(G))] #doesn't work for zero value for G
    #NewInputPulseShape=[int(MaxFETValue*abs((1+Q*(G[ii]-M[ii]))*I[ii])) for ii in range(len(G))]
    NewInputPulseShape=np.clip([abs((StepQ*(G[ii]-M[ii]))+I[ii])*math.ceil(G[ii]) for ii in range(len(G))],0,1)#math.ceil(G) is a mask that disallows values outside the goal
#this helps with problems caused by ripples/noise in the OSC trace
    return NewInputPulseShape #list of 140 pts, maybe written about 10% closer to where we need to go

def UpdatingShapingAlgorithm2(DesiredOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape, StepQ):
    #consider edges where there may be undershoot on waveform, so may actually need to turn up to go down? (!!!)
    #StepQ=.1 #step this much of the way; .1 = 10%
    G, M, I = DesiredOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape
    StepQarr=[StepQ*((float(.75*(ii-51+1))/40)+.25) for ii in range(len(G))]
    #NewInputPulseShape=[int(MaxFETValue*abs(((1+Q)*G[ii])-(Q*M[ii]))*I[ii]/G[ii]) for ii in range(len(G))] #doesn't work for zero value for G
    #NewInputPulseShape=[int(MaxFETValue*abs((1+Q*(G[ii]-M[ii]))*I[ii])) for ii in range(len(G))]
    NewInputPulseShape=np.clip([abs((StepQarr[ii]*(G[ii]-M[ii]))+I[ii])*math.ceil(G[ii]) for ii in range(len(G))],0,1)#math.ceil(G) is a mask that disallows values outside the goal
#this helps with problems caused by ripples/noise in the OSC trace
    return NewInputPulseShape #list of 140 pts, maybe written about 10% closer to where we need to go
    
def UpdatingShapingAlgorithm3(GoalOutputPulseShape, MeasuredOutputPulseShape, InputPulseShape, PercQ):##percent of percent error step in direction
    #consider edges where there may be undershoot on waveform, so may actually need to turn up to go down? (!!!)
    #StepQ=.1 #step this much of the way; .1 = 10%
    G, M, I = DesiredShapeToList(GoalOutputPulseShape), MeasuredOutputPulseShape, InputPulseShape
    #NewInputPulseShape=[int(MaxFETValue*abs(((1+Q)*G[ii])-(Q*M[ii]))*I[ii]/G[ii]) for ii in range(len(G))] #doesn't work for zero value for G
    #NewInputPulseShape=[int(MaxFETValue*abs((1+Q*(G[ii]-M[ii]))*I[ii])) for ii in range(len(G))]
    NewInputPulseShape=np.clip([abs(((1-PercQ)*((G[ii]-M[ii])/(.000000001+G[ii])))*I[ii])*math.ceil(GoalOutputPulseShape[ii]) for ii in range(len(G))],0,1)#math.ceil(G) is a mask that disallows values outside the goal
#this helps with problems caused by ripples/noise in the OSC trace
    return NewInputPulseShape #list of 140 pts, maybe written about 10% closer to where we need to go
#wish list:
#LogEntry() writes the input and corresponding output (and error) to file
#SetScopeParameters() protects viability of PDFETMap result by fixing the window of the scope
#enforce nothing less than 0 or greater than something?
#watch optimization first, maybe on electronic waveform first

#
def rch(OChan,LSock):
    rawdataq = send_and_reply("C{}:WAVEFORM? ALL".format(OChan), LSock)
    return parsewf(rawdataq, False)['DATA']
#
def rchxy(OChan,LSock):
    rawdataq = send_and_reply("C{}:WAVEFORM? ALL".format(OChan), LSock)
    fullaq=parsewf(rawdataq, False);
    yvals=fullaq['DATA'];xvals=[fullaq['HORIZ_OFFSET'] + ii*fullaq['HORIZ_INTERVAL'] for ii in range(len(fullaq['DATA']))];
    return [xvals,yvals]
#    
def rchall(LSock):
    rchans=[]
    for OChan in range(1,5):
        rchans.append(rch(OChan,LSock))
    return rchans
#    
def rchallxy(LSock):
    rchans=[]
    for OChan in range(1,5):
        rchans.append(rchxy(OChan,LSock))
    return rchans
#
def sch(OChan,LSock,FileName):
    rawdataq = send_and_reply("C{}:WAVEFORM? ALL".format(OChan), LSock)
    parseddataq=parsewf(rawdataq, False)
    with open(psfilepath()+'data/'+str(FileName)+'.csv','w',newline='') as f:
        writer=csv.writer(f, delimiter='\n')
        writer.writerow(parseddataq['DATA'])
    with open(psfilepath()+'data/'+str(FileName)+'-h.csv','w',newline='') as f:
        writer=csv.DictWriter(f, parseddataq.keys())
        writer.writeheader()
        writer.writerow(parseddataq)
    return 
#    
def schall(LSock,FileName):
    for OChan in range(1,5):
        sch(OChan,LSock,FileName+'_ch'+str(OChan))
    return
#
def pch(OChan,LSock):
    rawdataq = send_and_reply("C{}:WAVEFORM? ALL".format(OChan), LSock)
    df1=plt.figure()
    plt.plot(parsewf(rawdataq, False)['DATA']);
    df1.show()
    return 
#    
def YFEtrace():
    try:
        SLA=LXOpen('A');time.sleep(0.15);
        pch(1,SLA);time.sleep(0.15);
        LXClose(SLA);time.sleep(0.15);
    except:
        LXClose(SLA);
        print('Failed to display trace')
        return False
#
def pchall(LSock,xlb='none',ylb='none'):
    df1=plt.figure()
    for OChan in range(1,5):
        rawdataq = send_and_reply("C{}:WAVEFORM? ALL".format(OChan), LSock)
        plt.plot(parsewf(rawdataq, False)['DATA']);
    if xlb != 'none':
        plt.xlabel(str(xlb))
    if ylb != 'none':
        plt.ylabel(str(ylb))
    df1.show()
    return 
#    
def epl(listq,xlb='none',ylb='none'):
    df1=plt.figure()
    plt.plot(listq);
    if xlb != 'none':
        plt.xlabel(str(xlb))
    if ylb != 'none':
        plt.ylabel(str(ylb))
    df1.show()
    return
#
def eplxy(listxq,listyq,xlb='none',ylb='none'):
    df1=plt.figure()
    plt.plot(listxq,listyq);
    if xlb != 'none':
        plt.xlabel(str(xlb))
    if ylb != 'none':
        plt.ylabel(str(ylb))
    df1.show()
    return
#
def eplxyloglog(listxq,listyq,xlb='none',ylb='none'):
    df1=plt.figure()
    plt.loglog(listxq,listyq);
    if xlb != 'none':
        plt.xlabel(str(xlb))
    if ylb != 'none':
        plt.ylabel(str(ylb))
    df1.show()
    return
#
def eplsav(listq,FileNameQ,blockdisplay=True,xlb='none',ylb='none'):
    df1=plt.figure()
    plt.plot(listq);
    if xlb != 'none':
        plt.xlabel(str(xlb))
    if ylb != 'none':
        plt.ylabel(str(ylb))
    df1.savefig(str(FileNameQ+'.png'))
    if blockdisplay:
        plt.close(df1)
    return
#    
def eplxysav(listxq,listyq,FileNameQ,abs_path=False,xlb='none',ylb='none'):
    df1=plt.figure()
    plt.plot(listxq,listyq);
    if abs_path:
        figfilename=FileNameQ;
    else:
        figfilename=str(psfilepath()+FileNameQ+'.png')
    if xlb != 'none':
        plt.xlabel(str(xlb))
    if ylb != 'none':
        plt.ylabel(str(ylb))
    df1.savefig(figfilename)        
    plt.close(df1)
    return
#    
def eplcomp(listq,goalq,Map,tMax,xlb='none',ylb='none'):
    formtra=[]
    formtra.append(TraceFormatting(listq,Map,tMax))
    formtra.append(goalq)
    epll(formtra,xlb=xlb,ylb=ylb)
    return
#
def eplcsv(CSVname):
    with open(psfilepath()+'data/'+CSVname+'.csv','r') as filehead:
        RawListQ=filehead.read()
        ListedValues=RawListQ.split('\n')
    epl(ListedValues[:-1])
    return 
#
def epllcsv(CSVHeadname):
    ListofListedValues=[]
    for ii in range(1,5):
        with open(psfilepath()+'data/'+CSVHeadname+'_ch'+str(ii)+'.csv','r') as filehead:
            RawListQ=filehead.read()
            ListedValues=RawListQ.split('\n')
        ListofListedValues.append(ListedValues[:-1])
    epll(ListofListedValues)
    return 
#
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
#
def epll(llist,xlb='none',ylb='none'):
    df1=plt.figure()
    for ii in range(len(llist)):
        plt.plot(llist[ii]);
    if xlb != 'none':
        plt.xlabel(str(xlb))
    if ylb != 'none':
        plt.ylabel(str(ylb))
    df1.show()
#    
def epllxy(llistxyq,xlb='none',ylb='none'):
    df1=plt.figure()
    for ii in range(len(llistxyq)):
        plt.plot(llistxyq[ii][0],llistxyq[ii][1]);
    if xlb != 'none':
        plt.xlabel(str(xlb))
    if ylb != 'none':
        plt.ylabel(str(ylb))
    df1.show()
    return
#
def epllxyloglog(llistxyq,xlb='none',ylb='none'):
    df1=plt.figure()
    for ii in range(len(llistxyq)):
        plt.loglog(llistxyq[ii][0],llistxyq[ii][1]);
    df1.show()
    return
#
def epllt(listq,Map,xlb='none',ylb='none'):
    formtra=[]
    for ii in range(len(listq)):
        formtra.append(TraceFormatting(listq[ii],Map[ii],1))
    epll(formtra,xlb=xlb,ylb=ylb)
    return
#    
def epllcomp(listq,goalq,Map,tMax,xlb='none',ylb='none'):
    formtra=[]
    formtra.append(TraceFormatting(listq[-1],Map,tMax))
    formtra.append(goalq)
    epll(formtra,xlb=xlb,ylb=ylb)
    return
#    
def sumch(chnoQ,sockQ):
    templistQ=rch(chnoQ,sockQ)
    return np.sum(templistQ)

def weichYFE1w(chnoQ,sockQ):
    templistQ=rch(chnoQ,sockQ)
    if np.abs(len(templistQ)-10000)<10:
        bkgrdbuffer=380
    elif np.abs(len(templistQ)-1000)<10:
        bkgrdbuffer=38
    else:
        print('Warning: unanticipated pulse shape array length of '+str(len(templistQ))+', bkgrd subtraction may be off...')
        bkgrdbuffer=1
    bkgrdQ=np.mean(templistQ[:bkgrdbuffer])
    ensampQ=EG1wYFE1in()[0]
    weightQ=ensampQ/np.sum(np.array(templistQ)-bkgrdQ)
    return np.array(weightQ*(np.array(templistQ)-bkgrdQ))
    
def weich1in1w(chnoQ,sockQ):
    templistQ=rch(chnoQ,sockQ)
    if np.abs(len(templistQ)-10000)<10:
        bkgrdbuffer=380
    elif np.abs(len(templistQ)-1000)<10:
        bkgrdbuffer=38
    else:
        print('Warning: unanticipated pulse shape array length of '+str(len(templistQ))+', bkgrd subtraction may be off...')
        bkgrdbuffer=1
    bkgrdQ=np.mean(templistQ[:bkgrdbuffer])
    ensampQ=EG1wYFE1in()[1]
    weightQ=ensampQ/np.sum(np.array(templistQ)-bkgrdQ)
    return np.array(weightQ*(np.array(templistQ)-bkgrdQ))

def weich2in1w(chnoQ,sockQ):
    templistQ=rch(chnoQ,sockQ)
    if np.abs(len(templistQ)-10000)<10:
        bkgrdbuffer=380
    elif np.abs(len(templistQ)-1000)<10:
        bkgrdbuffer=38
    else:
        print('Warning: unanticipated pulse shape array length of '+str(len(templistQ))+', bkgrd subtraction may be off...')
        bkgrdbuffer=1
    bkgrdQ=np.mean(templistQ[:bkgrdbuffer])
    ensampQ=EG1w2in()
    tpva=EpicsSignal(str('MEC:PFN:CH'+str(2*(chnoQ-1)+1)+':ENABLE_RBV'))
    tpvb=EpicsSignal(str('MEC:PFN:CH'+str(2*(chnoQ-1)+2)+':ENABLE_RBV'))
    if tpva.get()*tpvb.get():
        weightQ=ensampQ[0][chnoQ-1]/np.sum(np.array(templistQ)-bkgrdQ)
    else:
        weightQ=1
    return np.array(weightQ*(np.array(templistQ)-bkgrdQ))

def weich2in2w(chnoQ,sockQ): #was weich
    templistQ=rch(chnoQ,sockQ)
    if np.abs(len(templistQ)-10000)<10:
        bkgrdbuffer=380
    elif np.abs(len(templistQ)-1000)<10:
        bkgrdbuffer=38
    else:
        print('Warning: unanticipated pulse shape array length of '+str(len(templistQ))+', bkgrd subtraction may be off...')
        bkgrdbuffer=1
    bkgrdQ=np.mean(templistQ[:bkgrdbuffer])
    ensampQ=EG()
    tpva=EpicsSignal(str('MEC:PFN:CH'+str(2*(chnoQ-1)+1)+':ENABLE_RBV'))
    tpvb=EpicsSignal(str('MEC:PFN:CH'+str(2*(chnoQ-1)+2)+':ENABLE_RBV'))
    if tpva.get()*tpvb.get():
        weightQ=ensampQ[0][0][chnoQ-1]/np.sum(np.array(templistQ)-bkgrdQ)
    else:
        weightQ=1
    return np.array(weightQ*(np.array(templistQ)-bkgrdQ))
    
def weichToPowerVsTime(weiarr):
    return (np.linspace(-5,45,len(weiarr)), np.array(weiarr)/(50e-9/len(weiarr)))

def avgch(chnoQ,sockQ):
    chavg=np.array(rch(chnoQ,sockQ))
    for ii in range(49):
        chavg=chavg+np.array(rch(chnoQ,sockQ))
        time.sleep(0.1)
    chavg=chavg/50.
    return chavg
    
def avgchall(sockQ):
    templistQ=[]
    templistQ.append(avgch(1,sockQ))
    templistQ.append(avgch(2,sockQ))
    templistQ.append(avgch(3,sockQ))
    templistQ.append(avgch(4,sockQ))
    return templistQ

def sumchall(sockQ):
    templistQ=[]
    templistQ.append(rch(1,sockQ))
    templistQ.append(rch(2,sockQ))
    templistQ.append(rch(3,sockQ))
    templistQ.append(rch(4,sockQ))
    energylistQ=[np.sum(templistQ[ii]) for ii in range(len(templistQ))]
    return energylistQ
    
def sumchall2(sockQ):
    coeffQ=[1.855,1.543,1.957,1.850]#[2.178,1.434,1.789,1.228]#[2.440,1.543,1.848,1.621]
    templistQ=[]
    templistQ.append(rch(1,sockQ))
    templistQ.append(rch(2,sockQ))
    templistQ.append(rch(3,sockQ))
    templistQ.append(rch(4,sockQ))
    energylistQ=[np.sum(templistQ[ii])*coeffQ[ii] for ii in range(len(coeffQ))]
    return [energylistQ,np.sum(energylistQ[:2]),np.sum(energylistQ[2:]),np.sum(energylistQ)]

def comchall(sockQ):#combine weighted channels
    coeffQ=[1.855,1.543,1.957,1.850]#[2.178,1.434,1.789,1.228]#[2.440,1.543,1.848,1.621]
    templistQ=[]
    for ii in range(len(coeffQ)):
        templistQ.append(coeffQ[ii]*np.array(rch(ii+1,sockQ)))
    return np.sum(templistQ,0)
    
def weichall(sockQ):
    templistQ=0.0*np.array(rch(1,sockQ))
    ensampQ=EG()
    for ii in range(4):
        subtemplistQ=rch(ii+1,sockQ)
        bkgrdtempQ=np.mean(subtemplistQ[:380])
        tpvii1=EpicsSignal(str('MEC:PFN:CH'+str(2*(ii)+1)+':ENABLE_RBV'))
        tpvii2=EpicsSignal(str('MEC:PFN:CH'+str(2*(ii)+2)+':ENABLE_RBV'))
        if tpvii1.get()*tpvii2.get():
            weighttempQ=ensampQ[0][0][ii]/np.sum(np.array(subtemplistQ)-bkgrdtempQ)
            print(weighttempQ)
        else:
            weighttempQ=1
        templistQ+=np.array(weighttempQ*(np.array(subtemplistQ)-bkgrdtempQ))
    return templistQ

def pchallw(LSock):
    coeffQ=[1.855,1.543,1.957,1.850]#[2.178,1.434,1.789,1.228]#[2.440,1.543,1.848,1.621]
    df1=plt.figure()
    for OChan in range(1,5):
        rawdataq = send_and_reply("C{}:WAVEFORM? ALL".format(OChan), LSock)
        plt.plot(coeffQ[OChan-1]*np.array(parsewf(rawdataq, False)['DATA']));
    df1.show()
    return 

def rgen():
    tpvw=EpicsSignal('MEC:GENTEC:01:CH2:MEAS')
    tpve=EpicsSignal('MEC:GENTEC:01:CH1:MEAS')
    print('Power meter: WEST: ' + str(tpvw.get()) + ', EAST: ' + str(tpve.get()) + ', TOTAL: ' + str(tpvw.get() + tpve.get()))
    return
    
def FixEdges(WavF,DurationListQ,StartStopListQ):
    FirstPix=0#was 50
    DurListQ=np.cumsum([0]+DurationListQ)
    fWavF=WavF[:]
    fWavF[FirstPix+1]=fWavF[FirstPix+2]*1.05
    fWavF[FirstPix]=fWavF[FirstPix+1]*1.1
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
    fWavF[FirstPix+int(4*DurListQ[-1])-2-ContCount]=fWavF[FirstPix+int(4*DurListQ[-1])-3-ContCount]*1.05
    fWavF[FirstPix+int(4*DurListQ[-1])-1-ContCount]=fWavF[FirstPix+int(4*DurListQ[-1])-2-ContCount]*1.05
    for ii in range(len(fWavF)):
        if np.mean(fWavF) < 1:
            if fWavF[ii] > 1:
                fWavF[ii]=1.0
        else:
            if fWavF[ii]>28000:
                fWavF[ii]=28000
    return fWavF

def MBCmodecheck():
    pvMBCmode=EpicsSignal('MEC:LPL:MBC:01:RunningMode_RBV',write_pv='MEC:LPL:MBC:01:RunningMode')#AUTO=0,MAN=1
    return pvMBCmode.get()

def isMBCsafe():#re-write checks as individual functions
    pvMBCpower=EpicsSignal('MEC:64B:PWR:2:Outlet:8:SetControlAction')#WAS 'MEC:S60:PWR:01:Outlet:7:SetControlAction'#read AND write:1=ON,2=OFF
    pvMBCmode=EpicsSignal('MEC:LPL:MBC:01:RunningMode_RBV',write_pv='MEC:LPL:MBC:01:RunningMode')#AUTO=0,MAN=1
    pvMBCsetpt=EpicsSignal('MEC:LPL:MBC:01:AutoCalibration.VAL',write_pv='MEC:LPL:MBC:01:AutoCalibration') #QUAD=0,MIN=1,MAX=2
    pvMBCbias=EpicsSignal('MEC:LPL:MBC:01:BiasValue_RBV',write_pv='MEC:LPL:MBC:01:BiasValue')
    pvMBCfault=EpicsSignal('MEC:LPL:MBC:01:ErrorStatus',write_pv='MEC:LPL:MBC:01:ClearErrors')
    status = True
    print('Checking MBC status...')
    if pvMBCpower.get() != 1:
        print('MBC is not on!!')
        status*=False
    if pvMBCmode.get() != 0:
        print('MBC is not in AUTO mode!!')
        status*=False
    if pvMBCsetpt.get() != 1:
        print('MBC is not in MIN mode!!')
        status*=False
    if not -7000<pvMBCbias.get()<7000:
        print('MBC is out of range!')
        status*=False
    if pvMBCfault.get() != 0:
        print('MBC fault detected!')
        status*=False
    if status:
        biaschk=[]
        print('Checking MBC bias level...',end='',flush=True) 
        for ii in range(3):
            biaschk.append(pvMBCbias.get())
            time.sleep(1);print('..',end='',flush=True);time.sleep(1);print('..',end='',flush=True);
        print('*')
        if np.max(np.abs(np.diff(biaschk))) > 5:
            print('MBC bias level unstable!')
            return False
        else:
            return True
    else:
        return False
#
def resetHighland():
    pvAWG=EpicsSignal('MEC:64B:PWR:2:Outlet:1:SetControlAction')#read AND write:1=ON,2=OFF
    print('Powering off Highland AWG, waiting 3sec...',end='',flush=True);
    pvAWG.put(2);
    dotsleep(3);
    print('Rebooting Highland AWG, waiting 10sec...',end='',flush=True)
    pvAWG.put(1);
    dotsleep(10);
    try:
        S=HOpen();time.sleep(0.15);
        ReadStatus(S,0);time.sleep(0.15);
        HClose(S);time.sleep(0.15);
    except:
        HClose(S);time.sleep(0.15);

def resetMBC():#includes YFEOff() at beginning, but must explicitly ask later to turn YFE back on w/ YFEOn() or YFEsetall(True)
    pvMBCpower=EpicsSignal('MEC:64B:PWR:2:Outlet:8:SetControlAction')#WAS 'MEC:S60:PWR:01:Outlet:7:SetControlAction'#read AND write:1=ON,2=OFF
    pvMBCmode=EpicsSignal('MEC:LPL:MBC:01:RunningMode_RBV',write_pv='MEC:LPL:MBC:01:RunningMode')#AUTO=0,MAN=1
    pvMBCsetpt=EpicsSignal('MEC:LPL:MBC:01:AutoCalibration.VAL',write_pv='MEC:LPL:MBC:01:AutoCalibration') #QUAD=0,MIN=1,MAX=2
    pvMBCbias=EpicsSignal('MEC:LPL:MBC:01:BiasValue_RBV',write_pv='MEC:LPL:MBC:01:BiasValue')
    pvMBCfault=EpicsSignal('MEC:LPL:MBC:01:ErrorStatus',write_pv='MEC:LPL:MBC:01:ClearErrors')
    YFEsetall(False);#pvMBCmode.put(0);time.sleep(1);#make sure it's in AUTO not MAN when it wakes up... saves time
    print('Begin resetting the MBC...')
    #pvMBCpower.put(2);time.sleep(2);
    if pvMBCpower.get() != 1:
        print('Powering on MBC, starting scan...',end='',flush=True)
        pvMBCpower.put(1);time.sleep(1);print('.',end='',flush=True);pvMBCmode.put(0);
        dotsleep(8);
    if pvMBCfault.get() != 0:
        print('Attempting to reset MBC fault...',end='',flush=True)
        pvMBCfault.put(1)
        time.sleep(2);print('*');
    if pvMBCmode.get() != 0:
        print('Setting MBC to AUTO mode, starting scan...',end='',flush=True)
        pvMBCmode.put(0)
        dotsleep(8);
    if pvMBCsetpt.get() != 1:
        print('Setting MBC to MIN mode,starting scan...',end='',flush=True)
        pvMBCsetpt.put(1)
        time.sleep(2);print('*');
    inibias=pvMBCbias.get()
    if not -7000<inibias<7000:
        print('MBC is out of range! Aborting and power-cycling...');
        pvMBCbias.put((np.round(time.time()*1000)%2)*9000*np.sign(inibias));time.sleep(1);
        pvMBCpower.put(2);time.sleep(2);
        resetMBC();return
    biaschk=[]
    print('Checking the initial MBC bias level...',end='',flush=True)
    for ii in range(3):
        biaschk.append(pvMBCbias.get())
        time.sleep(1);print('..',end='',flush=True);time.sleep(1);print('..',end='',flush=True);
    print('*')
    waitloop=True;loopcnt=0;
    while waitloop:
        if np.sum(np.abs(np.diff(biaschk))) > 3:
            print('MBC bias level unstable... '+str(biaschk),end='',flush=True)
            biaschk=[]
            for ii in range(3):
                biaschk.append(pvMBCbias.get())
                time.sleep(1);print('..',end='',flush=True);time.sleep(1);print('..',end='',flush=True);
            print('')
            loopcnt+=1
            if loopcnt >= 15:
                print('MBC bias level stability fail. Aborting and power-cycling...')
                pvMBCbias.put((np.round(time.time()*1000)%2)*9000*np.sign(biaschk[-1]));time.sleep(1);
                pvMBCpower.put(2);time.sleep(2);
                resetMBC();return
        else:
            print('MBC bias level stabilized... '+str(biaschk))
            waitloop = False
    return

def YFEoncheck(display=True):
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

def YFEon():
    if YFEoncheck(display=False):
        print('YFE emission already enabled.');return
    pvPC=EpicsSignal('MEC:S60:PWR:01:Outlet:6:SetControlAction')
    if pvPC.get() != 1:
        pvPC.put(1)
    pvpump=EpicsSignal('MEC:S60:PWR:01:Outlet:7:SetControlAction')
    if pvpump.get() != 1:
        pvpump.put(1)
    pvPS1=EpicsSignal('MEC:S60:PWR:01:Outlet:1:SetControlAction')
    if pvPS1.get() != 1:
        pvPS1.put(1)
    steps=10
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
    if not isMBCsafe():
        YFEsetall(False,displayQ=False)
        print('MBC not configured properly!')
        resetMBC()
        YFEsetall(True,displayQ=False);
    else: #later add check to avoid over-energizing by reading power meter
        for ii in range(len(YFEamp)):
            tempsetcurrpv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[1])
            tempsetemispv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[5])
            tempsetcurrpv.put(0)#set current to 0
            tempsetemispv.put(1)#turn on emission
        print('Initializing eDrives...',end='',flush=True)
        dotsleep(10);
        emissstatlist=[emisspv.get() for emisspv in emisspvlist]
        if all(emissstatlist):
            print('Ramping up currents...')
            YFEsetall(True);
            print('YFE LASER ON')
        else:
            print('Turn on sequence failed. Check emission!')
            YFEoff();
            return False
    if pvPC.get() != 1:
        print('Failed to turn on Pockels cell!')
    if pvpump.get() != 1:
        print('Failed to turn on scroll pump!')
    if pvPS1.get() != 1:
        print('Failed to turn on YFE PS1!')
    return True

def YFEoff():
    pvPC=EpicsSignal('MEC:S60:PWR:01:Outlet:6:SetControlAction')
    pvPC.put(2)
    pvpump=EpicsSignal('MEC:S60:PWR:01:Outlet:7:SetControlAction')
    pvpump.put(2)
    pvPS1=EpicsSignal('MEC:S60:PWR:01:Outlet:1:SetControlAction')
    pvPS1.put(2)
    YFEadd='MEC:LPL:LCO:0'
    YFEamp=['2','3','5','6','1','4']
    YFEsuf=[':SensedCurrent',':ActiveCurrent',':PowerSupply',':Temperature',':Emission_RBV',':Emission',':FaultState.RVAL']
    for ii in range(len(YFEamp)):
        tempsetcurrpv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[1])
        tempsetemispv=EpicsSignal(YFEadd+YFEamp[ii]+YFEsuf[5])
        tempsetcurrpv.put(0)#set current to 0
        tempsetemispv.put(0)#turn on emission
    PFNmodepv=EpicsSignal('MEC:PFN:MODE')
    PFNmodepv.put(0)#makes sure glass chillers turn off when YFEoff... people have been leaving them on
    pvMBC=EpicsSignal('MEC:64B:PWR:2:Outlet:8:SetControlAction')
    pvMBC.put(2)
    #pvAWG=EpicsSignal('MEC:64B:PWR:2:Outlet:1:SetControlAction')
    print('YFE LASER OFF')
    return
    
def YFEget(display=True):
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

def YFEset(mmQ,currQ,display=True):
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

def YFEsetall(IOBool,displayQ=False):
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
    YFEget(display=displayQ);
    return
    
def wshift(wq,pq):
    wshQ=wq[:]
    return np.append(np.delete(wshQ,range(pq)),[0]*pq)

def wswap(wq,b1,b2):#b1 and b2 are segment boards 1 through 10, i.e. b1 and b2 range from 1-10 inclusive
    wsQ=wq[:]
    if 1<=b1<=10 and 1<=b2<=10:
        for ii in range(14):
            wsQ[14*(b1-1)+ii],wsQ[14*(b2-1)+ii]=wsQ[14*(b2-1)+ii],wsQ[14*(b1-1)+ii]
        return wsQ
    else:
        print('Board numbers must be 1-10')
        return
    
def boardswap(b1,b2):
    if 1<=b1<=10 and 1<=b2<=10:
        HSQ=HOpen()
        time.sleep(.15)
        rmcQ=ReadMiscellaneousCalibrations(HSQ,0);
        time.sleep(.15)
        wmcQ=rmcQ[:]
        for ii in range(2):
            wmcQ[2*(b1-1)+ii],wmcQ[2*(b2-1)+ii]=wmcQ[2*(b2-1)+ii],wmcQ[2*(b1-1)+ii]
        WriteMiscellaneousCalibrations(HSQ,0,wmcQ)
        time.sleep(.15)
        HClose(HSQ)
        return
    else:
        print('Board numbers must be 1-10')
        return

def FETsurvey(HSQ,LS1Q):
    #HSQ=HOpen()
    #time.sleep(.15)
    #LS1Q=LOpen()
    #time.sleep(.15)
    qdatalist=[]
    for ii in range(140):
        WritePulseHeights(HSQ,0,IndFETWave([ii+1],28000))#could improve by doing a few spread-out points at a time
        time.sleep(6)
        qpixnodata=readchan(1,LS1Q)['DATA']
        qpixnomax=max(qpixnodata)
        qpixnomaxindex=np.mean([i for i,j in enumerate(qpixnodata) if j == qpixnomax])##could improve by changing to abbreviated centroid around peak, avoiding tail-end bump
        qdatalist.append([qpixnomaxindex,qpixnomax])
    #time.sleep(.15)
    #HClose(HSQ)
    #time.sleep(.15)
    #LClose(LS1Q)
    #time.sleep(.15)
    return qdatalist

def FastFETSurvey(HSQ,LS1Q):
    #HSQ=HOpen()
    #time.sleep(.15)
    #LS1Q=LOpen()
    #time.sleep(.15)
    qdatalist=[]
    for ii in range(10):
        WritePulseHeights(HSQ,0,IndFETWave([jj*10+ii+1 for jj in range(14)],28000))#could improve by doing a few spread-out points at a time
        time.sleep(6)
        qpixnodata=readchan(1,LS1Q)['DATA']
        qpeakcoords=signal.find_peaks_cwt(np.clip(qpixnodata,max(qpixnodata)/5,max(qpixnodata)),np.arange(180,200))#threshold->clip
        #this is 475ps/(2.5ps/pix)=190pix expected for S1 scope at max sampling; S2 scope needs different
        if len(qpeakcoords) != 15:
            print('\nWrong number of peaks detected!\n',len(qpeakcoords))
            return
        qdatalist.append([qpeakcoords[1:],qpixnodata[qpeakcoords[1:]]])
    #time.sleep(.15)
    #HClose(HSQ)
    #time.sleep(.15)
    #LClose(LS1Q)
    #time.sleep(.15)
    qdatalist2=np.array(qdatalist).transpose(2,0,1).reshape(140,2)
    qTimeErrInHpix=np.array([qdatalist2[ii,0]-np.mean(qdatalist2[:14,0])+100*(6.5-ii) for ii in range(140)])*2.5/1.8
    qTimeErrInHpixBoardAvg=np.array([np.mean(qTimeErrInHpix[14*ii:14*ii+14]) for ii in range(10)])
    epl(qTimeErrInHpix)
    epl(qTimeErrInHpixBoardAvg)
    return np.array([qTimeErrInHpix,qTimeErrInHpixBoardAvg])
    #qteHlist=[]
    #qteHlist.append(FastFETSurvey(S,S1)) #a few times
    
    #qteHlistavg=[[mean(np.array(list(np.array(qteHlist)[:,0]))[:,ii]) for ii in range(140)],[mean(np.array(list(np.array(qteHlist)[:,1]))[:,jj]) for jj in range(10)]]
    #qrmcnew=MiscCalBoardCorrection(S,qrmcnew,qteHlistavg[1])
    #qteHlistavg3=[[mean(np.array(list(np.array(qteHlist3)[:,0]))[:,ii]) for ii in range(140)],[mean(np.array(list(np.array(qteHlist3)[:,1]))[:,jj]) for jj in range(10)]]

def VeryFastFETSurvey(HSQ,LS1Q):#returns np.array([qTimeErrInHpix,qTimeErrInHpixBoardAvg])
    #HSQ=HOpen()
    #time.sleep(.15)
    #LS1Q=LOpen()
    #time.sleep(.15)
    qdatalist=[]
    qrfis=ReadFiducialImpulseSettings(HSQ,0)
    WriteFiducialImpulseSettings(HSQ,0,4000,45000)#try full-time 15000,45000
    time.sleep(.15)
    for ii in range(4):
        WritePulseHeights(HSQ,0,IndFETWave([jj*4+ii+1 for jj in range(35)],15000))#could improve by doing a few spread-out points at a time
        time.sleep(6)
        qpixnodata=readchan(1,LS1Q)['DATA']
        #qpeakcoords=signal.find_peaks_cwt(np.clip(qpixnodata,max(qpixnodata)/4,max(qpixnodata)),np.arange(90,110))#threshold->clip
        qpeakcoordspre=signal.find_peaks_cwt(np.clip(qpixnodata,max(qpixnodata)/4,max(qpixnodata))-max(qpixnodata)/4,np.arange(90,110));
        qpeakcoords=[pt for pt in qpeakcoordspre if qpixnodata[pt]>max(qpixnodata)/4]
        #this is 475ps/(2.5ps/pix)=190pix expected for S1 scope at max sampling; S2 scope needs different
        if len(qpeakcoords) != 36:
            print('\nWrong number of peaks detected!\n',len(qpeakcoords))
            #time.sleep(.15)
            #HClose(HSQ)
            #time.sleep(.15)
            #LClose(LS1Q)
            #time.sleep(.15)
            return
        #qdatalist.append([qpeakcoords[1:],qpixnodata[qpeakcoords[1:]]])
        #try post-pulse trigger instead
        qdatalist.append([qpeakcoords[:-1],qpixnodata[qpeakcoords[:-1]]])
    #time.sleep(.15)
    #HClose(HSQ)
    #time.sleep(.15)
    #LClose(LS1Q)
    #time.sleep(.15)
    qdatalist2=np.array(qdatalist).transpose(2,0,1).reshape(140,2)
    qTimeErrInHpix=np.array([qdatalist2[ii,0]-np.mean(qdatalist2[:14,0])+100*(6.5-ii) for ii in range(140)])*2.5/1.8
    qTimeErrInHpixBoardAvg=np.array([np.mean(qTimeErrInHpix[14*ii:14*ii+14]) for ii in range(10)])
    epl(qTimeErrInHpix)
    epl(qTimeErrInHpixBoardAvg)
    return np.array([qTimeErrInHpix,qTimeErrInHpixBoardAvg])

def MiscCalCorrection(HSQ,MiscCalOldQ,TimeErrInHpixBoardAvgQ):
    MiscCalNewQ=MiscCalOldQ[:]
    for ii in range(10):
        MiscCalNewQ[1+2*ii]=MiscCalNewQ[1+2*ii]-int(round(TimeErrInHpixBoardAvgQ[ii]))
    WriteMiscellaneousCalibrations(HSQ,0,MiscCalNewQ)
    time.sleep(.15)
    return MiscCalNewQ
    
def WaveTimeCalCorrection(HSQ,WaveTimeCalOldQ,TimeErrInHpixQ):
    WaveTimeCalNewQ=WaveTimeCalOldQ[:]
    for ii in range(140):
        WaveTimeCalNewQ[ii]=WaveTimeCalNewQ[ii]-int(round(TimeErrInHpixQ[ii]))
    WriteWaveTimeCalibrations(HSQ,0,WaveTimeCalNewQ)
    time.sleep(.15)
    return WaveTimeCalNewQ

def ScanAndShift(HSQ,LS1Q):
    #ideally: auto handling of sockets
    #need error handling/better stability of code first
    #HSQ=HOpen()
    #time.sleep(.15)
    #LS1Q=LOpen()
    #time.sleep(.15)
    #way to check scope settings? currently did this with the following:
    #5 mV/div, 10sweeps, 5ns/div, -33ns delay, 13ns deskew, -13.8mV offset
    ##### use on 20210209: on 5n/div on LeCroy1, trigger off another channel, set fiducial at 9th div
    #ideally: read YFE settings, turn down, turn back up after done
    #set fiducial and everything like that
    scanresults=[]
    #we don't care about historical values; we just want things fixed
    #so don't read in/pass in parameters; just get them straight from HSQ
    PulseHeightQ=ReadPulseHeights(HSQ,0)
    MiscCalQ=ReadMiscellaneousCalibrations(HSQ,0)
    time.sleep(.15)
    WaveTimeCalQ=ReadWaveTimeCalibrations(HSQ,0)
    time.sleep(.15)
    scanresults.append(VeryFastFETSurvey(HSQ,LS1Q))
    #test if need correction
    if any(abs(elem)>2.5 for elem in scanresults[-1][1]):
        print('Adjusting MiscCal\n')
        MiscCalQ=MiscCalCorrection(HSQ,MiscCalQ,scanresults[-1][1])
        time.sleep(.15)
        scanresults.append(VeryFastFETSurvey(HSQ,LS1Q))
    if any(abs(elem)>5.5 for elem in scanresults[-1][0]):
        #this is factor of 2 away from "bad" ("20ps"=11.1 Hpix of error)
        print('Adjusting WaveTimeCal\n')
        WaveTimeCalQ=WaveTimeCalCorrection(HSQ,WaveTimeCalQ,scanresults[-1][0])
        time.sleep(.15)
        scanresults.append(VeryFastFETSurvey(HSQ,LS1Q))
    if any(abs(elem)>2.5 for elem in scanresults[-1][1]) or any(abs(elem)>5.5 for elem in scanresults[-1][0]):
        print('Consider running a second iteration')
        #ideally: re-cal with for loops and iterate until corrected
    #time.sleep(.15)
    #HClose(HSQ)
    #time.sleep(.15)
    #LClose(LS1Q)
    #time.sleep(.15)
    WritePulseHeights(HSQ,0,PulseHeightQ)
    return


def FETsurveyfull():
    HSQ=HOpen()
    time.sleep(.15)
    LS1Q=LXOpen('1')
    time.sleep(.15)
    qdatalist=[]
    for ii in range(140):
        qptdatalist=[]
        for jj in range(5):
            WritePulseHeights(HSQ,0,IndFETWave([ii+1],int((jj+1)*65535/5)))#could improve by doing a few spread-out points at a time
            time.sleep(24)
            qpixnodata=readchan(1,LS1Q)['DATA'][2400:]
            qpixnomax=max(qpixnodata)
            qpixnomaxindex=np.mean([i for i,j in enumerate(qpixnodata) if j == qpixnomax])##could improve by changing to abbreviated centroid around peak, avoiding tail-end bump
            qptdatalist.append([2400+qpixnomaxindex,qpixnomax])
        qdatalist.append(qptdatalist)
    pickledump2(qdatalist,psfilepath()+'fullFETsurvey20181106.p')
    time.sleep(.15)
    HClose(HSQ)
    time.sleep(.15)
    LXClose(LS1Q)
    time.sleep(.15)
    return qdatalist
    
def HParamReset():##need to fix pickle... !!!
    HSQ=HOpen()
    time.sleep(.15)
    [qrphQ,qrwacQ,qrmcQ,qrwtcQ,qrwtQ]=pickle.load(open(psfilepath()+'HighlandParameterSnapshot20181116.p','rb'))#1108 original
    WriteFiducialImpulseSettings(HSQ,0,0,0) #turn off fiducial; was at (HSQ,0,65535,0)
    time.sleep(.15)
    WritePulseHeights(HSQ,0,qrphQ)
    #time.sleep(.15)
    #WriteWaveAmplitudeCalibrations(HSQ,0,qrwacQ)
    #time.sleep(.15)
    #WriteMiscellaneousCalibrations(HSQ,0,qrmcQ)
    #time.sleep(.15)
    #WriteWaveTimeCalibrations(HSQ,0,qrwtcQ)
    #time.sleep(.15)
    #WriteWalkTable(HSQ,0,qrwtQ)
    #time.sleep(.15)
    time.sleep(.15)
    HClose(HSQ)
    time.sleep(.15) 
    return

def EG():
    pveab=EpicsSignal('MEC:LAS:GENTEC:03:CH1:MEAS')
    eab=pveab.get()
    pveef=EpicsSignal('MEC:LAS:GENTEC:03:CH2:MEAS')
    eef=pveef.get()
    pvegh=EpicsSignal('MEC:LAS:GENTEC:04:CH1:MEAS')
    egh=pvegh.get()
    pveij=EpicsSignal('MEC:LAS:GENTEC:04:CH2:MEAS')
    eij=pveij.get()
    EAB,EEF,EGH,EIJ=round(1.05*1.17*eab/.00760/1.006/1.0412/1.0799/1.0478,4),round(1.13*.860*eef/.00686/1.006/.9634/0.8410/0.9517,4),round(.914*.897*egh/.00655/1.015/.9692/0.883/0.9650,4),round(1.07*1.25*eij/.00608/1.015/1.1232/1.075/1.0863,4)
    guessarray=[[EAB,EEF,EGH,EIJ],round(EAB+EEF,4),round(EGH+EIJ,4),round(EAB+EEF+EGH+EIJ,4)]
    #print(guessarray)
    try:
        eabefpv=EpicsSignal('MEC:GENTEC:01:CH2:MEAS')
        EABEF=eabefpv.get()
        eghijpv=EpicsSignal('MEC:GENTEC:01:CH1:MEAS')
        EGHIJ=eghijpv.get()
    except:
        EABEF=-1;EGHIJ=-1
    realarray=[EABEF,EGHIJ,EABEF+EGHIJ]
    #print(realarray)
    return [guessarray,realarray]

def EG1wYFE1in():
    pveyfe=EpicsSignal('MEC:LAS:LEM:03:A:CUR_DISP')
    eyfe=pveyfe.get()
    pve1in=EpicsSignal('MEC:LAS:LEM:03:B:CUR_DISP')
    e1in=pve1in.get()
    EYFE,E1IN=eyfe*0.3285-0.00039,e1in*0.5971#was 0.3578
    guessarray=[round(EYFE,4),round(E1IN,4)]
    return guessarray

def EG1w2in():
    pveab=EpicsSignal('MEC:LAS:GENTEC:02:CH1:MEAS')
    eab=pveab.get()
    pveef=EpicsSignal('MEC:LAS:GENTEC:02:CH2:MEAS')
    eef=pveef.get()
    pvegh=EpicsSignal('MEC:LAS:GENTEC:01:CH1:MEAS')
    egh=pvegh.get()
    pveij=EpicsSignal('MEC:LAS:GENTEC:01:CH2:MEAS')
    eij=pveij.get()
    EAB,EEF,EGH,EIJ=round(eab*224.0,4),round(eef*177.5,4),round(egh*307.4*0.849,4),round(eij*113.2,4)
    guessarray=[[EAB,EEF,EGH,EIJ],round(EAB+EEF,4),round(EGH+EIJ,4),round(EAB+EEF+EGH+EIJ,4)]
    return guessarray

def EGall(return_txt=False,chamber_meter_in=False):#add LasCoeff(); add total EAST/WEST/COMBINED energies
    [en1wYFE, en1w1in] =EG1wYFE1in()
    [en1wAB, en1wEF, en1wGH, en1wIJ] = EG1w2in()[0]
    [en2wAB, en2wEF, en2wGH, en2wIJ] = EG()[0][0]
    [enWEST, enEAST]=EG()[1][:2]
    tspv=EpicsSignal('MEC:PFN:SINCESHOT')
    [cAB,cEF,cGH,cIJ]=LasCoeff()
    ###
    headenb = LasCoeff()
    wpreadpvlist=['MEC:NS1:MMS:02.RBV','MEC:NS1:MMS:01.RBV','MEC:LAS:MMN:30.RBV','MEC:LAS:MMN:29.RBV']
    wpwritepvlist=['MEC:NS1:MMS:02.VAL','MEC:NS1:MMS:01.VAL','MEC:LAS:MMN:30.VAL','MEC:LAS:MMN:29.VAL']
    wppvlist=[EpicsSignal(wpreadpvlist[ii],write_pv=wpwritepvlist[ii]) for ii in range(4)]
    headstr='';wpstr='';headlist=['AB','EF','GH','IJ'];
    for ii in range(4):
        wpstr=wpstr+headlist[ii]+': '+str(round(wppvlist[ii].get(),3))+', '
    wpenlist=[1e-12 + headenb[ii]*np.cos((np.pi/180)*2*wppvlist[ii].get())**2 for ii in range(4)]
    ###
    strlist=[]
    strlist.append('Date: '+ datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
    print(strlist[-1]);strlist.append('\n');
    strlist.append('Time since last shot: '+tspv.get())
    print(strlist[-1]);strlist.append('\n');
    strlist.append('Current pulse target is: '+str(pickle.load(open(psfilepath()+'Psns.p','rb')))+' ns, '+str(pickle.load(open(psfilepath()+'SSs.p','rb')))+' % of max power.')
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




def EG2():
    print('use EG() instead')
    c1=.008508;c2=.0007500;c3=-.01210;
    d1=-.000676;d2=.005833;d3=.01347;
    e1=.005615;e2=-.0005531;e3=.01128;
    f1=.00004615;f2=.006142;f3=.001031;
    pveab=EpicsSignal('MEC:LAS:GENTEC:03:CH1:MEAS')
    sab=pveab.get()
    pveef=EpicsSignal('MEC:LAS:GENTEC:03:CH2:MEAS')
    sef=pveef.get()
    pvegh=EpicsSignal('MEC:LAS:GENTEC:04:CH1:MEAS')
    sgh=pvegh.get()
    pveij=EpicsSignal('MEC:LAS:GENTEC:04:CH2:MEAS')
    sij=pveij.get()
    EAB=(((sab-c3)*d2)-((sef-d3)*c2))/((c1*d2)-(c2*d1));
    EEF=(((sab-c3)*d1)-((sef-d3)*c1))/((c2*d1)-(c1*d2));
    EGH=(((sgh-e3)*f2)-((sij-f3)*e2))/((e1*f2)-(e2*f1));
    EIJ=(((sgh-e3)*f1)-((sij-f3)*e1))/((e2*f1)-(e1*f2));
    guessarray=[[EAB,EEF,EGH,EIJ],EAB+EEF,EGH+EIJ,EAB+EEF+EGH+EIJ]
    #print(guessarray)
    pvEABEF=EpicsSignal('MEC:GENTEC:01:CH2:MEAS')
    EABEF=pvEABEF.get()
    pvEGHIJ=EpicsSignal('MEC:GENTEC:01:CH1:MEAS')
    EGHIJ=pvEGHIJ.get()
    realarray=[EABEF,EGHIJ,EABEF+EGHIJ]
    #print(realarray)
    return [guessarray,realarray]

def EG3():
    print('use EG() instead')
    c1=.008508;c2=.0007500;c3=-.01210;
    d1=-.000676;d2=.005833;d3=.01347;
    e1=.005615;e2=-.0005531;e3=.01128;
    f1=.00004615;f2=.006142;f3=.001031;
    pveab=EpicsSignal('MEC:LAS:GENTEC:03:CH1:MEAS')
    sab=pveab.get()
    pveef=EpicsSignal('MEC:LAS:GENTEC:03:CH2:MEAS')
    sef=pveef.get()
    pvegh=EpicsSignal('MEC:LAS:GENTEC:04:CH1:MEAS')
    sgh=pvegh.get()
    pveij=EpicsSignal('MEC:LAS:GENTEC:04:CH2:MEAS')
    sij=pveij.get()
    EAB=(((sab-c3)*d2)-((sef-d3)*c2))/((c1*d2)-(c2*d1));
    EEF=(((sab-c3)*d1)-((sef-d3)*c1))/((c2*d1)-(c1*d2));
    EGH=(((sgh-e3)*f2)-((sij-f3)*e2))/((e1*f2)-(e2*f1));
    EIJ=(((sgh-e3)*f1)-((sij-f3)*e1))/((e2*f1)-(e1*f2));
    guessarray=[[EAB,EEF,EGH,EIJ],EAB+EEF,EGH+EIJ,EAB+EEF+EGH+EIJ]
    print('energy guesses:')
    print('AB: '+str(round(EAB,2))+' J, '+'EF: '+str(round(EEF,2))+' J, '+'GH: '+str(round(EGH,2))+' J, '+'IJ: '+str(round(EIJ,2))+' J')
    print('ABEF: '+str(round(EAB+EEF,2))+' J, '+'GHIJ: '+str(round(EGH+EIJ,2))+' J, '+'total: '+str(round(EAB+EEF+EGH+EIJ,2))+' J')
    #EABEF=pypsepics.get('MEC:GENTEC:01:CH2:MEAS')
    #EGHIJ=pypsepics.get('MEC:GENTEC:01:CH1:MEAS')
    #realarray=[EABEF,EGHIJ,EABEF+EGHIJ]
    #print(realarray)
    return

def EG3manual(sabQ,sefQ,sghQ,sijQ):
    c1=.008508;c2=.0007500;c3=-.01210;
    d1=-.000676;d2=.005833;d3=.01347;
    e1=.005615;e2=-.0005531;e3=.01128;
    f1=.00004615;f2=.006142;f3=.001031;
    sab=sabQ
    sef=sefQ
    sgh=sghQ
    sij=sijQ
    EAB=(((sab-c3)*d2)-((sef-d3)*c2))/((c1*d2)-(c2*d1));
    EEF=(((sab-c3)*d1)-((sef-d3)*c1))/((c2*d1)-(c1*d2));
    EGH=(((sgh-e3)*f2)-((sij-f3)*e2))/((e1*f2)-(e2*f1));
    EIJ=(((sgh-e3)*f1)-((sij-f3)*e1))/((e2*f1)-(e1*f2));
    guessarray=[[EAB,EEF,EGH,EIJ],EAB+EEF,EGH+EIJ,EAB+EEF+EGH+EIJ]
    print('energy guesses:')
    print('AB: '+str(round(EAB,2))+' J, '+'EF: '+str(round(EEF,2))+' J, '+'GH: '+str(round(EGH,2))+' J, '+'IJ: '+str(round(EIJ,2))+' J')
    print('ABEF: '+str(round(EAB+EEF,2))+' J, '+'GHIJ: '+str(round(EGH+EIJ,2))+' J, '+'total: '+str(round(EAB+EEF+EGH+EIJ,2))+' J')
    #EABEF=pypsepics.get('MEC:GENTEC:01:CH2:MEAS')
    #EGHIJ=pypsepics.get('MEC:GENTEC:01:CH1:MEAS')
    #realarray=[EABEF,EGHIJ,EABEF+EGHIJ]
    #print(realarray)
    return
    
def LasCoeff():
    cAB1pv=EpicsSignal('MEC:PFN:CH1:ENABLE_RBV')
    cAB2pv=EpicsSignal('MEC:PFN:CH2:ENABLE_RBV')
    cAB=cAB1pv.get()*cAB2pv.get()
    cEF1pv=EpicsSignal('MEC:PFN:CH3:ENABLE_RBV')
    cEF2pv=EpicsSignal('MEC:PFN:CH4:ENABLE_RBV')
    cEF=cEF1pv.get()*cEF2pv.get()
    cGH1pv=EpicsSignal('MEC:PFN:CH5:ENABLE_RBV')
    cGH2pv=EpicsSignal('MEC:PFN:CH6:ENABLE_RBV')
    cGH=cGH1pv.get()*cGH2pv.get()
    cIJ1pv=EpicsSignal('MEC:PFN:CH7:ENABLE_RBV')
    cIJ2pv=EpicsSignal('MEC:PFN:CH8:ENABLE_RBV')
    cIJ=cIJ1pv.get()*cIJ2pv.get()
    return np.array([cAB,cEF,cGH,cIJ])
    
def tempsave(FileNameQ):
    dumperQ=[weich(1,S2),weich(2,S2),weich(3,S2),weich(4,S2),LasCoeff()[0]*weich(1,S2)+LasCoeff()[1]*weich(2,S2)+LasCoeff()[2]*weich(3,S2)+LasCoeff()[3]*weich(4,S2),LasCoeff()*EG()[0][0]]
    pickledump2(dumperQ,str(FileNameQ)+'.p')
    print(FileNameQ)
    print(dumperQ[-1])
    epll(dumperQ[:4])
    return dumperQ


def abijswitch():
    PFNmodepv=EpicsSignal('MEC:PFN:MODE')
    PFNmodepv.put(0)
    PFNch1enpv=EpicsSignal('MEC:PFN:CH1:ENABLE')
    PFNch2enpv=EpicsSignal('MEC:PFN:CH2:ENABLE')
    PFNch7enpv=EpicsSignal('MEC:PFN:CH7:ENABLE')
    PFNch8enpv=EpicsSignal('MEC:PFN:CH8:ENABLE')
    NS1MMS02pv=EpicsSignal('MEC:NS1:MMS:02.VAL')
    LASMMN29pv=EpicsSignal('MEC:LAS:MMN:29.VAL')
    time.sleep(3)
    if pypsepics.get('MEC:PFN:CH1:ENABLE_RBV'):
        PFNch1enpv.put(0)
        PFNch2enpv.put(0)
        PFNch7enpv.put(1)
        PFNch8enpv.put(1)
        NS1MMS02pv.put(51)
        LASMMN29pv.put(27.5)
    else:
        PFNch1enpv.put(1)
        PFNch2enpv.put(1)
        PFNch7enpv.put(0)
        PFNch8enpv.put(0)
        NS1MMS02pv.put(96)
        LASMMN29pv.put(-17.5)
    time.sleep(3)
    PFNmodepv.put(1)
    time.sleep(3)
    PFNmodepv.put(2)
    return

def runrecall(RunNameQ):
    epl(np.genfromtxt('/reg/neh/operator/mecopr/mecpython/experiments/LT8917/scope/'+RunNameQ+'chsum.txt'))
    print('Pulse energies:')
    print(np.genfromtxt('/reg/neh/operator/mecopr/mecpython/experiments/LT8917/scope/'+RunNameQ+'energies.txt'))
    return
    
def findfid(TraceInQ): 
    #make it fast by narrowing area where we know peak should be
    TQ=TraceInQ[:]; maxTQ=max(TQ); minTQ=min(TQ);
    TQP=signal.find_peaks_cwt(np.clip(TQ,(maxTQ-.8*(maxTQ-minTQ)),maxTQ),np.arange(5,15));
    if (TQP[-1]-TQP[-2]) > 1000:
        return TQP[-1]
    else:
        print('check your answer...')
        return TQP[-1]

#returns the last pulse shape written
##originally for shapes based off old traces
#add condition to loop until PFN has been fired?
def YFEOptimizer(InputPulseShapeQ, DesiredOutputPulseShapeQ, ChannelNo, LSocket, HSocket, MaxIterations, AcceptableErrorCriterion, PDFETMap, StepQ, DurationListQ, StartStopListQ, verbose=False):
    IterationCounter=0
    InputPulseShape=InputPulseShapeQ[:]#140 pts scaled to 28000.
    DesiredOutputPulseShape=DesiredOutputPulseShapeQ[:]
    #WritePulseHeights(HSocket,0,InputPulseShape)
    #####
    #####epl(InputPulseShape)
    totplotQ=[]
    while True:
        ANewShotIsReady = (int(send_and_reply("INR?", LSocket).split()[1]) & 1) == 1
        if ANewShotIsReady:
            RawData = send_and_reply("C{}:WAVEFORM? ALL".format(ChannelNo), LSocket)
            PDTrace = parsewf(RawData, verbose)['DATA']
            FormattedMeasuredOutputPulseShape = TraceFormatting(PDTrace, PDFETMap, 1)#doing no scaling since no YFE energy meter...
            #hope for the best!!
            #####
            #epll([FormattedMeasuredOutputPulseShape,DesiredOutputPulseShape])
            totplotQ.append(FormattedMeasuredOutputPulseShape)
            if IterationCounter >= MaxIterations:
                print('Maximum number of iterations reached.')
                epll(totplotQ)
                return totplotQ #InputPulseShape
            PulseError = ErrorSignal(DesiredOutputPulseShape, FormattedMeasuredOutputPulseShape)
            if PulseError < AcceptableErrorCriterion:
                print('Error tolerance achieved.')
                return InputPulseShape
            InputPulseShape = 28000.0*np.array(FixEdges(UpdatingShapingAlgorithm0(DesiredOutputPulseShape, FormattedMeasuredOutputPulseShape, np.array(InputPulseShape)/28000.,DurationListQ,StartStopListQ,StepQ),DurationListQ,StartStopListQ))
            #####
            #####epl(InputPulseShape)
            #WritePulseHeights(HSocket,0,InputPulseShape)
            IterationCounter+=1

##test for shapes not based on previous traces
def YFEOptimizer2(InputPulseShapeQ, DesiredOutputPulseShapeQ, ChannelNo, LSocket, HSocket, MaxIterations, AcceptableErrorCriterion, PDFETMap, StepQ, DurationListQ, StartStopListQ, verbose=False):
    IterationCounter=0
    InputPulseShape=InputPulseShapeQ[:]#140 pts scaled to 28000.
    DesiredOutputPulseShape=DesiredOutputPulseShapeQ[:]
    #WritePulseHeights(HSocket,0,InputPulseShape)
    #####
    #####epl(InputPulseShape)
    while True:
        ANewShotIsReady = (int(send_and_reply("INR?", LSocket).split()[1]) & 1) == 1
        if ANewShotIsReady:
            RawData = send_and_reply("C{}:WAVEFORM? ALL".format(ChannelNo), LSocket)
            protoPDTrace = parsewf(RawData, verbose)['DATA']
            bkgrdQ=np.mean(protoPDTrace[:380])
            weightQ=1
            PDTrace=np.array(weightQ*(np.array(protoPDTrace)-bkgrdQ))
            #figure out how to check to make sure the fiducial is actually there; if not, then break
            #if PDTrace[26965]
            FormattedMeasuredOutputPulseShape = TraceFormatting(PDTrace, PDFETMap, 1)#doing no scaling since no YFE energy meter...
            #hope for the best!!
            #####
            #epll([FormattedMeasuredOutputPulseShape,DesiredOutputPulseShape])
            if IterationCounter >= MaxIterations:
                print('Maximum number of iterations reached.')
                return InputPulseShape
            PulseError = ErrorSignal(DesiredOutputPulseShape, FormattedMeasuredOutputPulseShape)
            if PulseError < AcceptableErrorCriterion:
                print('Error tolerance achieved.')
                return InputPulseShape
            InputPulseShape = 28000.0*np.array(FixEdges(UpdatingShapingAlgorithm0(DesiredOutputPulseShape, FormattedMeasuredOutputPulseShape, np.array(InputPulseShape)/28000.,DurationListQ,StartStopListQ,StepQ),DurationListQ,StartStopListQ))
            #####
            #####epl(InputPulseShape)
            WritePulseHeights(HSocket,0,InputPulseShape)
            IterationCounter+=1

def fidon():
    try:
        S=HOpen(); time.sleep(.15);
        WriteFiducialImpulseSettings(S,0,20000,45000); time.sleep(.15);
        HClose(S); time.sleep(.15);
    except:
        HClose(S);
        print('Error!')
        return False
    return

def fidoff():
    try:
        S=HOpen(); time.sleep(.15);
        WriteFiducialImpulseSettings(S,0,0,0); time.sleep(.15);
        HClose(S); time.sleep(.15);
    except:
        HClose(S);
        print('Error!')
        return False
    return

def HWPon(ArmStrQ,set_T=1):#fullON=1;fullOFF=0
    armlist=['AB','EF','GH','IJ'];
    if (set_T<0) or (set_T>1):
        set_T=1;print('Error: set_T must be between 0 and 1! Using set_T=1 instead...')
    if ArmStrQ.lower() == 'all':
        ArmStrQ = 'ABEFGHIJ'
    HWPABpv=EpicsSignal('MEC:NS1:MMS:02.VAL');HWPABpv2=EpicsSignal('MEC:NS1:MMS:02.RBV');
    HWPEFpv=EpicsSignal('MEC:NS1:MMS:01.VAL');HWPEFpv2=EpicsSignal('MEC:NS1:MMS:01.RBV');
    HWPGHpv=EpicsSignal('MEC:LAS:MMN:30.VAL');HWPGHpv2=EpicsSignal('MEC:LAS:MMN:30.RBV');
    HWPIJpv=EpicsSignal('MEC:LAS:MMN:29.VAL');HWPIJpv2=EpicsSignal('MEC:LAS:MMN:29.RBV');
    HWPpvlist=[HWPABpv,HWPEFpv,HWPGHpv,HWPIJpv];HWPpv2list=[HWPABpv2,HWPEFpv2,HWPGHpv2,HWPIJpv2];
    motor_angle=np.arccos(np.sqrt(set_T))*(180/np.pi)*0.5 
    print('Moving waveplates '+ArmStrQ+' to '+str(round(motor_angle,4))+'deg...');
    ABchk=('A' in ArmStrQ) or ('a' in ArmStrQ) or ('B' in ArmStrQ) or ('b' in ArmStrQ);
    EFchk=('E' in ArmStrQ) or ('e' in ArmStrQ) or ('F' in ArmStrQ) or ('f' in ArmStrQ);
    GHchk=('G' in ArmStrQ) or ('g' in ArmStrQ) or ('H' in ArmStrQ) or ('h' in ArmStrQ);
    IJchk=('I' in ArmStrQ) or ('i' in ArmStrQ) or ('J' in ArmStrQ) or ('j' in ArmStrQ);
    chklist=[ABchk,EFchk,GHchk,IJchk];
    HWPoldrbvlist=[HWPABpv2.get(),HWPEFpv2.get(),HWPGHpv2.get(),HWPIJpv2.get()];
    for ii in range(4):
        if chklist[ii]:
            HWPpvlist[ii].put(motor_angle)
    time.sleep(2)
    for ii in range(4):
        chklist[ii] = chklist[ii] and ((np.abs(HWPoldrbvlist[ii]-HWPpv2list[ii].get()) < 0.01) and (np.abs(motor_angle-HWPpv2list[ii].get()) > 0.01));
    retryii=0;
    while any(chklist):
        retryii+=1;
        for ii in range(4):
            if chklist[ii]:
                HWPpvlist[ii].put(motor_angle);print('HWP issue detected on '+armlist[ii]+'. Re-trying...');
            time.sleep(2)
            chklist[ii] = chklist[ii] and ((np.abs(HWPoldrbvlist[ii]-HWPpv2list[ii].get()) < 0.01) and (np.abs(motor_angle-HWPpv2list[ii].get()) > 0.01)) and (retryii<3);
    if any(chklist):
        for ii in range(4):
            if chklist[ii]:
                print('Re-try on '+armlist[ii]+' failed!')
    return
    
def PFNon(ArmStrQ):
    if ArmStrQ.lower() == 'all':
        ArmStrQ = 'ABCDEFGHIJ'
    PFNmodepv=EpicsSignal('MEC:PFN:MODE')
    PFNmodepv.put(0)
    time.sleep(2)
    PFNCDenpv=EpicsSignal('MEC:PFN:CH0:ENABLE')
    PFNAenpv=EpicsSignal('MEC:PFN:CH1:ENABLE')
    PFNBenpv=EpicsSignal('MEC:PFN:CH2:ENABLE')
    PFNEenpv=EpicsSignal('MEC:PFN:CH3:ENABLE')
    PFNFenpv=EpicsSignal('MEC:PFN:CH4:ENABLE')
    PFNGenpv=EpicsSignal('MEC:PFN:CH5:ENABLE')
    PFNHenpv=EpicsSignal('MEC:PFN:CH6:ENABLE')
    PFNIenpv=EpicsSignal('MEC:PFN:CH7:ENABLE')
    PFNJenpv=EpicsSignal('MEC:PFN:CH8:ENABLE')
    if ('A' in ArmStrQ) or ('a' in ArmStrQ):
        PFNAenpv.put(1)
    if ('B' in ArmStrQ) or ('b' in ArmStrQ):
        PFNBenpv.put(1)
    if ('CD' in ArmStrQ) or ('cd' in ArmStrQ):
        PFNCDenpv.put(1)
    if ('E' in ArmStrQ) or ('e' in ArmStrQ):
        PFNEenpv.put(1)
    if ('F' in ArmStrQ) or ('f' in ArmStrQ):
        PFNFenpv.put(1)
    if ('G' in ArmStrQ) or ('g' in ArmStrQ):
        PFNGenpv.put(1)
    if ('H' in ArmStrQ) or ('h' in ArmStrQ):
        PFNHenpv.put(1)
    if ('I' in ArmStrQ) or ('i' in ArmStrQ):
        PFNIenpv.put(1)
    if ('J' in ArmStrQ) or ('j' in ArmStrQ):
        PFNJenpv.put(1)
    time.sleep(2);PFNmodepv.put(1);time.sleep(3.5);PFNmodepv.put(2);
    return

def PFNoff(ArmStrQ):
    if ArmStrQ.lower() == 'all':
        ArmStrQ = 'ABCDEFGHIJ'
    PFNmodepv=EpicsSignal('MEC:PFN:MODE')
    PFNmodepv.put(0)
    time.sleep(2)
    PFNCDenpv=EpicsSignal('MEC:PFN:CH0:ENABLE')
    PFNAenpv=EpicsSignal('MEC:PFN:CH1:ENABLE')
    PFNBenpv=EpicsSignal('MEC:PFN:CH2:ENABLE')
    PFNEenpv=EpicsSignal('MEC:PFN:CH3:ENABLE')
    PFNFenpv=EpicsSignal('MEC:PFN:CH4:ENABLE')
    PFNGenpv=EpicsSignal('MEC:PFN:CH5:ENABLE')
    PFNHenpv=EpicsSignal('MEC:PFN:CH6:ENABLE')
    PFNIenpv=EpicsSignal('MEC:PFN:CH7:ENABLE')
    PFNJenpv=EpicsSignal('MEC:PFN:CH8:ENABLE')
    if ('A' in ArmStrQ) or ('a' in ArmStrQ):
        PFNAenpv.put(0)
    if ('B' in ArmStrQ) or ('b' in ArmStrQ):
        PFNBenpv.put(0)
    if ('CD' in ArmStrQ) or ('cd' in ArmStrQ):
        PFNCDenpv.put(0)
    if ('E' in ArmStrQ) or ('e' in ArmStrQ):
        PFNEenpv.put(0)
    if ('F' in ArmStrQ) or ('f' in ArmStrQ):
        PFNFenpv.put(0)
    if ('G' in ArmStrQ) or ('g' in ArmStrQ):
        PFNGenpv.put(0)
    if ('H' in ArmStrQ) or ('h' in ArmStrQ):
        PFNHenpv.put(0)
    if ('I' in ArmStrQ) or ('i' in ArmStrQ):
        PFNIenpv.put(0)
    if ('J' in ArmStrQ) or ('j' in ArmStrQ):
        PFNJenpv.put(0)
    time.sleep(2);PFNmodepv.put(1);time.sleep(3.5);PFNmodepv.put(2);
    return

def PFNonly(ArmStrQ):#delete and just use ARMonly?
    AllStrQ='ABEFGHIJ'
    AllStrq='abefghij'
    if ArmStrQ.lower() == 'all':
        ArmStrQ = AllStrQ
    PFNmodepv=EpicsSignal('MEC:PFN:MODE')
    PFNmodepv.put(0)
    time.sleep(2)
    for ii in range(len(AllStrQ)):
        if (AllStrQ[ii] in ArmStrQ) or (AllStrq[ii] in ArmStrQ):
            temppv=EpicsSignal(str('MEC:PFN:CH'+str(ii+1)+':ENABLE'))
            temppv.put(1)
        else:
            temppv=EpicsSignal(str('MEC:PFN:CH'+str(ii+1)+':ENABLE'))
            temppv.put(0)
    time.sleep(2);PFNmodepv.put(1);time.sleep(3.5);PFNmodepv.put(2);
    return


def ARMonly(ArmStrQ,set_T=1):
    HWPon(ArmStrQ,set_T=set_T)
    PFNonly(ArmStrQ)
    return

def TTL_shutter_status(display=True):
    pvstatuslist=['MEC:LAS:FLOAT:'+str(ii) for ii in range(14,21)];
    statuslist=[]
    for eapv in pvstatuslist:
        temppv=EpicsSignal(eapv);
        statuslist.append(int(temppv.get()))#0=open,1=closed
    if display:
        shutlist=['AB','EF','GH','IJ','WEST 527','EAST 527','REGEN']
        print('(0=out,1=in) '+', '.join([ea_shut+':'+str(ea_stat) for ea_shut,ea_stat in zip(shutlist,statuslist)]))
    return statuslist

def toggle_TTL_shutter(ArmStrQ,display=True):
    AllStrq='abefghijwwxxzz';#ww is WEST 527, xx is EAST 527, zz is REGEN SHUT
    ArmStrQ=ArmStrQ.lower().replace('all','abefghijwwxx')
    currTTLstate=TTL_shutter_status(display=False)[:-1];
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
    statuslist=TTL_shutter_status(display=display)
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
    statuslist=TTL_shutter_status(display=display)
    return statuslist

def shotpush():
    temppv1=EpicsSignal('EVR:MEC:USR01:TRIG7:TEC')
    temppv2=EpicsSignal('ECS:SYS0:6:PLYCTL')
    temppv1.put(182)
    temppv2.put(1)

def YFEwatch(ShotNumQ):
    plt.ion()
    fig=plt.figure()
    ax=fig.add_subplot(211)
    line1, = ax.plot(range(len(pwtF[:26])),pwtF[:26],'b-')
    line2, = ax.plot(range(len(pwtF[:26])),ops0F[:26],'r-')
    ax2=fig.add_subplot(212)
    xq=list();yq=list();
    try:
        SLA=LXOpen('A');time.sleep(.15);
        for ii in range(ShotNumQ):
            newops0F=TraceFormatting(rch(1,SLA),[25,500],1)[:26]
            line2.set_ydata(newops0F)
            meanerr=np.sum(np.abs(pwtF[:26]-newops0F[:26])/pwtF[:26])/len(pwtF[:26])
            xq.append(ii+1);yq.append(meanerr);
            ax2.scatter(xq[-1],yq[-1]);
            fig.canvas.draw()
            time.sleep(.1)
        LXClose(SLA);time.sleep(.15);
    except:
        LXClose(SLA);
        print('Error!')
        return False

def YFEshot():
    try:
        SLA=LXOpen('A');time.sleep(.15);
        for ii in range(100):
            newops0F=TraceFormatting(rch(1,SLA),[25,500],1)[:26]
            meanerr=np.sum(np.abs(pwtF[:26]-newops0F[:26])/pwtF[:26])/len(pwtF[:26])
            if meanerr < 0.04: #.15
                time.sleep(0.1);shotpush();time.sleep(0.1);
                break
            time.sleep(.1)
        newops0F=TraceFormatting(rch(1,SLA),[25,500],1)[:26]
        meanerr=np.sum(np.abs(pwtF[:26]-newops0F[:26])/pwtF[:26])/len(pwtF[:26])
        print(meanerr)
        epll([pwtF,newops0F])
        LXClose(SLA);time.sleep(.15);
        psacqx()
        psefc()
    except:
        LXClose(SLA);
        print('Error!')
        return False

def YFEerr():
    xq=list();eq=list();
    try:
        SLA=LXOpen('A');
        ii=0
        time.sleep(.15);
        while True:
            ii+=1
            if ii%1000 == 0:
                print(str('Iter: '+str(ii)))
            NewTraceQ=readchan(1,SLA)['DATA']
            newops0F=TraceFormatting(NewTraceQ,[25,500],1)[:26]
            meanerr=np.sum(np.abs(pwtF[:26]-newops0F[:26])/pwtF[:26])/len(pwtF[:26])
            xq.append(meanerr);
            eq.append(EG1wYFE1in()[0])
    except KeyboardInterrupt:
        LXClose(SLA);time.sleep(.15);
        return xq,eq
    finally:
        LXClose(SLA);
        print('Error!')
        return False

def efft(errlistQ,time_stepQ):
    #time_step1=0.1
    freqs1=np.fft.fftfreq(np.array(errlistQ).size, time_stepQ)
    idx1=np.argsort(freqs1)
    fftd1=np.fft.fft(errlistQ)
    ps1=np.abs(fftd1/max(fftd1))**2
    eplxyloglog(freqs1[idx1],ps1[idx1])
    return [freqs1[idx1],ps1[idx1]]








#end algo2.py

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################



def PulseGoal(DurationListQ,StartStopListQ):#140 pt list with max at 1
    BeginPix=1#was 51
    #PulseGoal([3,5],[[15,25],[75,100]],['l','l'])
    #PulseGoal([3,5],[[15,25],[75,100]],['l','p'])
    DurListQ=np.cumsum([0]+DurationListQ)
    SSListQ=StartStopListQ[:]
    SegTotQ=len(DurListQ)-1
    if (-1+len(DurListQ))!=len(SSListQ):
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
    for ii in range((-1+len(DurListQ))):
        #SegmentsQ.append(LinearWave(51+(DurListQ[ii]*4),int(20000.*SSListQ[ii][0]/100.),51+(DurListQ[ii+1]*4),int(20000.*SSListQ[ii][1]/100.)))
        SegmentsQ.append(LinearWave(int(BeginPix+(DurListQ[ii]*4)),int(20000.*SSListQ[ii][0]/100.),int(BeginPix+(DurListQ[ii+1]*4)-1),int(20000.*SSListQ[ii][1]/100.)))
    return np.append(np.delete(np.array(ComboWave(SegmentsQ)),np.array(DelListQ).astype(int)),[0]*len(DelListQ))

def PulseMax(DurationListQ,StartStopListQ,zzJQ):#fixed normalization part
    #Get amplitude setting for segmented, arbitrary PulseGoal
    return (1.*StartStopListQ[-1][-1]/100.)*(50.*zzJQ/(5.*500.*np.sum(PulseGoal(DurationListQ,StartStopListQ))))

def wIter2(sQ,wQ,DurationListQ,StartStopListQ,zzJQ,mapnowQ,stepQQ):#same as above but with TraceFormatting2, which has less averaging
    DurListQ=np.cumsum([0]+DurationListQ)
    w1,w2=0,int(DurListQ[-1]*4)+5 # 50-5, 50+int(DurListQ[-1]*4)+5
    PGQ=PulseGoal(DurationListQ,StartStopListQ)
    if np.abs(len(sQ)-10000)<10:
        PMQcorr=1
    elif np.abs(len(sQ)-1000)<10:
        PMQcorr=10
    else:
        print('Warning: unanticipated pulse shape array length of '+str(len(templistQ))+', PulseMax scaling may be off...')
        PMQcorr=1
    PMQ=PulseMax(DurationListQ,StartStopListQ,zzJQ)*PMQcorr
    wnew2=FixEdges(UpdatingShapingAlgorithm(PGQ,TraceFormatting2(sQ,mapnowQ,PMQ),wQ,stepQQ),DurationListQ,StartStopListQ)
    #epll([wQ[w1:w2],wnew2[w1:w2],np.array(TraceFormatting2(sQ,mapnowQ,PMQ))[w1:w2]*.6,np.array(PGQ)[w1:w2]*.6])
    epll([0.*np.array(wnew2[w1:w2]),np.array(wnew2[w1:w2])-np.array(wQ[w1:w2]),np.array(TraceFormatting2(sQ,mapnowQ,PMQ))[w1:w2]*.6,np.array(PGQ)[w1:w2]*.6])
    return wnew2








def pshostcheck():
    try:
        hostname=socket.gethostname()
        if (hostname != 'mec-monitor') and (hostname != 'mec-daq') and (hostname != 'mec-laser'):
            print('Host must be mec-monitor, mec-daq, or mec-laser! Current host: '+hostname)
            raise Exception
    except Exception:
        print('EXIT')
        #exit
        #sys.exit();
    try:
        curruser=getpass.getuser()
        if curruser != 'mecopr':
            print('Warning: you are logged in as '+curruser+', not as mecopr. Beware of permissions issues... You may even unwittingly cause some yourself!')
    except:
        print('Failed: could not ID current user!')
    return

def LMap():
    return [50,1000]

def DateString():
    qdate=date.today()
    return qdate.strftime('%Y%m%d')

def psfilepath():
    return '/reg/neh/operator/mecopr/mecpython/pulseshaping/'

def get_curr_exp(timeout=15): 
    script=SCRIPTS.format('mec','get_curr_exp') 
    exp=cache_script(script,timeout=timeout) 
    curr_exp=exp.lower().strip('\n')
    try:
        temppv=EpicsSignal('MEC:LAS:FLOAT:11.DESC')
        temppv.put(curr_exp)
    except:
        temppv.put('mecxx####')
    return curr_exp


def psheaders():
    mapnow=LMap()#line up fiducials of +45ns (40ns after edge) at pix no 9000#was [5,153]
    DateStr=DateString()
    psfpQ=psfilepath()
    try:
        globals()['w'+DateStr] = pickle.load(open(psfpQ+'w'+DateStr+'.p','rb'))
        globals()['y'+DateStr] = pickle.load(open(psfpQ+'y'+DateStr+'.p','rb'))
        globals()['s1in1w'+DateStr] = pickle.load(open(psfpQ+'s1in1w'+DateStr+'.p','rb'))
        globals()['s42in1w'+DateStr] = pickle.load(open(psfpQ+'s42in1w'+DateStr+'.p','rb'))
        globals()['s42in2w'+DateStr] = pickle.load(open(psfpQ+'s42in2w'+DateStr+'.p','rb'))
        globals()['s'+DateStr] = pickle.load(open(psfpQ+'s'+DateStr+'.p','rb'))
    except:
        print('No laser file found -- probably the first shot of the day.')
        try:
            globals()['s'+DateStr]
        except KeyError:#was NameError
            globals()['w'+DateStr]=[]
            globals()['y'+DateStr]=[]
            globals()['s1in1w'+DateStr]=[]
            globals()['s42in1w'+DateStr]=[]
            globals()['s42in2w'+DateStr]=[]
            globals()['s'+DateStr]=[]
    
    DateStrArr=[]
    exec('DateStrArr.append(w'+DateStr+')')
    exec('DateStrArr.append(y'+DateStr+')')
    exec('DateStrArr.append(s1in1w'+DateStr+')')
    exec('DateStrArr.append(s42in1w'+DateStr+')')
    exec('DateStrArr.append(s42in2w'+DateStr+')')
    exec('DateStrArr.append(s'+DateStr+')')
    globals()['wtoday']=DateStrArr[0]
    globals()['ytoday']=DateStrArr[1]
    globals()['s1in1wtoday']=DateStrArr[2]
    globals()['s42in1wtoday']=DateStrArr[3]
    globals()['s42in2wtoday']=DateStrArr[4]
    globals()['stoday']=DateStrArr[5]
    try:
        #globals()['RunNum']
        globals()['RunNum']=get_run_number(hutch='mec',timeout=10)
    except:
        #print('No RunNum given, set to 9000')
        globals()['RunNum']=9000
    try:
        globals()['RunFilePath']
    except:
        #print('No RunFilePath given, set to '+psfpQ+'temp/')
        globals()['RunFilePath']=psfpQ+'temp/'
    try:
        temppv=EpicsSignal('MEC:LAS:FLOAT:11')
        temppv.put(globals()['RunNum'])
    except:
        temppv.put(-1)
    return DateStr
        


def psacqx(save_flag=True, RunNumQQ=9000):#was psacqx_noLecroyA()
    pshostcheck()
    DateStr=psheaders()
    psfpQ=psfilepath()
    if save_flag:
        try:
            globals()['RunNum']=get_run_number(hutch='mec',timeout=10)
        except:
            print('Failed to find the run number')
    if RunNumQQ!=9000:
        globals()['RunNum']=RunNumQQ
    try:
        RunNuQ=globals()['RunNum']
    except:
        #print('No RunNum given, set to 9000')
        globals()['RunNum']=9000
        RunNuQ=globals()['RunNum']
    try:
        temppv=EpicsSignal('MEC:LAS:FLOAT:11')
        temppv.put(globals()['RunNum'])
    except:
        temppv.put(-1)
    try:
        fpQ=globals()['RunFilePath']
    except:
        #print('No RunFilePath given, set to '+psfpQ+'temp/')
        globals()['RunFilePath']=psfpQ+'temp/'
        fpQ=globals()['RunFilePath']
    RunNumStr=str(RunNuQ).zfill(4)#changed to 4, just in case! 20210810
    RunName='run'+str(RunNumStr)+'_'# ##add 'test' at beginning for test
    #fpQ='/reg/neh/operator/mecopr/experiments/optical_beam/ns_laser/pp/spare/'
    #fpQ='/reg/neh/operator/mecopr/experiments/mecx37917/lecroy/'
    #fpQ='/reg/neh/operator/mecopr/experiments/meclu6517/lecroy/'
    #fpQ='/reg/neh/operator/mecopr/experiments/meck05717/lecroy/'
    #fpQ='/reg/neh/operator/mecopr/experiments/meclv8018/lecroy/'
    #ExpName='meclv8018'
    if save_flag:
        ExpName=get_curr_exp()
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
    
    [cAB,cEF,cGH,cIJ] = LasCoeff()
##    cAB=pypsepics.get('MEC:PFN:CH1:ENABLE_RBV')*pypsepics.get('MEC:PFN:CH2:ENABLE_RBV')
##    cEF=pypsepics.get('MEC:PFN:CH3:ENABLE_RBV')*pypsepics.get('MEC:PFN:CH4:ENABLE_RBV')
##    cGH=pypsepics.get('MEC:PFN:CH5:ENABLE_RBV')*pypsepics.get('MEC:PFN:CH6:ENABLE_RBV')
##    cIJ=pypsepics.get('MEC:PFN:CH7:ENABLE_RBV')*pypsepics.get('MEC:PFN:CH8:ENABLE_RBV')
    if cAB:
        RunName=RunName+'AB'
    if cEF:
        RunName=RunName+'EF'
    if cGH:
        RunName=RunName+'GH'
    if cIJ:
        RunName=RunName+'IJ'

    try:
        SLA=LXOpen('A'); time.sleep(.15);#usually LXOpen('A');
        weichYFE00=np.array(weichYFE1w(1,SLA));
        weich1in1wCD=np.array(weich1in1w(2,SLA));
        time.sleep(.15); LXClose(SLA); time.sleep(.15);
    except:
        try:
            LXClose(SLA);
        except:
            print('Could not close socket!')
        weichYFE00=np.zeros(5002);weich1in1wCD=np.zeros(5002);
        print('Error! SLA')

    try:
        SLB=LXOpen('B'); time.sleep(.15);#usually LXOpen('B');
        weich2in1wAB=np.array(weich2in1w(1,SLB)); weich2in1wEF=np.array(weich2in1w(2,SLB));
        weich2in1wGH=np.array(weich2in1w(3,SLB)); weich2in1wIJ=np.array(weich2in1w(4,SLB));
        time.sleep(.15); LXClose(SLB); time.sleep(.15);
    except:
        try:
            LXClose(SLB);
        except:
            print('Could not close socket!')
        weich2in1wAB=np.zeros(5002); weich2in1wEF=np.zeros(5002);weich2in1wGH=np.zeros(5002); weich2in1wIJ=np.zeros(5002);
        print('Error! SLB')

    #SL2=L2Open(); time.sleep(.15);
    #weich1=np.array(weich(1,SL2)); weich2=np.array(weich(2,SL2));
    #weich3=np.array(weich(3,SL2)); weich4=np.array(weich(4,SL2));
    #time.sleep(.15); L2Close(SL2); time.sleep(.15);

    try:
        SL2=LXOpen('2'); time.sleep(.15);#usually LXOpen('2');
        weich2in2wAB=np.array(weich2in2w(1,SL2)); weich2in2wEF=np.array(weich2in2w(2,SL2));
        weich2in2wGH=np.array(weich2in2w(3,SL2)); weich2in2wIJ=np.array(weich2in2w(4,SL2));
        time.sleep(.15); LXClose(SL2); time.sleep(.15);
    except:
        try:
            LXClose(SL2);
        except:
            print('Could not close socket!')
        print('Error! SL2')

    total_print=''
    #check so don't overwrite if someone forgets to change run number

    try:
        olden=np.genfromtxt(fpQ+RunName+'energies.txt')
        RunName=str(RunName+'DUPLICATE')
        print(str('This run number already exists; run name '+RunName+' used'))
        total_print+=str('This run number already exists; run name '+RunName+' used')
        total_print+='\n'
    except:
        print(str('Run name: '+RunName+', shot number: '+str(len(stoday)+1)))
        total_print+=str('Run name: '+RunName+', shot number: '+str(len(stoday)+1))
        total_print+='\n'

    print(datetime.now().strftime('%A, %d. %B %Y %I:%M%p'))
    total_print+=str(datetime.now().strftime('%A, %d. %B %Y %I:%M%p'))
    total_print+='\n'

    if not save_flag:
        print('This is a test run. Use save_flag=True if you want to save data.')

    if save_flag:
        np.savetxt(str(fpQ+RunName+'ch1.txt'),weichToPowerVsTime(weich2in2wAB))
        np.savetxt(str(fpQ+RunName+'ch2.txt'),weichToPowerVsTime(weich2in2wEF))
        np.savetxt(str(fpQ+RunName+'ch3.txt'),weichToPowerVsTime(weich2in2wGH))
        np.savetxt(str(fpQ+RunName+'ch4.txt'),weichToPowerVsTime(weich2in2wIJ))

    WeightedSum=cAB*weich2in2wAB+cEF*weich2in2wEF+cGH*weich2in2wGH+cIJ*weich2in2wIJ

    if save_flag:
        np.savetxt(str(fpQ+RunName+'chsum.txt'),weichToPowerVsTime(WeightedSum))

    globals()['PulseEnergies']=np.array(list(map(lambda x: round(x,2),EG()[0][0])))*np.array([cAB,cEF,cGH,cIJ])
    EnMess='***'
    EnMess2=[' AB: ',' EF: ',' GH: ',' IJ: ']
    for ii in range(len(PulseEnergies)):
        EnMess+=EnMess2[ii]
        EnMess+=str(PulseEnergies[ii])
        EnMess+=' J ***'
    EnMess+=str(' total: '+str(np.sum(PulseEnergies))+' J ***')
    #print(EnMess)
    total_print+=EnMess
    total_print+='\n'#start over...
    total_print='Run: '+str(RunNuQ)+'\n'
    
    headenb = LasCoeff()
    wpreadpvlist=['MEC:NS1:MMS:02.RBV','MEC:NS1:MMS:01.RBV','MEC:LAS:MMN:30.RBV','MEC:LAS:MMN:29.RBV']
    wpwritepvlist=['MEC:NS1:MMS:02.VAL','MEC:NS1:MMS:01.VAL','MEC:LAS:MMN:30.VAL','MEC:LAS:MMN:29.VAL']
    wppvlist=[EpicsSignal(wpreadpvlist[ii],write_pv=wpwritepvlist[ii]) for ii in range(4)]
    headstr='';wpstr='';headlist=['AB','EF','GH','IJ'];
    for ii in range(4):
        if headenb[ii]:
            headstr+=headlist[ii]
        wpstr=wpstr+headlist[ii]+': '+str(round(wppvlist[ii].get(),3))+', '
    wpen=np.mean([headenb[ii]*np.cos((np.pi/180)*2*wppvlist[ii].get())**2 for ii in range(4)])
    total_print+=str('The following heads are enabled: '+headstr+'\n')
    total_print+=str('The waveplate settings are: '+wpstr[:-2]+'\n')
    total_print+=str('This is ~'+str(round(100*wpen,3))+'% of max energy.'+'\n')

    total_print+=EGall(return_txt=True);


    #EGall();
    #epl(TraceFormatting(WeightedSum,mapnow,1))

    if save_flag:
        np.savetxt(str(fpQ+RunName+'energies.txt'),PulseEnergies)
        #def eplxysav(listxq,listyq,FileNameQ,abs_path=False,xlb='none',ylb='none'):
        eplxysav(*weichToPowerVsTime(WeightedSum),str(fpQ+RunName+'_'+str(int(round(np.sum(PulseEnergies))))+'J'),abs_path=True,xlb='Time (ns)',ylb='Power (W)')
#        try:
#            SL1=LOpen();time.sleep(0.15);
#            TCCDiodeTrace=rch(2,SL1);time.sleep(0.15);
#            LClose(SL1);time.sleep(0.15);
#            SLA=LAOpen();time.sleep(0.15);
#            fsDiodeTrace=rch(3,SLA);time.sleep(0.15);
#            LAClose(SLA);time.sleep(0.15);
#            np.savetxt(str(fpQ+RunName+'_TCC+YFE+fsdiode.txt'),(TCCDiodeTrace,weichYFE00,fsDiodeTrace))
#        except:
#            LClose(SL1);time.sleep(0.15);
#            LAClose(SLA);time.sleep(0.15);

    try:
        SH=HOpen(); time.sleep(.15);
        wtoday.append(ReadPulseHeights(SH,0));
        time.sleep(.15); HClose(SH); time.sleep(.15);
    except:
        HClose(SH);
        print('Error! SH')

    #SLA=LAOpen(); time.sleep(.15); 
    #htoday.append(rch(1,SLA))
    #ytoday.append(rch(1,SLA));
    #LAClose(SLA); time.sleep(.15);
    ytoday.append(weichYFE00)
    s1in1wtoday.append(weich1in1wCD)
    s42in1wtoday.append([weich2in1wAB,weich2in1wEF,weich2in1wGH,weich2in1wIJ])
    s42in2wtoday.append([weich2in2wAB,weich2in2wEF,weich2in2wGH,weich2in2wIJ])
    stoday.append(WeightedSum)

    pickledump2(wtoday,str(psfpQ+'w'+str(DateStr)+'.p'));

    #pickle.dump(htoday,open(str(psfpQ+'h'+str(DateStr)+'.p'),'wb'));
    pickledump2(ytoday,str(psfpQ+'y'+str(DateStr)+'.p'));
    pickledump2(s1in1wtoday,str(psfpQ+'s1in1w'+str(DateStr)+'.p'));
    pickledump2(s42in1wtoday,str(psfpQ+'s42in1w'+str(DateStr)+'.p'));
    pickledump2(s42in2wtoday,str(psfpQ+'s42in2w'+str(DateStr)+'.p'));
    pickledump2(stoday,str(psfpQ+'s'+str(DateStr)+'.p'));

    if save_flag:
        mecel = elog.ELog({'experiment':ExpName},user='mecopr',pw=pickle.load(open(psfilepath()+'elogauth.p','rb'))) 
        #mecElog.submit(total_print,file=str(fpQ+RunName+'_'+str(int(round(np.sum(PulseEnergies))))+'J.png'))
        if True: #RunNuQ != 9000: ##FIX LATER MAYBE IF DAQ WORKING AGAIN???
            try:
                mecel.post(total_print,attachments=[str(fpQ+RunName+'_'+str(int(round(np.sum(PulseEnergies))))+'J.png')],run=RunNuQ,tags=['laser'])
                print('Auto-saved to eLog with run '+str(RunNuQ)) 
            except:
                try:
                    mecel.post(total_print,attachments=[str(fpQ+RunName+'_'+str(int(round(np.sum(PulseEnergies))))+'J.png')],tags=['laser']) 
                    print('Auto-saved to eLog')
                except:
                    print('Failed to auto-save to eLog')
            

    #execfile('tektest.py')

    globals()['RunNum']+=1
    
    ##EXECUTE THIS FILE AFTER A SHOT TO SAVE YOUR SCOPE TRACE



def psefc(JreqQ=0,AQQ=0.0):
    pshostcheck()
    psheaders();
    mapnow=LMap()

    psfpQ=psfilepath()
    Psns=pickle.load(open(psfpQ+'Psns.p','rb'))
    SSs=pickle.load(open(psfpQ+'SSs.p','rb'))
    Jscale=1

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
            mapnow=LMap()
        elif len(stoday[-1]) == 1002:
            mapnow=[5,100]
        else:
            print('Unanticipated pulse shot array length: '+str(len(stoday[-1])));
            print('Aborting...');return
        wupd=wIter2(stoday[-1],np.array(wtoday[-1])/28000.,Psns,SSs,Jreq,mapnow,AQ)
    else:
        print('No shots yet today; beginning with pre-loaded shape')
        try:
            SH=HOpen(); time.sleep(.15);
            wupd=ReadPulseHeights(SH,0); time.sleep(.15);
            HClose(SH); time.sleep(.15);
        except:
            HClose(SH);
            print('Error! SH')
    return wupd
    ##EXECUTE THIS FILE FIRST TO DETERMINE THE UPDATED WAVEFORM

def psefc10Hz(pwt=0,numIterQ=50,AQQ=0.03,displayPlot=False,reloopPrompt=False):#was B; replaced when LeCroyA fixed
    evrpv=EpicsSignal('EVR:MEC:USR01:TRIG7:TEC')
    evrpv.put(43)
    psfpQ=psfilepath()   
    SSs=pickle.load(open(psfpQ+'SSs.p','rb'))
    Psns=pickle.load(open(psfpQ+'Psns.p','rb'))
    pwtF=np.array(TraceFormatting2(pwt,[25,500],1))


    try:
        SLA=LXOpen('A');time.sleep(.15);S=HOpen();time.sleep(.15);#replaced w/LeCroyA
        ops00=rch(1,SLA);time.sleep(0.1);

        meanerr=[]
        meanerr.append(np.sum(np.abs(pwtF[:26]-TraceFormatting2(ops00,[5,100],1)[:26])/pwtF[:26])/len(pwtF[:26]))
        ops00F=TraceFormatting2(ops00,[5,100],1)
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
                ops0=rch(1,SLA);time.sleep(0.025);####added.215 when 200mV/div instead of 100mV/div
                if all(ops0 == ops00):
                    print('No scope update detected... no feedback applied!')
                else:
                    rph=ReadPulseHeights(S,0);time.sleep(0.025);
                    #pwtF=np.array(TraceFormatting2(pwt,[25,500],1))
                    ##ops0F=TraceFormatting(ops0,[25,500],1)
                    ops0F=TraceFormatting2(ops0,[5,100],1)
                    #epll([pwtF,ops0F])
                    meanerr.append(np.sum(np.abs(pwtF[:26]-ops0F[:26])/pwtF[:26])/len(pwtF[:26]));
                    if displayPlot:
                        axss[0].set_data(xdat,ops0F); axs[0].relim(); axs[0].autoscale_view(True,True,True);
                        axss[1].set_data(list(range(len(meanerr))),meanerr); axs[1].relim(); axs[1].autoscale_view(True,True,True);
                        fig.canvas.draw_idle(); plt.pause(0.01);
                    usa0=UpdatingShapingAlgorithm(pwtF,ops0F,np.array(rph)/28000.,AQQ)#.075#.25
                    usa0FE=FixEdges(usa0,Psns,SSs)
                    #usa0FE=FixEdges(usa0,[3,4.25],[[.98*100/8.0,100/8.0],[98,100]])
                    #epll([rph,usa0FE*28000.])
                    WritePulseHeights(S,0,usa0FE*28000.);time.sleep(0.05);
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
        HClose(S);time.sleep(.15);LXClose(SLA);time.sleep(.15)#replace w/ LeCroyA
    except:
        print('Failed')
        HClose(S);time.sleep(.15);LXClose(SLA);time.sleep(.15)#replace w/ LeCroyA

def psupd(newwavQ):
    ###EXECUTE THIS FILE ONCE YOU'RE SATISFIED WITH THE WAVEFORM UPDATE
    pshostcheck()
    wupdt=newwavQ[:]
    if max(wupdt) < 1.5:
        wupdt=28000.0*np.array(wupdt)
    try:
        SH=HOpen(); time.sleep(.15);
        WritePulseHeights(SH,0,wupdt); time.sleep(.15);
        HClose(SH); time.sleep(.15);
    except:
        HClose(SH);
        print('Error!')

def psall(upd_loop=True):
    pshostcheck()
    try:
        #RunNuQ=globals()['RunNum']
        psacqx()
        wupdQ=psefc()
        if upd_loop:
            psupd(wupdQ)
    except:
        try:
            globals()['RunNum']
        except:
            print('Current run number has not yet been entered! Please enter the current run number in this terminal like this (for example):\n \n RunNum=64')


def psloadwvfm(RecipeStrQ,WvGoal10HzHint=False):
    pshostcheck()
    print('Loading timestamp: '+datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
    try:
        [Psns,SSs,YFE02mmCurr,YFE06mmCurr,YFE10mmCurr,NewWvfm,WvGoal10Hz] = pickle.load(open(psfilepath()+'recipes/load'+RecipeStrQ+'.p','rb'))
    except:
        print('Recipe file '+psfilepath()+'recipes/load'+RecipeStrQ+'.p\' not found.')
        return
    pickledump2(Psns,psfilepath()+'Psns.p')
    pickledump2(SSs,psfilepath()+'SSs.p')
    YFEset(2,YFE02mmCurr);time.sleep(.15);
    YFEset(6,YFE06mmCurr);time.sleep(.15);
    YFEset(10,YFE10mmCurr);time.sleep(.15);
    if WvGoal10HzHint:
        if len(WvGoal10Hz) > 5:
            print('Hint for 10Hz waveform: '+WvGoal10Hz)
        else:
            print('No hint for 10Hz waveform available. MSG: \''+WvGoal10Hz+'\'')
    try:
        psupd(NewWvfm)
        print('New waveform loaded! ')
    except:
        print('Failed to load new waveform.')
    



def pssavewvfm(PsnsQ=0,SSsQ=0,YFEgetQ=0,TargetwlistDateQ='curr',TargetwindexQ=0,RecipeStrQ=0,WvGoal10HzQ='none'):
    pshostcheck()
    print('Saving timestamp: '+datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
    if TargetwlistDateQ == 'curr':
        print('Using current Highland waveform...')
        try:
            S=HOpen(); time.sleep(0.15);
            NewWvfmQ=ReadPulseHeights(S,0);
            HClose(S); time.sleep(0.15);
        except:
            HClose(S);
        try:
            wlastarr=pickle.load(open(psfilepath()+'w'+DateString()+'.p','rb'))
            if wlastarr[-1] == NewWvfmQ:
                print('Pulse looks equivalent to the most recent pulse, w'+DateString()+'['+str(len(wlastarr)-1)+']')
                WvGoal10HzQ=WvGoal10HzQ+';; w'+DateString()+'['+str(len(wlastarr)-1)+']'
            else:
                WvGoal10HzQ=WvGoal10HzQ+';; sometime after most recent w'+DateString()+'['+str(len(wlastarr)-1)+']'
        except:
            print('Failed to load most recent amplified shot.')
    else:
        wavehistQ=pickle.load(open(psfilepath()+'w'+TargetwlistDateQ+'.p','rb'))
        NewWvfmQ=wavehistQ[TargetwindexQ][:]
        WvGoal10HzQ=WvGoal10HzQ+', w'+TargetwlistDateQ+'['+str(TargetwindexQ)+']'

    if PsnsQ == 0:
        PsnsQ = pickle.load(open(psfilepath()+'Psns.p','rb'))
    if SSsQ == 0:
        SSsQ = pickle.load(open(psfilepath()+'SSs.p','rb'))
    if YFEgetQ == 0:
        YFEgetQ=YFEget(display=False)
    [YFE02mmCurrQ,YFE02mmCurrQ,YFE02mmCurrQ,YFE02mmCurrQ,YFE06mmCurrQ,YFE10mmCurrQ]=YFEgetQ
    if RecipeStrQ == 0:
        #learn formatting for pulse from Psns and SSs
        pass
    oldfilefound=True

    try:
        dummyQ = pickle.load(open(psfilepath()+'recipes/load'+RecipeStrQ+'.p','rb'))
        print('Old recipe found with same name: '+psfilepath()+'recipes/load'+RecipeStrQ+'.p')
    except:
        oldfilefound=False
    iiQ=0
    while oldfilefound:
        iiQ=iiQ+1
        try:
            dummyQ2 = pickle.load(open(psfilepath()+'recipes/load'+RecipeStrQ+'_'+str(iiQ).zfill(2)+'.p','rb'))
        except:
            oldfilefound=False
            dummyQ[-1]='**replaced on '+DateString()+'** '+dummyQ[-1]
            pickledump2(dummyQ,psfilepath()+'recipes/load'+RecipeStrQ+'_'+str(iiQ).zfill(2)+'.p')
            print('Saved old recipe as '+psfilepath()+'recipes/load'+RecipeStrQ+'_'+str(iiQ).zfill(2)+'.p')
    try:
        pickledump2([PsnsQ,SSsQ,YFE02mmCurrQ,YFE06mmCurrQ,YFE10mmCurrQ,NewWvfmQ,WvGoal10HzQ],psfilepath()+'recipes/load'+RecipeStrQ+'.p')
        print('Saved new recipe as '+psfilepath()+'recipes/load'+RecipeStrQ+'.p')
    except:
        print('Failed to save new recipe.')



def psviewwvfm(RecipeStrQ='none',TargetwlistDateQ='curr',TargetwindexQ=0,WvGoal10HzHint=False):
    pshostcheck()
    if RecipeStrQ == 'none':
        if TargetwlistDateQ == 'curr':
            foundlastshot=False
            iidQ=0
            while not foundlastshot:
                try:
                    wlastarr=pickle.load(open(psfilepath()+'w'+str(int(DateString())-iidQ)+'.p','rb'))
                    slastarr=pickle.load(open(psfilepath()+'s'+str(int(DateString())-iidQ)+'.p','rb'))
                    wlast=wlastarr[-1][:]
                    slast=slastarr[-1][:]
                    print('Retrieving most recent shot: w'+str(int(DateString())-iidQ)+'['+str(len(wlastarr)-1)+']')
                    foundlastshot=True
                except:
                    iidQ=iidQ+1
        else:
            try:
                wlastarr=pickle.load(open(psfilepath()+'w'+TargetwlistDateQ+'.p','rb'))
                slastarr=pickle.load(open(psfilepath()+'s'+TargetwlistDateQ+'.p','rb'))
                wlast=wlastarr[TargetwindexQ]
                slast=slastarr[TargetwindexQ]
            except:
                print('Failed to load at given date and index: '+TargetwlistDateQ+', '+str(TargetwindexQ))
        epl(wlast)
        epl(slast)
    else:
        try:
            [Psns,SSs,YFE02mmCurr,YFE06mmCurr,YFE10mmCurr,NewWvfm,WvGoal10Hz] = pickle.load(open(psfilepath()+'recipes/load'+RecipeStrQ+'.p','rb'))
        except:
            print('Recipe file '+psfilepath()+'recipes/load'+RecipeStrQ+'.p\' not found.')
            return
        print('Retrieved recipe: load'+RecipeStrQ)
        print('Psns: '+str(Psns)+', SSs: '+str(SSs))
        print('YFEcurr:: 2mm: '+'{:5.1f}'.format(YFE02mmCurr)+', 6mm: '+'{:5.1f}'.format(YFE06mmCurr)+', 10mm: '+'{:5.1f}'.format(YFE10mmCurr))
        print('Extra info: '+WvGoal10Hz)
        epl(NewWvfm)
        try:
            tempstr=WvGoal10Hz[-18:]
            wstrinx=tempstr.index('w')
            wstrinx0=wstrinx-len(tempstr)
            loaddate=tempstr[wstrinx0+1:wstrinx0+9]
            loadindx=int(tempstr[wstrinx0+10:-1])
            saskarr=pickle.load(open(psfilepath()+'s'+loaddate+'.p','rb'))
            epl(saskarr[loadindx])
        except:
            print('Failed to load 2w waveform for display.')
        return NewWvfm


def psrefrwvfm(RecipeStrQ,numStepsQ=50,stepSizeQ=0.25,displayPlotQ=False,reloopPrompt=False):
    pshostcheck()
    pvMBCpower=EpicsSignal('MEC:64B:PWR:2:Outlet:8:SetControlAction')#WAS 'MEC:S60:PWR:01:Outlet:7:SetControlAction'#read AND write:1=ON,2=OFF
    pvMBCmode=EpicsSignal('MEC:LPL:MBC:01:RunningMode_RBV',write_pv='MEC:LPL:MBC:01:RunningMode')#AUTO=0,MAN=1
    pvMBCsetpt=EpicsSignal('MEC:LPL:MBC:01:AutoCalibration.VAL',write_pv='MEC:LPL:MBC:01:AutoCalibration') #QUAD=0,MIN=1,MAX=2
    pvMBCbias=EpicsSignal('MEC:LPL:MBC:01:BiasValue_RBV',write_pv='MEC:LPL:MBC:01:BiasValue')
    pvMBCfault=EpicsSignal('MEC:LPL:MBC:01:ErrorStatus',write_pv='MEC:LPL:MBC:01:ClearErrors')
    pvlampEC=EpicsSignal('EVR:MEC:USR01:TRIG6:EC_RBV',write_pv='EVR:MEC:USR01:TRIG6:TEC')#lamp event code; needs 182
    pvslicerEC=EpicsSignal('EVR:MEC:USR01:TRIG7:EC_RBV',write_pv='EVR:MEC:USR01:TRIG7:TEC')#slicer event code; needs 182
    pvlampenable=EpicsSignal('EVR:MEC:USR01:TRIG6:TCTL') #lamp enable;
    pvslicerenable=EpicsSignal('EVR:MEC:USR01:TRIG7:TCTL') #slicer enable;
    print('Loading timestamp: '+datetime.now().strftime('%A, %d. %B %Y %I:%M:%S%p'))
    #load and extract the pulse target from the desired recipe
    try:
        [Psns,SSs,YFE02mmCurr,YFE06mmCurr,YFE10mmCurr,NewWvfm,WvGoal10Hz] = pickle.load(open(psfilepath()+'recipes/load'+RecipeStrQ+'.p','rb'))
    except:
        print('Recipe file '+psfilepath()+'recipes/load'+RecipeStrQ+'.p\' not found.')
        return
    print('Hint text: '+WvGoal10Hz)
    pseparams=np.array(re.findall('ExponentialWave2\((\d+),(\d+\.\d+|\.\d+),(\d+),(\d+\.\d+|\.\d+),0,5002\)',WvGoal10Hz),dtype=np.float32);
    pslparams=np.array(re.findall('LinearWave2\((\d+),(\d+\.\d+|\.\d+),(\d+),(\d+\.\d+|\.\d+),0,5002\)',WvGoal10Hz),dtype=np.float32);  
    yfegoal=LinearWave2(500,0,1025,0,0,5002);
    if (len(pslparams) > 0) or (len(pseparams) > 0):
        for ii in range(len(pslparams)): 
            yfegoal+=LinearWave2(pslparams[ii][0],pslparams[ii][1],pslparams[ii][2],pslparams[ii][3],0,5002)
        for ii in range(len(pseparams)): 
            yfegoal+=ExponentialWave2(pseparams[ii][0],pseparams[ii][1],pseparams[ii][2],pseparams[ii][3],0,5002)
    else:
        print('No wave extracted: '+WvGoal10Hz)
        return
    #close the shutters
    #print('NOT closing the shutters... hope you\'re protecting your sample!')
    print('Closing all shutters...')
    toggle_TTL_shutter('closeall',display=False)#close all the shutters
    time.sleep(4)
    try:#used to make sure shutters re-open even in case of error or KeyboardInterrupt
        #check if laser is on
        if np.sum(YFEget(display=False)) < 400:
            if np.sum(YFEget(display=False)) < 20:
                print('WARNING: YFE seems to be off... Attempting to turn on YFE...')
                YFEon();
            else:
                print('WARNING: eDrive currents seem low...')
                prechk = False
        #turn off bias dither; AUTO is 0; MAN is 1
        if pvMBCmode.get() != 1:
            print('Turning off MBC bias dither...',end='',flush=True);
            currbias=pvMBCbias.get();
            pvMBCmode.put(1);#set bias mode to MAN
            dotsleep(3);
            pvMBCbias.put(currbias+10);
            print('Testing bias responsivity...',end='',flush=True);
            for ii in range(2):
               time.sleep(1);print('..',end='',flush=True);
            if np.abs(pvMBCbias.get() - (currbias+10)) < 3:
                pvMBCbias.put(currbias);
                dotsleep(2);
            else:
                print('*')
                print('WARNING: MBC not responding!!')
        #else:
        #    print('MBC is not safe! Resetting the MBC...')#problem is here...
        #    resetMBC();
        #set and enable 10Hz output
        if pvslicerEC.get() != 43:
            pvslicerenable.put(0)
            pvslicerEC.put(43)
        if pvslicerenable.get() != 1:
            pvslicerenable.put(1)
        #run the update code
        print('Refreshing the YFE wavefront...')
        psefc10Hz(pwt=yfegoal,numIterQ=numStepsQ,AQQ=stepSizeQ,displayPlot=displayPlotQ,reloopPrompt=reloopPrompt)
        #reset to single shot on pulse picker
        pvslicerenable.put(0);time.sleep(0.5);
        pvslicerEC.put(182);time.sleep(0.5);
        pvslicerenable.put(1)
        #re-open shutters
        #print('Opening all shutters...')
        #toggle_TTL_shutter('openall',display=False);#open all the shutters
        resetMBC();
        YFEsetall(True,displayQ=False);
    except:#used to make sure shutters re-open even in case of error or KeyboardInterrupt
        #reset to single shot on pulse picker
        pvslicerenable.put(0);time.sleep(0.5);
        pvslicerEC.put(182);time.sleep(0.5);
        pvslicerenable.put(1)
        #re-open shutters
        #print('Opening all shutters...')
        #toggle_TTL_shutter('openall',display=False);#open all the shutters
        resetMBC();
        YFEsetall(True,displayQ=False);
        
def psrecipes():
    allrec=glob.glob(psfilepath()+'recipes/*.p');
    oldrec=glob.glob(psfilepath()+'recipes/*_*.p')
    currec=[ext[60:-2] for ext in allrec if ext not in oldrec];currec.sort();
    return currec

def pspreshot():
    pshostcheck()
    pvMBCpower=EpicsSignal('MEC:64B:PWR:2:Outlet:8:SetControlAction')#WAS 'MEC:S60:PWR:01:Outlet:7:SetControlAction'#read AND write:1=ON,2=OFF
    pvMBCmode=EpicsSignal('MEC:LPL:MBC:01:RunningMode_RBV',write_pv='MEC:LPL:MBC:01:RunningMode')#AUTO=0,MAN=1
    pvMBCsetpt=EpicsSignal('MEC:LPL:MBC:01:AutoCalibration.VAL',write_pv='MEC:LPL:MBC:01:AutoCalibration') #QUAD=0,MIN=1,MAX=2
    pvMBCbias=EpicsSignal('MEC:LPL:MBC:01:BiasValue_RBV',write_pv='MEC:LPL:MBC:01:BiasValue')
    pvMBCfault=EpicsSignal('MEC:LPL:MBC:01:ErrorStatus',write_pv='MEC:LPL:MBC:01:ClearErrors')
    pvlampEC=EpicsSignal('EVR:MEC:USR01:TRIG6:EC_RBV',write_pv='EVR:MEC:USR01:TRIG6:TEC')#lamp event code; needs 182
    pvslicerEC=EpicsSignal('EVR:MEC:USR01:TRIG7:EC_RBV',write_pv='EVR:MEC:USR01:TRIG7:TEC')#slicer event code; needs 182
    pvlampenable=EpicsSignal('EVR:MEC:USR01:TRIG6:TCTL') #lamp enable;
    pvslicerenable=EpicsSignal('EVR:MEC:USR01:TRIG7:TCTL') #slicer enable;
    prechk=True
    if not YFEoncheck(display=False):
        print('WARNING: YFE seems to be off... Attempting to turn on YFE...')
        YFEon();
    if np.sum(YFEget(display=False)) < 550:
        if np.sum(YFEget(display=False)) < 20:
            print('WARNING: YFE seems to be turned down. Attempting to turn up YFE...')
            YFEsetall(True,displayQ=True)
    if isMBCsafe():
        print('Turning off MBC bias dither...',end='',flush=True);
        currbias=pvMBCbias.get();
        pvMBCmode.put(1);#set bias mode to MAN
        dotsleep(3);
        pvMBCbias.put(currbias+10);
        print('Testing bias responsivity...',end='',flush=True);
        for ii in range(2):
           time.sleep(1);print('..',end='',flush=True);
        if np.abs(pvMBCbias.get() - (currbias+10)) < 3:
            pvMBCbias.put(currbias);
            dotsleep(2);
        else:
            print('*')
            print('WARNING: MBC not responding!!')
            prechk = False
    else:
        print('MBC is not safe! Resetting the MBC...')
        resetMBC();##up above, check emission AND check currents???
        YFEsetall(True,displayQ=False)
        pspreshot()
        prechk = False
    if pvlampEC.get() != 182:
        pvlampenable.put(0)
        pvlampEC.put(182)
    if pvslicerEC.get() != 182:
        pvslicerenable.put(0)
        pvslicerEC.put(182)
    if pvlampenable.get() != 1:
        pvlampenable.put(1)
    if pvslicerenable.get() != 1:
        pvslicerenable.put(1)
    headenb = LasCoeff()
    wpreadpvlist=['MEC:NS1:MMS:02.RBV','MEC:NS1:MMS:01.RBV','MEC:LAS:MMN:30.RBV','MEC:LAS:MMN:29.RBV']
    wpwritepvlist=['MEC:NS1:MMS:02.VAL','MEC:NS1:MMS:01.VAL','MEC:LAS:MMN:30.VAL','MEC:LAS:MMN:29.VAL']
    wppvlist=[EpicsSignal(wpreadpvlist[ii],write_pv=wpwritepvlist[ii]) for ii in range(4)]
    headstr='';wpstr='';headlist=['AB','EF','GH','IJ'];
    for ii in range(4):
        if headenb[ii]:
            headstr+=headlist[ii]
        wpstr=wpstr+headlist[ii]+': '+str(round(wppvlist[ii].get(),3))+', '
    wpen=np.mean([headenb[ii]*np.cos((np.pi/180)*2*wppvlist[ii].get())**2 for ii in range(4)])
    print('The following heads are enabled: '+headstr)
    print('The waveplate settings are: '+wpstr[:-2])
    print('This is ~'+str(round(100*wpen,3))+'% of max energy.')
    print('Current pulse target is: '+str(pickle.load(open(psfilepath()+'Psns.p','rb')))+' ns, '+str(pickle.load(open(psfilepath()+'SSs.p','rb')))+' % of max power.')
    not_charging='';not_enabled='';yes_enabled='';headchanlist=['CD','A','B','E','F','G','H','I','J'];
    for ii in range(9):
        temppv1=EpicsSignal('MEC:PFN:CH'+str(ii)+':ENABLE_RBV');temppv2=EpicsSignal('MEC:PFN:CH'+str(ii)+':CHARGE_STATE');
        if temppv1.get() == 0:
            not_enabled+=headchanlist[ii]
        if temppv1.get() == 1:
            yes_enabled+=headchanlist[ii]
        if (temppv1.get() == 1) and (temppv2.get() == 0):
            not_charging+=headchanlist[ii]
    if len(not_charging)>0:
        print('** WARNING: The following heads are enabled but NOT charging: '+not_charging) 

    toggle_TTL_shutter('open'+yes_enabled+'wwxx',display=False);time.sleep(1);#make sure all shutters are open...
    toggle_TTL_shutter('close'+not_enabled,display=False)#close shutters that aren't enabled
    return prechk
    #waveform pre-check? verify shutters are open?

def pspostshot(save_flag_q=True,RunNumQQ=9000):
    pvlampenable=EpicsSignal('EVR:MEC:USR01:TRIG6:TCTL') #lamp enable;
    pvslicerenable=EpicsSignal('EVR:MEC:USR01:TRIG7:TCTL') #slicer enable;
    pvlampenable.put(0)
    pvslicerenable.put(0)
    pshostcheck()
    psacqx(save_flag=save_flag_q,RunNumQQ=RunNumQQ)#took out _noLecroyA
    #psefc();
    toggle_TTL_shutter('openall',display=False);#make sure all shutters are open again...
    print('Resetting bias tracking...')
    resetMBC();
    YFEsetall(True);
    E_synth_refresh();


def save_scope_to_eLog(chan_to_eLog=2):
    ExpName=get_curr_exp()
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
        mecel = elog.ELog({'experiment':ExpName},user='mecopr',pw=pickle.load(open(psfilepath()+'elogauth.p','rb')))
        #print('Connected to eLog of current MEC experiment, which is: '+ExpName)
    except:
        print('Failed to connect to eLog!')
    try:
        S1=LXOpen('1');time.sleep(.15);
        timestamp=datetime.now().strftime('%Y%m%d_%H%M%S')
        chdata=rchallxy(S1);
        LXClose(S1);time.sleep(.15);
    except:
        print('Failed to read out data!')
        LXClose(S1);time.sleep(.15);
    for ii in range(4):
        np.savetxt(fpQ+'lecroy1_ch'+str(ii+1)+'_'+timestamp+'.dat',tuple(chdata[ii]))
    eplxysav(chdata[chan_to_eLog-1][0],chdata[chan_to_eLog-1][1],fpQ+'lecroy1_ch'+str(chan_to_eLog)+'_'+timestamp+'.png',abs_path=True)
    fullmsg=str('Scope trace data for all 4 channels saved to '+fpQ+' with time stamp '+timestamp+'. Attached are the data and plot files for channel '+str(chan_to_eLog)+'.')
    #eplxy(chdata[chan_to_eLog-1][0],chdata[chan_to_eLog-1][1])
    try:
        mecel.post(fullmsg,attachments=[fpQ+'lecroy1_ch'+str(chan_to_eLog)+'_'+timestamp+'.dat', fpQ+'lecroy1_ch'+str(chan_to_eLog)+'_'+timestamp+'.png'], tags=['scope_trace'])
        print('Auto-saved to eLog.') 
    except:
        print('Failed to auto-save to eLog!')
    
def SHG_opt(armsQ='ABEFGHIJ'):#check for trace height;#All shutters must start in the open state... 
    print('Running this routine requires ALL TTL shutters to begin in the open state! The YFE must be on with the bias dither initially enabled!')
    if np.sum(TTL_shutter_status(display=False)[-1]) > 0:
        print('Warning! The shutters don\'t all appear to be open! ',end='',flush=True);TTL_shutter_status(display=True);
    else:
        print('(Shutters seem OK...)')
    if not YFEoncheck(display=False):
        print('Warning! The YFE doesn\'t appear to be on! ',end='',flush=True);YFEoncheck(display=True);
    else:
        print('(YFE emission seems OK...)')
    if np.sum(YFEget(display=False)) < 550:
        print('Warning! The YFE doesn\'t appear to be turned up! ');YFEget(display=True);
    else:
        print('(YFE current seems OK...)')
    if MBCmodecheck() != 0:
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
    HWPon('all',set_T=1)
    pvslicerEC=EpicsSignal('EVR:MEC:USR01:TRIG7:EC_RBV',write_pv='EVR:MEC:USR01:TRIG7:TEC')#slicer event code; needs 43
    pvslicerenable=EpicsSignal('EVR:MEC:USR01:TRIG7:TCTL') #slicer enable; 0=off,1=on
    pvMBCmode=EpicsSignal('MEC:LPL:MBC:01:RunningMode_RBV',write_pv='MEC:LPL:MBC:01:RunningMode')#AUTO=0,MAN=1
    armlist=['AB','EF','GH','IJ']
    #YFEoff();YFEon();
    pvMBCmode.put(1)#set MAN mode on MBC
    if np.sum(YFEget(display=False)) < 100:
        print('Check YFE before optimizing!')
    optwvfm=pickle.load(open(psfilepath()+'opttrace.p','rb'));
    try:
        S=HOpen();time.sleep(.15);oldwvfm=ReadPulseHeights(S,0);time.sleep(.15);
        WritePulseHeights(S,0,optwvfm);time.sleep(.15);HClose(S);time.sleep(.15);
    except:
        HClose(S);
    motprefix='MEC:LAS:MMN:';motnamelist=['22','24','17','18'];#VAL/RBV
    SHGpvlist=[EpicsSignal(motprefix+motname+'.RBV',write_pv=motprefix+motname+'.VAL') for motname in motnamelist];
    print('Closing all shutters...')
    toggle_TTL_shutter('closeall',display=False);#close all the shutters
    time.sleep(4)
    pvslicerEC.put(43);pvslicerenable.put(1);#enable these...
    try:
        SLA=LXOpen('A');time.sleep(.15);
        tempchk1=rch(1,SLA);time.sleep(.15);tempchk2=rch(1,SLA);
        if np.sum(np.abs(tempchk1-tempchk2))<1e-6:
            print('Warning: scope trace doesn\'t appear to be updating, please check scope! Abort? [enter y/n]')
            checkprompt=input();
            if checkprompt.lower() != 'y':
                print('Try again later then!');
                return
            else:
                print('OK, I hope you know what you\'re doing!')
        LXClose(SLA);time.sleep(.15);
    except:
        print('Scope error, check scope status! Aborting...')
        LXClose(SLA);time.sleep(.15);
        return
    startposlist=[SHGrbv.get() for SHGrbv in SHGpvlist];
    newposlist=startposlist[:]
    #alphEFGH=0.5;
    alphlist=[1,0.5,0.5,1];
    for ii in range(4):
        if armlist[ii] in armsQ:#only prep the stage if it's going to be used
            #if armlist[ii] in ['EF','GH']:
            #    alph=alphEFGH
            #else:
            #    alph=1
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
        SLA=LXOpen('A');time.sleep(.15);#changed to LecroyA since repair
        for ii in range(4):
            if armlist[ii] in armsQ:
                #if armlist[ii] in ['EF','GH']:
                #    alph=alphEFGH;
                #else:
                #    alph=1;
                print('Begin optimizing '+armlist[ii]+'... ',end='',flush=True);
                shgarmdatax,shgarmdatay=[],[]
                toggle_TTL_shutter('open'+armlist[ii],display=False);currentshutter=ii;time.sleep(4);print('Shutter opened!');#open one shutter
                for jj in range(11):
                    print('.',end='',flush=True)
                    SHGpvlist[ii].put(startposlist[ii]+alphlist[ii]*(-.1+.02*(jj)));time.sleep(2.5);#step to new position
                    curr_x=SHGpvlist[ii].get();curr_y=np.max(rch(3,SLA));time.sleep(.15);#in testing, max is more stable than sum
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
                toggle_TTL_shutter('close'+armlist[ii],display=False);currentshutter=0;#close that shutter;
                xpq=np.arange(startposlist[ii]+alphlist[ii]*(-.1+.02*(-1)),startposlist[ii]+alphlist[ii]*(-.1+.02*(11)),.0001);
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
        LXClose(SLA);time.sleep(.15);#changed to LeCroyA
    except:
        print('Failed! Restoring original values and attempting to re-open most-recent shutter... you should verify!')
        LXClose(SLA);time.sleep(.15);#changed to LeCroyA
        if currentshutter > 0:
            toggle_TTL_shutter('open'+armlist[currentshutter],display=False);
        for ii in range(4):
            SHGpvlist[ii].put(startposlist[ii]);newposlist[ii]=startposlist[ii];
    time.sleep(2);#need time so that last shutter trigger ends before trying to open IJ
    try:
        S=HOpen();time.sleep(.15);WritePulseHeights(S,0,oldwvfm);time.sleep(.15);HClose(S);time.sleep(.15);
    except:
        HClose(S);
        print('Error! Check waveform!')
    pvslicerenable.put(0);#disable PC before re-opening shutters
    datestamp=int(datetime.now().strftime('%Y%m%d%H%M%S'))
    SHGlog=pickle.load(open(psfilepath()+'SHG_opt_log.p','rb'))
    SHGlog.append([datestamp,[newposlist[ii] for ii in range(4)]])
    pickledump2(SHGlog,psfilepath()+'SHG_opt_log.p')
    toggle_TTL_shutter('openall',display=False);#open all the shutters
    resetMBC();YFEsetall(True);#reset bias...

def HWP_opt(armsQ='ABEFGHIJ'):#check for trace height;#All shutters must start in the open state... 
    print('Running this routine requires ALL TTL shutters to begin in the open state! The YFE must be on with the bias dither initially enabled!')
    if np.sum(TTL_shutter_status(display=False)[-1]) > 0:
        print('Warning! The shutters don\'t all appear to be open! ',end='',flush=True);TTL_shutter_status(display=True);
    else:
        print('(Shutters seem OK...)')
    if not YFEoncheck(display=False):
        print('Warning! The YFE doesn\'t appear to be on! ',end='',flush=True);YFEoncheck(display=True);
    else:
        print('(YFE emission seems OK...)')
    if np.sum(YFEget(display=False)) < 550:
        print('Warning! The YFE doesn\'t appear to be turned up! ');YFEget(display=True);
    else:
        print('(YFE current seems OK...)')
    if MBCmodecheck() != 0:
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
    HWPon('all',set_T=1)
    pvslicerEC=EpicsSignal('EVR:MEC:USR01:TRIG7:EC_RBV',write_pv='EVR:MEC:USR01:TRIG7:TEC')#slicer event code; needs 43
    pvslicerenable=EpicsSignal('EVR:MEC:USR01:TRIG7:TCTL') #slicer enable; 0=off,1=on
    pvMBCmode=EpicsSignal('MEC:LPL:MBC:01:RunningMode_RBV',write_pv='MEC:LPL:MBC:01:RunningMode')#AUTO=0,MAN=1
    armlist=['AB','EF','GH','IJ']
    #YFEoff();YFEon();
    pvMBCmode.put(1)#set MAN mode on MBC
    if np.sum(YFEget(display=False)) < 100:
        print('Check YFE before optimizing!')
    optwvfm=pickle.load(open(psfilepath()+'opttrace.p','rb'));
    try:
        S=HOpen();time.sleep(.15);oldwvfm=ReadPulseHeights(S,0);time.sleep(.15);
        WritePulseHeights(S,0,optwvfm);time.sleep(.15);HClose(S);time.sleep(.15);
    except:
        HClose(S);
    motnamelist=['MEC:NS1:MMS:02','MEC:NS1:MMS:01','MEC:LAS:MMN:30','MEC:LAS:MMN:29'];#VAL/RBV
    HWPpvlist=[EpicsSignal(motname+'.RBV',write_pv=motname+'.VAL') for motname in motnamelist];
    print('Closing all shutters...')
    toggle_TTL_shutter('closeall',display=False);#close all the shutters
    time.sleep(4)
    pvslicerEC.put(43);pvslicerenable.put(1);#enable these...
    startposlist=[HWPrbv.get() for HWPrbv in HWPpvlist];
    newposlist=startposlist[:]
    for ii in range(4):
        if armlist[ii] in armsQ:#only prep the stage if it's going to be used
            HWPpvlist[ii].put(startposlist[ii]+(-20.0+4.0*0))
    currentshutter=0;#trying to re-open a shutter in case of failure...
    stepQ=1.0;rangeQ=20.0;
    try:
        SLA=LXOpen('A');time.sleep(.15);#changed to LecroyA since repair
        for ii in range(4):
            if armlist[ii] in armsQ:
                print('Begin optimizing '+armlist[ii]+'... ',end='',flush=True);
                hwparmdatax,hwparmdatay=[],[]
                toggle_TTL_shutter('open'+armlist[ii],display=False);currentshutter=ii;time.sleep(4);print('Shutter opened!');#open one shutter
                for jj in range(int(1+(2*rangeQ/stepQ))):
                    print('.',end='',flush=True)
                    HWPpvlist[ii].put(startposlist[ii]+(-rangeQ+stepQ*(jj)));time.sleep(4);#step to new position#was 2.5
                    curr_x=HWPpvlist[ii].get();curr_y=np.max(rch(3,SLA));time.sleep(.15);#in testing, max is more stable than sum
                    if curr_y > 0.005:#threshold so don't skew fit with noise; max is ~~10x this
                        hwparmdatax.append(curr_x);hwparmdatay.append(curr_y);#save x and y
                    print('.',end='',flush=True)
                print('*')
                qfit=np.polyfit(hwparmdatax,hwparmdatay,2);newpos=qfit[1]/(-2*qfit[0]);#find fit and new max
                if np.abs(startposlist[ii]-newpos)<0.85*rangeQ:
                    HWPpvlist[ii].put(newpos);newposlist[ii]=newpos;
                    print('HWP position on arm '+armlist[ii]+' changed from '+str(round(startposlist[ii],4))+' to '+str(round(newpos,4)))
                else:
                    print('Failed! New HWP position on arm '+armlist[ii]+' seems too far off... '+str(round(newpos,4))+' from '+str(round(startposlist[ii],4))+'... Restoring...')
                    HWPpvlist[ii].put(startposlist[ii])
                toggle_TTL_shutter('close'+armlist[ii],display=False);currentshutter=0;#close that shutter;
                xpq=np.arange(startposlist[ii]+(-rangeQ+stepQ*(-1)),startposlist[ii]+(-rangeQ+stepQ*int(1+(2*rangeQ/stepQ))),.1);
                qfitp=np.poly1d(qfit);
                epllxy([[hwparmdatax,hwparmdatay],[xpq,qfitp(xpq)]],xlb=armlist[ii])
            else:
                print('Skipping '+armlist[ii]+'...')
                pass
        LXClose(SLA);time.sleep(.15);#changed to LeCroyA
    except:
        print('Failed! Restoring original values and attempting to re-open most-recent shutter... you should verify!')
        LXClose(SLA);time.sleep(.15);#changed to LeCroyA
        if currentshutter > 0:
            toggle_TTL_shutter('open'+armlist[currentshutter],display=False);
        for ii in range(4):
            HWPpvlist[ii].put(startposlist[ii]);newposlist[ii]=startposlist[ii];
    time.sleep(2);#need time so that last shutter trigger ends before trying to open IJ
    try:
        S=HOpen();time.sleep(.15);WritePulseHeights(S,0,oldwvfm);time.sleep(.15);HClose(S);time.sleep(.15);
    except:
        HClose(S);
        print('Error! Check waveform!')
    pvslicerenable.put(0);#disable PC before re-opening shutters
    datestamp=int(datetime.now().strftime('%Y%m%d%H%M%S'))
    HWPlog=pickle.load(open(psfilepath()+'HWP_opt_log.p','rb'))
    HWPlog.append([datestamp,[newposlist[ii] for ii in range(4)]])
    pickledump2(HWPlog,psfilepath()+'HWP_opt_log.p')
    toggle_TTL_shutter('openall',display=False);#open all the shutters
    resetMBC();YFEsetall(True);#reset bias...
    motnamelist=['MEC:NS1:MMS:02','MEC:NS1:MMS:01','MEC:LAS:MMN:30','MEC:LAS:MMN:29']#.OFF
    #add adjustment to offset automatically....
    for ii in range(4):
        temppv=EpicsSignal(motnamelist[ii]+'.OFF');tempval=temppv.get();
        temppv.put(tempval-newposlist[ii]);

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

def smooth_wvfm(wvfm_in):
    wvfm_out=wvfm_in[:]
    for ii in range(len(wvfm_in)-2):
        wvfm_out[ii+1]=.25*wvfm_in[ii]+.5*wvfm_in[ii+1]+.25*wvfm_in[ii+2]
    return wvfm_out

def E_coeff_refresh():
    pvlist=['MEC:LAS:FLOAT:'+str(ii) for ii in range(31,41)];
    inddesclist=['YFE','CD1w','AB1w','EF1w','GH1w','IJ1w','AB2w','EF2w','GH2w','IJ2w']
    desclist=['E_coeff_'+inddesc for inddesc in inddesclist]
    valulist=[.3578,0.5971,224.0,177.5,307.4*0.849,113.2,111.0*1.17,187.9*0.860,182.1*0.897,123.5*1.25]
    for jj in range(len(pvlist)):
        temppv1=EpicsSignal(str(pvlist[jj]+'.DESC'));temppv2=EpicsSignal(pvlist[jj]);
        temppv1.put(desclist[jj]);temppv2.put(valulist[jj]);

def E_synth_refresh():
    pvlist=['MEC:LAS:FLOAT:'+str(ii).zfill(2) for ii in range(1,11)];
    inddesclist=['YFE','CD1w','AB1w','EF1w','GH1w','IJ1w','AB2w','EF2w','GH2w','IJ2w']
    desclist=['E_synth_'+inddesc for inddesc in inddesclist]

    pveyfe=EpicsSignal('MEC:LAS:LEM:03:A:CUR_DISP');eyfe=pveyfe.get();
    pve1in=EpicsSignal('MEC:LAS:LEM:03:B:CUR_DISP');e1in=pve1in.get();
    pveab1w=EpicsSignal('MEC:LAS:GENTEC:02:CH1:MEAS');eab1w=pveab1w.get();
    pveef1w=EpicsSignal('MEC:LAS:GENTEC:02:CH2:MEAS');eef1w=pveef1w.get();
    pvegh1w=EpicsSignal('MEC:LAS:GENTEC:01:CH1:MEAS');egh1w=pvegh1w.get();
    pveij1w=EpicsSignal('MEC:LAS:GENTEC:01:CH2:MEAS');eij1w=pveij1w.get();
    pveab2w=EpicsSignal('MEC:LAS:GENTEC:03:CH1:MEAS');eab2w=pveab2w.get();
    pveef2w=EpicsSignal('MEC:LAS:GENTEC:03:CH2:MEAS');eef2w=pveef2w.get();
    pvegh2w=EpicsSignal('MEC:LAS:GENTEC:04:CH1:MEAS');egh2w=pvegh2w.get();
    pveij2w=EpicsSignal('MEC:LAS:GENTEC:04:CH2:MEAS');eij2w=pveij2w.get();
    energyarr=np.array([eyfe,e1in,eab1w,eef1w,egh1w,eij1w,eab2w,eef2w,egh2w,eij2w])
    
    coefflist=[]
    for ii in range(31,41):
        temppv=EpicsSignal('MEC:LAS:FLOAT:'+str(ii));
        coefflist.append(temppv.get());
    
    valulist=energyarr*np.array(coefflist)
    for jj in range(len(pvlist)):
        temppv1=EpicsSignal(str(pvlist[jj]+'.DESC'));temppv2=EpicsSignal(pvlist[jj]);
        temppv1.put(desclist[jj]);temppv2.put(valulist[jj]);

def TTL_shutter_refresh():
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

def DG_refresh():
    pvlist=[['MEC:LAS:DDG:0'+str(numii)+':'+chii+'DelaySI.DESC' for chii in ['a','c','e','g']] for numii in [1,2,6,8]];
    desclist=[['A:PS LATE','C:INH PS EARLY','E:unused','G: unused'],['A:PS EARLY','C:unused','E:EvoHE1','G:EvoHE2'],['A:GaiaQSW','C:GaiaLamp','E:INH UNI','G:GigE TRIG IN'],['A:BIG UNI','C:small UNI','E:INH GigE','G:unused']]
    for eachboxii in range(len(pvlist)):
        for eachentryii in range(len(eachboxii)):
            temppv=EpicsSignal(pvlist[eachboxii][eachentryii])
            temppv.put(desclist[eachboxii][eachentryii])

def SPLEG():
    reen=EpicsSignal('MEC:LAS:GENTEC:07:CH1:MEAS').get() 
    toen=EpicsSignal('MEC:LAS:GENTEC:07:CH2:MEAS').get() 
    m1en=EpicsSignal('MEC:LAS:GENTEC:06:CH1:MEAS').get() 
    m2en=EpicsSignal('MEC:LAS:GENTEC:06:CH2:MEAS').get() 
    regenen=1.64e5*reen + 1.03156061e-01 
    topasen=3.48e7*toen - 1.63e1 
    mpa1en=1.81e5*m1en - 0.301 #1.76e5, -2.22e0
    mpa2en=1.05e5*m2en - 1.39e-1 
    return np.round([regenen,topasen,mpa1en,mpa2en],2)

def vinko_temp():
    dat = np.array(np.genfromtxt('/cds/home/e/efcunn/Documents/lw36_vinko.txt'))
    segs2=np.append((dat.T[1][1:]/0.178383/.01)[0::10],100);
    segs3=np.array([0.1,4,8,12,18.1,25,40,65,78,87,100]);
    rPsns=9*[1.75]+[1.75]
    rSSs=[[segs2[ii],segs2[ii+1]] for ii in range(len(segs2)-1)]
    rSSs3=[[segs3[ii],segs3[ii+1]] for ii in range(len(segs3)-1)]
    plt.plot(0.178383*PulseGoal(9*[1.75]+[1.75],[[segs2[ii],segs2[ii+1]] for ii in range(len(segs2)-1)]))
    pwtin=0*ExponentialWave2(500,.03,1000,.03,0,5002)
    r10SSs=[[.08*(((ii+1)/11.0)**1.5)*rSSs3[ii][0]/100,.08*(((ii+2)/11.0)**1.5)*rSSs3[ii][1]/100] for ii in range(len(rSSs3))]
    r10Psns=[[500+(100*1.5*ii),500+(100*1.5*(ii+1))] for ii in range(len(rPsns))]
    pwtinstr=''
    for ii in range(len(r10Psns)):
        pwtin+=ExponentialWave2(r10Psns[ii][0],r10SSs[ii][0],r10Psns[ii][1],r10SSs[ii][1],0,5002)
        pwtinstr+='ExponentialWave2('+str(int(r10Psns[ii][0]))+','+str(np.round(r10SSs[ii][0],5))+','+str(int(r10Psns[ii][1]))+','+str(np.round(r10SSs[ii][1],5))+',0,5002)'
        if ii < len(r10Psns)-1:
            pwtinstr+='+'
    pickledump2(rPsns,'Psns.p','wb')
    pickledump2(rSSs,'SSs.p','wb')

def pickledump2(objQ,fullFileNameQ):
    pickle.dump(objQ,open(fullFileNameQ,'wb'));
    os.chmod(fullFileNameQ,stat.S_IRUSR|stat.S_IWUSR|stat.S_IRGRP|stat.S_IWGRP|stat.S_IROTH|stat.S_IWOTH);#
    #os.chmod(fullFileNameQ,stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO);#
    return

def plzchkpv(inpvnam):
    try:
        qp=EpicsSignal(inpvnam);
        currval=qp.get()
        msgout = str(currval)+' vs oldval'
    except TimeoutError as err:
        msgout = 'Timeout!'
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
    msgout=['{:<34}'.format(pvnam[0]),'{:<26}'.format(iocnam[0]),'{:<19}'.format(hostnam[0]),'{:<25}'.format(locstat[0]),cstr('{:<6}'.format(str('true ' if netstat==0 else 'false')),str('blink,r' if netstat!=0 else ''))]
    return msgout
    
def pv_checker(las='lpl'):
    if las.lower() == 'lpl':
        qqpl=np.genfromtxt(psfilepath()+'_ps_pvlist_lpl.txt',delimiter='\n',dtype=str);pvmsg=[];
    elif las.lower() == 'spl':
        qqpl=np.genfromtxt(psfilepath()+'_ps_pvlist_spl.txt',delimiter='\n',dtype=str);pvmsg=[];
    else:
        print('No such laser!');return
    for eapv in qqpl:
        pvmsg.append(plzchkpv(eapv))
    with multiprocessing.Pool() as pool:
        srvmsg=pool.map(plzchksrv,qqpl)
    print(''.join(['{:<34}'.format('PV name'),'{:<26}'.format('IOC name'),'{:<19}'.format('Host name'),'{:<25}'.format('Host location'),'{:<6}'.format('Ping?'),'{:<15}'.format('PV value?')]))
    for ii in range(len(pvmsg)):
        print(''.join(srvmsg[ii])+pvmsg[ii])
    return #datq

def plzchkcmp(cmpnam):
    netstat=os.system('ping -c 1 -w2 '+cmpnam[0]+' > /dev/null 2>&1')
    msgout=['{:<15}'.format(cmpnam[2])+'{:<28}'.format(cmpnam[0])+cstr('{:<6}'.format(str('true ' if netstat==0 else 'false')),str('blink,r' if netstat!=0 else ''))]
    return msgout

def plzchkcmp2(cmpnam):
    netstat=os.system('ping -c 1 -w2 '+cmpnam[0]+' > /dev/null 2>&1')
    msgout=['{:<15}'.format(cmpnam[1])+'{:<28}'.format(cmpnam[0])+cstr('{:<6}'.format(str('true ' if netstat==0 else 'false')),str('blink,r' if netstat!=0 else ''))]
    return msgout

def MECcompylist():
    qqip=['172.21.46.147','172.21.46.148','172.21.46.146','172.21.46.60','172.21.46.128','172.21.46.100', '172.21.46.120','172.21.46.159', '172.21.46.197','172.21.46.142','172.21.46.70','172.21.46.78', '172.21.46.71','172.21.46.88','172.21.46.198','172.21.46.213','172.21.46.215','172.21.46.136', '172.21.46.218','172.21.46.219','172.21.46.182','172.21.46.144'];
    qqn=['evo1','evo2','gaia','lecroy1','lecroy2','lecroya','lecroyb','PIMikroMove','spider','spectrometer', 'tundra','topas','visar1','visar2','vitara','rga','emp','phasicslaptop','phasics1','phasics2','dacage','legend']
    nmlist=['mec-las-laptop06','mec-las-laptop07','mec-las-laptop05','scope-ics-mectc1-1','scope-ics-meclas-lecroy01','scope-ics-meclas-lecroy-a','scope-ics-meclas-lecroy-b','mec-las-laptop09','mec-las-laptop11','mec-las-laptop01','win-ics-mec-tundra','mec-las-laptop12','win-ics-mec-visar1','win-ics-mec-visar2','mec-las-vitara','mec-rga-laptop','scope-ics-mec-tektronix','mec-phasics-laptop01','win-ics-mec-phasics01','win-ics-mec-phasics02','mec-visar-cage','mec-las-laptop03']
    return list(zip(nmlist,qqip,qqn))

def MECcompylist2():
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

def cmp_checker():
    cmpmsg=[];
    #qqpl=[eacmp[0] for eacmp in MECcompylist()]
    qqpl=MECcompylist2()
    with multiprocessing.Pool() as pool:
        cmpmsg=pool.map(plzchkcmp,qqpl)
    print('{:<15}'.format('Computer name')+'{:<28}'.format('IP shorthand')+'{:<6}'.format('Ping?'))
    for ii in range(len(cmpmsg)):
        print(''.join(cmpmsg[ii]))
    return

def cmp2_checker():
    cmpmsg=[];
    #qqpl=[eacmp[0] for eacmp in MECcompylist()]
    qqpl=MECcompylist2()
    with multiprocessing.Pool() as pool:
        cmpmsg=pool.map(plzchkcmp2,qqpl)
    print('{:<15}'.format('Computer name')+'{:<28}'.format('IP shorthand')+'{:<6}'.format('Ping?'))
    for ii in range(len(cmpmsg)):
        print(''.join(cmpmsg[ii]))
    return

def tc():
    colors=['ENDC','BLINK','K','R','G','Y','B','M','C','W','BK','BR','BG','BY','BB','BM','BC','BW','BRK','BRR','BRG','BRY','BRB','BRM','BRC','BRW','BBRK','BBRR','BBRG','BBRY','BBRB','BBRM','BBRC','BBRW']#B- for background, BR+ for 
    colorcodes=['\033['+str(ii)+'m' for ii in [0,5,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47,90,91,92,93,94,95,96,97,100,101,102,103,104,105,106,107]];
    return dict(zip(colors,colorcodes))
    
def cprint(strQ,paramsQ):#delim with commas
    prargs=''
    if len(paramsQ) == 0:
        paramsQ='ENDC'
    for eaarg in paramsQ.split(','):
        prargs+=tc()[eaarg.upper()]
    print(f"{prargs}"+strQ+f"{tc()['ENDC']}")
    return

def cstr(strQ,paramsQ):
    prargs=''
    if len(paramsQ) == 0:
        paramsQ='ENDC'
    for eaarg in paramsQ.split(','):
        prargs+=tc()[eaarg.upper()]
    return f"{prargs}"+strQ+f"{tc()['ENDC']}"
    
def keybd():
    return dict(zip(['key_Enter','key_Esc','key_Up','key_Dn','key_Rt','key_Lt'],[13,27,'\033[A','\033[B','\033[C','\033[D']))

def getch():
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
    
def pldemo():
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
        qq=getch()
        if qq==keybd()['key_Dn']:
            ypos-=ydel
        elif qq==keybd()['key_Up']:
            ypos+=ydel
        elif qq==keybd()['key_Rt']:
            xpos+=xdel
        elif qq==keybd()['key_Lt']:
            xpos-=xdel
        elif qq=='w':
            ydel=ydel*2
        elif qq=='s':
            ydel=ydel/2
        elif qq=='a':
            xdel=xdel/2
        elif qq=='d':
            xdel=xdel*2
        elif qq==keybd()['key_Esc']:
            cont=False
        else:
            pass
        Z=[[np.exp(-((ii-xpos)/10)**2-((jj-ypos)/7.5)**2) for ii in range(100)] for jj in range(100)]
        print('['+str(xpos)+','+str(ypos)+']')
    plt.ioff()
    
def HWPclear_start():
    pvabef=['MEC:NS1:MMS:0'+str(num)+':SEQ_SELN' for num in [2,1]]#clear start for mforce chassis
    try:
        for eapv in pvabef:
            temppv=EpicsSignal(eapv);
            temppv.put(1)
    except:
        print('Failed!')
        
def NewportRestart():
    #XPS3 MEC:64A:PWR:2:Outlet:5:SetControlAction
    #XPS4 MEC:64B:PWR:1:Outlet:1:SetControlAction
    xps3and4=['1724','2532']
    ipvlist=['MEC:LAS:MMN_'+unitno+'.IALL' for unitno in xps3and4]
    rpvlist=['MEC:LAS:MMN_'+unitno+'.RALL' for unitno in xps3and4]
    try:
        for eapv in ipvlist:
            temppv=EpicsSignal(eapv);
            temppv.put(1)
        dotsleep(10)
        for eapv in rpvlist:
            temppv=EpicsSignal(eapv);
            temppv.put(1)
    except:
        print('Failed!')
#
def dotsleep(tSEC):
    for ii in range(tSEC):
        print('.',end='',flush=True);time.sleep(1);
    print('*')
    return
    
def ssleep():
    time.sleep(0.15);
    return

def reloadchk():
    print('Reload check: 20211025')


