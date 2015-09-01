
print("import common lib");
import glob;
import math;
import copy;
import enum;
print("import numpy");
import numpy as np;
print("import scipy");
from scipy import optimize;
print("import matplotlib");
import matplotlib.pyplot as plt;

# data type of struct
class DataType(enum.Enum):
    timeVoltage = 0;
    freqVoltage = 1;
    freqPower = 2;

# rawdata file type
class FileType(enum.Enum):
    TimeVoltage = 0;
    FreqVoltage = 1;
    FreqPower = 2;
    ROOscilloscope = 3;
    TektronixOscilloscope = 4;
    AnritsuSpectrumAnalyzer = 5;

# fitting function
class FittingType(enum.Enum):
    Lorentzian = 0;
    Gaussian = 1;
    ExpDecay = 2;

# data container
class DataSet:
    def __init__(self,_dataType,_startValue,_delta,_originalFileName,_array = np.empty(0),_impedance = 50):
        self.startValue = _startValue;
        self.delta = _delta;
        self.array= _array;
        self.impedance = _impedance;
        self.dataType = _dataType;

        self.originalFileName = _originalFileName;
        self.lastFileName = "";

        self.fitFlag = False;
        self.fitFunc = None;
        self.fitParam = None;
        self.fitError = None;
        self.fitArray = None;
        return;
    
    def getXArray(self):
        return np.arange(0,self.array.size)*self.delta+self.startValue;
    def getYArray(self):
        return self.array;
    def getYFitArray(self):
        return self.fitArray;

    def getStartValue(self):
        return self.startValue;
    def getEndValue(self):
        return self.startValue + self.array.size*self.delta;

    def getXString(self):
        if(self.dataType == DataType.timeVoltage):
            return "Time";
        else:
            return "Freq";
    def getYString(self):
        if(self.dataType == DataType.freqPower):
            return "Power Spectral Density";
        else:
            return "Voltage";
    def getXUnit(self):
        if(self.dataType == DataType.timeVoltage):
            return "s";
        else:
            return "Hz";
    def getYUnit(self):
        if(self.dataType == DataType.freqPower):
            return "W/Hz";
        else:
            return "V";

    def scale(self,coef):
        self.array*=coef;
    def shift(self,offset):
        self.array+=offset;

    def save(self,fileName):
        print("save imported "+self.originalFileName+" to "+fileName);

        # exit if try to overwrite original file
        if(fileName == self.originalFileName):
            print("cannot overwrite original file");
            return;

        self.lastFileName = fileName;
        fileHandle = open(fileName,"w");
        xvals = self.getXArray();
        yvals = self.getYArray();
        # write real value array
        if(self.dataType == DataType.timeVoltage or self.dataType == DataType.freqPower):
            for ind in range(self.array.size):
                fileHandle.write(str(xvals[ind])+" "+str(yvals[ind])+"\n");
        # write imaginary value array
        if(self.dataType == DataType.freqVoltage):
            for ind in range(len(data.array)):
                fileHandle.write(str(xvals[ind])+" "+str(yvals[ind].real)+" "+str(yvals[ind].imag)+"\n");
        fileHandle.close();

    def plot(self,fileName="",logFlag=False,show=True):
        plt.figure(figsize=(8,5));
        if(logFlag):
            plt.yscale("log");
        colors = ["b","g","r","c","m","y","k","w"];

        xdata = self.getXArray();
        if(self.dataType == DataType.freqVoltage):
            lineType = "b-";
            ydata = np.real(self.getYArray());
            plt.plot(xdata,ydata,lineType,label=self.originalFileName+" real");
            lineType = "r-";
            ydata = np.imag(self.getYArray());
            plt.plot(xdata,ydata,lineType,label=self.originalFileName+" imag");
        else:
            lineType = "b-";
            ydata = self.getYArray();
            plt.plot(xdata,ydata,lineType,label=self.originalFileName);

        if(self.fitFlag):
            yfdata = self.getYFitArray();
            lineType = "r-";
            plt.plot(xdata,yfdata,lineType,label=self.originalFileName+" fit");

        plt.xlabel(self.getXString());
        plt.ylabel(self.getYString());
        plt.legend(loc='best',fancybox=True, shadow=True);
        plt.grid(True);

        if(fileName != ""):
            plt.savefig(fileName,pid=100);
        if(show):
            plt.show();
        return;

    # fitting
    def fitting(self,fittingType,logResidual=False):
        print("fitting : "+self.originalFileName);
        xdata = self.getXArray();
        ydata = self.getYArray();

        if(not hasattr(self,"peakIndex")):
            extractPeak(self);

        peakPosition = self.startValue + self.peakIndex*self.delta; 
        peakAmplitude = self.array[self.peakIndex]-self.array[0];
        peakWidth = self.peakWidth*self.delta;
        peakOffset = self.array[0];
        param0 = [peakPosition,peakAmplitude,peakOffset,peakWidth];

        if(fittingType == FittingType.Lorentzian):
            func = lorentzian;
        if(fittingType == FittingType.Gaussian):
            func = gaussian;
        if(fittingType == FittingType.ExpDecay):
            func = expdecay;
        
        result = optimize.leastsq(getResidual,param0,args=(xdata,ydata,func,logResidual),full_output=True);
        x0 = result[0][0];
        a = result[0][1];
        o = result[0][2];
        w = result[0][3];
        ex0 = math.sqrt(np.abs(result[1][0][0]));
        ea = math.sqrt(np.abs(result[1][1][1]));
        eo = math.sqrt(np.abs(result[1][2][2]));
        ew = math.sqrt(np.abs(result[1][3][3]));
        pfx = self.getXUnit();
        pfy = self.getYUnit();

        print("reult  ");
        print(" * PeakPosition      : "+adjustUnit(x0,pfx)[2]+"   +-  "+str(myRound(ex0/x0*100,5))+" %");
        print(" * Amplitude         : "+adjustUnit(a,pfy)[2]+"  +-  "+str(myRound(ea/a*100,5))+" %");
        print(" * FullLineWidth     : "+adjustUnit(w*2,pfx)[2]+"  +-  "+str(myRound(ew/w*100,5))+" %");
        print(" * Offset            : "+adjustUnit(o,pfy)[2]+"  +-  "+str(myRound(eo/o*100,5))+" %");

        self.fitFlag=True;
        self.fitFunc = func;
        self.fitParam = [x0,a,o,w];
        self.fitError = [ex0,ea,eo,ew];
        self.fitArray = func(self.fitParam,xdata);

    def fittingEcho(self,fittingType,logResidual=False,echoTail=50):
        print("fitting echo : "+self.originalFileName);
        xdata = self.getXArray();
        ydata = self.getYArray();

        if(not hasattr(self,"peakIndex")):
            extractPeak(self);

        peakWidth = self.peakWidth*self.delta;
        peakDecay = peakWidth/10;
        peakAmplitude = self.array[self.peakIndex]-self.array[0];
        peakOffset = self.array[0];
        peakPosition = 0; 
        #param0 = [peakPositionFromLeft,peakAmplitude,peakOffset,peakWidth,peakDecay];
        param0 = [peakPosition,peakAmplitude,peakWidth,peakDecay,peakOffset];

        func = convolveWithLorentzian;

        result = optimize.leastsq(getEchoResidual,param0,args=(xdata,ydata,func,echoTail,self.peakIndex),full_output=True);
        
        x0 = self.startValue + self.peakIndex*self.delta + result[0][0];
        a = result[0][1];
        w = result[0][2];
        d = result[0][3];
        o = result[0][4];

        pfx = self.getXUnit();
        pfy = self.getYUnit();

        print("");
        if(result[1] is not None):
            ex0 = np.sqrt(result[1][0][0]);
            ea = np.sqrt(result[1][1][1]);
            ew = np.sqrt(result[1][2][2]);
            ed = np.sqrt(result[1][3][3]);
            eo = np.sqrt(result[1][4][4]);
            print("result  ");
            print(" * PeakPosition      : "+adjustUnit(x0,pfx)[2]+"   +-  "+str(myRound(ex0/x0*100,5))+" %");
            print(" * Amplitude         : "+adjustUnit(a,pfy)[2]+"  +-  "+str(myRound(ea/a*100,5))+" %");
            print(" * FullLineWidth     : "+adjustUnit(w*2,pfx)[2]+"  +-  "+str(myRound(ew/w*100,5))+" %");
            print(" * DecayWidth        : "+adjustUnit(d*2,pfx)[2]+"  +-  "+str(myRound(ed/d*100,5))+" %");
            print(" * Offset            : "+adjustUnit(o,pfy)[2]+"  +-  "+str(myRound(eo/o*100,5))+" %");
        else:
            ex0 = 0;
            ea = 0;
            ew = 0;
            ed = 0;
            eo = 0;
            print("result  ");
            print(" * Peak frequency    : "+adjustUnit(x0,pfx)[2]);
            print(" * Amplitude         : "+adjustUnit(a,pfy)[2]);
            print(" * FullLineWidth     : "+adjustUnit(w*2,pfx)[2]);
            print(" * DecayWidth        : "+adjustUnit(d*2,pfx)[2]);
            print(" * Offset            : "+adjustUnit(o,pfy)[2]);

#        plt.plot(xdata,ydata);
#        plt.plot(xdata,convolveWithLorentzian(param0,xdata,echoTail,self.peakIndex));
#        plt.plot(xdata,convolveWithLorentzian(result[0],xdata,echoTail,self.peakIndex));
#        plt.show();
        self.fitFlag=True;
        self.fitFunc = func;
        self.fitParam = [x0,a,o,w,d];
        self.fitError = [ex0,ea,eo,ew,ed];
        self.fitArray = convolveWithLorentzian(result[0],xdata,echoTail,self.peakIndex);

def getResidual(param,xdata,ydata,func,logResidual):
    ret = func(param,xdata);
    if(logResidual):
        return np.log(ydata)-np.log(ret);
    else:
        return ydata-ret;

def getEchoResidual(param,xdata,ydata,func,echoTail,peakIndex):
    ret = func(param,xdata,echoTail,peakIndex);
    print(".",end="",flush=True);
#    print(str(param)+" "+str(((ydata-ret)**2).sum()));
#    plt.plot(xdata,ret);
    return ydata-ret;

def convolveWithLorentzian(param,xdata,echoTail,peakIndex):
    x0 = param[0];
    amp = param[1];
    wid = np.abs(param[2]);
    dec = np.abs(param[3]);
    ofs = param[4];
    delta = xdata[1]-xdata[0];

    edlen = int(np.max([dec/delta,1])*echoTail);
    edx = np.arange(edlen)*delta;
    edy = expdecay([0,1,0,dec],edx);
    edy /= edy.sum()*delta;

#    lolen = np.max([int(np.max([wid/delta,1])*echoTail*2),xdata.size*2]);
    lolen = xdata.size*2;
    lox = np.arange(lolen)*delta;
    loy = lorentzian([lolen/2*delta-x0,1,0,wid],lox);
    cvy = np.convolve(edy,loy)*delta*amp+ofs;

    startIndex = np.max([lolen/2 - peakIndex,0]);
    cvy = cvy[startIndex:startIndex+xdata.size];
    return cvy;

# function for fitting
def gaussian(param,x):
    x0 = param[0];a = param[1];ofs = param[2];s = param[3];
    return a*np.exp(-(x-x0)**2/(2*s**2))+ofs;
def lorentzian(param,x):
    x0 = param[0];a = param[1];ofs = param[2];s = param[3];
    return a*(s**2)/((x-x0)**2 + s**2)+ofs;
def expdecay(param,x):
    x0 = param[0];a = param[1];ofs = param[2];s = param[3];
    return a*np.exp(-(x-x0)/s)+ofs;


def plotMultiple(dataArray,fileName="",logFlag=False,show=True):
    plt.figure(figsize=(8,5));
    if(logFlag):
        plt.yscale("log");
    colors = ["b","g","r","c","m","y","k","w"];
    cc = 0;
    for fileIndex in range(len(dataArray)):
        data = dataArray[fileIndex];
        xdata = data.getXArray();
        ydata = data.getYArray();
        lineType = colors[cc%len(colors)]+"o";
        cc+=1;
        plt.plot(xdata,ydata,lineType,label="data"+str(fileIndex));

        if(data.fitFlag):
            yfdata = data.getYFitArray();
            lineType = colors[cc%len(colors)]+"-";
            plt.plot(xdata,yfdata,lineType,label="data"+str(fileIndex)+" fit");
            cc+=1;

    plt.xlabel(dataArray[0].getXString());
    plt.ylabel(dataArray[0].getYString());
    plt.legend(loc='best',fancybox=True, shadow=True);
    plt.grid(True);

    if(fileName != ""):
        plt.savefig(fileName,pid=100);
    if(show):
        plt.show();
    return;

# time-voltage to freq-voltage
def fft(data):
    if(data.dataType != DataType.timeVoltage):
        print("cannot fft : " + str(data.dataType));
        return;
    sampleCount = len(data.array);
    print("start FFT : "+data.originalFileName + " with " + str(sampleCount) + " samples");
    fftlist = np.fft.fft(data.array);
    frqlist = np.fft.fftfreq(sampleCount,data.delta);
    print("end FFT");
    ndata = DataSet(DataType.freqVoltage,frqlist[0],abs(frqlist[1]-frqlist[0]),data.originalFileName,fftlist,data.impedance);
    return ndata;

# freq-voltage to time-voltage
def ifft(data):
    if(data.dataType != DataType.freqVoltage):
        print("cannot ifft : " + str(data.dataType));
        return;
    sampleCount = len(data.array);
    print("start IFFT : "+data.originalFileName + " with " + str(sampleCount) + " samples");
    ifftlist = np.fft.ifft(data.array);
    shiftlist = np.fft.fftfreq(sampleCount,data.delta);
    print("end IFFT");
    ndata = DataSet(DataType.timeVoltage,shiftlist[0],abs(shiftlist[1]-shiftlist[0]),data.originalFileName,ifftlist,data.impedance);
    return ndata;

# freq-voltage to freq-power
def calcPSD(data):
    if(data.dataType != DataType.freqVoltage):
        print("cannot calc power : " + str(data.dataType));
        return;
    print("convert : " + data.originalFileName + " to Freq-Power");
    rat = 1.0/data.impedance/len(data.array)/2/data.delta;
    ndata = DataSet(
        DataType.freqPower,
        data.startValue,
        data.delta,
        data.originalFileName,
        abs(data.array)**2*rat,
        data.impedance);
    return ndata;

# freq-power to integrated power
def calcPower(data):
    if(data.dataType != DataType.freqPower):
        print("cannot calc sum power");
        return;
    sum = 0.;
    for ind in range(0,len(data.array)):
        sum+= data.array[ind]*data.delta;
    rms = math.sqrt(sum*data.impedance);
    print("Total Power = "+str(sum)+" W  ->  V_RMS = "+str(rms)+"V");
    return rms;


# find peak value and width between (leftLimit,rightLimit)
#  and extract data from leftSigma to rightSigma
def extractPeak(data,leftSigma=-1,rightSigma=-1,leftLimit=-1,rightLimit=-1):
    print("extract data around peak : " + data.originalFileName);

    if(leftLimit < 0):
        leftLimit = 0;
    if(rightLimit < 0 ):
        rightLimit = len(data.array);

    maxIndex=max(range(leftLimit,rightLimit),key = lambda x:data.array[x]);
    maxValue=data.array[maxIndex];
    leftWidth = 0;
    rightWidth = 0;
    leftIndex = 0;
    rightIndex = len(data.array);

    for ind in list(reversed(range(0,maxIndex))):
        if(data.array[ind]<maxValue/2):
            leftWidth = (maxIndex-ind);
            if(leftSigma==-1):
                leftIndex = 0;
            else:
                leftIndex = max(0,maxIndex - leftWidth*leftSigma);
            break;

    for ind in range(maxIndex+1,len(data.array)):
        if(data.array[ind]<maxValue/2):
            rightWidth = (ind-maxIndex);
            if(rightSigma==-1):
                rightIndex = len(data.array);
            else:
                rightIndex = min(len(data.array),maxIndex + rightWidth*rightSigma+1);
            break;

    data.peakWidth = max(leftWidth,rightWidth);
    data.peakIndex = maxIndex-leftIndex;

    ndata = copy.copy(data);
    ndata.array = data.array[leftIndex:rightIndex];
    ndata.startValue += ndata.delta*leftIndex;
    return ndata;


# decrease number of elements in data
def combData(data,dataCount):
    print("comb : "+data.originalFileName);
    ndata = copy.copy(data);
    ndata.delta = (data.delta*(len(data.array)-1)/dataCount);
    ndata.array = [];
    ndataTime = 0;
    ndataInd = 0;
    for ind in range(0,dataCount):
        indPos = (ndata.delta*ind/data.delta);
        leftInd = math.floor(indPos);
        rightInd = math.ceil(indPos);
        if(leftInd == rightInd):
            rightInd +=1;
        leftValue = data.array[leftInd];
        rightValue = data.array[rightInd];
        ndata.array.append(leftValue + (rightValue-leftValue)*(indPos-leftInd));
    ndata.array = np.array(ndata.array);
    return ndata;

# extract data from leftValue to rightValue
def limitData(data,leftValue=None,rightValue=None):
    print("limiting : "+data.originalFileName);
    leftIndex = 0;
    rightIndex = len(data.array);
    if(leftValue != None):
        leftIndex = math.ceil((leftValue-data.startValue)/data.delta);
    if(rightValue != None):
        rightIndex = math.floor((rightValue-data.startValue)/data.delta);

    ndata = copy.copy(data);
    ndata.array = data.array[leftIndex:rightIndex];
    ndata.startValue += ndata.delta*leftIndex;
    return ndata;

# find preferable unit prefix and return [ratio, unitName, value+unitName  ]
def adjustUnit(refval,unitName):
    unitShift = 0;
    unitPrefixSmall = ["","m","u","n","p","f","a"];
    unitPrefixLarge = ["","k","M","G","T","P","E"];

    unitShift += math.floor(math.log10(abs(refval))/3);

    if(unitShift > 6 ):
        unitShift=6;
    if(unitShift < -6 ):
        unitShift=-6;

    ratio = math.pow(10,-3*unitShift);
    prefix = "";
    if(unitShift>=0):
        prefix = unitPrefixLarge[unitShift];
    else:
        prefix = unitPrefixSmall[-unitShift];
    return [ratio,prefix+unitName,str(myRound(refval*ratio,5))+" "+prefix+unitName];

# return rounded value with effective digit = prec
def myRound(val,prec):
    dec = math.floor(math.log10(abs(val)));
    return round(val,prec-dec-1);

# load datafile
def importDataFile(fileType , waveFileName):
    if(fileType == FileType.TimeVoltage):
        print("load TimeVoltage file : " +waveFileName);

        waveFileHandle = open(waveFileName,"r");
        lineCount=0;
        array = [];
        startValue = 0;
        delta = 0;
        for line in waveFileHandle:
            TimeVoltage = list(map(lambda x: float(x.strip()),line.split(" ")));
            array.append(TimeVoltage[1]);
            if(lineCount==0):
                startValue = TimeVoltage[0];
            if(lineCount==1):
                delta = TimeVoltage[0]-startValue;
            lineCount += 1;
        waveFileHandle.close();
        data = DataSet(DataType.timeVoltage,startValue,delta,waveFileName,np.array(array));
        return data;

    if(fileType == FileType.ROOscilloscope):
        print("load R&S Oscilloscope wave file : " +waveFileName);
        configFileName = waveFileName.replace(".Wfm.csv",".csv");
        config = {};
        configFileHandle = open(configFileName,"r");
        for line in configFileHandle:
            keyvalue = list(map(lambda x: x.strip(),line.split(":")));
            if(len(keyvalue)>=2):
                config[keyvalue[0]] = keyvalue[1];
        configFileHandle.close();
        delta = float(config["Resolution"]);
        startValue = float(config["HardwareXStart"]);

        array = [];
        waveFileHandle = open(waveFileName,"r");        
        for line in waveFileHandle:
            array.append(float(line.strip()));
        waveFileHandle.close();
        data = DataSet(DataType.timeVoltage,startValue,delta,waveFileName,np.array(array));
        return data;

    if(fileType == FileType.TektronixOscilloscope):
        print("load Tektronix Oscilloscope wave file : " +waveFileName);
        configLine = 18;
        lineCount = 0;
        deltaLogged = False;
        waveFileHandle = open(waveFileName,"r");
        lastTimeData = 1;
        lastValueData = 1;
        startValue = 0;
        delta = 1;
        array = [];
        for line in waveFileHandle:
            if(lineCount>=configLine):
                keyvalue = list(map(lambda x: x.strip(),line.split(",")));
                timeData = float(keyvalue[3]);
                valueData = float(keyvalue[4]);
                array.append(valueData);
                if(lineCount == configLine):
                    startValue = timeData;
                if(lineCount == configLine+1):
                    delta = timeData-startValue;
            lineCount += 1;
        waveFileHandle.close();
        data = DataSet(DataType.timeVoltage,startValue,delta,waveFileName,np.array(array));
        return data;

    if(fileType == FileType.AnritsuSpectrumAnalyzer):
        print("load Anritsu Spectrum Analyzer : " + waveFileName);
        waveFileHandle = open(waveFileName,"r");
        lineCount = 1;
        startFreq = 0;
        endFreq = 0;
        referenceLevel = 0;
        array = [];
        for line in waveFileHandle:
            if(lineCount == 2):
                values = list(map(lambda x:x.strip(),line.split(",")));
                startFreq = float(values[1]);
                endFreq = float(values[2]);
            if(lineCount == 11):
                values = list(map(lambda x:float(x.strip()),line.split(",")));
                referenceLevel = (1e-3) * math.pow(10.0,values[0]/10.0);
            if(lineCount>=14):
                value = float(line.strip());
                array.append(referenceLevel* math.pow(10,(value/10)));
            lineCount += 1;
        delta = (endFreq-startFreq)/(len(array)-1);
        array = list(map(lambda x:x/delta , array));
        waveFileHandle.close();

        data = DataSet(DataType.freqPower,startFreq,delta,waveFileName,np.array(array));
        return data;
