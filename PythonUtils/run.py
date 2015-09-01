
import Analyzer as a;
import numpy as np;
import matplotlib.pyplot as plt;
import glob;

delta = 1e-4;
sam = 20000;
dx = np.arange(sam)*delta;
dy = a.lorentzian([sam*delta/2-delta*4000 , 1.0 , 0.01 , 1000*delta ],dx) + np.random.rand(dx.size)/5;
d0 = a.DataSet(a.DataType.timeVoltage,-1.0,delta,"hoge.txt",dy);


#d0 = a.importDataFile(a.FileType.TimeVoltage,"10s0.5v_fb.txt");
#d0.scale(-1);
#d0 = a.limitData(d0,1,3);
#d0 = a.combData(d0,3000);

d0.plot();
d0.fittingEcho(None);
d0.plot();
#d0.fitting(a.FittingType.Gaussian);

#d0 = a.importDataFile(a.FileType.TimeVoltage,"10s0.5v_nfb.txt");
#d0.scale(-1);
#d0 = a.limitData(d0,2,4);
#d0 = a.combData(d0,3000);
#d0.plot();
#d0.fittingEcho(None);
