from ROOT import *
import numpy
from array import array
import sys
from searchClusters import searchClusters

   
def calibAmpl(fname,ofname):
    fileOut=TFile(ofname+".root","recreate");
    file=TFile(fname)

    #SET HERE SOME REASONABLE VALUES
    amplONE=10.
    sigmaONE=2.
    nSPECTRUM=5.
    
    hhL0=[]
    hhL1=[]
    gainL0=[1]*10
    gainL1=[1]*10

    ###L0### 
    #Get histograms
    for ii in range(0,10):
        h=file.Get("hBDXMiniVetoSIPMA_L0_C"+str(ii+1))
        hhL0.append(h)
        

    for ii in range(0,10):
        if (ii==2):
            continue
        
        h=hhL0[ii].Clone("hBDXMiniVetoSIPMA_L0_C"+str(ii+1)+"clone")
        min=0;
        max=0;
        sigma=0;
        ampl=0;
        if (ii<8):
            min=amplONE*1.5
            max=amplONE*nSPECTRUM
            sigma=sigmaONE
            ampl=amplONE
        else:
            min=(amplONE*1.5)/2
            max=(amplONE*nSPECTRUM)/2
            sigma=sigmaONE/2
            ampl=amplONE/2

        h.Smooth()
        h.GetXaxis().SetRangeUser(min,max);
        s=TSpectrum();
        n_found=s.Search(h,2,"",0.1);
        pm = h.GetListOfFunctions().FindObject("TPolyMarker")
        hhL0[ii].GetListOfFunctions().Add(pm)

        xP=s.GetPositionX();
        yP=s.GetPositionY();
        xxP=[]
        xxP=searchClusters(n_found,xP,yP,2*sigma)
        n_found=len(xxP)

      

        
        
        xPeak=array('d')
        xPeakErr=array('d')
        nPeak=array('d')
        nPeakErr=array('d')
        
        f=TF1("f_"+str(ii)+"_first","gaus",ampl-2.5*sigma,ampl+2.5*sigma)
        hhL0[ii].Fit(f,"R+","",ampl-2.5*sigma,ampl+2.5*sigma)
        xPeak.append(f.GetParameter(1))
        xPeakErr.append(f.GetParError(1))
        nPeak.append(1)
        nPeakErr.append(0)
        for jj in range(n_found):        
            x=xxP[jj]
            #y=yP[jj]
         
            f=TF1("f_"+str(ii)+"_"+str(jj),"gaus",x-2.5*sigma,x+2.5*sigma)
            hhL0[ii].Fit(f,"R+","",)
            xPeak.append(f.GetParameter(1))
            xPeakErr.append(f.GetParError(1))
            nPeak.append(jj+2)
            nPeakErr.append(0)
              
        g=TGraphErrors(n_found+1,nPeak,xPeak,nPeakErr,xPeakErr)
        g.SetName("g_L0_"+str(ii+1))
        g.Fit("pol1")
        g.SetMarkerColor(2)
        g.SetMarkerStyle(20)
        gainL0[ii]=g.GetFunction("pol1").GetParameter(1)
        
        fileOut.cd()
        hhL0[ii].Write();
        g.Write()


    ###L1### 
    #Get histograms
    for ii in range(0,10):
        h=file.Get("hBDXMiniVetoSIPMA_L1_C"+str(ii+1))
        hhL1.append(h)
        

    for ii in range(0,10):
        
        h=hhL1[ii].Clone("hBDXMiniVetoSIPMA_L1_C"+str(ii+1)+"clone")
        min=0;
        max=0;
        sigma=0;
        ampl=0;
        if (ii<8):
            min=amplONE*1.5
            max=amplONE*nSPECTRUM
            sigma=sigmaONE
            ampl=amplONE
        else:
            min=(amplONE*1.5)/2
            max=(amplONE*nSPECTRUM)/2
            sigma=sigmaONE/2
            ampl=amplONE/2

        h.Smooth()
        h.GetXaxis().SetRangeUser(min,max);
        s=TSpectrum();
        n_found=s.Search(h,2,"",0.1);
        pm = h.GetListOfFunctions().FindObject("TPolyMarker")
        hhL1[ii].GetListOfFunctions().Add(pm)
    
        xP=s.GetPositionX();
        yP=s.GetPositionY();
        xxP=[]
        xxP=searchClusters(n_found,xP,yP,2*sigma)
        n_found=len(xxP)
        
        
        print n_found
        print xxP
        
        xPeak=array('d')
        xPeakErr=array('d')
        nPeak=array('d')
        nPeakErr=array('d')
        
        f=TF1("f_"+str(ii)+"_first","gaus",ampl-2.5*sigma,ampl+2.5*sigma)
        hhL1[ii].Fit(f,"R+","",ampl-2.5*sigma,ampl+2.5*sigma)
        xPeak.append(f.GetParameter(1))
        xPeakErr.append(f.GetParError(1))
        nPeak.append(1)
        nPeakErr.append(0)
        for jj in range(n_found):        
            x=xxP[jj]
            #y=yP[jj]
            print("fit")
            f=TF1("f_"+str(ii)+"_"+str(jj),"gaus",x-2.5*sigma,x+2.5*sigma)
            hhL1[ii].Fit(f,"R+","",)
            xPeak.append(f.GetParameter(1))
            xPeakErr.append(f.GetParError(1))
            nPeak.append(jj+2)
            nPeakErr.append(0)
              
        g=TGraphErrors(n_found+1,nPeak,xPeak,nPeakErr,xPeakErr)
        g.SetName("g_L1_"+str(ii+1))
        g.Fit("pol1")

        gainL1[ii]=g.GetFunction("pol1").GetParameter(1)
        
        fileOut.cd()
        hhL1[ii].Write();
        g.SetMarkerStyle(20);
        g.SetMarkerColor(2);
        g.Write()
        
    fileOut.Write();
    fileOut.Close();

    ocalib = open(ofname,"w+")
    ocalib.write("#sector layer component readout ampl\n")
    for isector in range(0,20):
        for ilayer in range(0,4):
            for icomponent in range(0,20):
                for ireadout in range (1,6):
                    amp=1.;
                    if (isector==0)and(ilayer==1)and(icomponent>=1)and(icomponent<=10)and(ireadout==1):
                        amp=gainL1[icomponent-1]
                    if (isector==0)and(ilayer==0)and(icomponent>=1)and(icomponent<=10)and(ireadout==1):
                        amp=gainL0[icomponent-1]
                    
                    ocalib.write(str(isector)+" "+str(ilayer)+" "+str(icomponent)+" "+str(ireadout)+" "+str(amp)+"\n")
    ocalib.close()


calibAmpl("data/1338.root","data/1338.calibA")
