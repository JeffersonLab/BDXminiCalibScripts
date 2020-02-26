from ROOT import *
import numpy
import math
from array import array
from searchClusters import searchClusters



def calibCharge(fname,ofname):
    fileOut=TFile(ofname+".root","recreate");
    file=TFile(fname)

    #SET HERE SOME REASONABLE VALUES
    amplONE=90.
    sigmaONE=20.
    nSPECTRUM=4.
    nSIGMAL=1.5
    nSIGMAR=2.2

    hhL0=[]
    hhL1=[]
    gainL0=[1]*10
    gainL1=[1]*10
    pedL0=[0]*10
    pedL1=[0]*10
    ###L0### 
    #Get histograms
    for ii in range(0,10):
        h=file.Get("hBDXMiniVetoSIPMQ_L0_C"+str(ii+1))
        hhL0.append(h)
        

    for ii in range(0,10):
        if (ii==2):
            continue
        
        h=hhL0[ii].Clone("hBDXMiniVetoSIPMQ_L0_C"+str(ii+1)+"clone")
        min=0;
        max=0;
        sigma=0;
        ampl=0;

        if (ii<8):
            sigma=sigmaONE
            ampl=amplONE
        else:
            sigma=sigmaONE/2
            ampl=amplONE/2

        #PRE-FIT for sigma
        max=h.GetBinCenter(h.GetMaximumBin())
        h.GetXaxis().SetRangeUser(max-2*sigma,max+2*sigma)
        h.Fit("gaus")
        ampl=h.GetFunction("gaus").GetParameter(1)
        sigma=h.GetFunction("gaus").GetParameter(2)

        min=ampl*1.8;
        max=ampl*nSPECTRUM;


#        h.Smooth()
        h.GetXaxis().SetRangeUser(min,max);
        s=TSpectrum();
        n_found=s.Search(h,2,"",0.1);
        h.GetXaxis().UnZoom();
        pm = h.GetListOfFunctions().FindObject("TPolyMarker")
        hhL0[ii].GetListOfFunctions().Add(pm)

        xP=s.GetPositionX();
        yP=s.GetPositionY();

        xxP=[]
        xxP=searchClusters(n_found,xP,yP,3*sigma)
        n_found=len(xxP)
   
        xPeak=array('d')
        xPeakErr=array('d')
        nPeak=array('d')
        nPeakErr=array('d')
        
        hhL0[ii].Fit("gaus","LR+","",ampl-nSIGMAL*sigma,ampl+nSIGMAR*sigma)
        xPeak.append(hhL0[ii].GetFunction("gaus").GetParameter(1))
        xPeakErr.append(hhL0[ii].GetFunction("gaus").GetParError(1))
        nPeak.append(1)
        nPeakErr.append(0)
        for jj in range(n_found):  
            msigma=sigma*math.sqrt(jj+1)
            x=xxP[jj]
            #y=yP[jj]
            print("fit")
            f=TF1("f_"+str(ii)+"_"+str(jj),"gaus",x-nSIGMAL*msigma,x+nSIGMAR*msigma)
            hhL0[ii].Fit(f,"LR+","",)
            xPeak.append(f.GetParameter(1))
            xPeakErr.append(f.GetParError(1))
            nPeak.append(jj+2)
            nPeakErr.append(0)
              
        g=TGraphErrors(n_found+1,nPeak,xPeak,nPeakErr,xPeakErr)
        g.SetName("g_L0_"+str(ii+1))
        g.Fit("pol1")
        g.SetMarkerStyle(20)
        g.SetMarkerColor(2)

        gainL0[ii]=g.GetFunction("pol1").GetParameter(1)
        pedL0[ii]=g.GetFunction("pol1").GetParameter(0)
        
        fileOut.cd()
        hhL0[ii].Write();
        g.Write()


    ###L1### 
    #Get histograms
    for ii in range(0,10):
        h=file.Get("hBDXMiniVetoSIPMQ_L1_C"+str(ii+1))
        hhL1.append(h)
        

    for ii in range(0,10):
        
        h=hhL1[ii].Clone("hBDXMiniVetoSIPMQ_L1_C"+str(ii+1)+"clone")
        min=0;
        max=0;
        sigma=0;
        ampl=0;
        if (ii<8):
            sigma=sigmaONE
            ampl=amplONE
        else:
            sigma=sigmaONE/2
            ampl=amplONE/2

        #PRE-FIT for sigma
        max=h.GetBinCenter(h.GetMaximumBin())
        h.GetXaxis().SetRangeUser(max-2*sigma,max+2*sigma)
        h.Fit("gaus")
        ampl=h.GetFunction("gaus").GetParameter(1)
        sigma=h.GetFunction("gaus").GetParameter(2)
        
        min=ampl*1.8;
        max=ampl*nSPECTRUM;

#        h.Smooth()
        h.GetXaxis().SetRangeUser(min,max);
        s=TSpectrum();
        n_found=s.Search(h,2,"",0.1);
        functions = h.GetListOfFunctions();
        pm = functions.FindObject("TPolyMarker");
        hhL1[ii].Draw()
        hhL1[ii].GetListOfFunctions().Add(pm)
        pm.Draw("SAMES")
     
        
        xP=s.GetPositionX();
        yP=s.GetPositionY();
   
        xxP=[]
        xxP=searchClusters(n_found,xP,yP,3*sigma)
        n_found=len(xxP)
        
        xPeak=array('d')
        xPeakErr=array('d')
        nPeak=array('d')
        nPeakErr=array('d')

       
        hhL1[ii].Fit("gaus","LR+","",ampl-nSIGMAL*sigma,ampl+nSIGMAR*sigma)
        xPeak.append(hhL1[ii].GetFunction("gaus").GetParameter(1))
        xPeakErr.append(hhL1[ii].GetFunction("gaus").GetParError(1))
        nPeak.append(1)
        nPeakErr.append(0)
        for jj in range(n_found): 
            msigma=sigma*math.sqrt(jj+1)
            x=xxP[jj]
            #y=yP[jj]
            print("fit")
            f=TF1("f_"+str(ii)+"_"+str(jj),"gaus",x-nSIGMAL*msigma,x+nSIGMAR*msigma)
            hhL1[ii].Fit(f,"LR+","",)
            xPeak.append(f.GetParameter(1))
            xPeakErr.append(f.GetParError(1))
            nPeak.append(jj+2)
            nPeakErr.append(0)
              
        g=TGraphErrors(n_found+1,nPeak,xPeak,nPeakErr,xPeakErr)
        g.SetName("g_L1_"+str(ii+1))
        g.Fit("pol1")
        g.SetMarkerStyle(20)
        g.SetMarkerColor(2)

        gainL1[ii]=g.GetFunction("pol1").GetParameter(1)
        pedL1[ii]=g.GetFunction("pol1").GetParameter(0)
        
        fileOut.cd()
        hhL1[ii].Write();
        g.Write()
        
    fileOut.Write();
    fileOut.Close();

    ocalib = open(ofname,"w+")
    ocalib.write("#sector layer component readout gain ped\n")
    for isector in range(0,20):
        for ilayer in range(0,4):
            for icomponent in range(0,20):
                for ireadout in range (1,6):
                    amp=1.;
                    ped=0.;
                    if (isector==0)and(ilayer==1)and(icomponent>=1)and(icomponent<=10)and(ireadout==1):
                        amp=gainL1[icomponent-1]
                        ped=pedL1[icomponent-1]
                    if (isector==0)and(ilayer==0)and(icomponent>=1)and(icomponent<=10)and(ireadout==1):
                        amp=gainL0[icomponent-1]
                        ped=pedL0[icomponent-1]
                    ocalib.write(str(isector)+" "+str(ilayer)+" "+str(icomponent)+" "+str(ireadout)+" "+str(amp)+" "+str(ped)+"\n")
    ocalib.close()

