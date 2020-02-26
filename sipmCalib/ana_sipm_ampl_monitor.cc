
//This function, given the x and y positions of found peaks, determine clusters of nearby peaks and returns the SORTED array of peaks. 
//The minimum distance to merge is thr.
vector<double> searchClusters(int n,double *x,double *y,double thr){
  //First, sort by peak height in descending orde
  double *mx=new double[n];
  map<double, double, greater < double > > mymap;
  map<double, double>::iterator it; 
  int idx=0;
  vector<double> clusters;
  

  for (int ii=0;ii<n;ii++){
    mymap[y[ii]]=x[ii];
  }
  
  for (it=mymap.begin();it!=mymap.end();it++){
    mx[idx++]=(*it).second;
  }
  
  for (int ii=0;ii<n;ii++){
    bool flag=false;
    double x1=mx[ii];
    for (int jj=0;jj<clusters.size();jj++){
      double x2=clusters[jj];
      if (fabs(x2-x1)<thr){
	flag=true;
	break;
      }
    }
    if (flag==false){
      clusters.push_back(x1);
    }
  }

  sort(clusters.begin(),clusters.end());
  return clusters;
}


void ana_sipm_ampl_monitor(string fname,string ofname){

  TFile *fileOut=new TFile(Form("%s.root",ofname.c_str()),"recreate");
  TFile *file=new TFile(fname.c_str());

  //SET HERE SOME REASONABLE VALUES
  double amplONE=10.;
  double sigmaONE=2.;
  int nSPECTRUM=5.;
    
  TH1D *hhL0[10];
  TH1D *hhL1[10];
  
  double gainL0[10];
  double gainL1[10];

 
  double *xPeak,*xPeakErr,*nPeak,*nPeakErr;
  double min=0;
  double max=0;
  double sigma=0;
  double ampl=0;
  double x;
  TH1D *h;
  TSpectrum *s;
  int n_found;
  TPolyMarker *pm;
  
  double *xP,*yP;
  vector<double> xxP;

  TGraphErrors *g;
  TF1 *f;


  //Get Histograms
  for (int ii=0;ii<10;ii++){
    h=(TH1D*)file->Get(Form("hBDXMiniVetoSIPMA_L0_C%i",ii+1));
    hhL0[ii]=h;
    h=(TH1D*)file->Get(Form("hBDXMiniVetoSIPMA_L1_C%i",ii+1));
    hhL1[ii]=h;
    gainL0[ii]=1.;
    gainL1[ii]=1.;
  }
  




  for (int ii=0;ii<10;ii++){
    if (ii==2) continue;

    h=(TH1D*)hhL0[ii]->Clone(Form("hBDXMiniVetoSIPMA_L0_C%i_clone",ii+1));
   
    if (ii<8){
      min=amplONE*1.5;
      max=amplONE*nSPECTRUM;
      sigma=sigmaONE;
      ampl=amplONE;
    }
    else{
      min=(amplONE*1.5)/2;
      max=(amplONE*nSPECTRUM)/2;
      sigma=sigmaONE/2;
      ampl=amplONE/2;
    }

    h->Smooth();
    h->GetXaxis()->SetRangeUser(min,max);
    s=new TSpectrum();
    n_found=s->Search(h,2,"",0.1);
    pm = (TPolyMarker*)h->GetListOfFunctions()->FindObject("TPolyMarker");
    hhL0[ii]->GetListOfFunctions()->Add(pm);

    xP=s->GetPositionX();
    yP=s->GetPositionY();
      
    xxP=searchClusters(n_found,xP,yP,2*sigma);
    n_found=xxP.size();

    xPeak=new double[n_found+1];
    xPeakErr=new double[n_found+1];
    nPeak=new double[n_found+1];
    nPeakErr=new double[n_found+1];

    f=new TF1(Form("f0_%i_first",ii),"gaus",ampl-2.5*sigma,ampl+2.5*sigma);
    hhL0[ii]->Fit(f,"R+","",ampl-2.5*sigma,ampl+2.5*sigma);
    xPeak[0]=f->GetParameter(1);
    xPeakErr[0]=f->GetParError(1);
    nPeak[0]=1.;
    nPeakErr[0]=0.;
  
    for (int jj=0;jj<n_found;jj++){
      x=xxP[jj];
      f=new TF1(Form("f0_%i_%i",ii,jj),"gaus",x-2.5*sigma,x+2.5*sigma);
      hhL0[ii]->Fit(f,"R+","");
      xPeak[jj+1]=f->GetParameter(1);
      xPeakErr[jj+1]=f->GetParError(1);
      nPeak[jj+1]=jj+2;
      nPeakErr[jj+1]=0.0;
    }

    g=new TGraphErrors(n_found+1,nPeak,xPeak,nPeakErr,xPeakErr);
    g->SetName(Form("g_L0_%i",ii+1));
    g->Fit("pol1");
    g->SetMarkerColor(2);
    g->SetMarkerStyle(20);
    gainL0[ii]=g->GetFunction("pol1")->GetParameter(1);
    
   

    fileOut->cd();
    hhL0[ii]->Write();
    g->Write();

  }

  for (int ii=0;ii<10;ii++){
   
    h=(TH1D*)hhL1[ii]->Clone(Form("hBDXMiniVetoSIPMA_L1_C%i_clone",ii+1));
   
    if (ii<8){
      min=amplONE*1.5;
      max=amplONE*nSPECTRUM;
      sigma=sigmaONE;
      ampl=amplONE;
    }
    else{
      min=(amplONE*1.5)/2;
      max=(amplONE*nSPECTRUM)/2;
      sigma=sigmaONE/2;
      ampl=amplONE/2;
    }

    h->Smooth();
    h->GetXaxis()->SetRangeUser(min,max);
    s=new TSpectrum();
    n_found=s->Search(h,2,"",0.1);
    pm = (TPolyMarker*)h->GetListOfFunctions()->FindObject("TPolyMarker");
    hhL1[ii]->GetListOfFunctions()->Add(pm);

    xP=s->GetPositionX();
    yP=s->GetPositionY();
      
    xxP=searchClusters(n_found,xP,yP,2*sigma);
    n_found=xxP.size();

    xPeak=new double[n_found+1];
    xPeakErr=new double[n_found+1];
    nPeak=new double[n_found+1];
    nPeakErr=new double[n_found+1];

    f=new TF1(Form("f1_%i_first",ii),"gaus",ampl-2.5*sigma,ampl+2.5*sigma);
    hhL1[ii]->Fit(f,"R+","",ampl-2.5*sigma,ampl+2.5*sigma);
    xPeak[0]=f->GetParameter(1);
    xPeakErr[0]=f->GetParError(1);
    nPeak[0]=1.;
    nPeakErr[0]=0.0;
  
    for (int jj=0;jj<n_found;jj++){
      x=xxP[jj];
      f=new TF1(Form("f1_%i_%i",ii,jj),"gaus",x-2.5*sigma,x+2.5*sigma);
      hhL1[ii]->Fit(f,"R+","");
      xPeak[jj+1]=f->GetParameter(1);
      xPeakErr[jj+1]=f->GetParError(1);
      nPeak[jj+1]=jj+2;
      nPeakErr[jj+1]=0.0;
    }

    g=new TGraphErrors(n_found+1,nPeak,xPeak,nPeakErr,xPeakErr);
    g->SetName(Form("g_L1_%i",ii+1));
    g->Fit("pol1");
    g->SetMarkerColor(2);
    g->SetMarkerStyle(20);
    gainL0[ii]=g->GetFunction("pol1")->GetParameter(1);
        
    fileOut->cd();
    hhL1[ii]->Write();
    g->Write();
  }

  fileOut->Write();
 
  ofstream ocalib(ofname.c_str());
  ocalib<<"#sector layer component readout ampl"<<endl;

  for (int isector=0;isector<20;isector++){
    for (int ilayer=0;ilayer<4;ilayer++){
      for (int icomponent=0;icomponent<20;icomponent++){
	for (int ireadout=1;ireadout<6;ireadout++){
	  double amp=1;

	  cout<<isector<<" "<<ilayer<<" "<<icomponent<<" "<<ireadout<<endl;
	  if ((isector==0)and(ilayer==0)and(icomponent>=1)and(icomponent<=10)and(ireadout==1)){
	    amp=gainL0[icomponent-1];
	  }	  
	  if ((isector==0)and(ilayer==1)and(icomponent>=1)and(icomponent<=10)and(ireadout==1)){
	    amp=gainL1[icomponent-1];
	  }
	}
      }
    }
  }
  ocalib.close();
}
