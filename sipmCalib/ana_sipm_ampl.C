#include <iostream>
#include <map>
#include <TH1D.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>
#include <fstream>
#include <TSpectrum.h>
using namespace std;

typedef struct{
  int sector;
  int layer;
  int component;
  int readout;
} INT_VETO_IndexLight_t;



bool operator<(const INT_VETO_IndexLight_t &a, const INT_VETO_IndexLight_t &b) {
	if (a.sector < b.sector) return true;
	if (a.sector > b.sector) return false;
	if (a.layer < b.layer) return true;
	if (a.layer > b.layer) return false;
	if (a.component < b.component) return true;
	if (a.component > b.component) return false;
	if (a.readout < b.readout) return true;
	if (a.readout > b.readout) return false;

	return false;
}

std::map<INT_VETO_IndexLight_t,TH1D *> h_map;
std::map<INT_VETO_IndexLight_t,TH1D *>::iterator h_map_it;

std::map<INT_VETO_IndexLight_t,TF1 *> f_map;
std::map<INT_VETO_IndexLight_t,TF1 *>::iterator f_map_it;

double three_gaus(double *xx,double *par){
  double x=xx[0];
  
  double ped=par[0];
  double singlephe=par[1];
  double twophe=par[2];
  double A1=par[3];
  double A2=par[4];
  double A3=par[5];

  double sig1=par[6];
  double sig2=par[7];
  double sig3=par[8];

  double g1=A1*TMath::Gaus(x,ped,sig1);
  double g2=A2*TMath::Gaus(x,ped+singlephe,sig2);
  double g3=A3*TMath::Gaus(x,ped+twophe,sig3);
  
  double bck=par[9];
  
  return g1+g2+g3+bck;
}

void ana_sipm_ampl(string fname,string oname){
 
  static const int NBINS = 230.;
  static const double XLOW = 0.;
  static const double XUP = 63.;

  double flow,fup;
  double singlephe,twophe;
  string outname=oname;
  outname+=".root";
  string outCalibname=oname; 
  outCalibname=oname;
  string outCalibname_ideal=oname;
  outCalibname_ideal+=".ideal";
  
  ofstream ocalib(outCalibname.c_str());
  ofstream ocalib_ideal(outCalibname_ideal.c_str());

  TFile *file1=new TFile(outname.c_str(),"recreate");
  TFile *file=new TFile(fname.c_str());  
  TTree *t=(TTree*)file->Get("IntVeto_SipmCalib");

  INT_VETO_IndexLight_t index;
  double Q,T,Ampl;
  int n_found;
  TH1D *h;
  TF1 *f;
  TSpectrum *s;
  Double_t *xP,*yP;



   uint tWord;
  double A,ped,xlow,xup,single,two,sig1,sig2,sig3,A1,A2,A3;

  t->SetBranchAddress("sector",&index.sector);
  t->SetBranchAddress("layer",&index.layer);
  t->SetBranchAddress("component",&index.component);
  t->SetBranchAddress("readout",&index.readout);
  t->SetBranchAddress("A",&Ampl);
  t->SetBranchAddress("T",&T);
 t->SetBranchAddress("tWord",&tWord);
  int N=t->GetEntries();
  cout<<"There are: "<<N<<" entries " <<endl;
  
   for (int ii=0;ii<N;ii++){
    t->GetEntry(ii);
  //select random trg only
    if (((tWord>>31)&0x1)==0) continue;
    
    if ((index.sector==0)&&(index.layer==0)&&(index.component==3)) continue;

    if ((Ampl>XUP)||(Ampl<XLOW)) continue;

    if ((T<30)||(T>2000)) continue;
    h_map_it=h_map.find(index);
    if (h_map_it==h_map.end()){
      h_map.insert(std::make_pair(index,new TH1D(Form("h_%i_%i_%i_%i",index.sector,index.layer,index.component,index.readout),Form("h_%i_%i_%i_%i",index.sector,index.layer,index.component,index.readout),NBINS,XLOW,XUP)));
      f_map.insert(std::make_pair(index,new TF1(Form("f_%i_%i_%i_%i",index.sector,index.layer,index.component,index.readout),three_gaus,XLOW,XUP,10)));

    }
    h=h_map[index];
    h->Fill(Ampl);
  }
  
  file1->cd();
  for (h_map_it=h_map.begin();h_map_it!=h_map.end();h_map_it++){
    h=(h_map_it)->second;
    f=f_map.at(h_map_it->first);
    index=(h_map_it)->first;

    //First fit to determine A1, sigma1, ped
    TF1 g1("g1","gaus");
    A1=h->GetMaximum();
    ped=h->GetBinCenter(h->GetMaximumBin());
    g1.SetParameter(0,A1);
    g1.SetParameter(1,ped);
    h->GetXaxis()->SetRangeUser(ped-4,ped+4);
    sig1=h->GetRMS();
    h->GetXaxis()->SetRangeUser(XLOW,XUP);
    g1.SetParameter(2,sig1);
    h->Fit("g1","RQ","",ped-3*sig1,ped+3*sig1);
    A1=h->GetFunction("g1")->GetParameter(0);
    ped=h->GetFunction("g1")->GetParameter(1);
    sig1=h->GetFunction("g1")->GetParameter(2);
    
    /*now use TSpectrum to find peaks, but go away from the pedestal*/
    xup=XUP;
    xlow=ped+3*sig1;
   
    h->GetXaxis()->SetRangeUser(xlow,xup);
    /*Use TSpectrum to find peaks*/
    s=new TSpectrum(2); //max 2 peaks
    n_found=s->Search(h,2,"",0.05);
    h->GetXaxis()->SetRangeUser(XLOW,XUP);
    cout<<"Found "<<n_found<<" peaks "<<endl;
    xP=s->GetPositionX();
    yP=s->GetPositionY();
    
    if (n_found<=0){
      cout<<"Error "<<endl;
      continue;
    }
    
    else if (n_found>=1){ //best-case scenario
      /*Fit again with gaus the first peak*/
      sig2=sig1;
      TF1 g2("g2","gaus",xP[0]-3*sig2,xP[0]+3*sig2);
      h->Fit("g2","","RQ",xP[0]-3*sig2,xP[0]+3*sig2);
      A2=h->GetFunction("g2")->GetParameter(0);
      single=h->GetFunction("g2")->GetParameter(1);
      single-=ped;
      sig2=h->GetFunction("g2")->GetParameter(2);
      
      /*The third peak position is not reliable*/
      two=2*single+ped;
      h->GetXaxis()->SetRangeUser(two-2*sig2,two+5*sig2);
      A3=h->GetMaximum();
      sig3=sig2;
      two=h->GetBinCenter(h->GetMaximumBin());
      h->GetXaxis()->SetRangeUser(XLOW,XUP);
      cout<<"g3Fit: "<<two-3*sig3<<" "<<two+3*sig3<<endl;
      TF1 g3("g3","gaus",two-3*sig3,two+3*sig3);
      h->Fit("g3","","RQ",two-3*sig3,two+3*sig3);
      A3=h->GetFunction("g3")->GetParameter(0);
      two=h->GetFunction("g3")->GetParameter(1);
      two-=ped;
      sig3=h->GetFunction("g3")->GetParameter(2);      
      flow=ped-5*sig1;
      fup=two+ped+3*sig3;
      cout<<"DATA: "<<index.component<<" "<<index.readout<<" "<<flow<<" "<<fup<<" "<<ped<<" "<<single<<" "<<two<<" "<<sig1<<" "<<sig2<<" "<<sig3<<endl;
    }








 f->SetParameter(0,ped);
 f->SetParameter(1,single);
 f->SetParameter(2,two);
 f->SetParameter(3,A1);      //ampl1
 f->SetParameter(4,A2);   //ampl2
 f->SetParameter(5,A3);   //ampl3
 f->SetParameter(6,sig1);    //sigma1
 f->SetParameter(7,sig2);    //sigma2
 f->SetParameter(8,sig3);    //sigma3
    
    f->SetParameter(9,0);
    f->SetNpx(1000);
    f->SetRange(flow,fup);
    h->Fit(f->GetName(),"R","",flow,fup);
    h->Write();
  }
  
  ocalib<<"#sector layer component readout ampl"<<endl;
  double ampl;
  //  for (h_map_it=h_map.begin();h_map_it!=h_map.end();h_map_it++){
  for (int isector=0;isector<20;isector++){
    for (int ilayer=0;ilayer<4;ilayer++){
      for (int icomponent=0;icomponent<20;icomponent++){
	for (int ireadout=1;ireadout<=5;ireadout++){ //readout=0 is meaningless!

	  
	  index.sector=isector;
	  index.layer=ilayer;
	  index.component=icomponent;
	  index.readout=ireadout;
	  f_map_it=f_map.find(index);
	  if (f_map_it==f_map.end()) ampl=0;
	  else{
	    f=f_map_it->second;
	    ampl=fabs(f->GetParameter(2)-f->GetParameter(1));
	  }
	  ocalib<<isector<<" "<<ilayer<<" "<<icomponent<<" "<<ireadout<<" "<<ampl<<endl;
	  ampl=0;
	  ocalib_ideal<<isector<<" "<<ilayer<<" "<<icomponent<<" "<<ireadout<<" "<<ampl<<endl;
	}
      }
    }
  }
  ocalib.close();
  
  file1->Write();
  file1->Close();
}

