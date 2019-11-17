#include <iostream>
#include <map>
#include <TH1D.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <fstream>
#include <TLinearFitter.h>

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


std::map<INT_VETO_IndexLight_t,TF1 *> f0_map;
std::map<INT_VETO_IndexLight_t,TF1 *>::iterator f0_map_it;

std::map<INT_VETO_IndexLight_t,TF1 *> f_map;
std::map<INT_VETO_IndexLight_t,TF1 *>::iterator f_map_it;

std::map<INT_VETO_IndexLight_t,double> gain_map;
std::map<INT_VETO_IndexLight_t,double>::iterator gain_map_it;
std::map<INT_VETO_IndexLight_t,double> ped_map;
std::map<INT_VETO_IndexLight_t,double>::iterator ped_map_it;



double three_gaus(double *xx,double *par){
  double x=xx[0];
  
  double one=par[0];
  double two=par[1];
  double three=par[2];
  double A1=par[3];
  double A2=par[4];
  double A3=par[5];

  double sig1=par[6];
  double sig2=par[7];
  double sig3=par[8];

  double g1=A1*TMath::Gaus(x,one,sig1);
  double g2=A2*TMath::Gaus(x,two,sig2);
  double g3=A3*TMath::Gaus(x,three,sig3);

  return g1+g2+g3;
}

void ana_sipm_charge(string fname,string oname){
 
  static const int NBINS = 320;
  static const double XLOW = -0.5E3;
  static const double XUP = 1600;

  //static const int NBINS = 130;
  //static const double XLOW = -0.1E3;
  //static const double XUP = 500;

  double TMIN=200;
  double TMAX=1800;

  double flow,fup;
  double singlephe;
  string outname=oname; //".root"
  outname+=".root";

  string outCalibname=oname; 
  string outCalibname_ideal=oname;
  outCalibname_ideal+=".ideal";
  
  ofstream ocalib(outCalibname.c_str());
  ofstream ocalib_ideal(outCalibname_ideal.c_str());

  TFile *file1=new TFile(outname.c_str(),"recreate");
  TFile *file=new TFile(fname.c_str());  
  TTree *t=(TTree*)file->Get("IntVeto_SipmCalib");

  INT_VETO_IndexLight_t index;
  double Q,T;
 
  TH1D *h;
  TF1 *f,*f0;
  TLinearFitter *lf=new TLinearFitter();lf->SetFormula("1 ++ x");
  double xTMP[1];
  TSpectrum *s;
  Double_t *xP,*yP;
  int n_found;

  double A,sig,ped;
  
  double A1,A2,A3,sig1,sig2,sig3,one,two,three;
  double xup,xlow;
  double gain;

  uint tWord;

  t->SetBranchAddress("sector",&index.sector);
  t->SetBranchAddress("layer",&index.layer);
  t->SetBranchAddress("component",&index.component);
  t->SetBranchAddress("readout",&index.readout);
  t->SetBranchAddress("Qraw",&Q);
  t->SetBranchAddress("T",&T);
  t->SetBranchAddress("tWord",&tWord);
  int N=t->GetEntries();
  cout<<"There are: "<<N<<" entries " <<endl;
  
   for (int ii=0;ii<N;ii++){
    t->GetEntry(ii);
    
    //select random trg only
    //  if (((tWord>>31)&0x1)==0) continue;

    if ((T<TMIN)||(T>TMAX)) continue;
    h_map_it=h_map.find(index);
    if (h_map_it==h_map.end()){
      h_map.insert(std::make_pair(index,new TH1D(Form("h_%i_%i_%i_%i",index.sector,index.layer,index.component,index.readout),Form("h_%i_%i_%i_%i",index.sector,index.layer,index.component,index.readout),NBINS,XLOW,XUP)));
      f_map.insert(std::make_pair(index,new TF1(Form("f_%i_%i_%i_%i",index.sector,index.layer,index.component,index.readout),three_gaus,XLOW,XUP,9)));
      f0_map.insert(std::make_pair(index,new TF1(Form("gaus_%i_%i_%i_%i",index.sector,index.layer,index.component,index.readout),"gaus",XLOW,XUP)));     
    }
    //    cout<<index.component<<" "<<index.readout<<" "<<Q<<endl;
    //    cin.get();
    h=h_map[index];
    h->Fill(Q);
  }
  
  file1->cd();
  for (h_map_it=h_map.begin();h_map_it!=h_map.end();h_map_it++){
    h=(h_map_it)->second;
    f=f_map.at(h_map_it->first);
    f0=f0_map.at(h_map_it->first);
    h->Draw();
    /*First a 1D gaussian fit*/		     
    f0->SetParameter(1,h->GetBinCenter(h->GetMaximumBin()));
    f0->SetParameter(2,h->GetRMS());
    f0->SetParameter(0,h->GetMaximum());
    h->Fit(f0,"","",h->GetBinCenter(h->GetMaximumBin())-100,h->GetBinCenter(h->GetMaximumBin())+50);
    f0->SetLineColor(3);
    f0->Draw("SAME");

    one=f0->GetParameter(1);
    A1=f0->GetParameter(0);
    sig1=f0->GetParameter(2);

    /*Then set histo range out of 1 phe*/
    index=(h_map_it)->first;
    xlow=one+3.3*sig1;
    if (index.component==0){
      // if (index.readout==2) xup=250;
      //   else xup=550;
      xup=450;
    }
    else if (index.component==1){
      //if (index.readout==2) xup=250;
      //	else xup=550;
      xup=450;
    }
    else if (index.component==2){ 
      //if (index.readout==2) xup=250;
      //	else xup=550;    
      xup=450;    
    }
    else if (index.component==3){ 
      //  if (index.readout==1) xup=ped+1200;
      //  else xup=ped+1100;
      xup=450;
    }
    else if (index.component==4){ 
       xup=450;    
    }
    else if (index.component==5){ 
      xup=450;    
    }
    cout<<"Range: "<<index.component<<" "<<index.readout<<" "<<xlow<<" "<<xup<<endl;
    h->GetXaxis()->SetRangeUser(xlow,xup);

    /*Use TSpectrum to find peaks*/
    s=new TSpectrum(2); //max 3 peaks
    n_found=s->Search(h,2,"nobackground new",0.05);
   
    h->GetXaxis()->SetRangeUser(XLOW,XUP);
    cout<<"Found "<<n_found<<" peaks "<<endl;
    xP=s->GetPositionX();
    yP=s->GetPositionY();

    if (n_found<=0){
      cout<<"Error "<<endl;
      continue;
    }
    else if (n_found==1){
	 A2=yP[0];
	 A3=A2/2;
	
	 two=xP[0];
	 three=one+two;

	 
	 sig2=1.1*sig1;
	 sig3=1.2*sig1;
	 
	 flow=0.5*one;
	 fup= three+(two-one)*0.5;

    }
    else if (n_found==2){ //best-case scenario

      if ((yP[1]-yP[0])>50){
     
	A2=yP[0];
	A3=yP[1];

     
	two=xP[0];
	three=xP[1];
	 
      
	sig2=1.1*sig1;
	sig3=1.2*sig1;
      
	flow=0.5*one;
	fup= three+(two-one)*.75;

      }
      else{
	 A2=yP[0];
	 A3=A2/2;
	
	 two=xP[0];
	 three=xP[1];

	 
	 sig2=1.1*sig1;
	 sig3=1.2*sig1;

	 flow=0.5*one;
	 fup= three+(two-one)*.5;
      }

    }
    delete s;
    
    
    f->SetParameter(0,one);
    f->SetParameter(1,two);
    f->SetParameter(2,three);

    f->SetParameter(3,A1);
    f->SetParameter(4,A2);
    f->SetParameter(5,A3);

    f->SetParameter(6,sig1);
    f->SetParameter(7,sig2);
    f->SetParameter(8,sig3);
    
  

    h->Fit(f->GetName(),"","",flow,fup);
    f0->Draw("SAME");
    //   f0->Write();
    h->Write();
    
    xTMP[0]=0.;lf->AddPoint(xTMP,f->GetParameter(0));
    xTMP[0]=1.;lf->AddPoint(xTMP,f->GetParameter(1));
    xTMP[0]=2.;lf->AddPoint(xTMP,f->GetParameter(2));

    lf->Eval();
    cout<<"PED: "<<lf->GetParameter(0)<<" PHE: "<<lf->GetParameter(1)<<endl;

    gain_map[index]=lf->GetParameter(1);
    ped_map[index]=lf->GetParameter(0);
    lf->ClearPoints();

  }
  
  ocalib<<"#sector layer component readout gain ped"<<endl;

  //  for (h_map_it=h_map.begin();h_map_it!=h_map.end();h_map_it++){
  for (int isector=0;isector<20;isector++){
    for (int ilayer=0;ilayer<4;ilayer++){
      for (int icomponent=0;icomponent<20;icomponent++){
	for (int ireadout=1;ireadout<=5;ireadout++){ //readout=0 is meaningless!

	  
	  index.sector=isector;
	  index.layer=ilayer;
	  index.component=icomponent;
	  index.readout=ireadout;
	  gain_map_it=gain_map.find(index);
	  if (gain_map_it==gain_map.end()){ gain=1; ped=0;}
	  else{
	    gain=gain_map[index];
	    ped=ped_map[index];
	  }
	  ocalib<<isector<<" "<<ilayer<<" "<<icomponent<<" "<<ireadout<<" "<<gain<<" "<<ped<<endl;
	  gain=1;
	  ped=0;
	  ocalib_ideal<<isector<<" "<<ilayer<<" "<<icomponent<<" "<<ireadout<<" "<<gain<<" "<<ped<<endl;
	}
      }
    }
  }
  ocalib.close();
  
  file1->Write();
  file1->Close();
}

