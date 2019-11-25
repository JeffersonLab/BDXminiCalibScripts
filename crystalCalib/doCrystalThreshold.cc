#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"


#include <string>
#include <iostream>
#include <fstream>

using namespace std;

double fitFunE(double *x,double *par){
  double E=x[0];
  double thr=par[0];
  double sigma=par[1];
  double eff=par[2];
  double sig=0;

  sig=eff/(1+exp(-(E-thr)/sigma));
  
  return sig;

}

double fitFunQ(double *x,double *par){
  double Q=x[0];
  double thr=par[0];
  double sigma=par[1];
  double eff=par[2];
  double sig=0;

  sig=eff/(1+exp(-(Q-thr)/sigma));
  
  return sig;

}
double fitFunA(double *x,double *par){
  double A=x[0];
  double thr=par[0];
  double sigma=par[1];
  double eff=par[2];
  double sig=0;

  sig=eff/(1+exp(-(A-thr)/sigma));
  
  return sig;

}



map<pair<int,int>,int> geometry; //The map (x,y) -> ID
map<pair<int,int>,double> pTOP;  //The map (x,y) -> cal.constant K, where E=Q/K
map<pair<int,int>,double> pBOT;  //The map (x,y) -> cal.constant K, where E=Q/K
map<pair<int, int>, int>::iterator geometry_it;

typedef struct{
  int sector;
  int x;
  int y;
  int readout;
} CALO_IndexLight_t;


bool operator<(const CALO_IndexLight_t &a,const CALO_IndexLight_t &b){ //A.C. for the maps
  if (a.sector > b.sector)
    return true;
  if (a.sector < b.sector)
    return false;
  if (a.x > b.x)
    return true;
  if (a.x < b.x)
    return false;
  if (a.y > b.y)
    return true;
  if (a.y < b.y)
    return false;
  if (a.readout > b.readout)
    return true;
  if (a.readout < b.readout)
    return false;
  return false;
}



std::map<CALO_IndexLight_t,double> calo_map_thr_val,calo_map_thr_width;
std::map<CALO_IndexLight_t,double>::iterator calo_map_it;











void loadGlobal(){

  geometry[make_pair(-2, -2)] = 1;
  geometry[make_pair(-1, -2)] = 2;
  geometry[make_pair(0, -2)] = 3;
  geometry[make_pair(1, -2)] = 4;
  geometry[make_pair(2, -2)] = 5;
  geometry[make_pair(0, -1)] = 6;
  geometry[make_pair(-2, 0)] = 7;
  geometry[make_pair(-1, 0)] = 8;
  geometry[make_pair(0, 0)] = 9;
  geometry[make_pair(1, 0)] = 10;
  geometry[make_pair(2, 0)] = 11;
  geometry[make_pair(-2, 1)] = 12;
  geometry[make_pair(-1, 1)] = 13;
  geometry[make_pair(0, 1)] = 14;
  geometry[make_pair(1, 1)] = 15;
  geometry[make_pair(2, 1)] = 16;
  geometry[make_pair(-2, 2)] = 17;
  geometry[make_pair(-1, 2)] = 18;
  geometry[make_pair(1, 2)] = 19;
  geometry[make_pair(2, 2)] = 20;
  geometry[make_pair(-1, 3)] = 21;
  geometry[make_pair(1, 3)] = 22;

  pTOP[make_pair(-2, -2)] = 3.7;
  pTOP[make_pair(-1, -2)] = 3.2;
  pTOP[make_pair(0, -2)] = 3.5;
  pTOP[make_pair(1, -2)] = 3.5;
  pTOP[make_pair(2, -2)] = 3.1;
  pTOP[make_pair(0, -1)] = 3.6;
  pTOP[make_pair(-2, 0)] = 1.8;
  pTOP[make_pair(-1, 0)] = 1.7;
  pTOP[make_pair(0, 0)] = 3.7;
  pTOP[make_pair(1, 0)] = 1.8;
  pTOP[make_pair(2, 0)] = 2.0;
  pTOP[make_pair(-2, 1)] = 1.8;
  pTOP[make_pair(-1, 1)] = 1.8;
  pTOP[make_pair(0, 1)] = 3.0;
  pTOP[make_pair(1, 1)] = 2.2;
  pTOP[make_pair(2, 1)] = 1.8;
  pTOP[make_pair(-2, 2)] = 2.9;
  pTOP[make_pair(-1, 2)] = 2.3;
  pTOP[make_pair(1, 2)] = 1.9;
  pTOP[make_pair(2, 2)] = 1.9;
  pTOP[make_pair(-1, 3)] = 2.1;
  pTOP[make_pair(1, 3)] = 3.2;

  pBOT[make_pair(-2, -2)] = 3.6;
  pBOT[make_pair(-1, -2)] = 3.2;
  pBOT[make_pair(0, -2)] = 3.6;
  pBOT[make_pair(1, -2)] = 3.6;
  pBOT[make_pair(2, -2)] = 3.7;
  pBOT[make_pair(0, -1)] = 3.8;
  pBOT[make_pair(-2, 0)] = 2.6;
  pBOT[make_pair(-1, 0)] = 1.9;
  pBOT[make_pair(0, 0)] = 3.3;
  pBOT[make_pair(1, 0)] = 2.0;
  pBOT[make_pair(2, 0)] = 2.1;
  pBOT[make_pair(-2, 1)] = 3.2;
  pBOT[make_pair(-1, 1)] = 1.9;
  pBOT[make_pair(0, 1)] = 3.3;
  pBOT[make_pair(1, 1)] = 2.3;
  pBOT[make_pair(2, 1)] = 1.7;
  pBOT[make_pair(-2, 2)] = 2.2;
  pBOT[make_pair(-1, 2)] = 3.6;
  pBOT[make_pair(1, 2)] = 2.0;
  pBOT[make_pair(2, 2)] = 2.0;
  pBOT[make_pair(-1, 3)] = 3.2;
  pBOT[make_pair(1, 3)] = 1.6;

}





void doCrystalThreshold(string fname){
  string fname_simple = fname;
  fname_simple.erase(fname_simple.find_last_of("."), string::npos);

  loadGlobal();

 
  


  /*Here go the histograms*/
  const int nTOT = 22;
  TH1D *hBDXMiniCalibE_TOP[nTOT] = { 0 };
  TH1D *hBDXMiniCalibE_BOT[nTOT] = { 0 };
  TH1D *hBDXMiniCalibE_thr_TOP[nTOT] = { 0 };
  TH1D *hBDXMiniCalibE_thr_BOT[nTOT] = { 0 };

  TH1D *hBDXMiniCalibQ_TOP[nTOT] = { 0 };
  TH1D *hBDXMiniCalibQ_BOT[nTOT] = { 0 };
  TH1D *hBDXMiniCalibQ_thr_TOP[nTOT] = { 0 };
  TH1D *hBDXMiniCalibQ_thr_BOT[nTOT] = { 0 };
  
  TH1D *hBDXMiniCalibA_TOP[nTOT] = { 0 };
  TH1D *hBDXMiniCalibA_BOT[nTOT] = { 0 };
  TH1D *hBDXMiniCalibA_thr_TOP[nTOT] = { 0 };
  TH1D *hBDXMiniCalibA_thr_BOT[nTOT] = { 0 };


 
  TF1 *funTOP_E[nTOT]={0};
  TF1 *funBOT_E[nTOT]={0};
   
  TF1 *funTOP_Q[nTOT]={0};
  TF1 *funBOT_Q[nTOT]={0};
  
  TF1 *funTOP_A[nTOT]={0};
  TF1 *funBOT_A[nTOT]={0};
  

  TFile *f;
  f=new TFile(fname.c_str());

  TFile *fout=new TFile(Form("%s.CrystalThreshold.root",fname_simple.c_str()),"recreate");
  
  bool isFirst=true;

  map<pair<int, int>, int>::iterator geometry_it;
  int iX,iY,id;
  double Emin,Emax;
  double Qmin,Qmax;
  double Amin,Amax;
  

  Emin=3.;
  Emax=200.;
  
  Amin=0;

  Amax=200;
  //E=Q/p

 



 
  double s=1; 
  CALO_IndexLight_t index;

  double thrVal;

  TH1D *hThrA_val=new TH1D("hThrA_val","hThrA_val",2*geometry.size(),0.5,2*geometry.size()+0.5);
  TH1D *hThrQ_val=new TH1D("hThrQ_val","hThrQ_val",2*geometry.size(),0.5,2*geometry.size()+0.5);
  TH1D *hThrE_val=new TH1D("hThrE_val","hThrE_val",2*geometry.size(),0.5,2*geometry.size()+0.5);
  
  TH1D *hThrA_width=new TH1D("hThrA_width","hThrA_width",2*geometry.size(),0.5,2*geometry.size()+0.5);
  TH1D *hThrQ_width=new TH1D("hThrQ_width","hThrQ_width",2*geometry.size(),0.5,2*geometry.size()+0.5);
  TH1D *hThrE_width=new TH1D("hThrE_width","hThrE_width",2*geometry.size(),0.5,2*geometry.size()+0.5);
  

  for (geometry_it = geometry.begin(); geometry_it != geometry.end(); geometry_it++) {
    
    	iX = (geometry_it->first).first;
	iY = (geometry_it->first).second;
	id = (geometry_it->second);
	
	id = id - 1;

	//Get histograms
	hBDXMiniCalibE_TOP[id] = (TH1D*)f->Get(Form("hBDXMiniCalibE_s0_x%i_y%i", iX, iY));
	hBDXMiniCalibE_BOT[id] = (TH1D*)f->Get(Form("hBDXMiniCalibE_s1_x%i_y%i", iX, iY));	
	hBDXMiniCalibE_thr_TOP[id] = (TH1D*)f->Get(Form("hBDXMiniCalibE_thr_s0_x%i_y%i", iX, iY));
	hBDXMiniCalibE_thr_BOT[id] = (TH1D*)f->Get(Form("hBDXMiniCalibE_thr_s1_x%i_y%i", iX, iY));

	hBDXMiniCalibQ_TOP[id] = (TH1D*)f->Get(Form("hBDXMiniCalibQ_s0_x%i_y%i", iX, iY));
	hBDXMiniCalibQ_BOT[id] = (TH1D*)f->Get(Form("hBDXMiniCalibQ_s1_x%i_y%i", iX, iY));    
	hBDXMiniCalibQ_thr_TOP[id] = (TH1D*)f->Get(Form("hBDXMiniCalibQ_thr_s0_x%i_y%i", iX, iY));
	hBDXMiniCalibQ_thr_BOT[id] = (TH1D*)f->Get(Form("hBDXMiniCalibQ_thr_s1_x%i_y%i", iX, iY));

	hBDXMiniCalibA_TOP[id] = (TH1D*)f->Get(Form("hBDXMiniCalibA_s0_x%i_y%i", iX, iY));
	hBDXMiniCalibA_BOT[id] = (TH1D*)f->Get(Form("hBDXMiniCalibA_s1_x%i_y%i", iX, iY));    
	hBDXMiniCalibA_thr_TOP[id] = (TH1D*)f->Get(Form("hBDXMiniCalibA_thr_s0_x%i_y%i", iX, iY));
	hBDXMiniCalibA_thr_BOT[id] = (TH1D*)f->Get(Form("hBDXMiniCalibA_thr_s1_x%i_y%i", iX, iY));

	//Divide
	hBDXMiniCalibE_thr_TOP[id]->Divide(hBDXMiniCalibE_TOP[id]);
	hBDXMiniCalibE_thr_BOT[id]->Divide(hBDXMiniCalibE_BOT[id]);
	hBDXMiniCalibQ_thr_TOP[id]->Divide(hBDXMiniCalibQ_TOP[id]);
	hBDXMiniCalibQ_thr_BOT[id]->Divide(hBDXMiniCalibQ_BOT[id]);
	hBDXMiniCalibA_thr_TOP[id]->Divide(hBDXMiniCalibA_TOP[id]);
	hBDXMiniCalibA_thr_BOT[id]->Divide(hBDXMiniCalibA_BOT[id]);
	
	funTOP_E[id]=new TF1(Form("funTOP_E_x%i_y%i",iX,iY),fitFunE,Emin,Emax,3);
	funBOT_E[id]=new TF1(Form("funTOP_E_x%i_y%i",iX,iY),fitFunE,Emin,Emax,3);	

	thrVal=	hBDXMiniCalibE_thr_TOP[id]->GetBinCenter(hBDXMiniCalibE_thr_TOP[id]->GetMaximumBin());
	funTOP_E[id]->SetParameters(thrVal,2.,1.);

	thrVal=	hBDXMiniCalibE_thr_BOT[id]->GetBinCenter(hBDXMiniCalibE_thr_BOT[id]->GetMaximumBin());
	funBOT_E[id]->SetParameters(thrVal,2.,1.);

	funTOP_E[id]->SetNpx(1000);
	funBOT_E[id]->SetNpx(1000);
	hBDXMiniCalibE_thr_TOP[id]->Fit(funTOP_E[id],"R","",Emin,Emax);	
	hBDXMiniCalibE_thr_BOT[id]->Fit(funBOT_E[id],"R","",Emin,Emax);
       
	s=pTOP[make_pair(iX,iY)];
	Qmin=Emin*s;
	Qmax=Emax*s;
	funTOP_Q[id]=new TF1(Form("funTOP_Q_x%i_y%i",iX,iY),fitFunQ,Qmin,Qmax,3);
	thrVal=	hBDXMiniCalibQ_thr_TOP[id]->GetBinCenter(hBDXMiniCalibQ_thr_TOP[id]->GetMaximumBin());
	funTOP_Q[id]->SetParameters(thrVal,2.*s,1.);
	funTOP_Q[id]->SetNpx(1000);
	hBDXMiniCalibQ_thr_TOP[id]->Fit(funTOP_Q[id],"R","",Qmin,Qmax);	

	s=pBOT[make_pair(iX,iY)];
	Qmin=Emin*pBOT[make_pair(iX,iY)];
	Qmax=Emax*pBOT[make_pair(iX,iY)];
	funBOT_Q[id]=new TF1(Form("funTOP_Q_x%i_y%i",iX,iY),fitFunQ,Qmin,Qmax,3);
	thrVal=	hBDXMiniCalibQ_thr_BOT[id]->GetBinCenter(hBDXMiniCalibQ_thr_BOT[id]->GetMaximumBin());
	funBOT_Q[id]->SetParameters(thrVal,2.*s,1.);
	funBOT_Q[id]->SetNpx(1000);
	hBDXMiniCalibQ_thr_BOT[id]->Fit(funBOT_Q[id],"R","",Qmin,Qmax);


	funTOP_A[id]=new TF1(Form("funTOP_A_x%i_y%i",iX,iY),fitFunA,Amin,Amax,3);
	funBOT_A[id]=new TF1(Form("funTOP_A_x%i_y%i",iX,iY),fitFunA,Amin,Amax,3);
	thrVal=	hBDXMiniCalibA_thr_TOP[id]->GetBinCenter(hBDXMiniCalibA_thr_TOP[id]->GetMaximumBin());
	funTOP_A[id]->SetParameters(thrVal,2.,1.);
	thrVal=	hBDXMiniCalibA_thr_BOT[id]->GetBinCenter(hBDXMiniCalibA_thr_BOT[id]->GetMaximumBin());
	funBOT_A[id]->SetParameters(thrVal,2.,1.);
	funTOP_A[id]->SetNpx(1000);
	funBOT_A[id]->SetNpx(1000);
	hBDXMiniCalibA_thr_TOP[id]->Fit(funTOP_A[id],"R","",Amin,Amax);	
	hBDXMiniCalibA_thr_BOT[id]->Fit(funBOT_A[id],"R","",Amin,Amax);


	hThrA_val->SetBinContent(id+1,funTOP_A[id]->GetParameter(0));
	hThrA_val->SetBinError(id+1,funTOP_A[id]->GetParError(0));
	hThrA_val->SetBinContent(id+geometry.size()+1,funBOT_A[id]->GetParameter(0));
	hThrA_val->SetBinError(id+geometry.size()+1,funBOT_A[id]->GetParError(0));

	hThrA_width->SetBinContent(id+1,funTOP_A[id]->GetParameter(1));
	hThrA_width->SetBinError(id+1,funTOP_A[id]->GetParError(1));
	hThrA_width->SetBinContent(id+geometry.size()+1,funBOT_A[id]->GetParameter(1));
	hThrA_width->SetBinError(id+geometry.size()+1,funBOT_A[id]->GetParError(1));

	hThrQ_val->SetBinContent(id+1,funTOP_Q[id]->GetParameter(0));
	hThrQ_val->SetBinError(id+1,funTOP_Q[id]->GetParError(0));
	hThrQ_val->SetBinContent(id+geometry.size()+1,funBOT_Q[id]->GetParameter(0));
	hThrQ_val->SetBinError(id+geometry.size()+1,funBOT_Q[id]->GetParError(0));

	hThrQ_width->SetBinContent(id+1,funTOP_Q[id]->GetParameter(1));
	hThrQ_width->SetBinError(id+1,funTOP_Q[id]->GetParError(1));
	hThrQ_width->SetBinContent(id+geometry.size()+1,funBOT_Q[id]->GetParameter(1));
	hThrQ_width->SetBinError(id+geometry.size()+1,funBOT_Q[id]->GetParError(1));

	hThrE_val->SetBinContent(id+1,funTOP_E[id]->GetParameter(0));
	hThrE_val->SetBinError(id+1,funTOP_E[id]->GetParError(0));
	hThrE_val->SetBinContent(id+geometry.size()+1,funBOT_E[id]->GetParameter(0));
	hThrE_val->SetBinError(id+geometry.size()+1,funBOT_E[id]->GetParError(0));

	hThrE_width->SetBinContent(id+1,funTOP_E[id]->GetParameter(1));
	hThrE_width->SetBinError(id+1,1,funTOP_E[id]->GetParError(1));
	hThrE_width->SetBinContent(id+geometry.size()+1,funBOT_E[id]->GetParameter(1));
	hThrE_width->SetBinError(id+geometry.size()+1,funBOT_E[id]->GetParError(1));
	
	index.sector=0;
	index.x=iX;
	index.y=iY;
	index.readout=0;
	calo_map_thr_val[index]=funTOP_Q[id]->GetParameter(0);
	calo_map_thr_width[index]=funTOP_Q[id]->GetParameter(1);

	index.sector=1;
	index.x=iX;
	index.y=iY;
	index.readout=0;
	calo_map_thr_val[index]=funBOT_Q[id]->GetParameter(0);
	calo_map_thr_width[index]=funBOT_Q[id]->GetParameter(1);


  }
 

 
  TCanvas *c=new TCanvas("c","c");
  c->Divide(2,3);
  gStyle->SetOptFit(1111);
  for (geometry_it = geometry.begin(); geometry_it != geometry.end(); geometry_it++) {
    
    iX = (geometry_it->first).first;
    iY = (geometry_it->first).second;
    id = (geometry_it->second);
    
    id = id - 1;
    
    c->cd(1);
    hBDXMiniCalibE_thr_TOP[id]->Draw();
    c->cd(2);
    hBDXMiniCalibE_thr_BOT[id]->Draw();
    c->cd(3);
    hBDXMiniCalibQ_thr_TOP[id]->Draw();
    c->cd(4);
    hBDXMiniCalibQ_thr_BOT[id]->Draw();
    c->cd(5);
    hBDXMiniCalibA_thr_TOP[id]->Draw();
    c->cd(6);
    hBDXMiniCalibA_thr_BOT[id]->Draw();

  if (isFirst){
      isFirst=false;
      c->Print(Form("%s.CrystalThreshold.pdf(",fname_simple.c_str()));	    
    }
    else c->Print(Form("%s.CrystalThreshold.pdf",fname_simple.c_str()));	    
  }


  TCanvas *c2=new TCanvas("c2","c2");
  c2->Divide(2,3);
  c2->cd(1);
  hThrA_val->SetLineWidth(2);
  hThrA_val->Draw();
  c2->cd(2);
  hThrA_width->SetLineWidth(2);
  hThrA_width->Draw();
  c2->cd(3);
  hThrQ_val->SetLineWidth(2);
  hThrQ_val->Draw();
  c2->cd(4);
  hThrQ_width->SetLineWidth(2);
  hThrQ_width->Draw();
  c2->cd(5);
  hThrE_val->SetLineWidth(2);
  hThrE_val->Draw();
  c2->cd(6);
  hThrE_width->SetLineWidth(2);
  hThrE_width->Draw();

  c2->Print(Form("%s.CrystalThreshold.pdf)",fname_simple.c_str()));
  
  fout->cd();

 for (geometry_it = geometry.begin(); geometry_it != geometry.end(); geometry_it++) {
    
    iX = (geometry_it->first).first;
    iY = (geometry_it->first).second;
    id = (geometry_it->second);
    
    id = id - 1;

    hBDXMiniCalibE_thr_TOP[id]->Write();
    hBDXMiniCalibE_thr_BOT[id]->Write();
    hBDXMiniCalibQ_thr_TOP[id]->Write();
    hBDXMiniCalibQ_thr_BOT[id]->Write();
    hBDXMiniCalibA_thr_TOP[id]->Write();
    hBDXMiniCalibA_thr_BOT[id]->Write();
  }



  hThrA_val->Write();
  hThrA_width->Write();
  hThrQ_val->Write();
  hThrQ_width->Write();
  hThrE_val->Write();
  hThrE_width->Write();
  fout->Write();
  fout->Close();

  ofstream fCalib(Form("%s.CrystalThreshold.dat",fname_simple.c_str()));
  double val,width;
  fCalib<<"#sector X Y readout thrVal thrWidth"<<endl;
  for (int isector=0;isector<20;isector++){
    for (int ix=-10;ix<10;ix++){
      for (int iy=-10;iy<10;iy++){
	index.sector=isector;
	index.x=ix;
	index.y=iy;
	index.readout=0;
	
	calo_map_it=calo_map_thr_val.find(index);
	if (calo_map_it==calo_map_thr_val.end()){
	  val=0;
	  width=0;
	}
	else{ 
	  val=calo_map_it->second;
	  width=calo_map_thr_width[index];
	}
	fCalib<<isector<<" "<<ix<<" "<<iy<<" 0 "<<val<<" "<<width<<endl;
      }
    }
  }   
  

  fCalib.close();
  

}
