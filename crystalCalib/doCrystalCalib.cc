#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"

#include "RooRealVar.h"

using namespace std;
using namespace RooFit;


//Global variables
map<pair<int,int>,int> geometry; //The map (x,y) -> ID
map<pair<int,int>,double> pTOP;  //The map (x,y) -> cal.constant K, where E=Q/K
map<pair<int,int>,double> pBOT;  //The map (x,y) -> cal.constant K, where E=Q/K
map<pair<int, int>, int>::iterator geometry_it;

int isFirst=1;




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



std::map<CALO_IndexLight_t,double> calo_map;
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

//s0 is the initial scale parameter in the relation E=Q/s0  -> Q=E*s0
double doFit(TFile *fData,TFile *fMC,int sector,int iX,int iY,double Emin,double Emax){

  //Load histograms
  TH1D *hData=(TH1D*)fData->Get(Form("hBDXMiniCalibQ_s%i_x%i_y%i",sector, iX, iY));
  TH1D *hMC=(TH1D*)fMC->Get(Form("hBDXMiniCalibE_s%i_x%i_y%i",sector, iX, iY));

  hMC->Smooth();
  hMC->Rebin(4);
  hData->Rebin(2);
  
  
  
  //Get the initial parameter
  double s0;
  if (sector==0){
    s0=pTOP[make_pair(iX,iY)];
  }else{
    s0=pBOT[make_pair(iX,iY)];
  }
 
  //Load the data
  RooRealVar Q("Q","Q",Emin*s0,Emax*s0);
  RooDataHist data("data","data",Q,hData);


  
  

  //Import the PDF from the hMC
  RooRealVar E("E","E",0,500);
  RooDataHist histdata("histdata","histdata",E,hMC);
 
  s0=1./s0;
  RooRealVar scale("scale","scale",s0,0,1);  //ODG: s~2..3, 1/s is within 0 and 1
  RooRealVar p0("p0","p0",0.);
  RooPolyVar Qf("Qf","Qf",Q,RooArgSet(p0,scale)); //E=Q/s -> Q=E*s
  
  RooHistPdf histpdf("histpdf","histpdf",Qf,E,histdata,0);
  
  //Define bg function
  RooPolynomial pol0("pol0","pol0",Qf,RooArgList());

  RooRealVar f0("f0","f0",1.,0.8,1.0);
  
  RooAddPdf model("model","model",RooArgList(histpdf,f0),RooArgList(f0));


  
  //do the fit
  model.fitTo(data);


  //Plotting
  RooRealVar QPlot("QPlot","QPlot",0,800);
  RooDataHist dataPlot("dataPlot","data",QPlot,hData);
  
  RooPlot* frame = Q.frame();
  data.plotOn(frame,MarkerColor(kRed));
  histpdf.plotOn(frame);

  RooPlot* frame2 = QPlot.frame();
  dataPlot.plotOn(frame2);
  

  

  TCanvas *c=new TCanvas("c","c");
  frame2->Draw();
  frame->Draw("SAME");

  if (isFirst){
    c->Print(Form("%s.calib.pdf(",fData->GetName()));
    isFirst=0;
  }else{
    c->Print(Form("%s.calib.pdf",fData->GetName()));
  }

  s0=scale.getValV();
  s0=1/s0;

  c->Modified();
  c->Update();

  
  //delete c;

  return s0;
  
  
   
}




void doCrystalCalib(string fname){


  TFile *fData=new TFile(fname.c_str()); 
  TFile *fMC=new TFile("MC.root");

  double Emin=14;
  double Emax=200;

  TH1D *hData;
  TH1D *hMC;

  int iX,iY,id;

  loadGlobal();

  double ret;
  CALO_IndexLight_t index;

  //  doFit(fData,fMC,0,-2,-2,Emin,Emax);

  
  
  for (geometry_it = geometry.begin(); geometry_it != geometry.end(); geometry_it++) {
    
    	iX = (geometry_it->first).first;
	iY = (geometry_it->first).second;
	id = (geometry_it->second);
	
	id = id - 1;

	index.sector=0;
	index.x=iX;
	index.y=iY;
	index.readout=0;

	

	
	ret=doFit(fData,fMC,0,iX,iY,Emin,Emax);
	cout<<" 0 "<<iX<<" "<<iY<<" "<<ret<<" "<<pTOP[make_pair(iX,iY)]<<endl;


	calo_map[index]=ret;
	

  }

  for (geometry_it = geometry.begin(); geometry_it != geometry.end(); geometry_it++) {
    
    	iX = (geometry_it->first).first;
	iY = (geometry_it->first).second;
	id = (geometry_it->second);
	id = id - 1;

	index.sector=1;
	index.x=iX;
	index.y=iY;
	index.readout=0;
	
	ret=doFit(fData,fMC,1,iX,iY,Emin,Emax);
	cout<<" 1 "<<iX<<" "<<iY<<" "<<ret<<" "<<pBOT[make_pair(iX,iY)]<<endl;

	calo_map[index]=ret;
	}
  
  TCanvas *cFake=new TCanvas();
  cFake->Print(Form("%s.calib.pdf)",fData->GetName()));

  /*Now write*/
  ofstream ocalib(Form("%s.calib",fname.c_str()));
  ocalib<<"#sector X Y readout calib offset"<<endl;

  for (int isector=0;isector<20;isector++){
    for (int ix=-10;ix<10;ix++){
      for (int iy=-10;iy<10;iy++){
	index.sector=isector;
	index.x=ix;
	index.y=iy;
	index.readout=0;

	calo_map_it=calo_map.find(index);
	if (calo_map_it==calo_map.end())ret=0;
	else ret=calo_map_it->second;

	ocalib<<isector<<" "<<ix<<" "<<iy<<" 0 "<<ret<<" 0 "<<endl;
       
      }
    }
  }
  ocalib.close();
  
  
}
  
 
 

  



