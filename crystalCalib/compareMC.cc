std::map<std::pair<int,int>,std::vector<std::pair<double,double> > > geometry;
std::map<std::pair<int,int>,std::vector<std::pair<double,double> > >::iterator geometry_it;

std::map<std::pair<int,int>,std::pair<double,double> >  geometry2;
std::map<std::pair<int,int>,std::pair<double,double> >::iterator  geometry2_it;


TH2Poly *hTOP,*hBOT;

void makeGeometry(){
  
  //CRS
  geometry[make_pair(0,-2)].push_back(make_pair(-7.5,-60));
  geometry[make_pair(0,-2)].push_back(make_pair(7.5,-60));
  geometry[make_pair(0,-2)].push_back(make_pair(7.5,-30));
  geometry[make_pair(0,-2)].push_back(make_pair(-7.5,-30));
  geometry2[make_pair(0,-2)]=make_pair(0,-45);

  //CRS
  geometry[make_pair(1,-2)].push_back(make_pair(7.5,-60));
  geometry[make_pair(1,-2)].push_back(make_pair(22.5,-60));
  geometry[make_pair(1,-2)].push_back(make_pair(22.5,-30));
  geometry[make_pair(1,-2)].push_back(make_pair(7.5,-30));
  geometry2[make_pair(1,-2)]=make_pair(15,-45);

  //CRS
  geometry[make_pair(-1,-2)].push_back(make_pair(-22.5,-60));
  geometry[make_pair(-1,-2)].push_back(make_pair(-7.5,-60));
  geometry[make_pair(-1,-2)].push_back(make_pair(-7.5,-30));
  geometry[make_pair(-1,-2)].push_back(make_pair(-22.5,-30));
  geometry2[make_pair(-1,-2)]=make_pair(-15,-45);

  //CRS
  geometry[make_pair(2,-2)].push_back(make_pair(22.5,-50));
  geometry[make_pair(2,-2)].push_back(make_pair(37.5,-50));
  geometry[make_pair(2,-2)].push_back(make_pair(37.5,-20));
  geometry[make_pair(2,-2)].push_back(make_pair(22.5,-20));
  geometry2[make_pair(2,-2)]=make_pair(30,-35);
  
  //CRS
  geometry[make_pair(-2,-2)].push_back(make_pair(-22.5,-50));
  geometry[make_pair(-2,-2)].push_back(make_pair(-37.5,-50));
  geometry[make_pair(-2,-2)].push_back(make_pair(-37.5,-20));
  geometry[make_pair(-2,-2)].push_back(make_pair(-22.5,-20));
  geometry2[make_pair(-2,-2)]=make_pair(-30,-35);

  //CRS
  geometry[make_pair(0,-1)].push_back(make_pair(-15,-30));
  geometry[make_pair(0,-1)].push_back(make_pair(15,-30));
  geometry[make_pair(0,-1)].push_back(make_pair(15,-15));
  geometry[make_pair(0,-1)].push_back(make_pair(-15,-15));
  geometry2[make_pair(0,-1)]=make_pair(0,-22.5);
  
  //CRS
  geometry[make_pair(0,0)].push_back(make_pair(-15,-15));
  geometry[make_pair(0,0)].push_back(make_pair(15,-15));
  geometry[make_pair(0,0)].push_back(make_pair(15,0));
  geometry[make_pair(0,0)].push_back(make_pair(-15,0));
  geometry2[make_pair(0,0)]=make_pair(0,-7.5);

  //CRS
  geometry[make_pair(0,1)].push_back(make_pair(-15,0));
  geometry[make_pair(0,1)].push_back(make_pair(15,0));
  geometry[make_pair(0,1)].push_back(make_pair(15,15));
  geometry[make_pair(0,1)].push_back(make_pair(-15,15)); 
  geometry2[make_pair(0,1)]=make_pair(0,7.5);

  //CRS
  geometry[make_pair(1,0)].push_back(make_pair(15,-20));
  geometry[make_pair(1,0)].push_back(make_pair(35,-20));
  geometry[make_pair(1,0)].push_back(make_pair(35,0));
  geometry[make_pair(1,0)].push_back(make_pair(15,0));
  geometry2[make_pair(1,0)]=make_pair(25,-10);

  //CRS
  geometry[make_pair(2,0)].push_back(make_pair(35,-20));
  geometry[make_pair(2,0)].push_back(make_pair(55,-20));
  geometry[make_pair(2,0)].push_back(make_pair(55,0));
  geometry[make_pair(2,0)].push_back(make_pair(35,0));
  geometry2[make_pair(2,0)]=make_pair(45,-10);

  //CRS
  geometry[make_pair(-1,0)].push_back(make_pair(-15,-20));
  geometry[make_pair(-1,0)].push_back(make_pair(-35,-20));
  geometry[make_pair(-1,0)].push_back(make_pair(-35,0));
  geometry[make_pair(-1,0)].push_back(make_pair(-15,0));
  geometry2[make_pair(-1,0)]=make_pair(-25,-10);
  
  //CRS
  geometry[make_pair(-2,0)].push_back(make_pair(-35,-20));
  geometry[make_pair(-2,0)].push_back(make_pair(-55,-20));
  geometry[make_pair(-2,0)].push_back(make_pair(-55,0));
  geometry[make_pair(-2,0)].push_back(make_pair(-35,0));
  geometry2[make_pair(-2,0)]=make_pair(-45,-10);

  //CRS
  geometry[make_pair(1,1)].push_back(make_pair(15,20));
  geometry[make_pair(1,1)].push_back(make_pair(35,20));
  geometry[make_pair(1,1)].push_back(make_pair(35,0));
  geometry[make_pair(1,1)].push_back(make_pair(15,0));
  geometry2[make_pair(1,1)]=make_pair(25,10);
 
  //CRS
  geometry[make_pair(2,1)].push_back(make_pair(35,20));
  geometry[make_pair(2,1)].push_back(make_pair(55,20));
  geometry[make_pair(2,1)].push_back(make_pair(55,0));
  geometry[make_pair(2,1)].push_back(make_pair(35,0));
  geometry2[make_pair(2,1)]=make_pair(45,10);

  //CRS
  geometry[make_pair(-1,1)].push_back(make_pair(-15,20));
  geometry[make_pair(-1,1)].push_back(make_pair(-35,20));
  geometry[make_pair(-1,1)].push_back(make_pair(-35,0));
  geometry[make_pair(-1,1)].push_back(make_pair(-15,0));
  geometry2[make_pair(-1,1)]=make_pair(-25,10);
 
  //CRS
  geometry[make_pair(-2,1)].push_back(make_pair(-35,20));
  geometry[make_pair(-2,1)].push_back(make_pair(-55,20));
  geometry[make_pair(-2,1)].push_back(make_pair(-55,0));
  geometry[make_pair(-2,1)].push_back(make_pair(-35,0));
  geometry2[make_pair(-2,1)]=make_pair(-45,10);

  //CRS
  geometry[make_pair(-2,2)].push_back(make_pair(-40,20));
  geometry[make_pair(-2,2)].push_back(make_pair(-20,20));
  geometry[make_pair(-2,2)].push_back(make_pair(-20,40));
  geometry[make_pair(-2,2)].push_back(make_pair(-40,40));
  geometry2[make_pair(-2,2)]=make_pair(-30,30);

  //CRS
  geometry[make_pair(-1,2)].push_back(make_pair(-20,20));
  geometry[make_pair(-1,2)].push_back(make_pair(0,20));
  geometry[make_pair(-1,2)].push_back(make_pair(0,40));
  geometry[make_pair(-1,2)].push_back(make_pair(-20,40));
  geometry2[make_pair(-1,2)]=make_pair(-10,30);

  //CRS
  geometry[make_pair(1,2)].push_back(make_pair(20,20));
  geometry[make_pair(1,2)].push_back(make_pair(0,20));
  geometry[make_pair(1,2)].push_back(make_pair(0,40));
  geometry[make_pair(1,2)].push_back(make_pair(20,40));
  geometry2[make_pair(1,2)]=make_pair(10,30);
 
  //CRS
  geometry[make_pair(2,2)].push_back(make_pair(20,20));
  geometry[make_pair(2,2)].push_back(make_pair(40,20));
  geometry[make_pair(2,2)].push_back(make_pair(40,40));
  geometry[make_pair(2,2)].push_back(make_pair(20,40));
  geometry2[make_pair(2,2)]=make_pair(30,30);
 
  //CRS
  geometry[make_pair(-1,3)].push_back(make_pair(-20,40));
  geometry[make_pair(-1,3)].push_back(make_pair(0,40));
  geometry[make_pair(-1,3)].push_back(make_pair(0,60));
  geometry[make_pair(-1,3)].push_back(make_pair(-20,60));
  geometry2[make_pair(-1,3)]=make_pair(-10,50);
 
  //CRS
  geometry[make_pair(1,3)].push_back(make_pair(20,40));
  geometry[make_pair(1,3)].push_back(make_pair(0,40));
  geometry[make_pair(1,3)].push_back(make_pair(0,60));
  geometry[make_pair(1,3)].push_back(make_pair(20,60));
  geometry2[make_pair(1,3)]=make_pair(10,50);
  
  hTOP=new TH2Poly();
  hTOP->SetName("hTOP");
  hBOT=new TH2Poly();
  hBOT->SetName("hBOT");
  int N;
  double *x;
  double *y;
  
  for (geometry_it = geometry.begin();geometry_it !=geometry.end();geometry_it++){
    N=geometry_it->second.size();
    x=new double[N];
    y=new double[N];
    for (int ii=0;ii<N;ii++){
      x[ii]=(geometry_it->second)[ii].first;
      y[ii]=(geometry_it->second)[ii].second;
    }
    hTOP->AddBin(N,x,y);
    hBOT->AddBin(N,x,y);    
  }

  

}

void compareMC(string fname){

  makeGeometry();

  string fname_simple = fname;
  fname_simple.erase(fname_simple.find_last_of("."), string::npos);
  
  TFile *fout=new TFile(Form("%s.compareMC.root",fname_simple.c_str()),"recreate");


  
  TFile *fData=new TFile(fname.c_str());
  TFile *fMC=new TFile("MC.root");
  TTree *runInfo=(TTree*)fData->Get("RunInfo");

  int runN;
  Long64_t dT;
  runInfo->SetBranchAddress("runN",&runN);
  runInfo->SetBranchAddress("dT",&dT);
  runInfo->GetEntry(0);
 
  
   /*Now, start looping over all objects*/
  TIter next(fData->GetListOfKeys());
  TKey *key;
  TH1D *hData,*hMC;

 
  TCanvas *c=new TCanvas("cout");
  c->Print(Form("%s.compareMC.pdf(",fname_simple.c_str()));
  

  double integralData,integralMC;

  int iX,iY;
  for (geometry2_it = geometry2.begin();geometry2_it !=geometry2.end();geometry2_it++){
    iX=geometry2_it->first.first;
    iY=geometry2_it->first.second;

    hData=(TH1D*)fData->Get(Form("hBDXMiniCalibE_s0_x%i_y%i",iX,iY));
    hMC=(TH1D*)fMC->Get(hData->GetName());
    hData->Rebin(4);
    hData->Scale(1./dT,"width");
      
    c->cd();
    hData->Draw("HIST");
    hData->GetXaxis()->SetRangeUser(10,200);
    hMC->SetLineColor(2);
    hMC->SetMarkerColor(2);
    hMC->Draw("HISTSAME");
    
    c->Print(Form("%s.compareMC.pdf",fname_simple.c_str()));
  
    integralData=hData->Integral(hData->FindBin(20.),hData->FindBin(150.),"width");
    integralMC=hMC->Integral(hMC->FindBin(20.),hMC->FindBin(150.),"width");

    cout<<"TOP "<<iX<<" "<<iY<<" "<<integralData<<" "<<integralMC<<endl;
    hTOP->Fill(geometry2_it->second.first,geometry2_it->second.second,integralData/integralMC);
    

    hData=(TH1D*)fData->Get(Form("hBDXMiniCalibE_s1_x%i_y%i",iX,iY));
    hMC=(TH1D*)fMC->Get(hData->GetName());
    hData->Rebin(4);
    hData->Scale(1./dT,"width");
      
    c->cd();
    hData->Draw("HIST");
    hData->GetXaxis()->SetRangeUser(10,200);
    hMC->SetLineColor(2);
    hMC->SetMarkerColor(2);
    hMC->Draw("HISTSAME");
     c->Print(Form("%s.compareMC.pdf",fname_simple.c_str()));

 

    integralData=hData->Integral(hData->FindBin(20.),hData->FindBin(150.),"width");
    integralMC=hMC->Integral(hMC->FindBin(20.),hMC->FindBin(150.),"width");
    cout<<"BOT "<<iX<<" "<<iY<<" "<<integralData<<" "<<integralMC<<endl;
    hBOT->Fill(geometry2_it->second.first,geometry2_it->second.second,integralData/integralMC);






  }





  TCanvas *c2=new TCanvas("c2","c2");
  c2->cd();
  c2->Divide(2,1);
  c2->cd(1);
  hTOP->Draw("colz");
  c2->cd(2);
  hBOT->Draw("colz");
  c->Print(Form("%s.compareMC.pdf)",fname_simple.c_str()));

  fout->cd();
  hTOP->Write();
  hBOT->Write();
  fout->Write();
  fout->Close();

}
