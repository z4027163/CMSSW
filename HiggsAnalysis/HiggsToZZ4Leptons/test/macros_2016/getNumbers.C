{


  //TFile *_file0 = TFile::Open("plots/htotalfinal_hM4l_7_4mu.root");
  //TFile *_file0 = TFile::Open("plots/htotalfinal_hM4l_7_4e.root");
  //TFile *_file0 = TFile::Open("plots/htotalfinal_hM4l_7_2e2mu.root");
  TFile *_file0 = TFile::Open("plots/htotalfinal_hPFMET_3_2e2mu.root");  

  //TH1F *Histo = (TH1F*)_file0->Get("htotaldata");
  TH1F *Histo = (TH1F*)_file0->Get("htotalHisto");
  
  TAxis *axis = Histo->GetXaxis();

  float xmin=0.,xmax=10000.;
 //float xmin=118.,xmax=130.;
 //float xmin=70.,xmax=10000.;
 //float xmin=400.,xmax=10000.;
  //float xmin=70.,xmax=10000.;

 int bmin = axis->FindBin(xmin); 
 int bmax = axis->FindBin(xmax); 
 double integral = Histo->Integral(bmin,bmax);
 cout << "Integral before= " << integral << endl;
 integral-=Histo->GetBinContent(bmin)*(xmin-(axis->GetBinLowEdge(bmin)))/(axis->GetBinWidth(bmin));
 integral-=Histo->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/(axis->GetBinWidth(bmax));
 cout << "Integral after= " << integral << endl;

 double entry=0;
 for (int i=bmin;i<bmax;i++){
   entry+=Histo->GetBinContent(i)*Histo->GetEntries()/Histo->Integral();
   // cout << "Entry= " << entry << endl;
 }
 double error=sqrt(entry)*Histo->Integral()/Histo->GetEntries();
 cout << "Entry final= "<< entry << " " << integral << " weighted statistical error= " << error << " weighted systematic error (20%)= " << integral*0.2<< endl;
}
