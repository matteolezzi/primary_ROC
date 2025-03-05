#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <iterator>



#include <math.h>
#include <time.h>
#include <stdlib.h>
//#include <sys/types.h>
#include <unistd.h>
//#include <sys/time.h>






// ROOT



#include <TMath.h>
#include <TObject.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TRandom3.h>



using namespace std;

double separazione(){
	//creation of file
	TFile *pr = new TFile("C:/Users/Matteo Lezzi/Desktop/Matteo/Universita/Magistrale/laboratorio_di_analisi_dati/separazione massa/conex_eposlhc_1301560476_100.root");
	TFile *ir = new TFile("C:/Users/Matteo Lezzi/Desktop/Matteo/Universita/Magistrale/laboratorio_di_analisi_dati/separazione massa/conex_eposlhc_625923668_5600.root");
	TFile *ph = new TFile("C:/Users/Matteo Lezzi/Desktop/Matteo/Universita/Magistrale/laboratorio_di_analisi_dati/separazione massa/conex_eposlhc_2086832343_0.root");

	TFile *pr_ene = new TFile("C:/Users/Matteo Lezzi/Desktop/Matteo/Universita/Magistrale/laboratorio_di_analisi_dati/separazione massa/protoni_10000.root");

	int bins=300;
	// define tree 
	TTree *ShowerPR = (TTree*)pr->Get("Shower");
	TTree *ShowerIR = (TTree*)ir->Get("Shower");
	TTree *ShowerPH = (TTree*)ph->Get("Shower");
	
	//create graph
	TH1F *hpr = new TH1F("hpr","hpr",bins,500,2000);
	TH1F *hir = new TH1F("hir","hir",bins,500,2000);
	TH1F *hph= new TH1F("hph","hph",bins,500,2000);
	
	TH1F *hpr1 = new TH1F("hpr1","hpr1;x_max",bins,500,2000);
	TH1F *hir1 = new TH1F("hir1","hir1;x_max",bins,500,2000);
	TH1F *hph1 = new TH1F("hph1","hph1;x_max",bins,500,2000);
	
	TH1F *hpr2 = new TH1F("hpr2","hpr2;x_max",bins,500,2000);
	TH1F *hir2 = new TH1F("hir2","hir2;x_max",bins,500,2000);
	TH1F *hph2 = new TH1F("hph2","hph2;x_max",bins,500,2000);
	
	
	#define variables
	float x_max_ir=0.;
	float x_max_pr=0.;
	float x_max_ph=0.;
	float x_max_ir1=0.;
	float x_max_pr1=0.;
	float x_max_ph1=0.;
	
	float x_max_ir2=0.;
	float x_max_pr2=0.;
	float x_max_ph2=0.;
	
	float sigma1 = 20.0;
	float sigma2 = 40.0;
	
	//define tree branches
	ShowerPR->SetBranchAddress("Xmax", &x_max_pr);
	ShowerIR->SetBranchAddress("Xmax", &x_max_ir);
	ShowerPH->SetBranchAddress("Xmax", &x_max_ph);
	//load values from branches
	for(int i=0; i<ShowerPR->GetEntries(); i++){
		ShowerPR->GetEntry(i);
		ShowerIR->GetEntry(i);
		ShowerPH->GetEntry(i);
		x_max_ir1 = gRandom->Gaus(x_max_ir,sigma1);
		x_max_pr1 = gRandom->Gaus(x_max_pr,sigma1);
		x_max_ph1 = gRandom->Gaus(x_max_ph,sigma1);
		
		x_max_ir2 = gRandom->Gaus(x_max_ir,sigma2);
		x_max_pr2 = gRandom->Gaus(x_max_pr,sigma2);
		x_max_ph2 = gRandom->Gaus(x_max_ph,sigma2);
	
		hpr->Fill(x_max_pr);
		hir->Fill(x_max_ir);
		hph->Fill(x_max_ph);
		
		hpr1->Fill(x_max_pr1);
		hir1->Fill(x_max_ir1);
		hph1->Fill(x_max_ph1);
		
		hpr2->Fill(x_max_pr2);
		hir2->Fill(x_max_ir2);
		hph2->Fill(x_max_ph2);
		
		
		}
	
	TCanvas *separazione = new TCanvas("separazione","separazione",1000,600);
	
	hpr->Draw();
	hir->Draw("same");
	hph->Draw("same");
	
	hpr1->SetLineColor(1);
	hir1->SetLineColor(1);
	hph1->SetLineColor(1);
	
	hpr1->Draw("same");
	hir1->Draw("same");
	hph1->Draw("same");
	
	hpr2->SetLineColor(2);
	hir2->SetLineColor(2);
	hph2->SetLineColor(2);
	hpr2->Draw("same");
	hir2->Draw("same");
	hph2->Draw("same");
	
	
	float MF_prot_fe=0.; //Merit Factor
	float MF_prot_gamma=0.;
    float MF_fe_gamma=0.;
    
	float MF_prot_fe_1=0.; //Merit Factor
	float MF_prot_gamma_1=0.;
    float MF_fe_gamma_1=0.;
    
	float MF_prot_fe_2=0.; //Merit Factor
	float MF_prot_gamma_2=0.;
    float MF_fe_gamma_2=0.;
    
	
	float mean_fe=0.;
    float mean_pr=0.;
    float mean_ph=0.;
    
    float mean_fe_1=0.;
    float mean_pr_1=0.;
    float mean_ph_1=0.;
    
    float mean_fe_2=0.;
    float mean_pr_2=0.;
    float mean_ph_2=0.;
    
    float rms_pr=0.;
    float rms_fe=0.;
    float rms_ph=0.;
    float rms_pr_1=0.;
    float rms_fe_1=0.;
    float rms_ph_1=0.;
    float rms_pr_2=0.;
    float rms_fe_2=0.;
    float rms_ph_2=0.;
    
    
    
	mean_pr = hpr->GetMean();
    mean_fe = hir->GetMean();
    mean_ph = hph->GetMean();
    
	
	mean_pr_1 = hpr1->GetMean();
    mean_fe_1 = hir1->GetMean();
    mean_ph_1 = hph1->GetMean();
	
	mean_pr_2 = hpr2->GetMean();
    mean_fe_2 = hir2->GetMean();
    mean_ph_2 = hph2->GetMean();
	

	rms_pr = hpr->GetRMS();
    rms_fe = hir->GetRMS();
    rms_ph = hph->GetRMS();
    
    rms_pr_1 = hpr1->GetRMS();
    rms_fe_1 = hir1->GetRMS();
    rms_ph_1 = hph1->GetRMS();
    
    rms_pr_2 = hpr2->GetRMS();
    rms_fe_2 = hir2->GetRMS();
    rms_ph_2 = hph2->GetRMS();
    //calculate cumulative efficiency
    double * eff_ir = hir->GetIntegral(); 
    double * eff_pr = hpr->GetIntegral(); 
    double * eff_ph = hph->GetIntegral();
    
    double * eff_ir_2 = hir2->GetIntegral();
    double * eff_pr_2 = hpr2->GetIntegral();
    double * eff_ph_2 = hph2->GetIntegral();
    double  inv_eff_pr[300] = {0.};
    double  inv_eff_ph[300] = {0.};
	double  inv_eff_pr_2[300] = {0.};
    double  inv_eff_ph_2[300] = {0.};
	//calculate inverse cumulative efficency for inverse ROC curve
	for(int i=0; i<bins;i++){
    inv_eff_pr[i] = 1. - eff_pr[i];
	inv_eff_ph[i] = 1. - eff_ph[i];	
    
	inv_eff_pr_2[i] = 1. - eff_pr_2[i];
	inv_eff_ph_2[i] = 1. - eff_ph_2[i];	
	}
    
    double xq[3],yq[3];
    xq[0] = 0.16;
    xq[1] = 0.50;
    xq[2] = 0.84;
    hph->GetQuantiles(3,yq,xq);
    double x_max_16 = yq[2];
    double x_max_50 = yq[1];
    double x_max_84 = yq[0];
    
    cout<<"cut fotoni per efficienza 16% "<<x_max_16<<endl;
    cout<<"cut fotoni per efficienza 50% "<<x_max_50<<endl;
    cout<<"cut fotoni per efficienza 84% "<<x_max_84<<endl;
    
    hph2->GetQuantiles(3,yq,xq);
    double x_max_16_2 = yq[2];
    double x_max_50_2 = yq[1];
    double x_max_84_2 = yq[0];
    
    cout<<"cut fotoni per efficienza 16% reale "<<x_max_16_2<<endl;
    cout<<"cut fotoni per efficienza 50% reale "<<x_max_50_2<<endl;
    cout<<"cut fotoni per efficienza 84% reale "<<x_max_84_2<<endl;
    

    
    TGraph * cont_pr_fe = new TGraph(bins,eff_pr,eff_ir);
    TGraph * cont_ph_pr = new TGraph(bins,eff_ph,eff_pr);
    TGraph * cont_ph_fe = new TGraph(bins,eff_ph,eff_ir);
    
    
    TGraph * cont_pr_fe_2 = new TGraph(bins,eff_pr_2,eff_ir_2);
    TGraph * cont_ph_pr_2 = new TGraph(bins,eff_ph_2,eff_pr_2);
    TGraph * cont_ph_fe_2 = new TGraph(bins,eff_ph_2,eff_ir_2);
    
    TGraph * inv_cont_pr_ph = new TGraph(300,inv_eff_pr,inv_eff_ph);
    
    TGraph * inv_cont_pr_ph_2 = new TGraph(bins,inv_eff_pr_2,inv_eff_ph_2);
    
    float MuPR[300] ={0.};
    float MuIR[300] ={0.};
    float MuPH[300] ={0.};
    //5000 events
	float MuPR_ground[5000]={0.};
    float MuIR_ground[5000]={0.};
    float MuPH_ground[5000]={0.};
    float MuPR_ground_nolog[5000]={0.};
    
    float XPR_max[5000]={0.};
    float XIR_max[5000]={0.};
    float XPH_max[5000]={0.};
	
    float XPR_max[5000]={0.};
    float XIR_max[5000]={0.};
    float XPH_max[5000]={0.};
    
    float zenithPR=0;
    float zenithIR=0;
    float zenithPH=0;
    
	   
    int nXPR=0;
    int nXIR=0;
    int nXPH=0;
    //define branches
    ShowerPR->SetBranchAddress("Mu",&MuPR);
	ShowerIR->SetBranchAddress("Mu",&MuIR);
	ShowerPH->SetBranchAddress("Mu",&MuPH);
	
	ShowerPR->SetBranchAddress("nX",&nXPR);
	ShowerIR->SetBranchAddress("nX",&nXIR);
	ShowerPH->SetBranchAddress("nX",&nXPH);
	
	ShowerPR->SetBranchAddress("zenith",&zenithPR);
	ShowerIR->SetBranchAddress("zenith",&zenithIR);
	ShowerPH->SetBranchAddress("zenith",&zenithPH);
	
	
	
    gpr->SetMarkerColor(2);
    gpr->SetMarkerStyle(20);

    gir->SetMarkerColor(1);
    gir->SetMarkerStyle(20);
    
    gph->SetMarkerColor(4);
    gph->SetMarkerStyle(20);
    
    gph->Draw("AP");
    gpr->Draw("Psame");
    gir->Draw("Psame");
    
    
    
    MF_prot_fe = (mean_pr - mean_fe)/sqrt( pow(rms_pr,2) + pow(rms_fe,2)) ;
    MF_prot_gamma = -1*(mean_pr - mean_ph)/sqrt( pow(rms_pr,2) + pow(rms_ph,2)) ;
    MF_fe_gamma = -1*(mean_fe - mean_ph)/sqrt( pow(rms_fe,2) + pow(rms_ph,2)) ;
	
	MF_prot_fe_1 = (mean_pr_1 - mean_fe_1)/sqrt( pow(rms_pr_1,2) + pow(rms_fe_1,2)) ;
    MF_prot_gamma_1 = -1*(mean_pr_1 - mean_ph_1)/sqrt( pow(rms_pr_1,2) + pow(rms_ph_1,2)) ;
    MF_fe_gamma_1 = -1*(mean_fe_1 - mean_ph_1)/sqrt( pow(rms_fe_1,2) + pow(rms_ph_1,2)) ;
	
	MF_prot_fe_2 = (mean_pr_2 - mean_fe_2)/sqrt( pow(rms_pr_2,2) + pow(rms_fe_2,2)) ;
    MF_prot_gamma_2 = -1*(mean_pr_2 - mean_ph_2)/sqrt( pow(rms_pr_2,2) + pow(rms_ph_2,2)) ;
    MF_fe_gamma_2 = -1*(mean_fe_2 - mean_ph_2)/sqrt( pow(rms_fe_2,2) + pow(rms_ph_2,2)) ;
	
	
	
	cout<<"the merit factor between proton and iron is "<<MF_prot_fe<<endl;    
    cout<<"the merit factor between proton and photon "<<MF_prot_gamma<<endl;    
	cout<<"the merit factor between iron and photon is "<<MF_fe_gamma<<endl; 

	cout<<"the merit factor between proton and iron with smeared data (20 g/cm^2) is "<<MF_prot_fe_1<<endl;    
    cout<<"the merit factor between proton and photon with smeared data (20 g/cm^2) is "<<MF_prot_gamma_1<<endl;    
	cout<<"the merit factor between iron and iron photon smeared data (20 g/cm^2) is "<<MF_fe_gamma_1<<endl; 

	cout<<"the merit factor between proton and iron with smeared data (40 g/cm^2) is "<<MF_prot_fe_2<<endl;    
    cout<<"the merit factor between proton and photon with smeared data (40 g/cm^2) is  "<<MF_prot_gamma_2<<endl;    
	cout<<"the merit factor between iron and photon with smeared data (20 g/cm^2) is "<<MF_fe_gamma_2<<endl; 
	
	TCanvas *contamination = new TCanvas("contamination","contamination",1000,600);
	cont_pr_fe->SetTitle("ROC Contamination; Eff X ; Eff Y");
	cont_pr_fe->Draw();
	cont_ph_pr->SetLineColor(2);
	cont_ph_pr->Draw("same");
	cont_ph_fe->SetLineColor(3);
	cont_ph_fe->Draw("same");	
	
	TCanvas *contamination_40 = new TCanvas("contamination 40 ","contamination 40",1000,600);
	cont_pr_fe_2->Draw();
	cont_ph_pr_2->SetLineColor(2);
	cont_ph_pr_2->Draw("same");
	cont_ph_fe_2->SetLineColor(3);
	cont_ph_fe_2->Draw("same");	
	TCanvas *contamination_confront = new TCanvas("contamination confront","contamination confront",1000,600);
	cont_pr_fe_2->SetTitle("ROC Contaminazione Confron proton in iron; Eff X ; Eff Y");
	cont_pr_fe_2->SetLineStyle(2);
	cont_pr_fe_2->Draw();
	cont_pr_fe->Draw("same");
	
	TCanvas *inv_contamination_confront = new TCanvas("contamination confront inv","contamination confront inv",1000,600);
	inv_cont_pr_ph->Draw();
	inv_cont_pr_ph_2->SetLineStyle(2);
	inv_cont_pr_ph_2->Draw("same");
	
	
    TGraph * eval_cont_ph_pr_2 = new TGraph(bins,inv_eff_ph_2,inv_eff_pr_2);
    cout<<"real contamination 16%: "<<eval_cont_ph_pr_2->Eval(0.16)<<endl;
    cout<<"real contamination 50%: "<<eval_cont_ph_pr_2->Eval(0.50)<<endl;
	cout<<"real contamination 84%: "<<eval_cont_ph_pr_2->Eval(0.84)<<endl;

	TGraph * eval_cont_ph_pr = new TGraph(bins,inv_eff_ph,inv_eff_pr);
    cout<<"ideal contamination 16%: "<<eval_cont_ph_pr->Eval(0.16)<<endl;
    cout<<"ideal contamination 50%: "<<eval_cont_ph_pr->Eval(0.50)<<endl;
	cout<<"ideal contamination 84%: "<<eval_cont_ph_pr->Eval(0.84)<<endl;	
	
	
}
