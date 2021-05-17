#include "Riostream.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

const int nmult = 13;
const double multsx[] = {2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,14.5,19.5,29.5,39.5,49.5}; //range sx della molteplicità
const double multdx[] = {3.5,4.5,5.5,6.5,7.5,8.5,9.5,14.5,19.5,29.5,39.5,49.5,59.5}; //range dx
const double hrange[] = {1000,900,900,900,700,700,700,500,500,400,400,300,300}; //range istogrammi residui

//const int nmult = 19;
//const double multsx[] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,14.5,19.5,29.5,39.5,49.5,59.5,69.5,79.5,89.5}; //range sx della molteplicità
//const double multdx[] = {1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,14.5,19.5,29.5,39.5,49.5,59.5,69.5,79.5,89.5,99.5}; //range dx
//const double hrange[] = {1000,1000,1000,900,900,900,700,700,700,500,500,400,400,300,300,300,300,300,300}; //range istogrammi residui

const int nvert = 12;
const double vertsx[]={-200,-150,-110,-80,-50,-30,0,30,50,80,110,150}; //range dei residui nel vertice
const double vertdx[]={-150,-110,-80,-50,-30,0,30,50,80,110,150,200};

//const int nvert = 14;
//const double vertsx[]={-300,-200,-150,-110,-80,-50,-30,0,30,50,80,110,150,200}; //range dei residui nel vertice
//const double vertdx[]={-200,-150,-110,-80,-50,-30,0,30,50,80,110,150,200,300};

void treeAnalysis(TString drawopt = "on") {

    ifstream input("histoAnalysis.root");
    if(!input) {
        cout<<" Il file histoAnalysis.root NON ESISTE \n";
        return; //esce dal void
    }
    TFile histoan("histoAnalysis.root"); //file in lettura
    
    typedef struct {
        double zrec;
        double ztrue;
        int mult;} 
    VRec; //struttura con parametri della ricostruzione
    static VRec recvert;

    TTree *tree = (TTree*)histoan.Get("newtree"); //tree in lettura
    TBranch *b1 = tree->GetBranch("RecVert");
    b1->SetAddress(&recvert.zrec); 
    int kEntries = tree->GetEntries(); 

    TFile* histdraw = new TFile("histdraw.root","recreate"); 

    //residui
    TH1F* hrestot = new TH1F("hrestot","Residui a qualsiasi Molteplicita' ",100,-500,500); //residui a tutte le molteplicità
    hrestot->GetXaxis()->SetTitle("Z_{rec} - Z_{true}  [#mum]");
    hrestot->GetYaxis()->SetTitle(" Eventi ");

    char htitle[50]; 
    char hname[50];
    TH1F* hresm[nmult]; //residui divisi in molteplicità
    TH1F* hresv[nvert]; //residui divisi nel vertice

    for(int i=0;i<nmult;i++) { //riempiamo l'array di residui in molteplicita
        sprintf(hname,"hres_mult%d",i); //nome -> hres_mult0,hres_mult1,...
        sprintf(htitle," Residui da Molteplicita' %f a %f ",multsx[i],multdx[i]);
        hresm[i] = new TH1F(hname,htitle,hrange[i]/10,-hrange[i],hrange[i]);
        hresm[i]->GetXaxis()->SetTitle("Z_{rec} - Z_{true}  [#mum]");
        hresm[i]->GetYaxis()->SetTitle(" Eventi ");
    }

    for(int i=0;i<nvert;i++) { //riempiamo l'array di residui nel vertice
        sprintf(hname,"hres_vert%d",i); //nome -> hres_vert0,hres_vert1,...
        sprintf(htitle," Residui da Vertice %f a %f ",vertsx[i],vertdx[i]);
        hresv[i] = new TH1F(hname,htitle,100,-500,500);
        hresv[i]->GetXaxis()->SetTitle("Z_{rec} - Z_{true}  [#mum]");
        hresv[i]->GetYaxis()->SetTitle(" Eventi ");
    }

    //risoluzione nella molteplicità e nel vertice
    double binm[nmult+1];
    double binv[nvert+1];
    for(int i=0;i<nmult;i++) binm[i] = multsx[i];
    binm[nmult] = multdx[nmult-1];
    for(int i=0;i<nvert;i++) binv[i] = vertsx[i];
    binv[nvert] = vertdx[nvert-1];
    TH1F* hrism = new TH1F("hrism"," Risoluzione vs Molteplicita' ",nmult,binm);
    TH1F* hrisv = new TH1F("hrisv"," Risoluzione vs Vertice ",nvert,binv);
    hrisv->GetXaxis()->SetTitle(" Z_{true} [mm] ");
    hrisv->GetYaxis()->SetTitle(" Risoluzione #sigma [#mum] ");
    hrism->GetXaxis()->SetTitle(" Molteplicita' ");
    hrism->GetYaxis()->SetTitle(" Risoluzione #sigma [#mum] ");

    //efficienza nella molteplicità e vertice
    TH1F* heffm = new TH1F("heffm"," Efficienza vs Molteplicita' ",nmult,binm);
    TH1F* hefftotm = new TH1F("hefftotm"," Denominatore Efficienza in Molteplicita' ",nmult,binm);
    TH1F* heffv = new TH1F("heffv"," Efficienza vs Vertice ",nvert,binv);
    TH1F* hefftotv = new TH1F("hefftotv"," Denominatore Efficienza nel Vertice ",nvert,binv);
    heffv->GetXaxis()->SetTitle(" Z_{true} [mm] ");
    heffv->GetYaxis()->SetTitle(" Efficienza ");
    heffm->GetXaxis()->SetTitle(" Molteplicita' ");
    heffm->GetYaxis()->SetTitle(" Efficienza ");

    double dz;

    for(int i=0;i<kEntries;i++) {
        tree->GetEvent(i);
        hefftotm->Fill(recvert.mult);
        hefftotv->Fill(recvert.ztrue);
        if (recvert.zrec == -500) continue;
        else {
            heffm->Fill(recvert.mult);
            heffv->Fill(recvert.ztrue);
            dz = recvert.zrec - recvert.ztrue;
            dz = dz*1000; //in micrometri
            hrestot->Fill(dz); //residui totali in ogni molteplicità
            for(int j=0;j<nmult;j++) { //residui in molteplicità
                if (recvert.mult>=multsx[j] && recvert.mult<=multdx[j]) hresm[j]->Fill(dz); //residuo corrispondente in molteplicita
            }
            for(int j=0;j<nvert;j++) { //residui in vertice
                if (recvert.ztrue>=vertsx[j] && recvert.ztrue<=vertdx[j]) hresv[j]->Fill(dz);
            }
        }
    }

    for(int j=0;j<nmult;j++) { //risoluzione in molteplicità
        hrism->SetBinContent(j+1,hresm[j]->GetRMS(1));
        hrism->SetBinError(j+1,hresm[j]->GetRMSError(1));
    }

    for(int j=0;j<nvert;j++) { //risoluzione nel vertice
        hrisv->SetBinContent(j+1,hresv[j]->GetRMS(1));
        hrisv->SetBinError(j+1,hresv[j]->GetRMSError(1));
    }

    for(int j=0;j<nmult;j++) { //efficienza ed errori binomiali nella molteplicita
        if (hefftotm->GetBinContent(j+1)==0) { //non viene disegnato se non ci sono entries
            heffm->SetBinContent(j+1,0);
            heffm->SetBinError(j+1,0);
        }
        else {
            heffm->SetBinContent(j+1,heffm->GetBinContent(j+1)/hefftotm->GetBinContent(j+1));
            double deps = heffm->GetBinContent(j+1)*(1 - heffm->GetBinContent(j+1)); //epsilon per 1-eps
            deps = deps/hefftotm->GetBinContent(j+1); //diviso il numero di prove
            deps = TMath::Sqrt(deps);
            heffm->SetBinError(j+1,deps);
        }
    }

    for(int j=0;j<nvert;j++) { //efficienza nel vertice
        if (hefftotv->GetBinContent(j+1)==0) { //non viene disegnato se non ci sono entries
            heffv->SetBinContent(j+1,0);
            heffv->SetBinError(j+1,0);
        }
        else {
            heffv->SetBinContent(j+1,heffv->GetBinContent(j+1)/hefftotv->GetBinContent(j+1));
            double deps = heffv->GetBinContent(j+1)*(1 - heffv->GetBinContent(j+1)); //epsilon per 1-eps
            deps = deps/hefftotv->GetBinContent(j+1); //diviso il numero di prove
            deps = TMath::Sqrt(deps);
            heffv->SetBinError(j+1,deps);
        }
    }

    //heffv->Divide(heffv,hefftotv,1,1,"B"); //calcolo l'efficienza sul vertice -> B calcola errori binomiali
    //heffm->Divide(heffm,hefftotm,1,1,"B"); //calcola l'efficienza sulla molteplicità

    if (drawopt.Contains("on") || drawopt.Contains("ON")) {

        gROOT->SetBatch(kTRUE); //no pop up

        ifstream direc("Plot/hrestot.png"); //se esiste già è stato chiamata la macro con draw attivo
        if (!direc) {
            gSystem->Exec("mkdir Plot");
        }

        gStyle->SetOptFit(1);
        TCanvas *c = new TCanvas();
        if (hrestot->Integral()!=0) hrestot->Fit("gaus","Q"); //fit gaussiano a tutte le molteplicità -> NON è gaussiano
        hrestot->SetMarkerStyle(8);
        hrestot->SetMarkerSize(0.5);
        hrestot->Draw("PE");
        c->SaveAs("Plot/hrestot.png");
        c->Close();

        for(int i=0;i<nmult;i++) {
            c = new TCanvas();
            if (hresm[i]->Integral()!=0) hresm[i]->Fit("gaus","Q");
            hresm[i]->SetMarkerStyle(8);
            hresm[i]->SetMarkerColor(kBlue);
            hresm[i]->SetMarkerSize(0.5);
            hresm[i]->Draw("PE");
            sprintf(hname,"Plot/hresm_%d.png",i);
            c->SaveAs(hname);
            c->Close();
        }

        for(int i=0;i<nvert;i++) {
            c = new TCanvas();
            if (hresv[i]->Integral()!=0) hresv[i]->Fit("gaus","Q");
            hresv[i]->SetMarkerStyle(8);
            hresv[i]->SetMarkerColor(kBlue);
            hresv[i]->SetMarkerSize(0.5);
            hresv[i]->Draw("PE");
            sprintf(hname,"Plot/hresv_%d.png",i);
            c->SaveAs(hname);
            c->Close();
        }

        gStyle->SetOptStat(0);
        c = new TCanvas();
        heffv->SetMarkerStyle(8);
        heffv->SetMarkerSize(0.5);
        heffv->Draw("PE");
        c->SaveAs("Plot/heffv.png");
        //heffv->SaveAs("Plot/heffv.root");
        heffv->Write();
        c->Close();

        c = new TCanvas();
        heffm->SetMarkerStyle(8);
        heffm->SetMarkerSize(0.5);
        heffm->Draw("PE");
        c->SaveAs("Plot/heffm.png");
        //heffm->SaveAs("Plot/heffm.root");
        heffm->Write();
        c->Close();

        c = new TCanvas();
        hrism->SetMarkerStyle(8);
        hrism->SetMarkerSize(0.5);
        hrism->Draw("PE");
        c->SaveAs("Plot/hrism.png");
        //hrism->SaveAs("Plot/hrism.root");
        hrism->Write();
        c->Close();

        c = new TCanvas();
        hrisv->SetMarkerStyle(8);
        hrisv->SetMarkerSize(0.5);
        hrisv->Draw("PE");
        c->SaveAs("Plot/hrisv.png");
        //hrisv->SaveAs("Plot/hrisv.root");
        hrisv->Write();
        c->Close();

    } //fine if su drawopt

    histdraw->Write(); //scrive anche gli histo in un file root 
    histoan.Close();

}