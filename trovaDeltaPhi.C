#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TStopwatch.h"

//dichiarazione funzione per lo smeraring
void Smearing(double &phi, double &z, double r);

#ifndef HCOORD_H
#define HCOORD_H
//struttura con coordinate cilindriche degli hit
    struct HCoord : public TObject{
        double z;
        double cPhi;
        int label;
        ClassDef(HCoord,1);
    };
#endif

//funzione da eseguire
void trovaDeltaPhi(unsigned int seed = 567){

    TStopwatch timer;
    timer.Start();

    gRandom->SetSeed(seed);

    //definizione struct e TClonesArray analoghi a quelli salvati nel TTree che si vuole leggere
    typedef struct {
        double Z;
        int mult;} 
    VMult;
    static VMult vertex; 

    TClonesArray *hits1 = new TClonesArray("HCoord",100);
    TClonesArray *hits2 = new TClonesArray("HCoord",100);

    //apertura file di input
    ifstream input("simfile.root");
    if(!input) {
        cout<<"\nATTENZIONE: Non esiste il file simfile.root \n";
        return;
    }
    TFile hfile("simfile.root");

    //lettura TTree  e TBranch
    TTree *tree = (TTree*)hfile.Get("tree");
    TBranch *bv=tree->GetBranch("VertMult");
    TBranch *b1=tree->GetBranch("HitsL1");
    TBranch *b2=tree->GetBranch("HitsL2");

    //definizione degli indirizzi per la lettura dei dati su ttree
    bv->SetAddress(&vertex.Z);
    b1->SetAddress(&hits1);
    b2->SetAddress(&hits2);

    //dichiarazione dell'istogramma della distribuzione dei deltaPhi
    TH1D *h = new TH1D("h","Istogramma dei #Delta#varphi",100,-0.006,0.006);
    
    //loop sugli ingressi nel TTree
    for(int ev=0; ev<tree->GetEntries(); ev++){
        tree->GetEvent(ev);
        for (int j=0; j<hits1->GetEntries(); j++){ //loop su layer 1
            HCoord *t1 = (HCoord*)hits1->At(j);
            if(t1->cPhi==0) cout<<"Attenzione -> Phi Cilindrica è Zero "<<endl;
            Smearing(t1->cPhi,t1->z,40);
            for(int k=0; k<hits2->GetEntries(); k++){ //loop su layer 2
                HCoord *t2 = (HCoord*)hits2->At(k);
                if(t1->label == t2->label){ //si selezionano hit con lo stesso label, cioè dovuti alla stessa particella
                    if(t2->cPhi==0) cout<<"Attenzione -> Phi Cilindrica è Zero "<<endl;
                    Smearing(t2->cPhi,t2->z,70);
                    h->Fill(t2->cPhi - t1->cPhi); //riempimento dell'istogramma con il deltaphi 
                    break;
                } //fine dell'if sulle particelle del L1 e L2 con stesso label
            } //fine loop su layer 2
        }//fine loop su layer 1
    }//fine loop su eventi 

    h->SaveAs("istoPhi.root","recreate");
    double var = h->GetRMS();
    var = var*3; //la moltiplica per 3

    cout<<endl<<"Dispersione dei deltaPhi (3 sigma) = "<<var<<" rad"<<endl;
    cout<<" \n PERFORMANCE: \n";

    timer.Stop();
    timer.Print();
}

void Smearing(double &phi, double &z, double r){
    phi = phi + (gRandom->Gaus(0.,0.03))/r;
    z = z + gRandom->Gaus(0.,0.120);
}