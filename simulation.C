#include "SimPoint.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TStopwatch.h"
#include "Riostream.h"
#include "TMath.h"
#include "TRandom3.h"

//costanti globali
#ifndef COSTANTI
#define COSTANTI
    const double kR0 = 30.; //raggio beam pipe
    const double kR1 = 40.; //raggio layer 1
    const double kR2 = 70.; //raggio layer 2
    const double kHalfLength = 135; //metà lunghezza layer 1 e layer 2
#endif

//struct per salvare le hit in coordinate cilindriche
#ifndef HCOORD_H
#define HCOORD_H
    struct HCoord : public TObject{
        double z;
        double cPhi;
        int label;
        ClassDef(HCoord,1);
    };
#endif

//----------------------
//DICHIARAZIONE FUNZIONI
//----------------------

int multsimB(); //estrazione della molteplicità estraendola da una distribuzione uniforme
int multsimC(); //estrazione della molteplicità estraendola da una molteplicità assegnata
double etasimC(); //estrazione della pseudorapidità in una distribuzione uniforme
double vertsimB(); //estrazione del vertice entro un numero sigma
double vertsimC(); //estrazione del vertice uniforme
double thetap (); //angolo theta simulato per lo scattering multiplo
double phip (); //angolo phi simulato per lo scattering multiplo
TH1F* newheta (TH1F*); //funzione che restituisce la distribuzione di eta fra -2 e 2


void simulation (TString scattopt = "A", TString multopt = "A", TString etaopt = "A",TString vertopt = "A" ,int kEvents = 1000000, unsigned int seed = 5678) {

    //----------------------------------
    //SCELTE INIZIALI PER LA SIMULAZIONE
    //----------------------------------

    TStopwatch timer; 
    timer.Start();

    gRandom->SetSeed(seed); 

    bool scatt; //scattering acceso o spento
    if (scattopt.Contains("A")) scatt=true; 
    if (scattopt.Contains("B")) scatt=false;

    bool checkcross = true; //utilizzata nell'intersezione come controllo
    //serve a non chiamare due volte di fila l'intersezione senza salvare
    bool firstl = false; //diventa true se c'è una hit sul layer1
    bool secondl = false; //secondo layer 
    int both = 0; //conta hit su entrambi i layer
    int cnt1 = 0; //hit solo sul primo layer
    int cnt2 = 0; //hit solo sul secondo layer
    int tot = 0; //numero di volte totali

    //istogrammi delle distribuzioni assegnate
    TFile* distr = new TFile(); //sono inizializzati vuoti -> vengono riempiti solo se necessario
    TH1F* multhisto = new TH1F();
    TH1F* etahisto = new TH1F();
    TH1F* neweta = new TH1F();

    if (etaopt.Contains("A") || etaopt.Contains("B") || multopt.Contains("A")) {
        ifstream input("kinem.root");
        if(!input) {
            cout<<"\n ATTENZIONE: Il file kinem.root NON esiste \n";
            return; //esce dalla simulazione
        }
        distr = new TFile("kinem.root"); //file delle distribuzione assegnate 
        if (multopt.Contains("A")) multhisto = (TH1F*)distr->Get("hmul"); //distribuzione molteplicità
        etahisto = (TH1F*)distr->Get("heta"); //distribuzione pseudorapidità
        if (etaopt.Contains("A")) neweta = newheta(etahisto); //distribuzione pseudorapidità tra -2 e 2
    }

    //------------------
    //INIZIO SIMULAZIONE
    //------------------

    TFile simfile("simfile.root","recreate"); //file in scrittura
    TTree *tree = new TTree ("tree","Vertex and Hits");

    //def struct con coordinata z del vertice e molteplicità
    typedef struct {
        double Z;
        int mult;} 
    VMult;

    static VMult vertex; //struct con coordinata z del vertice e molteplicità
    static HCoord hitPoint; //struct con coordinate cilindriche delle hit e label della particella

    //il tree viene riempito con una struct che contiene vertice e molteplicità e 2 TClonesArray che contengono le hits sui layer 1 e 2
    TClonesArray *ptrhits1 = new TClonesArray("HCoord",100); //hits sul layer 1
    TClonesArray *ptrhits2 = new TClonesArray("HCoord",100); //hits sul layer 2
    TClonesArray &hits1 = *ptrhits1;
    TClonesArray &hits2 = *ptrhits2;
    
    tree->Branch("VertMult",&vertex.Z,"Z/D:mult/I");
    tree->Branch("HitsL1",&ptrhits1);
    tree->Branch("HitsL2",&ptrhits2);

    SimPoint* hitpt = new SimPoint(); //hitpt è l'oggetto che rappresenta le particelle
    //viene qui allocata la memoria per tutte le particelle, che saranno create per copia da questo indirizzo
    int multfix = 0; //molteplicità fissa, se richiesta dall'utente

    cout<<" -------------------------------------------------------------- \n SIMULAZIONE EFFETTUATA PER "<<kEvents<<" EVENTI CON UN SEED PARI A "<<seed<<endl<<endl;
    cout<<"\n DATI: \n";
    cout<<" - RAGGIO BEAM PIPE: "<<kR0<<" mm \n";
    cout<<" - RAGGIO LAYER 1: "<<kR1<<" mm \n";
    cout<<" - RAGGIO LAYER 2: "<<kR2<<" mm \n";
    cout<<" - ESTENSIONE LAYER 1 e 2: 270 mm \n \n -------------------------------------------------------------- \n ";

    for(int i=0;i<kEvents;i++){ //loop su eventi

        //estrazione della molteplicità
        if (multopt.Contains("A")) vertex.mult = (int)(multhisto->GetRandom()+0.5); //serve per non prendere 2 dall'istogramma (che parte da 2.5)
        else if (multopt.Contains("B")) vertex.mult=multsimB();
        else {
            if (i==0) multfix = multsimC();
            vertex.mult=multfix;
        }

        //estrazione della posizione vertice
        double Xv=gRandom->Gaus(0,0.1); //coordinata x del vertice, non è nella struct perché non verrà salvata nel tree
        double Yv=gRandom->Gaus(0,0.1); //coordinata y del vertice
        if (vertopt.Contains("A")) vertex.Z=gRandom->Gaus(0,53); //estrazione del vertice da una gaussiana con mu=0 e sigma=53mm
        else if (vertopt.Contains("B")) vertex.Z=vertsimB(); //estrazione entro un numero di sigma
        else if (vertopt.Contains("C")) vertex.Z=vertsimC(); //estrazione uniforme

        int pos1=0; //indici di hits1 e hits2 in L1 e L2
        int pos2=0; //In questo modo il TClonesArray è riempito senza buchi di memoria negli hit non accettabili geometricamente
        
        //j indica il numero della particella
        for(int j=0;j<vertex.mult;j++) { //loop sulla molteplicità
            firstl=false;
            secondl=false;
            tot++;

            //si estrae la direzione phi da 0 a 2Pi per ogni particella
            double phi = 2*TMath::Pi()*gRandom->Rndm();
            //estrazione della pseudorapidità
            double eta;
            if (etaopt.Contains("B")) eta = etahisto->GetRandom();
            else if (etaopt.Contains("A")) eta = neweta->GetRandom();
            else eta = etasimC();

            new(hitpt) SimPoint(eta,phi,j); //costruiamo sull'indirizzo di hitpt senza allocare memoria
            //intersezione con la beam pipe
            hitpt->Crossing(Xv,Yv,vertex.Z,kR0,checkcross); //prende la posizione iniziale del vertice e il raggio della beam pipe 
            //si effettua lo scattering solo se scatt=true
            if(scatt) hitpt->Scattering(thetap(),phip()); 
            //intersezione con il layer 1
            checkcross=true; //occorre rimettere questa booleana a true perché viene messa a false dopo l'intersezione della particella con la beam pipe
            hitpt->Crossing(hitpt->GetX(),hitpt->GetY(),hitpt->GetZ(),kR1,checkcross); //intersezione della particella con L1
            if (TMath::Abs(hitpt->GetZ()) > kHalfLength) scatt = false; //hit non accettabili su L1, ma potrebbero essere accettate su L2
            if(scatt) hitpt->Scattering(thetap(),phip()); //scattering su L1
            if (TMath::Abs(hitpt->GetZ()) <= kHalfLength) {
                firstl = true;
                hitpt->UpdateCyl(); //aggiornamento delle coordiante cilindriche
                //la coord radiale non viene salvata, ma si controlla che "coincida" con il raggio del rivelatore
                if( TMath::Abs(hitpt->GetRaggio()-kR1)>pow(10,-12)) cout<<"problema!"<<endl; 
                //aggiornamento dei data member della struct con le coordinate delle hit
                hitPoint.z=hitpt->GetZ();
                hitPoint.cPhi=hitpt->GetPhiPos();
                hitPoint.label=hitpt->GetParticella();
                new(hits1[pos1++]) HCoord(hitPoint); //salvataggio della struct nei TClonesArray
            }//vengono salvate nel TClonesArray solo le hit che possono essere geometricamente accettate

            //intersezione con il layer 2
            checkcross=true; 
            hitpt->Crossing(hitpt->GetX(),hitpt->GetY(),hitpt->GetZ(),kR2,checkcross);
            if (TMath::Abs(hitpt->GetZ()) <= kHalfLength) {
                secondl=true;
                hitpt->UpdateCyl();
                if( TMath::Abs(hitpt->GetRaggio()-kR2)>pow(10,-12)) cout<<"problema!"<<endl;
                hitPoint.z=hitpt->GetZ();
                hitPoint.cPhi=hitpt->GetPhiPos();
                hitPoint.label=hitpt->GetParticella();
                new(hits2[pos2++]) HCoord(hitPoint);
            }
            checkcross=true;
            if (firstl && secondl) both++;
            else if (firstl && !secondl) cnt1++;
            else if (!firstl && secondl) cnt2++;
        } //fine loop su molteplicità

    if (kEvents<10) { //se si hanno meno di 10 eventi vengono stampati a video i risultati di tutti
        printf(" \n Simulazione Evento %d con molteplicità: %d\n",i,vertex.mult);
        printf(" VERTICE: X=%f ; Y=%f ; Z=%f \n",Xv,Yv,vertex.Z);
        printf(" Hits nel Layer 1: %d\n",ptrhits1->GetEntries());
        printf(" Hits nel Layer 2: %d\n",ptrhits2->GetEntries());
    }

    else if (i%(kEvents/10)==0) { //se si hanno più di 10 eventi vengono stampati a video solo 10 eventi
        printf(" \n Simulazione Evento %d con molteplicità: %d\n",i,vertex.mult);
        printf(" VERTICE: X=%f ; Y=%f ; Z=%f \n",Xv,Yv,vertex.Z);
        printf(" Hits nel Layer 1: %d\n",ptrhits1->GetEntries());
        printf(" Hits nel Layer 2: %d\n",ptrhits2->GetEntries());
    }

    //viene riempito il tree e ripulita la memoria su cui si salvano le hit per riutilizzare la stessa porzione di memoria per la particella successiva
    tree->Fill();
    ptrhits1->Clear();
    ptrhits2->Clear();

    }//fine loop su eventi

    distr->Close(); //si chiude il file sulle distribuzioni assegnate
    //il file contente il tree con le informazioni sulla simulazione viene riempito e successivamente chiuso
    simfile.Write();
    simfile.Close();
    
    //vengono stampati a video i dati relativi all'accettanza relativi alla simulazione appena effettuata
    cout<<"\n ---------------------------------- \n ACCETTANZA \n";
    cout<<" Percentuale di Hit su L1 e L2 : "<<(double)both/tot*100<<"%\n";
    cout<<" Percentuale di Hit SOLO su L1: "<<(double)cnt1/tot*100<<"%\n";
    cout<<" Percentuale di Hit SOLO su L2: "<<(double)cnt2/tot*100<<"%\n";

    cout<<" \n ------------------------------------------ \n PERFORMANCE \n ";
    timer.Stop();
    timer.Print();

}

//------------------------
//IMPLEMENTAZIONE FUNZIONI
//------------------------

int multsimB(){
    static int multcall = 0;
    static int a,b;
    if (multcall==0) {
        cout<<" \n RANGE MOLTEPLICITA' UNIFORME \n";
        do{
            cout<<" Inserire numeri positivi e crescenti\n";
            cout<<" Min: ";
            cin>>a;
            cout<<" Max: ";
            cin>>b;
        } while (a<=0 || b<=0 || a>b);
    }
    multcall++; //in questo modo si scelgono gli estremi una volta sola 
    return (int)(a+(b-a)*gRandom->Rndm());
}  

int multsimC() {
    cout<<" \n MOLTEPLICITA' FISSA = ";
    int n;
    do {
        cin>>n;
    } while (n<=0);
    return n;
}

double etasimC() {
    static int etacall = 0;
    static double a,b;
    if (etacall==0) {
        cout<<"\n RANGE PSEUDORAPIDITA UNIFORME \n";
        do{
            cout<<" Inserire numeri fra -6 e 6 in modo crescente \n";
            cout<<" Min: ";
            cin>>a;
            cout<<" Max: ";
            cin>>b;
        } while (a<-6 || b>6 || a>b);
    }
    etacall++; //in questo modo si scelgono gli estremi una volta sola 
    return a+(b-a)*(gRandom->Rndm());
}

double vertsimB() { //simulazione del vertice entro un numero sigma
    static int vertcall=0;
    static int nsig; //numero di sigma entro il quale si vuole il vertice
    double vert; //vertice restituito
    if(vertcall==0) {
        cout<<"\n NUMERO DI SIGMA ENTRO IL QUALE SI VUOLE IL VERTICE \n";
        do {
            cout<<" Inserire un intero positivo: ";
            cin>>nsig;
        } while (nsig<0);
    } 
    do {
        vert=gRandom->Gaus(0,53);
    } while (TMath::Abs(vert)>(double)nsig*53); //finche non è entro nsig volte la sigma, cioè 53 mm
    vertcall++;
    return vert;
}

double vertsimC() {
    static int vertcall = 0;
    static double vertmin,vertmax;
    if (vertcall==0) {
        cout<<"\n RANGE VERTICE UNIFORME \n";
        do{
            cout<<" Inserire numeri [mm] in modo crescente \n";
            cout<<" Min: ";
            cin>>vertmin;
            cout<<" Max: ";
            cin>>vertmax;

        } while (vertmin>vertmax);
    }
    vertcall++; //in questo modo si scelgono gli estremi una volta sola 
    return vertmin+(vertmax-vertmin)*(gRandom->Rndm());
}

double phip () {
    return 2*TMath::Pi()*gRandom->Rndm();
}

double thetap () {
    return gRandom->Gaus(0,0.001);
}

TH1F* newheta (TH1F* oldheta) {  

    TAxis *xaxis=oldheta->GetXaxis();
    double step = xaxis->GetBinWidth(1);
    int b1=xaxis->FindBin(-2.);
    int b2=xaxis->FindBin(2.);
    double xlow=xaxis->GetBinLowEdge(b1);
    double xhig=xaxis->GetBinUpEdge(b2);
    int nbins=b2-b1+1;
    double step2 = (xhig-xlow)/nbins;

    TH1F* newheta = new TH1F("newheta","#eta Distribution [-2,2] ",nbins,xlow,xhig);
    int j=1;
    for(int i=b1;i<=b2;i++) newheta->SetBinContent(j++,oldheta->GetBinContent(i));  

    return newheta;

}







