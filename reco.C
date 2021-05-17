#include "TRandom3.h"
#include "TH1F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "Riostream.h"


//costanti globali
const double kDPhiMax = 0.00412939; //dispersione dei deltaPhi, calcolata come 3 Dev Std della distrib dei deltaPhi

#ifndef COSTANTI
#define COSTANTI
    const double kR0 = 30.; //raggio beam pipe
    const double kR1 = 40.; //raggio layer 1
    const double kR2 = 70.; //raggio layer 2
    const double kHalfLength = 135; //metà lunghezza layer 1 e layer 2
#endif

//struct con coordinate cilindriche degli hit
#ifndef HCOORD_H
#define HCOORD_H
    struct HCoord : public TObject{
        //double raggio;
        double z;
        double cPhi;
        int label;
        ClassDef(HCoord,1);
    };
#endif

//------------------------
//DICHIARAZIONE FUNZIONI
//------------------------

double tracklet (double z1, double z2, double r1, double r2); //retta passante per due punti, restituisce l'intersezione con l'asse z
double findzvec(std::vector<double> &v); //fornisce il valore del vertice, se ci sono state le condizioni per riconoscerlo, altrimenti restituisce -500
double findzhist (TH1F* hist , std::vector<double> &v); //implementazione dell'algoritmo con istogramma e vector
void printvec(std::vector<double> v); //permette di stampare a video il vector dei possibili vertici
void createNoise(double &z, double &phi); //modifica i valori di z e phi che gli vengono passati creando il niose
void smearing(double &phi, double &z, double r); //modifica il valori di z e phi che gli vengono passati per simulare lo smearing del rivelatore
void hidePartLabel(int &label); //nasconde le informazioni sul label della particella, così da essere sicuri di trattare hit simulate come fossero hit reali


//------------------------
//MAIN DELLA RICOSTRUZIONE
//------------------------

void reco (unsigned int seed = 675,int noise = 10) {
    
    //----------------------------------
    //SCELTE INIZIALI PER LA RICOSTRUZIONE
    //----------------------------------

    TStopwatch timer;
    timer.Start();

    int win = 0; //calcola il numero di volte che è stato possibile ricavare il vertice

    gRandom->SetSeed(seed);

    //vertice e molteplicità
    typedef struct {
        double Z;
        int mult;} 
    VMult;

    static VMult vertex; //struct con coordinata z del vertice e molteplicità
    static HCoord point; //struct con coordinate cilindriche delle hit e label della particella

    //TClonesArray per la lettura
    TClonesArray *ptrhits1 = new TClonesArray("HCoord",100); //hits sul layer 1
    TClonesArray *ptrhits2 = new TClonesArray("HCoord",100); //hits su layer 2
    TClonesArray &hits1 = *ptrhits1;
    TClonesArray &hits2 = *ptrhits2;

    //apertura in lettura del file di simulazione
    ifstream input("simfile.root");
    if(!input) {
        cout<<"\nATTENZIONE: Non esiste il file simfile.root \n";
        return;
    }
    TFile simfile("simfile.root");

    //tree in lettura e dichiarazione branch
    TTree *oldtree = (TTree*)simfile.Get("tree");
    TBranch *b1 = oldtree->GetBranch("VertMult");
    TBranch *b3 = oldtree->GetBranch("HitsL1");
    TBranch *b4 = oldtree->GetBranch("HitsL2");

    //indirizzi di lettura
    b1->SetAddress(&vertex.Z); 
    b3->SetAddress(&ptrhits1);
    b4->SetAddress(&ptrhits2);

    //numero di eventi
    int kEvents = oldtree->GetEntries();

    //istogramma e vector per la ricostruzione
    double binsz = 5; //5mm di bin size nell'istogramma
    int bins = (int)400/binsz; //numero di bin
    TH1F* zhisto = new TH1F(" zhisto ", " Calcolo Vertice ",bins,-200,200); //istogramma per il calcolo del vertice

    std::vector<double> zvec; //vector delle possibili posizioni del vertice

    //file in scrittura per l'analisi 
    TFile* histoan = new TFile("histoAnalysis.root","recreate"); 
    TTree *newtree = new TTree ("newtree","Reconstruction");
    typedef struct { //salviamo per ogni evento una struct
        double zrec; //z ricostruita
        double ztrue; //z vera
        int mult;} //molteplicità
    VRec; //struttura con parametri della ricostruzione
    static VRec recvert; 
    newtree->Branch("RecVert",&recvert.zrec,"zrec/D:ztrue:mult/I");

    cout<<" ------------------------- \n PARAMETRI RICOSTRUZIONE: \n";
    cout<<" Seed: "<<seed<<endl;
    cout<<" Numero di eventi: "<<kEvents<<endl;
    cout<<" Risoluzione Phi: "<<kDPhiMax<<" mrad\n";
    cout<<" Numero di Punti Spuri: "<<noise<<endl;
    cout<<" ------------------------- \n";

    //si popolano i TClonesArray e vertex
    for(int n=0;n<kEvents;n++) { //loop su eventi
        oldtree->GetEvent(n);
        zvec.reserve(vertex.mult); //prepariamo un numero di posti in memoria nel vector
        int entries1 = ptrhits1->GetEntries(); //numero vero di hit nel layer 1
        int entries2 = ptrhits2->GetEntries(); //layer 2
        //Si aggiungono i punti spuri al layer 1 e al layer 2 nei TClonesArray di questo evento
        for(int i=0;i<noise;i++) {
            createNoise(point.z, point.cPhi);
            new(hits1[entries1+i]) HCoord(point);
            createNoise(point.z, point.cPhi);
            new(hits2[entries2+i]) HCoord(point);        
        }
        
        for(int i=0;i<(entries1+noise);i++) { //loop su layer 1
            HCoord *hitrc1 = (HCoord*)ptrhits1->At(i); //prendo l'i-esimo elemento del TClonesArray
            if (i<entries1) { 
                hidePartLabel(hitrc1->label); //nascondo il label della particella
                smearing(hitrc1->cPhi,hitrc1->z,kR1); //effettuo lo smearing sulle hit del layer 1
            } 
            for (int j=0;j<(entries2+noise);j++) { //loop su particelle layer 2
                HCoord* hitrc2 = (HCoord*)ptrhits2->At(j); //prendo il j-esimo elemento del TClonesArray 
                if (j<entries2 && i==0) { //lo smearing si effettua solo una volta per le particelle del layer 2 e non sui punti spuri che sono già random
                    hidePartLabel(hitrc2->label);
                    smearing(hitrc2->cPhi,hitrc2->z,kR2);
                }
                double dphi = TMath::Abs(hitrc1->cPhi - hitrc2->cPhi); //modulo della differenza negli angoli phi
                if (dphi < kDPhiMax) {
                    double z = tracklet(hitrc2->z,hitrc1->z,kR2,kR1); //restituisce intersezione con asse Z
                    zvec.push_back(z);
                    zhisto->Fill(z);
                }
            } //fine loop su layer 2
        } //fine loop sul layer 1

        double zvert = findzvec(zvec); //implementazione con vector
        //double zvert = findzhist(zhisto,zvec); //implementazione con histo e vector
        if (zvert!=-500) win++;

        recvert.zrec = zvert; //-500 per i non riusciti
        recvert.ztrue = vertex.Z;
        recvert.mult = vertex.mult;
        newtree->Fill(); //fill del tree in scrittura

        if (kEvents<10) {
            cout<<" \nRicostruzione Evento numero: "<<n<<endl;
            cout<<" Molteplicità: "<<vertex.mult<<endl;
            if (zvert!=-500) cout<<" Vertice Ricostruito: "<<zvert<<"\n Vertice Vero: "<<vertex.Z<<"\n";
            else cout<<" Non è stato possibile costruire il vertice che si trovava in "<<vertex.Z<<"\n";
        }

        else if (n%(int)(kEvents/10)==0) {
            cout<<" \n Ricostruzione Evento numero: "<<n<<endl;
            cout<<" Molteplicità: "<<vertex.mult<<endl;
            if (zvert!=-500) cout<<" Vertice Ricostruito: "<<zvert<<"\n Vertice Vero: "<<vertex.Z<<"\n";
            else cout<<" Non è stato possibile costruire il vertice che si trovava in "<<vertex.Z<<"\n";
        }
        
        zvec.clear(); //reset zvec
        zhisto->Reset(); //anche l'istogramma
    } //finisce loop su eventi

    histoan->Write(); //scriviamo il file di analisi
    simfile.Close(); //chiudiamo il file di lettura
    cout<<" \n EFFICIENZA = "<<(double)win/kEvents<<endl; 
    cout<<" \n PERFORMANCE: \n";
    timer.Stop();
    timer.Print();
} //Fine del Main


//------------------------
//IMPLEMENTAZIONE FUNZIONI
//------------------------

double tracklet (double z1, double z2, double r1, double r2) {
    double k = (r1*z2-r2*z1)/(r1-r2);
    return k; //intersezione con asse
}

double findzvec (std::vector<double> &v) {

    if(v.empty()) return -500;
    std::sort(v.begin(),v.end()); //riordina da crescente a decrescente
    int sz = (int)v.size(); //nuova size del vector

    double zmax = 0; //primo massimo in media
    int cntmax = 0; //queste sono le entries
    int cntmax2 = 0; //queste sono le entries del secondo massimo
    double mean = 0.; //iteratori che si aggiornano a ogni posizione della running window
    int cnt = 0; //media ed entries
    double left = v[0]; //la running window parte dal minimo
    double right = left+5; //ha larghezza di 5 mm
    unsigned int iMin = 0; //parto dal primo elemento del vector

    while (left<=v[sz-1]) { //spostamento della running window 
        cnt = 0; //ogni volta resetto media e conteggio
        mean = 0; 
        bool cambia = true; //ogni volta che esce dal for sugli elementi del vector, cambia torna true
        for(unsigned int i=iMin; i<v.size(); i++) { //partiamo da iMin
            if (v[i]>=left && v[i]<=right) {
                cnt++; //calcolo le entries e la media degli elementi fra left e right (estremi della running window)
                mean = mean+v[i];
                if(cambia){ //parto a guardare il vettore dal primo elemento preso nel ciclo prima invece di iniziare tutte le volte da zero
                    iMin = i; //il primo elemento del vector che soddisfa i requisiti della running window
                    cambia = false; //cambia diventa false e finché v[i] è nel range della running window, iMin non viene più aggiornato
                }
            }    
            else if (v[i]>right) { //visto che il vettore v è ordinato se supera l'estremo destro della running window non può più avere ingressi
                break; //esce dal for
            }
        }

        if (cnt>cntmax) { //se nella running window il conteggio supera il massimo, allora si aggiorna il massimo e la media corrispondente
            cntmax = cnt;
            zmax = mean/cnt;
        }
        else if (cnt>cntmax2 && TMath::Abs(mean/cnt-zmax)>5.) cntmax2 = cnt;

        left = left + 1; //sposto la running window a destra
        right = right + 1; //sposto a destra 
    } //fine running window
    
    if(cntmax==0 && cntmax2==0) return -500; //no entries
    else if (cntmax>2*cntmax2) return zmax; //il picco massimo domina sull'altro se è più del doppio
    else return -500;
}

void printvec(std::vector<double> v) {
    cout<<" Z = {";
    for(int i=0;i<(int)v.size()-1;i++){
        cout<<v[i]<<", ";
    }
    cout<<v[v.size()-1]<<"}\n";
}

double findzhist (TH1F* hist, std::vector<double> &zvec) { 

    if(zvec.empty()) return -500; //non ci sono accoppiamenti
    int imax = hist->GetMaximumBin(); //restituisce l'indice del bin con il valore massimo
    double binmax = hist->GetBinContent(imax); //salvo il valore corrispondente al bin del massimo
    if (binmax==0) return -500; //non ci sono entries nell'istogramma
    hist->SetBinContent(imax,0.); //lo metto a 0
    int imax2 = hist->GetMaximumBin(); //secondo massimo
    double binmax2 = hist->GetBinContent(imax2); //valore corrispondente al secondo massimo
    double zctr; //zcentrale
    double width; //larghezza per media nel vector
    if (binmax>2*binmax2) { //il massimo domina sul secondo picco
        zctr = hist->GetXaxis()->GetBinCenter(imax); //caso in cui c'è un solo picco nella cella imax o picco dominante
        width = 2.5;
    }
    else if (TMath::Abs(imax-imax2)==1 && binmax2!=0) { //i due picchi si riferiscono allo stesso z
        if (imax>imax2) zctr=hist->GetXaxis()->GetBinLowEdge(imax2)+hist->GetXaxis()->GetBinUpEdge(imax);
        else zctr=hist->GetXaxis()->GetBinLowEdge(imax)+hist->GetXaxis()->GetBinUpEdge(imax2);
        zctr = zctr/2.; //restituisco il valore centrale come media
        width = 5.;
    }
    else return -500; //ci sono due picchi ben distinti
    double mn = 0; //numeratore della media
    int zcnt = 0; //contatore delle Z nel vector
    double zvert; //z finale del vertice
    for(int i=0;i<(int)zvec.size();i++) {
        if (TMath::Abs(zvec[i]-zctr)<=width) {
            zcnt++;
            mn=mn+zvec[i];
        }
    }
    if (zcnt==0) {
        cout<<" Errore -> No Entries nel Vector \n";
        return -500;
    }
    zvert = (double)mn/zcnt;
    return zvert;
}

void createNoise(double &z, double &phi){
    //si simulano i punti spuri aggiornando le coodinate passate con dei valori estratti uniformente sulla superficie dei layer
    z = -kHalfLength+2*kHalfLength*(gRandom->Rndm());
    phi = 2*TMath::Pi()*(gRandom->Rndm());
    
}

void smearing(double &phi, double &z, double r){
    phi = phi + (gRandom->Gaus(0.,0.03))/r;
    z = z + gRandom->Gaus(0.,0.120);
}

void hidePartLabel(int &label){
    label = 0;
}