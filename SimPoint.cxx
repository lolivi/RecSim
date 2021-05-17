#include "Riostream.h"
#include "TMath.h"
#include "SimPoint.h"


ClassImp(SimPoint)

//funzione arcotangente sui quattro quadranti che verrà usata nelle member function successive
double ATan4(double y, double x){
    if (x==0) {
        if (y>=0) return (TMath::Pi())/2.; //semiasse positivo di y
        if (y<0) return TMath::Pi()*3/2; //semiasse negativo di y
    }

    if (x>0 && y>=0) {
        return TMath::ATan(y/x);
    } //primo quadrante, l'arctotangente è già tra 0 e pi/2

    if (x<0 && y>=0) {
        return TMath::ATan(y/x) + TMath::Pi();
    } //secondo quadrante, l'arcotangente è fra -pi/2 e 0, va aggiunto pi

    if (x<0 && y<=0) {
        return TMath::ATan(y/x) + TMath::Pi();
    } //terzo quadrante, l'arcotangente è fra 0 e pi/2, va aggiunto pi 

    if (x>0 && y<=0) {
        return TMath::ATan(y/x) +2*TMath::Pi();
    } //quarto quadrante, l'arcotangente è fra -pi/2 e 0, vanno aggiunti 2pi

    else{
        cout <<"C'è qualcosa di errato!"<<endl;
        return -500;
    }
}

//implementazione del costruttore di default
SimPoint :: SimPoint() : TObject(),
    fX(0.),
    fY(0.),
    fZ(0.),
    fEta(0.),
    fTheta_dir(0.),
    fPhi_dir(0.),
    fPhi_pos(0.),
    fR(0.),
    fParticella(0),
    fLayer(0)
    {}


//implementazione del costruttore standard
SimPoint :: SimPoint(double eta, double phiDir, int particella) : TObject(),
    fX(0.),
    fY(0.),
    fZ(0.),
    fEta(eta),
    fTheta_dir( 2*TMath::ATan(TMath::Exp(-eta)) ),
    fPhi_dir(phiDir),
    fPhi_pos(0.),
    fR(0.),
    fParticella(particella),
    fLayer(-1)
    {}


//implementazione del distruttore
SimPoint :: ~SimPoint(){}


//implementazione della member function Crossing
//Crossing assegna/cambia i valori a fX, fY, fZ usando le intersezioni delle rette della direzione con i vari strati
//viene anche incrementato il valore di fLayer
//la funzione prende come argomenti le coordinate del punto di intersezione sul layer precendente (alla prima chiamata del vertice),
//del raggio del layer su cui si vuole effettuare l'intersezione, ed una booleana per autorizzare o meno lo scattering
void SimPoint :: Crossing(double x0, double y0, double z0, double raggio, bool& acceso){
    
    //i valori di fX, fY, fZ vengono cambiati solo se acceso=true
    if(acceso){

        //definizione di quantità che servono per calcolare l'intersezione
        double c1 = TMath::Sin(fTheta_dir) * TMath::Cos(fPhi_dir);
        double c2 = TMath::Sin(fTheta_dir) * TMath::Sin(fPhi_dir);
        double c3 = TMath::Cos(fTheta_dir);

        double delta = (x0*c1 + y0*c2)*(x0*c1 + y0*c2) - (c1*c1 + c2*c2)*(x0*x0 + y0*y0 - raggio*raggio);

        //si controlla che che il discriminante delta sia positivo
        if(delta<0){
            cout<<"Attenzione! Discriminante negativo!"<<endl;
        }

        //vengono calcolate le 2 soluzioni e viene scelta quella positiva
        double t = (-(x0*c1 + y0*c2) + TMath::Sqrt(delta) )/(c1*c1 + c2*c2);
        if(t<0){
            t = (-(x0*c1 + y0*c2)-TMath::Sqrt(delta))/(c1*c1 + c2*c2);
        }
        if(t<0){
            cout<<"Attenzione! 2 soluzioni negative!"<<endl;
            cout<<t<<endl;
        }

        //assegnazione a fX, fY, fZ dei valori degli hit
        fX = x0 + c1*t;
        fY = y0 + c2*t;
        fZ = z0 + c3*t;

        //viene segnato che la particella ha cambiato layer
        fLayer++;

    }

    //la booleana viene messa a false dopo il crossing per evitare problemi dovuti a chiamate consecutive non volute
    acceso = false;
}


//implementazione della member function Scattering
//Scattering cambia i valori di fTheta_dir e fPhi_dir in seguito allo scattering multiplo
//la funzione prende come argomenti i nuovo valori di direzione nel sistema di riferimento della particella
void SimPoint :: Scattering(double thp, double php){

    double mr[3][3];
    mr[0][0] = -TMath::Sin(fPhi_dir);
    mr[1][0] = TMath::Cos(fPhi_dir);
    mr[2][0] = 0.;
    mr[0][1] = -TMath::Cos(fTheta_dir) * TMath::Cos(fPhi_dir);
    mr[1][1] = -TMath::Cos(fTheta_dir) * TMath::Sin(fPhi_dir);
    mr[2][1] = TMath::Sin(fTheta_dir);
    mr[0][2] = TMath::Sin(fTheta_dir) * TMath::Cos(fPhi_dir);
    mr[1][2] = TMath::Sin(fTheta_dir) * TMath::Sin(fPhi_dir);
    mr[2][2] = TMath::Cos(fTheta_dir);

    double mp[3];
    mp[0] = TMath::Sin(thp) * TMath::Cos(php);
    mp[1] = TMath::Sin(thp) * TMath::Sin(php);
    mp[2] = TMath::Cos(thp);

    double cd[3];
    for(int i=0; i<3; i++){
        cd[i] = 0.;
        for(int j=0; j<3; j++){
            cd[i] += mr[i][j]*mp[j];
        }
    }

    double cnorm = cd[0]*cd[0]+ cd[1]*cd[1]+cd[2]*cd[2];
    cnorm = TMath::Sqrt(cnorm);

    fTheta_dir = TMath::ACos(cd[2]/cnorm);
    fPhi_dir = ATan4(cd[1],cd[0]);  

}

//implementazione di PrinStatus
//questa funzione stampa lo stato della particella al momento della chiamata della funzione
//IMPORTANTE! non tutti i data members vengono aggiornati in automatico durante il camminio della particella
void SimPoint :: PrintStatus(){
    cout<<" Particella n. = "<<fParticella<<endl;
    cout<<" {x,y,z} = {"<<fX<<", "<<fY<<", "<<fZ<<"}"<<endl;
    cout<<" Theta direzione = "<<fTheta_dir<<endl;
    cout<<" Phi direzione = "<<fPhi_dir<<endl;
    cout<<" {R, phi, z} = {"<<fR<<", "<<fPhi_pos<<", "<<fZ<<"}"<<endl;
    cout<<" La particella è sul layer n. = "<<fLayer<<endl;
    
}

//iplementazione di UpdateCyl
//UpdateCyl aggiorna le coordinate cilindriche, usate per salvare le hit
void SimPoint :: UpdateCyl(){
    fR = TMath::Sqrt(fX*fX + fY*fY);
    fPhi_pos = ATan4(fY,fX);
}

