#ifndef SIMPOINT_H
#define SIMPOINT_H

#include "TObject.h"

class SimPoint : public TObject{
//----------------------------------------------
//classe che contiene i punti e le funzioni per manipolarli
//----------------------------------------------

private:

    //data members
    double fX; //coordinata x
    double fY;
    double fZ;
    double fEta; //pseudorapidità
    double fTheta_dir; //angolo theta della direzione della particella
    double fPhi_dir; //angolo phi della direzione della particella
    double fPhi_pos; //coordinata phi della posizione della particella 
    double fR; //raggio del layer su cui è la particella
    int fParticella; //label identificativo del particella
    int fLayer; //layer su cui si trova la particella


public:

    //costruttori
    SimPoint();
    SimPoint(double eta, double phiDir, int particella);
    
    //distruttore
    virtual ~SimPoint();

    //getter
    double GetX() const {return fX;}
    double GetY() const {return fY;}
    double GetZ() const {return fZ;}
    double GetEta() const {return fEta;}
    double GetTheta() const {return fTheta_dir;}
    double GetPhiDir() const {return fPhi_dir;}
    double GetPhiPos() const {return fPhi_pos;}
    double GetRaggio() const {return fR;}
    int GetParticella() const {return fParticella;}
    int GetLayer() const {return fLayer;}

    //funzioni custom
    void Crossing(double x0, double y0, double z0, double raggio, bool& acceso);
    void Scattering(double thp, double php);
    void UpdateCyl();
    void PrintStatus();

ClassDef(SimPoint,1)

};

#endif