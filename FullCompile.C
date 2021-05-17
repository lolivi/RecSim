
void FullCompile (TString myopt = "fast") {

    TString opt; //compila o con il + (fast) o con il ++ (force)
    if (myopt.Contains("force")) {
        opt = "kfg";
    }
    else {
        opt = "kg";
    }


    //-----------------------
    //SIMULAZIONE
    //-----------------------
    
    //Opzioni Standard
    TString scattopt = "A"; //opzione std per lo scattering (acceso)
    TString multopt = "A"; //opzione per la molteplicità 
    TString etaopt = "A"; //opzione per la pseudorapidità
    TString vertopt = "A"; //opzione per il vertice
    TString esecopt = "A"; //opzione di esecuzione -> A corrisponde alla preimpostata, in B si possono modificare i parametri
    int kEvents = 1000000; //numero di eventi standard

    //Si permette all'user di scegliere tra un'esecuzione veloce, dove è tutto preimpostato e un'esecuzione attiva 
    cout<<" \n -> Si scriva A per un'esecuzione veloce e preimpostata \n";
    cout<<" \n -> Si scriva B per modificare i parametri \n";
    cout<<" \n Scelta: ";
    
    do {
        cin>>esecopt;
    } while (!esecopt.Contains("a") && !esecopt.Contains("A") && !esecopt.Contains("b") && !esecopt.Contains("B"));

    if (esecopt.Contains("a")) esecopt = "A";
    if (esecopt.Contains("b")) esecopt = "B";

    if (esecopt.Contains("B")) {

        cout<<" \n \n SCELTA SCATTERING MULTIPLO: \n";
        cout<<" A -> Acceso \n";
        cout<<" B -> Spento \n";
        cout<<" Scelta: ";
        do {
            cin>>scattopt;
        } while (!scattopt.Contains("a") && !scattopt.Contains("A") && !scattopt.Contains("b") && !scattopt.Contains("B"));

        if (scattopt.Contains("a")) scattopt = "A";
        if (scattopt.Contains("b")) scattopt = "B";
        
        cout<<" \n \n SCELTA MOLTEPLICITA': \n";
        cout<<" A -> Distribuzione ASSEGNATA \n";
        cout<<" B -> Distribuzione UNIFORME in un RANGE a scelta \n";
        cout<<" C -> Valore FISSO \n";
        cout<<" Scelta: ";
        do {
            cin>>multopt;
        } while (!multopt.Contains("a") && !multopt.Contains("A") && !multopt.Contains("b") && !multopt.Contains("B") && !multopt.Contains("C") && !multopt.Contains("c"));

        if (multopt.Contains("a")) multopt = "A";
        if (multopt.Contains("b")) multopt = "B";
        if (multopt.Contains("c")) multopt = "C";

        cout<<" \n \n SCELTA PSEUDORAPIDITA': \n";
        cout<<" A -> Distribuzione ASSEGNATA tra -2 e 2 \n";
        cout<<" B -> Distribuzione ASSEGNATA tra -6 e 6 \n";
        cout<<" C -> Distribuzione UNIFORME in un RANGE a scelta \n";
        cout<<" Scelta: ";
        do {
            cin>>etaopt;
        } while (!etaopt.Contains("a") && !etaopt.Contains("A") && !etaopt.Contains("b") && !etaopt.Contains("B") && !etaopt.Contains("C") && !etaopt.Contains("c"));

        if (etaopt.Contains("a")) etaopt = "A";
        if (etaopt.Contains("b")) etaopt = "B";
        if (etaopt.Contains("c")) etaopt = "C";

        cout<<" \n \n SCELTA NUMERO DI EVENTI: \n";
        cout<<" Scelta: ";
        cin>>kEvents;

        cout<<" \n \n SCELTA VERTICE: \n";
        cout<<" A -> Z Vertice gaussiana di sigma 53 mm \n";
        cout<<" B -> Z Vertice gaussiana ENTRO un NUMERO di sigma a scelta \n";
        cout<<" C -> Z Vertice Uniforme \n";
        cout<<" Scelta: ";
        do {
            cin>>vertopt;
        } while (!vertopt.Contains("a") && !vertopt.Contains("A") && !vertopt.Contains("b") && !vertopt.Contains("B") && !vertopt.Contains("c") && !vertopt.Contains("C"));

        if (vertopt.Contains("a")) vertopt = "A";
        if (vertopt.Contains("b")) vertopt = "B";
        if (vertopt.Contains("c")) vertopt = "C";

    }

    char command[100];
    sprintf(command,"simulation(\""+scattopt+"\",\""+multopt+"\",\""+etaopt+"\",\""+vertopt+"\",%d,%d)",kEvents,0);
    cout<<" \n COMPILAZIONE CLASSE \n ";
    gSystem->CompileMacro("SimPoint.cxx",opt.Data());
    cout<<" \n COMPILAZIONE  MACRO SIMULAZIONE \n ";
    gSystem->CompileMacro("simulation.C",opt.Data());
    cout<<" \n ESECUZIONE MACRO SIMULAZIONE \n ";
    gROOT->ProcessLine(command);
    cout<<"\n";

    //------------------------------
    //RICOSTRUZIONE
    //------------------------------

    cout<<" \n ------------------------------- \n RICOSTRUZIONE \n ------------------------------- \n ";

    int noise=10; //Scelta std

    do {
        cout<<" -> Si scelga un NUMERO di PUNTI SPURI POSITIVO: ";
        cin>>noise;
    } while(noise<0);

    cout<<" \n COMPILAZIONE MACRO RICOSTRUZIONE \n ";
    gSystem->CompileMacro("reco.C",opt.Data());
    char rec[100];
    sprintf(rec,"reco(%d,%d)",0,noise);
    cout<<" \n ESECUZIONE MACRO RICOSTRUZIONE \n ";
    gROOT->ProcessLine(rec);

    cout<<" \n ------------------------------- \n ANALISI BONTA' RICOSTRUZIONE \n ------------------------------- \n ";
    cout<<"\n ESECUZIONE MACRO ANALISI \n";
    gROOT->ProcessLine(".x treeAnalysis.C");
    gROOT->ProcessLine(".q");


}