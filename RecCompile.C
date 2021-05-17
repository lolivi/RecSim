
void RecCompile (TString myopt = "fast") {

    TString opt; //compila o con il + (fast) o con il ++ (force)
    if (myopt.Contains("force")) {
        opt = "kfg";
    }
    else {
        opt = "kg";
    }

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
    cout<<" \n ESECUZIONE MACRO RICOSTRUZIONE \n";
    gROOT->ProcessLine(rec);

    cout<<" \n ------------------------------- \n ANALISI BONTA' RICOSTRUZIONE \n ------------------------------- \n ";
    cout<<"\n ESECUZIONE MACRO ANALISI \n";
    gROOT->ProcessLine(".x treeAnalysis.C");
    gROOT->ProcessLine(".q");
    


}