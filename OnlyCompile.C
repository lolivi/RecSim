
void OnlyCompile(TString myopt = "fast"){
    TString opt; //compila o con il + (fast) o con il ++ (force)
    if (myopt.Contains("force")) {
        opt = "kfg";
    }
    else {
        opt = "kg";
    }

    cout<<endl<<"COMPILAZIONE DELLA CLASSE"<<endl;
    gSystem->CompileMacro("SimPoint.cxx",opt.Data());

    cout<<endl<<"COMPILAZIONE DELLA MACRO DI SIMULAZIONE"<<endl;
    gSystem->CompileMacro("simulation.C",opt.Data());
    
    cout<<endl<<"COMPILAZIONE DELLA MACRO DI RICOSTRUZIONE"<<endl;
    gSystem->CompileMacro("reco.C",opt.Data());
    cout<<endl;
}