  if (strstr(options, "cutwmass") != NULL){
    sscanf(options, "cutwmass=%lf", &cutwmass);
    if (print==1) cout << "  cutwmass = " <<  cutwmass << endl;
  }
  if (strstr(options, "leptonSep") != NULL){
    sscanf(options, "leptonSep=%d", &leptonSep);
    if (print==1) cout << "  leptonSep = " << leptonSep << endl;
  }
        
  if (strstr(options, "divideSample") != NULL){  //samor 14/06/2010
    sscanf(options, "divideSample=%d", &divideSample);
    if (print==1) cout << "  divideSample = " << divideSample << endl;
  }
  if (strstr(options, "lepSample") != NULL){
    sscanf(options, "lepSample=%d", &lepSample);
    if (print==1) cout << "  lepSample = " << lepSample << endl;
  }
  if (strstr(options, "alldata_ele") != NULL){
    sscanf(options, "alldata_ele=%d", &alldata_ele);
    if (print==1) cout << "  alldata_ele = " << alldata_ele << endl;
  }
  if (strstr(options, "alldata_muo") != NULL){
    sscanf(options, "alldata_muo=%d", &alldata_muo);
    if (print==1) cout << "  alldata_muo = " << alldata_muo << endl;
  }
  if (strstr(options, "pileup") != NULL){
    sscanf(options, "pileup=%d", &pileup);
    if (print==1) cout << "  pileup = " << pileup << endl;
  }
  if (strstr(options, "CutTriggerEle") != NULL){
    sscanf(options, "CutTriggerEle=%d", &CutTriggerEle);
    if (print==1) cout << "  CutTriggerEle = " <<  CutTriggerEle << endl;
  }
  if (strstr(options, "CutTriggerMuo") != NULL){
    sscanf(options, "CutTriggerMuo=%d", &CutTriggerMuo);
    if (print==1) cout << "  CutTriggerMuo = " <<  CutTriggerMuo << endl;
  }
  if (strstr(options, "CutTriggerEleMuo") != NULL){
    sscanf(options, "CutTriggerEleMuo=%d", &CutTriggerEleMuo);
    if (print==1) cout << "  CutTriggerEleMuo = " <<  CutTriggerEleMuo << endl;
  }
  if (strstr(options, "cuttopmass") != NULL){
    sscanf(options, "cuttopmass=%lf", &cuttopmass);
    if (print==1) cout << "  cuttopmass = " << cuttopmass << endl;
  }
  if (strstr(options, "cutptjet") != NULL){
    sscanf(options, "cutptjet=%lf", &cutptjet);
    if (print==1) cout << "  cutptjet = " <<  cutptjet << endl;
  }
  if (strstr(options, "signalpdf") != NULL){
    sscanf(options, "signalpdf=%s", signalpdf);
    if (print==1) cout << "  signalpdf = " <<  signalpdf << endl;
  }
  if (strstr(options, "corrlum") != NULL){
    sscanf(options, "corrlum=%lf", &corrlum);
    corrlum=1./corrlum;
    if (print==1) cout << "  corrlum = " <<  corrlum << endl;
  }
  if (strstr(options, "corrnorm_wjets") != NULL){
    sscanf(options, "corrnorm_wjets=%lf", &corrnorm_wjets);
    if (print==1) cout << "  corrnorm_wjets = " <<  corrnorm_wjets << endl;
  }

   if (strstr(options, "corrnorm_tt") != NULL){
    sscanf(options, "corrnorm_tt=%lf", &corrnorm_tt);
    if (print==1) cout << "  corrnorm_tt = " <<  corrnorm_tt << endl;
  }

  if (strstr(options, "corrnorm_wz") != NULL){
    sscanf(options, "corrnorm_wz=%lf", &corrnorm_wz);
    if (print==1) cout << "  corrnorm_wz = " <<  corrnorm_wz << endl;
  }
  if (strstr(options, "corrnorm_st") != NULL){
    sscanf(options, "corrnorm_st=%lf", &corrnorm_st);
    if (print==1) cout << "  corrnorm_st = " <<  corrnorm_st << endl;
  }
  if (strstr(options, "cutchi2") != NULL){
    sscanf(options, "cutchi2=%lf", &cutchi2);
    if (print==1) cout << "  cutchi2 = " <<  cutchi2 << endl;
  }
  if (strstr(options, "ref_weight") != NULL){
    sscanf(options, "ref_weight=%lf", &ref_weight);
    if (print==1) cout << "  ref_weight = " <<  ref_weight << endl;
  }

