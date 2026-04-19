int ene = 70;
void ana()
{
  // ------------------------ CALIBRATION ------------------------ 
  double siCal[4][4] = {0};
  double csiCal[4][4] = {0};

  for( int i=0; i<4; i++ )
  {
    for( int j=0; j<4; j++ )
    {
      if( i==0 && j==1 ) siCal[i][j] = 0.27091;
      if( i==3 && j==3 ) siCal[i][j] = 0.26171;
      if( i==0 && j==0 ) siCal[i][j] = 0.25190;
      if( i==0 && j==2 ) siCal[i][j] = 0.30027;
      if( i==2 && j==1 ) siCal[i][j] = 0.28011;

      if( i==0 && j==1 ) csiCal[i][j] = 2.644;
      if( i==3 && j==3 ) csiCal[i][j] = 3.219;
      if( i==0 && j==0 ) csiCal[i][j] = 2.618;
      if( i==0 && j==2 ) csiCal[i][j] = 2.851;
      if( i==2 && j==1 ) csiCal[i][j] = 2.174;
    }
  }

  TF1* fcsi[4][4];
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      fcsi[i][j] = new TF1( Form("fcsi_%d_%d", i, j), "[0] * (x*(1. - [1]/x * TMath::Log(1.+ x/[1])) + [3]*[1]*TMath::Log( (x+[1])/([2]+[1]) ))", 0, 80 );

  for( int i=0; i<4; i++ ) for( int j=0; j<4; j++ ) fcsi[i][j]->SetNpx( 1000 );

  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      fcsi[i][j]->SetParameters( csiCal[i][j], 0.25, 3.10, 0.27 );
  // ------------------------ CALIBRATION ------------------------ 




  // ------------------------ SETUP ------------------------ 
  //auto file = new TFile( Form("data_%dMeV_1.root", ene), "READ" );
  auto file = new TFile( Form("data_redo_%dMeV.root", ene), "READ" );

  TF1* f[4][4];
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      f[i][j] = (TF1*) file->Get( Form("fTG_q%d_t%d", i+1, j+1) );

  TF1* ff[4][4];
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
    {
      ff[i][j] = new TF1( Form("ff_q%d_t%d", i, j),
          "pow( pow([0]*[4],[1]+1) + pow([2],[1]+1)*pow([3],2)*pow(x,[1]) , 1./([1]+1) ) - [0]*[4] + [5]",
          0.1, 10 );

      ff[i][j]->SetNpx( 1000 );
      ff[i][j]->SetParameter( 0, f[i][j]->GetParameter(0) );
      ff[i][j]->SetParameter( 1, f[i][j]->GetParameter(1) );
      ff[i][j]->SetParameter( 2, f[i][j]->GetParameter(2) );
      ff[i][j]->SetParameter( 3, f[i][j]->GetParameter(3) );
      //ff[i][j]->SetParameter( 4, f[i][j]->GetParameter(4) );
      ff[i][j]->SetParameter( 5, f[i][j]->GetParameter(5) );
    }

  //for( int i=0; i<4; i++ ) for( int j=0; j<4; j++ ) f[i][j]->SetParameter( 4, 0.5 );



  TChain chain("faziatree");
  int startRun = 0; int endRun = 0;
  if( ene==70 ) { startRun = 340; endRun = 1032; } // 70MeV /day1
  if( ene==60 ) { startRun = 1033; endRun = 1699; } // 60MeV /day1

  for (int run=startRun; run<=endRun; run++) 
    chain.Add(Form("/Volumes/T7/converted_data/run_%d.root", run));



  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("Mtot",   1);
  chain.SetBranchStatus("fQua",   1);
  chain.SetBranchStatus("fTel",   1);
  chain.SetBranchStatus("sQ3max", 1);
  chain.SetBranchStatus("sQ2max", 1);

  chain.SetCacheSize(80*1024*1024);
  chain.AddBranchToCache("Mtot",   true);
  chain.AddBranchToCache("fQua",   true);
  chain.AddBranchToCache("fTel",   true);
  chain.AddBranchToCache("sQ3max", true);
  chain.AddBranchToCache("sQ2max", true);

  UShort_t mtot = 0;
  UShort_t fQua[100] = {0};
  UShort_t fTel[100] = {0};
  float    sQ3max[100] = {0};
  float    sQ2max[100] = {0};

  chain.SetBranchAddress("Mtot",   &mtot);
  chain.SetBranchAddress("fQua",   fQua);
  chain.SetBranchAddress("fTel",   fTel);
  chain.SetBranchAddress("sQ3max", sQ3max);
  chain.SetBranchAddress("sQ2max", sQ2max);


  // ------------------------ SETUP ------------------------ 


TH1D* histA[4][4];
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      histA[i][j] = new TH1D( Form("histA_q%d_t%d", i+1, j+1), "", 100, 0, 10 );

TH1D* histP[4][4];
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      histP[i][j] = new TH1D( Form("histP_q%d_t%d", i+1, j+1), "", 300, -30, 30 );


TH2D* histProj[4][4];
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      histProj[i][j] = new TH2D( Form("histProj_q%d_t%d", i+1, j+1), "", 400, 0, 100, 300, -30, 30 );


  // ------------------------ ANALYSIS ------------------------ 
  Long64_t nEntries = chain.GetEntries();
  //nEntries = nEntries/10;
  std::cout << "Total chained entries = " << nEntries << std::endl;

  for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) 
  {
    chain.GetEntry(iEntry);

    int nHit = (int)mtot;
    if (16 < nHit) continue;

    for (int i = 0; i < nHit; i++) 
    {

      if( !( (fQua[i]==1 && fTel[i]==2)
            || (fQua[i]==4 && fTel[i]==4)
            || (fQua[i]==1 && fTel[i]==1)
            || (fQua[i]==1 && fTel[i]==3)
            || (fQua[i]==3 && fTel[i]==2)
           ) ) continue;

        int q = fQua[i] - 1;
        int t = fTel[i] - 1;


        double cal1 = siCal[q][t];

        double x = 0;
        if( !(1<sQ3max[i]) ) continue;

        x = fcsi[q][t]->GetX( sQ3max[i], 1e-3, sQ3max[i] );
        double y = sQ2max[i] * cal1;


        if( q==0 && t==0 && x<10 ) continue;
        if( q==0 && t==1 && x<4 ) continue;
        if( q==0 && t==2 && x<8 ) continue;
        if( q==2 && t==1 && x<10 ) continue;
        if( q==3 && t==3 && x<8 ) continue;


        f[q][t]->SetParameter( 4, 0.5 );
        if( y < f[q][t]->Eval(x) ) continue;


        f[q][t]->SetParameter( 4, 1. );
        double diff = y - f[q][t]->Eval( x );
        histProj[q][t]->Fill( x, diff );
        histP[q][t]->Fill( diff );






//        continue;
        ff[q][t]->SetParameter( 4, x );

        double aa = ff[q][t]->GetX( y );

        histA[q][t]->Fill( aa );

    }
  }
  // ------------------------ ANALYSIS ------------------------ 


  /*
     TH2D* hist[4][4];
     for( int i=0; i<4; i++ )
     for( int j=0; j<4; j++ )
     hist[i][j] = (TH2D*) file->Get( Form("histCorr_q%d_t%d", i+1, j+1) );


   */

     auto outputFile = new TFile( Form( "res_%dMeV.root", ene), "RECREATE" );
     for( int i=0; i<4; i++ )
       for( int j=0; j<4; j++ )
         histProj[i][j]->Write();
     for( int i=0; i<4; i++ )
       for( int j=0; j<4; j++ )
         histP[i][j]->Write();
     for( int i=0; i<4; i++ )
       for( int j=0; j<4; j++ )
         histA[i][j]->Write();/f


         histA[3][3]->Draw();

}
