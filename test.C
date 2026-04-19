void test()
{
  int qSel = 3;
  int tSel = 2;

  gStyle->SetOptStat(0);

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




  const int    nbin = 400;
  const double xmin = 0.0, xmax = 100.0;
  const double ymin = 0.0, ymax = 20.0;

  // u slice settings
  const double xMinArc = 0.0;
  const double xMaxArc = 200.0;
  const double sliceW  = 10.0;
  //int nSlice = (int)((xMaxArc - xMinArc)/sliceW); // 20
  const int nSlice = 20;

  TChain chain("faziatree");
  //int startRun = 340, endRun = 1032; // 70MeV /day1
  int startRun = 1033, endRun = 1699; // 60MeV /day1

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

  TH2D* histCorr[4][4];
  TH2D* histdEE[4][4];
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
    {
      qSel = i+1;
      tSel = j+1;
      histCorr[i][j] = new TH2D(Form("histCorr_q%d_t%d", qSel, tSel),
          Form("Qua %d Tel %d; E_{CsI} (MeV); E_{Si2} (MeV)", qSel, tSel),
          nbin, xmin, xmax,
          nbin, ymin, ymax);


      histdEE[i][j] = new TH2D(Form("histdEE_q%d_t%d", qSel, tSel),
          Form("Qua %d Tel %d; E_{CsI} + E_{Si2} (MeV); E_{Si2} (MeV)", qSel, tSel),
          nbin, xmin, xmax,
          nbin, ymin, ymax);
    }


  // 5) ONE PASS over chain: fill 2D + store points
  Long64_t nEntries = chain.GetEntries();
  std::cout << "Total chained entries = " << nEntries << std::endl;

  for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
    chain.GetEntry(iEntry);

    int nHit = (int)mtot;
    if (nHit > 100) nHit = 100;

    for (int i = 0; i < nHit; i++) {

      if( !( (fQua[i]==1 && fTel[i]==2)
          || (fQua[i]==4 && fTel[i]==4)
          || (fQua[i]==1 && fTel[i]==1)
          || (fQua[i]==1 && fTel[i]==3)
          || (fQua[i]==3 && fTel[i]==2)
        ) ) continue;

      int q = fQua[i] - 1;
      int t = fTel[i] - 1;

      double cal1 = siCal[q][t];

      float x = 0;
      if( 1<sQ3max[i] ) x = fcsi[q][t]->GetX( sQ3max[i], 1e-3, sQ3max[i] );
      float y = sQ2max[i] * cal1;

      histCorr[q][t]->Fill(x, y);
      //histdEE[q][t]->Fill(x+y, y);
    }
  }


  auto file = new TFile( "histogram_corr_60MeV.root", "RECREATE" );
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      histCorr[i][j]->Write();
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      histdEE[i][j]->Write();
}
