void redo_PID()
{
  gStyle->SetOptStat(0);

  int En = 60;
  //int Energies[] = {70, 60, 50, 40};
  int qSel = 1;
  int tSel = 2;

  auto file = new TFile( Form("data_%dMeV_1.root", En), "READ" );

  TF1* fTG[4][4];
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      fTG[i][j] = (TF1*) file->Get( Form("fTG_q%d_t%d", i+1, j+1) );

  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      fTG[i][j]->SetNpx(1000);


  TGraph* g[4][4];
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
      g[i][j] = new TGraph();

  double sig = 1.;
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
    {
      qSel = i+1;
      tSel = j+1;

      if( !( (i==0 && j==1)
          || (i==3 && j==3)
          || (i==0 && j==0)
          || (i==0 && j==2)
          || (i==2 && j==1)
          ) ) continue;

      auto h = (TH2D*) file->Get( Form("histCorr_q%d_t%d", i+1, j+1) );

      int nx = h->GetNbinsX();
      int ny = h->GetNbinsY();

      int np = 20;
      for( int ip=0; ip<np; ip++ )
      {
        int binMin = ip*np + 1;
        int binMax = ip*np + np;

        double xx = h->GetXaxis()->GetBinUpEdge( ip*np + np/2 );
        if( En < xx ) continue;
        double fakey = fTG[i][j]->Eval( xx );

        auto hp = h->ProjectionY( Form("hproj_%d_%d_%d", i, j, ip), binMin, binMax );
        hp->Fit( "gaus", "", "", fakey-sig, fakey+sig );
        auto ff = (TF1*) gROOT->GetFunction( "gaus" );
        double yy = ff->GetParameter( 1 );

        g[i][j]->SetPoint( g[i][j]->GetN(), xx, yy );
      }

      if( 0&& i==0 && j==1 )
      {
        h->Draw( "colz" );
        g[i][j]->SetMarkerStyle( 29 );
        g[i][j]->SetMarkerSize( 2 );
        g[i][j]->SetMarkerColor( kRed );
        g[i][j]->Draw( "Psame" );
        g[i][j]->Fit( fTG[i][j], "", "", 0, En );
        return;
      }
    }

    auto outputFile = new TFile( Form("data_redo_%dMeV.root", En), "RECREATE" );
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
    {
      fTG[i][j]->Write();
      auto h2 = (TH2D*) file->Get( Form("histCorr_q%d_t%d", i+1, j+1) );
      h2->Write();
    }
}
