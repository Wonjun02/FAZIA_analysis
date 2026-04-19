void fit()
{
  gStyle->SetOptStat( 0 );

  auto file = new TFile( "res_60MeV.root", "READ" );
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
    {
      auto hist = (TH1D*) file->Get( Form("histA_q%d_t%d", i+1, j+1) );
      if( hist->GetEntries() < 1 ) continue;

      hist->SetTitle( Form("Qua: %d Tel: %d", i+1, j+1) );
      hist->GetXaxis()->SetTitle( "A" );
      hist->GetYaxis()->SetTitle( "dN/dA" );

      hist->Scale( 1./hist->GetEntries()/hist->GetBinWidth(1) );

      auto c = new TCanvas();
      hist->Draw( "HIST" );

      auto f = new TF1( "f", "gaus(0) + gaus(3)", -5, 5 );
      f->SetParameter( 0, 1 );
      f->SetParameter( 1, 0 );
      f->SetParameter( 2, 1 );
      f->SetParameter( 3, 1 );
      f->SetParameter( 4, 2 );
      f->SetParameter( 5, 1 );

      f->SetParLimits( 0, 0, 2 );
      f->SetParLimits( 3, 0, 2 );

      f->SetParLimits( 4, 1.9, 2.2 );

      // Q1 T3 don't need fit
      if( i==0 && j==2 )
        continue;

      if( i==0 && j==1 )
      {
        f->SetParLimits( 0, 0, 2 );
        f->SetParLimits( 3, 0, 1 );
      }

      if( i==2 && j==1 )
      {
        f->SetParLimits( 0, 0, 5 );
        f->SetParLimits( 3, 0, 5 );
      }


      hist->Fit( f, "", "", 0.5, 2.5 );
      f->Draw( "same" );

      double mean1 = f->GetParameter( 1 );
      double mean2 = f->GetParameter( 4 );
      double dev1 = std::abs(f->GetParameter( 2 ));// * 2.355;
      double dev2 = std::abs(f->GetParameter( 5 ));// * 2.355;

      double diff = std::abs( mean1 - mean2 );
      double fom = diff / (dev1+dev2);

      cout << "Qua: " << i+1 << " Tel: " << j+1 << " fom: " << fom << endl;

      c->SaveAs( Form("plots/hist_A_%d_%d.png", i+1, j+1) );
    }

}
