void draw()
{
  int ene = 70;

    int padMap[16][2] =
  {
    {1,1}, {1,2}, {2,1}, {2,2}, // row1
    {1,4}, {1,3}, {2,4}, {2,3}, // row2
    {4,1}, {4,2}, {3,1}, {3,2}, // row3
    {4,4}, {4,3}, {3,4}, {3,3} // row4
  };


  auto c = new TCanvas( "c", "", 1600, 1600 );
  c->Divide( 4, 4 );
  auto file = new TFile( Form("histogram_corr_%dMeV.root", ene), "READ" );
  for (int pad = 0; pad <= 15; pad++)//
  {
    int q = padMap[pad][0]-1;
    int t = padMap[pad][1]-1;
    c->cd( pad+1 )->SetLogz();
    auto h = (TH2D*) file->Get( Form( "histCorr_q%d_t%d", q+1, t+1) );
    h->Draw( "colz" );
  }

}
