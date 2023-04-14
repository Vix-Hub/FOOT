{
  gROOT->ProcessLine(".except");
  //gROOT->ProcessLine(".L EMU_Track_c.so");
    gROOT->ProcessLine(".L merge_trackings_SmallSample.C+");
    gROOT->ProcessLine("merge_trackings()");
}
