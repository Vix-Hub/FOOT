#using fedra methods to build the EdbTrackP list from the file
import ROOT
import fedrarootlogon
import sys
import time
#usage: python -i VerteTrackDisplay.py -f vertexfile -t tracksfile -nt trackIDS -nv vertexIDS


def close():
  '''Close canvas, allowing program to exit without crashing ROOT'''
  global ds
  ds = 0 
from argparse import ArgumentParser 

dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()
scancond = ROOT.EdbScanCond()
gAli.SetScanCond(scancond)

fedratrackslist = []
#vertexnumberlist = [10, 20]
isolatedtrackcolors = [ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kMagenta, ROOT.kMagenta, ROOT.kMagenta, ROOT.kMagenta] #so we can set different colors for different tracks
vertextrackcolors = [ROOT.kRed,ROOT.kBlue,ROOT.kRed, ROOT.kGreen, ROOT.kCyan, ROOT.kGray] #so we can set different colors for different tracks
#list of possible options
parser = ArgumentParser()
parser.add_argument("-t", "--tracks", dest="tracksfilename", help="file with fedra tracks",
                    required=True)
parser.add_argument("-fsegplates", "--firstsegsplates", dest="seg_plates", help="first seg plates", required=True, nargs='+')
parser.add_argument("-fsegids", "--firstsegids", dest="seg_ids", help="first seg ids ", required = True, nargs='+')

options = parser.parse_args()
tracksfilename = options.tracksfilename
seg_plates = options.seg_plates
seg_ids = options.seg_ids

inputfile = ROOT.TFile.Open(tracksfilename,"READ")
tracktree = inputfile.Get("tracks")
tracktree.SetAlias("trk","t.") #points create confusion to python

tracktree.BuildIndex("s[0].Plate()", "s[0].ID()")
tracks = ROOT.TObjArray(100)
for i in range(len(seg_plates)):
 tracktree.GetEntryWithIndex(int(seg_plates[i]), int(seg_ids[i]))
 #temporary object for reading the file and building EdbTrackP
 temptrack = ROOT.EdbTrackP()
 temptrack.Copy(ROOT.EdbTrackP(tracktree.trk))
 segments = tracktree.s
 #fittedsegments = tracktree.sf
 #loop on segments associated to the track
 for seg in segments:
     seg.SetDZ(300)
     seg.SetW(90)
     temptrack.AddSegment(seg)
     #temptrack .AddSegmentF(segf)
     temptrack.SetSegmentsTrack(temptrack.ID())
     temptrack.SetCounters()
     print(" Added segment " + str(seg.ID()) + " " + str(seg.Plate()))
 mytrack = ROOT.EdbTrackP()
 mytrack.Copy(temptrack)
 print(" Adding EdbTrackP with N " + str(mytrack.N()) + " segments ")
 tracks.Add(temptrack)


ROOT.gStyle.SetPalette(1);

dsname="FOOT FEDRA display"
ds = ROOT.EdbDisplay.EdbDisplayExist(dsname);
if not ds:  
  ds=ROOT.EdbDisplay(dsname,-50000.,50000.,-50000.,50000.,-4000.,80000.)


ds.SetDrawTracks(4)
ds.SetArrTr(tracks)
ds.Draw()
time.sleep(100)
close()


