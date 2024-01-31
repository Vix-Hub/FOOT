Instructions for CONDOR usage

Copy dm_tracks2.dm.root to /eos/user/v/viboccia/PRIN22/PLATE
python3 save_ausiliary_info.py from /eos/
note the total number of views

Change number of views in /afs/.../condor_merge_mts.sh
From /afs/cern.ch/work/v/viboccia/public/condor_submissions/PRIN22/PLATE -> condor_submit condor_merge_mts.sub

From /eos/... root -l merge_files.C
