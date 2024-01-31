#!/bin/bash

ProcId=$2
VIEWID=$3


#source /afs/cern.ch/work/v/viboccia/public/fedra/setup_new.sh	

echo  "go into reconstruction folder "
cd /eos/user/v/viboccia/PRIN22/P018_bot
echo " merge mts from view  " $VIEWID

source merge_mts_view.sh $VIEWID


#to avoid crashes put these in sub
#requirements = Machine =!= LastRemoteHost
#on_exit_remove          = (ExitBySignal == False) && ((ExitCode == 1) || (ExitCode == 0))
#max_retries             = 3