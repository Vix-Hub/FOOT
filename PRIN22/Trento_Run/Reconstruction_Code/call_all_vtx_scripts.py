import subprocess

# Script 1: remove_duplicates_server_new.py
subprocess.call(["python", "merge_mts_new.py"])

# Script 2: add_Nincoming.py
subprocess.call(["python", "merge_mts_post.py"])

# Script 3: add_charge_fast.py
subprocess.call(["python", "find_interaction_candidates_new.py"])

# Script 4: recover_missing_tracks_new.py (step 1)
subprocess.call(["python", "check_candidates.py"])

