#!/bin/bash

westlog=$1
iter=$(cat $westlog | grep "Beginning iteration" | tail -n 1 | awk '{print $3}')
numsegs=$(cat $westlog | grep "segments remain in iteration $iter" | awk '{print $1}') 
iterdir=$(printf "%0*d\n" 6 $iter)

progress=0

for i in $( eval echo {0..$((numsegs-1))} ); do
	segdir=$(printf "%0*d\n" 6 $i)
	logpath="traj_segs/$iterdir/$segdir/seg.log"
	if [ -f "$logpath" ]; then
	       if grep -q "Final Performance Info" $logpath; then
		       progress=$(expr $progress + 1)
	       fi
	else
		:
	fi
done

if [ $progress -eq 0 ]; then
    echo "Still waiting for the first segment to finish... please check back in a bit."
    exit 1
fi

current_time=$(date | awk '{print $4}')
start_time=$(cat $westlog | grep "2022" | tail -n 1 | awk '{print $4}')
difference=$(( $(date -d "$current_time" "+%s") - $(date -d "$start_time" "+%s") ))

rate=$(bc -l <<< "scale=2; $difference/$progress")

total_time=$(bc -l <<< "scale=2; $rate*$numsegs")

final_time=$(date '+%H:%M:%S' --date="$start_time + $total_time seconds")

waittime=$(( $(date -d "$final_time" "+%s") - $(date -d "$current_time" "+%s") ))

waittime=$(bc -l <<< "scale=1; $waittime/60")

echo "$waittime minutes remain in iteration $iter ($progress/$numsegs segments completed)"
