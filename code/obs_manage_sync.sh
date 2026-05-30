#!/bin/bash
# obs_manage sync + filter (lunar distance & elevation angle)
# Scheduled: 07:20 UTC (15:20 Beijing time) daily

export PATH="/home/liangrd/.nvm/versions/node/v24.13.1/bin:$PATH"

LOG_DIR="/home/liangrd/Follow_up/code/log"
mkdir -p "$LOG_DIR"
LOGFILE="$LOG_DIR/$(date +%Y%m%d)_obs_manage_sync.log"
REPORT_FILE="/tmp/obs_report_$(date +%Y%m%d).md"

{
    echo "=== $(date) === obs_manage sync run ==="

    # Step 1: sync cloud sheet → Candidates.csv
    /home/liangrd/anaconda3/envs/autophot/bin/python \
        /home/liangrd/Follow_up/code/obs_manage.py --sync
    echo "sync exit: $?"

    # Step 2: filter + generate report
    /home/liangrd/anaconda3/envs/autophot/bin/python \
        /home/liangrd/Follow_up/code/obs_manage.py --filter --no-push > "$REPORT_FILE" 2>/dev/null
    echo "filter exit: $?"
    echo "report written to: $REPORT_FILE"

    echo "=== $(date) done ==="
} >> "$LOGFILE" 2>&1
