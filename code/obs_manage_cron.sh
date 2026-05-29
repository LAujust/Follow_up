#!/bin/bash
# obs_manage daily cron wrapper — sync + agent-relayed push
# Scheduled: 07:30 UTC (15:30 Beijing time) daily

LOG_DIR="/home/liangrd/Follow_up/code/log"
mkdir -p "$LOG_DIR"
LOGFILE="$LOG_DIR/$(date +%Y%m%d)_obs_manage.log"
REPORT_FILE="/tmp/obs_report_$(date +%Y%m%d).md"

{
    echo "=== $(date) === obs_manage daily run ==="

    # Step 1: sync cloud sheet → Candidates.csv
    /home/liangrd/anaconda3/envs/autophot/bin/python \
        /home/liangrd/Follow_up/code/obs_manage.py --sync

    # Step 2: filter + generate report (no push via lark-cli)
    /home/liangrd/anaconda3/envs/autophot/bin/python \
        /home/liangrd/Follow_up/code/obs_manage.py --filter --no-push > "$REPORT_FILE" 2>/dev/null

    # Step 3: delegate push to obs_assistant agent
    openclaw agent \
        --agent obs_assistant \
        --message "请读取文件 $REPORT_FILE 的内容，原样发送到群聊 oc_34de0318721ab375e02e2afd4ccd495b" \
        --timeout 120000

    echo "=== $(date) done ==="
} >> "$LOGFILE" 2>&1
