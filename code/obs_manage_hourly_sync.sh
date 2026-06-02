#!/bin/bash
# obs_manage_hourly_sync.sh — 每小时将云文档 Targets 表格同步到 Candidates.csv
# Scheduled: 0 * * * * (UTC, 每小时整点)

LOG_DIR="/home/liangrd/Follow_up/code/log"
mkdir -p "$LOG_DIR"
LOGFILE="$LOG_DIR/$(date +%Y%m%d)_obs_manage_hourly_sync.log"

{
    echo "=== $(date -u) UTC ==="

    /home/liangrd/anaconda3/envs/autophot/bin/python \
        /home/liangrd/Follow_up/code/obs_manage.py --sync
    echo "sync exit: $?"

    echo "=== $(date -u) UTC done ==="
} >> "$LOGFILE" 2>&1
