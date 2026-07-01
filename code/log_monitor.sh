#!/bin/bash
# log_monitor.sh — 监视 /home/liangrd/Follow_up/code/log 的新增错误
# 每小时一次。首次运行初始化状态（不告警），后续只扫描增量变化。
set -uo pipefail

export PATH="/home/liangrd/.nvm/versions/node/v24.13.1/bin:$PATH"

LOG_DIR="/home/liangrd/Follow_up/code/log"
STATE_DIR="/home/liangrd/Follow_up/code/.log_state"
CHAT_ID="oc_34de0318721ab375e02e2afd4ccd495b"
MONITOR_LOG="$LOG_DIR/monitor.log"

mkdir -p "$STATE_DIR" "$LOG_DIR"

# ── 发送告警 ──────────────────────────────────────────────
ALERTS_TMP=$(mktemp)

send_alerts() {
    if [ ! -s "$ALERTS_TMP" ]; then
        return 0
    fi
    local msg
    msg=$(cat "$ALERTS_TMP")
    lark-cli im +messages-send --as user --chat-id "$CHAT_ID" --text "$msg" 2>/dev/null || true
}

add_alert() {
    local fname="$1"
    local errors="$2"
    {
        echo "📄 **\`$fname\`**"
        echo "\`\`\`"
        printf '%s\n' "$errors"
        echo "\`\`\`"
        echo ""
    } >> "$ALERTS_TMP"
}

# ── 提取关键错误行 ────────────────────────────────────────
extract_errors() {
    local text="$1"
    echo "$text" | grep -iE '(ERROR|Traceback|FileNotFoundError|EmptyDataError|ModuleNotFoundError|KeyError|error:|Exception|exit: [1-9]|Failed)' | head -10 || true
}

# ── 检查单个文件的新增内容 ────────────────────────────────
check_file() {
    local fname="$1"
    local fpath="$LOG_DIR/$fname"
    local state_file="$STATE_DIR/$fname.size"

    [ ! -f "$fpath" ] && return 0
    [ "$fname" = "monitor.log" ] && return 0

    if [ ! -f "$state_file" ]; then
        # 新文件：初始化状态（不告警）
        stat -c%s "$fpath" 2>/dev/null > "$state_file" || echo 0 > "$state_file"
        echo "  [$fname] NEW — initialized (${cur_size:-0}B)"
        return 0
    fi

    local cur_size prev_size
    cur_size=$(stat -c%s "$fpath" 2>/dev/null || echo 0)
    prev_size=$(cat "$state_file" 2>/dev/null || echo 0)

    if [ "$cur_size" -gt "$prev_size" ] 2>/dev/null; then
        local delta=$((cur_size - prev_size))
        local new_content
        new_content=$(tail -c "$delta" "$fpath" 2>/dev/null || true)

        if [ -n "$new_content" ]; then
            local errors
            errors=$(extract_errors "$new_content")
            if [ -n "$errors" ]; then
                echo "  [$fname] +${delta}B ⚠ errors"
                add_alert "$fname" "$errors"
            else
                echo "  [$fname] +${delta}B ✅ clean"
            fi
        fi

        echo "$cur_size" > "$state_file"
    fi
}

# ── 主流程 ────────────────────────────────────────────────
{
    echo "=== $(date -u) UTC ==="

    init_count=0
    while IFS= read -r -d '' fpath; do
        fname=$(basename "$fpath")
        [ "$fname" = "monitor.log" ] && continue
        state_file="$STATE_DIR/$fname.size"
        if [ ! -f "$state_file" ]; then
            stat -c%s "$fpath" 2>/dev/null > "$state_file" || echo 0 > "$state_file"
            init_count=$((init_count + 1))
        fi
    done < <(find "$LOG_DIR" -maxdepth 1 -type f -name "*.log" -print0 2>/dev/null | sort -z)
    [ "$init_count" -gt 0 ] && echo "  Initialized $init_count new log(s)"

    # Phase 2: 检查增量变化
    while IFS= read -r -d '' fpath; do
        fname=$(basename "$fpath")
        check_file "$fname"
    done < <(find "$LOG_DIR" -maxdepth 1 -type f -name "*.log" -print0 2>/dev/null | sort -z)

    send_alerts
    rm -f "$ALERTS_TMP"
    echo "=== $(date -u) UTC done ==="
} >> "$MONITOR_LOG" 2>&1