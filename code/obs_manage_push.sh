#!/bin/bash
# obs_manage_push.sh - 推送观测目标报告到群聊（飞书卡片 + 表格格式）
# Scheduled: 07:30 UTC (15:30 Beijing time) daily

LOG_FILE="/home/liangrd/Follow_up/code/log/$(date +%Y%m%d)_obs_manage_push.log"
REPORT_FILE="/tmp/obs_report_$(date +%Y%m%d).md"
CHAT_ID="oc_34de0318721ab375e02e2afd4ccd495b"
APP_ID="cli_a9138fa6e4b8dcc2"
APP_SECRET="KmODEuflvo0eJYN0BgIQWcwKeEDXg7lL"

{
    echo "=== $(date) === obs_manage push run ==="

    if [ ! -f "$REPORT_FILE" ]; then
        echo "ERROR: report file not found: $REPORT_FILE"
        echo "=== $(date) failed ==="
        exit 1
    fi

    TOKEN=$(curl -s -X POST https://open.feishu.cn/open-apis/auth/v3/tenant_access_token/internal \
        -H "Content-Type: application/json" \
        -d "{\"app_id\":\"$APP_ID\",\"app_secret\":\"$APP_SECRET\"}" | \
        python3 -c "import sys,json; print(json.load(sys.stdin).get('tenant_access_token',''))")

    if [ -z "$TOKEN" ]; then
        echo "ERROR: Failed to get tenant_access_token"
        echo "=== $(date) failed ==="
        exit 1
    fi

    export TOKEN REPORT_FILE CHAT_ID

    python3 << 'PYEOF'
import json, os, re, subprocess, sys

token = os.environ["TOKEN"]
report_file = os.environ["REPORT_FILE"]
chat_id = os.environ["CHAT_ID"]

with open(report_file, "r") as f:
    content = f.read()

lines = content.strip().split("\n")

# 提取头部信息
info = {}
for line in lines:
    m = re.match(r"📅\s+(.+)", line)
    if m: info["date"] = m.group(1).strip()
    m = re.match(r"🌙\s+月相:\s*(\S+)\s*\|\s*月距阈值:\s*(\S+)\s*\|\s*夜间最大高度角阈值:\s*(\S+)", line)
    if m: info["moon"] = m.group(1); info["dist_thresh"] = m.group(2); info["elev_thresh"] = m.group(3)
    m = re.match(r"共\s*\*{0,2}(\d+)\*{0,2}\s*个可观测目标", line)
    if m: info["count"] = m.group(1)

# 解析表格数据
rows = []
in_table = False
for line in lines:
    if line.startswith("| EP Name"):
        in_table = True; continue
    if in_table and line.startswith("|---"):
        continue
    if in_table and line.startswith("| "):
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) == 8:
            rows.append({
                "name": cells[0].replace(" 🔗", ""),
                "obs": cells[1], "sx": cells[2], "pri": cells[3],
                "moon": cells[4],
                "tnot": cells[5], "sitian": cells[6], "whut": cells[7],
            })

# ---------- 辅助函数 ----------
# 构建卡片目标条目 (无表格)
def build_target_elements(rows):
    elems = []
    for i, r in enumerate(rows):
        badge = "🔴" if r["pri"] == "3" else ("🟡" if r["pri"] == "2" else "🟢")
        elems.append({
            "tag": "div",
            "text": {
                "tag": "lark_md",
                "content": f"{badge} **{r['name']}**  |  Sx: {r['sx']}  |  Priority: {r['pri']}"
            }
        })
        elems.append({
            "tag": "div",
            "text": {
                "tag": "lark_md",
                "content": f"　　月距: {r['moon']}  |  TNOT: {r['tnot']}  |  Sitian: {r['sitian']}  |  WHUT: {r['whut']}"
            }
        })
        if i < len(rows) - 1:
            elems.append({"tag": "hr"})
    return elems

# ---------- 构建卡片 ----------
moon_phase = float(info.get("moon", "0").replace("%", ""))
moon_icon = "🌕" if moon_phase > 90 else ("🌔" if moon_phase > 65 else ("🌓" if moon_phase > 35 else ("🌒" if moon_phase > 10 else "🌑")))

card = {
    "config": {"wide_screen_mode": True},
    "header": {
        "title": {"tag": "plain_text", "content": f"🔭 观测目标报告  |  {info.get('date', '')}"},
        "template": "indigo"
    },
    "elements": [
        {
            "tag": "div",
            "text": {
                "tag": "lark_md",
                "content": f"{moon_icon} 月相 **{info.get('moon', '?')}**  |  月距阈值 **{info.get('dist_thresh', '?')}**  |  高度角阈值 **{info.get('elev_thresh', '?')}**"
            }
        },
        {"tag": "hr"},
        *build_target_elements(rows),
        {"tag": "hr"},
    ]
}
card["elements"].append({
    "tag": "note",
    "elements": [{"tag": "plain_text", "content": f"共 {len(rows)} 个目标  ·  🤖 obs_assistant"}]
})

# ---------- 发送 ----------
payload = {
    "receive_id": chat_id,
    "msg_type": "interactive",
    "content": json.dumps(card, ensure_ascii=False)
}

result = subprocess.run([
    "curl", "-s", "-X", "POST",
    "https://open.feishu.cn/open-apis/im/v1/messages?receive_id_type=chat_id",
    "-H", f"Authorization: Bearer {token}",
    "-H", "Content-Type: application/json",
    "-d", json.dumps(payload, ensure_ascii=False)
], capture_output=True, text=True)

resp = json.loads(result.stdout) if result.stdout else {}
code = resp.get("code", -1)

if code == 0:
    msg_id = resp.get("data", {}).get("message_id", "?")
    print(f"Push success: message_id={msg_id}")
else:
    print(f"Push failed: code={code}, msg={resp.get('msg','unknown')}")
    print(f"Response: {result.stdout[:500]}")
    sys.exit(1)
PYEOF

    echo "push exit: $?"
    echo "=== $(date) done ==="
} >> "$LOG_FILE" 2>&1
