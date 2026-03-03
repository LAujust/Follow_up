#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <github_username> <repo_name> <github_token>"
  exit 1
fi

USER_NAME="$1"
REPO_NAME="$2"
TOKEN="$3"

API="https://api.github.com/user/repos"

curl -sS -H "Authorization: token ${TOKEN}" -H "Accept: application/vnd.github+json" \
  -d "{\"name\":\"${REPO_NAME}\",\"private\":false}" \
  "$API" >/tmp/create_repo_resp.json

if rg -q '"full_name"' /tmp/create_repo_resp.json; then
  echo "GitHub repo created: ${USER_NAME}/${REPO_NAME}"
else
  echo "Repo creation response:" && cat /tmp/create_repo_resp.json
fi

if git remote get-url origin >/dev/null 2>&1; then
  git remote remove origin
fi

git remote add origin "https://${TOKEN}@github.com/${USER_NAME}/${REPO_NAME}.git"
git push -u origin main

echo "Pushed to https://github.com/${USER_NAME}/${REPO_NAME}"
