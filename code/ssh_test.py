import subprocess

fname = '/Volumes/T7/Shared_Files/EP/Results/Follow_up/wxtsource/wxt_candidates_api_20251208.csv'
password = 'nsqh.800@59726355'

subprocess.run([
    "sshpass", "-p", password,
    "scp", "-P", "5905",
    fname,
    "tnot@119.78.162.172:/home/tnot/EP/plans"
], check=True)