import os
import time
import re
import requests
from pathlib import Path
from typing import Optional
import pandas as pd
from playwright.sync_api import sync_playwright, TimeoutError as PWTimeout
from astropy.time import Time
import numpy as np
from astropy.table import Table, Column
import astropy.units as u
import requests
import subprocess


# ======== 配置 =========
STATE_PATH = "ep_oauth_state.json"     # 复用登录态（建议用 storage_state）
DOWNLOAD_ROOT = "/Volumes/T7/Shared_Files/EP/Results/Follow_up/wxtsource"  # 你的下载根目录
PROTECTED_URL = "https://ep.bao.ac.cn/ep/data_center/source_candidate_list"
API_URL = "https://ep.bao.ac.cn/ep/data_center/wxt_source_candidate_list_api"
HEADLESS = True
PER_PAGE_TIMEOUT = 5 * 60 * 1000       # 每个页面操作最大 5 分钟
ASYNC_MAX_WAIT = 10 * 60               # 等待异步打包最多 10 分钟（秒）
RETRY_TIMES = 2                        # 出错时每行重试次数
T_NOW = Time.now()
TIME_LAG = 36                          #小时
T_START = T_NOW - TIME_LAG * u.hour
yyyy, mm, dd, _, _, _ = T_NOW.ymdhms

transient_path = DOWNLOAD_ROOT

def mkrow_dir(destination_path: Path, obs_id: str, detnam: str, version: str) -> Path:
    d = destination_path / f"{obs_id}_{detnam}_{version}"
    d.mkdir(parents=True, exist_ok=True)
    return d

def wait_async_link(page) -> Optional[str]:
    """
    点击“Download User data”后，等待页面把下载链接写入 #async-download-link a 的 href。
    返回 href（download_url），若超时则返回 None。
    """
    # 触发开始按钮；按钮在源码里是 <button ... onclick="startAsyncDownload(...)" >
    page.get_by_role("button", name="Download User data").click()
    page.wait_for_selector("#async-download-status", state="visible", timeout=30_000)

    # 页面会每 3 秒自动轮询；我们也自定义一个“最长等待”
    t0 = time.time()
    while time.time() - t0 < ASYNC_MAX_WAIT:
        try:
            # 当下载就绪时，容器显示且 a 标签有 href
            page.wait_for_function(
                """() => {
                    const box = document.querySelector('#async-download-link');
                    const a = box ? box.querySelector('a') : null;
                    return box && getComputedStyle(box).display !== 'none' && a && a.href;
                }""",
                timeout=5_000
            )
            href = page.eval_on_selector("#async-download-link a", "el => el.href")
            return href
        except PWTimeout:
            # 也可读一眼状态文本，便于日志排查
            try:
                status = page.inner_text("#async-status-text")
                print("  [async] status =", status)
            except Exception:
                pass
            continue
    return None

def download_user_zip(page, out_dir: Path, obs_id: str, detnam: str, version: str) -> Path:
    """执行异步打包并下载 zip，返回保存路径。"""
    target = out_dir / f"{obs_id}_{detnam}_{version}_user.zip"
    if target.exists():
        print(f"  [skip] user zip exists: {target.name}")
        return target

    # 等待下载链接出现
    href = wait_async_link(page)
    if not href:
        raise RuntimeError("Async user-data link did not appear within time limit.")

    # 下载（a 标签 target=_blank；使用 expect_download 捕获）
    with page.expect_download(timeout=PER_PAGE_TIMEOUT) as dl_info:
        page.click("#async-download-link a")
    dl = dl_info.value
    dl.save_as(target.as_posix())
    print(f"  [ok] user zip -> {target.name}")
    return target

def download_sources_csv(page, out_dir: Path, obs_id: str, detnam: str, version: str) -> Path:
    """点击 Export CSV 下载表格，保存为定制文件名。"""
    target = out_dir / f"{obs_id}_{detnam}_{version}_sources.csv"
    if target.exists():
        print(f"  [skip] sources csv exists: {target.name}")
        return target

    # 按钮 id 固定为 download_src
    page.wait_for_selector("#download_src", timeout=30_000)
    with page.expect_download(timeout=PER_PAGE_TIMEOUT) as dl_info:
        page.click("#download_src")
    dl = dl_info.value
    dl.save_as(target.as_posix())
    print(f"  [ok] sources csv -> {target.name}")
    return target

def process_one_row(page, row: pd.Series, destination_path: Path):
    """处理单行：进入详情页，下载用户数据 zip + 源表 CSV。"""
    obs_id = str(row["obs_id"]).strip()
    detnam = str(row["detnam"]).strip()
    version = str(row.get("version", "02")).strip()  # 你 df 如果有 version 列，就直接取

    url = f"https://ep.bao.ac.cn/ep/data_center/fxt_obs_detail/{obs_id}/{detnam}/{version}"
    print(f"\n==> {obs_id} {detnam} {version}")
    page.goto(url, wait_until="domcontentloaded", timeout=PER_PAGE_TIMEOUT)

    # 确认页面加载到了期望模块（有“Get Data Files”和“Sources in the observation”）
    page.wait_for_selector("text=Get Data Files", timeout=30_000)
    page.wait_for_selector("text=Sources in the observation", timeout=30_000)

    out_dir = mkrow_dir(destination_path, obs_id, detnam, version)

    # 1) Download User data (异步打包 + 下载)
    try:
        download_user_zip(page, out_dir, obs_id, detnam, version)
    except Exception as e:
        print(f"  [warn] user zip failed: {e}")

    # 2) Export CSV（源表）
    try:
        download_sources_csv(page, out_dir, obs_id, detnam, version)
    except Exception as e:
        print(f"  [warn] sources csv failed: {e}")

def run_batch(df: pd.DataFrame, destination_path: Path):
    df_sorted = df

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=HEADLESS)
        ctx = browser.new_context(
            storage_state=STATE_PATH,
            accept_downloads=True
        )
        page = ctx.new_page()

        for i, row in df_sorted.iterrows():
            for attempt in range(1, RETRY_TIMES + 2):
                try:
                    process_one_row(page, row, destination_path)
                    break
                except Exception as e:
                    print(f"  [err] row {i} attempt {attempt} failed: {e}")
                    if attempt <= RETRY_TIMES:
                        time.sleep(3)
                        continue
                    else:
                        print("  [give up] move on.")
                        break

        browser.close()


def login_with_playwright(username, password, headless=True, state_path=STATE_PATH):
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        ctx = browser.new_context()
        page = ctx.new_page()

        # 1) 从受限页开始（若已登录就直接复用）
        page.goto(PROTECTED_URL, wait_until="domcontentloaded")
        if "oauth.china-vo.org" not in page.url:
            ctx.storage_state(path=state_path)
            cookies = ctx.cookies()
            browser.close()
            s = requests.Session()
            for c in cookies:
                s.cookies.set(c["name"], c["value"], domain=c.get("domain"), path=c.get("path", "/"))
            return s

        # 2) OAuth 页填写账号密码（可能会出现验证码）
        page.wait_for_selector('input[name="username"]', timeout=60_000)
        page.fill('input[name="username"]', username)
        page.fill('input[name="password"]', password)
        if page.is_visible('input[name="captcha"]'):
            page.screenshot(path="login_need_captcha.png")
            code = input("出现验证码，请查看 login_need_captcha.png 并输入: ").strip()
            page.fill('input[name="captcha"]', code)

        # 3) 点击登录
        page.get_by_role("button", name=re.compile("login", re.I)).click()

        # 4) 等待 URL 回到 ep 域（兼容所有版本）
        try:
            page.wait_for_url(re.compile(r"^https://ep\.bao\.ac\.cn/.*"), timeout=60_000)
        except PWTimeout:
            # 兜底：等一个“已登录”的特征元素
            try:
                page.wait_for_selector('a[href="/ep/user/logout"], img.avatar-xs, text=My Home',
                                       timeout=60_000)
                #page.locator('a[href="/ep/user/logout"], img.avatar-xs, text=My Home').first.wait_for()
            except PWTimeout:
                raise RuntimeError(f"登录可能失败，当前URL: {page.url}")

        # 5) 保存登录态 & 导出 cookies
        ctx.storage_state(path=state_path)
        cookies = ctx.cookies()
        browser.close()

    s = requests.Session()
    for c in cookies:
        s.cookies.set(c["name"], c["value"], domain=c.get("domain"), path=c.get("path", "/"))
    return s


session = login_with_playwright('aujust', 'Liang@981127')

#WXT
params = {
        # "start_datetime": f'{yyyy}-{mm:02d}-{dd-1:02d} 00:00:00',
        # "end_datetime":   f'{yyyy}-{mm:02d}-{dd:02d} 00:00:00',
        "start_datetime": T_START.iso[:-4],
        "end_datetime":   T_NOW.iso[:-4],
        "ra": "",
        "dec": "",
        "radius": ""
    }
headers = {
    "User-Agent": "Mozilla/5.0"
}


response = session.get(API_URL,params=params,headers=headers)
if response.status_code == 200:
    data = response.json()
else:
    print("Request failed with status code:", response.status_code)
response.raise_for_status()
#data = response.json()  # 列表，里面每个元素就是一行记录

df = pd.DataFrame(data)

if not df.empty:
    wxt_source = Table.from_pandas(df)
    notice_type = ['Unverified', 'Unclassified', 'Unclassfied']
    mask = np.isin(wxt_source['category'], notice_type)
    filtered_wxt = wxt_source[mask]
    filtered_wxt.sort('significance',reverse=True)
    filtered_wxt['tags'] = [str(tags) for tags in filtered_wxt['tags']]  #将Tags中的list转换为str
    front_cols = ['obs_id','simbad_name','significance','ra','dec','category','classification','obs_time','url']
    other_cols = [c for c in filtered_wxt.colnames if c not in front_cols]
    filtered_wxt = filtered_wxt[front_cols + other_cols]
    index_col = Column(data=range(1,len(filtered_wxt)+1), name='index')
    filtered_wxt.add_column(index_col, index=0)
    
    fname = f'wxt_candidates_api_{yyyy}{mm:02d}{dd:02d}.csv'
    filtered_wxt.write(os.path.join(DOWNLOAD_ROOT,fname),format='csv',overwrite=True)
    print(f'{fname} is saved to {DOWNLOAD_ROOT}')

else:
    print(response.json())