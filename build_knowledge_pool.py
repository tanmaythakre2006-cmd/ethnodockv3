import os
import urllib.request
import zipfile
import json
import time
import requests
from bs4 import BeautifulSoup
import traceback
import sys

# Ensure stdout handles Chinese characters properly on Windows
sys.stdout.reconfigure(encoding='utf-8')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
GITHUB_DIR = os.path.join(BASE_DIR, "05_github_corpora")
MANIFEST_FILE = os.path.join(BASE_DIR, "manifest.json")
ERROR_LOG = os.path.join(BASE_DIR, "download_errors.log")

def log_error(msg):
    with open(ERROR_LOG, "a", encoding="utf-8") as f:
        f.write(msg + "\n")
    print(f"ERROR: {msg}")

def update_manifest(entry):
    manifest = []
    if os.path.exists(MANIFEST_FILE):
        with open(MANIFEST_FILE, "r", encoding="utf-8") as f:
            try:
                manifest = json.load(f)
            except:
                pass
    manifest.append(entry)
    with open(MANIFEST_FILE, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=4, ensure_ascii=False)

def download_github_repos():
    repos = [
        {"name": "TCM-Sage", "url_main": "https://github.com/AndyZHENG0715/TCM-Sage/archive/refs/heads/main.zip", "url_master": "https://github.com/AndyZHENG0715/TCM-Sage/archive/refs/heads/master.zip"},
        {"name": "TCM-literature-corpus", "url_main": "https://github.com/yunzhangwww/TCM-literature-corpus/archive/refs/heads/main.zip", "url_master": "https://github.com/yunzhangwww/TCM-literature-corpus/archive/refs/heads/master.zip"},
        {"name": "ZhongJing-OMNI", "url_main": "https://github.com/pariskang/ZhongJing-OMNI/archive/refs/heads/main.zip", "url_master": "https://github.com/pariskang/ZhongJing-OMNI/archive/refs/heads/master.zip"}
    ]
    
    for repo in repos:
        zip_path = os.path.join(GITHUB_DIR, f"{repo['name']}.zip")
        extract_path = os.path.join(GITHUB_DIR, repo['name'])
        if os.path.exists(extract_path):
            print(f"Skipping {repo['name']}, already exists.")
            continue
            
        print(f"Downloading {repo['name']}...")
        try:
            urllib.request.urlretrieve(repo['url_main'], zip_path)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                try:
                    urllib.request.urlretrieve(repo['url_master'], zip_path)
                except Exception as e2:
                    log_error(f"Failed to download {repo['name']} from master: {e2}")
                    continue
            else:
                log_error(f"Failed to download {repo['name']} from main: {e}")
                continue
        except Exception as e:
            log_error(f"Failed to download {repo['name']}: {e}")
            continue
            
        print(f"Extracting {repo['name']}...")
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_path)
            os.remove(zip_path)
            
            update_manifest({
                "filename": f"{repo['name']}/",
                "book_title": repo['name'],
                "author": "Various",
                "dynasty": "Various",
                "era_folder": "05_github_corpora",
                "source_url": repo['url_main'].split('/archive')[0],
                "language": "Classical Chinese"
            })
            
        except Exception as e:
            log_error(f"Failed to extract {repo['name']}: {e}")

def scrape_wikisource():
    categories = ["Category:中醫學", "Category:本草", "Category:黃帝內經", "Category:傷寒論"]
    api_url = "https://zh.wikisource.org/w/api.php"
    
    headers = {"User-Agent": "EthnoDockBot/1.0 (bot@ethnodock.com)"}
    
    for category in categories:
        print(f"Scraping category: {category}")
        params = {
            "action": "query",
            "list": "categorymembers",
            "cmtitle": category,
            "cmlimit": "max",
            "format": "json"
        }
        
        try:
            r = requests.get(api_url, params=params, headers=headers)
            r.raise_for_status()
            data = r.json()
            members = data.get("query", {}).get("categorymembers", [])
            
            for member in members:
                title = member["title"]
                if member["ns"] != 0: # Only main namespace (articles)
                    continue
                    
                print(f"Fetching Wikisource: {title}")
                page_params = {
                    "action": "query",
                    "prop": "extracts",
                    "explaintext": True,
                    "titles": title,
                    "format": "json"
                }
                page_r = requests.get(api_url, params=page_params, headers=headers)
                page_r.raise_for_status()
                page_data = page_r.json()
                pages = page_data.get("query", {}).get("pages", {})
                for page_id, page_info in pages.items():
                    content = page_info.get("extract", "")
                    if content:
                        folder = "04_ming_qing_dynasties"
                        if "漢" in content[:100] or "汉" in content[:100]: folder = "01_han_dynasty"
                        elif "唐" in content[:100] or "隋" in content[:100]: folder = "02_sui_tang_dynasties"
                        elif "宋" in content[:100] or "元" in content[:100]: folder = "03_song_yuan_dynasties"
                        
                        safe_title = "".join([c for c in title if c.isalpha() or c.isdigit() or c==' ']).rstrip()
                        if not safe_title:
                            safe_title = str(page_id)
                            
                        filename = f"{safe_title}.txt"
                        filepath = os.path.join(BASE_DIR, folder, filename)
                        
                        with open(filepath, "w", encoding="utf-8") as f:
                            f.write(content)
                            
                        update_manifest({
                            "filename": filename,
                            "book_title": title,
                            "author": "Unknown",
                            "dynasty": "Various",
                            "era_folder": folder,
                            "source_url": f"https://zh.wikisource.org/wiki/{title}",
                            "language": "Classical Chinese"
                        })
                time.sleep(0.5) # Gentle delay
                
        except Exception as e:
            log_error(f"Wikisource scraping failed for {category}: {traceback.format_exc()}")

if __name__ == '__main__':
    print("Starting Knowledge Pool Build...")
    if not os.path.exists(MANIFEST_FILE):
        with open(MANIFEST_FILE, "w", encoding="utf-8") as f:
            f.write("[]")
            
    download_github_repos()
    print("GitHub corpora downloaded.")
    
    scrape_wikisource()
    print("Wikisource scraped.")
    
    print("Build complete.")
