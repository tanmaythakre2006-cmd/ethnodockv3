import os
import json
import time
import requests
import re
import sys

sys.stdout.reconfigure(encoding='utf-8')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MANIFEST_FILE = os.path.join(BASE_DIR, "manifest.json")
ERROR_LOG = os.path.join(BASE_DIR, "expansion_errors.log")

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

def strip_tei_xml(xml_content):
    text = re.sub(r'<[^>]+>', '', xml_content)
    text = re.sub(r'\n\s*\n', '\n\n', text)
    return text.strip()

def identify_era(title, repo_id):
    # Rudimentary heuristic to map Kanripo medicine texts to era folders
    title = title.lower()
    if repo_id in ["KR5a0001", "KR5a0002", "KR5a0010"]: # Neijing, Lingshu, Shanghan
        return "01_han_dynasty", "Han Dynasty"
    elif repo_id in ["KR5a0018", "KR5a0020"]: # Qianjin Fang, Waitai Miyao
        return "02_sui_tang_dynasties", "Sui / Tang Dynasty"
    elif repo_id in ["KR5a0033", "KR5a0035"]: # Taiping Huimin, etc
        return "03_song_yuan_dynasties", "Song / Yuan Dynasty"
    else:
        return "04_ming_qing_dynasties", "Ming / Qing / Unknown"

def fetch_kanripo():
    print("Starting Kanripo Mass Ingestion...")
    # List of known KR5 medical repositories to attempt fetching
    repos = [
        {"id": "KR5a0001", "title": "Huangdi Neijing Suwen"},
        {"id": "KR5a0002", "title": "Lingshu Jing"},
        {"id": "KR5a0010", "title": "Shanghan Lun"},
        {"id": "KR5a0018", "title": "Beiji Qianjin Yaofang"},
        {"id": "KR5a0033", "title": "Taiping Huimin Heji Jufang"}
    ]
    
    headers = {"User-Agent": "EthnoDockBot/1.0"}
    
    for repo in repos:
        repo_id = repo['id']
        title = repo['title']
        print(f"Querying GitHub API for Kanripo repo: {repo_id}...")
        
        api_url = f"https://api.github.com/repos/kanripo/{repo_id}/contents"
        try:
            r = requests.get(api_url, headers=headers)
            if r.status_code == 200:
                files = r.json()
                # Find .txt files first (Kanripo usually provides plain text exports)
                txt_files = [f for f in files if f['name'].endswith('.txt')]
                if not txt_files:
                    # Fallback to .xml files if no txt
                    txt_files = [f for f in files if f['name'].endswith('.xml')]
                
                for file_info in txt_files:
                    file_name = file_info['name']
                    download_url = file_info['download_url']
                    
                    print(f"Downloading {file_name} from {repo_id}...")
                    dl_r = requests.get(download_url, headers=headers)
                    if dl_r.status_code == 200:
                        content = dl_r.text
                        if file_name.endswith('.xml'):
                            content = strip_tei_xml(content)
                            
                        era_folder, dynasty = identify_era(title, repo_id)
                        save_name = f"{repo_id}_{file_name.replace('.xml', '.txt')}"
                        filepath = os.path.join(BASE_DIR, era_folder, save_name)
                        
                        with open(filepath, "w", encoding="utf-8") as out_f:
                            out_f.write(content)
                            
                        update_manifest({
                            "filename": save_name,
                            "book_title": title,
                            "author": "Unknown / Various",
                            "dynasty": dynasty,
                            "era_folder": era_folder,
                            "source_url": download_url,
                            "language": "Classical Chinese",
                            "provenance_chain": "Kanripo (GitHub API)",
                            "authenticity_confidence": "Highest (TEI-XML Academic)"
                        })
                        print(f"Successfully saved {save_name} to {era_folder}")
                        break # Just download the primary file to save time in this POC
                    else:
                        log_error(f"Failed to download file {file_name} for {repo_id}")
            else:
                log_error(f"GitHub API returned {r.status_code} for {repo_id}")
        except Exception as e:
            log_error(f"Error processing {repo_id}: {e}")
            
        time.sleep(1) # Respect GitHub API rate limits

if __name__ == '__main__':
    fetch_kanripo()
    print("Phase 4 Kanripo Expansion Complete.")
