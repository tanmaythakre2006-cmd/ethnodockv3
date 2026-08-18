import os
import json
import time
import requests
import re
import sys
import zipfile
import io

sys.stdout.reconfigure(encoding='utf-8')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MANIFEST_FILE = os.path.join(BASE_DIR, "manifest.json")
ERROR_LOG = os.path.join(BASE_DIR, "expansion_errors.log")

def log_error(msg):
    with open(ERROR_LOG, "a", encoding="utf-8") as f:
        f.write(msg + "\n")
    print(f"ERROR: {msg}")

def update_manifest(entries):
    manifest = []
    if os.path.exists(MANIFEST_FILE):
        with open(MANIFEST_FILE, "r", encoding="utf-8") as f:
            try:
                manifest = json.load(f)
            except:
                pass
    manifest.extend(entries)
    with open(MANIFEST_FILE, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=4, ensure_ascii=False)

def strip_tei_xml(xml_content):
    text = re.sub(r'<[^>]+>', '', xml_content)
    text = re.sub(r'\n\s*\n', '\n\n', text)
    return text.strip()

def identify_era(title, repo_id):
    title = title.lower()
    if repo_id in ["KR5a0001", "KR5a0002", "KR5a0010"]:
        return "01_han_dynasty", "Han Dynasty"
    elif repo_id in ["KR5a0017", "KR5a0018", "KR5a0020"]:
        return "02_sui_tang_dynasties", "Sui / Tang Dynasty"
    elif repo_id in ["KR5a0033", "KR5a0035", "KR5a0038"]:
        return "03_song_yuan_dynasties", "Song / Yuan Dynasty"
    else:
        return "04_ming_qing_dynasties", "Ming / Qing / Unknown"

def fetch_kanripo_zips():
    print("Starting Massive Data Lake Expansion (Direct ZIP Download)...")
    
    # Target large compendiums to massively expand text volume
    repos = [
        {"id": "KR5a0016", "title": "Bencao Gangmu (Compendium of Materia Medica)"},
        {"id": "KR5a0017", "title": "Bencao Gangmu Shiyi"},
        {"id": "KR5a0035", "title": "Shengji Zonglu"},
        {"id": "KR5a0038", "title": "Bencao Pinhui Jingyao"},
        {"id": "KR5a0068", "title": "Yixue Zhengchuan"}
    ]
    
    headers = {"User-Agent": "EthnoDockBot/1.0"}
    all_manifest_entries = []
    
    for repo in repos:
        repo_id = repo['id']
        title = repo['title']
        zip_url = f"https://github.com/kanripo/{repo_id}/archive/refs/heads/master.zip"
        
        print(f"Downloading massive archive: {zip_url}...")
        try:
            r = requests.get(zip_url, headers=headers, stream=True)
            if r.status_code == 200:
                print(f"Successfully downloaded archive for {repo_id}. Unpacking...")
                z = zipfile.ZipFile(io.BytesIO(r.content))
                
                extracted_count = 0
                for filename in z.namelist():
                    if filename.endswith(".xml") or filename.endswith(".txt"):
                        if "README" in filename.upper():
                            continue
                            
                        content = z.read(filename).decode('utf-8', errors='ignore')
                        if filename.endswith('.xml'):
                            content = strip_tei_xml(content)
                            
                        # If file is too small, skip it (often just meta-tags)
                        if len(content) < 100:
                            continue
                            
                        era_folder, dynasty = identify_era(title, repo_id)
                        safe_filename = os.path.basename(filename).replace('.xml', '.txt')
                        save_name = f"{repo_id}_{safe_filename}"
                        filepath = os.path.join(BASE_DIR, era_folder, save_name)
                        
                        # Ensure dir exists
                        os.makedirs(os.path.join(BASE_DIR, era_folder), exist_ok=True)
                        
                        with open(filepath, "w", encoding="utf-8") as out_f:
                            out_f.write(content)
                            
                        all_manifest_entries.append({
                            "filename": save_name,
                            "book_title": title,
                            "author": "Unknown / Various",
                            "dynasty": dynasty,
                            "era_folder": era_folder,
                            "source_url": zip_url,
                            "language": "Classical Chinese",
                            "provenance_chain": "Kanripo (Direct ZIP)",
                            "authenticity_confidence": "Highest (TEI-XML Academic)"
                        })
                        extracted_count += 1
                        
                print(f"Extracted and sorted {extracted_count} volumes from {repo_id}!")
            else:
                log_error(f"Failed to fetch ZIP for {repo_id}. HTTP {r.status_code}")
        except Exception as e:
            log_error(f"Error processing {repo_id}: {e}")
            
        time.sleep(2) # Prevent triggering GitHub anti-DDOS
        
    print(f"\nWriting {len(all_manifest_entries)} new volumes to manifest...")
    update_manifest(all_manifest_entries)
    print("Massive Expansion Complete!")

if __name__ == '__main__':
    fetch_kanripo_zips()
