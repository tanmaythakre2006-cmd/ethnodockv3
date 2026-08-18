import os
import json
import argparse
import sys

sys.stdout.reconfigure(encoding='utf-8')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MANIFEST_FILE = os.path.join(BASE_DIR, "manifest.json")

def verify_pool():
    if not os.path.exists(MANIFEST_FILE):
        print(f"Manifest file not found: {MANIFEST_FILE}")
        return
        
    with open(MANIFEST_FILE, "r", encoding="utf-8") as f:
        try:
            manifest = json.load(f)
        except json.JSONDecodeError:
            print("Manifest is not a valid JSON.")
            return
            
    print(f"Loaded {len(manifest)} entries from manifest.")
    
    missing_files = []
    for entry in manifest:
        folder = entry.get("era_folder", "")
        filename = entry.get("filename", "")
        filepath = os.path.join(BASE_DIR, folder, filename)
        
        if filename.endswith('/'):
            if not os.path.exists(filepath) and not os.path.exists(filepath[:-1]):
                missing_files.append(filepath)
        else:
            if not os.path.exists(filepath):
                missing_files.append(filepath)
                
    if missing_files:
        print(f"Verification Failed. Missing {len(missing_files)} files:")
        for mf in missing_files[:10]:
            print(f"  - {mf}")
        if len(missing_files) > 10:
            print(f"  ... and {len(missing_files) - 10} more.")
    else:
        print("Verification Passed. All files in manifest exist.")

def search_pool(keyword):
    print(f"Searching for '{keyword}'...")
    folders = [
        "01_han_dynasty",
        "02_sui_tang_dynasties",
        "03_song_yuan_dynasties",
        "04_ming_qing_dynasties",
        "05_github_corpora",
        "06_english_translations"
    ]
    
    match_count = 0
    for folder in folders:
        folder_path = os.path.join(BASE_DIR, folder)
        if not os.path.exists(folder_path):
            continue
            
        for root, _, files in os.walk(folder_path):
            for file in files:
                if file.endswith('.txt') or file.endswith('.md'):
                    filepath = os.path.join(root, file)
                    try:
                        with open(filepath, "r", encoding="utf-8") as f:
                            content = f.read()
                            if keyword in content:
                                match_count += 1
                                print(f"Match found in: {os.path.relpath(filepath, BASE_DIR)}")
                    except Exception:
                        pass
                        
    print(f"Total files containing '{keyword}': {match_count}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Search and verify EthnoDock Knowledge Pool")
    parser.add_argument("--verify", action="store_true", help="Verify the integrity of the data lake against the manifest")
    parser.add_argument("--search", type=str, help="Search for a keyword across all raw text files")
    
    args = parser.parse_args()
    
    if args.verify:
        verify_pool()
    elif args.search:
        search_pool(args.search)
    else:
        parser.print_help()
