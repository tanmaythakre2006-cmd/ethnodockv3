import os
import sys
import re
import json
import yaml
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OKF_INDEX_FILE = os.path.join(BASE_DIR, "okf_index.json")

def load_known_okf_entities():
    if not os.path.exists(OKF_INDEX_FILE):
        return set()
    with open(OKF_INDEX_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)
        return set(data.keys())

def discover_okf_candidates():
    print("[OKF Discovery Scanner] Initiating deterministic scan across classical corpora...")
    known_entities = load_known_okf_entities()
    
    candidates_counter = Counter()
    pattern = re.compile(r'([\u4e00-\u9fa5]{1,2}(?:草|根|叶|皮|子|仁|花|枝))')

    folders_to_scan = [
        "01_han_dynasty", "02_sui_tang_dynasties", "03_song_yuan_dynasties",
        "04_ming_qing_dynasties", "05_github_corpora", "06_english_translations"
    ]

    total_scanned_files = 0
    for folder in folders_to_scan:
        folder_path = os.path.join(BASE_DIR, folder)
        if not os.path.exists(folder_path):
            continue

        for root, _, files in os.walk(folder_path):
            for file in files:
                if file.endswith('.txt') or file.endswith('.md'):
                    total_scanned_files += 1
                    filepath = os.path.join(root, file)
                    try:
                        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                            content = f.read()
                            matches = pattern.findall(content)
                            for match in matches:
                                if match not in known_entities:
                                    candidates_counter[match] += 1
                    except Exception:
                        pass

    top_candidates = candidates_counter.most_common(20)
    print(f"[OKF Scanner] Complete. Scanned {total_scanned_files} files.")
    print(f"Found {len(candidates_counter)} unknown candidate terms.\n")
    print("--- TOP DISCOVERED ENTITIES FORMATTED AS OKF TEMPLATES ---")
    
    okf_templates = []
    for word, count in top_candidates:
        slug = f"herb_candidate_{word}"
        okf_template = {
            "okf_version": "1.0",
            "entity_id": slug,
            "type": "ethnobotanical_species",
            "title": f"Candidate: {word}",
            "taxonomy": {
                "chinese_name": word,
                "binomial_name": "Pending Botanical Verification"
            },
            "discovery_frequency": count,
            "status": "candidate_for_ingestion"
        }
        okf_templates.append(okf_template)

    print(yaml.dump(okf_templates[:3], sort_keys=False, allow_unicode=True))
    print("----------------------------------------------------------")

if __name__ == '__main__':
    discover_okf_candidates()
