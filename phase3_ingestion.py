import os
import json
import csv
import urllib.request
import sys

sys.stdout.reconfigure(encoding='utf-8')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MANIFEST_FILE = os.path.join(BASE_DIR, "manifest.json")
ERROR_LOG = os.path.join(BASE_DIR, "expansion_errors.log")

PHASE3_MD_DIR = os.path.join(BASE_DIR, "08_systems_pharmacology")
RAW_DATA_DIR = os.path.join(BASE_DIR, "10_raw_structured_data")

os.makedirs(PHASE3_MD_DIR, exist_ok=True)
os.makedirs(RAW_DATA_DIR, exist_ok=True)

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

def fetch_and_convert_csv():
    # Using herb_feature.csv from a public ML TCM repo as our structured data
    csv_url = "https://raw.githubusercontent.com/herb-medicne/Meridian-prediction/master/herb_feature.csv"
    raw_csv_path = os.path.join(RAW_DATA_DIR, "herb_feature.csv")
    
    print(f"Downloading raw structured data from {csv_url}...")
    try:
        urllib.request.urlretrieve(csv_url, raw_csv_path)
        print("Raw CSV saved successfully.")
        
        # Log the raw file in manifest
        update_manifest({
            "filename": "herb_feature.csv",
            "book_title": "TCM Herb ML Features",
            "author": "herb-medicne/Meridian-prediction",
            "dynasty": "Modern",
            "era_folder": "10_raw_structured_data",
            "source_url": csv_url,
            "language": "English",
            "provenance_chain": "GitHub (Raw CSV Extract)",
            "scientific_validation_level": "High (Computational Models)",
            "authenticity_confidence": "High"
        })
        
        print("Converting CSV to Markdown with YAML Frontmatter...")
        with open(raw_csv_path, "r", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f)
            count = 0
            for row in reader:
                # We'll just process the first 100 rows for the POC to avoid flooding the directory
                if count >= 100:
                    break
                    
                herb_id = row.get("herb_id", f"UNKNOWN_{count}")
                herb_name = row.get("herb_name", "Unknown Herb").strip()
                
                # Build MD content
                md_content = f"""---
source_db: "herb-medicne/Meridian-prediction"
herb_id: "{herb_id}"
herb_name: "{herb_name}"
scientific_validation_level: "High (Computational Models)"
---

# {herb_name} (Systems Pharmacology Features)

## Extracted Raw Features
"""
                for key, value in row.items():
                    if key not in ["herb_id", "herb_name"] and value:
                        md_content += f"* **{key}**: {value}\n"
                
                safe_name = "".join([c for c in herb_name if c.isalnum() or c.isspace()]).strip().replace(" ", "_")
                filename = f"{herb_id}_{safe_name}.md"
                filepath = os.path.join(PHASE3_MD_DIR, filename)
                
                with open(filepath, "w", encoding="utf-8") as out_f:
                    out_f.write(md_content)
                
                update_manifest({
                    "filename": filename,
                    "book_title": f"{herb_name} Pharmacology Report",
                    "author": "Automated Conversion",
                    "dynasty": "Modern",
                    "era_folder": "08_systems_pharmacology",
                    "source_url": "herb_feature.csv",
                    "language": "English",
                    "provenance_chain": "CSV -> MD Conversion",
                    "authenticity_confidence": "High"
                })
                
                count += 1
        print(f"Successfully converted {count} herb records into Markdown.")
        
    except Exception as e:
        log_error(f"Failed to fetch/convert CSV: {e}")

if __name__ == '__main__':
    print("Starting Phase 3 Ingestion...")
    fetch_and_convert_csv()
    print("Phase 3 Complete.")
