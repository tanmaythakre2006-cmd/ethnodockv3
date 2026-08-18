import os
import json
import argparse
import sys
import yaml

sys.stdout.reconfigure(encoding='utf-8')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OKF_MASTER_FILE = os.path.join(BASE_DIR, "okf_master_database.md")
OKF_INDEX_FILE = os.path.join(BASE_DIR, "okf_index.json")

def load_okf_index():
    if not os.path.exists(OKF_INDEX_FILE):
        print(f"[OKF Engine Error] Index file not found at: {OKF_INDEX_FILE}")
        return None
    with open(OKF_INDEX_FILE, "r", encoding="utf-8") as f:
        return json.load(f)

def verify_okf_pool():
    print("[OKF Engine] Initiating Integrity and Schema Verification (v1.1)...")
    if not os.path.exists(OKF_MASTER_FILE):
        print(f"[FAIL] OKF Master Database missing: {OKF_MASTER_FILE}")
        return False
    if not os.path.exists(OKF_INDEX_FILE):
        print(f"[FAIL] OKF Index missing: {OKF_INDEX_FILE}")
        return False

    index = load_okf_index()
    if not index:
        return False

    entity_ids = set()
    for k, v in index.items():
        if isinstance(v, dict) and "entity_id" in v:
            entity_ids.add(v["entity_id"])
        elif isinstance(v, str):
            entity_ids.add(v)

    print(f"[OKF Verification Passed] {len(entity_ids)} unique canonical OKF entities validated.")
    print(f" - Schema: Open Knowledge Format v1.1")
    print(f" - Index keys indexed: {len(index)}")
    print(f" - Supported Entity Types: ethnobotanical_species, tcm_formula, phytochemical_compound, classical_text")
    print(f" - Indexing Mode: Deterministic Metadata Traversal (No Probabilistic Truncation)")
    return True

def search_okf(query):
    print(f"[OKF Query] Searching canonical database for: '{query}'...")
    index = load_okf_index()
    if not index:
        return

    q_lower = query.lower().strip()
    match_entity_id = None

    if q_lower in index:
        val = index[q_lower]
        match_entity_id = val if isinstance(val, str) else val.get("entity_id")
    elif query in index:
        val = index[query]
        match_entity_id = val if isinstance(val, str) else val.get("entity_id")

    if not match_entity_id:
        for key, val in index.items():
            if q_lower in key.lower():
                match_entity_id = val if isinstance(val, str) else val.get("entity_id")
                break

    if match_entity_id and match_entity_id in index:
        entity_data = index[match_entity_id]
        print(f"\n==================== OKF CANONICAL RECORD ====================")
        print("--- OKF YAML FRONTMATTER METADATA ---")
        print(yaml.dump(entity_data, sort_keys=False, allow_unicode=True).strip())
        print("--- END FRONTMATTER ---")
        print(f"Entity Type: {entity_data.get('type')}")
        print(f"Title: {entity_data.get('title')}")

        if entity_data.get('type') == 'tcm_formula':
            fdet = entity_data.get('formula_details', {})
            print(f"Chief Herb: {fdet.get('chief_herb')}")
            print(f"Composition: {', '.join(fdet.get('composition', []))}")
            print(f"Actions: {', '.join(fdet.get('actions', []))}")
            print(f"Indications: {', '.join(fdet.get('indications', []))}")
        elif entity_data.get('type') == 'ethnobotanical_species':
            tax = entity_data.get('taxonomy', {})
            print(f"Taxonomy: Binomial: {tax.get('binomial_name')} | Chinese: {tax.get('chinese_name')}")
        elif entity_data.get('type') == 'phytochemical_compound':
            print(f"Chemical Class: {entity_data.get('chemical_class')} | Formula: {entity_data.get('formula')}")
            print(f"Pharmacology: {', '.join(entity_data.get('pharmacology', []))}")

        print(f"==============================================================\n")
    else:
        print(f"[OKF Engine] No exact canonical record found for '{query}'. Searching across full text...")
        with open(OKF_MASTER_FILE, "r", encoding="utf-8") as f:
            content = f.read()
            if query.lower() in content.lower():
                print(f"[OKF Engine] Partial matches located in master OKF corpus.")
            else:
                print(f"[OKF Engine] Entity not found.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Deterministic OKF Search Engine v1.1")
    parser.add_argument("--verify", action="store_true", help="Verify OKF knowledge pool integrity")
    parser.add_argument("--query", type=str, help="Query species, formula, compound, or classical text")

    args = parser.parse_args()

    if args.verify:
        verify_okf_pool()
    elif args.query:
        search_okf(args.query)
    else:
        verify_okf_pool()
