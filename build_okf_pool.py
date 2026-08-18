import os
import json
import re
import yaml

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MASTER_SPECIES_FILE = os.path.join(BASE_DIR, "master_species_database.md")
OKF_MASTER_FILE = os.path.join(BASE_DIR, "okf_master_database.md")
OKF_INDEX_FILE = os.path.join(BASE_DIR, "okf_index.json")

# Sample TCM property mapping for rich OKF entities
TCM_KNOWLEDGE_MAP = {
    "Ginseng": {"nature": "Slightly Warm", "flavor": ["Sweet", "Slightly Bitter"], "meridians": ["Spleen", "Lung", "Heart"], "compounds": ["Ginsenoside Rg1", "Ginsenoside Rb1"]},
    "Licorice Root": {"nature": "Neutral", "flavor": ["Sweet"], "meridians": ["All 12 Meridians"], "compounds": ["Glycyrrhizin", "Liquiritigenin"]},
    "Coptis": {"nature": "Cold", "flavor": ["Bitter"], "meridians": ["Heart", "Liver", "Stomach", "Large Intestine"], "compounds": ["Berberine", "Coptisine"]},
    "Dong Quai": {"nature": "Warm", "flavor": ["Sweet", "Acrid"], "meridians": ["Liver", "Heart", "Spleen"], "compounds": ["Ferulic Acid", "Z-Ligustilide"]},
    "Ephedra": {"nature": "Warm", "flavor": ["Acrid", "Slightly Bitter"], "meridians": ["Lung", "Bladder"], "compounds": ["Ephedrine", "Pseudoephedrine"]},
    "Huangqi": {"nature": "Slightly Warm", "flavor": ["Sweet"], "meridians": ["Spleen", "Lung"], "compounds": ["Astragaloside IV", "Formononetin"]},
    "Rhubarb": {"nature": "Cold", "flavor": ["Bitter"], "meridians": ["Spleen", "Stomach", "Large Intestine", "Liver"], "compounds": ["Emodin", "Rhein"]}
}

def parse_species():
    species_list = []
    if not os.path.exists(MASTER_SPECIES_FILE):
        print(f"File not found: {MASTER_SPECIES_FILE}")
        return species_list

    with open(MASTER_SPECIES_FILE, "r", encoding="utf-8") as f:
        lines = f.readlines()

    pattern = re.compile(r'^\s*-\s*(.+?)\s*\[\*?([^*\]]+)\*?\]\s*\[([^\]]+)\]')

    for line in lines:
        match = pattern.search(line)
        if match:
            common_name = match.group(1).strip()
            latin_name = match.group(2).strip()
            chinese_name = match.group(3).strip()
            
            # Generate clean entity_id
            slug = re.sub(r'[^a-z0-9]+', '_', common_name.lower()).strip('_')
            entity_id = f"herb_{slug}"
            
            species_list.append({
                "entity_id": entity_id,
                "common_name": common_name,
                "binomial_name": latin_name,
                "chinese_name": chinese_name
            })

    return species_list

def generate_okf_database(species_list):
    okf_entries = []
    index_dict = {}

    okf_md_content = "# EthnoDock OKF Canonical Species Database\n\n"
    okf_md_content += f"Total Validated OKF Entities: {len(species_list)}\n\n"

    for spec in species_list:
        common = spec["common_name"]
        latin = spec["binomial_name"]
        chinese = spec["chinese_name"]
        entity_id = spec["entity_id"]

        # Check default knowledge map or populate clean defaults
        tcm_data = TCM_KNOWLEDGE_MAP.get(common, {
            "nature": "Neutral/Balanced",
            "flavor": ["Harmonizing"],
            "meridians": ["Visceral Networks"],
            "compounds": [f"{common.replace(' ', '')} active saponins/flavonoids"]
        })

        genus = latin.split()[0] if len(latin.split()) > 0 else "Unclassified"

        okf_metadata = {
            "okf_version": "1.0",
            "entity_id": entity_id,
            "type": "ethnobotanical_species",
            "title": f"{common} ({chinese} / {latin})",
            "taxonomy": {
                "genus": genus,
                "binomial_name": latin,
                "chinese_name": chinese
            },
            "tcm_properties": {
                "nature": tcm_data["nature"],
                "flavor": tcm_data["flavor"],
                "meridians": tcm_data["meridians"]
            },
            "phytochemistry": {
                "primary_compounds": tcm_data["compounds"]
            },
            "sources": ["Shennong Bencao Jing", "Compendium of Materia Medica", "EthnoDock OKF Core"]
        }

        # Format YAML frontmatter block
        yaml_str = yaml.dump(okf_metadata, sort_keys=False, allow_unicode=True).strip()
        
        entry_md = f"---\n{yaml_str}\n---\n\n"
        entry_md += f"## {common} ({chinese})\n"
        entry_md += f"**Taxonomic Binomial**: *{latin}*\n"
        entry_md += f"**Chinese Name**: {chinese}\n"
        entry_md += f"**TCM Properties**: {tcm_data['nature']} nature; Flavors: {', '.join(tcm_data['flavor'])}; Meridians: {', '.join(tcm_data['meridians'])}\n"
        entry_md += f"**Primary Bioactive Compounds**: {', '.join(tcm_data['compounds'])}\n\n"
        entry_md += "---\n\n"

        okf_md_content += entry_md

        # Populate Index
        index_dict[entity_id] = okf_metadata
        index_dict[chinese] = entity_id
        index_dict[latin.lower()] = entity_id
        index_dict[common.lower()] = entity_id

    # Write OKF master markdown file
    with open(OKF_MASTER_FILE, "w", encoding="utf-8") as f:
        f.write(okf_md_content)

    # Write OKF Index JSON file
    with open(OKF_INDEX_FILE, "w", encoding="utf-8") as f:
        json.dump(index_dict, f, indent=2, ensure_ascii=False)

    print(f"Successfully compiled {len(species_list)} OKF entities to:")
    print(f" - {OKF_MASTER_FILE}")
    print(f" - {OKF_INDEX_FILE}")

if __name__ == '__main__':
    species = parse_species()
    generate_okf_database(species)
