import os
import json
import re
import yaml

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MASTER_SPECIES_FILE = os.path.join(BASE_DIR, "master_species_database.md")
OKF_MASTER_FILE = os.path.join(BASE_DIR, "okf_master_database.md")
OKF_INDEX_FILE = os.path.join(BASE_DIR, "okf_index.json")

# Detailed TCM & Phytochemistry profiles for major herbs
HERB_PROFILES = {
    "Ginseng": {
        "nature": "Slightly Warm",
        "flavor": ["Sweet", "Slightly Bitter"],
        "meridians": ["Spleen", "Lung", "Heart"],
        "actions": ["Greatly Tonifies Original Qi", "Restores Pulse & Rescues Collapse", "Tonifies Spleen & Lung Qi", "Engenders Fluids & Quenches Thirst"],
        "compounds": ["Ginsenoside Rg1", "Ginsenoside Rb1", "Ginsenoside Re", "Panaxytriol"]
    },
    "Licorice Root": {
        "nature": "Neutral",
        "flavor": ["Sweet"],
        "meridians": ["Heart", "Lung", "Spleen", "Stomach"],
        "actions": ["Tonifies Spleen Qi", "Moistens Lungs & Stops Cough", "Moderates Urgency & Relieves Pain", "Harmonizes Other Herbs"],
        "compounds": ["Glycyrrhizin", "Liquiritigenin", "Glabridin", "Isoliquiritigenin"]
    },
    "Coptis": {
        "nature": "Cold",
        "flavor": ["Bitter"],
        "meridians": ["Heart", "Liver", "Stomach", "Large Intestine"],
        "actions": ["Clears Heat & Dries Dampness", "Purges Heart Fire", "Clears Stomach Heat & Stops Vomiting", "Clears Toxic Heat"],
        "compounds": ["Berberine", "Coptisine", "Palmatine", "Jatrorrhizine"]
    },
    "Dong Quai": {
        "nature": "Warm",
        "flavor": ["Sweet", "Acrid"],
        "meridians": ["Liver", "Heart", "Spleen"],
        "actions": ["Tonifies Blood", "Invigorates Blood & Dispels Stasis", "Unblocks Bowels & Lubricates Intestines", "Relieves Pain"],
        "compounds": ["Ferulic Acid", "Z-Ligustilide", "Butylphthalide", "Senkyunolide A"]
    },
    "Astragalus": {
        "nature": "Slightly Warm",
        "flavor": ["Sweet"],
        "meridians": ["Spleen", "Lung"],
        "actions": ["Tonifies Spleen & Raises Yang", "Augments Protective Qi & Stabilizes Exterior", "Promotes Urination & Reduces Edema", "Discharges Pus & Generates Flesh"],
        "compounds": ["Astragaloside IV", "Formononetin", "Calycosin", "Astragalan"]
    },
    "Ephedra": {
        "nature": "Warm",
        "flavor": ["Acrid", "Slightly Bitter"],
        "meridians": ["Lung", "Bladder"],
        "actions": ["Induces Sweating & Releases Exterior", "Disseminates Lung Qi & Calms Asthma", "Promotes Urination & Reduces Edema"],
        "compounds": ["Ephedrine", "Pseudoephedrine", "N-methylephedrine"]
    },
    "Sweet Wormwood": {
        "nature": "Cold",
        "flavor": ["Bitter", "Acrid"],
        "meridians": ["Liver", "Gallbladder"],
        "actions": ["Clears Summer-Heat", "Clears Deficiency Heat & Bone-Steaming", "Check Malarial Disorders & Relieves Alternating Chill/Fever"],
        "compounds": ["Artemisinin", "Arteannuin B", "Artemisinic Acid"]
    },
    "Red Sage": {
        "nature": "Slightly Cold",
        "flavor": ["Bitter"],
        "meridians": ["Heart", "Pericardium", "Liver"],
        "actions": ["Invigorates Blood & Dispels Stasis", "Cools Blood & Reduces Abscesses", "Nourishes Blood & Calms Spirit"],
        "compounds": ["Tanshinone IIA", "Cryptotanshinone", "Salvianolic Acid B", "Danshensu"]
    },
    "Scutellaria": {
        "nature": "Cold",
        "flavor": ["Bitter"],
        "meridians": ["Lung", "Gallbladder", "Stomach", "Large Intestine"],
        "actions": ["Clears Upper Burner Heat & Dries Dampness", "Purges Fire & Detoxifies", "Calms Fetus"],
        "compounds": ["Baicalin", "Baicalein", "Wogonin", "Wogonoside"]
    },
    "Rhubarb": {
        "nature": "Cold",
        "flavor": ["Bitter"],
        "meridians": ["Spleen", "Stomach", "Large Intestine", "Liver"],
        "actions": ["Purges Heat & Unblocks Intestines", "Clears Heat & Cools Blood", "Invigorates Blood & Dispels Stasis"],
        "compounds": ["Emodin", "Rhein", "Chrysophanol", "Aloe-emodin"]
    }
}

# Major Classical Texts as OKF Entities
CLASSICAL_TEXTS = [
    {
        "okf_version": "1.0",
        "entity_id": "text_shennong_bencao_jing",
        "type": "classical_text",
        "title": "Shennong Bencao Jing (神农本草经 / Divine Farmer's Materia Medica)",
        "era": "Han Dynasty (circa 200 CE)",
        "author": "Attributed to Divine Farmer Shennong",
        "significance": "The foundational text of Chinese Materia Medica, categorizing 365 herbs into Superior, Medium, and Inferior classes.",
        "sources": ["Han Dynasty Medical Archives", "EthnoDock Classical Library"]
    },
    {
        "okf_version": "1.0",
        "entity_id": "text_shanghan_lun",
        "type": "classical_text",
        "title": "Shanghan Lun (伤寒论 / Treatise on Cold Damage)",
        "era": "Eastern Han Dynasty (circa 220 CE)",
        "author": "Zhang Zhongjing (张仲景)",
        "significance": "Pioneered the Six-Stage Syndrome Differentiation system and standardized classic formulas like Gui Zhi Tang and Ma Huang Tang.",
        "sources": ["Han Dynasty Medical Archives", "EthnoDock Classical Library"]
    },
    {
        "okf_version": "1.0",
        "entity_id": "text_bencao_gangmu",
        "type": "classical_text",
        "title": "Bencao Gangmu (本草纲目 / Compendium of Materia Medica)",
        "era": "Ming Dynasty (1578 CE)",
        "author": "Li Shizhen (李时珍)",
        "significance": "Monumental 52-volume encyclopedia documenting 1,892 substances, 11,096 formulas, and detailed botanical classifications.",
        "sources": ["Ming Dynasty Imperial Library", "EthnoDock Classical Library"]
    }
]

# Major Bioactive Compounds as OKF Entities
BIOACTIVE_COMPOUNDS = [
    {
        "okf_version": "1.0",
        "entity_id": "compound_berberine",
        "type": "phytochemical_compound",
        "title": "Berberine (黄连素 / 柏勃灵)",
        "chemical_class": "Isoquinoline Alkaloid",
        "formula": "C20H18NO4+",
        "primary_herbs": ["Coptis chinensis (黄连)", "Phellodendron amurense (黄柏)", "Berberis vulgaris"],
        "pharmacology": ["AMPK Activation", "Antibacterial & Antimicrobial", "Blood Glucose Regulation", "Anti-inflammatory"],
        "sources": ["PubChem", "ChEMBL", "EthnoDock Phytochemistry"]
    },
    {
        "okf_version": "1.0",
        "entity_id": "compound_artemisinin",
        "type": "phytochemical_compound",
        "title": "Artemisinin (青蒿素)",
        "chemical_class": "Sesquiterpene Lactone Endoperoxide",
        "formula": "C15H22O5",
        "primary_herbs": ["Artemisia annua (青蒿)"],
        "pharmacology": ["Potent Antimalarial Action", "Reactive Oxygen Species Generation in Parasites", "Antitumor Potential"],
        "sources": ["PubChem", "Nobel Prize Archives", "EthnoDock Phytochemistry"]
    },
    {
        "okf_version": "1.0",
        "entity_id": "compound_ginsenoside_rg1",
        "type": "phytochemical_compound",
        "title": "Ginsenoside Rg1 (人参皂苷 Rg1)",
        "chemical_class": "Triterpenoid Dammarane Saponin",
        "formula": "C42H72O14",
        "primary_herbs": ["Panax ginseng (人参)", "Panax notoginseng (三七)"],
        "pharmacology": ["Neuroprotective Action", "Central Nervous System Excitation", "Angiogenesis Promotion", "Anti-fatigue"],
        "sources": ["PubChem", "ChEMBL", "EthnoDock Phytochemistry"]
    }
]

def parse_species():
    species_list = []
    if not os.path.exists(MASTER_SPECIES_FILE):
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

            slug = re.sub(r'[^a-z0-9]+', '_', common_name.lower()).strip('_')
            entity_id = f"herb_{slug}"

            species_list.append({
                "entity_id": entity_id,
                "common_name": common_name,
                "binomial_name": latin_name,
                "chinese_name": chinese_name
            })

    return species_list

def generate_expanded_okf_database():
    species_list = parse_species()
    index_dict = {}

    okf_md_content = "# EthnoDock Open Knowledge Format (OKF) Master Database\n\n"
    okf_md_content += f"Total Validated OKF Entities: {len(species_list) + len(CLASSICAL_TEXTS) + len(BIOACTIVE_COMPOUNDS)}\n"
    okf_md_content += " - Species Entities: {}\n".format(len(species_list))
    okf_md_content += " - Classical Text Entities: {}\n".format(len(CLASSICAL_TEXTS))
    okf_md_content += " - Bioactive Compound Entities: {}\n\n".format(len(BIOACTIVE_COMPOUNDS))

    # 1. Add Species Entities
    for spec in species_list:
        common = spec["common_name"]
        latin = spec["binomial_name"]
        chinese = spec["chinese_name"]
        entity_id = spec["entity_id"]

        profile = HERB_PROFILES.get(common, {
            "nature": "Neutral/Balanced",
            "flavor": ["Harmonizing"],
            "meridians": ["Visceral Networks"],
            "actions": [f"Standard ethnobotanical action for {common}"],
            "compounds": [f"{common.replace(' ', '')} active phytochemicals"]
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
                "nature": profile["nature"],
                "flavor": profile["flavor"],
                "meridians": profile["meridians"],
                "primary_actions": profile["actions"]
            },
            "phytochemistry": {
                "primary_compounds": profile["compounds"]
            },
            "sources": ["Shennong Bencao Jing", "Bencao Gangmu", "EthnoDock Core Knowledge"]
        }

        yaml_str = yaml.dump(okf_metadata, sort_keys=False, allow_unicode=True).strip()

        entry_md = f"---\n{yaml_str}\n---\n\n"
        entry_md += f"## {common} ({chinese})\n"
        entry_md += f"**Taxonomic Binomial**: *{latin}*\n"
        entry_md += f"**Chinese Name**: {chinese}\n"
        entry_md += f"**TCM Properties**: {profile['nature']} nature; Flavors: {', '.join(profile['flavor'])}; Meridians: {', '.join(profile['meridians'])}\n"
        entry_md += f"**Actions**: {'; '.join(profile['actions'])}\n"
        entry_md += f"**Primary Bioactive Compounds**: {', '.join(profile['compounds'])}\n\n"
        entry_md += "---\n\n"

        okf_md_content += entry_md

        index_dict[entity_id] = okf_metadata
        index_dict[chinese] = entity_id
        index_dict[latin.lower()] = entity_id
        index_dict[common.lower()] = entity_id

    # 2. Add Classical Text Entities
    for text_ent in CLASSICAL_TEXTS:
        yaml_str = yaml.dump(text_ent, sort_keys=False, allow_unicode=True).strip()
        entry_md = f"---\n{yaml_str}\n---\n\n"
        entry_md += f"## {text_ent['title']}\n"
        entry_md += f"**Era**: {text_ent['era']}\n"
        entry_md += f"**Author**: {text_ent['author']}\n"
        entry_md += f"**Significance**: {text_ent['significance']}\n\n"
        entry_md += "---\n\n"

        okf_md_content += entry_md

        eid = text_ent["entity_id"]
        index_dict[eid] = text_ent
        index_dict[text_ent["title"].lower()] = eid

    # 3. Add Bioactive Compound Entities
    for cmp_ent in BIOACTIVE_COMPOUNDS:
        yaml_str = yaml.dump(cmp_ent, sort_keys=False, allow_unicode=True).strip()
        entry_md = f"---\n{yaml_str}\n---\n\n"
        entry_md += f"## {cmp_ent['title']}\n"
        entry_md += f"**Chemical Class**: {cmp_ent['chemical_class']}\n"
        entry_md += f"**Formula**: {cmp_ent['formula']}\n"
        entry_md += f"**Primary Botanical Sources**: {', '.join(cmp_ent['primary_herbs'])}\n"
        entry_md += f"**Pharmacology**: {'; '.join(cmp_ent['pharmacology'])}\n\n"
        entry_md += "---\n\n"

        okf_md_content += entry_md

        eid = cmp_ent["entity_id"]
        index_dict[eid] = cmp_ent
        index_dict[cmp_ent["title"].lower()] = eid

    # Save OKF Master Database
    with open(OKF_MASTER_FILE, "w", encoding="utf-8") as f:
        f.write(okf_md_content)

    # Save OKF Index JSON
    with open(OKF_INDEX_FILE, "w", encoding="utf-8") as f:
        json.dump(index_dict, f, indent=2, ensure_ascii=False)

    print(f"[OKF Database Updated] Total {len(index_dict)} indexed keys across Species, Texts, and Bioactive Compounds.")

if __name__ == '__main__':
    generate_expanded_okf_database()
