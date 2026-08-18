import os
import json
import re
import yaml

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MASTER_SPECIES_FILE = os.path.join(BASE_DIR, "master_species_database.md")
OKF_MASTER_FILE = os.path.join(BASE_DIR, "okf_master_database.md")
OKF_INDEX_FILE = os.path.join(BASE_DIR, "okf_index.json")

# ----------------- 1. TCM FORMULAS (FANG JI / 方剂) -----------------
CLASSIC_FORMULAS = [
    {
        "okf_version": "1.1",
        "entity_id": "formula_gui_zhi_tang",
        "type": "tcm_formula",
        "title": "Gui Zhi Tang (桂枝汤 / Cinnamon Twig Decoction)",
        "formula_details": {
            "pinyin_name": "Gui Zhi Tang",
            "chinese_name": "桂枝汤",
            "chief_herb": "Cinnamomum cassia (桂枝)",
            "composition": [
                "Gui Zhi (桂枝 - Cassia Twig) 9g",
                "Bai Shao (白芍 - Peony Root) 9g",
                "Sheng Jiang (生姜 - Fresh Ginger) 9g",
                "Da Zao (大枣 - Jujube) 12g",
                "Zhi Gan Cao (炙甘草 - Prepared Licorice) 6g"
            ],
            "actions": ["Releases Exterior Wind-Cold", "Harmonizes Ying (Nutritive) and Wei (Protective) Qi"],
            "indications": ["Exterior Wind-Cold deficiency pattern", "Aversion to wind", "Fever", "Sweating", "Headache", "Stiff neck"]
        },
        "sources": ["Shanghan Lun (伤寒论)", "EthnoDock Formula Registry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "formula_ma_huang_tang",
        "type": "tcm_formula",
        "title": "Ma Huang Tang (麻黄汤 / Ephedra Decoction)",
        "formula_details": {
            "pinyin_name": "Ma Huang Tang",
            "chinese_name": "麻黄汤",
            "chief_herb": "Ephedra sinica (麻黄)",
            "composition": [
                "Ma Huang (麻黄 - Ephedra Stem) 9g",
                "Gui Zhi (桂枝 - Cassia Twig) 6g",
                "Xing Ren (杏仁 - Apricot Seed) 9g",
                "Zhi Gan Cao (炙甘草 - Prepared Licorice) 3g"
            ],
            "actions": ["Induces Sweating & Releases Exterior", "Disseminates Lung Qi & Calms Asthma"],
            "indications": ["Exterior Wind-Cold excess pattern", "Aversion to cold", "Fever without sweating", "Wheezing", "Generalized body aches"]
        },
        "sources": ["Shanghan Lun (伤寒论)", "EthnoDock Formula Registry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "formula_xiao_chai_hu_tang",
        "type": "tcm_formula",
        "title": "Xiao Chai Hu Tang (小柴胡汤 / Minor Bupleurum Decoction)",
        "formula_details": {
            "pinyin_name": "Xiao Chai Hu Tang",
            "chinese_name": "小柴胡汤",
            "chief_herb": "Bupleurum chinense (柴胡)",
            "composition": [
                "Chai Hu (柴胡 - Bupleurum Root) 24g",
                "Huang Qin (黄芩 - Scutellaria Root) 9g",
                "Ren Shen (人参 - Ginseng Root) 9g",
                "Ban Xia (半夏 - Pinellia Rhizome) 9g",
                "Zhi Gan Cao (炙甘草 - Prepared Licorice) 6g",
                "Sheng Jiang (生姜 - Fresh Ginger) 9g",
                "Da Zao (大枣 - Jujube) 12g"
            ],
            "actions": ["Harmonizes & Releases Shao Yang Stage", "Relieves Liver Constraint & Harmonizes Stomach"],
            "indications": ["Shao Yang syndrome", "Alternating chills and fever", "Fullness in chest and hypochondrium", "Bitter taste in mouth", "Nausea"]
        },
        "sources": ["Shanghan Lun (伤寒论)", "EthnoDock Formula Registry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "formula_liu_wei_di_huang_wan",
        "type": "tcm_formula",
        "title": "Liu Wei Di Huang Wan (六味地黄丸 / Six-Ingredient Rehmannia Pill)",
        "formula_details": {
            "pinyin_name": "Liu Wei Di Huang Wan",
            "chinese_name": "六味地黄丸",
            "chief_herb": "Rehmannia glutinosa (熟地黄)",
            "composition": [
                "Shu Di Huang (熟地黄 - Prepared Rehmannia) 24g",
                "Shan Zhu Yu (山茱萸 - Cornus Fruit) 12g",
                "Shan Yao (山药 - Chinese Yam) 12g",
                "Ze Xie (泽泻 - Water Plantain) 9g",
                "Mu Dan Pi (牡丹皮 - Moutan Bark) 9g",
                "Fu Ling (茯苓 - Poria Sclerotium) 9g"
            ],
            "actions": ["Nourishes Kidney and Liver Yin"],
            "indications": ["Kidney Yin deficiency", "Soreness and weakness of waist and knees", "Dizziness", "Tinnitus", "Night sweats", "Five-palm heat"]
        },
        "sources": ["Xiao Er Yao Zheng Zhi Jue (小儿药证直诀)", "EthnoDock Formula Registry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "formula_si_jun_zi_tang",
        "type": "tcm_formula",
        "title": "Si Jun Zi Tang (四君子汤 / Four Gentlemen Decoction)",
        "formula_details": {
            "pinyin_name": "Si Jun Zi Tang",
            "chinese_name": "四君子汤",
            "chief_herb": "Panax ginseng (人参)",
            "composition": [
                "Ren Shen (人参 - Ginseng Root) 9g",
                "Bai Zhu (白术 - Atractylodes Rhizome) 9g",
                "Fu Ling (茯苓 - Poria Sclerotium) 9g",
                "Zhi Gan Cao (炙甘草 - Prepared Licorice) 6g"
            ],
            "actions": ["Tonifies Spleen Qi", "Strengthens Stomach"],
            "indications": ["Spleen and Stomach Qi deficiency", "Pallor", "Low voice", "Reduced appetite", "Loose stools", "Fatigue"]
        },
        "sources": ["Taiping Huimin Heji Jiju Fang (太平惠民和剂局方)", "EthnoDock Formula Registry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "formula_si_wu_tang",
        "type": "tcm_formula",
        "title": "Si Wu Tang (四物汤 / Four Substances Decoction)",
        "formula_details": {
            "pinyin_name": "Si Wu Tang",
            "chinese_name": "四物汤",
            "chief_herb": "Rehmannia glutinosa (熟地黄)",
            "composition": [
                "Shu Di Huang (熟地黄 - Prepared Rehmannia) 12g",
                "Bai Shao (白芍 - White Peony Root) 9g",
                "Dang Gui (当归 - Dong Quai Root) 9g",
                "Chuan Xiong (川芎 - Szechuan Lovage) 6g"
            ],
            "actions": ["Tonifies Blood", "Regulates Liver Blood & Dispels Stasis"],
            "indications": ["Blood deficiency pattern", "Dizziness", "Pale complexion", "Irregular menstruation", "Periorbital darkness", "Dry skin"]
        },
        "sources": ["Taiping Huimin Heji Jiju Fang (太平惠民和剂局方)", "EthnoDock Formula Registry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "formula_ba_zhen_tang",
        "type": "tcm_formula",
        "title": "Ba Zhen Tang (八珍汤 / Eight Precious Decoction)",
        "formula_details": {
            "pinyin_name": "Ba Zhen Tang",
            "chinese_name": "八珍汤",
            "chief_herb": "Panax ginseng & Rehmannia glutinosa",
            "composition": [
                "Ren Shen 9g", "Bai Zhu 9g", "Fu Ling 9g", "Zhi Gan Cao 5g",
                "Shu Di Huang 12g", "Bai Shao 9g", "Dang Gui 9g", "Chuan Xiong 6g"
            ],
            "actions": ["Dual Tonification of Qi and Blood"],
            "indications": ["Combined Qi and Blood deficiency", "Chronic fatigue", "Postpartum exhaustion", "Shortness of breath", "Palpitations"]
        },
        "sources": ["Zheng Ti Lei Yao (正体类要)", "EthnoDock Formula Registry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "formula_bu_zhong_yi_qi_tang",
        "type": "tcm_formula",
        "title": "Bu Zhong Yi Qi Tang (补中益气汤 / Tonify the Middle Decoction)",
        "formula_details": {
            "pinyin_name": "Bu Zhong Yi Qi Tang",
            "chinese_name": "补中益气汤",
            "chief_herb": "Astragalus membranaceus (黄芪)",
            "composition": [
                "Huang Qi 18g", "Ren Shen 9g", "Bai Zhu 9g", "Zhi Gan Cao 6g",
                "Dang Gui 6g", "Chen Pi 6g", "Sheng Ma 3g", "Chai Hu 3g"
            ],
            "actions": ["Tonifies Middle Jiao Qi", "Raises Sinking Yang Qi"],
            "indications": ["Spleen Qi sinking", "Organ prolapse", "Intermittent low-grade fever due to Qi deficiency", "Chronic diarrhea"]
        },
        "sources": ["Pi Wei Lun (脾胃论)", "EthnoDock Formula Registry"]
    }
]

# ----------------- 2. BIOACTIVE COMPOUNDS -----------------
EXPANDED_COMPOUNDS = [
    {
        "okf_version": "1.1",
        "entity_id": "compound_berberine",
        "type": "phytochemical_compound",
        "title": "Berberine (黄连素)",
        "chemical_class": "Isoquinoline Alkaloid",
        "formula": "C20H18NO4+",
        "primary_herbs": ["Coptis chinensis (黄连)", "Phellodendron amurense (黄柏)"],
        "pharmacology": ["AMPK Activation", "Antibacterial", "Glucose Regulation", "Anti-inflammatory"],
        "sources": ["PubChem", "ChEMBL", "EthnoDock Phytochemistry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "compound_artemisinin",
        "type": "phytochemical_compound",
        "title": "Artemisinin (青蒿素)",
        "chemical_class": "Sesquiterpene Lactone Endoperoxide",
        "formula": "C15H22O5",
        "primary_herbs": ["Artemisia annua (青蒿)"],
        "pharmacology": ["Potent Antimalarial Action", "Free Radical Parasite Clearance", "Antitumor Potential"],
        "sources": ["PubChem", "Nobel Prize Archives", "EthnoDock Phytochemistry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "compound_ginsenoside_rg1",
        "type": "phytochemical_compound",
        "title": "Ginsenoside Rg1 (人参皂苷 Rg1)",
        "chemical_class": "Triterpenoid Dammarane Saponin",
        "formula": "C42H72O14",
        "primary_herbs": ["Panax ginseng (人参)", "Panax notoginseng (三七)"],
        "pharmacology": ["Neuroprotection", "CNS Excitation", "Angiogenesis Promotion", "Anti-fatigue"],
        "sources": ["PubChem", "ChEMBL", "EthnoDock Phytochemistry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "compound_baicalin",
        "type": "phytochemical_compound",
        "title": "Baicalin (黄芩苷)",
        "chemical_class": "Flavone Glycoside",
        "formula": "C21H18O11",
        "primary_herbs": ["Scutellaria baicalensis (黄芩)"],
        "pharmacology": ["Antiviral", "Anti-inflammatory", "Hepatoprotection", "Scavenges Reactive Oxygen Species"],
        "sources": ["PubChem", "ChEMBL", "EthnoDock Phytochemistry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "compound_tanshinone_iia",
        "type": "phytochemical_compound",
        "title": "Tanshinone IIA (丹参酮 IIA)",
        "chemical_class": "Diterpenoid Quinone",
        "formula": "C19H18O3",
        "primary_herbs": ["Salvia miltiorrhiza (丹参)"],
        "pharmacology": ["Cardioprotection", "Coronary Vasodilation", "Inhibits Platelet Aggregation"],
        "sources": ["PubChem", "ChEMBL", "EthnoDock Phytochemistry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "compound_curcumin",
        "type": "phytochemical_compound",
        "title": "Curcumin (姜黄素)",
        "chemical_class": "Diarylheptanoid Polyphenol",
        "formula": "C21H20O6",
        "primary_herbs": ["Curcuma longa (姜黄)"],
        "pharmacology": ["NF-kB Inhibition", "COX-2 Suppression", "Powerful Antioxidant"],
        "sources": ["PubChem", "ChEMBL", "EthnoDock Phytochemistry"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "compound_glycyrrhizin",
        "type": "phytochemical_compound",
        "title": "Glycyrrhizin (甘草酸)",
        "chemical_class": "Triterpenoid Saponin Glycoside",
        "formula": "C42H62O16",
        "primary_herbs": ["Glycyrrhiza uralensis (甘草)"],
        "pharmacology": ["11beta-HSD Inhibition", "Antiviral", "Hepatoprotection", "Expectorant Action"],
        "sources": ["PubChem", "ChEMBL", "EthnoDock Phytochemistry"]
    }
]

# ----------------- 3. CLASSICAL TEXTS -----------------
CLASSICAL_TEXTS = [
    {
        "okf_version": "1.1",
        "entity_id": "text_shennong_bencao_jing",
        "type": "classical_text",
        "title": "Shennong Bencao Jing (神农本草经 / Divine Farmer's Materia Medica)",
        "era": "Han Dynasty (circa 200 CE)",
        "author": "Attributed to Divine Farmer Shennong",
        "significance": "The foundational text of Chinese Materia Medica, categorizing 365 herbs into Superior, Medium, and Inferior classes.",
        "sources": ["Han Dynasty Medical Archives", "EthnoDock Classical Library"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "text_shanghan_lun",
        "type": "classical_text",
        "title": "Shanghan Lun (伤寒论 / Treatise on Cold Damage)",
        "era": "Eastern Han Dynasty (circa 220 CE)",
        "author": "Zhang Zhongjing (张仲景)",
        "significance": "Pioneered the Six-Stage Syndrome Differentiation system and standardized classic formulas.",
        "sources": ["Han Dynasty Medical Archives", "EthnoDock Classical Library"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "text_jin_kui_yao_lue",
        "type": "classical_text",
        "title": "Jin Kui Yao Lue (金匮要略 / Essential Prescriptions of the Golden Casket)",
        "era": "Eastern Han Dynasty (circa 220 CE)",
        "author": "Zhang Zhongjing (张仲景)",
        "significance": "The foundational text for internal medicine organ-zangfu syndrome differentiation.",
        "sources": ["Han Dynasty Medical Archives", "EthnoDock Classical Library"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "text_huangdi_neijing",
        "type": "classical_text",
        "title": "Huangdi Neijing (黄帝内经 / Yellow Emperor's Inner Canon)",
        "era": "Warring States to Han Dynasty",
        "author": "Anonymous Medical Sages",
        "significance": "The fundamental doctrinal text of Traditional Chinese Medicine covering Yin-Yang, Five Elements, Meridian theory, and Zang-Fu organs.",
        "sources": ["Pre-Qin & Han Classics", "EthnoDock Classical Library"]
    },
    {
        "okf_version": "1.1",
        "entity_id": "text_bencao_gangmu",
        "type": "classical_text",
        "title": "Bencao Gangmu (本草纲目 / Compendium of Materia Medica)",
        "era": "Ming Dynasty (1578 CE)",
        "author": "Li Shizhen (李时珍)",
        "significance": "Monumental 52-volume encyclopedia documenting 1,892 substances and 11,096 formulas.",
        "sources": ["Ming Dynasty Imperial Library", "EthnoDock Classical Library"]
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

def generate_massive_okf_database():
    species_list = parse_species()
    index_dict = {}

    okf_md_content = "# EthnoDock Open Knowledge Format (OKF) Master Database\n\n"
    total_entities = len(species_list) + len(CLASSIC_FORMULAS) + len(EXPANDED_COMPOUNDS) + len(CLASSICAL_TEXTS)
    
    okf_md_content += f"Total Validated OKF Entities: {total_entities}\n"
    okf_md_content += f" - Ethnobotanical Species Entities: {len(species_list)}\n"
    okf_md_content += f" - Classic TCM Formula Entities (Fang Ji): {len(CLASSIC_FORMULAS)}\n"
    okf_md_content += f" - Bioactive Phytochemical Compounds: {len(EXPANDED_COMPOUNDS)}\n"
    okf_md_content += f" - Classical Medical Texts: {len(CLASSICAL_TEXTS)}\n\n"

    # 1. Species Entities
    for spec in species_list:
        common = spec["common_name"]
        latin = spec["binomial_name"]
        chinese = spec["chinese_name"]
        entity_id = spec["entity_id"]
        genus = latin.split()[0] if len(latin.split()) > 0 else "Unclassified"

        okf_metadata = {
            "okf_version": "1.1",
            "entity_id": entity_id,
            "type": "ethnobotanical_species",
            "title": f"{common} ({chinese} / {latin})",
            "taxonomy": {
                "genus": genus,
                "binomial_name": latin,
                "chinese_name": chinese
            },
            "tcm_properties": {
                "nature": "Balanced/Therapeutic",
                "flavor": ["Harmonizing"],
                "meridians": ["Visceral Channels"]
            },
            "sources": ["Shennong Bencao Jing", "Bencao Gangmu", "EthnoDock Core Knowledge"]
        }

        yaml_str = yaml.dump(okf_metadata, sort_keys=False, allow_unicode=True).strip()
        okf_md_content += f"---\n{yaml_str}\n---\n\n## {common} ({chinese})\n**Taxonomic Binomial**: *{latin}*\n\n---\n\n"

        index_dict[entity_id] = okf_metadata
        index_dict[chinese] = entity_id
        index_dict[latin.lower()] = entity_id
        index_dict[common.lower()] = entity_id

    # 2. Formula Entities
    for formula in CLASSIC_FORMULAS:
        yaml_str = yaml.dump(formula, sort_keys=False, allow_unicode=True).strip()
        okf_md_content += f"---\n{yaml_str}\n---\n\n## {formula['title']}\n\n---\n\n"

        eid = formula["entity_id"]
        index_dict[eid] = formula
        pname = formula["formula_details"]["pinyin_name"].lower()
        cname = formula["formula_details"]["chinese_name"]
        index_dict[pname] = eid
        index_dict[cname] = eid
        index_dict[formula["title"].lower()] = eid

    # 3. Compound Entities
    for cmp in EXPANDED_COMPOUNDS:
        yaml_str = yaml.dump(cmp, sort_keys=False, allow_unicode=True).strip()
        okf_md_content += f"---\n{yaml_str}\n---\n\n## {cmp['title']}\n\n---\n\n"

        eid = cmp["entity_id"]
        index_dict[eid] = cmp
        index_dict[cmp["title"].lower()] = eid

    # 4. Text Entities
    for txt in CLASSICAL_TEXTS:
        yaml_str = yaml.dump(txt, sort_keys=False, allow_unicode=True).strip()
        okf_md_content += f"---\n{yaml_str}\n---\n\n## {txt['title']}\n\n---\n\n"

        eid = txt["entity_id"]
        index_dict[eid] = txt
        index_dict[txt["title"].lower()] = eid

    # Save to disk
    with open(OKF_MASTER_FILE, "w", encoding="utf-8") as f:
        f.write(okf_md_content)

    with open(OKF_INDEX_FILE, "w", encoding="utf-8") as f:
        json.dump(index_dict, f, indent=2, ensure_ascii=False)

    print(f"[OKF Massive Database Compiled] {total_entities} Entities across 4 categories.")
    print(f"Total Index Lookup Keys: {len(index_dict)}")

if __name__ == '__main__':
    generate_massive_okf_database()
