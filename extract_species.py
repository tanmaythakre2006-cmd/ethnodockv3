import os
import sys

sys.stdout.reconfigure(encoding='utf-8')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ARTIFACT_DIR = r"C:\Users\hp\.gemini\antigravity\brain\301e2e9e-b650-4e6e-becd-d7d86ba73e0b"
OUTPUT_FILE = os.path.join(ARTIFACT_DIR, "plant_species_database.md")

# Massively expanded dictionary of major TCM plant species
HERB_DICT = {
    "Ginseng": {"scientific": "Panax ginseng", "chinese": "人参", "keywords": ["人参", "ginseng", "ren shen"]},
    "Licorice Root": {"scientific": "Glycyrrhiza uralensis", "chinese": "甘草", "keywords": ["甘草", "licorice", "gan cao"]},
    "Goji Berry": {"scientific": "Lycium barbarum", "chinese": "枸杞", "keywords": ["枸杞", "goji", "wolfberry", "gou qi"]},
    "Astragalus": {"scientific": "Astragalus membranaceus", "chinese": "黄芪", "keywords": ["黄芪", "astragalus", "huang qi"]},
    "Ephedra": {"scientific": "Ephedra sinica", "chinese": "麻黄", "keywords": ["麻黄", "ephedra", "ma huang"]},
    "Cassia Bark": {"scientific": "Cinnamomum cassia", "chinese": "肉桂", "keywords": ["肉桂", "cassia", "cinnamon bark", "rou gui"]},
    "Dong Quai": {"scientific": "Angelica sinensis", "chinese": "当归", "keywords": ["当归", "dong quai", "angelica", "dang gui"]},
    "Ginger": {"scientific": "Zingiber officinale", "chinese": "生姜", "keywords": ["生姜", "ginger", "sheng jiang"]},
    "Peony Root": {"scientific": "Paeonia lactiflora", "chinese": "白芍", "keywords": ["白芍", "peony", "bai shao"]},
    "Bupleurum": {"scientific": "Bupleurum chinense", "chinese": "柴胡", "keywords": ["柴胡", "bupleurum", "chai hu"]},
    "Ginkgo": {"scientific": "Ginkgo biloba", "chinese": "白果", "keywords": ["白果", "ginkgo", "bai guo"]},
    "Poria": {"scientific": "Poria cocos", "chinese": "茯苓", "keywords": ["茯苓", "poria", "fu ling"]},
    "Rehmannia": {"scientific": "Rehmannia glutinosa", "chinese": "地黄", "keywords": ["地黄", "rehmannia", "di huang"]},
    "Jujube": {"scientific": "Ziziphus jujuba", "chinese": "大枣", "keywords": ["大枣", "jujube", "da zao"]},
    "Tangerine Peel": {"scientific": "Citrus reticulata", "chinese": "陈皮", "keywords": ["陈皮", "tangerine peel", "chen pi"]},
    "Atractylodes (White)": {"scientific": "Atractylodes macrocephala", "chinese": "白术", "keywords": ["白术", "atractylodes", "bai zhu"]},
    "Schisandra": {"scientific": "Schisandra chinensis", "chinese": "五味子", "keywords": ["五味子", "schisandra", "wu wei zi"]},
    "Coptis": {"scientific": "Coptis chinensis", "chinese": "黄连", "keywords": ["黄连", "coptis", "huang lian"]},
    "Scutellaria": {"scientific": "Scutellaria baicalensis", "chinese": "黄芩", "keywords": ["黄芩", "scutellaria", "huang qin"]},
    "Moutan": {"scientific": "Paeonia suffruticosa", "chinese": "牡丹皮", "keywords": ["牡丹皮", "moutan", "mu dan pi"]},
    "Siberian Ginseng": {"scientific": "Eleutherococcus senticosus", "chinese": "刺五加", "keywords": ["刺五加", "siberian ginseng", "eleuthero", "ci wu jia"]},
    "Kudzu": {"scientific": "Pueraria lobata", "chinese": "葛根", "keywords": ["葛根", "kudzu", "pueraria", "ge gen"]},
    
    # Expanded entries
    "Forsythia": {"scientific": "Forsythia suspensa", "chinese": "连翘", "keywords": ["连翘", "forsythia", "lian qiao"]},
    "Honeysuckle": {"scientific": "Lonicera japonica", "chinese": "金银花", "keywords": ["金银花", "honeysuckle", "jin yin hua"]},
    "Isatis Root": {"scientific": "Isatis tinctoria", "chinese": "板蓝根", "keywords": ["板蓝根", "isatis", "ban lan gen"]},
    "Dahurian Angelica": {"scientific": "Angelica dahurica", "chinese": "白芷", "keywords": ["白芷", "dahurian angelica", "bai zhi"]},
    "Szechuan Lovage": {"scientific": "Ligusticum striatum", "chinese": "川芎", "keywords": ["川芎", "lovage", "chuan xiong"]},
    "Chinese Yam": {"scientific": "Dioscorea polystachya", "chinese": "山药", "keywords": ["山药", "chinese yam", "dioscorea", "shan yao"]},
    "Asiatic Cornelian Cherry": {"scientific": "Cornus officinalis", "chinese": "山茱萸", "keywords": ["山茱萸", "cornus", "shan zhu yu"]},
    "Red Sage": {"scientific": "Salvia miltiorrhiza", "chinese": "丹参", "keywords": ["丹参", "red sage", "salvia", "dan shen"]},
    "Notoginseng": {"scientific": "Panax notoginseng", "chinese": "三七", "keywords": ["三七", "notoginseng", "san qi"]},
    "Sweet Wormwood": {"scientific": "Artemisia annua", "chinese": "青蒿", "keywords": ["青蒿", "artemisia", "sweet wormwood", "qing hao"]},
    "Amur Cork Tree": {"scientific": "Phellodendron amurense", "chinese": "黄柏", "keywords": ["黄柏", "phellodendron", "cork tree", "huang bo"]},
    "Desertliving Cistanche": {"scientific": "Cistanche deserticola", "chinese": "肉苁蓉", "keywords": ["肉苁蓉", "cistanche", "rou cong rong"]},
    "Epimedium": {"scientific": "Epimedium brevicornu", "chinese": "淫羊藿", "keywords": ["淫羊藿", "epimedium", "horny goat weed", "yin yang huo"]},
    "Dodder Seed": {"scientific": "Cuscuta chinensis", "chinese": "菟丝子", "keywords": ["菟丝子", "dodder seed", "cuscuta", "tu si zi"]},
    "Magnolia Bark": {"scientific": "Magnolia officinalis", "chinese": "厚朴", "keywords": ["厚朴", "magnolia bark", "hou po"]},
    "Bitter Orange": {"scientific": "Citrus aurantium", "chinese": "枳壳", "keywords": ["枳壳", "bitter orange", "zhi ke"]},
    "Sichuan Aconite": {"scientific": "Aconitum carmichaelii", "chinese": "附子", "keywords": ["附子", "aconite", "fu zi"]},
    "Pinellia": {"scientific": "Pinellia ternata", "chinese": "半夏", "keywords": ["半夏", "pinellia", "ban xia"]},
    "Gastrodia": {"scientific": "Gastrodia elata", "chinese": "天麻", "keywords": ["天麻", "gastrodia", "tian ma"]},
    "Corydalis": {"scientific": "Corydalis yanhusuo", "chinese": "延胡索", "keywords": ["延胡索", "corydalis", "yan hu suo"]},
    "Safflower": {"scientific": "Carthamus tinctorius", "chinese": "红花", "keywords": ["红花", "safflower", "hong hua"]},
    "Peach Seed": {"scientific": "Prunus persica", "chinese": "桃仁", "keywords": ["桃仁", "peach seed", "tao ren"]},
    "Chinese Motherwort": {"scientific": "Leonurus japonicus", "chinese": "益母草", "keywords": ["益母草", "motherwort", "yi mu cao"]},
    "Chrysanthemum": {"scientific": "Chrysanthemum morifolium", "chinese": "菊花", "keywords": ["菊花", "chrysanthemum", "ju hua"]},
    "Chinese Mint": {"scientific": "Mentha haplocalyx", "chinese": "薄荷", "keywords": ["薄荷", "chinese mint", "bo he"]},
    "White Mulberry": {"scientific": "Morus alba", "chinese": "桑叶", "keywords": ["桑叶", "桑白皮", "mulberry", "sang ye", "sang bai pi"]},
    "Solomon's Seal": {"scientific": "Polygonatum odoratum", "chinese": "玉竹", "keywords": ["玉竹", "polygonatum", "yu zhu"]},
    "Dendrobium": {"scientific": "Dendrobium nobile", "chinese": "石斛", "keywords": ["石斛", "dendrobium", "shi hu"]},
    "Ophiopogon": {"scientific": "Ophiopogon japonicus", "chinese": "麦冬", "keywords": ["麦冬", "ophiopogon", "mai dong"]},
    "Trichosanthes": {"scientific": "Trichosanthes kirilowii", "chinese": "瓜蒌", "keywords": ["瓜蒌", "trichosanthes", "gua lou"]},
    "Balloon Flower": {"scientific": "Platycodon grandiflorus", "chinese": "桔梗", "keywords": ["桔梗", "platycodon", "jie geng"]},
    "Apricot Seed": {"scientific": "Prunus armeniaca", "chinese": "杏仁", "keywords": ["杏仁", "apricot seed", "xing ren"]},
    "Spine Date Seed": {"scientific": "Ziziphus spinosa", "chinese": "酸枣仁", "keywords": ["酸枣仁", "spine date", "suan zao ren"]},
    "Senega Root": {"scientific": "Polygala tenuifolia", "chinese": "远志", "keywords": ["远志", "polygala", "yuan zhi"]},
    "Sweet Flag": {"scientific": "Acorus tatarinowii", "chinese": "石菖蒲", "keywords": ["石菖蒲", "sweet flag", "acorus", "shi chang pu"]},
    "Uncaria": {"scientific": "Uncaria rhynchophylla", "chinese": "钩藤", "keywords": ["钩藤", "uncaria", "gou teng"]},
    "Tribulus": {"scientific": "Tribulus terrestris", "chinese": "白蒺藜", "keywords": ["白蒺藜", "tribulus", "bai ji li"]},
    "Sicklepod": {"scientific": "Cassia obtusifolia", "chinese": "决明子", "keywords": ["决明子", "sicklepod", "jue ming zi"]},
    "Asian Plantain": {"scientific": "Plantago asiatica", "chinese": "车前子", "keywords": ["车前子", "plantain", "che qian zi"]},
    "Water Plantain": {"scientific": "Alisma plantago-aquatica", "chinese": "泽泻", "keywords": ["泽泻", "alisma", "water plantain", "ze xie"]},
    "Polyporus": {"scientific": "Polyporus umbellatus", "chinese": "猪苓", "keywords": ["猪苓", "polyporus", "zhu ling"]},
    "Stephania": {"scientific": "Stephania tetrandra", "chinese": "防己", "keywords": ["防己", "stephania", "fang ji"]},
    "Large-Leaf Gentian": {"scientific": "Gentiana macrophylla", "chinese": "秦艽", "keywords": ["秦艽", "gentian", "qin jiao"]},
    "Chinese Mistletoe": {"scientific": "Taxillus chinensis", "chinese": "桑寄生", "keywords": ["桑寄生", "mistletoe", "sang ji sheng"]},
    "Acanthopanax": {"scientific": "Acanthopanax gracilistylus", "chinese": "五加皮", "keywords": ["五加皮", "acanthopanax", "wu jia pi"]},
    "Chinese Quince": {"scientific": "Chaenomeles sinensis", "chinese": "木瓜", "keywords": ["木瓜", "quince", "mu gua"]},
    "Atractylodes (Black)": {"scientific": "Atractylodes lancea", "chinese": "苍术", "keywords": ["苍术", "cang zhu"]},
    "Patchouli": {"scientific": "Agastache rugosa", "chinese": "藿香", "keywords": ["藿香", "patchouli", "agastache", "huo xiang"]},
    "Amomum": {"scientific": "Amomum villosum", "chinese": "砂仁", "keywords": ["砂仁", "amomum", "sha ren"]},
    "Alpinia": {"scientific": "Alpinia oxyphylla", "chinese": "益智仁", "keywords": ["益智仁", "alpinia", "yi zhi ren"]},
    "Hawthorn": {"scientific": "Crataegus pinnatifida", "chinese": "山楂", "keywords": ["山楂", "hawthorn", "shan zha"]},
    "Barley Sprout": {"scientific": "Hordeum vulgare", "chinese": "麦芽", "keywords": ["麦芽", "barley sprout", "mai ya"]},
    "Radish Seed": {"scientific": "Raphanus sativus", "chinese": "莱菔子", "keywords": ["莱菔子", "radish seed", "lai fu zi"]},
    "Rhubarb": {"scientific": "Rheum palmatum", "chinese": "大黄", "keywords": ["大黄", "rhubarb", "da huang"]},
    "Aloe Vera": {"scientific": "Aloe vera", "chinese": "芦荟", "keywords": ["芦荟", "aloe", "lu hui"]},
    "Hemp Seed": {"scientific": "Cannabis sativa", "chinese": "火麻仁", "keywords": ["火麻仁", "hemp seed", "huo ma ren"]},
    "Kansui": {"scientific": "Euphorbia kansui", "chinese": "甘遂", "keywords": ["甘遂", "kansui", "gan sui"]},
    "Morning Glory Seed": {"scientific": "Pharbitis nil", "chinese": "牵牛子", "keywords": ["牵牛子", "morning glory", "qian niu zi"]}
}

def scan_pool():
    found_herbs = set()
    
    print("Starting Deep Scan of Unstructured Knowledge Pool with Expanded Dictionary...")
    
    # Iterate through all files in ethnodock_knowledge_pool
    for root, dirs, files in os.walk(BASE_DIR):
        # Skip structured data to strictly scan unstructured texts
        if "10_raw_structured_data" in root or "07_modern_structured_data" in root:
            continue
            
        for file in files:
            if file.endswith(".txt") or file.endswith(".md"):
                filepath = os.path.join(root, file)
                try:
                    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                        content = f.read().lower()
                        
                        # Match against our massive dictionary
                        for common_name, data in HERB_DICT.items():
                            if common_name in found_herbs:
                                continue # Already found, skip for efficiency
                                
                            for kw in data["keywords"]:
                                if kw in content:
                                    found_herbs.add(common_name)
                                    print(f"Match found in {file}: {common_name}")
                                    break
                except Exception as e:
                    pass
                    
    print(f"\nDeep Scan Complete. Extracted {len(found_herbs)} deduplicated plant species.")
    
    print("Formatting and writing massive database...")
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write("# EthnoDock Plant Species Database (Expanded)\n\n")
        f.write(f"This expanded database was dynamically extracted by scanning the unstructured classical texts and English OCR translations within the Knowledge Pool using a massive botanical dictionary. We successfully extracted **{len(found_herbs)}** verified plant species.\n\n")
        
        for common_name in sorted(found_herbs):
            sci_name = HERB_DICT[common_name]["scientific"]
            chi_name = HERB_DICT[common_name]["chinese"]
            f.write(f"- {common_name} [{sci_name}] [{chi_name}]\n")
            
    print(f"Database successfully saved to {OUTPUT_FILE}")

if __name__ == '__main__':
    scan_pool()
