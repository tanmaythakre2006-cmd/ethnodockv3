import os
import sys
import re
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SCRATCH_DIR = r"C:\Users\hp\.gemini\antigravity\brain\301e2e9e-b650-4e6e-becd-d7d86ba73e0b\scratch"
os.makedirs(SCRATCH_DIR, exist_ok=True)
CANDIDATES_FILE = os.path.join(SCRATCH_DIR, "candidates.txt")

# The original 80 species to filter out
KNOWN_HERBS = {
    "人参", "甘草", "枸杞", "黄芪", "麻黄", "肉桂", "当归", "生姜", "白芍", "柴胡", 
    "白果", "茯苓", "地黄", "大枣", "陈皮", "白术", "五味子", "黄连", "黄芩", "牡丹皮", 
    "刺五加", "葛根", "连翘", "金银花", "板蓝根", "白芷", "川芎", "山药", "山茱萸", 
    "丹参", "三七", "青蒿", "黄柏", "肉苁蓉", "淫羊藿", "菟丝子", "厚朴", "枳壳", "附子", 
    "半夏", "天麻", "延胡索", "红花", "桃仁", "益母草", "菊花", "薄荷", "桑叶", "玉竹", 
    "石斛", "麦冬", "瓜蒌", "桔梗", "杏仁", "酸枣仁", "远志", "石菖蒲", "钩藤", "白蒺藜", 
    "决明子", "车前子", "泽泻", "猪苓", "防己", "秦艽", "桑寄生", "五加皮", "木瓜", "苍术", 
    "藿香", "砂仁", "益智仁", "山楂", "麦芽", "莱菔子", "大黄", "芦荟", "火麻仁", "甘遂", "牵牛子"
}

def discover_unknowns():
    print("Initiating RAG Discovery Scanner...")
    candidates_counter = Counter()
    
    # Regex to capture 1 or 2 Chinese characters followed by a botanical suffix
    # e.g., 夏枯草 (Xia Ku Cao), 蒲公英 (wait, suffix doesn't match here, but it captures many)
    # Suffixes: 草(herb), 根(root), 叶(leaf), 皮(bark/peel), 子(seed), 仁(kernel), 花(flower), 枝(branch)
    pattern = re.compile(r'([\u4e00-\u9fa5]{1,2}(?:草|根|叶|皮|子|仁|花|枝))')
    
    for root, dirs, files in os.walk(BASE_DIR):
        if "10_raw_structured_data" in root or "07_modern_structured_data" in root:
            continue
            
        for file in files:
            if file.endswith(".txt") or file.endswith(".md"):
                filepath = os.path.join(root, file)
                try:
                    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                        content = f.read()
                        matches = pattern.findall(content)
                        for match in matches:
                            if match not in KNOWN_HERBS:
                                candidates_counter[match] += 1
                except Exception as e:
                    pass
                    
    # Filter out extremely common non-herb words that match this pattern
    # e.g., 儿子 (son), 银子 (silver), 妻子 (wife), 树根 (tree root), 草根 (grass root), 落叶 (fallen leaves)
    # The RAG agent will do the final filtering, but we can drop obvious noise.
    
    top_candidates = candidates_counter.most_common(100)
    print(f"Extraction complete. Found {len(candidates_counter)} unique linguistic candidates.")
    
    with open(CANDIDATES_FILE, "w", encoding="utf-8") as f:
        f.write("TOP 100 RAG CANDIDATES (UNKNOWN SPECIES)\n")
        f.write("========================================\n")
        for word, count in top_candidates:
            f.write(f"{word}: {count}\n")
            
    print(f"Candidates successfully dumped to {CANDIDATES_FILE} for RAG Evaluation.")

if __name__ == '__main__':
    discover_unknowns()
