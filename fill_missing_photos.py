import os
import json
import urllib.request
import urllib.parse

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
GALLERY_PATH = os.path.join(BASE_DIR, "species_image_gallery.json")
OUTPUT_DIR = os.path.join(BASE_DIR, "species_images")

os.makedirs(OUTPUT_DIR, exist_ok=True)

HEADERS = {
    "User-Agent": "EthnoDockBot/5.0 (TCM Research Project; contact@ethnodock.org) Python/3.11"
}

def fetch_search_photo(query_term):
    api_url = f"https://en.wikipedia.org/w/api.php?action=query&generator=search&gsrsearch={urllib.parse.quote(query_term)}&gsrlimit=1&prop=pageimages&piprop=thumbnail&pithumbsize=800&format=json"
    try:
        req = urllib.request.Request(api_url, headers=HEADERS)
        with urllib.request.urlopen(req, timeout=12) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            pages = data.get('query', {}).get('pages', {})
            for page in pages.values():
                thumb = page.get('thumbnail', {}).get('source')
                if thumb:
                    return thumb
    except Exception:
        pass
    return None

def download_file(url, output_path):
    try:
        req = urllib.request.Request(url, headers=HEADERS)
        with urllib.request.urlopen(req, timeout=15) as resp:
            if resp.status == 200:
                data = resp.read()
                if len(data) > 500:
                    with open(output_path, "wb") as f:
                        f.write(data)
                    return True
    except Exception:
        pass
    return False

def main():
    if not os.path.exists(GALLERY_PATH):
        print("Gallery path missing.")
        return

    with open(GALLERY_PATH, "r", encoding="utf-8") as f:
        gallery_data = json.load(f)

    missing_tasks = []
    
    for key, item in gallery_data.items():
        local_photos = item.get("local_photos", [])
        has_local = any(os.path.exists(p) and os.path.getsize(p) > 1000 for p in local_photos)
        
        if not has_local:
            scientific = item.get("scientific_name", "")
            common = item.get("common_name", "")
            code = item.get("specimen_code", "TCM").replace("/", "_")
            missing_tasks.append((key, scientific, common, code))

    print(f"Searching Wikipedia photos for {len(missing_tasks)} missing species...")
    
    success_count = 0

    for key, scientific, common, code in missing_tasks:
        search_query = scientific if scientific and scientific != "N/A" else common
        photo_url = fetch_search_photo(search_query)
        
        if not photo_url and common:
            photo_url = fetch_search_photo(common)

        if photo_url:
            ext = ".jpg"
            if photo_url.lower().endswith(".png"):
                ext = ".png"
            filename = f"{code}_missing{ext}"
            filepath = os.path.join(OUTPUT_DIR, filename)
            
            if download_file(photo_url, filepath):
                gallery_data[key]["local_photos"] = [filepath]
                gallery_data[key]["photos"] = [photo_url]
                success_count += 1
                print(f"[OK] Downloaded photo for: {common}")
            else:
                print(f"[FAIL] Download failed for: {common}")
        else:
            print(f"[FAIL] No image found for: {common}")

    with open(GALLERY_PATH, "w", encoding="utf-8") as f:
        json.dump(gallery_data, f, indent=2, ensure_ascii=False)

    print(f"\nFinished! Added photos for {success_count} missing species.")

if __name__ == "__main__":
    main()
