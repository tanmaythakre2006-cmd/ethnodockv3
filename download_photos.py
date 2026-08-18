import os
import json
import urllib.request
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
GALLERY_PATH = os.path.join(BASE_DIR, "species_image_gallery.json")
OUTPUT_DIR = os.path.join(BASE_DIR, "species_images")

os.makedirs(OUTPUT_DIR, exist_ok=True)

HEADERS = {
    "User-Agent": "EthnoDockKnowledgeBot/4.0 (Computational Apothecary Project; contact@ethnodock.org) Python/3.11"
}

def resolve_and_download(filename, raw_url, output_path):
    if os.path.exists(output_path) and os.path.getsize(output_path) > 1000:
        return True, output_path

    # Try raw URL first
    try:
        req = urllib.request.Request(raw_url, headers=HEADERS)
        with urllib.request.urlopen(req, timeout=12) as resp:
            if resp.status == 200:
                data = resp.read()
                if len(data) > 500:
                    with open(output_path, "wb") as f:
                        f.write(data)
                    return True, output_path
    except Exception:
        pass

    # Extract filename for Wikipedia API resolution
    file_title = raw_url.split("/")[-1]
    file_title = urllib.parse.unquote(file_title)
    
    api_url = f"https://en.wikipedia.org/w/api.php?action=query&titles=File:{urllib.parse.quote(file_title)}&prop=imageinfo&iiprop=url&iiurlwidth=800&format=json"
    
    try:
        req = urllib.request.Request(api_url, headers=HEADERS)
        with urllib.request.urlopen(req, timeout=12) as resp:
            res_data = json.loads(resp.read().decode('utf-8'))
            pages = res_data.get('query', {}).get('pages', {})
            for page in pages.values():
                imageinfo = page.get('imageinfo', [])
                if imageinfo:
                    thumb_url = imageinfo[0].get('thumburl') or imageinfo[0].get('url')
                    if thumb_url:
                        img_req = urllib.request.Request(thumb_url, headers=HEADERS)
                        with urllib.request.urlopen(img_req, timeout=15) as img_resp:
                            if img_resp.status == 200:
                                data = img_resp.read()
                                if len(data) > 500:
                                    with open(output_path, "wb") as f:
                                        f.write(data)
                                    return True, output_path
    except Exception:
        pass

    return False, output_path

def main():
    if not os.path.exists(GALLERY_PATH):
        print(f"Gallery file not found: {GALLERY_PATH}")
        return

    with open(GALLERY_PATH, "r", encoding="utf-8") as f:
        gallery_data = json.load(f)

    download_tasks = []
    
    for species_key, item in gallery_data.items():
        specimen_code = item.get("specimen_code", "SPECIMEN").replace("/", "_")
        photos = item.get("photos", [])
        
        local_photos = []
        for idx, photo_url in enumerate(photos):
            ext = ".jpg"
            if photo_url.lower().endswith(".png"):
                ext = ".png"
            elif photo_url.lower().endswith(".jpeg"):
                ext = ".jpeg"
            
            filename = f"{specimen_code}_{idx+1}{ext}"
            filepath = os.path.join(OUTPUT_DIR, filename)
            local_photos.append(filepath)
            
            if not (os.path.exists(filepath) and os.path.getsize(filepath) > 1000):
                download_tasks.append((filename, photo_url, filepath))
            
        item["local_photos"] = local_photos

    if download_tasks:
        print(f"API Resolution Download Queue: {len(download_tasks)} photos remaining...")
        success_count = 0
        fail_count = 0

        with ThreadPoolExecutor(max_workers=5) as executor:
            futures = {
                executor.submit(resolve_and_download, task[0], task[1], task[2]): task 
                for task in download_tasks
            }
            
            for future in as_completed(futures):
                success, res_path = future.result()
                if success:
                    success_count += 1
                else:
                    fail_count += 1

                total_done = success_count + fail_count
                if total_done % 15 == 0 or total_done == len(download_tasks):
                    print(f"Progress: {total_done}/{len(download_tasks)} | Downloaded: {success_count} | Failed: {fail_count}")

    saved_files = [f for f in os.listdir(OUTPUT_DIR) if os.path.getsize(os.path.join(OUTPUT_DIR, f)) > 1000]

    with open(GALLERY_PATH, "w", encoding="utf-8") as f:
        json.dump(gallery_data, f, indent=2, ensure_ascii=False)

    print(f"Complete! Total local photos stored: {len(saved_files)} in '{OUTPUT_DIR}'.")

if __name__ == "__main__":
    main()
