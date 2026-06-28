from PIL import Image
from pathlib import Path
import sys

def crop_transparent(filepath):
    img = Image.open(filepath).convert("RGBA")
    bbox = img.getbbox()
    if bbox:
        img = img.crop(bbox)
        img.save(filepath)
        print(f"Cropped: {filepath}")
    else:
        print(f"Skipped (empty): {filepath}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        files = [Path(f) for f in sys.argv[1:]]
    else:
        files = sorted(Path(".").glob("*.png"))
    
    for f in files:
        crop_transparent(f)
