import requests
import os

os.makedirs("data", exist_ok=True)

MAP_KEY = "f9e70c64a2c21ce005f23551232d9315"  # paste your key from the FIRMS email

west, south, east, north = 79, 27, 88, 31
date = "2024-04-01"
day_range = 2

url = (
    f"https://firms.modaps.eosdis.nasa.gov/api/area/csv/"
    f"{MAP_KEY}/VIIRS_SNPP_NRT/"
    f"{west},{south},{east},{north}/"
    f"{day_range}/{date}"
)

print(f"Fetching FIRMS data ...")
response = requests.get(url, timeout=60)

if response.status_code == 200:
    with open("data/firms_viirs_active_fire.csv", "w") as f:
        f.write(response.text)
    lines = response.text.strip().split("\n")
    print(f"Downloaded {len(lines) - 1} fire detections")
    print("Done: data/firms_viirs_active_fire.csv")
else:
    print(f"Failed: {response.status_code}: {response.text}")