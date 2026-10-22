"""
download_ame2020_data.py
Downloads official AME2020 mass evaluation and NUBASE2020 files from IAEA NDS.
"""

import os
import sys
import urllib.request

FILES_TO_DOWNLOAD = {
    "mass_1.mas20.txt": "https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt",
    "nubase_4.mas20.txt": "https://www-nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt",
}


def download_file(url: str, output_path: str, chunk_size: int = 65536) -> None:
    """Streams a file download with size tracking and a text-based progress bar."""
    print(f"Connecting to: {url}")
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
    req = urllib.request.Request(url, headers=headers)

    with urllib.request.urlopen(req) as response:
        total_size = response.getheader("Content-Length")
        total_size = int(total_size) if total_size else None

        bytes_read = 0
        with open(output_path, "wb") as f:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                bytes_read += len(chunk)

                if total_size:
                    percent = (bytes_read / total_size) * 100
                    mb_read = bytes_read / (1024 * 1024)
                    mb_total = total_size / (1024 * 1024)
                    sys.stdout.write(
                        f"\rDownloading {os.path.basename(output_path)}: "
                        f"[{percent:5.1f}%] {mb_read:.2f} / {mb_total:.2f} MB"
                    )
                else:
                    mb_read = bytes_read / (1024 * 1024)
                    sys.stdout.write(
                        f"\rDownloading {os.path.basename(output_path)}: {mb_read:.2f} MB read"
                    )
                sys.stdout.flush()

    print(f"\n[✓] Saved to: {os.path.abspath(output_path)}\n")


def main(target_dir: str = "./ame2020_data") -> None:
    os.makedirs(target_dir, exist_ok=True)

    for filename, url in FILES_TO_DOWNLOAD.items():
        destination = os.path.join(target_dir, filename)
        if os.path.exists(destination):
            print(f"File already exists: {destination} (skipping)")
            continue
        try:
            download_file(url, destination)
        except Exception as exc:
            print(f"\n[✗] Failed to download {filename}: {exc}")


if __name__ == "__main__":
    main()
