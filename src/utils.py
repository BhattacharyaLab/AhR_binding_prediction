import requests
import shutil


# def download_file(url, dest):
#     with requests.get(url, stream=True) as r:
#         with open(dest, 'wb') as f:
#             shutil.copyfileobj(r.raw, f, length=16*1024*1024)

#     return dest


def download_file(url, dest):
    # keep track of how much of the file is downloaded
    total = 0
    last_percent = 0
    
    # get the size of the file
    info = requests.head(url)
    file_size = int(info.headers['Content-Length'])
    
    # download the file
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(dest, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024*1024):
                total += len(chunk)
                if (total*100 // (file_size)) > last_percent:
                    last_percent = total*100 // (file_size)
                    print(f'Downloaded {last_percent}%')
                f.write(chunk)
    return dest