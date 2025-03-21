from __future__ import annotations
import os
import shutil
from os import PathLike
from pathlib import Path
import requests
from bs4 import BeautifulSoup
from nplinker.utils import download_and_extract_archive
from nplinker.utils import extract_archive
from nplinker.utils import list_dirs
from nplinker.utils import list_files


# urls to be given to download antismash DB data
ANTISMASH_DB_DOWNLOAD_URL = "https://antismash-db.secondarymetabolites.org/output/{}/{}"
# The antiSMASH DBV2 is for the availability of the old version, better to keep it.
ANTISMASH_DBV2_DOWNLOAD_URL = "https://antismash-dbv2.secondarymetabolites.org/output/{}/{}"

# url to download antismash API data
ANTISMASH_API_DOWNLOAD_URL = "https://antismash.secondarymetabolites.org/upload/{}/{}"


def download_and_extract_from_antismash_api(
    job_id: str, download_root: str | PathLike, extract_root: str | PathLike
) -> None:
    """Downloads and extracts results from an antiSMASH API job for a given job ID.

    This function constructs the download URL using the provided job ID then
    downloads the results as a ZIP file and extracts its contents to the specified directories.

    Args:
        job_id (str): The job ID for the antiSMASH API job.
        download_root (str or PathLike): The root directory where the ZIP file will be downloaded.
        extract_root (str or PathLike): The root directory where the contents of the ZIP file will be extracted.

    Raises:
        requests.exceptions.RequestException: If there is an issue with the HTTP request.
        zipfile.BadZipFile: If the downloaded file is not a valid ZIP file.
        OSError: If there is an issue with file operations such as writing or extracting.
    """
    antismash_id = _get_antismash_id(job_id)
    url = ANTISMASH_API_DOWNLOAD_URL.format(job_id, antismash_id + ".zip")
    extract_path = download_and_extract_antismash_data(
        url, antismash_id, download_root, extract_root
    )
    return extract_path


def download_and_extract_from_antismash_db(
    refseq_acc: str, download_root: str | PathLike, extract_root: str | PathLike
) -> None:
    """Download and extract antiSMASH BGC archive for a specified genome.

    The antiSMASH database (https://antismash-db.secondarymetabolites.org/)
    is used to download the BGC archive. And antiSMASH use RefSeq assembly id
    of a genome as the id of the archive.

    Args:
        refseq_acc: The id used to download BGC archive from antiSMASH database.
            If the id is versioned (e.g., "GCF_004339725.1") please be sure to
            specify the version as well.
        download_root: Path to the directory to place downloaded archive in.
        extract_root: Path to the directory data files will be extracted to.
            Note that an `antismash` directory will be created in the specified `extract_root` if
            it doesn't exist. The files will be extracted to `<extract_root>/antismash/<antismash_id>` directory.

    Raises:
        ValueError: if `<extract_root>/antismash/<refseq_assembly_id>` dir is not empty.

    Examples:
        >>> download_and_extract_antismash_metadata("GCF_004339725.1", "/data/download", "/data/extracted")
    """
    for base_url in [ANTISMASH_DB_DOWNLOAD_URL, ANTISMASH_DBV2_DOWNLOAD_URL]:
        url = base_url.format(refseq_acc, refseq_acc + ".zip")
        if requests.head(url).status_code == 404:  # not found
            continue
        extract_path = download_and_extract_antismash_data(
            url, refseq_acc, download_root, extract_root
        )
        return extract_path
    raise ValueError(f"No results in antiSMASH DB for {refseq_acc}")


def download_and_extract_antismash_data(
    url: str, antimash_id: str, download_root: str | PathLike, extract_root: str | PathLike
) -> str | PathLike:
    download_root = Path(download_root)
    extract_path = Path(extract_root) / "antismash" / antimash_id
    try:
        _prepare_antismash_extract_path(extract_path)
        download_and_extract_archive(url, download_root, extract_path, antimash_id + ".zip")
        _cleanup_extracted_files(extract_path)
    except Exception as e:
        shutil.rmtree(extract_path)
        raise e
    return extract_path


def extract_antismash_data(
    archive: str, extract_root: str | PathLike, antimash_id: str
) -> str | PathLike:
    extract_path = Path(extract_root) / "antismash" / antimash_id
    try:
        _prepare_antismash_extract_path(extract_path)
        extract_archive(archive, extract_path, remove_finished=False)
        _cleanup_extracted_files(extract_path)

    except Exception as e:
        shutil.rmtree(extract_path)
        raise e


def _prepare_antismash_extract_path(extract_path: str | PathLike) -> None:
    # check if extract_path is empty
    if extract_path.exists() and any(extract_path.iterdir()):
        raise ValueError(f'Nonempty directory: "{extract_path}"')
    else:
        extract_path.mkdir(parents=True, exist_ok=True)


def _cleanup_extracted_files(extract_path: Path) -> None:
    # delete subdirs
    for subdir_path in list_dirs(extract_path):
        shutil.rmtree(subdir_path)

    # delete unnecessary files
    files_to_keep = list_files(extract_path, suffix=(".json", ".gbk"))
    for file in list_files(extract_path):
        if file not in files_to_keep:
            os.remove(file)


def _get_antismash_id(job_id: str):
    """Fetches the antiSMASH ID for a given job ID by scraping the antiSMASH website.

    Args:
        job_id (str): The antismash job ID.

    Returns:
        str: The antiSMASH ID if found, otherwise None.

    Raises:
        requests.RequestException: If there is an error while making the HTTP request.
    """
    url = f"https://antismash.secondarymetabolites.org/upload/{job_id}/index.html"
    response = requests.get(url)
    response.raise_for_status()

    soup = BeautifulSoup(response.text, "html.parser")
    link = soup.find("a", href=lambda href: href and href.endswith(".zip"))
    antismash_id = link["href"].rstrip(".zip")
    return antismash_id
