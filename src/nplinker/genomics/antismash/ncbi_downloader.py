import logging
import os
import shutil
from os import PathLike
from pathlib import Path
from typing import Optional
from nplinker.utils import check_md5
from nplinker.utils import download_url
from nplinker.utils import extract_archive


logger = logging.getLogger(__name__)


# TODO: Fix this function

def download_and_extract_ncbi_genome(
    refseq_id: str,
    download_root: str | PathLike,
    extract_root: str | PathLike,
    max_retries: int = 10,
) -> Optional[Path]:
    """Downloads and extracts an NCBI dataset for a given genome refseq ID.

    This function attempts to download a dataset from the NCBI database using
    the provided refseq ID. It retries the download and extraction process up
    to a maximum number of times if any errors occur. The function verifies
    the integrity of the downloaded files using MD5 checksums and moves and
    renames the GenBank files upon successful verification.

    Args:
        refseq_id (str): The refseq ID for the dataset to be downloaded.
        download_root (str or Path): The root directory where the dataset will be downloaded.
        extract_root (str or Path): The root directory where the dataset will be extracted.
        max_retries (int): The maximum number of times to retry downloading and extracting

    Returns:
        Path: The path to the extracted dataset if successful, otherwise None.

    Raises:
        Exception: If the maximum number of retries is reached and the dataset could
        not be successfully downloaded and extracted.
    """
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{refseq_id}/download?include_annotation_type=GENOME_GB"
    )

    download_root = Path(download_root)
    extract_path = Path(extract_root) / "ncbi_genomes"
    filename = f"ncbi_{refseq_id}.zip"

    extract_path.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, max_retries + 1):
        try:
            # logger.info(
            #     f"Attempt {attempt}: Downloading and extracting genome assembly"
            #     f"with RefSeq accession {refseq_id} from NCBI ..."
            # )

            download_url(url, download_root, filename)
            archive = download_root / filename
            break
        except Exception as e:
            logger.warning(f"Error occurred during attempt {attempt}: {e}")

        extract_archive(archive, extract_path)

        md5_ok = verify_ncbi_dataset_md5_sums(extract_path)
        if not md5_ok:
            continue

        # Move and rename GenBank file
        genbank_path = extract_path / "ncbi_dataset" / "data" / refseq_id / "genomic.gbff"
        new_genbank_path = extract_path / f"{refseq_id}.gbff"
        genbank_path.rename(new_genbank_path)

        # Delete unnecessary files
        shutil.rmtree(extract_path / "ncbi_dataset")
        os.remove(extract_path / "md5sum.txt")
        os.remove(extract_path / "README.md")

        return new_genbank_path

    logger.error("Maximum retries reached. Download and extraction failed.")
    return None


def verify_ncbi_dataset_md5_sums(extract_path: PathLike) -> bool:
    """Check MD5 checksums for files in the extraction path.

    Args:
        extract_path (Path): Path to the extraction directory.

    Returns:
        bool: True if all MD5 checksums are correct, False otherwise.
    """
    md5_ok = True
    with open(extract_path / "md5sum.txt", "r") as f:
        for line in f:
            md5sum, file_name = line.strip().split()
            file_path = extract_path / file_name
            if not check_md5(file_path, md5sum):
                logger.error(f"MD5 checksum mismatch for {file_path}")
                md5_ok = False
                break
    return md5_ok
