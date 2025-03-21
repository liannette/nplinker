from __future__ import annotations
import json
import logging
import re
import time
import warnings
from collections.abc import Mapping
from collections.abc import Sequence
from os import PathLike
from pathlib import Path
import httpx
from bs4 import BeautifulSoup
from jsonschema import validate
from nplinker.defaults import GENOME_STATUS_FILENAME
from nplinker.genomics.antismash import antismash_job_is_done
from nplinker.genomics.antismash import download_and_extract_from_antismash_api
from nplinker.genomics.antismash import download_and_extract_from_antismash_db
from nplinker.genomics.antismash import download_and_extract_ncbi_genome
from nplinker.genomics.antismash import extract_antismash_data
from nplinker.genomics.antismash import initiate_antismash_job

# TODO: update GENOME_STATUS_SCHEMA
from nplinker.schemas import GENOME_STATUS_SCHEMA


logger = logging.getLogger(__name__)

JGI_GENOME_LOOKUP_URL = (
    "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid={}"
)
USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:86.0) Gecko/20100101 Firefox/86.0"


class GenomeStatus:
    """Class to represent the status of a single genome.

    The status of genomes is tracked in the file
    [GENOME_STATUS_FILENAME][nplinker.defaults.GENOME_STATUS_FILENAME].
    """

    def __init__(
        self,
        original_id: str,
        resolved_genbank_id: str = "",
        resolved_refseq_id: str = "",
        resolve_attempted: bool = False,
        antismash_job_id: str = "",
        bgc_path: str = "",
    ):
        """Initialize a GenomeStatus object for the given genome.

        Args:
            original_id: The original ID of the genome.
            resolved_genbank_id: The resolved GenBank accession of the genome assembly.
                Defaults to "".
            resolved_refseq_id: The resolved RefSeq accession of the genome assembly.
                Defaults to "".
            resolve_attempted: A flag indicating whether an
                attempt to resolve the RefSeq ID has been made. Defaults to False.
            antismash_job_id: The ID of the antiSMASH job for the genome
                if a job has been submitted to the antiSMASH API. Defaults to "".
            bgc_path: The path to the downloaded BGC file for
                the genome. Defaults to "".
        """
        self.original_id = original_id
        self.resolved_genbank_id = "" if resolved_genbank_id == "None" else resolved_genbank_id
        self.resolved_refseq_id = "" if resolved_refseq_id == "None" else resolved_refseq_id
        self.resolve_attempted = resolve_attempted
        self.antismash_job_id = antismash_job_id
        self.bgc_path = bgc_path

    @staticmethod
    def read_json(file: str | PathLike) -> dict[str, "GenomeStatus"]:
        """Get a dict of GenomeStatus objects by loading given genome status file.

        Note that an empty dict is returned if the given file doesn't exist.

        Args:
            file: Path to genome status file.

        Returns:
            Dict keys are genome original id and values are GenomeStatus
                objects. An empty dict is returned if the given file doesn't exist.
        """
        genome_status_dict = {}
        if Path(file).exists():
            with open(file, "r") as f:
                data = json.load(f)

            # validate json data before using it
            validate(data, schema=GENOME_STATUS_SCHEMA)

            genome_status_dict = {
                gs["original_id"]: GenomeStatus(**gs) for gs in data["genome_status"]
            }
        return genome_status_dict

    @staticmethod
    def to_json(
        genome_status_dict: Mapping[str, "GenomeStatus"], file: str | PathLike | None = None
    ) -> str | None:
        """Convert the genome status dictionary to a JSON string.

        If a file path is provided, the JSON string is written to the file. If
        the file already exists, it is overwritten.

        Args:
            genome_status_dict: A dictionary of genome
                status objects. The keys are the original genome IDs and the values
                are GenomeStatus objects.
            file: The path to the output JSON file.
                If None, the JSON string is returned but not written to a file.

        Returns:
            The JSON string if `file` is None, otherwise None.
        """
        gs_list = [gs._to_dict() for gs in genome_status_dict.values()]
        json_data = {"genome_status": gs_list, "version": "1.0"}

        # validate json object before dumping
        validate(json_data, schema=GENOME_STATUS_SCHEMA)

        if file is not None:
            with open(file, "w") as f:
                json.dump(json_data, f)
            return None
        return json.dumps(json_data)

    def _to_dict(self) -> dict:
        """Convert the GenomeStatus object to a dict."""
        return {
            "original_id": self.original_id,
            "resolved_genbank_id": self.resolved_refseq_id,
            "resolved_refseq_id": self.resolved_refseq_id,
            "resolve_attempted": self.resolve_attempted,
            "antismash_job_id": self.antismash_job_id,
            "bgc_path": self.bgc_path,
        }


def podp_download_and_extract_antismash_data(
    genome_records: Sequence[Mapping[str, Mapping[str, str]]],
    project_download_root: str | PathLike,
    project_extract_root: str | PathLike,
):
    """Download and extract antiSMASH BGC archive for the given genome records.

    Args:
        genome_records: list of dicts representing genome records.

            The dict of each genome record contains a key of genome ID with a value
            of another dict containing information about genome type, label and
            accession ids (RefSeq, GenBank, and/or JGI).
        project_download_root: Path to the directory to place
            downloaded archive in.
        project_extract_root: Path to the directory downloaded archive will be extracted to.

            Note that an `antismash` directory will be created in the specified
            `extract_root` if it doesn't exist. The files will be extracted to
            `<extract_root>/antismash/<antismash_id>` directory.

    Warnings:
        UserWarning: when no antiSMASH data is found for some genomes.
    """
    if not Path(project_download_root).exists():
        # otherwise in case of failed first download, the folder doesn't exist and
        # genome_status_file can't be written
        Path(project_download_root).mkdir(parents=True, exist_ok=True)

    # Load the Genome Status from the json
    gs_file = Path(project_download_root, GENOME_STATUS_FILENAME)
    gs_dict = GenomeStatus.read_json(gs_file)

    # Make sure that every genome record has a genome status in the gs_dict
    logger.info(f"Checking {len(genome_records)} genome records from the PODP dataset")
    prepare_genome_status_for_all(genome_records, gs_dict)

    # Resolve for refseq_accessions
    add_refseq_accessions(genome_records, gs_dict)

    # Look for already downloaded antismash data
    logger.info("Checking for already downloaded antismash data")
    check_for_existing_antismash_data(gs_dict, project_extract_root)
    GenomeStatus.to_json(gs_dict, gs_file)

    # Check already submitted antismash jobs
    logger.info("Checking already submitted jobs to antiSMASH API")
    process_antismash_jobs(gs_dict, project_download_root, project_extract_root)
    GenomeStatus.to_json(gs_dict, gs_file)

    # Try to get antismash data from antiSMASH-DB
    logger.info("Checking for BGC data in antismash-DB")
    obtain_data_from_antismash_db(gs_dict, project_download_root, project_extract_root)
    GenomeStatus.to_json(gs_dict, gs_file)

    # Submit jobs to antiSMASH API
    logger.info(
        "Submitting jobs to antiSMASH API (maximum 5 concurrent jobs to avoid overloading the service)"
    )
    submit_antismash_jobs(gs_dict, project_download_root, project_extract_root)
    GenomeStatus.to_json(gs_dict, gs_file)

    # Process antiSMASH API jobs
    logger.info("Processing results from submitted antiSMASH API jobs")
    process_antismash_jobs(gs_dict, project_download_root, project_extract_root)
    GenomeStatus.to_json(gs_dict, gs_file)

    # raise and log warning for genomes with no antismash data
    failed_ids = [gs.original_id for gs in gs_dict.values() if not gs.bgc_path]
    if failed_ids:
        warning_message = f"No antiSMASH data for the following genome IDs: {failed_ids}"
        logger.warning(warning_message)
        warnings.warn(warning_message, UserWarning)

    # Raise error if no antismash data available at all
    if len(failed_ids) == len(genome_records):
        raise ValueError("No antiSMASH data found for any genome")


def get_best_available_genome_id(genome_id_data: Mapping[str, str]) -> str | None:
    """Get the best available ID from genome_id_data dict.

    Args:
        genome_id_data: dictionary containing information for each genome record present.

    Returns:
        ID for the genome, if present, otherwise None.
    """
    if "RefSeq_accession" in genome_id_data:
        best_id = genome_id_data["RefSeq_accession"]
    elif "GenBank_accession" in genome_id_data:
        best_id = genome_id_data["GenBank_accession"]
    elif "JGI_Genome_ID" in genome_id_data:
        best_id = genome_id_data["JGI_Genome_ID"]
    else:
        best_id = None

    if best_id is None or len(best_id) == 0:
        logger.warning(f"Failed to get valid genome ID in genome data: {genome_id_data}")
        return None
    return best_id


def _resolve_genbank_accession(genbank_id: str) -> str:
    """Try to get RefSeq assembly id through given GenBank assembly id.

    Note that GenBank assembly accession starts with "GCA_" and RefSeq assembly
    accession starts with "GCF_". For more info, see
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/troubleshooting/faq

    Args:
        genbank_id: ID for GenBank assembly accession.

    Raises:
        httpx.ReadTimeout: If the request times out.

    Returns:
        RefSeq assembly ID if the search is successful, otherwise an empty string.
    """
    logger.info(
        f"Attempting to resolve Genbank assembly accession {genbank_id} to RefSeq accession"
    )
    # NCBI Datasets API https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/
    # Note that there is a API rate limit of 5 requests per second without using an API key
    # For more info, see https://www.ncbi.nlm.nih.gov/datasets/docs/v2/troubleshooting/faq/

    # API for getting revision history of a genome assembly
    # For schema, see https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#get-/genome/accession/-accession-/revision_history
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{genbank_id}/revision_history"

    refseq_id = ""
    try:
        resp = httpx.get(
            url, headers={"User-Agent": USER_AGENT}, timeout=10.0, follow_redirects=True
        )
        resp.raise_for_status()

        data = resp.json()
        if not data:
            raise ValueError("No Assembly Revision data found")

        assembly_entries = [
            entry for entry in data["assembly_revisions"] if "refseq_accession" in entry
        ]
        if not assembly_entries:
            raise ValueError("No RefSeq assembly accession found")

        latest_entry = max(assembly_entries, key=lambda x: x["release_date"])
        refseq_id = latest_entry["refseq_accession"]

    except httpx.RequestError as exc:
        logger.warning(f"An error occurred while requesting {exc.request.url!r}: {exc}")
    except httpx.HTTPStatusError as exc:
        logger.warning(
            f"Error response {exc.response.status_code} while requesting {exc.request.url!r}"
        )
    except httpx.ReadTimeout:
        logger.warning("Timed out waiting for result of GenBank assembly lookup")
    except ValueError as exc:
        logger.warning(f"Error while resolving GenBank assembly accession {genbank_id}: {exc}")

    return refseq_id


def _resolve_jgi_accession(jgi_id: str) -> str:
    """Try to get RefSeq id through given JGI id.

    Args:
        jgi_id: JGI_Genome_ID for GenBank accession.

    Returns:
        RefSeq ID if search is successful, otherwise an empty string.
    """
    url = JGI_GENOME_LOOKUP_URL.format(jgi_id)
    logger.info(f"Attempting to resolve JGI_Genome_ID {jgi_id} to GenBank accession via {url}")
    # no User-Agent header produces a 403 Forbidden error on this site...
    try:
        resp = httpx.get(
            url, headers={"User-Agent": USER_AGENT}, timeout=10.0, follow_redirects=True
        )
    except httpx.ReadTimeout:
        logger.warning("Timed out waiting for result of JGI_Genome_ID lookup")
        return ""

    soup = BeautifulSoup(resp.content, "html.parser")
    # Find the table entry giving the "NCBI Assembly Accession" ID
    link = soup.find("a", href=re.compile("https://www.ncbi.nlm.nih.gov/datasets/genome/.*"))
    if link is None:
        return ""

    assembly_id = link.text
    # check if the assembly ID is already a RefSeq ID
    if assembly_id.startswith("GCF_"):
        return assembly_id  # type: ignore
    else:
        return _resolve_genbank_accession(assembly_id)


def _resolve_refseq_id(genome_id_data: Mapping[str, str]) -> str:
    """Get the RefSeq ID to which the genome accession is linked.

    Check https://pairedomicsdata.bioinformatics.nl/schema.json.

    Args:
        genome_id_data: dictionary containing information
        for each genome record present.

    Returns:
        RefSeq ID if present, otherwise an empty string.
    """
    if "RefSeq_accession" in genome_id_data:
        # best case, can use this directly
        return genome_id_data["RefSeq_accession"]
    if "GenBank_accession" in genome_id_data:
        # resolve via NCBI
        return _resolve_genbank_accession(genome_id_data["GenBank_accession"])
    if "JGI_Genome_ID" in genome_id_data:
        # resolve via JGI => NCBI
        return _resolve_jgi_accession(genome_id_data["JGI_Genome_ID"])

    logger.warning(f"Unable to resolve genome_ID: {genome_id_data}")
    return ""


def prepare_genome_status_for_all(genome_records, gs_dict):
    # Make sure that every genome record has a Genome Status
    for genome_record in genome_records:
        # get the best available ID from the dict
        genome_id_data = genome_record["genome_ID"]
        raw_genome_id = get_best_available_genome_id(genome_id_data)
        if not raw_genome_id:
            logger.warning(f'Invalid input genome record "{genome_record}"')
            continue

        if raw_genome_id not in gs_dict:
            gs_dict[raw_genome_id] = GenomeStatus(raw_genome_id)

    if len(genome_records) > len(gs_dict):
        logger.warning(
            "More genome records than genome status entries: "
            f"{len(genome_records) - len(gs_dict)} extra records "
        )
    if len(gs_dict) > len(genome_records):
        logger.warning(
            "More genome status entries than genome records: "
            f"{len(gs_dict) - len(genome_records)} extra entries"
        )


def add_refseq_accessions(genome_records, gs_dict):
    successful_count = 0
    unsuccessful_count = 0
    for genome_record in genome_records:
        # get genome id data and genome status object
        genome_id_data = genome_record["genome_ID"]
        raw_genome_id = get_best_available_genome_id(genome_id_data)
        if not raw_genome_id:
            continue
        gs_obj = gs_dict[raw_genome_id]

        # Do not resolve if not necessary
        if gs_obj.bgc_path:
            continue
        if gs_obj.resolved_refseq_id:
            continue
        if gs_obj.resolve_attempted:
            continue

        # Resolve for refseq assembly accession
        gs_obj.resolved_refseq_id = _resolve_refseq_id(genome_id_data)
        gs_obj.resolve_attempted = True

        # count success
        if gs_obj.resolved_refseq_id:
            successful_count += 1
        else:
            unsuccessful_count += 1

    if successful_count > 0:
        logger.info(
            f"Successfully resolved RefSeq genome assembly accessions for {successful_count} genome records"
        )
    if unsuccessful_count > 0:
        logger.warning(
            f"Failed to resolve RefSeq genome assembly accessions for {unsuccessful_count} genome records"
        )


def check_for_existing_antismash_data(gs_dict, project_extract_root):
    successful_count = 0
    for raw_genome_id, gs_obj in gs_dict.items():
        if not gs_obj.bgc_path:
            continue
        try:
            antimash_id = Path(gs_obj.bgc_path).stem
            extract_antismash_data(gs_obj.bgc_path, project_extract_root, antimash_id)
            output_path = Path(project_extract_root, "antismash", antimash_id)
            if output_path.exists():
                Path.touch(output_path / "completed", exist_ok=True)
            successful_count += 1
        except Exception as e:
            gs_obj.bgc_path = ""
            logger.warning(
                f"Failed to extract already downloaded antismash data for {raw_genome_id}. Error: {e}"
            )
    logger.info(f"Extracted already present antismash results for {successful_count} genomes")


def obtain_data_from_antismash_db(gs_dict, project_download_root, project_extract_root):
    successful_count = 0
    for raw_genome_id, gs_obj in gs_dict.items():
        if gs_obj.bgc_path:
            continue
        if not gs_obj.resolved_refseq_id:
            continue

        # if resolved id is valid, try to download and extract data from antismash DB
        # logger.info(f"Downloading and extracting BGC data of {raw_genome_id} from antiSMASH-DB.")
        try:
            # TODO: check if in antismash DB first before attempting to download and extract
            output_path = download_and_extract_from_antismash_db(
                gs_obj.resolved_refseq_id, project_download_root, project_extract_root
            )
            if output_path.exists():
                Path.touch(output_path / "completed", exist_ok=True)
            gs_obj.bgc_path = str(
                Path(project_download_root, gs_obj.resolved_refseq_id + ".zip").absolute()
            )
            successful_count += 1
        except Exception as e:
            logger.warning(
                f"Failed to obtain antiSMASH DB data for genome ID {raw_genome_id} found. Error: {e}"
            )
    logger.info(f"Successfully obtained BGC data from antiSMASH-DB for {successful_count} genomes")


def submit_antismash_jobs(gs_dict, project_download_root, project_extract_root):
    successful_count = 0
    running_jobs = {}
    for raw_genome_id, gs_obj in gs_dict.items():
        if gs_obj.bgc_path:
            continue
        if not gs_obj.resolved_refseq_id:
            continue
        if gs_obj.antismash_job_id != "":
            raise ValueError("AntiSMASH job ID is not empty")

        # Ensure that max 5 jobs are running in parallel for karma points
        while len(running_jobs) > 5:
            completed_jobs = set()
            for job_id in running_jobs:
                try:
                    if antismash_job_is_done(job_id):
                        completed_jobs.add(job_id)
                except Exception as e:
                    logger.warning(f"Error checking antismash job {job_id}: {e}")
                    completed_jobs.add(job_id)
                    gs_obj.antismash_job_id = ""

            if completed_jobs:
                running_jobs = running_jobs - completed_jobs
            else:
                time.sleep(10)  # Wait before polling again

        # Download and extract genome from ncbi, then submit to antismash API
        try:
            genbank_file_path = download_and_extract_ncbi_genome(
                gs_obj.resolved_refseq_id, project_download_root, project_extract_root
            )
            # logger.info(f"Genome assembly {raw_genome_id} downloaded to {genbank_file_path}")
            job_id = initiate_antismash_job(genbank_file_path)
            if not job_id:
                raise ValueError("No antiSMASH job ID")
            # logger.info(f"Successfully submitted antismash job for genome ID {raw_genome_id}")
            gs_obj.antismash_job_id = job_id
            running_jobs.add(job_id)
            successful_count += 1
        except Exception as e:
            logger.warning(
                f"Failed to submit antismash job for genome ID {raw_genome_id}. Error: {e}"
            )
    logger.info(f"Successfully submitted antiSMASH API jobs for {successful_count} genomes")


def process_antismash_jobs(gs_dict, project_download_root, project_extract_root):
    # Process jobs from the antismash API
    successful_count = 0
    unsuccessful_count = 0

    genomes_with_unprocessed_jobs = {
        raw_genome_id
        for raw_genome_id, gs_obj in gs_dict.items()
        if not gs_obj.bgc_path and gs_obj.antismash_job_id
    }

    while genomes_with_unprocessed_jobs:
        processed = set()

        for raw_genome_id in list(genomes_with_unprocessed_jobs):
            gs_obj = gs_dict[raw_genome_id]

            try:
                if not antismash_job_is_done(gs_obj.antismash_job_id):
                    continue

                output_path = download_and_extract_from_antismash_api(
                    gs_obj.antismash_job_id,
                    gs_obj.resolved_refseq_id,
                    project_download_root,
                    project_extract_root,
                )
                if output_path.exists():
                    Path.touch(output_path / "completed", exist_ok=True)
                gs_obj.bgc_path = str(
                    Path(project_download_root, gs_obj.resolved_refseq_id + ".zip").absolute()
                )
                processed.add(raw_genome_id)
                successful_count += 1

            except Exception as e:
                logger.warning(
                    f"Failed to obtain antismash API results for genome ID {raw_genome_id}. Error: {e}"
                )
                processed.add(raw_genome_id)
                gs_obj.antismash_job_id = ""
                unsuccessful_count += 1

        if processed:
            genomes_with_unprocessed_jobs = genomes_with_unprocessed_jobs - processed
        if genomes_with_unprocessed_jobs:
            time.sleep(30)
    logger.info(f"Obtained antismash results of submitted jobs for {successful_count} genomes")
