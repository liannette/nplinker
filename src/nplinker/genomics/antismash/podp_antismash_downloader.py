from __future__ import annotations
import json
import logging
import time
import warnings
from collections.abc import Mapping
from collections.abc import Sequence
from os import PathLike
from pathlib import Path
from jsonschema import validate
from nplinker.defaults import GENOME_STATUS_FILENAME
from nplinker.genomics.antismash import GenomeAccessionResolver
from nplinker.genomics.antismash import antismash_job_is_done
from nplinker.genomics.antismash import download_and_extract_from_antismash_api
from nplinker.genomics.antismash import download_and_extract_from_antismash_db
from nplinker.genomics.antismash import download_and_extract_ncbi_genome
from nplinker.genomics.antismash import extract_antismash_data
from nplinker.genomics.antismash import initiate_antismash_job

# TODO: update GENOME_STATUS_SCHEMA
from nplinker.schemas import GENOME_STATUS_SCHEMA


logger = logging.getLogger(__name__)


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
            "resolved_genbank_id": self.resolved_genbank_id,
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

    # Ensure a Genome Status object for each genome record
    gs_file = Path(project_download_root, GENOME_STATUS_FILENAME)
    gs_dict = GenomeStatus.read_json(gs_file)
    prepare_genome_records(genome_records, gs_dict)
    GenomeStatus.to_json(gs_dict, gs_file)

    # Resolve genome assembly accessions
    resolve_assembly_accessions(genome_records, gs_dict)
    GenomeStatus.to_json(gs_dict, gs_file)

    # Look for already downloaded antismash data
    check_for_existing_antismash_data(gs_dict, project_extract_root)
    GenomeStatus.to_json(gs_dict, gs_file)

    # Check already submitted antismash jobs
    process_antismash_jobs(gs_dict, project_download_root, project_extract_root)
    GenomeStatus.to_json(gs_dict, gs_file)

    # Try to get antismash data from antiSMASH-DB
    obtain_data_from_antismash_db(gs_dict, project_download_root, project_extract_root)
    GenomeStatus.to_json(gs_dict, gs_file)

    # Submit jobs to antiSMASH API
    submit_antismash_jobs(gs_dict, project_download_root, project_extract_root, gs_file)

    # Process antiSMASH API jobs
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


def prepare_genome_records(genome_records, gs_dict):
    logger.info(f"Total genome records from the PODP dataset to process: {len(genome_records)}")
    for genome_record in genome_records:
        genome_id_data = genome_record["genome_ID"]

        # get the best available genome ID (RefSeq > GenBank > JGI)
        raw_genome_id = get_best_available_genome_id(genome_id_data)
        if not raw_genome_id:
            logger.warning(f'Invalid input genome record "{genome_record}"')
            continue

        # create the GenomeStatus Object if not already there
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
        


def resolve_assembly_accessions(genome_records, gs_dict):
    genome_ids_to_check = [
        genome_id for genome_id, gs_obj in gs_dict.items() if not gs_obj.resolve_attempted
    ]
    if not genome_ids_to_check:
        return

    logger.info(
        f"Attempting to resolve genome assembly accessions for {len(genome_ids_to_check)} genomes."
    )
    successful_count = 0
    for genome_record in genome_records:
        genome_id_data = genome_record["genome_ID"]
        raw_genome_id = get_best_available_genome_id(genome_id_data)
        if not raw_genome_id:
            continue
        gs_obj = gs_dict[raw_genome_id]

        try:
            genbank_acc, refseq_acc = GenomeAccessionResolver().resolve(genome_id_data)
            gs_obj.resolved_genbank_id = genbank_acc
            gs_obj.resolved_refseq_id = refseq_acc
            successful_count += 1
        except Exception:
            logger.warning(
                f"Failed to resolve genome assembly accessions for genome ID data: {genome_id_data}"
            )
        finally:
            gs_obj.resolve_attempted = True

    if successful_count > 0:
        logger.info(
            f"Successfully resolved genome assembly accessions for {successful_count} genomes"
        )
    unsuccessful_count = len(genome_ids_to_check) - successful_count
    if unsuccessful_count > 0:
        logger.warning(
            f"Failed to resolve genome assembly accessions for {unsuccessful_count} genomes"
        )


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


def check_for_existing_antismash_data(gs_dict, project_extract_root):
    """Extract existing antiSMASH data.

    This function iterates over a dictionary of genome objects, checks if they have
    an associated BGC path, and attempts to extract the antiSMASH data from these paths.
    If the extraction is unsucessful, the associated BGC path is removed.
    It logs the number of successful and unsuccessful extractions.

    Args:
        gs_dict (dict): A dictionary where keys are genome IDs and values are genome
            objects containing BGC paths.
        project_extract_root (str): The root directory where the extracted antiSMASH
            data should be stored.

    Returns:
        None
    """
    genome_ids_to_check = [genome_id for genome_id, gs_obj in gs_dict.items() if gs_obj.bgc_path]
    if not genome_ids_to_check:
        return

    logger.info(
        f"Extracting already existing antiSMASH results for {len(genome_ids_to_check)} genomes."
    )

    successful_count = 0

    for genome_id in genome_ids_to_check:
        gs_obj = gs_dict[genome_id]
        try:
            antimash_id = Path(gs_obj.bgc_path).stem
            extract_antismash_data(gs_obj.bgc_path, project_extract_root, antimash_id)

            output_path = Path(project_extract_root, "antismash", antimash_id)
            if output_path.exists():
                Path.touch(output_path / "completed", exist_ok=True)
                successful_count += 1

        except Exception as e:
            logger.warning(
                f"Failed to extract already downloaded antismash data for {genome_id}. Error: {e}"
            )
            gs_obj.bgc_path = ""

    if successful_count > 0:
        logger.info(f"Successfully extracted {successful_count} already existing antismash results")
    unsuccessful_count = len(genome_ids_to_check) - successful_count
    if unsuccessful_count > 0:
        logger.warning(f"Failed to extract {unsuccessful_count} already existing antismash results")


def process_existing_bgc_data(gs_dict, project_extract_root):
    """
    Processes existing antiSMASH results.

    This function checks if there are any existing antiSMASH results for the 
    genomes provided in the `gs_dict`. If such data exists, it attempts to extract 
    the antiSMASH data into the specified project extraction root directory.
    It logs the number of successful, unsuccessful, and skipped extractions.

    Args:
        gs_dict (dict): A dictionary where keys are genome IDs and values are objects containing genome data.
                        Each object should have a `bgc_path` attribute indicating the path to the antiSMASH data.
        project_extract_root (str): The root directory where the extracted antiSMASH data should be stored.

    Returns:
        None
    """

    # Check if any genomes has already downloaded data
    if all(not gs_obj.bgc_path for gs_obj in gs_dict.values()):
        return

    logger.info(
        "Extracting already existing antiSMASH results..."
    )

    successful_cnt = 0
    unsuccessful_cnt = 0
    skipped_cnt = 0

    for genome_id, gs_obj in gs_dict.items():

        # skip if no downloaded antismash data
        if not gs_obj.bgc_path:
            skipped_cnt += 1
            continue

        # try to extract antismash data
        try:
            antimash_id = Path(gs_obj.bgc_path).stem
            extract_antismash_data(gs_obj.bgc_path, project_extract_root, antimash_id)
            output_path = Path(project_extract_root, "antismash", antimash_id)
            Path.touch(output_path / "completed", exist_ok=True)
            successful_count += 1
        except Exception as e:
            logger.warning(
                f"Failed to extract antismash results for {genome_id}. Error: {e}"
            )
            gs_obj.bgc_path = ""
            unsuccessful_cnt += 1

    logger.info(
        "Summary of existing antiSMASH data extraction: "
        f"{successful_cnt} successful, {unsuccessful_cnt} failed, {skipped_cnt} skipped."
        )

def obtain_data_from_antismash_db(gs_dict, project_download_root, project_extract_root):
    """Obtain BGC data from the antiSMASH database.

    This function attempts to download and extract BGC data from the antiSMASH database
    for genomes that do not yet have BGC data paths and have a resolved RefSeq ID.
    It logs the progress and any errors encountered during the process.

    Args:
        gs_dict (dict): A dictionary where keys are genome IDs and values are genome objects.
            Each genome object should have attributes `bgc_path` and `resolved_refseq_id`.
        project_download_root (str): The root directory where downloaded data will be stored.
        project_extract_root (str): The root directory where extracted data will be stored.

    Returns:
        None
    """
    genome_ids_to_check = [
        genome_id
        for genome_id, gs_obj in gs_dict.items()
        if not gs_obj.bgc_path and gs_obj.resolved_refseq_id
    ]
    if not genome_ids_to_check:
        return

    logger.info(
        f"Attempting to obtain BGC data from antiSMASH-DB for {len(genome_ids_to_check)} genomes."
    )

    successful_count = 0

    for genome_id in genome_ids_to_check:
        gs_obj = gs_dict[genome_id]
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
        except Exception:
            continue
        # except Exception as e:
        #     logger.warning(
        #         f"Failed to obtain antiSMASH DB data for genome ID {genome_id} found. Error: {e}"
        #     )

    if successful_count > 0:
        logger.info(
            f"Successfully obtained BGC data from antiSMASH-DB for {successful_count} genomes"
        )
    unsuccessful_count = len(genome_ids_to_check) - successful_count
    if unsuccessful_count > 0:
        logger.warning(
            f"Failed to obtain BGC data from antiSMASH-DB for {unsuccessful_count} genomes"
        )


def process_antismash_jobs(gs_dict, project_download_root, project_extract_root):
    """Process submitted antiSMASH jobs and download the results.

    This function checks the status of submitted antiSMASH jobs, downloads and extracts
    the results for completed jobs, and updates the genome status dictionary accordingly.

    Args:
        gs_dict (dict): A dictionary where keys are genome IDs and values are GenomeStatus objects.
        project_download_root (str): The root directory where downloaded data will be stored.
        project_extract_root (str): The root directory where extracted data will be stored.

    Returns:
        None
    """
    # Get the genome IDs that have submitted antiSMASH jobs but no BGC path
    genome_ids_to_check = [
        genome_id
        for genome_id, gs_obj in gs_dict.items()
        if not gs_obj.bgc_path and gs_obj.antismash_job_id
    ]
    if not genome_ids_to_check:
        return

    logger.info(f"Downloading and extracting results from {len(genome_ids_to_check)} submitted antiSMASH jobs.")

    successful_count = 0
    pending_genome_ids = set(genome_ids_to_check)

    TODO: use pop()
    while pending_genome_ids:
        finished_genome_ids = set()

        for genome_id in list(pending_genome_ids):
            gs_obj = gs_dict[genome_id]

            try:
                # Check if the antiSMASH job is done
                if not antismash_job_is_done(gs_obj.antismash_job_id):
                    continue

                # Download and extract the results
                output_path = download_and_extract_from_antismash_api(
                    gs_obj.antismash_job_id,
                    project_download_root,
                    project_extract_root,
                )
                if output_path.exists():
                    Path.touch(output_path / "completed", exist_ok=True)

                antimash_id = Path(output_path).stem
                gs_obj.bgc_path = str(Path(project_download_root, antimash_id + ".zip").absolute())

                finished_genome_ids.add(genome_id)
                successful_count += 1

            except Exception as e:
                logger.warning(
                    f"Failed to obtain antiSMASH job results for genome {genome_id}. Error: {e}"
                )
                gs_obj.antismash_job_id = ""
                finished_genome_ids.add(genome_id)

        # Remove finished genome IDs from the pending set
        pending_genome_ids -= finished_genome_ids

        if pending_genome_ids:
            # At least one job is still running, wait before polling again
            time.sleep(30)

    if successful_count > 0:
        logger.info(f"Successfully downloaded and extracted {successful_count} antiSMASH jobs.")
    unsuccessful_count = len(genome_ids_to_check) - successful_count
    if unsuccessful_count > 0:
        logger.warning(f"Failed to download and extract {unsuccessful_count} antiSMASH jobs.")


def submit_antismash_jobs(gs_dict, project_download_root, project_extract_root, gs_file):
    genome_ids_to_check = [
        genome_id
        for genome_id, gs_obj in gs_dict.items()
        if not gs_obj.bgc_path and (gs_obj.resolved_refseq_id or gs_obj.resolved_genbank_id)
    ]
    if not genome_ids_to_check:
        return

    logger.info(
        f"Submitting {len(genome_ids_to_check)} jobs to antiSMASH API with a maximum of 5 "
        "concurrent jobs to avoid overloading the service."
    )
    
    # TODO: use list pop here
    successful_count = 0
    running_jobs = set()
    for genome_id in genome_ids_to_check:
        # Ensure that max 5 jobs are running in parallel for karma points
        while len(running_jobs) >= 5:
            completed_jobs = set()
            for job_id in running_jobs:
                try:
                    if antismash_job_is_done(job_id):
                        completed_jobs.add(job_id)
                except Exception as e:
                    logger.warning(
                        f"Error occurred while checking the status of antiSMASH job {job_id}: {e}"
                    )
                    completed_jobs.add(job_id)

            if completed_jobs:
                running_jobs = running_jobs - completed_jobs
            else:
                # Wait before polling again
                time.sleep(30)

        gs_obj = gs_dict[genome_id]
        if gs_obj.antismash_job_id != "":
            raise ValueError("AntiSMASH job ID is not empty")

        try:
            # First get genome from ncbi
            genome_acc = gs_obj.resolved_refseq_id or gs_obj.resolved_genbank_id
            genbank_file_path = download_and_extract_ncbi_genome(
                genome_acc, project_download_root, project_extract_root
            )
            # then submit genome it to antismash API
            job_id = initiate_antismash_job(genbank_file_path)

            if not job_id:
                raise RuntimeError("No job ID returned.")
            gs_obj.antismash_job_id = job_id
            running_jobs.add(job_id)
            successful_count += 1

            # Update GenomeStatus file to ensure no job IDs are lost in case of an interruption
            GenomeStatus.to_json(gs_dict, gs_file)

        except Exception as e:
            logger.warning(f"Failed to submit antismash job for genome ID {genome_id}. Error: {e}")

    if successful_count > 0:
        logger.info(f"Successfully submitted antiSMASH jobs for {successful_count} genomes")
    unsuccessful_count = len(genome_ids_to_check) - successful_count
    if unsuccessful_count > 0:
        logger.warning(f"Failed to submit antiSMASH jobs for {unsuccessful_count} genomes")
