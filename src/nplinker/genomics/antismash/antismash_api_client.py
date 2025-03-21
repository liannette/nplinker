from __future__ import annotations
import logging
from os import PathLike
from pathlib import Path
from typing import Optional
import requests


logger = logging.getLogger(__name__)


def initiate_antismash_job(genbank_filepath: str | PathLike) -> Optional[str]:
    """Submits an antiSMASH job using the provided GenBank file.

    This function submits a job to the antiSMASH API and
    returns the job ID if the submission is successful. If an HTTP error or any
    other exception occurs during the submission, it prints an error message and
    returns None.

    Args:
        genbank_filepath (Path): The path to the GenBank file to be submitted to the antiSMASH API.

    Returns:
        str: The job ID if successful, otherwise None.
    """
    url = "https://antismash.secondarymetabolites.org/api/v1.0/submit"

    genbank_filepath = Path(genbank_filepath)

    try:
        with open(genbank_filepath, "rb") as file:
            files = {"seq": file}
            data = {
                "knownclusterblast": "true",
                "cc_mibig": "true",
            }
            response = requests.post(url, files=files, data=data)
            response.raise_for_status()  # Raise an exception for HTTP errors

            return response.json().get("id")

    except requests.exceptions.RequestException as req_err:
        logger.error(f"Request failed: {req_err}")
    except ValueError as json_err:  # Handles JSON decoding errors
        logger.error(f"Invalid JSON response: {json_err}")
    except Exception as err:
        logger.error(f"Unexpected error: {err}")


def query_antismash_job(job_id: str) -> Optional[dict]:
    """Gets the status of an antiSMASH job.

    Args:
        job_id (str): The job ID to query.

    Returns:
        dict: The response JSON if successful, otherwise None.
    """
    url = f"https://antismash.secondarymetabolites.org/api/v1.0/status/{job_id}"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for HTTP errors
        return response.json()

    except requests.exceptions.RequestException as req_err:
        logger.error(f"Request failed for job_id {job_id}: {req_err}")
    except ValueError as json_err:  # Handles JSON decoding errors
        logger.error(f"Invalid JSON response for job_id {job_id}: {json_err}")
    except Exception as err:
        logger.error(f"Unexpected error while getting job state for job_id {job_id}: {err}")

    return None


def antismash_job_is_done(job_id: str) -> bool:
    """Checks if the antiSMASH job is complete by polling the job status.

    Args:
        job_id (str): The job ID to query.

    Returns:
        bool: True if the job is done, False if the job is still running.

    Raises:
        RuntimeError: If the job status could not be retrieved or if the job failed.
        ValueError: If the job state is missing or unexpected in the response.
    """
    response = query_antismash_job(job_id)
    if response is None:
        raise RuntimeError(f"Failed to retrieve job status for job_id {job_id}")

    job_state = response.get("state")
    if not job_state:
        raise ValueError(f"Job state missing in response for job_id: {job_id}")
    elif job_state == "failed":
        job_status = response.get("status", "No error message provided")
        raise RuntimeError(f"AntiSMASH job {job_id} failed with an error: {job_status}")
    elif job_state in ("running", "queued"):
        return False
    elif job_state == "done":
        return True
    else:
        raise ValueError(
            f"Unexpected job state for antismash job ID {job_id}. Job state: {job_state}"
        )
