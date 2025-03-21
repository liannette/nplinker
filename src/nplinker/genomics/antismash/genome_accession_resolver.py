import logging
import re
from typing import Dict
from typing import Literal
from typing import Tuple
import httpx
from bs4 import BeautifulSoup


JGI_GENOME_LOOKUP_URL = (
    "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid={}"
)
USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:86.0) Gecko/20100101 Firefox/86.0"
logger = logging.getLogger(__name__)


class GenomeAccessionResolver:
    """Resolves genome accessions from various identifier sources."""

    RESOLVER_PRIORITY = ["RefSeq_accession", "GenBank_accession", "JGI_Genome_ID"]

    def __init__(self):
        self._resolvers: Dict[str, callable] = {
            "RefSeq_accession": self._resolve_refseq,
            "GenBank_accession": self._resolve_genbank,
            "JGI_Genome_ID": self._resolve_jgi,
        }

    def resolve(self, genome_id_data: dict) -> Tuple[str, str]:
        """Resolves NCBI genome assembly accessions from various genome ids.

        Attempts to resolve the most recent GenBank and RefSeq accessions by trying different
        identifier types in priority order (RefSeq > GenBank > JGI).

        Args:
            genome_id_data (dict): Dictionary with genome identifiers.
                Expected keys: "RefSeq_accession", "GenBank_accession", or "JGI_Genome_ID"
                Values: Corresponding identifier strings

        Returns:
            tuple[str, str]: (genbank_acc, refseq_acc)
                - genbank_acc: GenBank assembly accession string
                - refseq_acc: RefSeq assembly accession string

        Raises:
            RuntimeError: If no valid assembly accessions could be resolved from the provided
                identifiers

        Example:
            >>> resolver = GenomeAccessionResolver()
            >>> gb_acc, rs_acc = resolver.resolve({"GenBank_accession": "GCA_000123456.1"})
        """
        for id_type in self.RESOLVER_PRIORITY:
            if id_type not in genome_id_data:
                continue

            resolver = self._resolvers.get(id_type)
            if not resolver:
                logger.warning(f"No resolver found for {id_type}")
                continue

            try:
                return resolver(genome_id_data[id_type].strip())
            except Exception as e:
                logger.warning(f"Failed to resolve {id_type}: {e}")

        raise RuntimeError("No valid assembly accessions found")

    def _resolve_refseq(self, acc: str) -> tuple[str, str]:
        self._validate_assembly_acc(acc, "refseq")
        return resolve_ncbi_assembly_accession(acc)

    def _resolve_genbank(self, acc: str) -> Tuple[str, str]:
        self._validate_assembly_acc(acc, "genbank")
        return resolve_ncbi_assembly_accession(acc)

    def _resolve_jgi(self, acc: str) -> Tuple[str, str]:
        return resolve_jgi_genome_id(acc)

    def _validate_assembly_acc(self, acc: str, acc_type: Literal["refseq", "genbank"]) -> None:
        """Validates NCBI genome assembly accession format."""
        prefix = "GCF_" if acc_type == "refseq" else "GCA_"
        if not acc.startswith(prefix):
            raise ValueError(f"Invalid {acc} assembly accession (must start with {prefix}): {acc}")

        if "." not in acc:
            raise ValueError(f"Invalid assembly accession (missing version number): {acc}")


def resolve_jgi_genome_id(jgi_genome_id: str) -> str:
    """Try to get the latest GenBank and RefSeq ID for a given JGI genome ID.

    This function queries the JGI genome lookup URL to find the corresponding
    NCBI Assembly Accession ID for the provided JGI genome ID. It then resolves
    this NCBI Assembly Accession ID to obtain the latest GenBank and RefSeq ID.

    Args:
        jgi_genome_id (str): The JGI genome ID to resolve.

    Returns:
        str: The resolved GenBank or RefSeq ID.

    Raises:
        RuntimeError: If no NCBI Assembly Accession ID can be found for the given JGI genome ID.
        httpx.HTTPStatusError: If the HTTP request to the JGI genome lookup URL ails.
    """
    url = JGI_GENOME_LOOKUP_URL.format(jgi_genome_id)
    # no User-Agent header produces a 403 Forbidden error on this site...
    resp = httpx.get(url, headers={"User-Agent": USER_AGENT}, timeout=10.0, follow_redirects=True)
    resp.raise_for_status()

    soup = BeautifulSoup(resp.content, "html.parser")
    # Find the table entry giving the "NCBI Assembly Accession" ID
    link = soup.find("a", href=re.compile("https://www.ncbi.nlm.nih.gov/datasets/genome/.*"))
    if not link:
        raise RuntimeError(
            f"Failed to find NCBI Assembly Accessions for JGI Genome ID: {jgi_genome_id}"
        )
    assembly_acc = link.text
    return resolve_ncbi_assembly_accession(assembly_acc)


def resolve_ncbi_assembly_accession(assembly_acc: str) -> Tuple[str, str]:
    """Try to get latest GenBank and RefSeq accession for a given NCBI genome assembly accession.

    This method takes an NCBI assembly accession and retrieves its revision history.
    It then determines the latest GenBank and RefSeq accessions from the revision history.

    Args:
        assembly_acc (str): The NCBI genome assembly accession to resolve (GenBank of RefSeq).

    Returns:
        Tuple[str, str]: A tuple containing the latest GenBank accession and the latest RefSeq accession.

    Raises:
        ValueError: If no assembly revision data is found for the given accession.
        ValueError: If no valid GenBank or RefSeq accession is found in the assembly revision history.
    """
    revision_history = _get_revision_history(assembly_acc)
    if not revision_history:
        raise ValueError(f"No Assembly Revision data found for {assembly_acc}")

    resolved_genbank_acc = _get_latest_accession(revision_history, "genbank_accession")
    resolved_refseq_acc = _get_latest_accession(revision_history, "refseq_accession")

    if resolved_genbank_acc or resolved_refseq_acc:
        return resolved_genbank_acc, resolved_refseq_acc

    raise ValueError("No valid GenBank or RefSeq accession found in assembly revision history")


def _get_revision_history(assembly_acc: str) -> dict:
    """Fetches revision history from NCBI Datasets API."""
    url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{assembly_acc}/revision_history"
    )

    resp = httpx.get(url, headers={"User-Agent": USER_AGENT}, timeout=10.0, follow_redirects=True)
    resp.raise_for_status()
    return resp.json()


def _get_latest_accession(
    revision_history: dict, accession_type: Literal["refseq_accession", "genbank_accession"]
) -> str:
    """Gets the latest accession for either RefSeq or GenBank."""
    revisions_with_acc = [
        entry for entry in revision_history["assembly_revisions"] if accession_type in entry
    ]

    if revisions_with_acc:
        latest_revision = max(revisions_with_acc, key=lambda x: x["release_date"])
        return latest_revision[accession_type]

    return ""
