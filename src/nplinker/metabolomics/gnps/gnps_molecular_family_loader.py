import csv
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.abc import MolecularFamilyLoaderBase
from nplinker.metabolomics.molecular_family import MolecularFamily
from nplinker.metabolomics.singleton_family import SingletonFamily


logger = LogConfig.getLogger(__file__)


class GNPSMolecularFamilyLoader(MolecularFamilyLoaderBase):
    def __init__(self, filename: str):
        self._families: list[MolecularFamily] = []
        
        for family_id, spectra_ids in _load_molecular_families(filename).items():
            if family_id == -1:
                for spectrum_id in spectra_ids:
                    family = SingletonFamily()
                    family.spectra_ids = set([spectrum_id])
                    self._families.append(family)
            else:
                family = MolecularFamily(family_id)
                family.spectra_ids = spectra_ids
                self._families.append(family)
    
    def families(self) -> list[MolecularFamily]:
        return self._families


def _load_molecular_families(filename: str) -> dict[int, set[int]]:
    """Load ids of molecular families and corresponding spectra from GNPS output file.

    Args:
        filename(str): Path to GNPS .pairsinfo file.

    Returns:
        dict[int, set[int]]: Mapping from molecular family/cluster id to the spectra ids.
    """
    logger.debug('loading edges file: %s', filename)

    families: dict = {}
    
    with open(filename, mode='rt', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        cid1_index, cid2_index, fam_index = _sniff_column_indices(filename, headers)

        for line in reader:
            spec1_id = int(line[cid1_index])
            spec2_id = int(line[cid2_index])
            family_id = int(line[fam_index])

            if families.get(family_id) is None:
                families[family_id] = set([spec1_id, spec2_id])
            else:
                families[family_id].add(spec1_id)
                families[family_id].add(spec2_id)

    return families

def _sniff_column_indices(filename: str, headers: list[str]) -> tuple[int, int, int]:
    """Get indices of required columns from the file.

    Args:
        filename(string): Path to the edges file.
        headers(string): Header line of the edges file.

    Raises:
        Exception: Raises exception if one of the required columns is not present.

    Returns:
        Tuple[int, int, int]: Tuple of indices for the required columns.
    """
    try:
        cid1_index = headers.index('CLUSTERID1')
        cid2_index = headers.index('CLUSTERID2')
        fam_index = headers.index('ComponentIndex')
    except ValueError as ve:
        message = f'Unknown or missing column(s) in edges file: {filename}'
        raise Exception(message) from ve
                
    return cid1_index,cid2_index,fam_index