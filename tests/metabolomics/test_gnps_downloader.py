import filecmp
from pathlib import Path
from tempfile import gettempdir
import zipfile
import pytest
from requests.exceptions import ReadTimeout
from typing_extensions import Self
from nplinker.metabolomics.gnps.gnps_downloader import GNPSDownloader
from .. import DATA_DIR


class GNPSDownloaderBuilder:
    def __init__(self):
        self._task_id = None
        self._download_root = gettempdir()

    def with_task_id(self, task_id: str) -> Self:
        self._task_id = task_id
        return self

    def with_download_root(self, download_root: Path) -> Self:
        self._download_root = download_root
        return self

    def build(self) -> GNPSDownloader:
        return GNPSDownloader(self._task_id, self._download_root)



def test_has_gnps_task_id():
    sut = GNPSDownloaderBuilder().with_task_id("c22f44b14a3d450eb836d607cb9521bb").build()
    assert sut.get_task_id() == "c22f44b14a3d450eb836d607cb9521bb"


def test_has_url():
    try:
        sut = GNPSDownloaderBuilder().with_task_id("c22f44b14a3d450eb836d607cb9521bb").build()
        assert sut.get_url() == 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=c22f44b14a3d450eb836d607cb9521bb&view=download_clustered_spectra'
    except ReadTimeout:
        pytest.skip("GNPS is down")


@pytest.mark.parametrize("task_id, filename_expected", [
    ["92036537c21b44c29e509291e53f6382", "ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip"],
    ["c22f44b14a3d450eb836d607cb9521bb", "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip"]
])
def test_downloads_file(tmp_path: Path, task_id, filename_expected):
    outpath = tmp_path.joinpath(task_id + ".zip")
    sut = GNPSDownloader(task_id, tmp_path)
    try:
        sut.download()
        actual = zipfile.ZipFile(outpath)

        expected = zipfile.ZipFile(DATA_DIR / filename_expected)

        actual_names = actual.namelist()
        expected_names = [x.filename for x in expected.filelist if x.compress_size > 0]
        assert all(item in actual_names for item in expected_names)
    except ReadTimeout:
        pytest.skip("GNPS is down")