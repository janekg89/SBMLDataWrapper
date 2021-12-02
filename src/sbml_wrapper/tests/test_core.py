from pathlib import Path
import pytest
from sbml_wrapper.apap_example import apapMapping
from sbml_wrapper.core.experiment_factory import ExperimentFactory

def test_model_creation(tmp_path: Path, module: str) -> None:
    """Test data from outputs."""
    PKDATA_ZIP_PATH = Path(
        "/home/janek/Dev/pkdb_models/pkdb_models/models/dextromethorphan/results/pkdata/pkdb_data.zip"
    )
    chen1996 = ExperimentFactory(
        sid="Chen1996", zip_path=PKDATA_ZIP_PATH, key_mapping=apapMapping
    )
    chen1996.create_timecourse_data()





