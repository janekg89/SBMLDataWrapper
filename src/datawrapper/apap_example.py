from pathlib import Path

from pkdb_analysis import PKData

from core.experiment_factory import ExperimentFactory
from core.objects import KeyMappings




intervention_mapping = {}

substance_mapping = {}

task_mapping = {}

route_mapping = {"oral": "PO"}

key_mapping = {
    # substances

    # groups

    # measurement
    "concentration": "concentration",
    "recovery": "recovery",
    # tissue
    "urine": "urine",
    "plasma": "plasma",
    "serum": "serum",
    # units
    "gram": "g",
    "milligram": "mg",
}

apapMapping = KeyMappings(
    intervention_mapping, substance_mapping, task_mapping
)
apapMapping.set_key_mapping(key_mapping)


if __name__ == "__main__":

    # FIXME: this is a hack to avoid dependency hell
    PKDATA_ZIP_PATH = Path(
        "/home/janosch/Coding/Work/Matthias/pkdb_models/pkdb_models/models/dextromethorphan/results/pkdata/pkdb_data.zip"
    )
    Chen1996 = ExperimentFactory(
        sid="Chen1996", zip_path=PKDATA_ZIP_PATH, key_mapping=apapMapping
    )

    print(Chen1996)
