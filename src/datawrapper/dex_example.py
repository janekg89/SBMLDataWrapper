from pathlib import Path

from pkdb_analysis import PKData

from datawrapper.core.experiment_factory import ExperimentFactory
from datawrapper.core.objects import KeyMappings


class DexKeyMapping(KeyMappings):

    intervention_mapping = {
        "dmthbr": "DEXHBr",
        "dmt": "DEX",
        "dextromethorphan": "DEX",
        "dextromethorphan hydrobromide": "DEXHBr",
        "quinidine": "Q",
        "quinidine sulphate": "Q",
        "qui-s": "Q",
    }

    substance_mapping = {"dmt": "dex", "dtf": "dor", "dtf-plus-dtfglu": "dor"}

    task_mapping = {
        "dmthbr": "dex",
        "dmt": "dex",
        "qui": "qui",
        "qui-s": "qui",
    }

    route_mapping = {"oral": "PO"}

    key_mapping = {
        # substances
        "dextromethorphan": "dex",
        "dextrorphan": "dor",
        "3-methoxymorphinan": "mom3",
        "3-hydroxymorphinan": "hom3",
        "3-hydroxymorphinan-glucuronide": "hom3-glu",
        "dextrorphan-glucuronide": "dor-glu",
        "hm3+dtf+dex+mom3": "dex_total",
        # groups
        "all": "all",
        "um": "UM",
        "em": "EM",
        "im": "IM",
        "pm": "PM",
        "UM": "UM",
        "EM": "EM",
        "IM": "IM",
        "PM": "PM",
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


if __name__ == "__main__":
    # from pkdb_models.models.dextromethorphan import PKDATA_ZIP_PATH
    # print(type(PKDATA_ZIP_PATH))

    # FIXME: this is a hack to avoid dependency hell
    PKDATA_ZIP_PATH = Path(
        "/home/mkoenig/git/pkdb_models/pkdb_models/models/dextromethorphan/results/pkdata/pkdb_data.zip"
    )
    Capon1996 = ExperimentFactory(
        sid="Capon1996", zip_path=PKDATA_ZIP_PATH, key_mapping=DexKeyMapping
    )

    # TODO: paracetamol example
    # TODO: outputs/correlation

    print(Capon1996)
