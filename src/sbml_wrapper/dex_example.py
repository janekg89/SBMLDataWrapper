from pathlib import Path

from pkdb_analysis import PKData

from core.experiment_factory import ExperimentFactory
from core.objects import KeyMappings

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
observables = {
    "[Cve_dex]"
}

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

dex_mapping = KeyMappings(
    intervention_mapping, substance_mapping, task_mapping
)
dex_mapping.set_key_mapping(key_mapping)

if __name__ == "__main__":
    #from pkdb_models.models.dextromethorphan import PKDATA_ZIP_PATH
    #print(type(PKDATA_ZIP_PATH))

    # FIXME: this is a hack to avoid dependency hell
    #PKDATA_ZIP_PATH = Path(
    #    "/home/janosch/Coding/Work/Matthias/pkdb_models/pkdb_models/models/dextromethorphan/results/pkdata/pkdb_data.zip"
    #)
    PKDATA_ZIP_PATH = Path(
                "/home/janek/Dev/pkdb_models/pkdb_models/models/dextromethorphan/results/pkdata/pkdb_data.zip"
    )
    #capon1996 = ExperimentFactory(
    #    sid="Capon1996", zip_path=PKDATA_ZIP_PATH, key_mapping=dex_mapping
    #)

    lopez2005 = ExperimentFactory(
        sid="Lopez2005", zip_path=PKDATA_ZIP_PATH, key_mapping=dex_mapping
    )
    print(lopez2005)

    # TODO: paracetamol example
    # TODO: outputs/correlation

