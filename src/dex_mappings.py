from src.key_mappings import KeyMappings


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

    substance_mapping = {
        "dmt": "dex",
        "dtf": "dor",
        "dtf-plus-dtfglu": "dor"
    }

    task_mapping = {
        "dmthbr": "dex",
        "dmt": "dex",
        "qui": "qui",
        "qui-s": "qui",
    }

    route_mapping = {
        "oral": "PO"
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