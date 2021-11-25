class KeyMappings:
    """Experiment specific information."""

    intervention_mapping = {}

    substance_mapping = {}

    task_mapping = {}

    route_mapping = {"oral": "PO"}

    key_mapping = {
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
