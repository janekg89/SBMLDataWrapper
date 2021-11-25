import ast
import json
import sys
from typing import Dict, List

import numpy as np
from pint import Quantity

# from pkdb_models.models.dextromethorphan.experiments.base_experiment import (
#     DexSimulationExperiment,
# )


from sbmlutils import log
from sbmlsim.data import DataSet
from sbmlsim.experiment import SimulationExperiment

# from src.experiment_factory import (
#     ExperimentFactory,
#     TimecourseMetaData,
#     Timecourse,
# )  # this is circular
# from src.key_mappings import KeyMappings

logger = log.get_logger(__name__)


# def metadata_to_key(data: TimecourseMetaData, mapping: KeyMappings) -> str:
#     """
#     Creates the yid for each dset.
#     Assumes only one substance and measurement_type per dset.
#     """
#
#     substance = data.substance
#     measurement = data.timecourse.measurement
#
#     # TODO: This is dex specific: generailze.
#     if substance in ["dtf-plus-dtfglu"]:
#         logger.warning(f"Ambiguous substance found: {substance}.")
#
#     if measurement == "concentration":
#         yid = f"[Cve_{mapping.substance_mapping[substance]}]"
#     elif measurement == "cumulative amount":
#         yid = f"Aurine_{mapping.substance_mapping[substance]}"
#     else:
#         raise ValueError(f"Unexpected measurement: {measurement}.")
#
#     return yid



# deprecated
# def min_max_time(experiment: SimulationExperiment) -> (float, float):
#     """
#     Return the minimum and maximum time point of all time-courses.
#     """
#     time = experiment.tcs["time"].values.copy()
#     time = list(map(ast.literal_eval, time))
#     time_array = np.array([np.array(xi) for xi in time], dtype=object)
#     min_time = sys.float_info.max
#     max_time = sys.float_info.min
#     for element in time_array:
#         if min_time > min(element):
#             min_time = min(element)
#         if max_time < max(element):
#             max_time = max(element)
#     if min_time > 0:
#         min_time = 0.0
#     return min_time, max_time





def intervention_details(
    experiment: SimulationExperiment, medication: Dict[str, str]
) -> (str, str, str):
    substance = medication["substance"]
    if medication["route"] == "oral":
        route = "PO"
    elif medication["route"] == "iv":
        route = "IV"
    else:
        raise ValueError("Unexpected route in interventions.")

    dose = experiment.Q_(medication["dose"], medication["unit"])
    if substance == "dexhbr":
        dose = experiment.f_dexhbr_to_dex()
    return substance, route, dose


def append_task_id(experiment: SimulationExperiment, intervention: str, task_id: str):
    for key, dset in experiment.dsets.items():
        if dset["interventions"][0] == str(intervention):
            dset["task"] = task_id


def base_meta_data_dict() -> Dict:
    meta_data_dict = {
        "diplotype": None,
        "tissue": None,
        "diplotypic_phenotype": None,
        "metabolic_phenotype": None,
        "quinidine": False,
        "inhibition": False,
    }
    return meta_data_dict


def fill_meta_dict(dset) -> Dict:
    """Generates the metadata dict."""
    meta_data_dict = base_meta_data_dict()
    for meta_data_entry, item in meta_data_dict.items():
        if meta_data_entry in dset.columns:
            assert len(dset[meta_data_entry].unique()) == 1
            meta_data_dict[meta_data_entry] = dset[meta_data_entry][0]
    for intervention in dset["interventions"].unique():
        if "qui" in intervention:
            meta_data_dict["quinidine"] = True
        # FIXME: for now inhibition is always false
    return meta_data_dict


def reporting_type(dset):
    """Returns the reporting-type of a measurement (mean or value) for a given DataSet."""
    if "mean" in dset.keys():
        measured_yid = "mean"
    elif "value" in dset.key():
        measured_yid = "value"
    else:
        raise ValueError("Unexpected reporting type in dset.")
    return measured_yid


def convert_units(experiment: DexSimulationExperiment, dset: DataSet):
    Q_ = experiment.Q_

    # TODO: encapsulate -> utils
    if "mean" in dset.columns:
        subject_type = "mean"
    else:
        subject_type = "value"

    assert len(dset.measurement.unique()) == 1
    assert len(dset.substance.unique()) == 1
    assert len(dset.intervention.unique()) == 1
    assert len(dset[f"{subject_type}_unit"].unique()) == 1
    measurement = dset.measurement.unique()[0]
    substance = dset.substance.unique()[0]
    intervention = dset.intervention.unique()[0]

    # get DEX dose for each intervention (assuming exactly one dex intervention)
    doses = interventions_from_json(experiment, intervention)
    for key, intervention_dict in doses.items():
        if intervention_dict["substance"] == "dextromethorphan hydrobromide":
            # FIXME: why is multiplication of quantity by quantity not possible in this context?
            dose = float(
                float(intervention_dict["dose"].magnitude)
                * experiment.f_dexhbr_to_dex()
            )  # DEX-HBr
            DOSE = Q_(dose, intervention_dict["dose"].units)
            continue
        elif intervention_dict["substance"] == "dextromethorphan":
            DOSE = intervention_dict["dose"]
            continue
    assert DOSE

    # % -> mMol
    if measurement == "recovery":
        if (
            not Q_(dset.uinfo[subject_type]).dimensionality
            == Q_(1, "mole").dimensionality
        ):
            if substance == "dextromethorphan":
                dset.unit_conversion(
                    subject_type, DOSE / experiment.Mr.dex / Q_(100, "percent")
                )
            elif substance == "dextrorphan":
                dset.unit_conversion(
                    subject_type, DOSE / experiment.Mr.dor / Q_(100, "percent")
                )
            elif substance == "dextrorphan-glucuronide":
                dset.unit_conversion(
                    subject_type, DOSE / experiment.Mr.dorglu / Q_(100, "percent")
                )
            else:
                logger.warning(f"skipped {substance} in {dset.index}")

    # mg/mL -> mMol/mL
    if measurement == "concentration":
        if (
            not Q_(dset.uinfo[subject_type]).dimensionality
            == Q_(1, "mole/l").dimensionality
        ):
            if substance == "dextromethorphan":
                dset.unit_conversion(subject_type, 1 / experiment.Mr.dex)
            elif substance == "dextrorphan":
                dset.unit_conversion(subject_type, 1 / experiment.Mr.dor)
            elif substance == "dextrorphan-glucuronide":
                dset.unit_conversion(subject_type, 1 / experiment.Mr.dorglu)
            else:
                logger.warning(f"skipped {substance} in {dset.index}")


def dose_from_json(
    experiment: DexSimulationExperiment, interventions
) -> Dict[str, Quantity]:
    Q_ = experiment.Q_
    doses = {}
    with open(f"{experiment.data_path[0]}/{experiment.sid}/study.json") as f:
        study = json.load(f)

    for element in study["interventionset"]["interventions"]:
        if element["name"] in interventions:
            # FIXME: this will fail if different doses of same substance are given
            doses[element["substance"]] = Q_(element["value"], element["unit"])

    return doses


def interventions_from_json(
    experiment: DexSimulationExperiment, interventions
) -> Dict[str, Dict[str, Quantity]]:
    Q_ = experiment.Q_
    interventions_dict = {}
    with open(f"{experiment.data_path[0]}/{experiment.sid}/study.json") as f:
        study = json.load(f)

    for element in study["interventionset"]["interventions"]:
        if element["name"] in interventions:
            interventions_dict[element["name"]] = {
                "substance": element["substance"],
                "dose": Q_(element["value"], element["unit"]),
            }
    return interventions_dict
