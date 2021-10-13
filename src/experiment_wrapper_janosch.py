import ast
import json
import logging
import sys
import warnings
from typing import Dict, List

import coloredlogs
import numpy as np
import pandas as pd
import utils
from pandas import DataFrame
from pint import Quantity
from sbmlsim.data import DataSet, load_pkdb_dataframes_by_substance
from sbmlsim.experiment import SimulationExperiment, ExperimentDict
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.simulation import TimecourseSim, Timecourse, AbstractSim
from sbmlsim.task import Task
from sbmlsim.units import UnitsInformation

from pkdb_models.models.dextromethorphan import PKDATA_ZIP_PATH
from pkdb_models.models.dextromethorphan.experiments.base_experiment import DexSimulationExperiment, DexMappingMetaData
from pkdb_analysis.data import PKData, PKDataFrame



coloredlogs.install(
    level="INFO",
    fmt='%(levelname)s %(pathname)s:%(lineno)d %(message)s'
)
logger = logging.getLogger(__name__)

intervention_mapping = {
    # interventions
    "dmthbr": "DEXHBr",
    "dmt": "DEX",
    "dextromethorphan": "DEX",
    "dextromethorphan hydrobromide": "DEXHBr",
    "quinidine": "Q",
    "quinidine sulphate": "Q",
    "qui-s": "Q",
}

substance_mapping = {
    # interventions
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


# TODO:     inheritance seems not well structured for me.
#           Maybe DexSimulationExperiment(PKDataSimulationExperiment)
#           and PKDataSimulationExperiment(SimulationExperiment would be better)
class PKDataSimulationExperiment(DexSimulationExperiment):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.sid = kwargs["sid"]
        print(kwargs)
        print(self)
        # load data
        self.pkdata_raw: PKData = PKData.from_archive(PKDATA_ZIP_PATH).filter({"outputs": {"study_name": self.sid}})
        self.tcs: PKDataFrame = self.pkdata_raw.timecourses
        self.ops: PKDataFrame = self.pkdata_raw.outputs
        self.groups: PKDataFrame = self.pkdata_raw.groups
        self.interventions: PKDataFrame = self.pkdata_raw.interventions

        # manipulations on tcs, ops
        self.add_measurement("cyp2d6 phenotype")
        self.add_group_data()
        self.add_interventions()

        print(self.ops.keys())

        self.dsets: Dict[str, DataSet] = self.datasets_from_pkdata()

        # manipulation on dsets
        self.add_yid()

        # simulation, fitting and evaluation
        self.simulation_tasks = self.set_simulation_tasks()

        # FIXME: docu says nono
        self.initialize()  # this should be redundant when each experiment is realsed as one single instance of the same class

        self.fit_mappings = self.set_fit_mappings()
        self.plots = {}

    def get_dsets(self):
        return self.dsets

    def get_simulation_tasks(self):
        return self.simulation_tasks

    def get_fit_mappings(self):
        return self.fit_mappings

    def get_plots(self):
        return self.plots

    def datasets_from_pkdata(self, ) -> Dict[str, DataSet]:
        """ creates datasets from pkdata object."""
        datasets = {}

        # TODO also for ops (only timed data like concentration, cumulative amount, recovery etc.)
        for index, tc in self.tcs.iterrows():
            name, dset = self.create_tc_dset(tc)
            datasets[name] = dset
        return datasets

    # TODO: maybe first create dset without phenotype, groupdate etc and add additional columns afterwards
    def create_tc_dset(self, tc) -> (str, DataSet):
        time = json.loads(tc['time'])
        # def __dict__
        # pd.merge
        dset = {
            "study": tc['study_name'],
            "label": tc['label'],
            "group": tc['group_name'],
            "phenotype": tc['cyp2d6 phenotype'],
            "count": tc['count'],
            "count_unit": tc['count_unit'],
            "measurement": tc['measurement_type'],
            "substance": tc['substance'],
            "tissue": tc['tissue'],
            "interventions": str(tc['interventions']),
            "time": time,
            "time_unit": tc['time_unit'],
            "mean": tc['mean'],
            "mean_sd": tc['sd'],
            "mean_unit": tc['unit'],
        }
        name = self.dset_to_name(tc)
        # TODO: adapt and enncapulate; all_udict as global parameter in separate file? see sbmlsim.data.py
        all_udict: Dict[str, str] = {}
        for key in dset.keys():
            print(key)
            key = str(key)
            # handle '*_unit columns'
            if key.endswith("_unit"):
                # parse the item and unit in dict
                units = dset[key].unique()
                if len(units) > 1:
                    logger.error(
                        f"Column '{key}' units are not unique: '{units}' in \n" f"{dset}"
                    )
                elif len(units) == 0:
                    logger.error(f"Column '{key}' units are missing: '{units}'")
                    print(dset.head())
                item_key = key[0:-5]
                if item_key not in dset.columns:
                    logger.error(
                        f"Missing * column '{item_key}' for unit " f"column: '{key}'"
                    )
                else:
                    all_udict[item_key] = units[0]
        d = DataSet(dset)
        d.uinfo = UnitsInformation(all_udict, ureg=self.ureg)  # TODO: add all_udict
        d.Q_ = d.uinfo.ureg.Quantity
        return name, d

    # TODO: calculate amount from recovery an replace measurement_type==recovery by cumulative amount
    def add_yid(self):
        """
        Creates the yid for each dset and adds it as a column. This yid can later be used for FitMappings.
        Assumes only one substance, measurment_type, tissue per dset.
        """
        for key, dset in self.dsets.items():
            substance = dset.substance.unique()[0]
            measurement = dset.measurement.unique()[0]
            tissue = dset.tissue.unique()[0]

            if substance in ["dtf-plus-dtfglu"]:
                logger.warning(f"Ambiguous substance found: {substance}.")

            if measurement == "concentration":
                yid = f"[Cve_{substance_mapping[substance]}]"
            elif measurement == "cumulative amount":
                yid = f"Aurine_{substance_mapping[substance]}"
            else:
                raise ValueError(f"Unexpected measurement: {measurement}.")

            dset["yid"] = yid

    def add_measurement(self, additional_measurement: str):
        """
            Gets measurement type from group definition and adds it to respective outputs.
        """
        for output_type in ["tcs", "ops"]:
            result = "NaN"
            for index, content in getattr(self, output_type).iterrows():
                for group_index, group_row in self.groups.iterrows():
                    if content["group_pk"] == group_row["group_pk"] and group_row[
                        'measurement_type'] == additional_measurement:
                        result = group_row["choice"]
                        break
                getattr(self, output_type).at[index, additional_measurement] = result

    def add_group_data(self):
        """
            Adds group name and count to each timecourse and each output
        """
        for output_type in ["tcs", "ops"]:
            for index, content in getattr(self, output_type).iterrows():
                for group_index, group_row in self.groups.iterrows():
                    if content["group_pk"] == group_row["group_pk"]:
                        getattr(self, output_type).at[index, "count"] = group_row["group_count"]
                        getattr(self, output_type).at[index, "count_unit"] = self.Q_(1, "dimensionless").units
                        getattr(self, output_type).at[index, "group_name"] = group_row["group_name"]

    # TODO: encapsulate functional parts
    def add_interventions(self):
        # create dictionary with intevention_pk as key that contains substance dose and unit of all substances in each intervention
        self.tcs["interventions"] = ""
        self.ops["interventions"] = ""
        interventions: Dict[str, List] = {}
        for intervention_index, intervention_row in self.interventions.iterrows():
            if intervention_row["intervention_pk"] in interventions.keys():
                interventions[intervention_row["intervention_pk"]].append(
                    {
                        "substance": intervention_row["substance"],
                        "dose": intervention_row["value"],
                        "unit": intervention_row["unit"],
                        "route": intervention_row["route"],
                    }
                )
            else:
                interventions[intervention_row["intervention_pk"]] = [
                    {
                        "substance": intervention_row["substance"],
                        "dose": intervention_row["value"],
                        "unit": intervention_row["unit"],
                        "route": intervention_row["route"],
                    }
                ]
        # add dictionary entry to each tc/op
        for index, tc in self.tcs.iterrows():
            self.tcs.at[index, "interventions"] = interventions[tc["intervention_pk"]]
        for index, op in self.ops.iterrows():
            self.ops.at[index, "interventions"] = interventions[op["intervention_pk"]]

    def dset_to_name(self, output) -> str:
        """
            Create name that contains intervention, substance, tissue
        """
        name = f"{output['group_name']}_{output['cyp2d6 phenotype']}_{output['substance']}_{output['tissue']}"
        for intervention in output["interventions"]:
            intervention_key = self.get_intervention_key(intervention=intervention)
            name += f"_{intervention_key}"
        return name

    def get_intervention_key(self, intervention) -> str:
        """
            Return a unique string for each intervention that contains information on the given substance and the
            respective amount in mg.
        """
        dose = self.Q_(intervention['dose'], intervention['unit']).m_as('mg')
        if intervention['substance'] in intervention_mapping:
            key = f"{intervention_mapping[intervention['substance']]}{int(dose)}mg"
        else:
            key = f"{intervention['substance']}{int(dose)}mg"
        return key

    # TODO: this method is too long -> encapsulate functional parts
    def set_simulation_tasks(self) -> Dict[str, TimecourseSim]:
        """
            Returns a dictionary with one simulation task for each unique intervention.
        """
        tcsims = {}

        # TODO: encapsulate
        interventions = []
        for intervention in self.tcs["interventions"]:
            if intervention not in interventions:
                interventions.append(intervention)

        # This finds the minimum and maximum time point of measurement
        # TODO: this has to be done simpler (or encapulate)
        time = self.tcs['time'].values.copy()
        time = list(map(ast.literal_eval, time))
        time_array = np.array([np.array(xi) for xi in time], dtype=object)
        min_time = sys.float_info.max
        max_time = sys.float_info.min
        for element in time_array:
            if min_time > min(element):
                min_time = min(element)
            if max_time < max(element):
                max_time = max(element)
        if min_time > 0:
            min_time = 0.0

        # this iterates over list  of lists
        for intervention in interventions:
            id = ""
            changes = {}
            for medication in intervention:
                id += self.get_intervention_key(medication)
                substance = medication["substance"]

                # TODO: encapsulate
                if medication['route'] == "oral":
                    route = "PO"
                elif medication['route'] == "iv":
                    route = "IV"
                else:
                    raise ValueError("Unexpected route in interventions.")

                dose = self.Q_(medication["dose"], medication["unit"])
                if substance == "dexhbr":
                    dose = self.f_dexhbr_to_dex()
                changes[f"{route}DOSE_{task_mapping[substance]}"] = dose

            task_dict = {
                "start": min_time,
                "end": max_time * 1.1,  # +10% simulation time
                "steps": int(60 * (max_time * 1.1 - min_time)),  # one step per min   TODO: step per min as param
                "changes": changes,
                # "model_changes": Dict[str, Quantity],
                # "model_manipulations": dict ,
                # "discard": bool,
            }

            tcsims[id] = TimecourseSim(timecourses=Timecourse(**task_dict))
        return tcsims

    def set_fit_mappings(self) -> Dict[str, FitMapping]:
        """
            Iterates over all outputs and returns a dictionary with one entry for each timecourse/output.
        """
        mappings = {}
        for key, dset in self.dsets.items():
            assert len(dset['tissue'].unique()) == 1
            tissue = dset['tissue']
            task_yid = dset['yid']

            if "mean" in dset.keys():
                measured_yid = "mean"
            elif "value" in dset.key():
                measured_yid = "value"
            else:
                raise ValueError("Unexpected reporting type in dset.")

            print(self)
            print(self.datasets().keys())
            print(key)

            print(type(self.dsets[key]))

            print(self.datasets()["PM_pm_dtf-plus-dtfglu_plasma_DEXHBr30mg"].uinfo)

            test = FitData(
                experiment=self,
                dataset=key,
                xid="time",
                yid=measured_yid,
                count="count"
            )

            print(test)

            mappings[f"fm_{key}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=key,
                    xid="time",
                    yid=measured_yid,
                    count="count"
                ),
                observable=FitData(
                    self,
                    task="task_po30",  # TODO: map from key, interventions ...?
                    xid="time",
                    yid=task_yid  # TODO: get from dset (create function: measurement type -> str)
                ),
                metadata=DexMappingMetaData(
                    diplotype=None,  # TODO: get from dset
                    tissue=tissue,  # TODO: get from dset
                    diplotypic_phenotype=None,  # TODO: get from dset
                    metabolic_phenotype="PM",  # TODO: get from dset
                    quinidine=False,  # TODO: get from dset
                    inhibition=False,  # TODO: get from dset
                )
            )
        return mappings


def from_csv(experiment: DexSimulationExperiment, fig_ids: List[str], **kwargs) -> Dict[
    str, DataSet]:
    # TODO:
    #   - Load info from study.json -> map intervention to dose
    #       - intervention: substance, dose [x]
    #       - group: phenotype [ ]

    # generate dicts for specified figures/tables (split by substance)
    sheets = {}
    for fig_id in fig_ids:
        dframes = load_pkdb_dataframes_by_substance(f"{experiment.sid}_{fig_id}", data_path=experiment.data_path)
        sheets[f"{fig_id}"] = dframes

    # iterate over sheets
    dsets = {}
    sheet: Dict[str, DataFrame]
    for sheet_id, sheet in sheets.items():
        # iterate over substances substance
        for substance, df in sheet.items():
            if substance in key_mapping.keys():
                d = {f"{sheet_id}_{key_mapping[substance]}": DataSet.from_df(df, experiment.ureg)}
            else:
                warnings.warn(f"{substance} not in mapping keys")
                continue
                # d = {f"{sheet_id}_{substance}": DataSet.from_df(df, experiment.ureg)}
            # split each data frame in single time-courses
            for split_key in ["group", "id", "measurement", "tissue", "intervention"]:
                d = split_dataset(experiment, d, split_key)
                # TODO: remove empty columns (including key)
            # TODO: remove empty datasets
            dsets.update(d)

    # convert units to molar
    dset: DataSet
    for key, dset in dsets.items():
        convert_units(experiment, dset)
    return dsets


def split_dataset(experiment, dsets: Dict[str, DataSet], split_key: str) -> Dict[str, DataSet]:
    dsets_split = {}
    for key, dset in dsets.items():
        if split_key in dset.columns.values:
            for splitter in getattr(dset, split_key).unique():
                name = dataset_name(experiment, splitter, split_key)
                dsets_split[f"{key}_{name}"] = DataSet.from_df(dset[getattr(dset, split_key) == splitter],
                                                               experiment.ureg)
        else:
            return dsets
    return dsets_split


def dataset_name(experiment, splitter: str, split_key: str) -> str:
    name = ""
    if split_key == "intervention":
        interventions = splitter.split(',')
        # split multiple interventions
        for intervention in interventions:
            intervention.replace(' ', '')
        doses = dose_from_json(experiment, splitter)
        for medication, dose in doses.items():
            if medication in ["quinidine", "quinidine sulphate"]:
                name += "Q"
            elif medication in ["dextromethorphan", "dextromethorphan hydrobromide"]:
                name += f"{intervention_mapping[medication]}{dose.magnitude}{key_mapping[str(dose.units)]}"
            else:
                warnings.warn("Non standard intervention.")
                name += f"{medication}{dose.magnitude}{key_mapping[str(dose.units)]}"
    else:
        if splitter in key_mapping.keys():
            name += key_mapping[splitter]
        else:
            warnings.warn(f"{splitter} not in mapping keys")
            name = splitter
    return name


def convert_units(experiment: DexSimulationExperiment, dset: DataSet):
    Q_ = experiment.Q_

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
            dose = float(float(intervention_dict["dose"].magnitude) * experiment.f_dexhbr_to_dex())  # DEX-HBr
            DOSE = Q_(dose, intervention_dict["dose"].units)
            continue
        elif intervention_dict["substance"] == "dextromethorphan":
            DOSE = intervention_dict["dose"]
            continue
    assert DOSE

    # % -> mMol
    if measurement == "recovery":
        if not Q_(dset.uinfo[subject_type]).dimensionality == Q_(1, "mole").dimensionality:
            if substance == "dextromethorphan":
                dset.unit_conversion(subject_type, DOSE / experiment.Mr.dex / Q_(100, "percent"))
            elif substance == "dextrorphan":
                dset.unit_conversion(subject_type, DOSE / experiment.Mr.dor / Q_(100, "percent"))
            elif substance == "dextrorphan-glucuronide":
                dset.unit_conversion(subject_type, DOSE / experiment.Mr.dorglu / Q_(100, "percent"))
            else:
                warnings.warn(f"skipped {substance} in {dset.index}")

    # mg/mL -> mMol/mL
    if measurement == "concentration":
        if not Q_(dset.uinfo[subject_type]).dimensionality == Q_(1, "mole/l").dimensionality:
            if substance == "dextromethorphan":
                dset.unit_conversion(subject_type, 1 / experiment.Mr.dex)
            elif substance == "dextrorphan":
                dset.unit_conversion(subject_type, 1 / experiment.Mr.dor)
            elif substance == "dextrorphan-glucuronide":
                dset.unit_conversion(subject_type, 1 / experiment.Mr.dorglu)
            else:
                warnings.warn(f"skipped {substance} in {dset.index}")


def dose_from_json(experiment: DexSimulationExperiment, interventions) -> Dict[str, Quantity]:
    Q_ = experiment.Q_
    doses = {}
    with open(f"{experiment.data_path[0]}/{experiment.sid}/study.json") as f:
        study = json.load(f)

    for element in study["interventionset"]["interventions"]:
        if element["name"] in interventions:
            # FIXME: this will fail if different doses of same substance are given
            doses[element["substance"]] = Q_(element["value"], element["unit"])

    return doses


def interventions_from_json(experiment: DexSimulationExperiment, interventions) -> Dict[str, Dict[str, Quantity]]:
    Q_ = experiment.Q_
    interventions_dict = {}
    with open(f"{experiment.data_path[0]}/{experiment.sid}/study.json") as f:
        study = json.load(f)

    for element in study["interventionset"]["interventions"]:
        if element["name"] in interventions:
            interventions_dict[element["name"]] = {
                "substance": element["substance"],
                "dose": Q_(element["value"], element["unit"])
            }
    return interventions_dict
