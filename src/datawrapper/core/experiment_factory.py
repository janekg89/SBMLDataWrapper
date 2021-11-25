"""
Factory which allows to create SimulationExperiments from PKDB information.

Necessary to provide the model and necessary changes.

"""
import json
from pathlib import Path
from typing import List, Dict, Any
from pint import Quantity
from pydantic import BaseModel

import numpy as np

from sbmlutils.factory import Q_

from pkdb_analysis.units import ureg
from pkdb_analysis.data import PKData, PKDataFrame

import utils

from sbmlsim.experiment import SimulationExperiment
from sbmlsim.simulation import TimecourseSim
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.simulation import Timecourse as Timecourse_sbmlsim

# from key_mappings import KeyMappings


class Timecourse:
    label: str  # never use this for anything
    unit: str
    time: np.ndarray
    value: np.ndarray
    mean: np.ndarray
    sd: np.ndarray
    se: np.ndarray
    median: np.ndarray
    count: int

    def __init__(self, tc):
        self.label = tc["label"]
        self.count = tc["count"]
        self.count_unit = tc["count_unit"]
        self.measurement = tc["measurement_type"]
        self.time = json.loads(tc["time"])
        self.time_unit = tc["time_unit"]
        self.mean = tc["mean"]
        self.sd = tc["sd"]
        self.unit = tc["unit"]


class Intervention:
    name: str
    substance: str
    dose: float
    unit: str
    route: str

    def __init__(self, intervention):
        self.name = intervention["substance"]
        self.substance = intervention["substance"]
        self.dose = intervention["dose"]
        self.unit = intervention["unit"]
        self.route = intervention["route"]

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self.name != other.name:
            return False
        if self.substance != other.substance:
            return False
        if self.route != other.route:
            return False
        if Q_(self.dose, self.unit) != Q_(other.dose, other.unit):
            return False
        return True

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """Get string representation."""
        return f"{self.substance}_{self.dose}{self.unit}_{self.route}"


class Group(BaseModel):
    name: str
    count: int

    def __init__(self, name, count, **data: Any):
        super().__init__(**data)
        self.name = name
        self.count = count

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """Get string representation."""
        return f"{self.name}_({self.count})"


class Individual(BaseModel):
    name: str
    group: Group

    def __init__(self, name, group, **data: Any):
        super().__init__(**data)
        self.name = name
        self.group = group

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """Get string representation."""
        return f"{self.name}_({self.group.name})"


class Task:
    interventions: List[Intervention]
    tcsim: TimecourseSim

    def __init__(
        self,
        intervention_set: List[Intervention],
        t_0=0,
        t_end=24,
        key_mapping=KeyMappings,
        steps_per_min=1,
    ):
        self.interventions = intervention_set
        changes = self.default_changes()
        for intervention in intervention_set:
            # TODO: add dmthbr -> dmt conversion somewhere
            DOSE = Q_(intervention.dose, intervention.unit)
            changes[
                f"{key_mapping.route_mapping[intervention.route]}DOSE_{key_mapping.task_mapping[intervention.substance]}"
            ] = DOSE

        steps = int(t_end * 3600 * steps_per_min)
        self.tcsim = TimecourseSim(
            timecourses=Timecourse_sbmlsim(
                start=t_0, end=t_end, steps=steps, changes=changes
            )
        )

    def default_changes(self):
        return {}

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """Get string representation."""
        return f"task_{self.interventions}"


class TimecourseMetaData:
    group: Group
    individual: Individual
    interventions: List[Intervention]
    tissue: str
    substance: str
    timecourse: Timecourse
    task: Task

    def __init__(self, tc):
        self.group = None
        self.individual = None
        if "group_name" in tc.keys():
            self.group = tc["group_name"]
        if "id_name" in tc.keys():
            self.individual = tc["id_name"]
        self.interventions = []
        for intervention in tc["interventions"]:
            self.interventions.append(Intervention(intervention))
        self.timecourse = Timecourse(tc)
        self.tissue = tc["tissue"]
        self.substance = tc["substance"]

        # self.phenotype = tc['cyp2d6 phenotype']

    def set_task(self, task: Task):
        self.task = task

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """Get string representation."""
        info = [
            "-" * 80,
            f"TimecourseMetaData:",
            "-" * 80,
            f"{'Group':20} {self.group}",
            f"{'Individual':20} {self.individual}",
            f"{'Interventions':20} {self.interventions}",
            f"{'Tissue':20} {self.tissue}",
            f"{'Substance':20} {self.substance}",
            f"{'Task':20} {self.task}",
        ]
        return "\n".join(info)


# -------------------------------------------------
# Here manual information is required;
# What subset of time courses?
# What parameters/conditions should be changed for what groups?
# Model specific MetaData for filtering/analysis (EM/PM), ...


class ModelDefinition:
    path: Path
    changes: Dict[str, Quantity]

    def __init__(self, path, changes):
        self.path = path
        self.changes = changes


class Observable:
    key: str  # e.g. [Cve_ome]
    model: ModelDefinition
    unit: str

    def __init__(self, key, unit, model=None):
        self.key = key
        self.model = model
        self.unit = unit

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """Get string representation."""
        return self.key


class Mapping:
    task: Task
    data: TimecourseMetaData
    observable: Observable
    mapping: FitMapping

    def __init__(self, data: TimecourseMetaData, mapping=KeyMappings()):
        self.data = data
        self.task = data.task
        self.observable = Observable(
            utils.metadata_to_key(data, mapping), data.timecourse.unit
        )
        # self.mapping = FitMapping() TODO for Matthias

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """Get string representation."""
        info = [
            "-" * 80,
            f"Mapping:",
            "-" * 80,
            f"{'Data':20} \n{self.data} \n",
            f"{'Observable':20} {self.observable}",
        ]
        return "\n".join(info)


class ExperimentFactoryOptions:
    pkdb_id: str
    map_substances_to_model_observable: Dict


class ExperimentFactory:
    data: List[TimecourseMetaData]
    tasks: List[Task]
    mappings: List[Mapping]

    def __init__(self, sid: str, zip_path, key_mapping, **kwargs):
        self.sid = sid
        self.zip_path = zip_path
        self.Q_ = ureg.Quantity
        self.pkdata_raw: PKData = PKData.from_archive(zip_path).filter(
            {"outputs": {"study_name": self.sid}}
        )

        # not sure if ok here
        self.tcs: PKDataFrame = self.pkdata_raw.timecourses
        self.ops: PKDataFrame = self.pkdata_raw.outputs
        self.groups: PKDataFrame = self.pkdata_raw.groups
        self.interventions: PKDataFrame = self.pkdata_raw.interventions
        self.key_mapping = key_mapping
        utils.add_group_data(self)
        utils.add_interventions(self)

        # this is needed for mapping and fitting
        self.data: List[TimecourseMetaData] = []
        self.tasks: List[Task] = []
        self.mappings: List[Mapping] = []

        for index, tc in self.tcs.iterrows():
            self.data.append(TimecourseMetaData(tc))

        # goes through all tcs (and later also outputs such as recovery)
        unique_interventions: List[List[Intervention]] = []
        for dset in self.data:
            # collects unique intervention-sets
            if dset.interventions not in unique_interventions:
                unique_interventions.append(dset.interventions)
                # creates a task fore each unique intervention
                self.tasks.append(Task(dset.interventions, key_mapping=key_mapping))
            for task in self.tasks:
                if task.interventions == dset.interventions:
                    # assigns a task to each timecourse (and later also output)
                    dset.set_task(task)
                    break
            # initialise Mapping objects
            self.mappings.append(Mapping(dset, mapping=key_mapping))

    def create_experiment(self) -> SimulationExperiment:
        """Uses the instance information to create a simulation experiment."""

        # FIXME: this requires instances of simulation experiments
        experiment = None

        return experiment

    def __repr__(self):
        return "ExperimentFactory()"

    def __str__(self):
        """Get string representation."""
        info = [
            "-" * 80,
            f"SimulationExperiment: {self.__class__.__name__}: {self.sid}",
            "-" * 80,
            f"{'Timecourses':20} \n {self.tcs}",
            f"{'Outputs':20} \n {self.ops}",
            f"{'Mappings':20}",
        ]
        for mapping in self.mappings:
            info.append(f"{mapping} \n")
        info.append(
            "-" * 80,
        )
        return "\n".join(info)
