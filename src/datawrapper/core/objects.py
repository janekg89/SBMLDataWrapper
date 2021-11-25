"""
Factory which allows to create SimulationExperiments from PKDB information.

Necessary to provide the model and necessary changes.

"""
import json
from pathlib import Path
from typing import List, Dict, Any

from pydantic import BaseModel

import numpy as np

from pint import Quantity

from pkdb_analysis.data import PKData, PKDataFrame

from sbmlsim.experiment import SimulationExperiment
from sbmlsim.simulation import TimecourseSim
from sbmlsim.fit import FitMapping
from sbmlsim.simulation import Timecourse as Timecourse_sbmlsim

from datawrapper.core.key_mappings import KeyMappings


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

    def t0(self) -> float:
        """Get initial time."""
        print(self.time)
        time = self["time"].values.copy()
        raise NotImplementedError

    def tend(self) -> float:
        """Get end time."""
        raise NotImplementedError


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
        if Quantity(self.dose, self.unit) != Quantity(other.dose, other.unit):
            return False
        return True

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

    def __str__(self):
        """Get string representation."""
        return f"{self.name}_({self.group.name})"


class Task:
    interventions: List[Intervention]
    tcsim: TimecourseSim

    def __init__(
        self,
        intervention_set: List[Intervention],
        t0: float = 0,
        t_end=24,
        key_mapping=None,
        steps_per_min=1,
    ):
        self.interventions = intervention_set
        changes = self.default_changes()
        for intervention in intervention_set:
            # TODO: add dmthbr -> dmt conversion somewhere
            DOSE = Quantity(intervention.dose, intervention.unit)
            changes[
                f"{key_mapping.route_mapping[intervention.route]}DOSE_{key_mapping.task_mapping[intervention.substance]}"
            ] = DOSE

        steps = int(t_end * 3600 * steps_per_min)
        self.tcsim = TimecourseSim(
            timecourses=Timecourse_sbmlsim(
                start=t0, end=t_end, steps=steps, changes=changes
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

    def __init__(self, data: TimecourseMetaData, mapping: KeyMappings = None):
        self.data = data
        self.task = data.task
        # self.observable = Observable(
        #     utils.metadata_to_key(data, mapping), data.timecourse.unit
        # )
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
            # f"{'Observable':20} {self.observable}",
        ]
        return "\n".join(info)


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
