"""
Factory which allows to create SimulationExperiments from PKDB information.

Necessary to provide the model and necessary changes.

"""
from __future__ import annotations
import json
from pathlib import Path
from typing import List, Dict, Any

from pydantic import BaseModel

import numpy as np

from pint import Quantity
from sbmlsim.simulation import TimecourseSim
from sbmlsim.fit import FitMapping
from sbmlsim.simulation import Timecourse as Timecourse_sbmlsim

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
        self.observable = Observable
        self.interventions = intervention_set
        self.tcsim = TimecourseSim
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


class MetaData:
    group: Group
    individual: Individual
    interventions: List[Intervention]
    tissue: str
    substance: str

    def __init__(
        self,
        group: Group,
        individual: Individual,
        interventions: List[Intervention],
    ):
        self.group = group
        self.individual = individual
        self.interventions = interventions

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
        ]
        return "\n".join(info)





class Timecourse:
    _keys: List[str] = [
        "label",
        "measurement_type",
        "tissue",
        "substance",
        "time",
        "time_unit",
        "mean",
        "sd",
        "se",
        "median",
        "cv",
        "unit",
    ]

    def __init__(self, tc):
        self.label: str = tc["label"]

        self.measurement_type: str = tc["measurement_type"]
        self.tissue = tc["tissue"]
        self.substance = tc["substance"]

        self.time = tc["time"]
        self.time_unit = tc["time_unit"]
        self.mean = tc["mean"]
        self.sd = tc["sd"]
        self.se = tc["se"]
        self.median = tc["median"]
        self.cv = tc["cv"]
        self.unit = tc["unit"]

    def __str__(self):
        info = [
            self.__class__.__name__,
        ]
        for key in self._keys:
            info.append(f"{key:20}{getattr(self, key)}")
        return "\n".join(info)

    def t0(self) -> float:
        """Get initial time."""
        print(self.time)
        time = self["time"].values.copy()
        raise NotImplementedError

    def tend(self) -> float:
        """Get end time."""
        raise NotImplementedError


class Data:
    _keys: List[str] = [

        "value",
        "value_err",
        "unit",
        "meta",
    ]
    def __init__(self, value: float, value_err: float, unit: str, meta: MetaData):

        self.value = value
        self.value_err = value_err
        self.unit = unit
        self.meta = meta



class TimecourseMetaData:
    group: Group
    individual: Individual
    interventions: List[Intervention]
    tissue: str
    substance: str
    timecourse: Timecourse

    def __init__(
        self,
        timecourse: Timecourse,
        group: Group,
        individual: Individual,
        interventions: List[Intervention],
    ):
        self.group = group
        self.individual = individual
        self.interventions = interventions
        self.timecourse = timecourse

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
            #f"{'Task':20} {self.task}",
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

    def __init__(self, data: Data, task: Task, observable: Observable):
        self.data = data
        self.task = task
        self.observable = observable

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
    intervention_mapping: Dict[str, str]
    substance_mapping: Dict[str, str]
    task_mapping: Dict[str, str]
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
    route_mapping: Dict[str, str] = {
        "oral": "PO"
    }

    def __init__(self, intervention_mapping, substance_mapping, task_mapping):
        self.intervention_mapping = intervention_mapping
        self.substance_mapping = substance_mapping
        self.task_mapping = task_mapping

    def set_route_mapping(self, route_mapping: Dict[str, str]):
        """ Overrides default route mappings."""
        self.route_mapping = route_mapping

    def set_key_mapping(self, key_mapping: Dict[str, str]):
        """ Overrides default key mappings."""
        self.key_mapping = key_mapping


