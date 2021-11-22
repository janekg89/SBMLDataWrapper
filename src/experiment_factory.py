"""
Factory which allows to create SimulationExperiments from PKDB information.

Necessary to provide the model and necessary changes.

"""
import json
from pathlib import Path
from typing import List, Dict, Any, Optional
from pkdb_analysis.data import PKData, PKDataFrame

import numpy as np
from pkdb_analysis.units import ureg
from sbmlutils.factory import Q_

import utils
from pint import Quantity
from pydantic import BaseModel
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.simulation import TimecourseSim
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.simulation import Timecourse as Timecourse_sbmlsim

from dex_mappings import DexKeyMapping

STEPS_PER_SEC = 1.0

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
    
    def __init__(self, tc, **data: Any):
        self.label = tc['label']
        self.count = tc['count']
        self.count_unit = tc['count_unit']
        self.measurement = tc['measurement_type']
        self.time = json.loads(tc['time'])
        self.time_unit = tc['time_unit']
        self.mean = tc['mean']
        self.sd = tc['sd']
        self.unit = tc['unit']


class Intervention:
    name: str
    substance: str
    dose: float
    unit: str
    route: str

    def __init__(self, intervention, **data: Any):
        self.name = intervention["substance"]
        self.substance = intervention["substance"]
        self.dose = intervention["dose"]
        self.unit = intervention["unit"]
        self.route = intervention["route"]

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


class Individual(BaseModel):
    name: str
    group: Group

    def __init__(self, name, group,**data: Any):
        super().__init__(**data)
        self.name = name
        self.group = group


class TimecourseMetaData:
    group: Group
    individual: Individual
    interventions: List[Intervention]
    tissue: str
    substance: str
    timecourse: Timecourse

    def __init__(self, tc):
        if "group_name" in tc.keys():
            self.group = tc["group_name"]
        if "id_name" in tc.keys():
            self.individual = tc["id_name"]
        self.interventions = []
        for intervention in tc["interventions"]:
            self.interventions.append(Intervention(intervention))
        self.timecourse = Timecourse(tc)
        self.tissue = tc['tissue']
        self.substance = tc['substance']

        # self.phenotype = tc['cyp2d6 phenotype']

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

class Mapping:
    key: str  # [ome_urine_em]
    data: TimecourseMetaData
    observable: Observable
    tcsim: TimecourseSim
    mapping: FitMapping

    def __init__(self, tc: Timecourse):
        self.data = TimecourseMetaData(tc)
        self.key = utils.metadata_to_key(self.data, DexKeyMapping())
        self.observable = Observable(self.key, self.data.timecourse.unit)
        self.tcsim = self.create_timecourse_simulation()

    def create_timecourse_simulation(self) -> TimecourseSim:
        """Based on Dosing"""
        tcsim = None
        for intervention in self.data.interventions:
            # FIXME: how to handle multiple interventions?
            DOSE = Q_(intervention.dose, intervention.unit)
            t_end = self.data.timecourse.time[-1]
            steps = int(t_end * 3600 * STEPS_PER_SEC)
            tcsim = TimecourseSim(timecourses=Timecourse_sbmlsim(
                    start=0, end=t_end, steps=steps,
                    changes={
                        **self.default_changes(),
                        f"{intervention.route}DOSE_{intervention.substance}": DOSE,
                    }
                )
            )
        return tcsim

    def create_fit_mapping(self) -> FitMapping:
        # TODO: FitMapping requiers name of task and name of data set. This should be changed.
        mapping = None
        return mapping

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """Get string representation."""
        info = [
            f"{'Timecourse':20} {self.data.timecourse.label}",
            f"{'Key':20} {self.key}",
        ]
        return "\n".join(info)

    def default_changes(self):
        return {}


class ExperimentFactoryOptions:
    pkdb_id: str
    map_substances_to_model_observable: Dict


class ExperimentFactory:
    mappings: List[Mapping]

    def __init__(self, sid: str, zip_path, **kwargs):
        self.sid = sid
        self.zip_path = zip_path
        self.Q_ = ureg.Quantity
        self.pkdata_raw: PKData = PKData.from_archive(zip_path).filter({"outputs": {"study_name": self.sid}})

        # not sure if ok here
        self.tcs: PKDataFrame = self.pkdata_raw.timecourses
        self.ops: PKDataFrame = self.pkdata_raw.outputs
        self.groups: PKDataFrame = self.pkdata_raw.groups
        self.interventions: PKDataFrame = self.pkdata_raw.interventions
        utils.add_group_data(self)
        utils.add_interventions(self)

        self.mappings: List[Mapping] = []

        for index, tc in self.tcs.iterrows():
            self.mappings.append(Mapping(tc))

    def create_experiment(self) -> SimulationExperiment:
        """Uses the instance information to create a simulation experiment."""

        # FIXME: this requires instances of simulation experiments
        experiment = None

        return experiment

    def __repr__(self):
        return "ExperimentFactory()"

    def __str__(self):
        """Get string representation."""
        info = ["-" * 80, f"SimulationExperiment: {self.__class__.__name__}: {self.sid}",
                "-" * 80,
                f"{'Timecourses':20} \n {self.tcs}",
                f"{'Outputs':20} \n {self.ops}",
                f"{'Mappings':20}",
                ]
        for mapping in self.mappings:
            info.append(f"{mapping} \n")
        info.append("-" * 80,)
        return "\n".join(info)


