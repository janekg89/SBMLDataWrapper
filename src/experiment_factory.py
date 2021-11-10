"""
Factory which allows to create SimulationExperiments from PKDB information.

Necessary to provide the model and necessary changes.

"""
from pathlib import Path
from typing import List, Dict

import numpy as np
from pint import Quantity
from pydantic import BaseModel
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.simulation import TimecourseSim
from sbmlsim.fit import FitMapping

# -------------------------------------------------
# This can be done generically, from the database
class Timecourse(BaseModel):
    _label: str  # never use this for anything
    substance: str
    unit: str
    time: np.ndarray
    value: np.ndarray
    mean: np.ndarray
    sd: np.ndarray
    se: np.ndarray
    median: np.ndarray
    count: int
    # tissue ?

class Intervention(BaseModel):
    pass

class Group(BaseModel):
    pass

class Individual(BaseModel):
    pass

class TimecourseMetaData(BaseModel):
    group: None
    individual: None
    intervention: None
    timecourse: Timecourse

# -------------------------------------------------
# Here manual information is required;
# What subset of time courses?
# What parameters/conditions should be changed for what groups?
# Model specific MetaData for filtering/analysis (EM/PM), ...

class ModelDefinition(BaseModel):
    path: Path
    changes: Dict[str, Quantity]

class Observable(BaseModel):
    key: str  # e.g. [Cve_ome]
    model: ModelDefinition
    unit: str


class Mapping(BaseModel):
    # key: str  # [ome_urine_em]
    data: TimecourseMetaData
    observable: Observable

    def create_timecourse_simulation(self) -> TimecourseSim:
        """Based on Dosing"""
        tcsim = None
        return tcsim

    def create_fit_mapping(self) -> FitMapping:
        mapping = None
        return mapping


class ExperimentFactoryOptions:
    pkdb_id: str
    map_substances_to_model_observable: Dict


class ExperimentFactory:

    def __init__(self, ...):

        self.mappings: List[Mapping] = None

    def create_experiment(self) -> SimulationExperiment:
        """Uses the instance information to create a simulation experiment."""

        # FIXME: this requires instances of simulation experiments
        experiment = None

        return experiment


# Capon1996
