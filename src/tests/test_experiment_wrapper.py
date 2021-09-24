from copy import deepcopy
from timeit import timeit
from unittest import TestCase
import pytest
from pkdb_analysis import PKData
from pydantic import ValidationError
from pytest import fixture

from experiment_validator import SimulationExperimentValidator
from experiment_wrapper import CustomSimulationExperiment
from tests import TEST_PKDATA, TEST_MODEL_BASE_PATH


class TestCustomSimulationExperiment(TestCase):
    """ Test automatic merging of data with simulation Experiment"""

    def setUp(self) -> None:
        c = CustomSimulationExperiment(TEST_PKDATA, base_path=TEST_MODEL_BASE_PATH)
        c.__class__.__name__ = "Nyunt2008"
        self.nyunt2008 = c


        class CustomSimulationExperiment2(CustomSimulationExperiment):
            pass
        c2 = CustomSimulationExperiment2(TEST_PKDATA, base_path=TEST_MODEL_BASE_PATH)
        c2.__class__.__name__ = "Unknown"
        self.unknown = c2

    def test_init_wrong_init(self):
        try:
            SimulationExperimentValidator(simulation_experiment=self.unknown)
            assert 0 == 1, "ValidationError expected"
        except ValidationError:
            pass

    def test_init(self):
        SimulationExperimentValidator(simulation_experiment=self.nyunt2008)
        study_name = self.nyunt2008.pkdata().studies["name"].iloc[0]
        assert study_name == self.nyunt2008.__class__.__name__

    def test_datasets_pkdata(self):
        c = self.nyunt2008
        c.datasets_pkdata()







