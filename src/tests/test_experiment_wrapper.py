from unittest import TestCase
import pytest
from pytest import fixture

from experiment_wrapper import CustomSimulationExperiment
from tests import TEST_PKDATA


class TestCustomSimulationExperiment(TestCase):
    """ Test automatic merging of data with simulation Experiment"""

    def setUp(self) -> None:
        c = CustomSimulationExperiment(TEST_PKDATA,)
        c.__class__.__name__ = "Nyunt2008"
        self.nyunt2008 = c

    def test_init(self):
        c = self.nyunt2008


    def test_datasets_pkdata(self):
        c = self.nyunt2008
        c.datasets_pkdata()





