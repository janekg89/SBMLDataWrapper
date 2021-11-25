"""
Factory which allows to create SimulationExperiments from PKDB information.

Necessary to provide the model and necessary changes.

"""
import json
from pathlib import Path
from typing import List, Dict, Any

from pint import Quantity

from pkdb_analysis.data import PKData, PKDataFrame


from sbmlsim.experiment import SimulationExperiment
from datawrapper.core.objects import (
    Timecourse,
    TimecourseMetaData,
    Task,
    Mapping,
    Intervention,
)

class ExperimentFactory:
    data: List[TimecourseMetaData]
    tasks: List[Task]
    mappings: List[Mapping]

    def __init__(self, sid: str, zip_path: Path, key_mapping):
        self.sid = sid
        self.pkdata: PKData = PKData.from_archive(zip_path).filter(
            {"outputs": {"study_name": self.sid}}
        )

        # not sure if ok here

        self.key_mapping = key_mapping
        self.data: List[TimecourseMetaData] = []

        self.tasks: List[Task] = []
        self.mappings: List[Mapping] = []

        # self.initialize()

    def initialize(self):
        self.tcs: PKDataFrame = self.pkdata.timecourses
        self.ops: PKDataFrame = self.pkdata.outputs
        self.groups: PKDataFrame = self.pkdata.groups
        self.interventions: PKDataFrame = self.pkdata.interventions

        self.add_group_data(self)
        self.add_interventions(self)

        for index, tc in self.tcs.iterrows():
            self.data.append(TimecourseMetaData(tc))

        # goes through all tcs (and later also outputs such as recovery)
        unique_interventions: List[List[Intervention]] = []
        for dset in self.data:
            # collects unique intervention-sets
            if dset.interventions not in unique_interventions:
                unique_interventions.append(dset.interventions)
                # creates a task fore each unique intervention
                self.tasks.append(Task(dset.interventions, key_mapping=self.key_mapping))
            for task in self.tasks:
                if task.interventions == dset.interventions:
                    # assigns a task to each timecourse (and later also output)
                    dset.set_task(task)
                    break
            # initialise Mapping objects
            self.mappings.append(Mapping(dset, mapping=self.key_mapping))


    def add_group_data(self):
        """Adds group name and count to each timecourse and each output."""
        for output_type in ["tcs", "ops"]:
            for index, content in getattr(self, output_type).iterrows():
                for group_index, group_row in self.groups.iterrows():
                    if content["group_pk"] == group_row["group_pk"]:
                        getattr(self, output_type).at[index, "count"] = group_row[
                            "group_count"
                        ]
                        getattr(self, output_type).at[
                            index, "count_unit"
                        ] = Quantity(1, "dimensionless").units
                        getattr(self, output_type).at[
                            index, "group_name"
                        ] = group_row["group_name"]

    def add_interventions(self):
        """
        Creates dictionary with intevention_pk as key that contains substance dose and unit
        of all substances in each intervention.
        :param experiment:
        :return:
        """
        self.tcs["interventions"] = ""
        self.ops["interventions"] = ""
        interventions: Dict[str, List] = {}
        for intervention_index, intervention_row in self.interventions.iterrows():
            if intervention_row["intervention_pk"] in interventions.keys():
                # TODO: DataClass instead of dict?
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
            self.tcs.at[index, "interventions"] = interventions[
                tc["intervention_pk"]]
        for index, op in self.ops.iterrows():
            self.ops.at[index, "interventions"] = interventions[
                op["intervention_pk"]]

    def create_experiment(self) -> SimulationExperiment:
        """Uses the instance information to create a simulation experiment."""

        # FIXME: this requires instances of simulation experiments
        experiment = None
        return experiment

    def __str__(self):
        """Get string representation."""
        info = [
            "-" * 80,
            f"{self.__class__.__name__}: {self.sid}",
            "-" * 80,
            f"{'Data':20}{self.data}",
            f"{'Tasks':20}{self.tasks}",
            f"{'Mappings':20}{self.mappings}",
        ]
        for mapping in self.mappings:
            info.append(f"{mapping} \n")
        info.append(
            "-" * 80,
        )
        return "\n".join(info)


