"""
Factory which allows to create SimulationExperiments from PKDB information.

Necessary to provide the model and necessary changes.

"""
from pathlib import Path
from typing import List, Dict

import numpy as np
import pandas as pd
from rich import print
from pint import Quantity

from pkdb_analysis.data import PKData, PKDataFrame

from sbmlutils import log
from sbmlsim.experiment import SimulationExperiment
from src.sbml_wrapper.core.objects import (
    Timecourse,
    TimecourseMetaData,
    Task,
    Mapping,
    Intervention,
    Individual,
    Group, Observable, Data, MetaData,
)

logger = log.get_logger(__name__)


class ValueType:
    VALUE = "value"
    MEAN = "mean"
    MEDIAN = "median"


class ErrType:
    SD = "sd"
    SE = "se"
    CV = "cv"

class ExperimentFactory:
    #data: List[TimecourseMetaData]
    #tasks: List[Task]
    #mappings: List[Mapping]

    def __init__(self, sid: str, zip_path: Path, key_mapping, observables):
        self.sid = sid
        self.pkdata: PKData = PKData.from_archive(zip_path).filter(
            {"outputs": {"study_name": self.sid}}
        )
        # not sure if ok here

        self.key_mapping = key_mapping
        self.task_intervention_mapping = {}

        self.data: List[Data] = self.create_data()

        self.tasks: List[Task] = []
        self.mappings: List[Mapping] = []
        self.observables: List[Observable] = observables


        #self.create_tasks()
        #self.create_mappings()
        # self.initialize()


    def create_data_single(self, output: pd.Series) -> Data:
        # data_dict = {"timecourses": self.create_timecourse_data()}
        #  data_dict
        meta_kwargs = {}
        subject_pk = output[["group_pk", "individual_pk"]].astype(int).idxmax()
        if subject_pk == "individual_pk":
            individual_data = self.pkdata.individuals[self.pkdata.individuals[subject_pk] == output[subject_pk]]

            meta_kwargs["individual"] = Individual(name=individual_data.individual_name.iloc[0])

            value = output.value
            value_type = ValueType.VALUE
            err = np.NAN
            err_type = None

        else:
            value_type = output[[ValueType.MEAN, ValueType.MEDIAN]].idxmax()
            value = output[value_type]

            err_type = output[[ErrType.SD, ErrType.SE, ErrType.CV]].idxmax()
            err = output[err_type]
            group_data = self.pkdata.groups[self.pkdata.groups[subject_pk] == output[subject_pk]]
            meta_kwargs["group"] = Group(name=group_data.name, count=group_data.count)
        interventions_data = self.pkdata.interventions[self.pkdata.interventions["intervention_pk"] == output["intervention_pk"]]
        interventions=[]
        for i, intervention_data in interventions_data.iterrows():
            intervention = Intervention(
                name=intervention_data["name"],
                substance=intervention_data["substance"],
                dose=intervention_data["value"],
                unit=intervention_data["unit"],
                route=intervention_data["route"]
            )
            interventions.append(intervention)

        meta_kwargs["interventions"] = interventions

        meta = MetaData(tissue=output["tissue"],substance=output["substance"], **meta_kwargs)

        return Data(
            value=value,
            value_type=value_type,
            err=err,
            err_type=err_type,
            unit=output.unit,
            meta=meta
            )

    def create_data(self) -> List[Data]:
        data  = []
        for index, output in self.pkdata.outputs.iterrows():
            data_single = self.create_data_single(output)
            data.append(data_single)
        return data

    def create_timecourse_data(self):
        """Creates the data objects."""

        for index, pkdb_tc in self.pkdata.timecourses.iterrows():
            #timecourse_outputs = self.pkdata.outputs.df[self.pkdata.outputs["output_pk"].isin(pkdb_tc.output_pk)]
            timecourse = Timecourse(pkdb_tc)
            logger.info(timecourse)

            # FIXME: use pkdb_data to read this
            individual: Individual = None
            group: Group = None
            interventions: List[Intervention] = []

            tc_metadata = TimecourseMetaData(
                timecourse=timecourse,
                individual=individual,
                group=group,
                interventions=interventions,
            )

            self.data.append(tc_metadata)

    def create_tasks(self) -> List[Task]:
        """Creates the task objects."""
        tasks = []
        for (intervention_pk, intervention_set) in self.pkdata.interventions.groupby("intervention_pk"):
            outputs = self.pkdata.outputs
            t_end = outputs.time.max()
            t_min = outputs.time.min()

            task = Task()
            self.task_interventionset_mapping[intervention_pk] = task.sid
            tasks.append(task)

        return tasks

    def create_task(self) -> Task:
        """Creates one task object."""
        raise NotImplementedError

    def create_mapping(self) -> Mapping:
        data = self.get_data()
        data = self.observable()
        return Mapping(data=data)

    def create_mappings(self):
        """Creates the mapping objects."""
        for data_single in self.data:
            self.create_data()
        task = self.get_task()


        pass

    def initialize(self):
        """"""
        # FIXME: this is create_tasks and mappings
        # goes through all tcs (and later also outputs such as recovery)
        unique_interventions: List[List[Intervention]] = []
        for dset in self.data:
            # collects unique intervention-sets
            if dset.interventions not in unique_interventions:
                unique_interventions.append(dset.interventions)
                # creates a task fore each unique intervention
                self.tasks.append(
                    Task(dset.interventions, key_mapping=self.key_mapping)
                )

            # for task in self.tasks:
            #     if task.interventions == dset.interventions:
            #         # assigns a task to each timecourse (and later also output)
            #         dset.set_task(task)
            #         break
            # initialise Mapping objects
            self.mappings.append(Mapping(dset, mapping=self.key_mapping))

    # def add_interventions(self):
    #     """
    #     Creates dictionary with intevention_pk as key that contains substance dose and unit
    #     of all substances in each intervention.
    #     :param experiment:
    #     :return:
    #     """
    #     self.tcs["interventions"] = ""
    #     self.ops["interventions"] = ""
    #     interventions: Dict[str, List] = {}
    #     for intervention_index, intervention_row in self.interventions.iterrows():
    #         if intervention_row["intervention_pk"] in interventions.keys():
    #             # TODO: DataClass instead of dict?
    #             interventions[intervention_row["intervention_pk"]].append(
    #                 {
    #                     "substance": intervention_row["substance"],
    #                     "dose": intervention_row["value"],
    #                     "unit": intervention_row["unit"],
    #                     "route": intervention_row["route"],
    #                 }
    #             )
    #         else:
    #             interventions[intervention_row["intervention_pk"]] = [
    #                 {
    #                     "substance": intervention_row["substance"],
    #                     "dose": intervention_row["value"],
    #                     "unit": intervention_row["unit"],
    #                     "route": intervention_row["route"],
    #                 }
    #             ]
    #     # add dictionary entry to each tc/op
    #     for index, tc in self.tcs.iterrows():
    #         self.tcs.at[index, "interventions"] = interventions[tc["intervention_pk"]]
    #     for index, op in self.ops.iterrows():
    #         self.ops.at[index, "interventions"] = interventions[op["intervention_pk"]]

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
