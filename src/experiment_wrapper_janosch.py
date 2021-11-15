import json
import logging
from typing import Dict, List

import coloredlogs
from pkdb_models.models.dextromethorphan.experiments.genotype_phenotype import colors_phenotype
from sbmlsim.plot import Figure, Axis

import utils
from sbmlsim.data import DataSet
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.simulation import TimecourseSim, Timecourse
from sbmlsim.units import UnitsInformation

from pkdb_models.models.dextromethorphan import PKDATA_ZIP_PATH
from pkdb_models.models.dextromethorphan.experiments.base_experiment import DexSimulationExperiment, DexMappingMetaData
from pkdb_analysis.data import PKData, PKDataFrame

coloredlogs.install(
    level="INFO",
    fmt='%(levelname)s %(pathname)s:%(lineno)d %(message)s'
)
logger = logging.getLogger(__name__)

STEPS_PER_MIN = 1.0

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


# TODO: - Inheritance seems not well structured for me.
#         Maybe DexSimulationExperiment(PKDataSimulationExperiment)
#         and PKDataSimulationExperiment(SimulationExperiment) would be better
#       - Add automatic unit conversion
class PKDataSimulationExperiment(DexSimulationExperiment):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.sid = kwargs["sid"]
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

        self.dsets: Dict[str, DataSet] = self.datasets_from_pkdata()

        # manipulation on dsets
        self.add_yid()

        # simulation, fitting and evaluation
        self.simulation_tasks = self.set_simulation_tasks()

        # FIXME: docu says no
        self.pre_initialize()  # this should be redundant when each experiment is realsed as one single instance of the same class

        print(self)

        self.fit_mappings = self.set_fit_mappings()
        self.plots = self.set_figures()

        self.finalize()

    def pre_initialize(self) -> None:
        """Pre-initialize SimulationExperiment.

        Initialization must be separated from object construction due to
        the parallel execution of the problem later on.
        Certain objects cannot be serialized and must be initialized.
        :return:
        """
        # process all information necessary to run the simulations, i.e.,
        # all data required from the model
        self._datasets = self.dsets  # storage of datasets

        # validation of information
        self._check_keys()
        self._check_types()

    def finalize(self) -> None:
        self.pre_initialize()  # update dsets etc
        self._simulations = self.simulation_tasks  # storage of simulation definition
        self._tasks = {}
        self._fit_mappings = {}  # type: Dict[str, FitMapping]
        self.datagenerators()  # definition of data accessed later on (sets self._data)
        self._fit_mappings = self.fit_mappings

        # validation of information
        self._check_keys()
        self._check_types()

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
            key = str(key)
            # handle '*_unit columns'
            if key.endswith("_unit"):
                # parse the item and unit in dict
                units = dset[key]
                item_key = key[0:-5]  # remove "_unit"
                if item_key not in dset.keys():
                    logger.error(
                        f"Missing * column '{item_key}' for unit " f"column: '{key}'"
                    )
                else:
                    all_udict[item_key] = str(units)
        d = DataSet(dset)
        d.uinfo = UnitsInformation(all_udict, ureg=self.ureg)  # TODO: add all_udict
        d.Q_ = d.uinfo.ureg.Quantity
        return name, d

    # TODO: calculate amount from recovery an replace measurement_type==recovery by cumulative amount
    def add_yid(self):
        """
        Creates the yid for each dset and adds it as a column. This yid can later be used for FitMappings.
        Assumes only one substance and measurement_type per dset.
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

    def set_simulation_tasks(self) -> Dict[str, TimecourseSim]:
        """
            Returns a dictionary with one simulation task for each unique intervention.
        """
        tcsims = {}

        interventions = []
        for intervention in self.tcs["interventions"]:
            if intervention not in interventions:
                interventions.append(intervention)

        min_time, max_time = utils.min_max_time(self)

        for intervention in interventions:
            task_id = ""                                     # initialize new task name
            changes = {}
            for medication in intervention:
                task_id += self.get_intervention_key(medication)
                substance, route, dose = utils.intervention_details(self, medication)
                changes[f"{route}DOSE_{task_mapping[substance]}"] = dose

            task_dict = {
                "start": min_time,
                "end": max_time * 1.1,  # +10% simulation time
                "steps": int(60 * (max_time * 1.1 - min_time) * STEPS_PER_MIN),
                "changes": changes,
                # "model_changes": Dict[str, Quantity],
                # "model_manipulations": dict ,
                # "discard": bool,
            }
            tcsims[task_id] = TimecourseSim(timecourses=Timecourse(**task_dict))
            utils.append_task_id(self, intervention, task_id)
        return tcsims

    def set_fit_mappings(self) -> Dict[str, FitMapping]:
        """
            Iterates over all outputs and returns a dictionary with one entry for each timecourse/output.
        """
        mappings = {}
        for key, dset in self.dsets.items():
            meta_data_dict = utils.fill_meta_dict(dset)
            measured_yid = utils.reporting_type(dset)

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
                    task=f"task_{dset['task']}",
                    xid="time",
                    yid=dset['yid']
                ),
                metadata=DexMappingMetaData(**meta_data_dict)
            )
        return mappings

    def set_figures(self) -> Dict[str, Figure]:

        fig = Figure(experiment=self, sid=self.sid, name=f"{self.sid}", num_rows=4)
        plots = fig.create_plots(xaxis=Axis("time", unit=self.unit_time), legend=True)

        # simulation data
        for task in self.simulation_tasks:
            for substance_index, substance in enumerate(["dex", "dor"]):
                for tissue_index, tissue in enumerate([f"[Cve_{substance}]", f"Aurine_{substance}"]):
                    index = substance_index * 2 + tissue_index
                    plots[index].add_data(task=task, xid="time", yid=tissue,
                                          label=tissue, color="black")


        # experimental data
        substance_mapping = {
            "[Cve_dex]": 0,
            "[Cve_dor]": 1,
            "Aurine_dex": 2,
            "Aurine_dor": 3,
        }

        # TODO: create dset to plot index mapping
        for key, dset in self.dsets.items():
            assert len(dset['yid'].unique()) == 1

            # TODO: get from utils
            if "mean" in dset.columns:
                subject_type = "mean"
            else:
                subject_type = "value"

            if "cyp2d6 phenotype" in dset.columns:
                color = colors_phenotype[dset['cyp2d6 phenotype']][0]
            else:
                color = "tab:grey"

            plots[substance_mapping[dset['yid'][0]]].add_data(dataset=key, xid="time", yid=subject_type,
                              count="count", label="PM", color=color)

        return {
            "fig1": fig,
        }



