from pathlib import Path
from typing import Dict, Optional
from dataclasses import dataclass
from pkdb_analysis.filter import pk_info
from pkdb_analysis.meta_analysis import MetaAnalysis
from sbmlsim.data import DataSet
from sbmlsim.experiment import SimulationExperiment
from pkdb_analysis.data import PKData
from sbmlsim.simulation import TimecourseSim, Timecourse


@dataclass
class Substance:
    """Substance details"""
    name: str
    mass: float
    bioavialability: float


def add_genotype(df):
    categorical_field = "cyp2d6 genotype"
    return pk_info(df.extra, categorical_field, ["choice"]).rename(
        columns={f"choice_{categorical_field}": categorical_field}
    )[categorical_field]


def add_phenotype(df):
    categorical_field = "cyp2d6 phenotype"
    return pk_info(df.extra, categorical_field, ["choice"]).rename(
        columns={f"choice_{categorical_field}": categorical_field}
    )[categorical_field]


class CustomSimulationExperiment(SimulationExperiment):

    def __init__(self, pkdata_path: Path, **kwargs) -> None:
        super().__init__(**kwargs)
        self.substance_details: Dict[str: Substance] = None
        self.tasks_pkdata: Optional[Dict] = None
        self.pkdata_path = pkdata_path


    def pkdata(self):
        pkdata = PKData.from_archive(self.pkdata_path).filter({"outputs": {"study_name": self.__class__.__name__}})
        pkdata.scatters = pkdata.scatters._emptify()
        return pkdata



    def datasets_pkdata(self) -> Dict[str, DataSet]:
        """ creates datasets from pkdata object."""
        datasets = {}
        pkdata = self.pkdata()
        if not self.substance_details:
            substance_details = {}
        else:
            substance_details = self.substance_details
        for label, df in pkdata.timecourses.df.groupby("label"):

            this_pkdata = pkdata.filter(
                {
                "timecourses": {"label": label},
                "outputs": {"measurement_type": "concentration"}})
            print(this_pkdata)


            meta_analysis = MetaAnalysis(this_pkdata, intervention_substances=set(substance_details.keys()))
            meta_analysis.create_results()
            results = meta_analysis.results
            for key, add_func in {"genotype": add_genotype, "phenotype": add_phenotype}.items():
                this_info = results.apply(add_func, axis=1)
                if not this_info.empty:
                    results[key] = this_info
                else:
                    results[key] = None

            print(results.columns)


            #results["time"] = results["time"].astype(float)


            time_delta = results.time.max() - results.time.min()
            this_task = {
                "start": results.time.min()*60,
                "end":  results.time.max()*60,
                "steps": int(time_delta),
            }
            results_interventions = results[results["intervention_substance"] == "dmt"]
            single_result = results.iloc[0]

            assert single_result["intervention_number"] == 1, ("Can not be considered currently.", single_result)
            changes = {}
            if single_result.intervention_route == "oral":
                Q_ = self.Q_
                changes["PODOSE_dex"] = Q_(single_result["intervention_value"], single_result["intervention_unit"])
            this_task["changes"] = changes
            self.tasks_pkdata[label] = this_task
            results = results.rename({"group_count": "count"})
            dset = DataSet.from_df(results, self.ureg)

            #print(df.substance)
            if df.substance.unique()[0] == "dmt":
                dset.unit_conversion("value", 1 / self.Mr.dex)
            elif df.substance.unique()[0] == "dor":
                dset.unit_conversion("value", 1 / self.Mr.dor)
            datasets[label] = dset
        return datasets

    def simulation_pkdata(self) -> Dict[str, TimecourseSim]:
        tcsims = {}
        for id,task_dict in self.tasks_pkdata.items():
            tcsims[id] = TimecourseSim(timecourses=Timecourse(**task_dict))
        return tcsims

