

# from pkdb_models.models.dextromethorphan.helpers import run_experiments
from pydantic.dataclasses import dataclass
from experiment_wrapper_janosch import PKDataSimulationExperiment
from experiment_factory import ExperimentFactory


@dataclass
class Experiment:
    sid: str


if __name__ == "__main__":
    from pkdb_models.models.dextromethorphan import PKDATA_ZIP_PATH
    print(type(PKDATA_ZIP_PATH))
    Capon1996 = ExperimentFactory(sid="Capon1996", zip_path=PKDATA_ZIP_PATH)
    print(Capon1996)
    # Capon1996: Experiment = Experiment("Capon1996")
    # print(Capon1996)
    # test: PKDataSimulationExperiment = PKDataSimulationExperiment(sid=Capon1996.sid)
    # print(test)
    # run_experiments(test, output_dir="test")
    
