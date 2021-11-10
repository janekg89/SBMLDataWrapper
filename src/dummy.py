

# from pkdb_models.models.dextromethorphan.helpers import run_experiments
from pydantic.dataclasses import dataclass
from experiment_wrapper_janosch import PKDataSimulationExperiment


@dataclass
class Experiment:
    sid: str


if __name__ == "__main__":
    Capon1996: Experiment = Experiment("Capon1996")
    print(Capon1996)
    test: PKDataSimulationExperiment = PKDataSimulationExperiment(sid=Capon1996.sid)
    print(test)
    #run_experiments(test, output_dir="test")
