from pydantic.dataclasses import dataclass
from experiment_wrapper_janosch import PKDataSimulationExperiment


@dataclass
class Experiment:
    sid: str

if __name__ == "__main__":
    Capon1996: Experiment = Experiment("Capon1996")
    print(Capon1996)
    test: PKDataSimulationExperiment = PKDataSimulationExperiment(sid=Capon1996.sid)

    # run_experiments(Capon1996, output_dir="test")
