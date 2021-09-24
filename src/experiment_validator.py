from pydantic import BaseModel, validator, ValidationError
from experiment_wrapper import CustomSimulationExperiment


class SimulationExperimentValidator(BaseModel):
    simulation_experiment: CustomSimulationExperiment

    class Config:
        arbitrary_types_allowed = True

    @validator('simulation_experiment')
    def pkdata_is_not_empty(cls, v):
        if len(v.pkdata().study_sids) == 0:
            raise ValidationError('PKdata does not contain your study.')
        return v
