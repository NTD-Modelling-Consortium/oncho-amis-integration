from epioncho_ibm.endgame_simulation import EndgameSimulation
from epioncho_ibm.state.params import EpionchoEndgameModel

from numpy.typing import NDArray
import numpy as np
import pandas as pd
import os 
from pathlib import Path
from multiprocessing import cpu_count
from tqdm.contrib.concurrent import process_map

seeds = int
gamma_distribution = float
bite_rate_per_person_per_year = float

prev = NDArray[np.float_]

        
id = os.getenv("SLURM_ARRAY_TASK_ID")

# Read in csv's
PATH_TO_MODEL_OUTPUT = Path(os.getenv("PATH_TO_MODEL_OUTPUT"))
mda_path = PATH_TO_MODEL_OUTPUT / f'InputMDA_MTP_{id}.csv'
mda_history = pd.read_csv(mda_path)
mda_history = mda_history.sort_values(by=['Year'])
mda_history = mda_history.reset_index(drop=True)

vc_path = PATH_TO_MODEL_OUTPUT / f'InputVC_MTP_{id}.csv'
vc_history = pd.read_csv(vc_path)
vc_history = vc_history.sort_values(by=['Year'])
vc_history = vc_history.reset_index(drop=True)

# Scenario for mda history
treatment_program = []
if mda_history.shape[0] > 0:
    for i in range((mda_history.shape[0])):
        treatment_program.append({
                "first_year": mda_history.Year[i],
                "last_year":  mda_history.Year[i],
                "interventions": {
                    "treatment_interval": mda_history.treatment_interval[i],
                    "total_population_coverage": mda_history.ModelledCoverage[i],
                    "correlation": mda_history.adherence_par[i],
                }
        })


def run_sim(params):
    changes_params = []
    changes_params.append({
        "year": 1970,
        "params": {
            "delta_time_days": 1
        }
    })

    if vc_history.shape[0] > 0:
        prev_bite_rate = None
        for i in range((vc_history.shape[0])):
          bite_rate = vc_history.abr_multiplier[i] * params[0][2]
          if prev_bite_rate is None or bite_rate != prev_bite_rate:
            changes_params.append({
                "year": vc_history.Year[i],
                "params":{
                    "blackfly": {
                        "bite_rate_per_person_per_year": bite_rate
                    }                    
                }
            })
            prev_bite_rate = bite_rate
             
    endgame_structure = {
        "parameters": {
            "initial": {
                "n_people": 400, 
                "seed": params[0][0], 
                "delta_time_days": 7,
                "gamma_distribution": params[0][1],
                "sequela_active": [
                    "HangingGroin",
                    "Atrophy",
                    "Blindness",
                    "APOD",
                    "CPOD",
                    "RSD",
                    "Depigmentation",
                    "SevereItching"
                ],                
                "blackfly": {
                    "bite_rate_per_person_per_year": params[0][2]
                }
            },
            "changes": changes_params
        },
        "programs": treatment_program
    }

    endgame_obj = EpionchoEndgameModel.parse_obj(endgame_structure)
    endgame_sim = EndgameSimulation(
        start_time=1895, endgame=endgame_obj, verbose=False, debug=True
    )

    prev = []
    save_trajectory = []
    for state in endgame_sim.iter_run(end_time=2019,sampling_interval=1):
       
        if state.current_time in [1975,2000,2018]: #timepoints saved should correspond with year of map samples
            prev.append((state.mf_prevalence_in_population()))
        
        if state.current_time in range(1975,2019,1): # for plotting the simulated trajectories
            save_trajectory.append((state.mf_prevalence_in_population()))

    return prev, save_trajectory

def wrapped_parameters(
    parameters: tuple[list[seeds],list[gamma_distribution], list[bite_rate_per_person_per_year]]
):
    restructured_params: list[tuple[seeds, gamma_distribution, bite_rate_per_person_per_year]] = list(zip(parameters)) #type:ignore

    if __name__ == "r_wrapper_endgame_fitting_multipletimepts":
        list_of_stats = process_map(
            run_sim, restructured_params, max_workers=cpu_count()
        )
        #print(np.mean(list_of_stats))
        #print(np.std(list_of_stats))

        return list_of_stats
