import pandas as pd
import h5py
import os
import time
from pathlib import Path
from multiprocessing import cpu_count
from tqdm.contrib.concurrent import process_map

from epioncho_ibm.endgame_simulation import (
    EndgameSimulation,
)
from epioncho_ibm.state.params import EpionchoEndgameModel
from epioncho_ibm.tools import Data, add_state_to_run_data, convert_data_to_pandas
from endgame_postprocessing.post_processing.single_file_post_processing import process_single_file


PATH_TO_MODEL_OUTPUT = Path(os.getenv("PATH_TO_MODEL_OUTPUT", "./model_output"))

start = time.time()

def get_endgame(seed,exp,abr,treatment_program,vc_history):
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
          bite_rate = vc_history.abr_multiplier[i] * abr
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
    return {
        "parameters": {
            "initial": {
                "n_people": 400, 
                "seed": seed, 
                "delta_time_days": 7,
                "gamma_distribution": exp,
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
                    "bite_rate_per_person_per_year": abr
                }
            },
            "changes": changes_params
        },
        "programs": treatment_program
    }

def run_sim(i, file_path,endgame_structure):
    
    endgame = EpionchoEndgameModel.parse_obj(endgame_structure)
    endgame_sim = EndgameSimulation(
        start_time=1895, endgame=endgame, verbose=False, debug=True
    )
    new_file=h5py.File(file_path,'a')

    run_data: Data = {}
    for state in endgame_sim.iter_run(end_time=2026,sampling_interval=1):
        #print(state.current_time)
        add_state_to_run_data(
            state,
            run_data=run_data,
            prevalence=True,
            intensity=False,
            number=True,
            mean_worm_burden=False,
            prevalence_OAE=True,
            n_treatments=False,
            achieved_coverage=False,
            with_age_groups=False,
            with_sequela=True,
            with_pnc=True,
        )
    
    grp = new_file.create_group(
        f"draw_{str(i)}"
    )
    endgame_sim.save(grp)

    return run_data
   
def post_processing_calculation(
    data: list[Data],
    iuName: str,
    csv_file: str,
    num_draws: int
) -> None:
    pandas_data = convert_data_to_pandas(data)
    df = process_single_file(pandas_data,
                            "model_fits", # change if you want to label the "scenario"
                            iuName,
                            post_processing_start_time=1895,  # change to be the year you start the model
                            post_processing_end_time=2026, # change to be the year you end the model
                            num_draws=num_draws # change to match the number of draws you are using
                            )
    df.to_csv(csv_file)

def process_single_simulation(args):
    idx, file_path, endgame_structure = args
    print(f"Running simulation {idx}")
    run_data = run_sim(idx, file_path, endgame_structure)
    return run_data


def wrapped_parameters(IU):

    started=str(IU) + ' started'
    print(started)

    file_path = PATH_TO_MODEL_OUTPUT / f'OutputVals_MTP_{IU}.hdf5'
    new_file = h5py.File(file_path, 'w')

    # Read in csv's
    mda_path = PATH_TO_MODEL_OUTPUT / f'InputMDA_MTP_proj_{IU}.csv'
    mda_history = pd.read_csv(mda_path)
    mda_history = mda_history.sort_values(by=['Year'])
    mda_history = mda_history.reset_index(drop=True)
    print(mda_history)

    vc_path = PATH_TO_MODEL_OUTPUT / f'InputVC_MTP_proj_{IU}.csv'
    vc_history = pd.read_csv(vc_path)
    vc_history = vc_history.sort_values(by=['Year'])
    vc_history = vc_history.reset_index(drop=True)
    print(vc_history)

    param_path = PATH_TO_MODEL_OUTPUT / f'InputPars_MTP_proj_{IU}.csv'
    params = pd.read_csv(param_path)
    print(params)

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


    endgame_structures = [
        get_endgame(params.loc[p,"seed"], params.loc[p,"individual_exposure"], params.loc[p,"bite_rate_per_person_per_year"],treatment_program, vc_history) for p in range(n_runs) 
    ]

    sim_args = [
        (i, file_path, endgame_structures[i]) for i in range(n_runs)
    ]

    results = process_map(process_single_simulation, sim_args, max_workers=min(cpu_count(), len(sim_args)), chunksize=1)    

    output_data: list[Data] = []
    for run_data in results:
        output_data.append(run_data)

    post_processing_calculation(output_data, str(IU), PATH_TO_MODEL_OUTPUT / f"model_output_MTP_{IU}.csv", n_runs)

    finished = str(IU) + ' completed'
    print(finished)
    
n_runs = 200

id = os.getenv("SLURM_ARRAY_TASK_ID")
IU_list = pd.read_csv(PATH_TO_MODEL_OUTPUT / f'IUs_MTP_proj_{id}.csv', header=None)
IU_list.columns = ["IU"]

# Run simulations 
if __name__ == "__main__":
    for iu in IU_list.IU:
        wrapped_parameters(iu)

end = time.time()
elapsed_time = end-start
print("elapsed time in seconds: " + str(elapsed_time) )

# # for testing on local computer                                                 
# #iu = "GHA0216121382"
# #wrapped_parameters(iu)

# id = "25"
# path = './model_output/IUs_MTP_proj_' + id + '.csv'
# IU_list = pd.read_csv(path, header=None)
# IU_list.columns = ["IU"]

# # Run simulations 
# if __name__ == "__main__":
#     list_of_simulations = process_map(
#         wrapped_parameters, IU_list.IU, max_workers=cpu_count()
#     )

# #for IU in IU_list.IU:
# #    wrapped_parameters(IU)

# end = time.time()
# elapsed_time = end-start
# print("elapsed time in seconds: " + str(elapsed_time) )