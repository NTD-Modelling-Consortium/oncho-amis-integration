#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import enum
from pathlib import Path


class Stage(enum.Enum):
    FITTING_PREP = "fitting-prep"
    FITTING = "fitting"
    PROJECTIONS_PREP = "projections-prep"
    NEARTERM_PROJECTIONS = "nearterm-projections"
    ALL = "all"
    SKIP_FITTING_PREP = "skip-fitting-prep"


def validate_environment():
    """Validate required environment variables are set and we're in the correct directory."""
    oncho_amis_dir = os.getenv("ONCHO_AMIS_DIR")
    if not oncho_amis_dir:
        raise ValueError("ONCHO_AMIS_DIR environment variable is not set")

    current_dir = os.getcwd()
    if current_dir != oncho_amis_dir:
        raise ValueError(
            f"This script must be run from '{oncho_amis_dir}', "
            f"but it is running from '{current_dir}'."
        )


def run_command(command, description=None):
    """Run a command and handle any errors."""
    if description:
        print(f"{description}...")

    try:
        result = subprocess.run(command, shell=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}", file=sys.stderr)
        return False


def get_args_for_fitting(args):
    """Prepare arguments for the fitting script based on provided command-line arguments."""
    r_args = []
    if args.id:
        r_args.append(f"--id={args.id}")
    elif args.failed_ids:
        r_args.append(f"--failed-ids={args.failed_ids}")

    if args.amis_sigma:
        r_args.append(f"--amis-sigma={args.amis_sigma}")

    if args.amis_target_ess:
        r_args.append(f"--amis-target-ess={args.amis_target_ess}")

    if args.amis_n_samples:
        r_args.append(f"--amis-n-samples={args.amis_n_samples}")

    if args.amis_max_iters:
        r_args.append(f"--amis-max-iters={args.amis_max_iters}")

    return r_args


def get_args_for_projections_prep(args):
    """Prepare arguments for the projections preparation script based on provided command-line arguments."""
    r_args = []
    if args.id:
        r_args.append(f"--id={args.id}")
    elif args.failed_ids:
        r_args.append(f"--failed-ids={args.failed_ids}")

    if args.amis_sigma:
        r_args.append(f"--amis-sigma={args.amis_sigma}")

    if args.amis_n_samples:
        r_args.append(f"--amis-n-samples={args.amis_n_samples}")

    if args.ess_threshold:
        r_args.append(f"--ess-threshold={args.ess_threshold}")

    return r_args


def get_args_for_nearterm_projections(args):
    """Prepare arguments for the near-term projections script based on provided command-line arguments."""
    r_args = []
    if args.id:
        r_args.append(f"--id={args.id}")

    if args.amis_n_samples:
        r_args.append(f"--amis-n-samples={args.amis_n_samples}")

    return r_args


def do_fitting_prep(args):
    """Run the fitting preparation stage of the pipeline."""
    PATH_TO_FITTING_PREP_SCRIPTS = Path(
        os.getenv("PATH_TO_FITTING_PREP_SCRIPTS", "fitting-prep/scripts")
    )

    if not run_command(
        f"bash {PATH_TO_FITTING_PREP_SCRIPTS / 'run_fitting_prep.sh'}",
        "Running oncho fitting preparation",
    ):
        return False

    return True


def do_fitting(args):
    """Run the fitting stage of the pipeline."""
    PATH_TO_FITTING_SCRIPTS = Path(
        os.getenv("PATH_TO_FITTING_SCRIPTS", "fitting/scripts")
    )
    r_args = get_args_for_fitting(args)

    if not run_command(
        f"bash {PATH_TO_FITTING_SCRIPTS / 'run_fitting.sh'} {' '.join(r_args)}",
        "Running oncho fitting",
    ):
        return False

    return True


def do_projections_prep(args):
    """Run the projections preparation stage of the pipeline."""
    PATH_TO_PROJECTIONS_PREP_SCRIPTS = Path(
        os.getenv("PATH_TO_PROJECTIONS_PREP_SCRIPTS", "projections-prep/scripts")
    )

    r_args = get_args_for_projections_prep(args)

    if not run_command(
        f"bash {PATH_TO_PROJECTIONS_PREP_SCRIPTS / 'run_projections_inputs.sh'} {' '.join(r_args)}",
        "Preparing inputs for near-term projections",
    ):
        return False

    return True


def do_nearterm_projections(args):
    """Run the near-term projections (to 2026) stage of the pipeline."""
    PATH_TO_PROJECTIONS_SCRIPTS = Path(
        os.getenv("PATH_TO_PROJECTIONS_SCRIPTS", "projections/scripts")
    )

    r_args = get_args_for_nearterm_projections(args)

    if not run_command(
        f"bash {PATH_TO_PROJECTIONS_SCRIPTS / 'run_projections_to_2026.sh'} {' '.join(r_args)}",
        "Running historic/near-term simulations",
    ):
        return False

    return True


STAGE_SEQUENCE_MAP = {
    Stage.FITTING_PREP.value: (do_fitting_prep,),
    Stage.FITTING.value: (do_fitting,),
    Stage.PROJECTIONS_PREP.value: (do_projections_prep,),
    Stage.NEARTERM_PROJECTIONS.value: (do_nearterm_projections,),
    Stage.ALL.value: (
        do_fitting_prep,
        do_fitting,
        do_projections_prep,
        do_nearterm_projections,
    ),
    Stage.SKIP_FITTING_PREP.value: (
        do_fitting,
        do_projections_prep,
        do_nearterm_projections,
    ),
}


def main():
    parser = argparse.ArgumentParser(
        description="Run the oncho AMIS pipeline end-to-end"
    )

    # Required arguments
    parser.add_argument(
        "-i",
        "--id",
        type=int,
        required=True,
        help="Batch/task ID to process",
    )

    # Optional arguments
    parser.add_argument(
        "--failed-ids",
        type=str,
        required=False,
        help="Comma-separated list ('id1,id2,id3...') of failed batch/task IDs to skip."
        " Only used when --id is not specified and in the Projections-Prep stage.",
    )

    parser.add_argument(
        "--stage",
        type=str,
        choices=[
            Stage.FITTING_PREP.value,
            Stage.FITTING.value,
            Stage.PROJECTIONS_PREP.value,
            Stage.NEARTERM_PROJECTIONS.value,
            Stage.ALL.value,
            Stage.SKIP_FITTING_PREP.value,
        ],
        required=False,
        default="skip-fitting-prep",
        help="Stage of the pipeline to run. "
        "Options: 'fitting-prep', 'fitting', 'projections-prep', "
        "'nearterm-projections', 'all', 'skip-fitting-prep'. "
        "Default is 'skip-fitting-prep'. "
        "If 'all' is specified, it runs all stages in order. "
        "If 'skip-fitting-prep' is specified, it skips the fitting preparation stage.",
    )

    # AMIS-related parameters
    parser.add_argument(
        "--amis-sigma",
        type=float,
        required=False,
        help="AMIS 'sigma' parameter (default: 0.0025)",
    )
    parser.add_argument(
        "--amis-target-ess",
        type=int,
        required=False,
        help="Target ESS parameter for AMIS (default: 500)",
    )
    parser.add_argument(
        "--amis-n-samples",
        type=int,
        required=False,
        help="Number of AMIS samples (default: 500)",
    )
    parser.add_argument(
        "--amis-max-iters",
        type=int,
        required=False,
        help="Maximum number of AMIS iterations (default: 50)",
    )
    parser.add_argument(
        "--ess-threshold",
        type=int,
        required=False,
        help="ESS threshold parameter (default: 200)",
    )

    args = parser.parse_args()

    try:
        # Validate environment
        validate_environment()

        # Validate arguments
        if args.id:
            print(f"Running pipeline for ID: {args.id}")
            os.environ["SLURM_ARRAY_TASK_ID"] = str(args.id)
            print(f"SLURM_ARRAY_TASK_ID set to {os.environ['SLURM_ARRAY_TASK_ID']}")
            if args.failed_ids:
                print("Warning: --failed-ids is ignored when --id is specified")

        for stage_function in STAGE_SEQUENCE_MAP[args.stage]:
            print(f"Running stage: {stage_function.__name__}")
            if not stage_function(args):
                print(f"Stage {stage_function.__name__} failed.", file=sys.stderr)
                return 1

        print("Pipeline completed successfully!")
        return 0

    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
