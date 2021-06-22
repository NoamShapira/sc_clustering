import os
from datetime import datetime
from pathlib import Path

import pytz

import config


def get_now_timestemp_as_string() -> str:
    return datetime.now(tz=pytz.timezone("Israel")).strftime("%Y_%m_%d__%H_%M_%S")


def create_experiment_dir_and_return_path(experiment_name, main_results_dir: Path = config.RESULTS_DIR) -> Path:
    experiment_name = f"{experiment_name}_{get_now_timestemp_as_string()}"
    exp_results_dir_path = Path(main_results_dir, experiment_name)
    os.mkdir(exp_results_dir_path)
    return exp_results_dir_path