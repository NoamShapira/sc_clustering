from datetime import datetime
from glob import glob
from os import listdir
from pathlib import Path

import anndata as ad
import numpy as np

import config
from utils import TIME_PATTERN_LEN, TIME_PATTERN


def get_name_of_last_experiment() -> str:
    all_experiment_names = listdir(config.RESULTS_DIR)
    all_experimnt_times = [datetime.strptime(exp_name[-TIME_PATTERN_LEN:], TIME_PATTERN) for exp_name in
                           all_experiment_names]
    return all_experiment_names[np.argmax(all_experimnt_times)]


def get_last_result_from_experiment(experiment_name: str) -> ad.AnnData:
    if experiment_name == "last":
        experiment_name = get_name_of_last_experiment()

    experiment_path = Path(config.RESULTS_DIR, experiment_name)
    all_results_paths = [Path(path) for path in glob(str(experiment_path) + "/final*")]
    all_results_times = [datetime.strptime(result_path.name.split(".")[0][-TIME_PATTERN_LEN:], TIME_PATTERN)
                         for result_path in all_results_paths]
    return ad.read(all_results_paths[np.argmax(all_results_times)])


def get_raw_data_from_experiment(experiment_name: str) -> ad.AnnData:
    if experiment_name == "last":
        experiment_name = get_name_of_last_experiment()

    return ad.read(Path(config.RESULTS_DIR, experiment_name, "loaded_data.h5ad"))
