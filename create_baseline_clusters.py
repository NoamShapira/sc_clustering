import os
from pathlib import Path

import config
import utils
from data.data_loading import get_all_experiments_in_one_anndata


def create_experiment_dir_and_return_path(experiment_name, main_results_dir: Path = config.RESULTS_DIR) -> Path:
    experiment_name = f"{experiment_name}_{utils.get_now_timestemp_as_string()}"
    exp_results_dir_path = Path(main_results_dir, experiment_name)
    os.mkdir(exp_results_dir_path)
    return exp_results_dir_path


if __name__ == '__main__':
    experiment_results_dir_path = create_experiment_dir_and_return_path("simple_clustering")

    adata = get_all_experiments_in_one_anndata()
    # create_clusters(adata, experiment_results_dir_path)
    adata.write(experiment_results_dir_path)
