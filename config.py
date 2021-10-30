from pathlib import PurePosixPath

# labs paths
from data.serono_data_loading_description import SeranoDataLoaderDescription

AmitLab_Path = PurePosixPath("/home/labs/amit")
UMI_DIR_PATH = PurePosixPath(AmitLab_Path, "eyald/sc_pipeline/scdb_v4_mouse/output/umi.tab/")
ASAF_SERANO_DIR_PATH = PurePosixPath(AmitLab_Path, 'weiner/Serono')
ASAF_SERANO_META_DATA_DIR_PATH = PurePosixPath(ASAF_SERANO_DIR_PATH, 'Serono11')
ASAF_META_CELL_ATLAS_DIR_PATH = PurePosixPath(
    ASAF_SERANO_DIR_PATH, "Serono14/clustering_results/simple_clustering_2021_09_05__15_04_43")
META_DATA_PATH = PurePosixPath(ASAF_SERANO_META_DATA_DIR_PATH, 'serono11_metadata.txt')
NEW_META_DATA_PATH = PurePosixPath(ASAF_SERANO_META_DATA_DIR_PATH, 'SeronoMetadata13.txt')
UPDATED_META_DATA_PATH = PurePosixPath(AmitLab_Path, 'noamsh/data/SeronoMetadata_v2_2021_08_28.xlsx')

# personal paths
NOAM_DIR_PATH = PurePosixPath(AmitLab_Path, 'noamsh')
RESULTS_DIR = PurePosixPath(NOAM_DIR_PATH, 'clustering_results')
ANNOTATION_DIR = PurePosixPath(NOAM_DIR_PATH, 'clustering_annotation')
META_CELL_PATH = PurePosixPath(NOAM_DIR_PATH, 'mc.csv')


# data loading
SERANO_DATA_LOADING_DESCRIPTION :SeranoDataLoaderDescription = SeranoDataLoaderDescription.ARMS_1_2_3_FROM_NOAMSH

# env consts
IO_N_WORKERS = 32

# pre_proccising
PP_MIN_NUMBER_OF_CELLS_PER_GENE = 20
PP_MIN_NUMBER_OF_GENES_PER_CELL = 200
PP_MAX_NUMBER_OF_GENES_PER_CELL = 7000
PP_MAX_PCT_OF_MT_GENES_PER_CELL = 20

PP_NORMELIZED_READS_PER_CELL = 1e4

PP_REGRESS_OUT_OBS_KEYS = ['total_counts', 'pct_counts_mt']
PP_SCALE_MAX_VALUE = 10

PP_DROP_UNVARIABLE_GENES = True
PP_REGRESS_OUT_TOTAL_CONT_AND_MT = False


# transform
TL_PCA_SVD_SOLVER = 'arpack'
TL_PCA_N_COMPS = 100

NEIGHBORHOOD_GRAPH_N_PCS = 100
NEIGHBORHOOD_GRAPH_N_NEIGHBORS = 15

TL_LOUVAIN_RESOLUTION = 1.0
TL_LEIDEN_RESOLUTION = 1.0


PP_HIGHLY_VARIATION_PARAMS = {
    "min_mean": 0.0125,
    "max_mean": 3,
    "min_disp": 0.5
}
PP_QUALITY_CONTROL_PARAMS = {
    "qc_vars": ['mt'],
    "percent_top": None,
    "log1p": False,
    "inplace": True
}

# run config
DEBUG_MODE = False
DEBUG_N_BATCHES = 10
CLUSTERING_METHOD = "leiden"