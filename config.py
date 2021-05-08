from pathlib import PurePosixPath

# labs paths
UMI_PATH = PurePosixPath("/home/labs/amit/eyald/sc_pipeline/scdb_v4_mouse/output/umi.tab/")
META_DATA_PATH = PurePosixPath('/home/labs/amit/weiner/Serono/Serono11/serono11_metadata.txt')

# personal paths
RESULTS_DIR = PurePosixPath('/home/labs/amit/noamsh/clustering_results')

# pre_proccising
PP_MIN_NUMBER_OF_CELLS_PER_GENE = 3
PP_MIN_NUMBER_OF_GENES_PER_CELL = 200
PP_MAX_NUMBER_OF_GENES_PER_CELL = 2500
PP_MAX_NUMBER_OF_MT_GENES_PER_CELL = 5

PP_NORMELIZED_READS_PER_CELL = 1e4

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
