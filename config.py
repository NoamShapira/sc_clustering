from pathlib import PurePosixPath

# labs paths
UMI_PATH = PurePosixPath("/home/labs/amit/eyald/sc_pipeline/scdb_v4_mouse/output/umi.tab/")
META_DATA_PATH = PurePosixPath('/home/labs/amit/weiner/Serono/Serono11/serono11_metadata.txt')

# personal paths
RESULTS_DIR = PurePosixPath('/home/labs/amit/noamsh/clustering_results')
META_CELL_PATH = PurePosixPath('/home/labs/amit/noamsh/mc.csv')

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

# transform
TL_PCA_SVD_SOLVER = 'arpack'
TL_PCA_N_COMPS = 100

NEIGHBORHOOD_GRAPH_N_PCS = 100
NEIGHBORHOOD_GRAPH_N_NEIGHBORS = 10

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
