FROM nvcr.io/nvidia/pytorch:21.05-py3

RUN pip install -r repos/sc_clusterinf/requirments
RUN export PYTHONPATH="${PYTHONPATH}:/home/labs/amit/noamsh/repos/sc_clustering"