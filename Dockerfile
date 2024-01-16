# DOCKER FOR LUMI DIFFDOCK
FROM koubic/lumi-imgs:pyg

RUN pip install --upgrade pip
RUN apt-get -y update
RUN apt-get -y install git

RUN pip install torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric==2.0.4 -f https://data.pyg.org/whl/torch-1.11.0+cu117.html
RUN python3 -m pip install PyYAML scipy "networkx[default]" biopython rdkit-pypi e3nn spyrmsd pandas biopandas

RUN pip install "fair-esm"
RUN pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
