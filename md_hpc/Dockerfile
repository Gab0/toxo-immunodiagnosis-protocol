FROM nvcr.io/hpc/gromacs:2021.3
#FROM gromacs/gromacs
#FROM nvidia/cuda

RUN uname
#RUN apt-get update
#RUN apt-get install -y wget

# Update nvidia keyring;
ARG DISTRO_ARCH=ubuntu2004/x86_64
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/${DISTRO_ARCH}/3bf863cc.pub
#RUN wget https://developer.download.nvidia.com/compute/cuda/repos/${DISTRO_ARCH}/cuda-keyring_1.0-1_all.deb

RUN apt-get update
RUN apt-get install -y apt-utils python git
RUN apt-get install -y python3-pip
RUN apt-get install -y pymol 

# Install these autoGromacs dependencies first
# for faster image rebuilding.
RUN pip3 install numpy==1.21.3
RUN pip3 install mdanalysis freesasa
RUN pip3 install correlationplus seaborn umap-learn


# Copy FORCE-FIELD files.
ARG FF_DIR=/usr/local/gromacs/avx2_256/share/gromacs/top
COPY FF/amber14sb_parmbsc1_cufix.ff ${FF_DIR}/amber14sb_parmbsc1_cufix.ff
COPY FF/charmm36_ljpme-jul2021.ff ${FF_DIR}/charmm36_ljpme-jul2021.ff
COPY FF/charmm36-jul2021.ff ${FF_DIR}/charmm36-jul2021.ff

# SETUP DSSP

# Install build dependencies;
RUN apt-get install -y wget automake 
RUN apt-get install -y libboost-all-dev

# Prepare files and directories;
ARG DSSP_VERSION=3.1.4
RUN wget https://github.com/cmbi/dssp/archive/${DSSP_VERSION}.tar.gz -O dssp.tar.gz
RUN tar -xvf dssp.tar.gz
WORKDIR dssp-${DSSP_VERSION}

# Build;
RUN ./autogen.sh
RUN ./configure --prefix=/usr/
RUN make

# Install;
RUN install -D -m0755 mkdssp /usr/bin
RUN install -D -m0755 test_readpdb /usr/bin/
RUN install -D -m0644 LICENSE_1_0.txt "/usr/share/licenses/dssp/LICENSE"


ARG CACHEBUST1=6
RUN git clone https://github.com/Gab0/epitope-protocol
WORKDIR epitope-protocol
RUN pip3 install -r requirements.txt
RUN pip3 install .



ARG CACHEBUST=9
RUN git clone https://github.com/Gab0/auto-gromacs
WORKDIR auto-gromacs
RUN pip3 install -r requirements.txt
RUN pip3 install .

# Address matplotlib backend dependencies and tests;
RUN apt-get install -y python3-cairo
RUN python3 -c "import matplotlib; matplotlib.use('Gtk3Agg')"

WORKDIR /home
COPY test test
RUN python3 test/test.py

RUN ln -s /usr/local/gromacs/avx2_256/bin/gmx /usr/local/bin/gmx

RUN autogromacs --help; exit 0


CMD ["/bin/bash"]
