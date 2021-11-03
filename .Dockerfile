FROM deeplearnphysics/larcv2:minimal

LABEL maintainer="drinkingkazu"
LABEL contact="contact@deeplearnphysics.org"
LABEL type="gpu"
LABEL version="ub18.04-cuda10.2-pytorch1.7.1-edepsim"

ARG DEBIAN_FRONTEND=noninteractive

# root
ENV ROOTSYS=/app/root
ENV PATH="${ROOTSYS}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${ROOTSYS}/lib:${LD_LIBRARY_PATH}"
ENV PYTHONPATH="${ROOTSYS}/lib:${PYTHONPATH}"


# geant4
ENV PATH="/app/geant4/bin:${PATH}"
ENV LD_LIBRARY_PATH="/app/geant4/lib:${LD_LIBRARY_PATH}"

#edepsim
ENV PATH="/app/edep/bin:${PATH}"
ENV LD_LIBRARY_PATH="/app/edep/lib:${LD_LIBRARY_PATH}"

# for geant4
ENV G4NEUTRONHPDATA="/app/geant4/share/Geant4-10.6.2/data/G4NDL4.6"
ENV G4LEDATA="/app/geant4/share/Geant4-10.6.2/data/G4EMLOW7.9.1"
ENV G4LEVELGAMMADATA="/app/geant4/share/Geant4-10.6.2/data/PhotonEvaporation5.5"
ENV G4RADIOACTIVEDATA="/app/geant4/share/Geant4-10.6.2/data/RadioactiveDecay5.4"
ENV G4PARTICLEXSDATA="/app/geant4/share/Geant4-10.6.2/data/G4PARTICLEXS2.1"
ENV G4PIIDATA="/app/geant4/share/Geant4-10.6.2/data/G4PII1.3"
ENV G4REALSURFACEDATA="/app/geant4/share/Geant4-10.6.2/data/RealSurface2.1.1"
ENV G4SAIDXSDATA="/app/geant4/share/Geant4-10.6.2/data/G4SAIDDATA2.0"
ENV G4ABLADATA="/app/geant4/share/Geant4-10.6.2/data/G4ABLA3.1"
ENV G4INCLDATA="/app/geant4/share/Geant4-10.6.2/data/G4INCL1.0"
ENV G4ENSDFSTATEDATA="/app/geant4/share/Geant4-10.6.2/data/G4ENSDFSTATE2.2"

# Geant4
RUN apt-get -y install libxerces-c-dev && \
    wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1vvPiZgeR2l_gtlaVNGySXhL3Z5zo2h52' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1vvPiZgeR2l_gtlaVNGySXhL3Z5zo2h52" -O geant4.tar.gz && \
    rm /tmp/cookies.txt && \
    tar -xzf geant4.tar.gz && \
    rm geant4.tar.gz && \
    mkdir -p /app && mv geant4-docker /app/geant4


LABEL org.opencontainers.image.authors="kterao@slac.stanford.edu"

# larcv build
ENV LARCV_BASEDIR=/app/larcv3
ENV LARCV_BUILDDIR="${LARCV_BASEDIR}/build"
#ENV LARCV_COREDIR="${LARCV_BASEDIR}/larcv/core"
#ENV LARCV_APPDIR="${LARCV_BASEDIR}/larcv/app"
ENV LARCV_LIBDIR="${LARCV_BUILDDIR}/lib"
#ENV LARCV_INCDIR="${LARCV_BUILDDIR}/include"
ENV LARCV_BINDIR="${LARCV_BUILDDIR}/bin"
#ENV LARCV_ROOT6=1
#ENV LARCV_CXX=g++

# without numpy
#ENV LARCV_NUMPY=0
#ENV LARCV_INCLUDES="-I${LARCV_INCDIR} "
#ENV LARCV_LIBS="-L${LARCV_LIBDIR} -llarcv"

# with numpy
#ENV LARCV_NUMPY=1
#ENV LARCV_INCLUDES="-I${LARCV_INCDIR} -I/usr/include/python2.7 -I/usr/include/x86_64-linux-gnu/python2.7"
#ENV LARCV_LIBS="-L/usr/lib/ -L/usr/lib/python2.7/config-x86_64-linux-gnu -L/usr/lib"
#ENV LARCV_LIBS="${LARCV_LIBS} -lpthread -ldl -lutil -lm -lpython2.7 -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions"
#ENV LARCV_LIBS="${LARCV_LIBS} -L${LARCV_LIBDIR} -llarcv"

# set bin and lib path
ENV PATH=${LARCV_BASEDIR}/bin:${LARCV_BINDIR}:${PATH}
ENV LD_LIBRARY_PATH=${LARCV_LIBDIR}:${LD_LIBRARY_PATH}:
ENV PYTHONPATH=${LARCV_BASEDIR}/python3:${PYTHONPATH}


# larcv
#RUN mkdir -p /app && \
#    cd /app && \
#    git clone https://github.com/DeepLearnPhysics/larcv2 && \
#    cd larcv2 && \
#    mkdir -p $LARCV_BUILDDIR && \
#    mkdir -p $LARCV_LIBDIR && \
#    mkdir -p $LARCV_BINDIR && \
#    make -j4

#python
RUN apt-get update -y
RUN apt-get install -y python3 python-pip libhdf5-dev libopenmpi-dev
RUN pip3 install --upgrade pip
RUN pip3 install cmake
RUN pip3 install h5py
RUN pip3 install scikit-build
RUN pip3 install pytest
RUN pip3 install mpi4py
RUN pip3 install numpy scipy matplotlib ipython jupyter pandas sympy nose

# ROOT
RUN apt-get -y install libx11-dev libxpm-dev libxft-dev libxext-dev libssl-dev
#RUN wget https://root.cern/download/root_v6.22.02.Linux-ubuntu18-x86_64-gcc7.5.tar.gz && \
#    tar -xzf root_v6.22.02.Linux-ubuntu18-x86_64-gcc7.5.tar.gz && \
#    rm root_v6.22.02.Linux-ubuntu18-x86_64-gcc7.5.tar.gz && \
#    mkdir -p /app && mv root /app/root
RUN wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1tZpuZTWN0sY8QA-OG2ACZAYLukHDL-Lb' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1tZpuZTWN0sY8QA-OG2ACZAYLukHDL-Lb" -O root.tar.gz && \
    rm /tmp/cookies.txt && \
    tar -xzf root.tar.gz && \
    rm root.tar.gz && \
    mkdir -p /app && mv root /app/root
RUN pip3 --no-cache-dir install rootpy root_numpy



# larcv3
RUN mkdir -p /app && \
    cd /app && \
    git clone --branch IO_separation https://github.com/zhulcher/larcv3.git && \
    cd larcv3 && \
    mkdir -p $LARCV_BUILDDIR && \
    mkdir -p $LARCV_LIBDIR && \
    mkdir -p $LARCV_BINDIR && \
    git submodule update --init &&\
    python3 setup.py build -j 2 && \
    python3 setup.py install

# edep-sim
#RUN apt-get -y install libxerces-c-dev && \
#    wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1ObLd_Z63jcCmfj4iUGB0n0nKO2kE00nQ' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1ObLd_Z63jcCmfj4iUGB0n0nKO2kE00nQ" -O edep.tar.gz && \
#    rm /tmp/cookies.txt && \
#    tar -xzf edep.tar.gz && \
#    rm edep.tar.gz && \
#    mkdir -p /app && mv edep /app/edep
RUN mkdir -p edep-build /app/edep /source && \
    git clone https://github.com/ClarkMcGrew/edep-sim && mv edep-sim /source/edep && \
    cd edep-build && \
    cmake -DCMAKE_INSTALL_PREFIX=/app/edep -DEDEPSIM_DISPLAY=OFF /source/edep && \
    make && make install && \
    cd .. && rm -r edep-build

# remove temporary cache
RUN apt-get autoremove -y && apt-get clean -y