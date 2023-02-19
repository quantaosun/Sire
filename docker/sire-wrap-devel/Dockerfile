# This image is used to create a container for building Sire Python wrappers.

FROM siremol/sire-devel:latest@sha256:e51f2dfe9d97e890006c1517d141b1f4b13855413ed093a599339894b584b1cc

# Configure environment
ENV SHELL=/bin/bash \
    FN_USER=sireuser \
    SIRE_SILENT_PHONEHOME=1 \
    SIRE_DONT_PHONEHOME=1 \
    HOME=/home/$FN_USER \
    PATH=$HOME/sire.app/bin:$PATH \
    LD_LIBRARY_PATH=$HOME/sire.app/lib:$HOME/sire.app/lib64:$LD_LIBRARY_PATH

# Switch to the non-root user.
USER $FN_USER

# Install pyplusplus and pygccxml.
RUN pip install pyplusplus pygccxml==1.8.5 fuzzywuzzy && \
    conda install -y -c conda-forge clang clangdev llvmdev gcc_linux-64 gxx_linux-64

# Download and install CastXML.
RUN cd $HOME && \
    git clone https://github.com/CastXML/CastXML.git && \
    cd CastXML && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/sire.app .. && \
    make -j 4 && \
    make install

# Clean up so that the Docker image size is minimised.
RUN conda clean -tipy

# Set the entry directory.
WORKDIR $HOME/Sire/wrapper

ENTRYPOINT ["bash"]
