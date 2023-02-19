
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-devel:latest@sha256:e51f2dfe9d97e890006c1517d141b1f4b13855413ed093a599339894b584b1cc

WORKDIR $HOME

USER $FN_USER

RUN $HOME/sire.app/bin/sire_test

ENTRYPOINT ["bash"]
