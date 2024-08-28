FROM shfiorenza/cylaks_base:latest
WORKDIR /build
RUN git clone --recursive --single-branch --branch angular_potential https://github.com/Betterton-Lab/CyLaKS.git &&\
    cd CyLaKS &&\
    ./install.sh -I
RUN rm -rf /build
WORKDIR /mnt
ENTRYPOINT []
