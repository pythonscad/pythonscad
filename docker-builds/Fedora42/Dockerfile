FROM fedora:42

COPY ../../ /opt/pythonscad

RUN rm -rf /opt/pythonscad/build

WORKDIR /opt/pythonscad

RUN scripts/uni-get-dependencies.sh

RUN dnf install -y mimalloc-devel python3-devel qt5-qtgamepad-devel \
    hidapi-devel libspnav-devel curl lib3mf-devel guile22-devel guile30-devel \
    python3-setuptools

RUN cmake -B build -DENABLE_LIBFIVE=ON

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest --output-on-failure -j$(nproc)

# Default command
CMD ["build/pythonscad"]
