FROM linuxmintd/mint22.2-amd64

COPY ../../ /opt/pythonscad

RUN rm -rf /opt/pythonscad/build

WORKDIR /opt/pythonscad

# Install Python dependencies if requirements.txt exists
RUN scripts/uni-get-dependencies.sh

WORKDIR /opt/pythonscad

RUN cmake -B build -DENABLE_LIBFIVE=ON

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest --output-on-failure -j$(nproc)

# Default command
CMD ["build/pythonscad"]
