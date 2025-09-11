FROM debian:12

COPY ../../ /opt/pythonscad

RUN rm -rf /opt/pythonscad/build

# Install dependencies
WORKDIR /opt/pythonscad

RUN scripts/uni-get-dependencies.sh

RUN cmake -B build

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest --output-on-failure -j$(nproc)

# Default command
CMD ["build/pythonscad"]
