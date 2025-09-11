FROM debian:12

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        git \
        python3 \
        python3-pip \
        build-essential

# Clone the pythonscad repository
RUN git clone --recursive https://github.com/pythonscad/pythonscad.git /opt/pythonscad

WORKDIR /opt/pythonscad

RUN scripts/uni-get-dependencies.sh

RUN cmake -B build

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest --output-on-failure -j$(nproc)

# Default command
CMD ["build/pythonscad"]
