FROM ubuntu:25.04@sha256:57665ab8178042ef197191fd77d21d8a2f7f535acd26ff7bd548b1f340f081d7

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        git \
        python3 \
        python3-pip \
        build-essential

# Clone the pythonscad repository
RUN git clone --recursive --branch "add-docker-files" https://github.com/pythonscad/pythonscad.git /opt/pythonscad

WORKDIR /opt/pythonscad

# Install Python dependencies if requirements.txt exists
RUN scripts/uni-get-dependencies.sh

RUN cmake -B build

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest --output-on-failure -j$(nproc)

# Default command
CMD ["build/pythonscad"]
