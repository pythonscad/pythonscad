FROM linuxmintd/mint22.2-amd64@sha256:f0fdbb7c9e8e0e1dcea5c76458215c45f4eda541a3d27731fbb5a6452fba39be

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        git \
        python3 \
        python3-pip \
        build-essential \
        wget \
        sudo

# Clone the pythonscad repository
RUN git clone --recursive --branch "add-docker-files" https://github.com/pythonscad/pythonscad.git /opt/pythonscad

WORKDIR /opt/pythonscad

# Install Python dependencies if requirements.txt exists
RUN scripts/uni-get-dependencies.sh

WORKDIR /opt/pythonscad

RUN cmake -B build

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest --output-on-failure -j$(nproc) --tests-regex "pythonscad_sdf" || true

# Default command
CMD ["build/pythonscad"]
