FROM ubuntu:22.04@sha256:d0afa9fbcf16134b776fbba4a04c31d476eece2d080c66c887fdd2608e4219a9

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        git \
        python3 \
        python3-pip \
        build-essential \
        wget \
        sudo

# Add Lib3MF repo
RUN wget -qO - https://files.openscad.org/OBS-Repository-Key.pub | sudo apt-key add - && \
  echo yes | sudo add-apt-repository "deb https://download.opensuse.org/repositories/home:/t-paul:/lib3mf/xUbuntu_22.04/ ./" && \
  echo yes | sudo add-apt-repository "deb https://download.opensuse.org/repositories/home:/t-paul:/cgal/xUbuntu_22.04/ ./" && \
  sudo apt-get update -qq && \
  sudo apt-get purge -qq fglrx || true

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
