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

# Install Python dependencies if requirements.txt exists
RUN scripts/uni-get-dependencies.sh

RUN apt-get install -y libmimalloc-dev python3-dev libqt5gamepad5-dev python3-setuptools

RUN cmake -B build -DENABLE_LIBFIVE=ON

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest --output-on-failure -j$(nproc)

# Default command
CMD ["build/pythonscad"]
