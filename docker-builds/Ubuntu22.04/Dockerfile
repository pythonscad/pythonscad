FROM ubuntu:22.04

COPY ../../ /opt/pythonscad

RUN rm -rf /opt/pythonscad/build

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

RUN cmake -B build -DENABLE_LIBFIVE=ON

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest --output-on-failure -j$(nproc)

# Default command
CMD ["build/pythonscad"]
