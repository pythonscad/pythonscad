FROM ubuntu:22.04

COPY ../../ /opt/pythonscad

RUN rm -rf /opt/pythonscad/build

# Add Lib3MF repo
RUN apt update && apt install -y curl gpg software-properties-common && \
  curl https://files.openscad.org/OBS-Repository-Key.pub | gpg --import && \
  echo yes | add-apt-repository "deb https://download.opensuse.org/repositories/home:/t-paul:/lib3mf/xUbuntu_22.04/ ./" && \
  echo yes | add-apt-repository "deb https://download.opensuse.org/repositories/home:/t-paul:/cgal/xUbuntu_22.04/ ./" && \

WORKDIR /opt/pythonscad

# Install Python dependencies if requirements.txt exists
RUN scripts/uni-get-dependencies.sh

RUN cmake -B build -DENABLE_LIBFIVE=ON

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest --output-on-failure -j$(nproc)

# Default command
CMD ["build/pythonscad"]
