FROM archlinux:latest

COPY ../../ /opt/pythonscad

RUN rm -rf /opt/pythonscad/build

WORKDIR /opt/pythonscad

RUN scripts/uni-get-dependencies.sh

RUN cmake -B build

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest -j$(nproc) --output-on-failure

# Default command
CMD ["build/pythonscad"]
