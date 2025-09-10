FROM archlinux:latest

COPY ../../ /opt/pythonscad

RUN rm -rf /opt/pythonscad/build

RUN pacman -Syu --noconfirm --needed \
    base-devel \
    boost \
    cairo \
    cgal \
    cmake \
    double-conversion \
    double-conversion \
    eigen \
    fontconfig \
    fontconfig \
    freetype2 \
    gcc-libs \
    ghostscript \
    glew \
    glib2 \
    glibc \
    glu \
    gmp \
    harfbuzz \
    harfbuzz \
    hicolor-icon-theme \
    hidapi \
    imagemagick \
    lib3mf \
    libglvnd \
    libspnav \
    libx11 \
    libxml2 \
    libzip \
    mimalloc \
    mpfr \
    nettle \
    opencsg \
    python \
    python-pip \
    python-setuptools \
    procps \
    qscintilla-qt5 \
    qt5-base \
    qt5-multimedia \
    qt5-svg \
    tbb \
    xorg-server-xvfb

WORKDIR /opt/pythonscad

RUN cmake -B build -DENABLE_LIBFIVE=ON

RUN cmake --build build -j$(nproc)

WORKDIR /opt/pythonscad/build

RUN ctest -j$(nproc) --output-on-failure

# Default command
CMD ["build/pythonscad"]
