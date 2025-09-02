#! /bin/sh
export PYTHON_VERSION=$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:2])))")
cd build
cmake  -DCMAKE_INSTALL_PREFIX=/usr -DEXPERIMENTAL=1 -DENABLE_PYTHON=1 -DPYTHON_VERSION=${PYTHON_VERSION} -DENABLE_LIBFIVE=1 ..
make -j3
if [ ! $? == 0 ] ; then exit 1 ; fi
rm -rf ../AppDir
make install DESTDIR=../AppDir
cd ..

export EXTRA_QT_PLUGINS=svg
/home/gsohler/git/linuxdeploy-plugin-python/linuxdeploy-plugin-python.sh --appdir AppDir
export PATH=/appimage/usr/local/bin:$PATH
rm AppDir/usr/lib/python${PYTHON_VERSION}/config-${PYTHON_VERSION}-x86_64-linux-gnu/python.o
linuxdeploy --plugin qt --output appimage --appdir AppDir
echo finished
