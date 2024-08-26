#!/usr/bin/env bash

if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    brew install boost
elif ldd --version 2>&1 | grep -q musl; then
    # musllinux
    apk add --no-cache boost-dev
else
    # manylinux
    yum -y install boost-devel
fi

# Common commands for both systems
pip install meson
meson wrap install eigen
meson wrap install fmt
