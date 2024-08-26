#!/usr/bin/env bash

if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    brew install boost
else
    # Linux (assuming yum-based system)
    yum -y install boost-devel
fi

# Common commands for both systems
pip install meson
meson wrap install eigen
meson wrap install fmt
