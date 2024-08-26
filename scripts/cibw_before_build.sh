#!/usr/bin/env bash

yum -y install boost-devel
pip install meson
meson wrap install eigen
meson wrap install fmt
