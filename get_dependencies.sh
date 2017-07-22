#!/bin/bash

command -v wget >/dev/null 2>&1 || { echo >&2 "Error, please install wget and rerun get_dependencies.sh."; exit 1; }
command -v git >/dev/null 2>&1 || { echo >&2 "Error, please install git and rerun get_dependencies.sh."; exit 1; }
command -v md5sum >/dev/null 2>&1 || command -v md5 >/dev/null 2>&1 || { echo >&2 "Error, please install md5 or md5sum and rerun get_dependencies.sh."; exit 1; }

./scripts_include/get_eigen.sh || { echo >&2 "Error, failed to install Eigen."; exit 1; }
