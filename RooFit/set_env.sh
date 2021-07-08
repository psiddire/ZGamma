#!/bin/bash

export ANALYSIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo ANALYSIS_DIR is set to  $ANALYSIS_DIR

export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ANALYSIS_DIR/lib
