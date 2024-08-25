#!/usr/bin/env bash
# Synchronize codes in developing directory to my mini server to test
if [ $# -eq 0 ]; then # copy codes to remote
    rsync -avh --exclude rst --exclude ".*" . amd:disk/xybng/test
elif [ $# -eq 2 ]; then # copy remote summary results to local
    if [ ! -d rst/$2 ]; then
        mkdir -p rst/$2
    fi
    scp amd:disk/xybng/$1/s* rst/$2
else
    echo Not defined
fi
