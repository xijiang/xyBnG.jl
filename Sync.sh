#!/usr/bin/env bash
# Synchronize codes in developing directory to my mini server to test
if [ $# -eq 0 ]; then # copy codes to remote
    rsync -avh --exclude rst --exclude ".*" . amd:/mnt/a/store/xybng/a
    rsync -avh --exclude rst --exclude ".*" . amd:/mnt/a/store/xybng/b
elif [ $# -eq 3 ]; then # copy remote summary results to local
    if [ ! -d rst/$3 ]; then
        mkdir -p rst/$3
    fi
    scp amd:/mnt/a/store/xybng/$1/rst/$2/s* rst/$3
else
    echo Not defined
fi
