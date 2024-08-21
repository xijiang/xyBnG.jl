#!/usr/bin/env bash
# Synchronize codes in developing directory to my mini server to test
if [ $# -eq 0 ]; then # copy codes to remote
    rsync -avh --exclude rst --exclude ".*" . amd:/mnt/a/store/xybng/a
    rsync -avh --exclude rst --exclude ".*" . amd:/mnt/a/store/xybng/b
elif [ $# -eq 2 ]; then # copy remote summary results to local
    scp amd:/mnt/a/store/xybng/a/rst/$1/s* rst/$2
else
    echo Not defined
fi
