#!/usr/bin/env bash
# Synchronize codes in developing directory to my mini server to test
if [ $# -eq 0 ]; then
    rsync -avh --exclude rst . amd:/mnt/a/store/xybng/a
    rsync -avh --exclude rst . amd:/mnt/a/store/xybng/b
elif [ $# -eq 2 ]; then
    scp amd:/mnt/a/store/xybng/a/rst/$1/s* $2
else
    echo Not defined
fi
