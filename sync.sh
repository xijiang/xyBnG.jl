#!/usr/bin/env bash
# Synchronize codes in developing directory to my mini server to test
rsync -avh --exclude rst . amd:/mnt/a/store/xybng/a
rsync -avh --exclude rst . amd:/mnt/a/store/xybng/b
