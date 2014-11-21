#!/bin/sh
if [ "${PWD##*/}" != "build" ] ; then
    echo "please run under build/"
    exit 1
fi

rm -rf CMakeCache.txt *.cmake CMakeFiles src test \
install_manifest.txt Makefile > /dev/null
