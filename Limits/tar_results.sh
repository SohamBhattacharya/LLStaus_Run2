#!/bin/bash -x

set -e -u

DIR=$1

find $DIR -name "*.json" | sort -V | tar -czvf ${DIR}.tar.gz -T -
