#!/bin/bash -x

set -e -u

DIR=$1

find $DIR | grep -v configs | grep -e yaml$ -e json$ -e pdf$ | sort -V | tar -czvf ${DIR}.tar.gz -T -
