#!/bin/bash -x

set -e -u

DIR=$1

find $DIR -type f -not -empty | grep -v configs | grep -e yaml$ -e json$ -e pdf$ -e postfit_s.root$ -e postfit_b.root$ | tar -czvf ${DIR}.tar.gz -T -

#find $DIR -type f -not -empty | grep -v configs | grep -e yaml$ -e json$ -e pdf$ -e postfit_s.root$ -e postfit_b.root$ | tar -cv -I 'zstd -10 -T15' -f ${DIR}.tar.xz -T -
