#!/bin/bash

cp -av configs/limits/llstau_maximally-mixed/* configs/limits/llstau_mass-degenerate/

for f in $(find configs/limits/llstau_mass-degenerate/ | grep "yaml$"); do
    echo "Converting file: $f"
    sed -i "s/maximally-mixed/mass-degenerate/g" $f
done
