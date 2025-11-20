grep -h exp0 `find results/limits_signal_v12/llstau_mass-degenerate/channels_BRT2/eras_all/ | grep json$ | grep -e MStau-100 -e MStau-200 -e MStau-300 -e MStau-400 | grep -e ctau-3mm -e ctau-30mm -e ctau-300mm | sort -V`

grep -h obs `find results/limits_signal_v12/llstau_mass-degenerate/channels_BRT2/eras_all/ | grep json$ | grep -e MStau-100 -e MStau-200 -e MStau-300 -e MStau-400 | grep -e ctau-3mm -e ctau-30mm -e ctau-300mm | sort -V`
