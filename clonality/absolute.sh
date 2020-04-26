#! usr/bin/env bash

IFS='_' read -ra SAMPLE_SET <<< $1
for i in $SAMPLE_SET; do echo $i; done &&
