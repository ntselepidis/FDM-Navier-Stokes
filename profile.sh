#!/usr/bin/env bash

make clean all DEBUG=0
#perf record --call-graph lbr -- ./main.x
perf record --call-graph fp -- ./main.x
perf report --sort symbol --call-graph fractal,5 --percent-limit 5 >> profile.txt
