# fio_graphs
Plot fio results with matplotlib.

## Notes

* tries to remain agnostic to your actual jobs
* aggregates on job name; i.e. aggregates job duplication with numjobs and
  aggregates across multiple client nodes if the sam job file is used with fio's
  client/server mode

## Requirements
Requires python libraries

* matplotlib
* pandas
* numpy

# Running
Run ``./fio_plots.py --help``
