#!/usr/bin/env python3

import argparse
import json
import os
import pprint

import pandas
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def get_arg_parser():
    p = argparse.ArgumentParser(
        description='Create graphs from various fio json outputs')
    p.add_argument('path', help='Source path for fio output')
    p.add_argument(
        '-d',
        '--dir', action="store_true",
        help='Read output files from a directory and consider files to be of the same run')
    p.add_argument('-o', '--output', help='output directory for graphs',
                   default='graphs')
    return p


class FioResults(object):

    def __init__(self, args):
        # two parsing modes: single file, dir with files to aggregate
        self.b_width = 0.35
        self.args = args
        self.data = {
            'results': [],
            'directory': self.args.dir
        }
        os.makedirs(self.args.output, exist_ok=True)

    def parse_data(self):
        if self.args.dir:
            self._parse_dir()
        else:
            self._parse_file(self.args.path)

    def _parse_dir(self):
        for f in os.listdir(self.args.path):
            path = '{}/{}'.format(self.args.path, f)
            if os.path.isfile(path):
                self._parse_file(path)

    def _parse_file(self, path):
        with open(path) as file_:
            try:
                d = json.load(file_)
                self.data['results'].append(d)
            except ValueError:
                print('IGNORING file {}, contains no valid JSON'.format(path))

    def get_aggregate_bw(self):
        d = {job['jobname']: {'read': 0, 'write': 0}
             for job in self.data['results'][0]['jobs']}
        for result in self.data['results']:
            for job in result['jobs']:
                d[job['jobname']]['read'] += job['read']['bw']
                d[job['jobname']]['write'] += job['write']['bw']
        return pandas.DataFrame(data={
            'name': [k for k in d.keys()],
            'read': [v['read'] for v in d.values()],
            'write': [v['write'] for v in d.values()]})

    def print(self):
        pprint.pprint(self.get_aggregate_bw())

    def aggregate_bw_graph(self):
        dframe = self.get_aggregate_bw()
        ind = np.arange(dframe.index.size)
        if max(dframe.read) + max(dframe.write) > 9900000:
            b1_data = dframe.read / 1024
            b2_data = dframe.write / 1024
            plt.ylabel('Bandwidth (MiB/s)')
        else:
            b1_data = dframe.read
            b2_data = dframe.write
            plt.ylabel('Bandwidth (KiB/s)')

        bar1 = plt.bar(ind, b1_data, self.b_width)
        bar2 = plt.bar(ind, b2_data, self.b_width, bottom=b1_data)
        plt.title('Aggregated bandwidth over {} clients'.format(
            len(self.data['results'])))
        # adjust xscale if stacked is > 1000000 or so
        plt.xticks(ind, dframe.name, rotation=45)
        plt.legend((bar1[0], bar2[0]),
                   ('read', 'write')).get_frame().set_facecolor('#FFFFFF')
        plt.savefig('{}/bw.png'.format(self.args.output), bbox_inches='tight')


def get_fio(path):
    return FioResults(argparse.Namespace(directory=True, path=path))


def main():
    a_parser = get_arg_parser()
    args = a_parser.parse_args()

    if args.dir:
        if not os.path.isdir(args.path):
            raise a_parser.ArgumentError(('-d was passed but path is not a ',
                                         'directory'))
    else:
        if os.path.isdir(args.path):
            raise a_parser.ArgumentError(('-d was not passed but path is a ',
                                         'directory'))

    results = FioResults(args)
    results.parse_data()
    results.print()
    results.aggregate_bw_graph()


if __name__ == "__main__":
    main()
