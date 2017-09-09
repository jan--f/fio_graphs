#!/usr/bin/env python3

import argparse
import json
import pprint
import os
import sys

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
        self.cache = {}

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

    def _aggregate_data(self):
        if not self.data['results']:
            print('ERROR...no data found.')
            sys.exit()

        if 'jobs' in self.data['results'][0]:
            result_key = 'jobs'

        if 'client_stats' in self.data['results'][0]:
            result_key = 'client_stats'

        d = {job['jobname']: {'read': 0,
                              'write': 0,
                              'r_iops': 0,
                              'w_iops': 0,
                              'lat_us': {},
                              'lat_ms': {},
                              'count': 0,
                             }
             for job in self.data['results'][0][result_key]}
        # remove 'All clients' if present
        d.pop('All clients', None)
        for result in self.data['results']:
            for job in result[result_key]:
                # Skip 'All clients' if present
                if job['jobname'] == 'All clients':
                    continue
                if job['error'] is not 0:
                    print('job {} reported an error...skipping'.format(
                        job['jobname']
                    ))
                    continue

                d[job['jobname']]['count'] += 1
                d[job['jobname']]['read'] += job['read']['bw']
                d[job['jobname']]['write'] += job['write']['bw']
                d[job['jobname']]['r_iops'] += job['read']['iops']
                d[job['jobname']]['w_iops'] += job['write']['iops']
                for k, v in job['latency_us'].items():
                    if k in d[job['jobname']]['lat_us']:
                        d[job['jobname']]['lat_us'][k] += job['latency_us'][k]
                    else:
                        d[job['jobname']]['lat_us'][k] = job['latency_us'][k]
                for k, v in job['latency_ms'].items():
                    if k in d[job['jobname']]['lat_ms']:
                        d[job['jobname']]['lat_ms'][k] += job['latency_ms'][k]
                    else:
                        d[job['jobname']]['lat_ms'][k] = job['latency_ms'][k]

        self.cache['bw'] = pandas.DataFrame(data={
            'name': [k for k in d.keys()],
            'read': [v['read'] for v in d.values()],
            'write': [v['write'] for v in d.values()]})
        self.cache['iops'] = pandas.DataFrame(data={
            'name': [k for k in d.keys()],
            'read': [v['r_iops'] for v in d.values()],
            'write': [v['w_iops'] for v in d.values()]})
        lat_data = {'lats': list(d[next(iter(d))]['lat_us'].keys())
                    + [k + '000' for k in d[next(iter(d))]['lat_ms'].keys()]}
        for name in d.keys():
            # TODO don't use values (eliminates duplicates?) but iterate over
            # the keys
            c = []
            for k in d[name]['lat_us'].keys():
                c.append(d[name]['lat_us'][k])
            for k in d[name]['lat_ms'].keys():
                c.append(d[name]['lat_ms'][k])
            lat_data[name] = c
        srt_fct = lambda x: int(str.lstrip(x, '>='))
        print(sorted(lat_data['lats'], key=srt_fct))
        self.cache['lat_dist'] = pandas.DataFrame(data=lat_data)

    def get_aggregate_bw(self):
        if not 'bw' in self.cache:
            self._aggregate_data()
        return self.cache['bw']

    def get_aggregate_iops(self):
        if not 'iops' in self.cache:
            self._aggregate_data()
        return self.cache['iops']

    def get_aggregate_lat_dist(self):
        if not 'lat_dist' in self.cache:
            self._aggregate_data()
        return self.cache['lat_dist']

    def print_(self):
        lats = self.get_aggregate_lat_dist()
        print('aggregate latency distribution')
        pprint.pprint(lats)
        print('aggregate bandwidth')
        pprint.pprint(self.get_aggregate_bw())
        print('aggregate iops')
        pprint.pprint(self.get_aggregate_iops())

    def aggregate_bw_graph(self):
        plt.clf()
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
        plt.legend((bar2[0], bar1[0]),
                   ('write', 'read')).get_frame().set_facecolor('#FFFFFF')
        plt.savefig('{}/bw_aggr.png'.format(self.args.output), bbox_inches='tight')

    def aggregate_iops_graph(self):
        plt.clf()
        dframe = self.get_aggregate_iops()
        ind = np.arange(dframe.index.size)
        if max(dframe.read) + max(dframe.write) > 9900000:
            b1_data = dframe.read / 1024
            b2_data = dframe.write / 1024
        else:
            b1_data = dframe.read
            b2_data = dframe.write
        plt.ylabel('IOPS')

        bar1 = plt.bar(ind, b1_data, self.b_width)
        bar2 = plt.bar(ind, b2_data, self.b_width, bottom=b1_data)
        plt.title('Aggregated IOPS over {} clients'.format(
            len(self.data['results'])))
        # adjust xscale if stacked is > 1000000 or so
        plt.xticks(ind, dframe.name, rotation=45)
        plt.legend((bar2[0], bar1[0]),
                   ('write', 'read')).get_frame().set_facecolor('#FFFFFF')
        plt.savefig('{}/iops_aggr.png'.format(self.args.output), bbox_inches='tight')

    def aggregate_lat_dist_graph(self):
        plt.clf()
        dframe = self.get_aggregate_lat_dist()
        ind = np.arange(dframe['lats'].size)
        # if max(dframe.read) + max(dframe.write) > 9900000:
        #     b1_data = dframe.read / 1024
        #     b2_data = dframe.write / 1024
        # else:
        #     b1_data = dframe.read
        #     b2_data = dframe.write
        # plt.ylabel('IOPS')


        plt.title('Aggregated latencyc distribution over {} clients'.format(
            len(self.data['results'])))
        def srt_fct (e):
            if not e[0].isdigit():
                return int(e.lstrip('>=')) + 1
            else:
                return int(e)
        dframe['sort'] = dframe['lats'].apply(srt_fct)
        d = dframe.sort_values(by='sort')
        # d = dframe.reindex(sorted(dframe.lats, key=srt_fct))
        pprint.pprint(d.iloc[:, :-2])
        for c in d.iloc[:, :-2]:
            plt.plot(ind, d[c])
        plt.xticks(ind, d.columns[-2], rotation=45)
        # plt.legend((bar2[0], bar1[0]),
        #            ('write', 'read')).get_frame().set_facecolor('#FFFFFF')
        plt.savefig('{}/lat_dist.png'.format(self.args.output), bbox_inches='tight')


def get_fio(path):
    return FioResults(argparse.Namespace(dir=True, path=path, output='graphs'))


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
    results.print_()
    results.aggregate_bw_graph()
    results.aggregate_iops_graph()
    results.aggregate_lat_dist_graph()


if __name__ == "__main__":
    main()
