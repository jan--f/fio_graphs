#!/usr/bin/env python3

import argparse
import json
import pprint
import os
import re
import sys

import pandas
import numpy as np
import matplotlib.colors as color
import matplotlib.pyplot as plt
# plt.style.use('ggplot')
plt.style.use('bmh')


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
        self.b_width = 0.3
        self.args = args
        self.data = {
            'results': [],
            'directory': self.args.dir
        }
        os.makedirs(self.args.output, exist_ok=True)
        self.cache = {}
        self.meta = {}

    @property
    def num_clients(self):
        if self.meta is {}:
            return 0
        # TODO fix dirty hack
        k = list(self.meta.keys())[0]
        return len(self.meta[k]['clients'])

    @property
    def num_threads(self):
        if self.meta is {}:
            return 0
        # TODO fix dirty hack
        k = list(self.meta.keys())[0]
        return self.meta[k]['count'] / len(self.meta[k]['clients'])

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

        d = {}
        for result in self.data['results']:
            if 'jobs' in result:
                result_key = 'jobs'
            elif 'client_stats' in result:
                result_key = 'client_stats'

            for job in result[result_key]:
                # Skip 'All clients' if present
                if job['jobname'] == 'All clients':
                    continue
                if job['error'] is not 0:
                    print('job {} reported an error...skipping'.format(
                        job['jobname']
                    ))
                    continue
                # Extract data from json
                if job['jobname'] not in d:
                    d[job['jobname']] = {'read': 0,
                                         'write': 0,
                                         'r_iops': 0,
                                         'w_iops': 0,
                                         'lat_us': {},
                                         'lat_ms': {},
                                         'clients': [],
                                         'options': {},
                                         'count': 0}
                    d[job['jobname']]['options'] = job['job options']

                d[job['jobname']]['count'] += 1
                if job['hostname'] not in d[job['jobname']]['clients']:
                    d[job['jobname']]['clients'].append(job['hostname'])
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

        # create data frames from extracted data
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
        self.cache['meta_clients'] = {k: v['count'] for k, v in d.items()}
        for name in d.keys():
            c = []
            for k in d[name]['lat_us'].keys():
                c.append(d[name]['lat_us'][k] / d[name]['count'])
            for k in d[name]['lat_ms'].keys():
                c.append(d[name]['lat_ms'][k] / d[name]['count'])
            lat_data[name] = c
        self.cache['lat_dist'] = pandas.DataFrame(data=lat_data)

        # collect some metadata about the jobs
        for name in d.keys():
            self.meta[name] = {
                'count': d[name]['count'],
                'clients': d[name]['clients'],
            }

    def get_aggregate_bw(self):
        if 'bw' not in self.cache:
            self._aggregate_data()
        return self.cache['bw']

    def get_aggregate_iops(self):
        if 'iops' not in self.cache:
            self._aggregate_data()
        return self.cache['iops']

    def get_aggregate_lat_dist(self):
        if 'lat_dist' not in self.cache:
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

        dframe['sort1'] = dframe['name'].apply(get_workers)
        dframe['sort2'] = dframe['name'].apply(get_op)
        dframe['sort3'] = dframe['name'].apply(get_bs)

        dframe = dframe.sort_values(by=['sort1', 'sort2', 'sort3'])

        pprint.pprint(dframe)

        if max(dframe.read) + max(dframe.write) > 1000000:
            b1_data = dframe.read / 1024
            b2_data = dframe.write / 1024
            plt.ylabel('Bandwidth (MiB/s)')
        else:
            b1_data = dframe.read
            b2_data = dframe.write
            plt.ylabel('Bandwidth (KiB/s)')

        randr_color = color.to_rgba('C4', 0.8)
        randw_color = color.to_rgba('C4', 1)
        seqr_color = color.to_rgba('C3', 0.8)
        seqw_color = color.to_rgba('C3', 1)
        colors1 = ([randr_color] * 5 + [seqr_color] * 5) * 4
        colors2 = ([randw_color] * 5 + [seqw_color] * 5) * 4
        bar1 = plt.bar(ind, b1_data, self.b_width, color=colors1)
        bar2 = plt.bar(ind, b2_data, self.b_width, bottom=b1_data,
                       color=colors2)
        plt.title('Aggregated bandwidth over {} clients'.format(
            self.num_clients))
        # adjust xscale if stacked is > 1000000 or so
        plt.xticks(ind, dframe.name, rotation=-45, ha='left',
                   rotation_mode='anchor')

        plt.legend((bar2[0], bar1[0], bar2[5], bar1[5]),
                   ('random write', 'random read', 'seq write','seq read')).get_frame().set_facecolor('#FFFFFF')
        fig = plt.gcf()
        fig.set_size_inches(24, 15)
        plt.savefig('{}/bw_aggr.png'.format(self.args.output), bbox_inches='tight')
        plt.savefig('{}/bw_aggr.svg'.format(self.args.output), bbox_inches='tight')

    def aggregate_iops_graph(self):
        plt.clf()
        dframe = self.get_aggregate_iops()
        ind = np.arange(dframe.index.size)
        plt.ylabel('IOPS')

        dframe['sort1'] = dframe['name'].apply(get_workers)
        dframe['sort2'] = dframe['name'].apply(get_op)
        dframe['sort3'] = dframe['name'].apply(get_bs)

        dframe = dframe.sort_values(by=['sort1', 'sort2', 'sort3'])


        if max(dframe.read) + max(dframe.write) > 9900000:
            b1_data = dframe.read / 1024
            b2_data = dframe.write / 1024
        else:
            b1_data = dframe.read
            b2_data = dframe.write

        randr_color = color.to_rgba('C4', 0.8)
        randw_color = color.to_rgba('C4', 1)
        seqr_color = color.to_rgba('C3', 0.8)
        seqw_color = color.to_rgba('C3', 1)
        colors1 = ([randr_color] * 5 + [seqr_color] * 5) * 4
        colors2 = ([randw_color] * 5 + [seqw_color] * 5) * 4
        bar1 = plt.bar(ind, b1_data, self.b_width, color=colors1)
        bar2 = plt.bar(ind, b2_data, self.b_width, bottom=b1_data,
                       color=colors2)
        plt.title('Aggregated IOPS over {} clients'.format(
            self.num_clients))
        plt.yscale('log')
        # adjust xscale if stacked is > 1000000 or so
        plt.xticks(ind, dframe.name, rotation=-45, ha='left',
                   rotation_mode='anchor')
        plt.legend(
            (bar2[0], bar1[0], bar2[5], bar1[5]),
            ('random write', 'random read', 'seq write', 'seq read')
        ).get_frame().set_facecolor('#FFFFFF')
        fig = plt.gcf()
        fig.set_size_inches(16, 9)
        plt.savefig('{}/iops_aggr.png'.format(self.args.output), bbox_inches='tight')
        plt.savefig('{}/iops_aggr.svg'.format(self.args.output), bbox_inches='tight')

    def aggregate_lat_dist_graph(self):
        plt.clf()
        dframe = self.get_aggregate_lat_dist()
        ind = np.arange(dframe['lats'].size)
        plt.ylabel('% of all IOPS')
        plt.xlabel('Completion latency bucket [μs]')


        plt.title('Aggregated latency distribution over {} clients'.format(
            self.num_clients))
        def strip_fct (e):
            if not e[0].isdigit():
                return int(e.lstrip('>=')) + 1
            else:
                return int(e)

        # TODO add colour gradient based on number of workes maybe using
        # colormaps? just darker -> more workers...maybe use num workers as
        # darkening factor or so
        def get_color(c):
            bs = get_bs(c)
            if bs == 4:
                return 'C1'
            elif bs == 64:
                return 'C2'
            elif bs == 1024:
                return 'C3'
            elif bs == 4096:
                return 'C4'
            else:
                return 'C5'

        def get_alpha(c):
            w = get_workers(c)
            if w == 1:
                return 0.4
            if w == 4:
                return 0.6
            if w == 8:
                return 0.8
            if w == 16:
                return 1

        dframe['sort'] = dframe['lats'].apply(strip_fct)
        d = dframe.sort_values(by='sort')
        pprint.pprint(d)
        legend = []
        for c in d.iloc[:, :-2]:
            a = get_alpha(c)
            col = color.to_rgba(get_color(c), a)
            line = plt.plot(ind, d[c].cumsum(), color=col)
            legend.append((line[0], c))
        plt.xticks(ind, d['lats'], rotation=45)

        legend = sorted(legend, key=lambda tup: get_workers(tup[1]))
        legend = sorted(legend, key=lambda tup: get_bs(tup[1]))

        plt.legend([l[0] for l in legend],
                   [l[1] for l in legend]).get_frame().set_facecolor('#FFFFFF')
        fig = plt.gcf()
        fig.set_size_inches(16, 9)
        plt.savefig('{}/lat_dist.png'.format(self.args.output), bbox_inches='tight')
        plt.savefig('{}/lat_dist.svg'.format(self.args.output), bbox_inches='tight')


def get_workers(val):
    return int(re.search('^\d+', val)[0])


def get_bs(val):
    bs = re.findall('\d+', val)[1]
    u = re.findall('[k,m]', val)[0]
    if u == 'm':
        return int(bs) * 1024
    return int(bs)


def get_op(val):
    return val.split('_')[-1]


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
