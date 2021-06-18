#!/usr/bin/env python3


########
#
#  get.combinedListOfBeads.py
#
#  Author: YU Hao (yuhao@genomics.cn)
#          BAI Yinqi (baiyinqi@genomics.cn)
#          LIU weiqing (liuweiqing@genomics.cn)
#
#  Date:   2021-01-20 (v1.0)
#  Date:   2021-03-04 (v1.1) To fix up a bug of BFS
#          2021-03-16 (v1.2) To add a function for completing out BFS network
#          2021-03-18 (v1.3) To add a limit for CPU using
#
########




from datatable import dt, f
from sys import argv




########
dt.options.nthreads = 4




########
class BFS():
    """"""
    def __init__(self, graph):
        """"""
        self.graph = graph
        self.visited_queue = []


    def visit(self):
        """"""
        while list(self.graph.keys()):
            _current_queue = [list(self.graph.keys()).pop(0)]
            _current_visited_queue = list()

            while _current_queue:
                _current_node = _current_queue.pop(0)

                if _current_node not in _current_visited_queue:
                    _current_visited_queue.append(_current_node)

                    if _current_node in self.graph:
                        for _node in self.graph[_current_node]:
                            if _node in self.graph:
                                _current_queue.append(_node)

                        del self.graph[_current_node]

                    else:
                        _current_queue.append(_current_node)

            self.visited_queue.append(_current_visited_queue)
########



#weiwqing D = dt.fread('DATA/' + argv[1] + '/Similarity.droplet.filtered.csv', header=False)

D = dt.fread(argv[1], header=False)

if D.ncols == 5:
    D.names = ['B1', 'B2', 'similarity', 'C0', 'C1']

if D.ncols == 3:
    D.names = ['B1', 'B2', 'similarity']




D0 = D[:, [f.B1, f.B2]]
D1 = D[:, [f.B2, f.B1]]

D1.names = ['B1', 'B2']
D1.key = ['B1', 'B2']

D = dt.rbind(D0, D1)




graph = {}

for i in zip(D[:, f.B1].to_list()[0], D[:, f.B2].to_list()[0]):
    try:
        graph[i[0]]

    except KeyError:
        graph.update({i[0]:[]})

    graph[i[0]].append(i[1])

G = BFS(graph)
G.visit()


#weiqing with open('DATA/' + argv[1] + '/' + argv[1] + '.beads_combined_list.txt', 'wt') as OU:
with open(argv[2], 'wt') as OU:
    n = 1
    for i in G.visited_queue:
        for j in i:
            print(j, '\t', 'CELL', n, '_N', len(i), sep='', file=OU)
        n += 1
