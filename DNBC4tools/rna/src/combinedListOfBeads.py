#!/usr/bin/env python3

from typing import Optional
import warnings
import argparse
from datatable import dt, f
from tarjan import tarjan

dt.options.nthreads = 4

class BFS:
    """
    Breadth-first search algorithm
    """
    def __init__(self, graph: dict):
        """
        Initialize BFS instance

        :param graph: dictionary representing graph
        """
        self.graph = graph
        self.visited_queue = []

    def visit(self):
        """
        Visit nodes in graph using BFS algorithm
        """
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

def similarity_droplet_file(
    SIMI_DROPLET_FNAME: str,
    BEADS_LIST_FNAME: str,
    COMB_LIST_FNAME: str,
    simi_thres: float = 0.2,
    max_merge_thres: int = 1,
    COUNT_MTX_FNAME: Optional[str] = None,
    BEADS_BARCODE_FNAME: Optional[str] = None
):
    """
    Perform clustering on droplet data

    :param SIMI_DROPLET_FNAME: file name of similarity droplet data
    :param BEADS_LIST_FNAME: file name to output list of bead barcodes
    :param COMB_LIST_FNAME: file name to output list of combined clusters
    :param simi_thres: threshold for similarity between droplets (default 0.2)
    :param max_merge_thres: maximum number of droplets that can be merged (default 1)
    :param COUNT_MTX_FNAME: file name of count matrix (default None)
    :param BEADS_BARCODE_FNAME: file name of bead barcodes (default None)
    """
    # If count matrix file name is specified and bead barcode file name is specified, output list of bead barcodes
    if COUNT_MTX_FNAME is not None and BEADS_BARCODE_FNAME is not None:
        df = dt.fread(COUNT_MTX_FNAME, header=False)
        umi = dt.fread(BEADS_BARCODE_FNAME, header=False)
        beads = list(set(umi.to_list()[0]).intersection(set(df[:,f.C2].to_list()[0])))
        with open(BEADS_LIST_FNAME, 'wt') as OUT:
            for bead in beads:
                OUT.write(f"{bead}\n")

    # Read similarity droplet data and filter by similarity threshold
    D = dt.fread(SIMI_DROPLET_FNAME, header=False)
    if D.ncols == 0:
        open(COMB_LIST_FNAME,'w').close()
    else:
        D = D[f.C2 > float(simi_thres),:]

        # Read list of bead barcodes and candidates
        beads = dt.fread(BEADS_LIST_FNAME, header=False)
        candidates = set(D[:, f.C0].to_list()[0]).union(set(D[:, f.C1].to_list()[0]))
        outers = list(set(beads[:, f.C0].to_list()[0]).difference(candidates))
        if D.ncols == 5:
            D.names = ['B1', 'B2', 'similarity', 'C0', 'C1']
        if D.ncols == 3:
            D.names = ['B1', 'B2', 'similarity']
        if D.ncols == 4:
            D.names = ['B1', 'B2', 'similarity', 'pval']
        
        # Create undirected and directed graphs
        graph_ud = {}
        graph_d = {}
        for i in zip(D[:, f.B1].to_list()[0], D[:, f.B2].to_list()[0]):
            try:
                _ = graph_ud[i[0]]
            except KeyError:
                graph_ud.update({i[0]:[]})
            try:
                _ = graph_d[i[0]]
            except KeyError:
                graph_d.update({i[0]:[]})
            if i[1] not in graph_ud[i[0]]: graph_ud[i[0]].append(i[1])
            if i[1] not in graph_d[i[0]]: graph_d[i[0]].append(i[1])
            try:
                _ = graph_ud[i[1]]
            except KeyError:
                graph_ud.update({i[1]:[]})
            if i[0] not in graph_ud[i[1]]: graph_ud[i[1]].append(i[0])

        G_BFS = BFS(graph_ud)
        G_BFS.visit()
        combine_list_BFS = G_BFS.visited_queue
        combine_list_Tarjan = tarjan(graph_d)
        combine_list = []
        combine_visit = {i:False for i in candidates}
        for i in combine_list_BFS:
            if len(i) <= max_merge_thres:
                combine_list.append(i)
                for j in i:
                    combine_visit[j] = True
        for i in combine_list_Tarjan:
            if not combine_visit[i[0]]:
                combine_list.append(i)
                for j in i:
                    combine_visit[j] = True
        if sum([i for i in combine_visit.values()]) != len(combine_visit):
            warnings("unexpected inconsistency")

        with open(COMB_LIST_FNAME, 'wt') as outfile:
            cell_number = 1
            for cell in sorted(combine_list, key=lambda x: x[0]):
                for bead in cell:
                    outfile.write(f"{bead}\tCELL{cell_number}_N{len(cell)}\n")
                cell_number += 1
            for bead in sorted(outers):
                outfile.write(f"{bead}\tCELL{cell_number}_N1\n")
                cell_number += 1

def parse_args():
    parser = argparse.ArgumentParser(description="combinedListOfBeads.py")
    parser.add_argument(
        "--similarity_droplet", 
        type=str,
        help="input filename for candidate similarity list, step1 output"
        )
    parser.add_argument(
        "--beads_list", 
        type=str,
        help="input filename for non-empty beads list, step1 output. if count_mtx and beads_barcode are provided, this become output filename of recalculated beads_list"
        )
    parser.add_argument(
        "--combined_list", 
        type=str,
        help="output filename for combined beads-cells list"
        )
    parser.add_argument(
        "--max_merge", 
        type=int, 
        default=1,
        help="when BFS is attempting to merge more than MAX_MERGE beads, Tarjan will be employed"
        )
    parser.add_argument(
        "--count_mtx", 
        type=str, 
        default=None,
        help="input filename for beads-m280 counts matrix"
        )
    parser.add_argument(
        "--beads_barcode", 
        type=str, 
        default=None,
        help="input filename for beads barcode"
        )
    parser.add_argument(
        "--simi_threshold", 
        type=float, 
        default=0.2,
        help="simi_thres"
        )
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    SIMI_DROPLET_FNAME = args.similarity_droplet
    BEADS_LIST_FNAME = args.beads_list
    COMB_LIST_FNAME = args.combined_list
    max_merge_thres = args.max_merge
    COUNT_MTX_FNAME = args.count_mtx
    BEADS_BARCODE_FNAME = args.beads_barcode
    simi_thres = args.simi_threshold
    similarity_droplet_file(SIMI_DROPLET_FNAME,BEADS_LIST_FNAME,COMB_LIST_FNAME,simi_thres,max_merge_thres,COUNT_MTX_FNAME,BEADS_BARCODE_FNAME)

if __name__=='__main__':
    main()
    