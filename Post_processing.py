import pandas as pd
import numpy as np
import scipy.sparse as scisp
import tqdm
import Bio.SeqIO as SeqIO
import os
import logging
import igraph as ig
import leidenalg
from utils import open_input, count_fasta_sequences

# package logger
logger = logging.getLogger(__name__)

class Postprocess:

    def __init__(self , path , checkm_result , cl):
        #######Do the recursive Leiden Algorithm##################
        cm = pd.read_csv(checkm_result , sep = ',' , header=None)
        ##ind.csv is the name of fasta file
        cm = cm.values[:,0]

        ##ref help you determine the location of contigs in the fasta file#####
        ref = {}
        numnode = 0
        rang = []
        dist_cluster={}
        self.path = path
        self.len = cl.len
        self.name = cl.name
        self.map = cl.seq_map
        self.binsize = cl.binsize

        for counter, value in enumerate(list(cl.name)):
            ref[value] = counter

        for k in range(cm.shape[0]):
            fi = str(cm[k]) + '.fa'
            logger.info('Handling fasta file {}'.format(fi))
            name = []
            with open_input(os.path.join(self.path, 'BIN' ,fi)) as multi_fasta:
                fasta_count = count_fasta_sequences(os.path.join(self.path, 'BIN' ,fi))
                for seqrec in tqdm.tqdm(SeqIO.parse(multi_fasta, 'fasta'), total=fasta_count):
                    name.append(seqrec.id) 

            index_sub = []
            for i in name:
                index_sub.append(ref[i])

            map_sub = self.map.tocsr()
            map_sub = map_sub[index_sub , :]
            map_sub = map_sub.tocsc()
            map_sub = map_sub[: , index_sub]
            map_sub = map_sub.tocoo()

            len_sub = self.len[index_sub]
            name_sub = self.name[index_sub]
            logger.debug('There are {} contigs with total length {} in {}'.format(len(len_sub),sum(len_sub),fi))

            vcount = map_sub.shape[0]
            sources = map_sub.row
            targets = map_sub.col
            wei = map_sub.data
            index = sources>targets
            sources = sources[index]
            targets = targets[index]
            wei = wei[index]
            edgelist = list(zip(sources, targets))
            g = ig.Graph(vcount, edgelist)

            part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition  , weights=wei , n_iterations=-1)
            part = list(part)

            # dict of communities
            for ci in range(len(part)):
                if np.sum(len_sub[part[ci]]) >= self.binsize:
                    rang.append(ci)
                    numnode = numnode+len(part[ci])
                    for id in part[ci]:
                        dist_cluster[name_sub[id]] = str(cm[k])+str(ci)
        logger.debug('There are {} contigs in {} sub bins'.format(numnode,len(rang)))

        ########create file for checkm################
        with open(os.path.join(self.path ,'cluster_sub.txt'),'w') as out:
            for key , value in dist_cluster.items():
                out.write(str(key)+ '\t' +str(value))
                out.write('\n')