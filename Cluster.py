#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scipy.sparse as scisp
from math import log,exp,sqrt
import logging
import igraph as ig
import leidenalg
from sklearn import metrics
import os


# package logger
logger = logging.getLogger(__name__)


class ClusterBin:
    def __init__(self, path , contig_info , seq_map , norm_result , min_signal , min_binsize ):
        '''
        perc: threshold of spurious contacts
        min_signal: minimum signal of acceptable contigs
        min_binsize: minimum bin size of output bins
        '''
        self.path = path
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.signal = min_signal
        self.binsize = min_binsize
        self.dist_cluster={}
        self.name = []
        self.site = []
        self.len = []
        self.cov = []
        self.tax = []

        for i in range(len(contig_info)):
            temp = contig_info[i]
            self.name.append(temp.name)
            self.site.append(temp.sites)
            self.len.append(temp.length)
            self.cov.append(temp.cov)
            self.tax.append(temp.tax)

        del contig_info

        self.name = np.array(self.name)
        self.site = np.array(self.site)
        self.len = np.array(self.len)
        self.cov = np.array(self.cov)
        self.tax= np.array(self.tax)
        
        self.cov[self.cov==0] = np.min(self.cov[self.cov!=0])
        
        self.norm()  
        del self.site, self.cov 
        logger.info('Run Leiden Algorithm')
        self.leiden()
        del self.tax
        self._write_cluster()
        del self.dist_cluster

    def norm(self):
        self.seq_map = self.seq_map.tocoo()
        _map_row = self.seq_map.row
        _map_col = self.seq_map.col
        _map_data = self.seq_map.data
        _map_coor = list(zip(_map_row , _map_col , _map_data))
        coeff = self.norm_result[0:4]

        self.seq_map = self.seq_map.tolil()
        self.seq_map = self.seq_map.astype(np.float)
        for x,y,d in _map_coor:

            s1 = self.site[x]
            if s1 == 0:
                s1 = 1
            s2 = self.site[y]
            if s2 == 0:
                s2 = 1
            s = (log(s1*s2)-self.norm_result[5])/self.norm_result[6]

            l1 =  self.len[x]
            l2 =  self.len[y]
            l = (log(l1*l2)-self.norm_result[7])/self.norm_result[8]

            c1 = self.cov[x]
            c2 = self.cov[y]
            c = (log(c1*c2)-self.norm_result[9])/self.norm_result[10]
            
            d_norm = d/exp(coeff[0] + coeff[1]  * s  + coeff[2] * l + coeff[3]* c)   
    
            if d_norm > self.norm_result[4]:
                self.seq_map[x , y] = d_norm
            else:
                self.seq_map[x , y] = 0
        del _map_row, _map_col, _map_data, _map_coor

    def leiden(self):
        #########Use Leiden Algorithm to do clustering########
        map_del = self.seq_map.tocoo()
        vcount = map_del.shape[0]
        sources = map_del.row
        targets = map_del.col
        wei = map_del.data
        index = sources>targets
        sources = sources[index]
        targets = targets[index]
        wei = wei[index]
        edgelist = list(zip(sources, targets))
        g = ig.Graph(vcount, edgelist)
        
        #############determine the best resolution parameter###########
        st = []
        res_option = np.arange(0,300,10)
        res_option[0] = 1
        for res in res_option:
            part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , weights=wei , resolution_parameter = res , n_iterations = -1)
            part = list(part)
            label_true = []
            label_pred = []
            for i in range(len(part)):
                for j in self.tax[part[i]]:
                    if j != 'Unassign':
                        label_true.append(j)
                        label_pred.append(i)
            ARI_score = metrics.adjusted_rand_score(label_true, label_pred)
            NMI_score = metrics.normalized_mutual_info_score(label_true, label_pred)
            st.append((ARI_score+NMI_score)/2) 

        ind = st.index(max(st))
        res_optimal = res_option[ind]
        part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , weights=wei , resolution_parameter = res_optimal, n_iterations = -1)
        part = list(part)

        # dict of communities
        numnode = 0
        rang = []
        for ci in range(len(part)):
            if np.sum(self.len[part[ci]]) >= self.binsize:
                rang.append(ci)
                numnode = numnode+len(part[ci])
                for id in part[ci]:
                    self.dist_cluster[self.name[id]] = 'group'+str(ci)

        logger.debug('The optimal resolution is {}'.format(res_optimal))
        logger.debug('There are {} contigs in {} bins'.format(numnode , len(rang)))
        del map_del, sources, targets, wei, index, edgelist, g, part, label_true, label_pred


    def _write_cluster(self):
        ########create file for checkm################
        with open(os.path.join(self.path ,'cluster.txt'),'w') as out:
            for key , value in self.dist_cluster.items():
                out.write(str(key)+ '\t' +str(value))
                out.write('\n') 



class ClusterBin_LC:
    def __init__(self, path , contig_info , seq_map , norm_result , min_signal , min_binsize ):
        '''
        perc: threshold of spurious contacts
        min_signal: minimum signal of acceptable contigs
        min_binsize: minimum bin size of output bins
        '''
        self.path = path
        self.seq_map = seq_map
        self.norm_result = norm_result
        self.signal = min_signal
        self.binsize = min_binsize
        self.dist_cluster={}
        self.name = []
        self.len = []
        self.cov = []
        self.tax = []

        for i in range(len(contig_info)):
            temp = contig_info[i]
            self.name.append(temp.name)
            self.len.append(temp.length)
            self.cov.append(temp.cov)
            self.tax.append(temp.tax)

        del contig_info

        self.name = np.array(self.name)
        self.len = np.array(self.len)
        self.cov = np.array(self.cov)
        self.tax= np.array(self.tax)        

        
        self.norm()  
        del self.cov 
        logger.info('Run Leiden Algorithm')
        self.leiden()
        del self.tax
        self._write_cluster()
        del self.dist_cluster

    def norm(self):
        self.seq_map = self.seq_map.tocoo()
        _map_row = self.seq_map.row
        _map_col = self.seq_map.col
        _map_data = self.seq_map.data
        _map_coor = list(zip(_map_row , _map_col , _map_data))
        coeff = self.norm_result[0:4]

        self.seq_map = self.seq_map.tolil()
        self.seq_map = self.seq_map.astype(np.float)
        for x,y,d in _map_coor:

            l1 =  self.len[x]
            l2 =  self.len[y]
            l = (log(l1*l2)-self.norm_result[4])/self.norm_result[5]

            c1 = self.cov[x]
            c2 = self.cov[y]
            c = (log(c1*c2)-self.norm_result[6])/self.norm_result[7]
            
            d_norm = d/exp(coeff[0] + coeff[1] * l + coeff[2]* c)   
    
            if d_norm > self.norm_result[3]:
                self.seq_map[x , y] = d_norm
            else:
                self.seq_map[x , y] = 0
        del _map_row, _map_col, _map_data, _map_coor

    def leiden(self):
        #########Use Leiden Algorithm to do clustering########
        map_del = self.seq_map.tocoo()
        vcount = map_del.shape[0]
        sources = map_del.row
        targets = map_del.col
        wei = map_del.data
        index = sources>targets
        sources = sources[index]
        targets = targets[index]
        wei = wei[index]
        edgelist = list(zip(sources, targets))
        g = ig.Graph(vcount, edgelist)
        
        #############determine the best resolution parameter###########
        st = []
        res_option = np.arange(0,300,10)
        res_option[0] = 1
        for res in res_option:
            part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , weights=wei , resolution_parameter = res , n_iterations = -1)
            part = list(part)
            label_true = []
            label_pred = []
            for i in range(len(part)):
                for j in self.tax[part[i]]:
                    if j != 'Unassign':
                        label_true.append(j)
                        label_pred.append(i)
            ARI_score = metrics.adjusted_rand_score(label_true, label_pred)
            NMI_score = metrics.normalized_mutual_info_score(label_true, label_pred)
            st.append((ARI_score+NMI_score)/2) 

        ind = st.index(max(st))
        res_optimal = res_option[ind]
        part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , weights=wei , resolution_parameter = res_optimal, n_iterations = -1)
        part = list(part)

        # dict of communities
        numnode = 0
        rang = []
        for ci in range(len(part)):
            if np.sum(self.len[part[ci]]) >= self.binsize:
                rang.append(ci)
                numnode = numnode+len(part[ci])
                for id in part[ci]:
                    self.dist_cluster[self.name[id]] = 'group'+str(ci)

        logger.debug('The optimal resolution is {}'.format(res_optimal))
        logger.debug('There are {} contigs in {} bins'.format(numnode , len(rang)))
        del map_del, sources, targets, wei, index, edgelist, g, part, label_true, label_pred

    def _write_cluster(self):
        ########create file for checkm################
        with open(os.path.join(self.path ,'cluster.txt'),'w') as out:
            for key , value in self.dist_cluster.items():
                out.write(str(key)+ '\t' +str(value))
                out.write('\n') 









