import bz2
import pickle
import gzip
import io
import os
import sys
import subprocess
import scipy.sparse as scisp
import numpy as np

def save_object(file_name, obj):
    """
    Serialize an object to a file with gzip compression. .gz will automatically be
    added if missing.

    :param file_name: output file name
    :param obj: object to serialize
    """
    with open_output(file_name, compress='gzip') as out_h:
        pickle.dump(obj, out_h)


def load_object(file_name):
    """
    Deserialize an object from a file with automatic support for compression.

    :param file_name: input file name
    :return: deserialzied object
    """
    with open_input(file_name) as in_h:
        return pickle.load(in_h)


def open_input(file_name):
    """
    Open a text file for input. The filename is used to indicate if it has been
    compressed. Recognising gzip and bz2.

    :param file_name: the name of the input file
    :return: open file handle, possibly wrapped in a decompressor
    """
    suffix = file_name.split('.')[-1].lower()
    if suffix == 'bz2':
        return bz2.BZ2File(file_name, 'r')
    elif suffix == 'gz':
        return gzip.GzipFile(file_name, 'r')
    else:
        return open(file_name, 'r')


def open_output(file_name, append=False, compress=None, gzlevel=6):
    """
    Open a text stream for reading or writing. Compression can be enabled
    with either 'bzip2' or 'gzip'. Additional option for gzip compression
    level. Compressed filenames are only appended with suffix if not included.

    :param file_name: file name of output
    :param append: append to any existing file
    :param compress: gzip, bzip2
    :param gzlevel: gzip level (default 6)
    :return:
    """

    mode = 'w' if not append else 'w+'

    if compress == 'bzip2':
        if not file_name.endswith('.bz2'):
            file_name += '.bz2'
        # bz2 missing method to be wrapped by BufferedWriter. Just directly
        # supply a buffer size
        return bz2.BZ2File(file_name, mode, buffering=65536)
    elif compress == 'gzip':
        if not file_name.endswith('.gz'):
            file_name += '.gz'
        return io.BufferedWriter(gzip.GzipFile(file_name, mode, compresslevel=gzlevel))
    else:
        return io.BufferedWriter(io.FileIO(file_name, mode))


def make_dir(path, exist_ok=False):
    """
    Convenience method for making directories with a standard logic.
    An exception is raised when the specified path exists and is not a directory.
    :param path: target path to create
    :param exist_ok: if true, an existing directory is ok. Existing files will still cause an exception
    """
    if not os.path.exists(path):
        os.mkdir(path)
    elif not exist_ok:
        raise IOError('output directory already exists!')
    elif os.path.isfile(path):
        raise IOError('output path already exists and is a file!')


def app_path(subdir, filename):
    """
    Return path to named executable in a subdirectory of the running application

    :param subdir: subdirectory of application path
    :param filename: name of file
    :return: absolute path
    """
    return os.path.join(sys.path[0], subdir, filename)



def count_fasta_sequences(file_name):
    """
    Estimate the number of fasta sequences in a file by counting headers. Decompression is automatically attempted
    for files ending in .gz. Counting and decompression is by why of subprocess calls to grep and gzip. Uncompressed
    files are also handled. This is about 8 times faster than parsing a file with BioPython and 6 times faster
    than reading all lines in Python.

    :param file_name: the fasta file to inspect
    :return: the estimated number of records
    """
    if file_name.endswith('.gz'):
        proc_uncomp = subprocess.Popen(['gzip', '-cd', file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_read = subprocess.Popen(['grep', r'^>'], stdin=proc_uncomp.stdout, stdout=subprocess.PIPE)
    else:
        proc_read = subprocess.Popen(['grep', r'^>', file_name], stdout=subprocess.PIPE)

    n = 0
    for _ in proc_read.stdout:
        n += 1
    return n


def gen_bins(fastafile,resultfile,outputdir):
    # read fasta file
    sequences={}
    if fastafile.endswith("gz"):
        with gzip.open(fastafile,'r') as f:
            for line in f:
                line=str(line,encoding="utf-8")
                if line.startswith(">"):
                    if " " in line:
                        seq,others=line.split(' ', 1)
                        sequences[seq] = ""
                    else :
                        seq=line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    else:
        with open(fastafile,'r') as f:
            for line in f:
                if line.startswith(">"):
                    if " " in line:
                        seq,others=line.split(' ', 1)
                        sequences[seq] = ""
                    else :
                        seq=line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    dic={}
    with open(resultfile,"r") as f:
        for line in f:
            contig_name,cluster_name=line.strip().split('\t')
            try:
                dic[cluster_name].append(contig_name)
            except:
                dic[cluster_name]=[]
                dic[cluster_name].append(contig_name)
    print("Writing bins in \t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    bin_name=0
    for _,cluster in dic.items():
        if bin_name < 10:
            bin = 'BIN'+ '000' + str(bin_name) + '.fa'
        elif bin_name >= 10 and bin_name < 100:
            bin = 'BIN'+ '00' + str(bin_name) + '.fa'
        elif bin_name >= 100 and bin_name < 1000:
            bin = 'BIN'+ '0' + str(bin_name) + '.fa'
        else:
            bin = 'BIN'+str(bin_name) + '.fa'
        binfile=os.path.join(outputdir,"{}".format(bin))
        with open(binfile,"w") as f:
            for contig_name in cluster:
                contig_name=">"+contig_name
                try:
                    sequence=sequences[contig_name]
                except:
                    continue
                f.write(contig_name+"\n")
                f.write(sequence+"\n")
        bin_name+=1


def gen_sub_bins(fastafile,resultfile,outputdir):
    # read fasta file
    sequences={}
    if fastafile.endswith("gz"):
        with gzip.open(fastafile,'r') as f:
            for line in f:
                line=str(line,encoding="utf-8")
                if line.startswith(">"):
                    if " " in line:
                        seq,others=line.split(' ', 1)
                        sequences[seq] = ""
                    else :
                        seq=line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    else:
        with open(fastafile,'r') as f:
            for line in f:
                if line.startswith(">"):
                    if " " in line:
                        seq,others=line.split(' ', 1)
                        sequences[seq] = ""
                    else :
                        seq=line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    dic={}
    with open(resultfile,"r") as f:
        for line in f:
            contig_name,cluster_name=line.strip().split('\t')
            try:
                dic[cluster_name].append(contig_name)
            except:
                dic[cluster_name]=[]
                dic[cluster_name].append(contig_name)
    print("Writing sub bins in \t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    bin_name=0
    for _,cluster in dic.items():
        if bin_name < 10:
            bin = 'SUB'+ '000' + str(bin_name) + '.fa'
        elif bin_name >= 10 and bin_name < 100:
            bin = 'SUB'+ '00' + str(bin_name) + '.fa'
        elif bin_name >= 100 and bin_name < 1000:
            bin = 'SUB'+ '0' + str(bin_name) + '.fa'
        else:
            bin = 'SUB'+str(bin_name) + '.fa'
        binfile=os.path.join(outputdir,"{}".format(bin))
        with open(binfile,"w") as f:
            for contig_name in cluster:
                contig_name=">"+contig_name
                try:
                    sequence=sequences[contig_name]
                except:
                    continue
                f.write(contig_name+"\n")
                f.write(sequence+"\n")
        bin_name+=1
