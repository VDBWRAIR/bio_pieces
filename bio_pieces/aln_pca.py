from __future__ import print_function

import itertools
import os
import sys

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch

'''
Builds pca 3d plots based on an identity matrix for a given fasta alignment file

Heavily based on this:
http://sebastianraschka.com/Articles/2014_pca_step_by_step.html#introduction

Also using code from this:
https://zulko.wordpress.com/2012/09/29/animate-your-3d-plots-with-pythons-matplotlib/

TODO:
    - Somehow detect different data sets(aka, each sample's data should be graphed
      as different color)
    - Make interactive? Seems like if we can launch gtk matplotlib then it should
      just work?
'''

def pairwise_identity(seq1, seq2, substitution_matrix=None):
    '''
    Compare two sequences by counting only bases that
    are identical between them

    Substitution matrix rules:
        - len of rows must equal len of cols
        - must contain a key '~' for rows and columns that will be used in case
          the substitution is not in the rest of the matrix. this is the 'all' key
        - keys can be upper and/or lower case as well as the sequences as the
          case is ignored

    :param str seq1: string of chars
    :param str seq2: string of chars
    :param mapping substitution_matrix: 2D mapping that allows lookup of
        substitutions at each position that returns the identity score for that
        substitution(by default, uses matrix of 1's for any match and 0 for 
        non-matches). If None is given, the simplest lookup table is used of 0's 
        and 1.
    :return: sum of same base positions
    '''
    missingallkey = ValueError('Substitution matrix is missing the ~ key that is ' \
        'used to fetch values not contained in the substitution matrix')
    # Lazy alias
    sm = substitution_matrix
    if len(seq1) != len(seq2):
        raise ValueError('Sequence lengths did not match')

    ident = 0
    for x,y in itertools.izip(seq1,seq2):
        if substitution_matrix is None:
            if x.lower() == y.lower():
                ident += 1
        else:
            # Ignore case complexity
            x_id = sm.get(
                x.upper(),
                sm.get(
                    x.lower(),
                    sm.get('~', None)
                )
            )
            if x_id is None:
                raise missingallkey
            sub_score = x_id.get(
                y.upper(),
                x_id.get(
                    y.lower(),
                    x_id.get('~', None)
                )
            )
            if sub_score is None:
                raise missingallkey
            ident += sub_score
    return ident

def index_fasta(aln_fh):
    '''
    Return a pandas.Series for the fasta sequence alignment in order
    of the sequences in the file

    :param file aln_fh: fasta file like iterator
    :return: pandas.Series indexed by Bio.SeqRecord.id
    '''
    # Build list of tuples (id,seq) and convert to series later
    seq_index = OrderedDict()
    for record in SeqIO.parse(aln_fh, 'fasta'):
        seq_index[record.id] =  str(record.seq)
    return pd.Series(seq_index)

def identity_matrix(aln, subst_matrix=None):
    '''
    Build an identity matrix from all the pairwise identities
    of all sequences in the supplied fasta index

    n^2 loop over sequences to generate all identities.
    Only does top right calculations of matrix and copies values
    into bottom left as they are identical

    If subst_matrix is provided then it iwll be used as the lookup
    to build the pairwise identity such that each nucleotide pair will be looked
    up in the matrix to get the value that will represent that identity and all
    of the identities will be sumed at the end.

    By default the matrix that is used is a simple pairwise as follows:
    * Any matching pairs yield 1
    * Any mismatch pairs yeild 0

    :param mapping aln: Aligned indexed fasta
    :param matrix subst_matrix: Single Nucleotide Substitution matrix as a 2D lookup
        table(dictionary like) or None to use the default.
    :return: pandas.DataFrame representing identity matrix
    '''
    id_matrix = np.empty([len(aln),len(aln)])
    for i in range(len(aln)):
        for j in range(len(aln)):
            #print(i,j)
            # Only compute top right of matrix
            if i > j:
                #print("Copying {0}{1} from {1}{0}".format(j,i))
                id_matrix[i][j] = id_matrix[j][i]
            else:
                _id = pairwise_identity(aln[i], aln[j], subst_matrix)
                #print("Pident of {0} and {1} is {2}".format(aln[i],aln[j],_id))
                id_matrix[i][j] = _id
    return pd.DataFrame(id_matrix, index=aln.keys(), columns=aln.keys())

class Arrow3D(FancyArrowPatch):
    '''
    This is an arrow class for matplotlib to draw 3d arrow..i think
    '''
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def build_pca_from_fasta(fastapath, outputfile, substitution_matrix):
    '''
    Prototype function to build pca graphics for aligned fasta file

    :param str fastapath: File path of fasta alignment
    :param str outputfile: Path to write output file image
    :param dict substitution_matrix: 2D dictionary like object that supports matrix
        lookup for Nucleotides or None to use default lookup table
    '''
    # First get identity matrix for fasta
    ifasta = index_fasta(open(fastapath))
    # Now parse substitution matrix if needed
    if substitution_matrix is not None:
        subs_matrix = pd.read_csv(
            substitution_matrix, header=0, delimiter='\s+', index_col=0)
        substitution_matrix = subs_matrix.to_dict()
    id_matrix = identity_matrix(ifasta, substitution_matrix)

    # Build vector of means for all sequences in id_matrix
    mean_vector = id_matrix.mean().as_matrix()
    # Get ndarray matrix instead of dataframe
    m = id_matrix.as_matrix()

    # builds scatter matrix(unsure what that is)
    scatter_matrix = np.zeros(m.shape)
    for i in range(m.shape[1]):
        scatter_matrix += (m[:,i].reshape(m.shape[0],1) - mean_vector).dot(
            (m[:,i].reshape(m.shape[0],1) - mean_vector).T)

    # Get eigen values and vectors(unsure what they are)
    eig_val_sc, eig_vec_sc = np.linalg.eig(scatter_matrix)

    # Get tuple of eigen val,vec so we can sort together
    eig_pairs = []
    for i in range(len(eig_val_sc)):
        a = np.abs(eig_val_sc[i])
        b = eig_vec_sc[:,i]
        eig_pairs.append((a,b))
        
    # Sort the (eigenvalue, eigenvector) tuples from high to low
    eig_pairs.sort(key=lambda x: x[0])
    eig_pairs.reverse()

    # Don't know what np.hstack does, but I think we are using only the
    # top 3 dimensions to build a matrix for transformation later
    matrix_w = np.hstack((
        eig_pairs[0][1].reshape(m.shape[0],1),
        eig_pairs[1][1].reshape(m.shape[0],1),
        eig_pairs[2][1].reshape(m.shape[0],1)
    ))

    # Now get transformed matrix from original against our matrix_w
    transformed = matrix_w.T.dot(m)

    # Get means for x,y and z of transformed
    mean_x = transformed[0,:].mean()
    mean_y = transformed[1,:].mean()
    mean_z = transformed[2,:].mean()

    # Start the 3D plot
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')
    plt.rcParams['legend.fontsize'] = 10
    # Plot the transformed matrix
    ax.plot(transformed[0,:], transformed[1,:], transformed[2,:],
            'o', markersize=8, color='blue', alpha=0.5, label='Seqs')
    # Plot the mean location(single dot)
    ax.plot([mean_x], [mean_y], [mean_z],
           'o', markersize=10, color='red', alpha=0.5, label='Mean')

    plt.title('Samples for class 1 and class 2')
    ax.legend(loc='upper right')

    # Now we will write a png file that contains the 6 rotated images
    # TODO: Make it interactive?
    angles = np.linspace(0,360,21)[:-1] # A list of 20 angles between 0 and 360
    angles = np.array([0,45,90])
    elevs = np.array([0,45,90])
     
    # create an animated gif (20ms between frames)
    rotanimate(ax, angles, outputfile, elevations=elevs, delay=100, width=10, heigh=10)

# I snipped this from the following site and made a few modififications to allow
# you to also rotate over x and y axis(elevation and angle
#https://zulko.wordpress.com/2012/09/29/animate-your-3d-plots-with-pythons-matplotlib/
##### TO CREATE A SERIES OF PICTURES
 
def make_views(ax, angles, elevations=None, width=4, height = 3,
                prefix='tmprot_',**kwargs):
    """
    Makes jpeg pictures of the given 3d ax, with different angles.
    Args:
        ax (3D axis): te ax
        angles (list): the list of angles (in degree) under which to
                       take the picture.
        width,height (float): size, in inches, of the output images.
        prefix (str): prefix for the files created. 
     
    Returns: the list of files created (for later removal)
    """
     
    files = []
    ax.figure.set_size_inches(width,height)
    ax1 = ax.figure.gca()
    ax1.axes.set_xticklabels([])
    ax1.axes.set_yticklabels([])
    ax1.axes.set_zticklabels([])
    ax1.axes.set_xlabel('X Axis')
    ax1.axes.set_ylabel('Y Axis')
    ax1.axes.set_zlabel('Z Axis')
     
    i = 0
    for elevation in elevations:
        for angle in angles:
            ax.view_init(elev=elevation, azim=angle)
            fname = '%s%03d.png'%(prefix,i)
            ax.figure.savefig(fname)
            files.append(fname)
            i += 1
     
    return files
 
##### TO TRANSFORM THE SERIES OF PICTURE INTO AN ANIMATION
def make_movie(files,output, fps=10,bitrate=1800,**kwargs):
    """
    Uses mencoder, produces a .mp4/.ogv/... movie from a list of
    picture files.
    """
     
    output_name, output_ext = os.path.splitext(output)
    command = { '.mp4' : 'mencoder "mf://%s" -mf fps=%d -o %s.mp4 -ovc lavc\
                         -lavcopts vcodec=msmpeg4v2:vbitrate=%d'
                         %(",".join(files),fps,output_name,bitrate)}
                          
    command['.ogv'] = command['.mp4'] + '; ffmpeg -i %s.mp4 -r %d %s'%(output_name,fps,output)
     
    print(command[output_ext])
    output_ext = os.path.splitext(output)[1]
    os.system(command[output_ext])
 
def make_gif(files,output,delay=100, repeat=True,**kwargs):
    """
    Uses imageMagick to produce an animated .gif from a list of
    picture files.
    """
     
    loop = -1 if repeat else 0
    os.system('convert -delay %d -loop %d %s %s'
              %(delay,loop," ".join(files),output))
 
def make_strip(files,output,**kwargs):
    """
    Uses imageMagick to produce a .jpeg strip from a list of
    picture files.
    """
     
    os.system('montage -tile 1x -geometry +0+0 %s %s'%(" ".join(files),output))
     
##### MAIN FUNCTION
def rotanimate(ax, angles, output, **kwargs):
    """
    Produces an animation (.mp4,.ogv,.gif,.jpeg,.png) from a 3D plot on
    a 3D ax
     
    Args:
        ax (3D axis): the ax containing the plot of interest
        angles (list): the list of angles (in degree) under which to
                       show the plot.
        output : name of the output file. The extension determines the
                 kind of animation used.
        **kwargs:
            - width : in inches
            - heigth: in inches
            - framerate : frames per second
            - delay : delay between frames in milliseconds
            - repeat : True or False (.gif only)
    """
         
    output_ext = os.path.splitext(output)[1]
 
    files = make_views(ax,angles, **kwargs)
     
    D = { '.mp4' : make_movie,
          '.ogv' : make_movie,
          '.gif': make_gif ,
          '.jpeg': make_strip,
          '.png':make_strip}
           
    D[output_ext](files,output,**kwargs)
     
    for f in files:
        os.remove(f)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'fasta',
        help='Aligned fasta file to build pca graphics from'
    )
    parser.add_argument(
        '--substitution-matrix',
        '-s',
        default=None,
        help='Substitution matrix file to use'
    )
    parser.add_argument(
        '--outfile',
        default='pca.png',
        help='Output path[Default: %(default)s]'
    )
    return parser.parse_args()

def main():
    args = parse_args()
    build_pca_from_fasta(args.fasta, args.outfile, args.substitution_matrix)

if __name__ == '__main__':
    main()
