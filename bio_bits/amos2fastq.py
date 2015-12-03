'''
From what I can tell, Ray will not use all the reads.
It will, however, load all the reads as RED entries. Then onlly some get used for TLE entries.
The following will work (except for that `Bio.SeqIO` drops the description when it writes out a file, so you would have to cut out all the non-matching "+" lines.
```bash
diff actual.sorted  expected.sorted | grep "^<" -E | wc -l
```
figure out parseargs by file extension and variable number.
'''
from functools import partial
from Bio import SeqIO
import itertools
import pandas as pd
from bio_bits import amos
''' Python3 compatibility '''
from past.builtins import map , filter


def series_contains_nan(series):
    return series.isnull().any()

def add_cumcount_index(df):
    '''
    Get a dataframe with a second index. This adds the cumulative count as a MultiIndex.
    The "cumulative count" is the number of times the (possibly non-unique) index entry has occurred so far up to that row in the index.
    This adds an additional index which will make each index entry "unique", as no two instances in
    the original index can have the same cumulative count.  Used to make a Unique index from a DataFrame with non-unique index values.
    :param pd.DataFrame df: A pandas DataFrame
    :return pd.DataFrame  with a MultiIndex, including the cumulative count of the index.
    '''
    return df.set_index(df.groupby(df.index).cumcount(), append=True)


def join_non_unique_dataframes(df1, df2):
    '''
    Sourced form: http://stackoverflow.com/questions/20297021/grouping-or-merging-dataframes-pandas-python
    Get two dataframes joined on their index, even if their indices hold non-unique values.
    The joined DataFrame will have a MultiIndex; the original shared index + a cumulative count index.
    :param pandas.DataFrame df1: A DataFrame with equivalent (but non-unique) index values to df2
    :param pandas.DataFrame df2: A DataFrame with equivalent (but non-unique) index values to df1
    :return pandas.Dataframe new joined DataFrame. Will also have a "cumulative count" index
    '''
    df_multi_index_1, df_multi_index_2 = map(add_cumcount_index, (df1, df2))
    return df_multi_index_1.join(df_multi_index_2)

def df_from_collection_attributes(attributes, collection):
    '''
    Make a pandas DataFrame out of a list of objects (collection) given a list of the attributes (attributes).
    The columns will be titled the same as their named attributes.
    Note: Will fail if every object in collection does not have every attribute in attributes.
    :param list attributes: A list of strings. Valid attributes of the objects in the collection. Also used as the column names.
    :param list collection: A list of arbitrary objects. Will fail if every object in collection does not have every attribute in attributes.
    :return pandas.DataFrame  A DataFrame, matrix of the attributes for all objects in the collection.
    '''
    #TODO: It would be nice to know that the collection obj actually have these values, but that requires accessing the list
    lambdas = [lambda obj, col=col: getattr(obj, col) for col in attributes]
    return collection_as_df(lambdas, attributes, collection)

def collection_as_df(lambdas, columns, collection):
    '''
    Create a pandas DataFrame by applying a series of functions to a collection.
    :param list lambdas: a list of functions which take exactly one argument (an objects in collection)
    :param list columns: (str) the column names, in order with lambdas
    :param list collection: a list of objects which are valid arguments for the lambdas.
    :return pandas.DataFrame the lambda results on the collection objects as a Matrix.
    '''
    assert len(lambdas) == len(columns), "lambdas must have same length as columns"
    '''use list here to force the evaluation of the functions. otherwise the lambda grabs the last obj evaluated from collection, as in a closure.'''
    values = (list( func(obj) for func in lambdas) for obj in collection)
    return pd.DataFrame(values, columns=columns)

def get_df_subset(df, iterable, key=None):
    '''
    Note: Uses the index if no key is passed.
    :param pandas.DataFrame df: a DataFrame where in df[key] matches type of iterable, or if no key is supplied, df.inidex matches type of iterable
    :param iterator iterable: iterator of values which may be in df[key]/df.index
    :return pandas.Series of those objects which match the objects in iterable
    '''
    return df[df.index.isin(iterable)] if not key else df[df[key].isin(iterable)]

amos_reds_as_df = partial(df_from_collection_attributes, attributes=['seq', 'iid'])
amos_reds_as_df.__doc__ = ''' See df_from_collection_attributes; Expects list of amos.RED objects.'''
bio_records_as_df = partial(collection_as_df, [lambda rec: rec.seq.tostring(), lambda rec: rec], ['seq', 'seq_obj'])
bio_records_as_df.__doc__ = ''' See df_from_collection_attributes; Expects list of Bio.SeqRecord.'''

def flatten_multiple_seq_files(filehandles, format):
    '''
    Get a flat iterator of SeqRecords from a list of opened files.
    :param iterable filehandles: collection of valid fastq/fasta `file` objects
    :param str format: Either "fastq" or "fasta"
    :return generator of all fastq/fasta Bio.SeqRecords as a flat iterator
    '''
    open_biofile = partial(SeqIO.parse, format=format)
    return itertools.chain(*map(open_biofile, filehandles))

def get_iids(contig):
    '''
    Given a amos.CTG object, get the iids for all the Reads which aligned to form that contig.
    :param amos.CTG contig object with a tlelist attribute
    :return [int] the TLE "src" attributes (int) for each TLE in the contig, e.g. the iid for all composing reads.
    '''
    return (tle.src for tle in contig.tlelist)

def extract_dfs_by_iids(df, iids_by_ctg):
    '''
    Get a list of "sub-frames" organized according to the organization of the iids.
    :param pandas.DataFrame df: DataFrame with a 'seq' column
    :param list iids_by_ctg: 2D list of iids (ints) organized by the contig they map to
    :return a list of sub-dataframes, where each sub-frame matches the iids in the 2D list iids_by_ctg
    '''
    get_df_subset_seqs = partial(get_df_subset, df, key='iid')
    dfs_by_ctg = map(get_df_subset_seqs, iids_by_ctg)
    return dfs_by_ctg

def get_not_nulls(series):
    return series[series.notnull()]

def get_seqs_by_ctg(fastq_records, reds, iids_by_ctg):
    '''
    Transforms the fastq records and reds into pandas.DataFrame objects and joins them on the sequence string column.
    Then, this DataFrame is sliced according to ``iids_by_ctg`` and returned. The result is a 2D list of SeqRecord objects
    organized by the contigs they mapped to.
    :param list fastq_records: a collection of Bio.SeqRecord objects
    :param list reds: a collection of amos.RED objects
    :param list iids_by_ctg: a 2D list of iids (organized by contig)
    :return A 2D list of bio.SeqRecord objects organized by the contig they map to.
    '''
    fastq_df = bio_records_as_df(fastq_records)
    reds_df = amos_reds_as_df(collection=reds)
    fastq_df, reds_df = fastq_df.set_index('seq'), reds_df.set_index('seq')
    assert reds_df.shape == fastq_df.shape, "should have the same number of columns (seqs, seq_obj/iid) and rows (fastq reads / AMOS REDs."
    reds_with_seqs_df = join_non_unique_dataframes(reds_df, fastq_df)
    dfs_by_ctg = extract_dfs_by_iids(reds_with_seqs_df, iids_by_ctg)
    unfiltered_seqs_by_ctg = [df['seq_obj'] for df in dfs_by_ctg]
    if filter(series_contains_nan, unfiltered_seqs_by_ctg):
        print("Warning: AMOS records without matching fastq records found.")
    matched_seqs_by_ctg = map(get_not_nulls, unfiltered_seqs_by_ctg)
    return matched_seqs_by_ctg


def make_fastqs_by_contigs(fastqs, amos_file, fformat='fastq'):
    '''
    Loads fastq records and amos object, get the sequences and write them to the file.
    :param list fastqs: list of valid fastq `file` objects
    :param file amos_file: single amos `file` (usually .afg extension)
    :return int 0 success code
    '''
    # Make list to end IO here (keep IO in main function, and fail immediately)
    fastq_records = list(flatten_multiple_seq_files(fastqs, fformat))
    amos_obj = amos.AMOS(amos_file)
    for f in fastqs + [amos_file]:
        f.close()
    reds = amos_obj.reds.values()
    contigs = amos_obj.ctgs.values()
    iids_by_ctg = map(get_iids, contigs)
    # Do the heavy-lifting
    seqs_by_ctg = get_seqs_by_ctg(fastq_records, reds, iids_by_ctg)
    write_to_file = partial(SeqIO.write, format=fformat)
    filenames = ("{0}.{1}".format(ctg.eid, fformat) for ctg in contigs)
    map(write_to_file, seqs_by_ctg, filenames)
    return 0



