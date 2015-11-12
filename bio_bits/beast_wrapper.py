import sys
import re
import functools

import sh

def hours_states_to_sec(line, chainlength):
    r'''
    >>> line = '20000   -61227.3011     -30325.1333     -30902.1678     932.046         1.60076E-4      0.35573         3487.00         3.39 hours/million states'
    >>> chainlength = 500000000
    >>> hours_states_to_sec(line, chainlength)
    6101755
    >>> hours_states_to_sec(line.replace('million','billion'), chainlength*1000)
    6101999
    '''
    p = '^(\d+).*?((\d+\.\d+)\shours/[mb]illion\sstates)$'
    m = re.search(p, line)
    if not m:
        raise ValueError(line + ' is not a valid beast output line')

    completed_states, stateline, states = m.groups()
    states = float(states)
    completed_states = float(completed_states)

    multiplier = 1000000
    if 'billion' in stateline:
        multiplier = multiplier * 1000

    states_to_complete = (float(chainlength) - completed_states) / multiplier
    sec = states * states_to_complete * 3600
    return int(sec) 

def sec_to_time(secs):
    '''
    >>> sec_to_time(1800)
    '0d 00:30:00'
    >>> sec_to_time(3600)
    '0d 01:00:00'
    >>> sec_to_time(86400)
    '1d 00:00:00'
    >>> sec_to_time(86400+3600+1800+60+10)
    '1d 01:31:10'
    '''
    days, hours = divmod(secs, 86400)
    hours, minutes = divmod(hours, 3600)
    minutes, seconds = divmod(minutes, 60)
    return "{0}d {1:02}:{2:02}:{3:02}".format(days,hours,minutes,seconds)

def get_xmlpath_from_argv(argv):
    r'''
    Get the beast input xml file from argv list

    >>> argv = ['beast', 'foo.xml', '-beagle_GPU']
    >>> get_xmlpath_from_argv(argv)
    'foo.xml'
    >>> argv[1] = '/path/to/foo.xml'
    >>> get_xmlpath_from_argv(argv)
    '/path/to/foo.xml'
    >>> get_xmlpath_from_argv(argv[:1]+argv[-1:])
    ValueError: Missing input xml file path
    '''
    try:
        xmlpath = list(filter(lambda x: x.endswith('.xml'), argv))[0]
        return xmlpath
    except IndexError:
        raise ValueError('Missing input xml file path')

def get_chainlength_from_xml(xml_fh):
    r'''
    Get chainLength from opened xml file handle

    >>> from StringIO import StringIO
    >>> xml = StringIO('foo bar baz\n"<\n chainLength="10000" bar\n baz')
    >>> get_chainlength_from_xml(xml)
    '10000'
    >>> get_chainlength_from_xml(StringIO('foo\nbar baz'))
    ValueError
    '''
    chainlength = None
    p = re.compile('chainLength="(\d+)"')
    for line in xml_fh:
        m = p.search(line)
        if m:
            chainlength = m.group(1)
    if chainlength is None:
        raise ValueError('Could not get chainLength from {0}'.format(xml_fh))
    return chainlength

def add_est_time_to_line(chainlength, line):
    '''
    Adds estimated time column to beast output if hour/ is found in line
    '''
    if 'hours/' in line:
        sec = hours_states_to_sec(line, chainlength)
        pretty_time = sec_to_time(sec)
        sys.stdout.write(line.rstrip() + '\t' +  pretty_time + '\n')
    else:
        sys.stdout.write(line)
    sys.stdout.flush()

def run_beast(*args, **kwargs):
    '''
    Simply run beast with same args supplied to beast_wrapper but add sec_to_time
    column for remaining time left

    kwarg['_out'] is passed to sh.beast and should be a function that accepts
        a single argument that is a beast output line

    args are a list of arguments that you would normally supply to beast
    '''
    # Use _out since _iter seems do deadlock probably because of GIL
    sh.beast(*args, _out=kwargs.get('_out', sys.stdout.write))

def beast_wrapper():
    xmlpath = get_xmlpath_from_argv(sys.argv)
    chainlength = get_chainlength_from_xml(open(xmlpath))
    process_line = functools.partial(add_est_time_to_line, chainlength)
    run_beast(*sys.argv[1:], _out=process_line)

def beast_est_time():
    states = sys.argv[1]
    line = ' '.join(sys.argv[2:])
    s = hours_states_to_sec(line, states)
    sys.stdout.write(sec_to_time(s) + '\n')

if __name__ == '__main__':
    main()
