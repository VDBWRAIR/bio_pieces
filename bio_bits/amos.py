import re
import importlib

def splitline(strvalue, delimiter, converter):
    '''
    Split a string on delimiter into 2 pieces and convert the
    right piece using converter

    :param str strvalue: String to split
    :param str delimiter: string delimiter to split with
    :param func converter: function to convert right piece of split with
    :return: tuple(left, converter(right))
    '''
    left, right = strvalue.split(delimiter)
    return (left, converter(right))

class AmosBlock(object):
    '''
    Base class for RED, CTG and TLE
    '''
    def format(self, fmt):
        return fmt.format(**self.__dict__)

    def __str__(self):
        return self.format(self.FMT)

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

class AMOS(object):
    def __init__(self, file_handle):
        self.reds = {}
        self.ctgs = {}
        self._parse_file(file_handle)

    def _parse_file(self, file_handle):
        '''
        Iterate over file_handle and parse the AMOS contents
        '''
        hdr = re.compile('^\{([A-Z]{3})$')
        # Build this up to be parsed
        curstr = ''
        curtype = None
        for line in file_handle:
            m = hdr.match(line)
            # New type found
            if m:
                _typename = m.group(1)
                # First iteration
                if curtype is None:
                    curstr = ''
                    curtype = getattr(importlib.import_module(__name__), _typename)
                elif _typename != 'TLE':
                    # Append parsed curtype to the correct list
                    # either self.reds or self.ctgs
                    p = curtype.parse(curstr)
                    classname = curtype.__name__.lower()
                    getattr(self, classname + 's')[p.iid] = p
                    # Get new type and reset curstr
                    curtype = getattr(importlib.import_module(__name__), _typename)
                    curstr = ''
            curstr += line
        # Catch the last CTG
        p = curtype.parse(curstr)
        self.ctgs[p.iid] = p

class CTG(object):
    def __init__(self, iid, eid, com, seq, qlt, tlelist):
        self.iid = iid
        self.eid = eid
        self.com = com
        self.seq = seq
        self.qlt = qlt
        self.tlelist = tlelist

    @classmethod
    def parse(klass, ctgstr):
        '''
        Parse an AMOS CTG block and return a CTG instance

        {CTG
        iid:1
        eid:contig-000001
        com:
        Ray software bla bla
        .
        seq:
        ATGC
        .
        qlt:
        DDDD
        .
        {TLE
        src:1
        off:1
        clr:0,4
        }
        }
        
        :param str ctgstr: CTG block string
        :param dict reddict: Dictionary of iid -> RED instances
        :return: CTG instance
        '''
        lines = ctgstr.splitlines()
        iid = splitline(lines[1], ':', int)[1]
        eid = splitline(lines[2], ':', str)[1]
        com = lines[4]
        seq = lines[7]
        qlt = lines[10]
        # Lines 12 -> end-1 are TLE blocks
        tlelist = []
        for i in range(12, len(lines)-1, 5):
            # The \n gets stripped from the splitlines above
            tlestr = '\n'.join(lines[i:i+5])
            tle = TLE.parse(tlestr)
            tlelist.append(tle)
        return CTG(iid, eid, com, seq, qlt, tlelist)

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

class TLE(AmosBlock):
    FMT = '{{TLE\nsrc:{src}\noff:{off}\nclr:{clr}\n}}'

    def __init__(self, src, off, clr):
        self.src = src
        self.off = off
        self.clr = clr

    @classmethod
    def parse(klass, tlestr):
        '''
        Parse an AMOS TLE block and return a TLE instance

        {TLE
        src:1
        off:1
        clr:0,10
        }

        :param str tlestr: TLE block string
        :return: TLE instance
        '''
        lines = tlestr.splitlines()
        l = len(lines)
        if l != 5:
            raise ValueError('Got {0} lines instead of 5'.format(l))
        src = splitline(lines[1], ':', int)[1]
        off = splitline(lines[2], ':', int)[1]
        clr = splitline(lines[3], ':', str)[1]
        return TLE(src, off, clr)

    def __repr__(self):
        return self.format(
            "TLE({src},{off},'{clr}')"
        )


class RED(AmosBlock):
    FMT = '{{RED\niid:{iid}\neid:{eid}\nseq:\n{seq}\n.\nqlt:\n{qlt}\n.\n}}'
    
    class MismatchedSequenceError(Exception):
        pass

    def __init__(self, iid, eid, seq, qlt):
        self.iid = iid
        self.eid = eid
        self.seq = seq
        self.qlt = qlt

    @classmethod
    def parse(klass, redstr):
        '''
        Parse an AMOS RED block and return a RED instance

        {RED
        iid:1
        eid:1
        seq:
        ATGC
        .
        qlt:
        DDDD
        .
        }

        :param str redstr: RED block string
        :return: RED instance
        '''
        lines = redstr.splitlines()
        l = len(lines)
        if l != 10:
            raise ValueError("Got {0} lines instead of 10".format(l))
        iid = splitline(lines[1], ':', int)[1]
        eid = splitline(lines[2], ':', int)[1]
        seq = lines[4]
        qlt = lines[7]
        return RED(iid, eid, seq, qlt)

    def set_from_seqrec(self, seqrec):
        '''
        Set seq, qlt based on a Bio.SeqRecord.SeqRecord
    
        :param Bio.SeqRecord.SeqRecord seqrec: SeqRecord to pull id, seq and qual from
        '''
        if str(seqrec.seq) != self.seq:
            raise RED.MismatchedSequenceError(
                "RED sequence and SeqRecord sequence do not match"
            )
        if 'phred_quality' not in seqrec._per_letter_annotations:
            raise ValueError("Given SeqRecord is missing phred_quality")
        self.eid = seqrec.id
        self.qlt = seqrec._per_letter_annotations['phred_quality']

    def __repr__(self):
        return self.format(
            "RED({iid},{eid},'{seq}','{qlt}')"
        )
