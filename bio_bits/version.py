import sys
from bio_bits import __release__

def get_release():
    return __release__

def print_release():
    sys.stdout.write("{0}\n".format(get_release()))

def main():
    print_release()
