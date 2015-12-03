# Imports all of bio_bits into bio_pieces namespace
import sys
import warnings
import bio_bits

warnings.warn(
    'bio_pieces will be removed in the future. Please use import bio_bits'
)
sys.modules['bio_pieces'] = bio_bits
