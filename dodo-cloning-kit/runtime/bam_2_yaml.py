#!/usr/bin/env python

import glob2
import yaml
import argparse
import re

parser = argparse.ArgumentParser(prog='bam_2_yaml.py', description='Create samples.yaml from a bam directory')
parser.add_argument('sample_file', nargs = '+')
args = parser.parse_args()

bam_files = []
all_bams = glob2.glob('bam/' + '*' + '.bam')
# for bam_file in all_bams:
#     if bam_file.startswith('bam/') and bam_file.endswith('.bam'):
#         bam_file = bam_file[4:]
#         bam_file = bam_file[:-4]
#         bam_files.append([bam_file])
        
# all_samples = []
# for sample in bam_files:
#     all_samples.append({'name': sample, 'tumor': [sample]})
# 
# yaml.dump(all_samples, args.sample_file,
#     default_style=None,
#     default_flow_style=None,
#     encoding='utf-8',
#     explicit_start=None,
#     explicit_end=None,
#     version=None,
#     tags=None,
#     canonical=None,
#     indent=None,
#     width=None,
#     allow_unicode=None,
#     line_break=None)
