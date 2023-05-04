'''
Created on May 4, 2023

@author: dan
'''

import sys
import os

from open_ephys.analysis.formats import BinaryRecording

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter


if __name__ == '__main__':
    
    # Setup argument parser

    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by djs on 2023-05-04.
  Copyright 2023 UC Davis. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc)

    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--root-dir", dest="rootdir", required=True, help='join(root,intermediate) contains struture.oebin')
    parser.add_argument("-i", "--intermediate-dir", dest="intdir", default="", required=False, help="intermediate folder, join(root,intermediate) should be the folder with structure.oebin")
    parser.add_argument("-s", "--sp2-dir", dest="sp2dir", required=True, help='folder with spike2 extract files *.bdx, *.inx')

    # Process arguments
    args = parser.parse_args()

    structure_folder = os.path.join(args.rootdir, args.intdir)
    print(f'Load open ephys events from {structure_folder}')
    
    br = BinaryRecording(structure_folder)
    br.load_events()
    print(br.events)

    
