#!/usr/bin/env python
import os
import subprocess
from argparse import ArgumentParser

TEMPLATE = """
Usage
=====

CONCOCT uses several command line options to control the clustering, here is a complete documentation of these. These can also be viewed by typing ``concoct -h`` on the command line.::

"""

def indent(s):
    return "\n".join(map(lambda w: '\t'+w, s.splitlines()))


def help_doc_rst(script):
    """ Fetch the --help info from the script, outputs rst formatted string """
    sp = subprocess.Popen([script, "--help"],
                          stdout=subprocess.PIPE)
    stdout,stderr = sp.communicate()
    
    # Add help message to template
    return "{0}\n\n".format(indent(stdout))
    

if __name__ == "__main__":
    # Argumentparser only to add --help option
    args = ArgumentParser(description=("Generates basic command line documentation for "
                                       "concoct. This is used since the hassle of getting "
                                       "readthedocs.org to install concoct within a virtualenv "
                                       "felt like the worse option.")).parse_args()

    file_path = os.path.dirname(os.path.realpath(__file__))
    script = "concoct"
    TEMPLATE += help_doc_rst(script)

    # Print the help message to a sphinx (restructured text) markup file
    docs_path = os.path.join(file_path, 'source', 'usage.rst')
    with open(docs_path,'w') as doc_f:
        doc_f.write(TEMPLATE)
    
