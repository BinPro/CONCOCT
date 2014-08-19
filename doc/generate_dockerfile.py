#!/usr/bin/env python
import os
from argparse import ArgumentParser
import jinja2
import sys

def main(args):
    with open(args.template) as tf:
        t = jinja2.Template(tf.read())
    
    with open(args.Dockerfile, 'w') as df:
        df.write(t.render(version=args.version))

if __name__ == "__main__":
    # Argumentparser only to add --help option
    parser = ArgumentParser(description=("Generates the Dockerfile for the given release by changing the "
                                       " version number in the Dockerfile template."))
    parser.add_argument("template", help="Path to the Dockerfile template")
    parser.add_argument("Dockerfile", help="Path to where new dockerfile will be printed")
    parser.add_argument("version", help=("Version number for current release, "
            " need to be present as a tag on github."))

    args = parser.parse_args()
   
    main(args)
