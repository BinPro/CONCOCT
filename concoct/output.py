# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 2013

@author: Johannes Alneberg
"""
import os
from os.path import isdir
from os import mkdir
import sys
from datetime import datetime

class Output(object):
    path = None
    @classmethod
    def create_dir(self,dir_name):
        if isdir(dir_name):
            dir_name2="{0}_{1}".format(dir_name,datetime.now().strftime("%Y-%m-%d-%H.%M.%S"))
            mkdir(dir_name2)
            
        else:
            mkdir(dir_name)
            self.path = dir_name
