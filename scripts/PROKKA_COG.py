#!/usr/bin/env python
# ***************************************************************
# Name:      PROKKA_COG.py
# Purpose:   This script integrates with PROKKA and generates Cogs assignments for the protein sequences.
# Version:   0.1
# Authors:   Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk)
#                 http://userweb.eng.gla.ac.uk/umer.ijaz
# Created:   2014-01-11
# License:   Copyright (c) 2014 Computational Microbial Genomics Group, University of Glasgow, UK
#
#            This program is free software: you can redistribute it and/or modify
#            it under the terms of the GNU General Public License as published by
#            the Free Software Foundation, either version 3 of the License, or
#            (at your option) any later version.
#
#            This program is distributed in the hope that it will be useful,
#            but WITHOUT ANY WARRANTY; without even the implied warranty of
#            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#            GNU General Public License for more details.
#
#            You should have received a copy of the GNU General Public License
#            along with this program.  If not, see <http://www.gnu.org/licenses/>.
# **************************************************************/
import sys, getopt, urllib
from xml.dom import minidom
from BCBio import GFF

def get_record_from_cdd(query):
    # We need CDD accession to COG accession mapping. For this we will use NCBI eutils and parse the returned XML
    # file. For example,
    #
    # 	http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=cdd&id=223855
    #
    # returns the following XML record
    #
    #	<eSummaryResult>
    #		<DocSum>
    #			<Id>223855</Id>
    #			<Item Name="Accession" Type="String">COG0784</Item>
    #			<Item Name="Title" Type="String">CheY</Item>
    #			<Item Name="Abstract" Type="String">FOG: CheY-like receiver [Signal transduction mechanisms]</Item>
    #			<Item Name="Status" Type="Integer">0</Item>
    #			<Item Name="LivePssmID" Type="Integer">0</Item>
    #		</DocSum>
    #	</eSummaryResult>

    params = {
        'db':'cdd',
    }

    params['id'] = query
    # get citation info:
    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?' + urllib.urlencode(params)
    data = urllib.urlopen(url).read()
    xmldoc = minidom.parseString(data)
    items=xmldoc.getElementsByTagName("Item")	
    r={}
    for i in range(items.length):	
	r[items[i].getAttribute('Name')]=items[i].firstChild.data
    return r

def usage():
    print '\n'.join([
	   'Usage:',
           '\t./PROKKA_COG.py -g <gfffile> -b <blastoutfile>',
           '',
	   'Optional parameters:',
	   '\t-s (--scovs-threshold)\t\tsubject coverage threshold (Default:60)',
           '\t-p (--pident-threshold)\t\tpident threshold (Default:0)',
	   '',
           'Example usage:',
	   '',   			
           '\tStep 1: Run PROKKA_XXXXXXXX.faa with rpsblast against the  Cog database',
	   '\twith following format:',    
           '\t\t\trpsblast -query PROKKA_XXXXXXXX.faa -db Cog -evalue 0.00001', 
           '\t\t\t-outfmt \"6 qseqid sseqid evalue pident score qstart qend', 
           '\t\t\tsstart send length slen\" -out blast_output.out',
           '',
           '\tStep 2: Run this script to generate COG anotations:',
           '\t\t\t./PROKKA_COG.py -g PROKKA_XXXXXXXX.gff -b blast_output.out', 
           '\t\t\t > annotation.cog',
	   '',	
           'Refer to rpsblast tutorial: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/rpsblast/'])

def main(argv):

   # = Parameters to set ============== #
   RPSBLAST_SCOVS_THRESHOLD=60.0
   RPSBLAST_PIDENT_THRESHOLD=0.0	
   RPSBLAST_QSEQID_FIELD=0
   RPSBLAST_SSEQID_FIELD=1
   RPSBLAST_EVALUE_FIELD=2
   RPSBLAST_PIDENT_FIELD=3
   RPSBLAST_SCORE_FIELD=4
   RPSBLAST_QSTART_FIELD=5
   RPSBLAST_QEND_FIELD=6
   RPSBLAST_SSTART_FIELD=7
   RPSBLAST_SEND_FIELD=8
   RPSBLAST_LENGTH_FIELD=9
   RPSBLAST_SLEN_FIELD=10
   # = /Parameters to set ============= #


   gfffile = ''
   blastoutfile=''
			
   try:
      opts, args = getopt.getopt(argv,"hg:b:s:p:",["gfffile=","blastoutfile=","scovs-threshold=","pident-threshold="])
   except getopt.GetoptError:
      usage()	
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
	 usage()
         sys.exit()
      elif opt in ("-g", "--gfffile"):
         gfffile = arg
      elif opt in ("-b", "--blastoutfile"):
         blastoutfile = arg
      elif opt in ("-s", "--scovs-threshold"):
         RPSBLAST_SCOVS_THRESHOLD = float(arg)	
      elif opt in ("-p", "--pident-threshold"):
         RPSBLAST_PIDENT_THRESHOLD = float(arg)

   if (gfffile =='' or blastoutfile == ''):
	usage()	
	sys.exit()

   featureid_locations={}
   limits=dict(gff_type=["gene","mRNA","CDS"])
   in_handle=open(gfffile)
   for rec in GFF.parse(in_handle,limit_info=limits):
	for feature in rec.features:
		if str(feature.location.strand)!="-1":
			featureid_locations[feature.id]=[rec.id,str(feature.location.start),str(feature.location.end),'+']
		else:
			featureid_locations[feature.id]=[rec.id,str(feature.location.start),str(feature.location.end),'-']
		
   in_handle.close()

   print  '#Query\tHit\tE-value\tIdentity\tScore\tQuery-start\tQuery-end\tHit-start\tHit-end\tHit-length\tDescription\tTitle\tClass-description\tComments'	

   in_handle=open(blastoutfile)
   for line in in_handle:
        record=line.split("\t")
        if (float(record[RPSBLAST_PIDENT_FIELD])>= RPSBLAST_PIDENT_THRESHOLD and ((float(abs(int(record[RPSBLAST_SEND_FIELD])-int(record[RPSBLAST_SSTART_FIELD]))+1)/float(record[RPSBLAST_SLEN_FIELD]))*100.0)>= RPSBLAST_SCOVS_THRESHOLD):
		cogrecord=get_record_from_cdd(record[RPSBLAST_SSEQID_FIELD].split('|')[2])
		featureidlocrecord=featureid_locations[record[RPSBLAST_QSEQID_FIELD]]
		print (	featureidlocrecord[0]+'_'+record[RPSBLAST_QSEQID_FIELD][7:]+'\t'+
			cogrecord['Accession']+'\t'+
			record[RPSBLAST_EVALUE_FIELD]+'\t'+
			record[RPSBLAST_PIDENT_FIELD]+'\t'+
			record[RPSBLAST_SCORE_FIELD]+'\t'+
			record[RPSBLAST_QSTART_FIELD]+'\t'+
			record[RPSBLAST_QEND_FIELD]+'\t'+
			record[RPSBLAST_SSTART_FIELD]+'\t'+
			record[RPSBLAST_SEND_FIELD]+'\t'+
			record[RPSBLAST_LENGTH_FIELD]+'\t'+
			cogrecord['Abstract'].split('[')[0].strip()+'\t'+
			cogrecord['Title']+'\t'+cogrecord['Abstract'].split('[')[1].strip()[:-1]+'\t'+
                        '['+featureidlocrecord[1]+','+featureidlocrecord[2]+']('+featureidlocrecord[3]+')'
			)
   in_handle.close()

if __name__ == "__main__":
   main(sys.argv[1:])
