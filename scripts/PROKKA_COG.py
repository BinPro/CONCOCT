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
import sys
from BCBio import GFF
import argparse
from Bio import Entrez

def get_records_from_cdd(queries, email):
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
    Entrez.email = email # Always tell ncbi who you are.
    search_result = Entrez.read(Entrez.epost("cdd", id=",".join(queries)))
    records = Entrez.read(Entrez.efetch(db="cdd",
            rettype='docsum',
            webenv=search_result['WebEnv'],
            query_key=search_result['QueryKey']))
    return records

def usage():
    print '\n'.join([
           'Example usage:',
	   '',   			
           '\tStep 1: Run PROKKA_XXXXXXXX.faa with rpsblast against the  Cog database',
	   '\twith following format:',    
           '\t\t\trpsblast -query PROKKA_XXXXXXXX.faa -db Cog -evalue 0.00001', 
           '\t\t\t-outfmt \"6 qseqid sseqid evalue pident score qstart qend', 
           '\t\t\tsstart send length slen\" -out blast_output.out',
           '',
           '\tStep 2: Run this script to generate COG anotations:',
           '\t\t\t./PROKKA_COG.py -g PROKKA_XXXXXXXX.gff -b blast_output.out -e mail@example.com',
           '\t\t\t > annotation.cog',
	   '',	
           'Refer to rpsblast tutorial: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/rpsblast/',
           ''])

def main(args):
   blastoutfile = args.blastoutfile
   gfffile = args.gfffile
   RPSBLAST_SCOVS_THRESHOLD = args.scovs_threshold
   RPSBLAST_PIDENT_THRESHOLD = args.pident_threshold

   # = Parameters to set ============== #
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


   featureid_locations={}
   limits=dict(gff_type=["gene","mRNA","CDS"])
   with open(gfffile) as in_handle:
       for rec in GFF.parse(in_handle, limit_info=limits):
           for feature in rec.features:
               l = [rec.id, str(feature.location.start), str(feature.location.end)]
               if feature.location.strand == 1:
                   l.append('+')
               else:
                   l.append('-')
               featureid_locations[feature.id] = l

   print  '#Query\tHit\tE-value\tIdentity\tScore\tQuery-start\tQuery-end\tHit-start\tHit-end\tHit-length\tDescription\tTitle\tClass-description\tComments'	

   sseq_ids = []
   with open(blastoutfile) as in_handle:
       for line in in_handle:
           sseq_ids.append(line.split("\t")[RPSBLAST_SSEQID_FIELD].split('|')[2])
   cogrecords_l = get_records_from_cdd(sseq_ids, args.email)
   cogrecords = {}
   for rec in cogrecords_l:
       cogrecords[rec['Id']] = rec

   in_handle=open(blastoutfile)
   for line in in_handle:
        record=line.split("\t")
        l_covered = (float(abs(int(record[RPSBLAST_SEND_FIELD])-int(record[RPSBLAST_SSTART_FIELD]))+1))

        if (float(record[RPSBLAST_PIDENT_FIELD])>= RPSBLAST_PIDENT_THRESHOLD and
                ((l_covered/float(record[RPSBLAST_SLEN_FIELD]))*100.0 >= RPSBLAST_SCOVS_THRESHOLD)):

            cogrecord = cogrecords[record[RPSBLAST_SSEQID_FIELD].split('|')[2]]
            featureidlocrecord=featureid_locations[record[RPSBLAST_QSEQID_FIELD]]
            print(featureidlocrecord[0]+'_'+record[RPSBLAST_QSEQID_FIELD][7:]+'\t'+
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
			cogrecord['Title']+'\t'+
                        cogrecord['Abstract'].split('[')[1].strip()[:-1]+'\t'+
                        '['+featureidlocrecord[1]+','+featureidlocrecord[2]+']('+featureidlocrecord[3]+')'
			)
   in_handle.close()

if __name__ == "__main__":
   parser = argparse.ArgumentParser(usage=usage())
   parser.add_argument('-g', '--gfffile', required=True,
           help='GFF file generated by e.g. prodigal')
   parser.add_argument('-b', '--blastoutfile', required=True,
           help='Output of rpsblast run')
   parser.add_argument('-s', '--scovs-threshold', type=float, default=60.0,
           help='Threshold covered in percent, default=60.0')
   parser.add_argument('-p', '--pident-threshold', type=float, default=0.0,
           help='Threshold identity in percent, default=0.0')
   parser.add_argument('-e', '--email',
           help='Email adress needed to fetch data through ncbi api')

   args = parser.parse_args()

   main(args)
