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
    return '\n'.join([
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

def read_blast_output(blastoutfile): 
    sseq_ids = []
    records = []
    with open(blastoutfile) as in_handle:
        for line in in_handle:
            line_items = line.split("\t")
            qseq = line_items[0]
            sseq = line_items[1]
            pident = line_items[3]
            send = line_items[7]
            sstart = line_items[8]
            slen = line_items[10]

            records.append({'qseqid': qseq,
                            'sseqid': sseq,
                            'pident': float(pident),
                            'send': float(send),
                            'sstart': float(sstart),
                            'slen': float(slen)})

            sseq_ids.append(sseq.split('|')[2])
    return records, sseq_ids

def read_gff_file(gfffile):
    featureid_locations={}
    limits=dict(gff_type=["gene","mRNA","CDS"])
    with open(gfffile) as in_handle:
        for rec in GFF.parse(in_handle, limit_info=limits):
            for feature in rec.features:
                featureid_locations[feature.id] = rec.id
    return featureid_locations

def main(args):
   RPSBLAST_SCOVS_THRESHOLD = args.scovs_threshold
   RPSBLAST_PIDENT_THRESHOLD = args.pident_threshold

   print  "\t".join(['#Query','Hit'])

   records, sseq_ids = read_blast_output(args.blastoutfile)

   # Retrieve the cog accession number from ncbi
   cogrecords_l = get_records_from_cdd(sseq_ids, args.email)
   cogrecords = {}
   for rec in cogrecords_l:
       cogrecords[rec['Id']] = rec

   featureid_locations = read_gff_file(args.gfffile)

   for record_d in records:
       pident_above_threshold = record_d['pident'] >= RPSBLAST_PIDENT_THRESHOLD

       # A certain fraction of the cog should be covered to avoid the same cog 
       # to be counted twice in the case when a cog is split across two or more contigs.
       alignment_length_in_subject = abs(record_d['send'] - record_d['sstart']) + 1
       percent_seq_covered = (alignment_length_in_subject / record_d['slen']) * 100.0
       seq_cov_above_threshold =  percent_seq_covered >= RPSBLAST_SCOVS_THRESHOLD
        
       if pident_above_threshold and seq_cov_above_threshold:
           cog_accession = cogrecords[record_d['sseqid'].split('|')[2]]['Accession']
           featureidlocrecord = featureid_locations[record_d['qseqid']]
           
           print(featureidlocrecord + '_' + record_d['qseqid'].split('_')[-1] + '\t' +
		cog_accession)

if __name__ == "__main__":
   parser = argparse.ArgumentParser(usage=usage())
   parser.add_argument('-g', '--gfffile', required=True,
           help='GFF file generated by e.g. prodigal')
   parser.add_argument('-b', '--blastoutfile', required=True,
           help=('Output of rpsblast run, assumed to be in tabular format whith '
               'columns: qseqid sseqid evalue pident score qstart qend sstart send length slen '))
   parser.add_argument('-s', '--scovs-threshold', type=float, default=60.0,
           help='Threshold covered in percent, default=60.0')
   parser.add_argument('-p', '--pident-threshold', type=float, default=0.0,
           help='Threshold identity in percent, default=0.0')
   parser.add_argument('-e', '--email',
           help='Email adress needed to fetch data through ncbi api')

   args = parser.parse_args()

   main(args)
