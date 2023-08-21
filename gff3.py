import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('--gff3')
parser.add_argument('--bed')
args = parser.parse_args()

class GFF3:
    def __init__(self, chromosome, source, type, start, end, score, strand, phase, name):
        self.chromosome = chromosome
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.name = name


    def write_file(out, filename):
        with open(filename, "a") as outfile:  # Open the file in append mode ("a")
            outfile.write(out + '\n')


    def parse_gff3 ( input_file ): ## opens gff3 file and creates new instance of the gff3 class

        gff3_instances = []
        
        with open ( input_file , 'r' ) as gff3:
            lines = gff3.readlines()

            for line in lines:
                line = line.strip('\n')
                comment = re.compile('#')
               
                if not re.search (comment, line): ## create a new gff3 instance from the files columns and append it into a list
                    cols = line.split ('\t')
                    gff3_instance = GFF3( cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8] )
                    gff3_instances.append(gff3_instance)

            return gff3_instances  # Return the list of instances

                   
    def gff2bed ( input_file, bed_file ): ## opens gff3 file, creates new instances of the gff3 class and outputs a bed file
        gff3_instances = GFF3.parse_gff3(input_file)  # Get the list of instances from parse_gff3

        for gff3_instance in gff3_instances:
        
            if re.search("mRNA",gff3_instance.type):
                gene = re.sub(".*Name=|;.*","", gff3_instance.name)
                gene = re.sub ("-00001","",gene)
               
                out = '\t'.join([gff3_instance.chromosome, gff3_instance.start, gff3_instance.end, gff3_instance.strand, gene] )
                GFF3.write_file ( out, args.bed )
                                
GFF3.gff2bed(args.gff3, args.bed)
