import argparse
import re
import os

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--gff3', type=str, help='GFF3 file that you would like to sort.')
parser.add_argument('-l', '--gene_list', type=str, help='List with gene IDs that you would like to extract.')
parser.add_argument('--extract', type=str, help='Option to extract range from GFF3 file.')
parser.add_argument('--chromosome', type=str, help='Chromosome whose range you would like to extract.')
parser.add_argument('--start', type=str, help='Start position of the range you would like to extract.')
parser.add_argument('--end', type=str, help='End position of the range you would like to extract.')
args = parser.parse_args()

# Define a class for GFF3 entries
class GFF3:
    __slots__ = ['chromosome', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'name']
    
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

    def __str__(self):
        return f"{self.chromosome}\t{self.source} \t {self.type} \t {self.start} \t {self.end} \t {self.score} \t {self.strand} \t {self.phase} \t {self.name}"

    # Method to parse GFF3 file and create instances
    @staticmethod
    def parse_gff3(input_file):
        gff3_instances = []

        with open(input_file, 'r') as gff3:
            lines = gff3.readlines()

            for line in lines:
                line = line.strip('\n')
                comment = re.compile('#')

                if not re.search(comment, line):
                    cols = line.split('\t')
                    gff3_instance = GFF3(cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8])
                    gff3_instances.append(gff3_instance)

        return gff3_instances

    # Method to convert GFF3 to BED format
    @staticmethod
    def gff2bed(input_file, bed_file, only_genes=False):
        gff3_instances = GFF3.parse_gff3(input_file)

        for gff3_instance in gff3_instances:
            out = str()
            gene = re.sub(".*Name=|;.*", "", gff3_instance.name)
            gene = re.sub("-00001", "", gene)
            out = '\t'.join([gff3_instance.chromosome, gff3_instance.start, gff3_instance.end, gene, gff3_instance.score, gff3_instance.strand])
            write_file(out, args.bed)

    # Method to extract range from GFF3 file
    @staticmethod
    def extract_range(input_file, chromosome, start, end):
        gff3_instances = GFF3.parse_gff3(input_file)
        subset = []

        for gff3_instance in gff3_instances:
            gene = re.sub(".*Name=|;.*", "", gff3_instance.name)
            gene = re.sub("-00001", "", gene)

            if re.search(chromosome, gff3_instance.chromosome):
                if gff3_instance.start >= start and gff3_instance.end <= end:
                    subset.append(gff3_instance)

        inp = re.sub(".gff3", "", input_file)
        with open(f"{inp}_extracted_{start}_{end}.gff3", "w") as f:
            for out in subset:
                print(out, file=f)

    # Method to extract gene coordinates from GFF3 file
    @staticmethod
    def extract_genes(input_file, gene_list):
        subset = []
        found = False
        gff3_instances = GFF3.parse_gff3(input_file)

        if os.path.isfile(gene_list):
            glist = []
            with open(gene_list, "r") as file:
                lines = file.readlines()
                for line in lines:
                    glist.append(line.strip())
        elif isinstance(gene_list, str):
            glists = gene_list

        for gff3_instance in gff3_instances:
            if re.search("gene", gff3_instance.type):
                geneid = re.sub(".*Name=|\.t\d|;.*", "", gff3_instance.name)

                if glist:
                    if geneid in glist:
                        found = True
                    else:
                        found = False
                else:
                    if geneid == gene_list:
                        found = True
                    else:
                        found = False

            if found == True:
                subset.append(gff3_instance)

        if len(subset) > 0:
            entry_type_order = {"gene": 0, "pseudogene": 1, "mRNA": 2, "transcript": 3, "exon": 4, "CDS": 5, "stop_codon_read_through": 6}
            def custom_sort_key(element):
                return (element.chromosome, element.start, entry_type_order[element.type])
            sorted_list = sorted(subset, key=custom_sort_key)
            inp = re.sub(".gff3", "", input_file)
            with open(f"{inp}_extracted_{gene_list}.gff3", "w") as f:
                for out in sorted_list:
                    print(out, file=f)
        else:
            print("Gene list is empty." )


           
## Implementation ##

#GFF3.extract_range ( args.gff3, args.chromosome, args.start, args.end  )
GFF3.extract_genes ( args.gff3, args.gene_list )

