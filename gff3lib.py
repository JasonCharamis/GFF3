import argparse
from natsort import natsorted
import re
import os


def check_isfile ( input_file ): ## Function defined outside of class to check if input is a file or a string
    if os.path.isfile(input_file):
        glist = []
        
        with open(input_file, "r") as file:
            lines = file.readlines()
            for line in lines:
                glist.append(line.strip())
        return glist
        
    elif isinstance(input_file, str):
        glists = input_file
        return glists


class GFF3: ## Define class for GFF3 object
    
    __slots__ = [ 'chromosome', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'name' ]
    
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

    
    def parse_gff3 (input_file ): ## opens gff3 file and creates new instance of the gff3 class
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

            
    def gff2bed ( input_file, bed_file, only_genes = False ): ## convert gff3 to bed format
        gff3_instances = GFF3.parse_gff3(input_file)  # Get the list of instances from parse_gff3

        for gff3_instance in gff3_instances:

            out = str()
            gene = re.sub(".*Name=|;.*","", gff3_instance.name)
            gene = re.sub ("-00001","",gene)
                        
            out = '\t'.join([gff3_instance.chromosome, gff3_instance.start, gff3_instance.end, gene, gff3_instance.score, gff3_instance.strand] )
                
            write_file ( out, args.bed )

            
    def extract_range (input_file, chromosome, start, end ): ## extract range from gff3 file
        gff3_instances = GFF3.parse_gff3(input_file)  # Get the list of instances from parse_gff3
        
        subset = []
        out = str ()

        for gff3_instance in gff3_instances:
            gene = re.sub(".*Name=|;.*","", gff3_instance.name)
            gene = re.sub ("-00001","",gene)

            if re.search (chromosome, gff3_instance.chromosome):

                if gff3_instance.start >= start and gff3_instance.end <= end:
                    subset.append ( gff3_instance )

        inp = re.sub(".gff3","",input_file)

        with open(f"{inp}_extracted_{start}_{end}.gff3", "w") as f:
            for out in subset:
                print ( out, file = f)

     
    def extract_genes ( input_file, gene_list, coords = True ): ## extract gene coordinates from gff3 file        
        subset = []
        found = False

        gff3_instances = GFF3.parse_gff3(input_file)  # Get the list of instances from parse_gff3


        for gff3_instance in gff3_instances:            
            if re.search("gene",gff3_instance.type):
                geneid = re.sub(".*Name=|\.t\d|;.*","", gff3_instance.name)

                if check_isfile ( gene_list ):
                    if geneid in check_isfile ( gene_list ):
                        found = True
                    else:
                        found = False

                else:
                    if geneid == gene_list:
                        found = True
                    else:
                        found = False

            if found == True:
                subset.append ( gff3_instance )

                
        if len(subset) > 0:
            
            # Define sorting order for entry types
            entry_type_order = {"gene":0,"pseudogene":1,"mRNA":2,"transcript":3,"exon": 4, "CDS": 5,"stop_codon_read_through":6}
        
            # Define a custom sorting key function
            def custom_sort_key(element):
                return (element.chromosome,element.start, entry_type_order[element.type])

            sorted_list = natsorted(subset, key=custom_sort_key)  # Use the custom sorting function

            inp = re.sub ( ".gff3", "", input_file )

            if coords == True:
                with open(f"{inp}_extracted_{gene_list}.gff3", "w") as f:
                    for out in sorted_list:
                        print ( out, file = f )
            else:
                with open(f"{inp}_extracted_{gene_list}.bed", "w") as f:
                    for out in sorted_list:
                        if re.search("gene", out.type ):
                            out.name = re.sub (".*Name=|\.t\d|;.*","",out.name)
                            bed = '\t'.join([out.name, out.chromosome, out.start, out.end, out.strand])
                            print ( bed, file = f )

        else:
            print ( "Gene list is empty or not provided." )

           
## Implementation ##

def parse_args(): # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating gff3 files.')
    parser.add_argument('--gff3', type=str, help='GFF3 file that you would like to sort.')
    parser.add_argument('--bed', type=str, help='BED output file.')
    parser.add_argument('-l', '--gene_list', type=str, help='List with gene IDs that you would like to extract.')
    parser.add_argument('--extract', type=str, help='Option to extract range from GFF3 file.')
    parser.add_argument('--coords', type=str, help='Option to extract bed-like features for selected genes.')
    parser.add_argument('--chromosome', type=str, help='Chromosome whose range you would like to extract.')
    parser.add_argument('--start', type=str, help='Start position of the range you would like to extract.')
    parser.add_argument('--end', type=str, help='End position of the range you would like to extract.')
    return parser.parse_args()

def main():
    args = parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")
        return

    if args.gff3:
        if args.extract:
            GFF3.extract_range(args.gff3, args.chromosome, args.start, args.end)
        elif args.gene_list:
            GFF3.extract_genes(args.gff3, args.gene_list)
        elif args.gene_list and args.coords:
            GFF3.extract_genes(args.gff3, args.gene_list, args.coords)
        elif args.gff2bed:
            GFF3.gff2bed(args.gff3, args.bed)
        else:
            print("Please provide either --extract, --bed or --gene_list option.")
    else:
        print("Please provide a gff3 file as input.")
            

if __name__ == "__main__":
    main()

