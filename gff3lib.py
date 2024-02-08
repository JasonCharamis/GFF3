
#!/usr/bin/env python3

from natsort import natsorted
import os, re
import argparse

def isfile(input_file, field=0): # Function to check if input is a file or a string

    if os.path.isfile(input_file):
        with open(input_file, "r") as file:
            # Check if single or multi-column file, if latter is true select the one with gene list (good for using in associative lists of genes)
            lines = file.readlines()
            glist = []

            for line in lines:
                num_columns = len(line.strip().split('\t'))  # Split the first line into columns based on the delimiter

                if num_columns == 1:
                    glist.append(line.strip('\n'))
                    
                elif num_columns > 1:
                    glist.append(line.split('\t')[field].strip('\n'))

            sorted_list = natsorted( glist )

            return sorted_list

    elif isinstance(input_file, str):
        glists = input_file

    return glists


class GFF3:
    
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
        return f"{self.chromosome}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{self.name}"
    

    def parse_gff3(input_file):
        gff3_instances = []
        
        with open(input_file, 'r') as gff3:
            lines = gff3.readlines()

            for line in lines:
                line = line.strip('\n')
                comment = re.compile('#')

                if not re.search(comment, line):  # Create a new GFF3 instance from the file's columns and append it to a list
                    cols = line.split('\t')
                    gff3_instance = GFF3(cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8])
                    gff3_instances.append(gff3_instance)

        return gff3_instances  # Return the list of instances

    
    def sort_gff3(input_file):
        gff3_instances = GFF3.parse_gff3(input_file)
        
        # Define sorting order for entry types
        entry_type_order = {
            "gene": 0,
            "pseudogene": 1,
            "mRNA": 2,
            "transcript": 3,
            "exon": 4,
            "CDS": 5,
            "non_canonical_five_prime_splice_site": 6,
            "non_canonical_three_prime_splice_site": 7,
            "stop_codon_read_through": 8
        }
        
        # Define a custom sorting key function
        def custom_sort_key(element):
            return (element.chromosome, element.start, entry_type_order.get(element.type, float('inf')))
        
        sorted_list = natsorted(gff3_instances, key=custom_sort_key)

        return sorted_list


    def gff2gtf ( input_file ): ## convert gff3 to bed format
        gff3_instances = GFF3.sort_gff3(input_file)

        outl = []

        for gff3_instance in gff3_instances:
            out = str()           
            out = '\t'.join( [gff3_instance.chromosome, gff3_instance.source, gff3_instance.start, gff3_instance.end, gff3_instance.score, gff3_instance.strand, gff3_instance.phase, gff3_instance.name] )
            outl.append(out)

        return ( outl )

           
    def gff2bed ( input_file ): ## convert gff3 to bed format
        gff3_instances = GFF3.sort_gff3(input_file)  # Get the list of instances from sort_gff3

        out = str()                        
        bed = []
        
        for gff3_instance in gff3_instances:

            if gff3_instance.type == "gene":
                gene = re.sub(".*Name=|;.*","", gff3_instance.name)
                out = '\t'.join([gff3_instance.chromosome, gff3_instance.start, gff3_instance.end, gene] )
                bed.append(out)

        return sort_gff3 ( bed )


            
    def extract_range (input_file, chromosome, start, end ): ## extract range from gff3 file
        gff3_instances = GFF3.sort_gff3(input_file)  # Get the list of instances from sort_gff3
        
        subset = []
        out = str ()

        for gff3_instance in gff3_instances:
            gene = re.sub(".*Name=|;.*","", gff3_instance.name)
            gene = re.sub ("-00001","",gene)

            if re.search (chromosome, gff3_instance.chromosome):

                if gff3_instance.start >= start and gff3_instance.end <= end:
                    subset.append ( gff3_instance )

        if len(subset) > 0:  
            return subset

        else:
            print ("No gff3 instances in this range.")

    
    def extract_genes ( input_file, gene_list, coords = True ): ## extract gene coordinates from gff3 file        
        subset = []
        found = False

        gff3_instances = GFF3.sort_gff3(input_file)  # Get the list of instances from sort_gff3

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

            inp = re.sub ( ".gff3", "", input_file )
            beds = []

            if coords == True:
                return sorted_list

            else:
                with open(f"{inp}_extracted_{gene_list}.bed", "w") as f:
                    for out in sorted_list:
                        if re.search("gene", out.type ):
                            out.name = re.sub (".*Name=|\.t\d|;.*","",out.name)
                            bed = '\t'.join([out.name, out.chromosome, out.start, out.end, out.strand])
                            beds.append(bed)

                return beds

        else:
            print ( "Gene list is empty or not provided." )

           
## Implementation ##
def parse_arguments():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating gff3 files.')
    parser.add_argument('--gff3', type=str, help='GFF3 file that you would like to sort.')
    parser.add_argument('--bed', action="store_true",help='BED output file.')
    parser.add_argument('-e','--extract', action="store_true", help='Option to extract range from GFF3 file.')
    parser.add_argument('-s','--sort', action="store_true", help='Option to sort the GFF3 file.')
    parser.add_argument('-gtf','--gtf', action="store_true", help='Option to convert GFF3 to GTF file format.')
    parser.add_argument('-l', '--gene_list', type=str, help='List with gene IDs that you would like to extract.')
    parser.add_argument('-c','--coords', type=str, help='Option to extract bed-like features for selected genes.')
    parser.add_argument('-r','--range', action="store_true", help='Option to extract range from GFF3 file.')
    parser.add_argument('-chr','--chromosome', type=str, help='Chromosome whose range you would like to extract.')
    parser.add_argument('-st','--start', type=str, help='Start position of the range you would like to extract.')
    parser.add_argument('-end','--end', type=str, help='End position of the range you would like to extract.')
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")

    return parser.parse_args()


def main():
    
    parser = argparse.ArgumentParser(description='Library for efficiently manipulating gff3 files.')
    args = parse_arguments()
    
    if args.gff3:

        inp = re.sub (".gff$|.gff3$","",args.gff3)
        
        if args.sort:
            with open(f"{inp}.sorted.gff3", "w") as f:
                for out in GFF3.sort_gff3(args.gff3):
                    print ( out, file = f )

        elif args.gtf:
            with open(f"{inp}.gtf", "w") as f:
                for out in GFF3.gff2gtf(args.gff3):
                    print ( out, file = f )                   
                   
        elif args.extract:
            with open(f"{inp}_extracted_{gene_list}.gff3", "w") as f:
                for out in GFF3.extract_genes(args.gff3, args.gene_list):
                    print ( out, file = f )

        elif args.range:
            with open(f"{inp}_extracted_{args.start}_{args.end}.gff3", "w") as f:
                for out in GFF3.extract_range(args.gff3, args.chromosome, args.start, args.end):
                    print ( out, file = f )
                    
        elif args.bed:
            with open(f"{inp}.bed", "w") as f:
                for out in GFF3.gff2bed(args.gff3):
                    print ( out, file = f )                  
        else:
            print("Please provide either --sort, --extract, --range, --bed or --gene_list option.")
            
    else:
        print ("Please provide a gff3 file as input.")
            

if __name__ == "__main__":
    main()
