#!/usr/bin/python

# Motif-Marker
# Tucker Bower
# 2/13/21

# This program searches a .fasta file for motifs of interest and produces a vector image graphic illustrating the locations of each given motif
# and exons for each read or gene of that .fasta file. 

import argparse
import cairo
import re

parser = argparse.ArgumentParser(description = 'Input: fasta file with exons capitalized and motifs file | Output: Motif location visualization ')
parser.add_argument('-f' , '--fasta' , type = str, nargs = 1, help = 'input .fasta file and file path (ex: ./example_data/Figure_1.fasta')
parser.add_argument('-m' , '--motifs' , type = str, nargs = 1, help = 'file containing motifs of interest, one per line')
parser.add_argument('-o' , '--output' , type = str, nargs = 1, help = 'Output directory', default = '.')
parser.add_argument('-c' , '--colors' , type = str, nargs = 1, help = 'Colors file', default = ['./example_data/pastels.txt'])
args = parser.parse_args()

FASTA_FILE = open(args.fasta[0], "r")
MOTIFS_FILE = open(args.motifs[0], "r")
COLORS_FILE = open(args.colors[0], "r")

OUTPUT_DIR = args.output[0]
FASTA_PREFIX = re.findall('/(.+)\..+', args.fasta[0])[0]

OUTPUT_FILE_NAME = './' + OUTPUT_DIR + FASTA_PREFIX + '.svg'

IUPAC_CODES_DICT = {'A':'[Aa]', 'T':'[TtUu]', 'U':'[TtUu]', 'C':'[Cc]', 'G':'[Gg]', 'C':'[Cc]', 'R':'[AaGg]', 'Y':'[CcTtUU]',
    'W':'[AaTtUu]', 'S':'[GgCc]', 'M':'[AaCc]', 'K':'[GgTtUu]', 'B':'[CcGgTtUu]', 'D':'[AaGgTtUu]', 'H':'[AaCcTtUu]',
    'V':'[AaCcGg]', 'N':'[GgCcAaTtUu]'}
#Dictionary where the keys are IUPAC codes and values are sequence bases

def fasta_string_inator(fasta_file):
    '''Creates a dictionary where
    keys: gene names
    values: DNA/RNA sequence as string'''
    gene_sequence_dict = {}
    while True:
        line = fasta_file.readline().rstrip()
        if line == '':
            break
            #End of file
        if line.startswith('>'):
            #Header
            gene_ID = re.findall('>(.+)\schr', line)[0]
            #Extracting gene ID
        else:
            if gene_ID not in gene_sequence_dict:
                gene_sequence_dict[gene_ID] = line
                #Initializing string in dictionary
            else:
                gene_sequence_dict[gene_ID] += line
                #Adding all lines of sequence together into a single continuous string 
    return gene_sequence_dict


GENE_SEQUENCE_DICT = fasta_string_inator(FASTA_FILE)


def iupac_regex_inator(motifs_file):
    '''Creates a dictionary where:
    keys: motifs 
    values: regex terms that can be used to locate that motif, accounting for IUPAC ambiguous nucleotide notation'''
    motif_regex_dict = {}
    while True:
        motif = motifs_file.readline().rstrip()
        motif_regex = ''
        if motif == '':
            break
        #EOF
        for character in motif:
            motif_regex += IUPAC_CODES_DICT[character.upper()]
        #Translating every single NUC character into a regex term for every possible character that COULD be that character
        motif_regex_dict[motif] = motif_regex
        #Storing those regex terms in a dictionary
    return motif_regex_dict

MOTIF_REGEX_DICT = iupac_regex_inator(MOTIFS_FILE)


def motif_coordinates_inator(motif, sequence):
    '''Uses regular expressions to search for motifs in a sequence of bases and returns the midpoint locations of each occurence of that motif in a list.'''
    motif_locations_list = []
    for motif_match in re.finditer(motif, sequence):
        motif_start = motif_match.start() + 1
        motif_end = motif_match.end()
        motif_midpoint =  (motif_start + motif_end) /2
        motif_locations_list.append(motif_midpoint)
    return motif_locations_list

def color_brew_inator():
    '''Loads the color palette in from a user given file'''
    colors_list = []
    while True:
        line =  COLORS_FILE.readline().rstrip()
        if line == '':
            break
        color_rgb_values = re.findall('[0-9]{1,3}', line)
        colors_list.append(color_rgb_values)
    return colors_list
    #grabbing all of the colors. 
    #   The rgb values of each color are a list. 
    #   All color lists are stored in an outer list (COLORS_LIST)'

COLORS_LIST = color_brew_inator()

def make_colors_ints_inator():
    '''Changes the rgb values in the colors list from strings to ints'''
    colors_ints_list = []
    for color in COLORS_LIST:
        color_ints = []
        for rgb_string in color:
            rgb_int = int(rgb_string)
            color_ints.append(rgb_int)
        colors_ints_list.append(color_ints)
    return colors_ints_list

COLORS_INTS_LIST = make_colors_ints_inator()

def pycairo_figure_inator():
    '''Uses pycairo to create the output figure'''
    number_of_genes = len(GENE_SEQUENCE_DICT.keys())
    number_of_motifs = 1
    #Number of genes and motifs needed for adjusting height of figure


    seq_lengths_list = []
    for seq in GENE_SEQUENCE_DICT.values():
        seq_lengths_list.append(len(seq))
    max_seq_length = max(seq_lengths_list)
    #Max seq length for adjusting length of figure

    figure_surface = cairo.SVGSurface(OUTPUT_FILE_NAME, max_seq_length + 100, number_of_genes * 75 + 100)
    context = cairo.Context(figure_surface)
    #initializing cairo surface

    gene_counter = 0
    for gene in GENE_SEQUENCE_DICT.keys():
        #Looping through each gene
        gene_counter +=1
        gene_sequence = GENE_SEQUENCE_DICT[gene]
        y_coordinate = gene_counter * 30
        #Spacing each gene out by 30 units
        context.set_source_rgba(0, 0, 0, 1)
        context.select_font_face('Arial', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        context.set_font_size(12)
        context.move_to(10,y_coordinate)
        context.show_text(gene)
        #Write the name of each gene
        
        context.set_line_width(1)
        context.move_to(80,y_coordinate)
        context.line_to(len(gene_sequence) + 80, y_coordinate)
        context.stroke()
        #Draw a line of length equal to the length of the gene

        exon_start_location = re.search('[A-Z]+', gene_sequence).start() +1
        exon_end_location = re.search('[A-Z]+', gene_sequence).end()
        width = exon_end_location - exon_start_location
        #Regex searching for the location of the uppercase exons

        context.set_line_width(20)
        context.move_to(80 + exon_start_location, y_coordinate)
        context.line_to(80 + exon_end_location, y_coordinate)
        context.stroke()
        #Drawing thick black rectangle at location of exons

        context.set_line_width(2)
        motif_counter = 0
        for motif in MOTIF_REGEX_DICT.values():
            motif_counter += 1
            context.set_source_rgba(COLORS_INTS_LIST[motif_counter-1][0]/255, COLORS_INTS_LIST[motif_counter-1][1]/255, COLORS_INTS_LIST[motif_counter-1][2]/255, 0.7)
            motif_locations_list = motif_coordinates_inator(motif, gene_sequence)
            
            for location in motif_locations_list:
                context.move_to(80 + location, y_coordinate + 10)
                context.line_to(80 + location, y_coordinate - 10)
                context.stroke()
        #Marking the motifs

    context.set_source_rgba(0,0,0,1)
    context.move_to(10,y_coordinate+30)
    context.set_font_size(8)
    context.show_text('Legend')
    #Writing the word legend

    legend_counter = 0
    for motif in MOTIF_REGEX_DICT.keys():
        legend_counter += 1
        context.move_to(10, y_coordinate + 40 + 15 * legend_counter)
        context.set_line_width(12)
        context.set_source_rgba(COLORS_INTS_LIST[legend_counter-1][0]/255,COLORS_INTS_LIST[legend_counter-1][1]/255,COLORS_INTS_LIST[legend_counter-1][2]/255,1)
        context.line_to(15, y_coordinate + 40 + 15 * legend_counter)
        context.stroke()
        context.move_to(20, y_coordinate + 40 + 15 * legend_counter)
        context.set_source_rgba(0,0,0,1)
        context.show_text(motif)
        context.stroke()
        


pycairo_figure_inator()
