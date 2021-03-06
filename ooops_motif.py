#!/usr/bin/env python

# imports
import argparse
import os
import copy
import random
import cairo
import argparse


# todo argparse input
# import arguments from command line
def get_args():
        parser = argparse.ArgumentParser(description='mark up motifs on a gene')
        parser.add_argument("-f", "--fasta_file_in", type=str, help='specifies input fasta file.')
        parser.add_argument("-m", "--motif_file_in", type=str, help='specifies input motif file.')



        return parser.parse_args()

# get arg parse arguments
parseArgs = get_args()

motif_file = parseArgs.motif_file_in
fasta_file = parseArgs.fasta_file_in
margin = 10

# svg filename 
svg_name = fasta_file.split(".fa")[0] + ".svg"



'''
GENNE CLASS
'''
class GENE:

    def __init__(self, gene_seq, gene_num):

        self.gene_length = len(gene_seq)
        self.gene_num    = gene_num

    def draw(self, context):
        context.set_source_rgb(0, 0, 0)
        context.set_line_width(3)
        context.move_to(margin, self.gene_num)
        context.line_to(self.gene_length + margin, self.gene_num)
        context.stroke()


'''
EXON CLASS
'''
class EXON:

    def __init__(self, gene_seq, gene_num):

        self.gene_seq = gene_seq
        self.Exons    = self.find_exons()
        self.gene_num    = gene_num

    def find_exons(self):

        # empty list for all exon coordinates
        all_exon_coordinates = []
        
        # check if first base is part of an exon or an intron
        if self.gene_seq[0] == self.gene_seq[0].upper:
            exon = True
        else:
            exon = False

        # initialize empty list to store first exon
        coordinate = []

        # iterate through the string 
        for i,base in enumerate(self.gene_seq):

            # find a nex exon
            if exon == False and base == base.upper() and len(coordinate) == 0:
                exon = True
                coordinate.append(i)

            # fall off an exon
            elif exon == True and base != base.upper() and len(coordinate) == 1:
                exon = False
                coordinate.append(i-1)
                all_exon_coordinates.append(coordinate)
                coordinate = []
        
        return all_exon_coordinates

    def draw(self, context):
        exon = self.Exons[0]
        context.set_source_rgb(0, 0, 0)
        context.set_line_width(10)
        context.move_to(exon[0] + margin, self.gene_num)
        context.line_to(exon[1] + margin, self.gene_num)
        context.stroke()

'''
MOTIF CLASS
'''

class MOTIF:

    def __init__(self, gene_seq, gene_num, motifs, motif_colors):

 
        self.iupac_dic          =            {
                        "A":["A","W","M","R","a","w","m","r"],
                        "C":["C","S","M","Y","c","s","m","y"],
                        "G":["G","S","K","R","g","s","k","r"],
                        "T":["U","T","W","K","Y","u","t","w","k","y"],
                        "U":["T","U","W","K","Y","t","u","w","k","y"],
                        "W":["A","T","W","a","t","w"        ],
                        "S":["S","C","G","s","c","g"    ],
                        "M":["M","A","C","m","a","c"    ],
                        "K":["K","G","T","k","g","t"    ],
                        "R":["A","R","G","a","r","g"    ],
                        "Y":["Y","C","T","y","c","t"    ],
                        "B":["B","C","G","T","b","c","g","t"],
                        "D":["A","D","G","T","a","d","g","t"],
                        "H":["A","C","H","T","a","c","h","t"],
                        "V":["A","C","G","V","a","c","g","v"],
                        "N":["A","C","G","T","N","a","c","g","t","n"],
                        "Z":["Z"            ],

                        "a":["A","W","M","R","a","w","m","r"],
                        "c":["C","S","M","Y","c","s","m","y"],
                        "g":["G","S","K","R","g","s","k","r"],
                        "t":["U","T","W","K","Y","u","t","w","k","y"],
                        "u":["T","U","W","K","Y","t","u","w","k","y"],
                        "w":["A","T","W","a","t","w"        ],
                        "s":["S","C","G","s","c","g"    ],
                        "m":["M","A","C","m","a","c"    ],
                        "k":["K","G","T","k","g","t"    ],
                        "r":["A","R","G","a","r","g"    ],
                        "y":["Y","C","T","y","c","t"    ],
                        "b":["B","C","G","T","b","c","g","t"],
                        "d":["A","D","G","T","a","d","g","t"],
                        "h":["A","C","H","T","a","c","h","t"],
                        "v":["A","C","G","V","a","c","g","v"],
                        "n":["A","C","G","T","N","a","c","g","t","n"],
                        "z":["Z"            ]                  
                                            }

        self.gene_seq = gene_seq
        self.motifs   = motifs
        self.motif_cords = self.find_motifs()
        self.gene_num    = gene_num
        self.motif_colors = motif_colors


    def find_motifs(self):

        # gene motifs
        gene_motifs = {motif:[] for motif in self.motifs.keys()}

        # iterate through sequence  
        for i in range(len(self.gene_seq)):

            # slice into each length chunk
            for motif in self.motifs.keys():

                # check to make sure slice won't fall off end of list
                # then slice and add coordinates if its a motif
                end_slice = i + self.motifs[motif]

                if end_slice <= len(self.gene_seq) + 1:
                    seq_slice = self.gene_seq[i:end_slice]

                    if self.compare_sequences(seq_slice, motif):
                        gene_motifs[motif].append((i, end_slice - 1))

        return gene_motifs

    # compare 2 sequences using iupac dictionary       
    def compare_sequences(self, seq1, seq2):

        for base1,base2 in zip(seq1,seq2):
            
            if base1 in self.iupac_dic[base2]:
                match = True

            else:
                match = False
                break

        return match 

    def draw(self, context):

        # draw motifs
        #context.scale (10000, 10000)
        for motif in self.motifs.keys():
                
                #context.set_source_rgb(52, 235, 207)
                red     = self.motif_colors[motif][0]
                green   = self.motif_colors[motif][1]
                blue    = self.motif_colors[motif][2]

                for coord in self.motif_cords[motif]:

                        # coord = (start, stop)
                        m0 = coord[0]
                        m1 = coord[1]
                        
                        # draw line at motif
                        context.set_line_width(15)
                        context.set_source_rgba(red, green, blue, .5)
                        context.move_to(m0+margin, self.gene_num)
                        context.line_to(m1+margin, self.gene_num)
                        context.stroke()



'''
FASTA HEADER CLASS
'''

class FASTA:

    def __init__(self, header, gene_num):

        self.header = header
        self.gene_num    = gene_num

    def draw(self, context):
        # draw gene title
        context.set_source_rgb(0, 0, 0)
        context.set_font_size(15)
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) 
        context.move_to(margin, self.gene_num-25)
        context.show_text(self.header)
        context.stroke()


'''
GENE GROUP CLASS
'''
class GENE_GROUP:

    def __init__(self, fasta_header, gene_obj, exon_obj, motif_obj):

        self.fasta  = fasta_header
        self.gene   = gene_obj
        self.exon   = exon_obj
        self.motif  = motif_obj

    def draw_all(self, context):

        self.fasta.draw(context)
        self.gene.draw(context)
        self.exon.draw(context)
        self.motif.draw(context)


'''
MAIN
'''

def main():

        # make motif dictionary 
        motif_dic = {}
        motif_fh = open(motif_file, "r")
        for line in motif_fh:
                line = line.strip()
                motif_dic[line] = len(line)

        # define motif colors
        motif_colors = {motif:(random.random(), random.random(), random.random()) for motif in motif_dic.keys()}

        # go through fasta file and send each read to fasta_process
        # open file, read and store header
        fasta_fh = open(fasta_file, "r")
        fasta_genes = []
        gene_num = 100



        # store the rest of the fasta read as a string
        fasta_seq = ""
        header = fasta_fh.readline().strip()
        fasta_header = FASTA(header, gene_num)

        for line in fasta_fh:
                #print(line)
                line = line.strip()

                if line[0] == ">":
                        # make 4 objects
                        gene_obj       = GENE(fasta_seq, gene_num)
                        exon_obj       = EXON(fasta_seq, gene_num)
                        motif_obj      = MOTIF(fasta_seq, gene_num, motif_dic, motif_colors)

                        # combine 5 objects
                        gene_group_obj = GENE_GROUP(fasta_header, gene_obj, exon_obj, motif_obj)
                        
                        # add object to list
                        fasta_genes.append(gene_group_obj)

                        # reset all
                        fasta_seq       = ""
                        gene_num += 100
                        fasta_header = FASTA(line, gene_num)
                
                else:
                        fasta_seq = fasta_seq + line

        # make 4 objects
        gene_obj       = GENE(fasta_seq, gene_num)
        exon_obj       = EXON(fasta_seq, gene_num)
        motif_obj      = MOTIF(fasta_seq, gene_num, motif_dic, motif_colors)

        # combine 5 objects
        gene_group_obj = GENE_GROUP(fasta_header, gene_obj, exon_obj, motif_obj)

        # add gene group to list
        fasta_genes.append(gene_group_obj)


        # draw genes
        # initialize surface
        surface = cairo.SVGSurface(svg_name, 1000, 1000)
        context = cairo.Context(surface)

        # make a figure title
        context.set_source_rgb(0, 0, 0)
        context.set_font_size(30)
        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) 
        context.move_to(margin, 40)
        context.show_text("Motif Finder Results")
        context.stroke()

        for gene in fasta_genes:
            gene.draw_all(context)
        
        surface.finish()


main()