#!/usr/bin/env python

# imports
import re
import gene_class
import argparse
import os
import copy
import random
import cairo
import argparse

Gene = gene_class.Gene

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

# svg filename 
svg_name = fasta_file.split(".fa")[0] + ".svg"

# draw a gene
def draw_a_gene(gene_objects, gene_nums, margin, motif_ext, motif_colors, svg_name):

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


        # iterate through all genes
        for gene_num,gene_obj in zip(gene_nums, gene_objects):
                # draw gene title
                context.set_source_rgb(0, 0, 0)
                context.set_font_size(15)
                context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) 
                context.move_to(margin, gene_num-25)
                context.show_text(gene_obj.fasta_header)
                context.stroke()

                # draw gene length
                context.set_line_width(3)
                context.move_to(margin, gene_num)
                context.line_to(gene_obj.end, gene_num)
                context.stroke()


                # draw exons
                context.set_line_width(15)
                for exon in gene_obj.exons:
                        e0 = exon[0]
                        e1 = exon[1]
                        context.move_to(e0 + margin, gene_num)
                        context.line_to(e1 + margin, gene_num)
                        context.stroke()

                # draw motifs
                #context.scale (10000, 10000)
                for motif in gene_obj.motif_coordinates.keys():
                        
                        #context.set_source_rgb(52, 235, 207)
                        red     = motif_colors[motif][0]
                        green   = motif_colors[motif][1]
                        blue    = motif_colors[motif][2]

                        for coord in gene_obj.motif_coordinates[motif]:

                                # coord = (start, stop)
                                m0 = coord[0]
                                m1 = coord[1]
                                
                                # draw line at motif
                                context.set_line_width(15)
                                context.set_source_rgba(red, green, blue, .5)
                                context.move_to(m0+margin, gene_num)
                                context.line_to(m1+margin, gene_num)
                                context.stroke()
        
        # make motif legend
        start_legend = max(gene_nums) + 40
        for motif in motif_colors.keys():

                # get rgb colors
                red     = motif_colors[motif][0]
                green   = motif_colors[motif][1]
                blue    = motif_colors[motif][2]

                #draw label
                context.set_source_rgb(0, 0, 0)
                context.set_font_size(15)
                context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) 
                context.move_to(margin+60, start_legend+3)
                context.show_text(motif)
                context.stroke()

                # draw line color
                # draw line at motif
                context.set_line_width(10)
                context.set_source_rgb(red, green, blue)
                context.move_to(margin, start_legend)
                context.line_to(margin + 50, start_legend)
                context.stroke()
                
                start_legend += 20


        surface.finish()



# main
'''
1. split fasta file up by fasta records
2. pass each 1-record fasta file into process_fasta() function
    i.      use regex to find all ORFs and make gene object from each
    ii.     within each gene object, use find_motif() method, which will
            create a dictionary of all motifs and store this as an
            attribute of this object
    iii.    now create vissualization information of fasta sequence which will 
            include annotations for: gene boundaries, introns, exons, and motifs 
3. use these different files to make svg image file     
'''
def main():

        # make motif dictionary 
        motif_dic = {}
        motif_fh = open(motif_file, "r")
        for line in motif_fh:
                line = line.strip()
                motif_dic[line] = len(line)

        # go through fasta file and send each read to fasta_process
        # open file, read and store header
        fasta_fh = open(fasta_file, "r")
        fasta_genes = []
        
        # store the rest of the fasta read as a string
        fasta_seq = ""
        header = fasta_fh.readline().strip()

        for line in fasta_fh:
                #print(line)
                line = line.strip()

                if line[0] == ">":
                        # self, ORF, start, end, motifs, header
                        gene_objects    = Gene(fasta_seq, motif_dic, header)
                        fasta_genes.append(gene_objects)
                        fasta_seq       = ""
                        header          = line
                
                else:
                        fasta_seq = fasta_seq + line
        
        # process last fasta read
        gene_objects    = Gene(fasta_seq, motif_dic, header)
        fasta_genes.append(gene_objects)

        # make pretty pictures from the gene objects  
        margin = 10
        motif_ext = 15
        gene_nums = [num for num in range(100, (len(fasta_genes)+1)*100, 100)]
        motif_colors = {motif:(random.random(), random.random(), random.random()) for motif in motif_dic.keys()}
        draw_a_gene(fasta_genes, gene_nums, margin, motif_ext, motif_colors, svg_name)
        



main()

