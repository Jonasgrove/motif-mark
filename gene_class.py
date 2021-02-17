#!/usr/bin/env python

class Gene:

    def __init__(self, ORF, motifs, header):
        
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
        self.motifs             = motifs
        self.fasta_header       = header                        # header sequence
        self.gene_seq           = ORF                           # gene sequence
        self.motif_coordinates  = self.find_motifs()            # dictionary of motif coordinates ex {motif_sequence:(start,stop),(start,stop)}
        self.end                = len(self.gene_seq) - 1        # ORF end position
        self.exons              = self.find_exons()             # exon coordinates


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

    # compare 2 sequences using iupac dictionary       
    def compare_sequences(self,seq1, seq2):

        for base1,base2 in zip(seq1,seq2):
            
            if base1 in self.iupac_dic[base2]:
                match = True

            else:
                match = False
                break

        return match 

'''
sequence = "ACGTYYYAC"
motif = {"CCC":3}

obj = Gene(sequence, motif, "hello")

print(obj.motif_coordinates)
'''