__author__ = 'Navratan Bagwan'


import pdb
import re
from Bio import SeqIO
import numpy as np
import math
from difflib import get_close_matches


def closeMatches(patterns, word):

    return(get_close_matches(word, patterns))


def trim_seq_window(window, length):
    if len(window) == length:
        pass
    else:
        limit = int((len(window) - length) / 2)
        window_trimmed = window[limit:-limit]

    return window_trimmed

def getList(fastaQ):
    listout = []
    for seq_record in SeqIO.parse(fastaQ, "fasta"):
        if "GN=" in seq_record.description:
            if seq_record.description.split("GN=")[1].split(" ")[0].strip() not in listout:
                listout.append(seq_record.description.split("GN=")[1].split(" ")[0].strip())

    return listout


def parseFasta(fastaQ):
    FastaDic = {}
    ## parsing fasta accroding to gene name
    for seq_record in SeqIO.parse(fastaQ, "fasta"):
        # print(seq_record.id)
        # print(seq_record.seq)
        sequence = str(seq_record.seq)
        if "GN=" in seq_record.description:
            if seq_record.description.split("GN=")[1].split(" ")[0].strip() not in FastaDic:
                FastaDic[seq_record.description.split("GN=")[1].split(" ")[0].strip()] = str(sequence)

    return FastaDic

def parseFasta_acc(fastaQ):
    ##parsing fasta file accroding to accesion(protein id)
    FastaDic = {}
    for seq_record in SeqIO.parse(fastaQ, "fasta"):
        # print(seq_record.id)
        # print(seq_record.seq)
        sequence = str(seq_record.seq)
        if "GN=" in seq_record.description:
            if seq_record.description.split("GN=")[1].split(" ")[0].strip() not in FastaDic:
                FastaDic[seq_record.description.split("GN=")[1].split(" ")[0].strip()] = seq_record.description.split("|")[1].strip()
    return FastaDic


def createDic(file, key, species, species_col):
    species_dict = {}
    with open(file) as dict_file:
        next(dict_file)
        for line in dict_file:
            if line != "\n":
                splits = line.split("\t")
                # pdb.set_trace()
                if splits[species_col - 1].strip() == species.strip():
                    if splits[key - 1].strip() not in species_dict:
                        species_dict[splits[key - 1].strip()] = line.strip()
    return species_dict


#### this is the main method, takes phosphosite plus file, human and mice fasta and the input file
#### exmaple of all these files are in the exmaple folder
def ortholog_peptides(psp_file, mouse_fasta, human_fasta, input_data):

    ### creating the dictionary of phosph site plus data data. once for the human and once for mouse
    mouse_psp_dict = createDic(file=psp_file, key=5, species="mouse", species_col=6)
    human_psp_dict = createDic(file=psp_file, key=5, species="human", species_col=6)

    ### dictionary for teh fasta file
    human_fasta_dict = parseFasta(fastaQ=human_fasta)
    mouse_fasta_dict = parseFasta(fastaQ=mouse_fasta)

    human_fasta_dict_acc = parseFasta_acc(fastaQ=human_fasta)
    human_fasta_geneList = getList(fastaQ=human_fasta)
    # pdb.set_trace()

    count = 0
    NOT_count = 0
    with open(input_data) as inputFile:
        header = next(inputFile)
        print("Source", "\t", "Fol changeC", "\t", header.strip(), "\t", "site (mouse)", "\t", "gene name (mouse)", "\t",
              "seq window 7+/-7 (mouse)", "\t", "Uniprot_ID (mouse)", "\t", "site (human)", "\t", "gene name (human)", "\t",
              "seq window 7+/-7 (human)", "\t", "Uniprot_ID (human)")
        for line in inputFile:
            if line != "\t":
                splits = line.split("\t")
                ## phospho site position in protein coloumn.
                position_mono = splits[2].split(";")[0].strip()
                ## protein id coloumn
                accesion = splits[9].split(";")[0].strip()
                ## gene name coloumn
                gene_name = splits[6].split(";")[0].strip()
                ## amino acid coloum
                amino_acid = splits[0].strip()
                ### creating a palin peptiode sequence by remoing the probability from the sequence
                plain_seq = re.sub("[^a-zA-Z]+", "", splits[8].strip())
                ### sequence window
                seq_window = trim_seq_window(window=splits[7].strip(), length=15)
                seq_window_short = trim_seq_window(window=splits[7].strip(), length=7)
                FC = math.pow(2, float(splits[10].strip()))

                check = False
                for site_grp_id in mouse_psp_dict:
                    # pdb.set_trace()
                    if mouse_psp_dict[site_grp_id].split("\t")[2].strip() == accesion.strip():
                        aa_postion_mouse_psp = mouse_psp_dict[site_grp_id].split("\t")[3].split("-")[0][0]
                        position_mouse_psp = re.sub('[^0-9]', '',  mouse_psp_dict[site_grp_id].split("\t")[3].split("-")[0])

                        if position_mouse_psp == position_mono:
                            # pdb.set_trace()
                            if site_grp_id in human_psp_dict:
                                count = count + 1
                                # print(line.strip(), "\t", mouse_psp_dict[site_grp_id].strip(), "\t", human_psp_dict[site_grp_id].strip(), "\t", "PSP")
                                print("PSP", "\t", FC, "\t", line.strip(),
                                      "\t", mouse_psp_dict[site_grp_id].strip().split("\t")[3].split("-")[0].strip(), ## site
                                      "\t", mouse_psp_dict[site_grp_id].strip().split("\t")[0].strip(), ## gene
                                      "\t", mouse_psp_dict[site_grp_id].strip().split("\t")[7].strip(), "\t", accesion, ### seq window ### accesion (mouse)
                                      "\t", human_psp_dict[site_grp_id].strip().split("\t")[3].split("-")[0].strip(), ### site human
                                      "\t", human_psp_dict[site_grp_id].strip().split("\t")[0].strip(), "\t", ### gene
                                      human_psp_dict[site_grp_id].strip().split("\t")[7].strip(), "\t", human_psp_dict[site_grp_id].strip().split("\t")[2].strip())  ## seq window ## accesion (human)

                                check = True
                                break

                if not check:
                    # pdb.set_trace()
                    if gene_name.upper() in human_fasta_dict:
                        # pdb.set_trace()
                        # check_gene = False
                        if plain_seq in human_fasta_dict[gene_name.upper()]:
                            position_in_peptide = int(position_mono) - mouse_fasta_dict[gene_name].index(plain_seq)
                            position_in_human = human_fasta_dict[gene_name.upper()].index(plain_seq)
                            human_length = human_fasta_dict[gene_name.upper()].index(plain_seq)
                            mouse_length = mouse_fasta_dict[gene_name].index(plain_seq)

                            if len(human_fasta_dict[gene_name.upper()]) == len(mouse_fasta_dict[gene_name]):
                                if human_fasta_dict[gene_name.upper()][int(position_mono) - 1] == amino_acid:
                                    print("1-Exact match: same seq & postions", "\t", FC, "\t", line.strip(), "\t", str(amino_acid) + str(position_mono), "\t", gene_name.strip(), "\t", seq_window, "\t", accesion,
                                          "\t", str(amino_acid) + str(position_mono), "\t", gene_name.upper(), "\t", seq_window, "\t", human_fasta_dict_acc[gene_name.upper()].strip())

                                elif human_fasta_dict[gene_name.upper()][(position_in_human + position_in_peptide) - 1] == amino_acid:
                                    print("2-Exact match: same seq & diff postions","\t", FC, "\t", line.strip(), "\t", str(amino_acid) + str(position_mono), "\t", gene_name.strip(), "\t", seq_window, "\t", accesion,
                                          "\t", str(amino_acid) + str(position_in_human + position_in_peptide), "\t", gene_name.upper(), "\t", seq_window,  "\t", human_fasta_dict_acc[gene_name.upper()].strip())

                            elif human_fasta_dict[gene_name.upper()][(position_in_human + position_in_peptide) - 1] == amino_acid:
                                print("3-Exact match: same seq & diff postions", "\t", FC, "\t", line.strip(), "\t", str(amino_acid) + str(position_mono), "\t", gene_name.strip(), "\t", seq_window, "\t", accesion,
                                      "\t", str(amino_acid) + str(position_in_human + position_in_peptide), "\t", gene_name.upper(), "\t", seq_window, "\t", human_fasta_dict_acc[gene_name.upper()].strip())

                        elif seq_window.strip() in human_fasta_dict[gene_name.upper()]:
                            position_in_peptide = int(position_mono) - mouse_fasta_dict[gene_name].index(plain_seq)
                            position_in_human = human_fasta_dict[gene_name.upper()].index(seq_window)
                            # human_length = human_fasta_dict[gene_name.upper()].index(plain_seq)
                            mouse_length = mouse_fasta_dict[gene_name].index(plain_seq)

                            if human_fasta_dict[gene_name.upper()][int(position_mono) - 1] == amino_acid:
                                print("4-seq-window match: same position", "\t", FC, "\t", line.strip(), "\t", str(amino_acid) + str(position_mono), "\t", gene_name.strip(), "\t", seq_window, "\t", accesion,
                                      "\t", str(amino_acid) + str(position_mono), "\t", gene_name.upper(), "\t", seq_window.strip(), "\t", human_fasta_dict_acc[gene_name.upper()].strip())

                            elif human_fasta_dict[gene_name.upper()][int(position_in_human + 7)] == amino_acid:
                                print("5-seq-window match: diff position", "\t", FC, "\t", line.strip(), "\t", str(amino_acid) + str(position_mono), "\t", gene_name.strip(), "\t", seq_window, "\t", accesion,
                                      "\t", str(amino_acid) + str(int(position_in_human + 7)), "\t", gene_name.upper(), "\t", seq_window.strip(), "\t", human_fasta_dict_acc[gene_name.upper()].strip())
                                # print(line.strip())

                        else:
                            print("No solutions found", "\t", FC, "\t", line.strip())

                    else:
                        possiblematch = closeMatches(patterns=human_fasta_geneList, word=gene_name.upper())
                        if len(possiblematch) == 0:
                            print("gene is not in humans/or no simmiar gene found", "\t", FC, "\t", line.strip())
                        else:
                            # pdb.set_trace()
                            if plain_seq in human_fasta_dict[possiblematch[0]] and gene_name in mouse_fasta_dict:

                                position_in_peptide = int(position_mono) - mouse_fasta_dict[gene_name].index(plain_seq)
                                position_in_human = human_fasta_dict[possiblematch[0]].index(plain_seq)
                                if len(human_fasta_dict[possiblematch[0]]) == len(mouse_fasta_dict[gene_name]):
                                    if human_fasta_dict[possiblematch[0]][int(position_mono) - 1] == amino_acid:
                                        print("NEW-1-Exact match: same seq & postions", "\t", FC, "\t", line.strip(), "\t", str(amino_acid) + str(position_mono), "\t", gene_name.strip(), "\t",
                                              seq_window, "\t", accesion, "\t", str(amino_acid) + str(position_mono), "\t", gene_name.upper(), "\t", seq_window, "\t", human_fasta_dict_acc[possiblematch[0]].strip())
                                    else:
                                        print("T-1")

                                else:
                                    print("T-2")

                            elif plain_seq in human_fasta_dict[possiblematch[0]] and possiblematch[0].capitalize() in mouse_fasta_dict:
                                position_in_peptide = int(position_mono) - mouse_fasta_dict[possiblematch[0].capitalize()].index(plain_seq)
                                position_in_human = human_fasta_dict[possiblematch[0]].index(plain_seq)

                                if len(human_fasta_dict[possiblematch[0]]) == len(mouse_fasta_dict[possiblematch[0].capitalize()]):
                                    if human_fasta_dict[possiblematch[0]][int(position_mono) - 1] == amino_acid:
                                        print("NEW-1-Exact match: same seq & postions", "\t", FC, "\t", line.strip(), "\t", str(amino_acid) + str(position_mono), "\t", gene_name.strip(), "\t",
                                              seq_window, "\t", accesion, "\t", str(amino_acid) + str(position_mono), "\t", gene_name.upper(), "\t", seq_window, "\t", human_fasta_dict_acc[possiblematch[0]].strip())

                                    elif human_fasta_dict[possiblematch[0]][ (position_in_human + position_in_peptide) - 1] == amino_acid:
                                        print("NEW-2-Exact match: same seq & diff postions", "\t", FC, "\t", line.strip(), "\t", str(amino_acid) + str(position_mono), "\t", gene_name.strip(), "\t", seq_window, "\t", accesion,
                                              "\t", str(amino_acid) + str(position_in_human + position_in_peptide), "\t", gene_name.upper(), "\t", seq_window, "\t",human_fasta_dict_acc[possiblematch[0]].strip())

                                    else:
                                        print("T-3")

                                elif human_fasta_dict[possiblematch[0]][(position_in_human + position_in_peptide) - 1] == amino_acid:
                                        print("NEW-3-Exact match: same seq & diff postions", "\t", FC, "\t", line.strip(), "\t", str(amino_acid) + str(position_mono), "\t", gene_name.strip(), "\t", seq_window, "\t", accesion,

                                              "\t", str(amino_acid) + str(position_in_human + position_in_peptide),"\t", gene_name.upper(), "\t", seq_window, "\t", human_fasta_dict_acc[possiblematch[0]].strip())

                                else:
                                     print("T-4")
                            else:
                                print("T-5", "\t", FC, "\t", line.strip())




#### these parameters needs to changed
#### provide paths to these file as you see below in these example
### for all these file there is an example file in teh example folder.
ortholog_peptides(
    psp_file=r"P:\ARVC-PKP2\Mice\Left_ventricle\Phospho_proteome\Data_analysis\files\LV-Set1\KO-sed vs KO-exc --New\KSEA\Phosphorylation_site_dataset_2020_10_23_human_mouse.tsv",
    mouse_fasta=r"P:\ARVC-PKP2\Mice\Left_ventricle\Phospho_proteome\Data_analysis\files\LV-Set1\KO-sed vs KO-exc --New\KSEA\2019_12_03_Uniprot_reviewed_Mouse.fasta",
    human_fasta=r"P:\ARVC-PKP2\Mice\Left_ventricle\Phospho_proteome\Data_analysis\files\LV-Set1\KO-sed vs KO-exc --New\KSEA\2020_01_06_uniprot_reviewed_human.fasta",
    input_data=r"P:\Insulin\Liver\files\4-vv-in-insulin\WO_i3\Ortholog-input.txt")