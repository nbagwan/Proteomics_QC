__author__ = 'Navratan Bagwan'

## this program takes the PTM-expand file nd removes duplicate and checks all the phospho site in mouse database and tries to extract the
## site not provides by Maxquant ####
## it also check every site in database and report the position in protein.
## if the program can not extract the all pospho position from a peptide, it extracts the evidence peptides from vidence file

import pdb
from collections import defaultdict
from Bio import SeqIO
import re
import QC_stats


### method takes 3 input, expand table you created in perseus, evidence table by MQ,
### the protein fasta file you used for the MQ search (for info check the README file)
def validate(phosphoEXPAND_table, evidenceTable, fastaDB):


    ### evidence being parsed and saved as dictionary
    ### structure:  multiplicity + protein id = all the peptides with phospho
    evidence_Dic = QC_stats.evidenceTableParse(evidenceFile=evidenceTable)

    ### parsing the fasta file and saved as dictionary
    fastDB_dic = QC_stats.parseFasta(fastaQ= fastaDB)
    # pdb.set_trace()
    baseDic = {}
    with open(phosphoEXPAND_table) as expandFile:
        header_names = str(next(expandFile))
        ### hear i am trying to access all the coloum which has the these (approximatly) names ####
        ### check find_col_index in the stasts file
        intensity_col = QC_stats.find_col_index(string= header_names, substring="Reporter intensity")
        reverse_col = QC_stats.find_col_index(string=header_names, substring="Reverse")
        contaminant_col = QC_stats.find_col_index(string=header_names, substring="contaminant")
        leading_protein_col = QC_stats.find_col_index(string=header_names, substring="Leading protein")
        leading_position_col = QC_stats.find_col_index(string=header_names, substring="Positions within proteins")
        multiplicity_col = QC_stats.find_col_index(string=header_names, substring="Multiplicity")
        prob_seq_col = QC_stats.find_col_index(string=header_names, substring="Phospho (STY) Probabilities")


        ###here i am creating an output header. two new coloumn plus the header in your input file
        header =  "Evidence Sequences" + "\t" + "Cross_Check_DB_sites" +  "\t" + header_names.strip()
        # pdb.set_trace()
        print(header.strip())

        for line in expandFile:
            if line != "\n":
                splits = line.split("\t")
                ### filtering the Contanminant and reverse peptide
                if splits[contaminant_col[0]].strip() != "+" and splits[reverse_col[0]].strip() != "+":
                    # pdb.set_trace()
                    #### filtering the peptides where all the intensities values are zero
                    Intensity_sum = sum(map(float, splits[:max(intensity_col) + 1]))
                    if Intensity_sum != 0.0:
                        if "REV__" not in splits[leading_protein_col[0]].split(";")[0].strip():
                            uniprot_accession = splits[leading_protein_col[0]].split(";")[0].strip()
                            Unique_ID = str(uniprot_accession.strip()) + "__" + str(Intensity_sum)
                            if Unique_ID not in baseDic:
                                baseDic[Unique_ID] = [[line.strip()],[splits[leading_position_col[0]].split(";")[0].strip()]]
                            else:
                                baseDic[Unique_ID][1].append(splits[leading_position_col[0]].split(";")[0].strip())
                                ### this basedic contains now all the peptide with intensities without intensity duplicate


    for i in baseDic:
        if "REV_" not in baseDic[i][0][0].split("\t")[leading_protein_col[0]].split(";")[0].strip(): #### leading protein
            accesion = baseDic[i][0][0].split("\t")[leading_protein_col[0]].split(";")[0].strip() ## leading protein
            prob_seq = baseDic[i][0][0].split("\t")[prob_seq_col[0]].strip()
            plain_seq = re.sub("[^a-zA-Z]+", "", prob_seq)

            prob_String = "".join(c for c in prob_seq if not c.isalpha())
            positions_list = QC_stats.positionList(stringInput=prob_String)
            multiplicity = baseDic[i][0][0].split("\t")[multiplicity_col[0]].strip()
            position_in_protein = fastDB_dic[accesion].index(plain_seq)
            PTM_pos_in_seq = prob_seq.index(str(max(positions_list)))
            # temp_pos = prob_seq.index(str(max(positions_list)))
            PTM_pos_in_seq_length = QC_stats.singleSTY_cout(stringInput=prob_seq[0:int(PTM_pos_in_seq)])

            ### for mono phospho peptides, we are just cross chekcing the position in fasta
            if multiplicity == "___1":

                print("NA", "\t", str(position_in_protein + PTM_pos_in_seq_length).strip(), "\t", baseDic[i][0][0].strip())

            elif multiplicity == "___2":
                positions_list_filtered = [n for n in positions_list if float(n) > 0.4]
                # pdb.set_trace()
                if len(positions_list_filtered) == 1:
                    PTM_pos_in_seq_1 = prob_seq.index(str(max(positions_list_filtered)))
                    PTM_pos_in_seq_length_1 = QC_stats.singleSTY_cout(stringInput=prob_seq[0:int(PTM_pos_in_seq_1)])
                    # pdb.set_trace()
                    evidence_seq = QC_stats.check_evidence_check(id="2_" + str(accesion), evi_dict=evidence_Dic, sequence=plain_seq)
                    # pdb.set_trace()
                    print(evidence_seq.strip(), "\t", str(position_in_protein + PTM_pos_in_seq_length_1).strip(), "\t", baseDic[i][0][0].strip())

                else:
                    post_list = QC_stats.peptide_position(substring=prob_seq, position_in_protein=position_in_protein)
                    print("NA", "\t", ";".join(post_list), "\t", baseDic[i][0][0].strip())


            else:
                post_list = QC_stats.peptide_position(substring=prob_seq, position_in_protein=position_in_protein)
                print("NA", "\t", ";".join(post_list), "\t", baseDic[i][0][0].strip())




validate(phosphoEXPAND_table=r"O:\CardiacProteomics\Scripts\NB_Scripts_and Data\Example_data_phospho\Phospho (STY)Sites-PTM-expand.txt",
         evidenceTable=r"O:\CardiacProteomics\Scripts\NB_Scripts_and Data\Example_data_phospho\evidence.txt",
         fastaDB=r"O:\CardiacProteomics\Scripts\NB_Scripts_and Data\Example_data_phospho\2019_12_03_Uniprot_reviewed_Mouse.fasta")