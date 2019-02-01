# program made by Nikola Vinko
# contact : unispieler@gmail.com
# _______________________________________________________________________#

# Clonal Rearrangement Visualisation Tool helps you to visualise IGH Fr1, IGH Fr2, IGK and TRG Clonal Rearrangements
#Copyright (C) GPL v3 

#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.


import matplotlib
matplotlib.use('Agg')
import os
import sys
import plotly
from plotly.offline import download_plotlyjs, init_notebook_mode, plot
from plotly.graph_objs import Bar, Figure, Layout
import pandas as pd
import plotly.graph_objs as go
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import pandas as pd
import pandas
import re


def C_graph_TCRG(folder_name, loc, loc_save, region, patient_id_1, patient_id_2, region_graph_1, region_graph_2, srr, src):
    print("Running Comparison TCRG graph for patient: " + patient_id_1 + " and " + patient_id_2)

	# Creating folder : "Comparison"
    if not os.path.exists(folder_name):  
        os.makedirs(folder_name)

	# Creating folder : "CSV_Files"
    if not os.path.exists(loc_save):  
        os.makedirs(loc_save)

	# error checking
    tsv_file_1 = loc + region + region_graph_1  
    if not os.path.exists(tsv_file_1):
        print("ERROR: file not found: " + tsv_file_1)

    print("Running TCRG for file_1: " + tsv_file_1)

    try:
        csv_table = pandas.read_table(tsv_file_1, sep="\t")
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")

	# saving the CSV file from file_1
    try:
        print("Saving data to: " + loc_save + 'index' + patient_id_1 + ".fastq_read_summary.csv")  

        csv_table.to_csv(loc_save + 'index' + patient_id_1 + ".fastq_read_summary.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    if not os.path.exists(loc_save):
        os.makedirs(loc_save)

    tsv_file_2 = loc + region + region_graph_2
    if not os.path.exists(tsv_file_2):
        print("ERROR: file not found: " + tsv_file_2)

    print("Running TCRG for file_2: " + tsv_file_2)

    try:
        csv_table = pandas.read_table(tsv_file_2, sep="\t")
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")
	
	# saving the CSV file from file_2
    try:
        print("Saving data to: " + loc_save + 'index' + patient_id_2 + ".fastq_read_summary.csv")

        csv_table.to_csv(loc_save + 'index' + patient_id_2 + ".fastq_read_summary.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

	# selecting the first csv file
    with open(loc_save + "index" + patient_id_1 + ".fastq_read_summary.csv") as f:  
        fastq_data = pd.read_csv(f, header=None)
		# selecting the %total_reads and fastq sequences column
        total_reads1 = fastq_data.iloc[srr:, src]
        seq1 = fastq_data.iloc[1:, 4]  

	# selecting the second csv file
    with open(loc_save + "index" + patient_id_2 + ".fastq_read_summary.csv") as f:  
        fastq_data = pd.read_csv(f, header=None)
		# selecting the %total_reads and fastq sequences column
        total_reads2 = fastq_data.iloc[srr:, src] 
        seq2 = fastq_data.iloc[1:, 4]

	# creating a dictionary that will contain both %total read values of merged fastqs
    dic1 = {}  
    for i in range(1, len(seq1) + 1):
		# adding the %total reads from first sample (1_1)
        dic1[seq1[i]] = [float(total_reads1[i])]  

    for i in range(1, len(seq2) + 1):
		# if the sequence from the second sample has a mathc in the first sample :
        if seq2[i] in dic1:  
			# then add that second %total read value  to that corresponding sequence key
            dic1[seq2[i]].append(float(total_reads2[i]))  
		
		# if it doesnt match, then create a new sequence key with list value : [0,value]
        else:
            dic1[seq2[i]] = [0, float(total_reads2[i])]  
	# appending to every key one more 0 (just to make sure that all keys have at minimum 2 values (x,y))
    for key in dic1:
        dic1[key].append(0)  

    x = []
    y = []

	# adding the first values from each dic1[key] to x that will represent x axis values
    for key in dic1:
        x.append(dic1[key][0])  

	# adding the second values from each dic1[key] to y that will represent y axis values
    for key in dic1:
        y.append(dic1[key][1])  

    plt.scatter(x, y, s=1)
    plt.title("Comparison graph TRG")
    plt.savefig("./plots/" + patient_id_1[0:3] + "_Comparison/" + patient_id_1 + "_vs_" + patient_id_2 + "_Comparison_TRG.png", dpi=200)

    print("TRG Comparison Done!" + "\n")
    plt.close()


def C_graph_IGH_FR1(folder_name, loc, loc_save, region, patient_id_1, patient_id_2, region_graph_1, region_graph_2, srr, src):
    print("Running Comparison TCRG graph for patient: " + patient_id_1 + " and " + patient_id_2)

	# Creating folder : "Comparison"
    if not os.path.exists(folder_name):  
        os.makedirs(folder_name)

	# Creating folder : "CSV_Files"
    if not os.path.exists(loc_save):
        os.makedirs(loc_save)

	# error checking
    tsv_file_1 = loc + region + region_graph_1
    if not os.path.exists(tsv_file_1):
        print("ERROR: file not found: " + tsv_file_1)

    print("Running IGH_FR1 for file: " + tsv_file_1)

    try:
        csv_table = pandas.read_table(tsv_file_1, sep="\t")
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")

	# saving the CSV file from file_1
    try:
        print("Saving data to: " + loc_save + 'index' + patient_id_1 + ".fastq_read_summary_family.csv")

        csv_table.to_csv(loc_save + 'index' + patient_id_1 + ".fastq_read_summary_family.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    if not os.path.exists(loc_save):
        os.makedirs(loc_save)

    tsv_file_2 = loc + region + region_graph_2
    if not os.path.exists(tsv_file_2):
        print("ERROR: file not found: " + tsv_file_2)

    print("Running IGK_FR1 for file: " + tsv_file_2)

    try:
        csv_table = pandas.read_table(tsv_file_2, sep="\t")
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")

	# saving the CSV file from file_2
    try:
        print("Saving data to: " + loc_save + 'index' + patient_id_2 + ".fastq_read_summary_family.csv")

        csv_table.to_csv(loc_save + 'index' + patient_id_2 + ".fastq_read_summary_family.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

	# opening the first csv file
    with open(loc_save + "index" + patient_id_1 + ".fastq_read_summary_family.csv") as f:
        fastq_data = pd.read_csv(f, header=None)
		# selecting the %total_reads and fastq sequences column
        total_reads1 = fastq_data.iloc[srr:, src]
        seq1 = fastq_data.iloc[1:, 4]
	# opening hte second csv file
    with open(loc_save + "index" + patient_id_2 + ".fastq_read_summary_family.csv") as f:  
        fastq_data = pd.read_csv(f, header=None)
		# selecting the %total_reads and fastq sequences column
        total_reads2 = fastq_data.iloc[srr:, src]
        seq2 = fastq_data.iloc[1:, 4]

	# creating a dictionary that will contain both %total read values of merged fastqs
    dic1 = {}

	# adding the %total reads from first sample (1_1)
    for i in range(1, len(seq1) + 1):
        dic1[seq1[i]] = [float(total_reads1[i])]  

    for i in range(1, len(seq2) + 1):
		# if the sequence from the second sample has a mathc in the first sample :
        if seq2[i] in dic1:
			# then add that second %total read value  to that corresponding sequence key
            dic1[seq2[i]].append(float(total_reads2[i]))
			# if it doesnt match, then create a new sequence key with list value : [0,value]
        else:
            dic1[seq2[i]] = [0, float(total_reads2[i])]
	# appending to every key one more 0 (just to make sure that all keys have at minimum 2 values (x,y))
    for key in dic1:
        dic1[key].append(0)  

    x = []
    y = []
	# adding the first values from each dic1[key] to x that will represent x axis values
    for key in dic1:
        x.append(dic1[key][0])
	# adding the second values from each dic1[key] to y that will represent y axis values
    for key in dic1:
        y.append(dic1[key][1])

    plt.scatter(x, y, s=1)
    plt.title("Comparison graph IGH_FR1")
    plt.savefig("./plots/" + patient_id_1[0:3] + "_Comparison/" + patient_id_1 + "_vs_" + patient_id_2 + "_Comparison_IGH_FR1.png", dpi=200)

    print("IGH_FR1 Comparison Done!" + "\n")
    plt.close()


def C_graph_IGH_FR2(folder_name, loc, loc_save, region, patient_id_1, patient_id_2, region_graph_1, region_graph_2, srr, src):
    print("Running Comparison IGH_FR2 graph for patient: " + patient_id_1 + " and " + patient_id_2)

	# Creating folder : "Comparison"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

	# Creating folder : "CSV_Files"
    if not os.path.exists(loc_save):
        os.makedirs(loc_save)

	# error checking
    tsv_file_1 = loc + region + region_graph_1
    if not os.path.exists(tsv_file_1):
        print("ERROR: file not found: " + tsv_file_1)

    print("Running IGH_FR2 for file: " + tsv_file_1)

    try:
        csv_table = pandas.read_table(tsv_file_1, sep="\t")
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")

    try:
        print("Saving data to: " + loc_save + 'index' + patient_id_1 + ".fastq_read_summary_family.csv")

        csv_table.to_csv(loc_save + 'index' + patient_id_1 + ".fastq_read_summary_family.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    if not os.path.exists(loc_save):
        os.makedirs(loc_save)

    tsv_file_2 = loc + region + region_graph_2
    if not os.path.exists(tsv_file_2):
        print("ERROR: file not found: " + tsv_file_2)

    print("Running IGH_FR2 for file: " + tsv_file_2)

    try:
        csv_table = pandas.read_table(tsv_file_2, sep="\t")
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")

    try:
        print("Saving data to: " + loc_save + 'index' + patient_id_2 + ".fastq_read_summary_family.csv")

        csv_table.to_csv(loc_save + 'index' + patient_id_2 + ".fastq_read_summary_family.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

    # opening the first csv file
    with open(loc_save + "index" + patient_id_1 + ".fastq_read_summary_family.csv") as f:
        fastq_data = pd.read_csv(f, header=None)
		# selecting the %total_reads and fastq sequences column
        total_reads1 = fastq_data.iloc[srr:, src]
        seq1 = fastq_data.iloc[1:, 4]
	# opening hte second csv file
    with open(loc_save + "index" + patient_id_2 + ".fastq_read_summary_family.csv") as f:  
        fastq_data = pd.read_csv(f, header=None)
		# selecting the %total_reads and fastq sequences column
        total_reads2 = fastq_data.iloc[srr:, src]
        seq2 = fastq_data.iloc[1:, 4]

	# creating a dictionary that will contain both %total read values of merged fastqs
    dic1 = {}

	# adding the %total reads from first sample (1_1)
    for i in range(1, len(seq1) + 1):
        dic1[seq1[i]] = [float(total_reads1[i])]  

    for i in range(1, len(seq2) + 1):
		# if the sequence from the second sample has a mathc in the first sample :
        if seq2[i] in dic1:
			# then add that second %total read value  to that corresponding sequence key
            dic1[seq2[i]].append(float(total_reads2[i]))
			# if it doesnt match, then create a new sequence key with list value : [0,value]
        else:
            dic1[seq2[i]] = [0, float(total_reads2[i])]
	# appending to every key one more 0 (just to make sure that all keys have at minimum 2 values (x,y))
    for key in dic1:
        dic1[key].append(0)  

    x = []
    y = []
	# adding the first values from each dic1[key] to x that will represent x axis values
    for key in dic1:
        x.append(dic1[key][0])
	# adding the second values from each dic1[key] to y that will represent y axis values
    for key in dic1:
        y.append(dic1[key][1])
		
    plt.scatter(x, y, s=1)
    plt.title("Comparison graph IGH_FR2")
    plt.savefig("./plots/" + patient_id_1[0:3] + "_Comparison/" + patient_id_1 + "_vs_" + patient_id_2 + "_Comparison_IGH_FR2.png", dpi=200)

    print("IGH_Fr2 Comparison Done!" + "\n")
    plt.close()


def C_graph_IGK(folder_name, loc, loc_save, region, patient_id_1, patient_id_2, region_graph_1, region_graph_2, srr, src):
    print("Running Comparison IGK graph for patient: " + patient_id_1 + " and " + patient_id_2)

    # Creating folder : "Comparison"
    if not os.path.exists(folder_name):  
        os.makedirs(folder_name)

	# Creating folder : "CSV_Files"
    if not os.path.exists(loc_save):  
        os.makedirs(loc_save)

	# error checking
    tsv_file_1 = loc + region + region_graph_1  
    if not os.path.exists(tsv_file_1):
        print("ERROR: file not found: " + tsv_file_1)

    print("Running IGK for file: " + tsv_file_1)

    try:
        csv_table = pandas.read_table(tsv_file_1, sep="\t")
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")

	# saving the CSV file from file_1
    try:
        print("Saving data to: " + loc_save + 'index' + patient_id_1 + ".fastq_read_summary_family.csv")

        csv_table.to_csv(loc_save + 'index' + patient_id_1 + ".fastq_read_summary_family.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    if not os.path.exists(loc_save):
        os.makedirs(loc_save)

    tsv_file_2 = loc + region + region_graph_2
    if not os.path.exists(tsv_file_2):
        print("ERROR: file not found: " + tsv_file_2)

    print("Running IGK for file: " + tsv_file_2)

	# saving the CSV file from file_2
    try:
        csv_table = pandas.read_table(tsv_file_2, sep="\t")
    except ValueError:
        print("Oops!  That was no valid number.  Try again...")
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data read done ")

    try:
        print("Saving data to: " + loc_save + 'index' + patient_id_2 + ".fastq_read_summary_family.csv")

        csv_table.to_csv(loc_save + 'index' + patient_id_2 + ".fastq_read_summary_family.csv", index=False)
    except:
        print("Oops! " + sys.exc_info()[0])

    print("Data save done !")

    # opening the first csv file
    with open(loc_save + "index" + patient_id_1 + ".fastq_read_summary_family.csv") as f:
        fastq_data = pd.read_csv(f, header=None)
		# selecting the %total_reads and fastq sequences column
        total_reads1 = fastq_data.iloc[srr:, src]
        seq1 = fastq_data.iloc[1:, 4]
	# opening hte second csv file
    with open(loc_save + "index" + patient_id_2 + ".fastq_read_summary_family.csv") as f:  
        fastq_data = pd.read_csv(f, header=None)
		# selecting the %total_reads and fastq sequences column
        total_reads2 = fastq_data.iloc[srr:, src]
        seq2 = fastq_data.iloc[1:, 4]

	# creating a dictionary that will contain both %total read values of merged fastqs
    dic1 = {}

	# adding the %total reads from first sample (1_1)
    for i in range(1, len(seq1) + 1):
        dic1[seq1[i]] = [float(total_reads1[i])]  

    for i in range(1, len(seq2) + 1):
		# if the sequence from the second sample has a mathc in the first sample :
        if seq2[i] in dic1:
			# then add that second %total read value  to that corresponding sequence key
            dic1[seq2[i]].append(float(total_reads2[i]))
			# if it doesnt match, then create a new sequence key with list value : [0,value]
        else:
            dic1[seq2[i]] = [0, float(total_reads2[i])]
	# appending to every key one more 0 (just to make sure that all keys have at minimum 2 values (x,y))
    for key in dic1:
        dic1[key].append(0)  

    x = []
    y = []
	# adding the first values from each dic1[key] to x that will represent x axis values
    for key in dic1:
        x.append(dic1[key][0])
	# adding the second values from each dic1[key] to y that will represent y axis values
    for key in dic1:
        y.append(dic1[key][1])
		
    plt.scatter(x, y, s=1)
    plt.title("Comparison graph IGK")
    plt.savefig("./plots/" + patient_id_1[0:3] + "_Comparison/" + patient_id_1 + "_vs_" + patient_id_2 + "_Comparison_IGK.png", dpi=200)

    print("IGK Comparison Done!" + "\n")
    plt.close()


def main(argv):
    basepath = os.getcwd()

    if not os.path.exists(basepath):
        print("ERROR: Path \"" + basepath + "\" does not exist!")
        exit(1)

    os.chdir(basepath)

    mypath = "./TCRG_output_v1.3.3_automated"
    if os.path.exists(mypath):
        onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

        # get all patient_ids in the folder
        lista_id = []
        for file in onlyfiles:
            match = re.search(r"index([0-9]+_[0-9]+)\.", file, re.IGNORECASE)
            if match:
                lista_id.append(match.group(1))

        lista_id = set(lista_id)

        for patient_id_1 in lista_id:
            for patient_id_2 in lista_id:
                print("Working on: " + patient_id_1 + " <-> " + patient_id_2)

                # PARAMETERS for TCRG
				
				# location of the tsv files
                loc = "." 
				# location where to save the csv files
                loc_save = "./CSV_files_TCRG/"  
                region = "/TCRG_output_v1.3.3_automated"
                folder_name = "./plots/" + patient_id_1[0:3] + "_Comparison"
                region_graph_1 = "/index" + patient_id_1 + ".fastq_read_summary.tsv"
                region_graph_2 = "/index" + patient_id_2 + ".fastq_read_summary.tsv"
                srr = 1
                src = 9

                C_graph_TCRG(folder_name, loc, loc_save, region, patient_id_1, patient_id_2, region_graph_1, region_graph_2, srr, src)

    mypath = "./IGH_FR1_output_v1.3.3_automated"
    if os.path.exists(mypath):
        onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

        # get all patient_ids in the folder
        lista_id = []
        for file in onlyfiles:
            match = re.search(r"index([0-9]+_[0-9]+)\.", file, re.IGNORECASE)
            if match:
                lista_id.append(match.group(1))

        lista_id = set(lista_id)

        for patient_id_1 in lista_id:
            for patient_id_2 in lista_id:
                print("Working on: " + patient_id_1 + " <-> " + patient_id_2)

                # PARAMETERS for IGH FR1
				
				# location of the tsv files
                loc = "."  
                loc_save = "./CSV_files_IGH_FR1/"
				# location where to save the csv files
                region = "/IGH_FR1_output_v1.3.3_automated"
                folder_name = "./plots/" + patient_id_1[0:3] + "_Comparison"
                region_graph_1 = "/index" + patient_id_1 + ".fastq_read_summary_family.tsv"
                region_graph_2 = "/index" + patient_id_2 + ".fastq_read_summary_family.tsv"
                srr = 1
                src = 9

                C_graph_IGH_FR1(folder_name, loc, loc_save, region, patient_id_1, patient_id_2, region_graph_1, region_graph_2, srr, src)

    mypath = "./IGH_FR2_output_v1.3.3_automated"
    if os.path.exists(mypath):
        onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

        # get all patient_ids in the folder
        lista_id = []
        for file in onlyfiles:
            match = re.search(r"index([0-9]+_[0-9]+)\.", file, re.IGNORECASE)
            if match:
                lista_id.append(match.group(1))

        lista_id = set(lista_id)

        for patient_id_1 in lista_id:
            for patient_id_2 in lista_id:
                print("Working on: " + patient_id_1 + " <-> " + patient_id_2)
                # PARAMETERS for IGH FR2

				# location of the tsv files
                loc = "."  
				# location where to save the csv files
                loc_save = "./CSV_files_IGH_FR2/"  
                region = "/IGH_FR2_output_v1.3.3_automated"
                folder_name = "./plots/" + patient_id_1[0:3] + "_Comparison"
                region_graph_1 = "/index" + patient_id_1 + ".fastq_read_summary_family.tsv"
                region_graph_2 = "/index" + patient_id_2 + ".fastq_read_summary_family.tsv"
                srr = 1
                src = 9

                C_graph_IGH_FR2(folder_name, loc, loc_save, region, patient_id_1, patient_id_2, region_graph_1, region_graph_2, srr, src)

    mypath = "./IGK_output_v1.3.3_automated"
    if os.path.exists(mypath):
        onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

        # get all patient_ids in the folder
        lista_id = []
        for file in onlyfiles:
            match = re.search(r"index([0-9]+_[0-9]+)\.", file, re.IGNORECASE)
            if match:
                lista_id.append(match.group(1))

        lista_id = set(lista_id)

        for patient_id_1 in lista_id:
            for patient_id_2 in lista_id:
                print("Working on: " + patient_id_1 + " <-> " + patient_id_2)

                # PARAMETERS for IGK

				# location of the tsv files
                loc = "."
				# location where to save the csv files
                loc_save = "./CSV_files_IGK/"
                region = "/IGK_output_v1.3.3_automated"
                folder_name = "./plots/" + patient_id_1[0:3] + "_Comparison"
                region_graph_1 = "/index" + patient_id_1 + ".fastq_read_summary_family.tsv"
                region_graph_2 = "/index" + patient_id_2 + ".fastq_read_summary_family.tsv"
                srr = 1
                src = 9

                C_graph_IGK(folder_name, loc, loc_save, region, patient_id_1, patient_id_2, region_graph_1, region_graph_2, srr, src)

                


if __name__ == "__main__":
    main(sys.argv[1:])
