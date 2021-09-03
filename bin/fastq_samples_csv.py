'''This file snipped for nextflow loads in all the gz compressed files and generates a csv file in the right format matching samples and sapmple names
need to make sure that the files are named as NAME.{1,2}.fastq.gz
'''

import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description="Convert to CSV file of paired fastq")
parser.add_argument("-f", dest = 'file_path', required=True, default='FN', help="The path to the file")
args = parser.parse_args()
file_in  = args.file_path
all_data ={}
# Convert this now to the list.
file_in=file_in.replace("[",'')
file_in=file_in.replace("]",'')
file_in=file_in.replace("\\",'')
files_in=file_in.split(", ")
files_in.sort()
for file in files_in:

    file_name = (file.split('/')[-1])
    Sample_ID = file_name.split(".")[0]
    if (Sample_ID[-2:]=='_2' or Sample_ID[-2:]=='_1'):
        # this may cause issues if users define the files in this particular way, but hopefully not
        fataq_id= Sample_ID.split("_")[-1]
        Sample_ID=Sample_ID[:-2]
        
    else:
        fataq_id= file_name.split(".")[1]

    try:
        all_data[Sample_ID][fataq_id]=file
    except:
        all_data[Sample_ID]={}
        all_data[Sample_ID][fataq_id]=file

#in case that there is no paired fastq files we check it by:
length_dict = {key: len(value) for key, value in all_data.items()}
listOfKeys = [key  for (key, value) in length_dict.items() if value == 2]
all_data_matched = dict((k, all_data[k]) for k in listOfKeys if k in all_data)

Dataset= pd.DataFrame(all_data_matched).T

# Need to make sure that the data is in correct order, if not it doesnt allign correctly.
Dataset = Dataset[['1', '2']]
Dataset.to_csv('Fastq_files.csv',header=False)