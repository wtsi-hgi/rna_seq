
"""This code snippet for nexflow performs the gz compression checks and if not compressed then performs this compression"""
import os
import argparse

def fastq_copy_or_compress(file_in,file_out):
    print('Assessing file extension and performing either copy or gz compression')
    if file_in.split('.')[-1]=='gz':
        file_out=file_out.replace('.gz.gz','.gz')
        os.system(f"cp {file_in} {file_out}")
    else:
        os.system(f"gzip -k --force {file_in} {file_out}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarise samtools idxstats output even further")
    parser.add_argument("-f", dest = 'file_path', required=True, default='FN', help="The path to the file")
    parser.add_argument("-o", dest = 'output_path', required=True, default='FN', help="The path to the file output")
    args = parser.parse_args()
    file_in  = args.file_path
    file_out  = args.output_path
    fastq_copy_or_compress(file_in,file_out)