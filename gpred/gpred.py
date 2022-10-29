import argparse
import sys
import os
import csv
import re

# Arnaud LY

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int,
                        default=50, help="Minimum gene length to consider")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int,
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes (shine box not included).")
    parser.add_argument('-p', dest='predicted_genes_file', type=str,
                        default=os.curdir + os.sep + "predict_genes.csv",
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=str,
                        default=os.curdir + os.sep + "genes.fna",
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file):
    """Extract the complete genome sequence as a single string
    """
    seq = ""
    with open(fasta_file, "r") as file:
        for line in file:
            if line[0] == ">":
                continue
            seq += line.strip()
    if seq.isupper():
        return seq
    else:
        print("Sequence containing lowercase characters")


def find_start(start_regex, sequence, start, stop):
    """Find the start codon
    """
    match = start_regex.search(sequence, start, stop)

    if match is not None:
        return match.start()
    else:
        print('No start')


def find_stop(stop_regex, sequence, start):  # Attention a return match.end !!!
    """Find the stop codon
    """
    for match in stop_regex.finditer(sequence, start):
        if match is not None:
            if (match.start(0) - start) % 3 == 0:
                return match.start(0)


def has_shine_dalgarno(shine_regex, sequence, start, max_shine_dalgarno_distance):
    """Find a shine dalgarno motif before the start codon
    """

    match = shine_regex.finditer(sequence, max(
        0, start - max_shine_dalgarno_distance), start)
    for i in match:
        if i.end(0) < start - 6:
            return True
        return False
    return False


def predict_genes(sequence, start_regex, stop_regex, shine_regex,
                  min_gene_len, max_shine_dalgarno_distance, min_gap):
    """Predict most probable genes
    """

    start = 0
    gene_list = []

    while (len(sequence) - start) >= min_gap:
        start = find_start(start_regex, sequence, start, len(sequence))
        if start is not None:
            stop = find_stop(stop_regex, sequence, start)
            if stop is not None:
                l = len(sequence[start:stop])
                if l >= min_gene_len:
                    if has_shine_dalgarno(shine_regex, sequence, start, max_shine_dalgarno_distance):
                        gene_list.append([start + 1, stop + 3])
                        start = min_gap + stop + 2  # Stop = position 1ère lettre du codon
                    start += 1
                start += 1
            start += 1
    return gene_list


def write_genes_pos(predicted_genes_file, probable_genes):
    """Write list of gene positions
    """
    try:
        with open(predicted_genes_file, "wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_genes(fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp):
    """Write gene sequence in fasta format
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i, gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i+1, os.linesep,
                    fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j, gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(kmer):
    """Get the reverse complement"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in kmer[::-1]])


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Initialize necessary data
    args = get_arguments()
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    seq = read_fasta(args.genome_file)
    seq = seq.replace("U", "T")
    
    probable_genes = predict_genes(seq, start_regex, stop_regex, shine_regex,
                                   args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap)

    # We reverse and complement
    sequence_rc = reverse_complement(seq)
    probable_genes_rc = predict_genes(sequence_rc, start_regex, stop_regex, shine_regex,
                                    args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap)
    probable_genes_comp = []
    for i in range(len(probable_genes_rc)):
        p1 = len(seq)-probable_genes_rc[i][1]
        p2 = len(seq)-probable_genes_rc[i][0]
        probable_genes_comp.append([p1, p2])

    # Call to output functions
    write_genes_pos(args.predicted_genes_file, probable_genes)
    write_genes(args.fasta_file, seq, probable_genes,
                sequence_rc, probable_genes_comp)


if __name__ == '__main__':
    main()
