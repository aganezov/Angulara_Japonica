# -*- coding: utf-8 -*-
from collections import defaultdict
import itertools

__author__ = "Sergey Aganezov"
__email__ = "aganezov(at)gwu.edu"
__status__ = "develop"

import os
import csv

DATA_ROOT = os.path.abspath(os.path.sep.join(["..", "data", "fish_genome"]))

if __name__ == "__main__":
    gtf_files = [file_name for file_name in os.listdir(DATA_ROOT) if file_name.endswith(".gtf")]
    gene_name_files = [file_name for file_name in os.listdir(DATA_ROOT) if
                       "gene_name" in file_name and not file_name.startswith(".")]
    print("GTF files:\n\t", "\n\t".join(gtf_files), "\n", sep="")
    print("Gene name files:\n\t", "\n\t".join(gene_name_files), "\n", sep="")

    scaffold_counter = defaultdict(set)
    genes_per_genome = defaultdict(set)
    genes_orth = defaultdict(dict)
    genomes = defaultdict(lambda: defaultdict(list))
    for gtf_files in gtf_files:
        genome_name = gtf_files.split(".")[0]
        full_file_name = os.path.sep.join([DATA_ROOT, gtf_files])
        with open(full_file_name, "rt") as source:
            reader = csv.reader(source, delimiter='\t', quotechar='"')
            for row in reader:
                scaffold_counter[genome_name].add(row[0])
                genes_per_genome[genome_name].add(row[8].split(" ")[1][1:-2])
                genomes[genome_name][row[0]].append((row[8].split(" ")[1][1:-2], row[3], row[4], row[6]))

    for genome_name in scaffold_counter:
        assert len(genomes[genome_name]) == len(scaffold_counter[genome_name])
    for gene_name_file in gene_name_files:
        genome_name = gene_name_file.split(".")[0]
        full_file_name = full_file_name = os.path.sep.join([DATA_ROOT, gene_name_file])
        full_file_name = os.path.sep.join([DATA_ROOT, gene_name_file])
        with open(full_file_name, "rt") as source:
            reader = csv.reader(source, delimiter='\t', quotechar='"')
            if genome_name != "Anguilla_japonica":
                next(reader, None)
            for row in reader:
                if len(row) == 2:
                    gene_name, annotation = row
                else:
                    gene_name, annotation = row[0], row[2]
                genes_orth[genome_name][gene_name] = annotation
    print("Scaffold counter:\n\t",
          "\n\t".join(genome_name + ": " + str(len(scaffolds)) for genome_name, scaffolds in scaffold_counter.items()),
          "\n", sep="")
    print("Gene counter:\n\t",
          "\n\t".join(genome_name + ": " + str(len(genes)) for genome_name, genes in genes_per_genome.items()),
          "\n\n", sep="")

    annotated = defaultdict(set)
    print("Gene annotations:")
    for genome_name in scaffold_counter:
        annotated[genome_name] = set(
            gene_name for gene_name in genes_per_genome[genome_name] if gene_name in genes_orth[genome_name])
        print("\t", genome_name + ":", len(annotated[genome_name]), " out of ", len(genes_per_genome[genome_name]))

    print("\n\nShared genes (as annotations):")
    for genome_1, genome_2 in itertools.combinations(scaffold_counter.keys(), 2):
        print("\t", genome_1, "&", genome_2, ":",
              len(set(genes_orth[genome_1].values()).intersection(set(genes_orth[genome_2].values()))))

    print("\n\nNon empty genome scaffolds count:")
    for genome_name in scaffold_counter:
        for scaffold_name in scaffold_counter[genome_name]:
            genomes[genome_name][scaffold_name] = [(gene_name, start, end, strand) for gene_name, start, end, strand in
                                                   genomes[genome_name][scaffold_name] if
                                                   gene_name in annotated[genome_name]]

            if len(genomes[genome_name][scaffold_name]) == 0:
                del genomes[genome_name][scaffold_name]
            else:
                genomes[genome_name][scaffold_name] = sorted(genomes[genome_name][scaffold_name],
                                                             key=lambda item: (int(item[1]), int(item[2]), item[0]))
        print("\t", genome_name, ":", len(genomes[genome_name]), "vs", len(scaffold_counter[genome_name]))

    # print("\n\n Gene parts check:")
    # for genome_name in scaffold_counter:
    #     for scaffold_name in scaffold_counter[genome_name]:
    #         if len(genomes[genome_name][scaffold_name]) > 1:
    #             for prev_gene, following_gene in zip(genomes[genome_name][scaffold_name][:-1],
    #                                                  genomes[genome_name][scaffold_name][1:]):
    #                 if int(prev_gene[2]) > int(following_gene[1]):
    #                     print("\t", genome_name, ":", scaffold_name, ": ", prev_gene, " -- ", following_gene)
