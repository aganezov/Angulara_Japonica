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

    with open("not_scaffold_consistent_gene_ids.txt", "w") as source:
        print("Not scaffold-consistent genes:", file=source)
        print("\n\nNon-consistent genes:")
        for genome in genomes:
            visited_genes = defaultdict(set)
            bad_genes = set()
            for scaffold_name, genes in genomes[genome].items():
                visited_genes_on_this_scaffold = set()
                current_gene, *genes = genes
                visited_genes_on_this_scaffold.add(current_gene[0])
                for gene in genes:
                    gene_name, start, finish, strand = gene
                    if gene_name in visited_genes and gene_name not in bad_genes:
                        bad_genes.add(gene_name)
                    visited_genes_on_this_scaffold.add(gene_name)
                for gene_name in visited_genes_on_this_scaffold:
                    visited_genes[gene_name].add(scaffold_name)
            print("\tGenome {genome_name} contains {bgid_cnt} scaffold inconsistent gene ids".format(genome_name=genome,
                                                                                                     bgid_cnt=len(
                                                                                                         bad_genes)))
            if len(bad_genes) > 0:
                print(genome, file=source)
            for gene_name in bad_genes:
                print("\t\tGene id {gene_id} was found on multiple scaffolds ({scaffold_list})"
                      "".format(gene_id=gene_name,
                                scaffold_list=",".join(visited_genes[gene_name])), file=source)

    shrunk_genomes = defaultdict(lambda: defaultdict(list))
    with open("non_strand_consistent_gene_ids.txt", "w") as source:
        print("Not strand-consistent gene_ids", file=source)
        print("\n\nNot strand-consistent gene_ids")
        for genome in genomes:
            bad_genes = set()
            for scaffold_name, genes in genomes[genome].items():
                current_gene_id, *genes = genes
                current_gene_id, current_strand = current_gene_id[0], current_gene_id[3]
                shrunk_genomes[genome][scaffold_name].append(current_gene_id)
                for gene in genes:
                    gene_id, start, finish, strand = gene
                    if gene_id != current_gene_id:
                        current_gene_id = gene_id
                        current_strand = strand
                        shrunk_genomes[genome][scaffold_name].append(gene)
                    elif strand != current_strand and gene_id not in bad_genes:
                        bad_genes.add(gene_id)
            print("\tGenome {genome_name} contains {bgid_cnt} strand inconsistent gene ids"
                  "".format(genome_name=genome, bgid_cnt=len(bad_genes)))
            if len(bad_genes) > 0:
                print(genome, file=source)
                print("\n".join(bad_genes), file=source)
