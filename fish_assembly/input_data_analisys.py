# -*- coding: utf-8 -*-
from collections import defaultdict, Counter
import itertools

__author__ = "Sergey Aganezov"
__email__ = "aganezov(at)gwu.edu"
__status__ = "develop"

import os
import csv

from bg import BreakpointGraph, GRIMMReader, BGVertex, Multicolor

RAW_DATA_ROOT = os.path.abspath(os.path.sep.join(["..", "data", "fish_genome"]))
GRIMM_DATA_ROOT = os.path.abspath(os.path.sep.join(["..", "data", "grimm"]))

if __name__ == "__main__":
    gtf_files = [file_name for file_name in os.listdir(RAW_DATA_ROOT) if file_name.endswith(".gtf")]
    gene_name_files = [file_name for file_name in os.listdir(RAW_DATA_ROOT) if
                       "gene_name" in file_name and not file_name.startswith(".")]
    print("GTF files:\n\t", "\n\t".join(gtf_files), "\n", sep="")
    print("Gene name files:\n\t", "\n\t".join(gene_name_files), "\n", sep="")

    scaffold_counter = defaultdict(set)
    genes_per_genome = defaultdict(set)
    genes_orth = defaultdict(dict)
    genomes = defaultdict(lambda: defaultdict(list))
    for gtf_files in gtf_files:
        genome_name = gtf_files.split(".")[0]
        full_file_name = os.path.sep.join([RAW_DATA_ROOT, gtf_files])
        with open(full_file_name, "rt") as source:
            reader = csv.reader(source, delimiter='\t', quotechar='"')
            for row in reader:
                scaffold_counter[genome_name].add(row[0])
                genes_per_genome[genome_name].add(row[8].split(" ")[1][1:-2])
                genomes[genome_name][row[0]].append((row[8].split(" ")[1][1:-2], int(row[3]), int(row[4]), row[6]))

    for genome_name in scaffold_counter:
        assert len(genomes[genome_name]) == len(scaffold_counter[genome_name])
    for gene_name_file in gene_name_files:
        genome_name = gene_name_file.split(".")[0]
        full_file_name = full_file_name = os.path.sep.join([RAW_DATA_ROOT, gene_name_file])
        full_file_name = os.path.sep.join([RAW_DATA_ROOT, gene_name_file])
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
                # else:
                # genomes[genome_name][scaffold_name] = sorted(genomes[genome_name][scaffold_name],
                # key=lambda item: (int(item[1]), int(item[2]), item[0]))
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
    with open("not_strand_consistent_gene_ids.txt", "w") as source:
        print("Not strand-consistent gene_ids", file=source)
        print("\n\nNot strand-consistent gene_ids")
        for genome in genomes:
            bad_genes = set()
            for scaffold_name, genes in genomes[genome].items():
                gene, *genes = genes
                current_gene_id, current_strand, current_start, current_finish = gene[0], gene[3], gene[1], gene[2]
                shrunk_genomes[genome][scaffold_name].append(gene)
                for gene in genes:
                    gene_id, start, finish, strand = gene
                    if gene_id != current_gene_id:
                        shrunk_genomes[genome][scaffold_name].append(
                            (current_gene_id, current_start, current_finish, current_strand))
                        current_gene_id = gene_id
                        current_strand = strand
                        current_start = start
                        current_finish = finish
                    else:
                        current_finish = finish
                        if strand != current_strand and gene_id not in bad_genes:
                            bad_genes.add(gene_id)
                            # new_data = tuple(item if cnt != 2 else finish for item, cnt in
                            # enumerate(shrunk_genomes[genome][scaffold_name][-1]))
                            # shrunk_genomes[genome][scaffold_name][-1] = new_data
            print("\tGenome {genome_name} contains {bgid_cnt} strand inconsistent gene ids"
                  "".format(genome_name=genome, bgid_cnt=len(bad_genes)))
            if len(bad_genes) > 0:
                print(genome, file=source)
                print("\n".join(bad_genes), file=source)

    print("\nGlued genomes gene cnt:")
    for genome in shrunk_genomes:
        cnt = sum(map(lambda scaffold: len(shrunk_genomes[genome][scaffold]), shrunk_genomes[genome]))
        print("\t{genome_name}: {g_cnt}".format(genome_name=genome,
                                                g_cnt=cnt))

    scaffold_lengths_per_genome = {}
    scaffold_length_counter_per_genome = {}
    print("\nScaffold lengths")
    for genome in shrunk_genomes:
        scaffold_lengths_per_genome[genome] = [len(shrunk_genomes[genome][scaffold]) for scaffold in
                                               shrunk_genomes[genome]]
        scaffold_length_counter_per_genome[genome] = Counter(scaffold_lengths_per_genome[genome])
        print("\t{genome_name}: {counter}".format(genome_name=genome,
                                                  counter=str(scaffold_length_counter_per_genome[genome])))

    # print(shrunk_genomes.keys())
    for genome in scaffold_counter:
        if genome in shrunk_genomes:
            for scaffold_name in scaffold_counter[genome]:
                shrunk_genomes[genome][scaffold_name] = [(gene_id, (start + finish) / 2, strand) for
                                                         gene_id, start, finish, strand in
                                                         shrunk_genomes[genome][scaffold_name]]
                shrunk_genomes[genome][scaffold_name] = sorted(shrunk_genomes[genome][scaffold_name],
                                                               key=lambda entry: entry[1])

    all_orth = set()
    for genome in shrunk_genomes:
        for key, value in genes_orth[genome].items():
            all_orth.add(value)
    orth_dict = {orth: number for number, orth in enumerate(sorted(all_orth))}

    print("\nOverall different orthologous: {orth_cnt}".format(orth_cnt=len(all_orth)))

    # bad_genomes = []
    bad_genomes = ["Atlantic_cod"]

    for genome_name in scaffold_counter:
        if genome_name in bad_genomes:
            del shrunk_genomes[genome_name]

    grimm_formatted_genomes = defaultdict(lambda: defaultdict(list))
    for genome in shrunk_genomes:
        for scaffold_name in shrunk_genomes[genome]:
            for gene_name, coordinate, strand in shrunk_genomes[genome][scaffold_name]:
                strand = "-" if strand == "-1" else "+"
                orth_gene_name = orth_dict[genes_orth[genome][gene_name]]
                grimm_formatted_genomes[genome][scaffold_name].append((strand, orth_gene_name))

    # print("\nGRIMM Formatted data")
    if not os.path.exists(GRIMM_DATA_ROOT):
        os.mkdir(GRIMM_DATA_ROOT)
    target_directory_rel_name = "8_genomes_no_Atlantic_cod"
    target_directory = os.path.join(GRIMM_DATA_ROOT, target_directory_rel_name)
    if not os.path.exists(target_directory):
        os.mkdir(target_directory)
    for genome in grimm_formatted_genomes:
        with open(os.path.join(target_directory, genome) + ".grimm", "w") as target:
            # print("\t>", genome)
            print(">", genome, file=target, sep="")
            for scaffold_name in grimm_formatted_genomes[genome]:
                # print("\t#", scaffold_name)
                print("#", scaffold_name, file=target)
                # print("\t", " ".join(
                #     strand + str(gene_name) for strand, gene_name in grimm_formatted_genomes[genome][scaffold_name]),
                #       "$")
                print(" ".join(
                    strand + str(gene_name) for strand, gene_name in grimm_formatted_genomes[genome][scaffold_name]),
                      "$", file=target)
        # print()

    bg_graphs = dict()
    for file_name in os.listdir(target_directory):
        file_name = os.path.join(target_directory, file_name)
        with open(file_name, "r") as source:
            bg = GRIMMReader.get_breakpoint_graph(source)
            bg_graphs[file_name.split(".")[0]] = bg
    bg = BreakpointGraph()
    for br_gr in bg_graphs.values():
        bg.update(br_gr, merge_edges=True)

    target_multicolor = Multicolor("Anguilla_japonica")

    print("Breakpoint graph stats:")
    print("\t", "non-infinity nodes count:", len(list(node for node in bg.nodes() if not BGVertex.is_infinity_vertex(node))))
    normal_edges, infinity_edges = [], []
    for edge in bg.edges():
        if edge.is_infinity_edge:
            infinity_edges.append(edge)
        else:
            normal_edges.append(edge)
    print("\t", "non-infinity edges count:", len(normal_edges))
    print("\t", "infinity edges count:", len(infinity_edges))
    print("\t", "connected component count:", len(list(bg.connected_components_subgraphs())))

    print()
    print()

    cc_cnt_with_exactly_two_targeted_infinity_edges = 0
    for cnt, cc in enumerate(bg.connected_components_subgraphs()):
        infinity_edges = defaultdict(int)
        targeted_infinity_edges = []
        targeted_infinity_edges_no_account_for_multiplicity = []
        for edge in cc.edges():
            if edge.is_infinity_edge:
                infinity_edges[tuple(sorted(edge.multicolor.colors))] += 1
                if edge.multicolor == target_multicolor:
                    targeted_infinity_edges.append(edge)
                elif edge.multicolor.colors == target_multicolor.colors:
                    targeted_infinity_edges_no_account_for_multiplicity.append(edge)
        if len(targeted_infinity_edges) > 1 or len(targeted_infinity_edges_no_account_for_multiplicity) > 1:
            print("\t", "processing connected component", cnt)
            print("\t\t", "overall number of infinity edges", sum(infinity_edges.values()))
            print("\t\t", "number of targeted infinity edges accounting for multiplicity", len(targeted_infinity_edges))
            print("\t\t", "number of targeted infinity edges no accounting for multiplicity", len(targeted_infinity_edges_no_account_for_multiplicity))
            if len(targeted_infinity_edges) == 2 or len(targeted_infinity_edges_no_account_for_multiplicity) == 2:
                cc_cnt_with_exactly_two_targeted_infinity_edges += 1
            print()

    print("possible assembly points: ", cc_cnt_with_exactly_two_targeted_infinity_edges)
