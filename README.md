This directory contains python scripts to take Target Finder output files and create corresponding Weka input files for a Naive Bayes machine learning classifier for distinguishing true and false positive Target Finder predictions (Cernadas et al, 2014; Wilkins et al, under review).

#License

All source code is available under an ISC license.

Copyright (c) 2015-2016, Katie Wilkins kxw116@gmail.com and Nick Booher njbooher@gmail.com.

Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted, provided that the above copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

#Description

The script `msu_gff_parser.py` creates a dictionary from a gff genome annotation file. Like all scripts in this directory, it was designed for the MSU7 rice genome annotation (Kawahara et al, 2013) and may not handle the quirks of other gff file formats.
 
The script `generate_msu7_promoters.py` uses a fasta-formatted genome sequence and the corresponding annotation to create a fasta-formatted file containing all promoters in that genome, where a promoter is defined as the 1000 base pairs upstream of the annotated transcriptional start site plus the 5' UTR if it is annotated. The script also creates pickled files containing all the genome features used by the machine learning classifier. The script depends on `msu_gff_parser.py`. The output file prefix should end with the beginning of a file name, not with a directory. The script can be run as follows:

>`python ./generate_msu7_promoters.py -a annotation_file.gff -g genome_file.fa -o output_file_prefix`

The script `get_msu7_Ypatch_TATABox.py` searches a fasta-formatted file of promoter sequences for TATA box and Y patch sequences as defined by Yamamoto et al (2007). The script can be run as follows:

>`python ./get_msu7_Ypatch_TATABox.py -p promoters.fa -o output_file_prefix`

The script `write_weka_commands_and_input_files` takes as input the text tab-delimited output from the Target Finder TAL effector binding site prediction tool, as well as the information about genomic features output by the previous two scripts, and outputs Weka input files for the Target Finder results and Weka commands for running the classifier on these input files. The script can be run as follows:

>`python ./write_weka_commands_and_input_files -i input_folder_containing_target_finder_output -g output_file_prefix_from_previous_scritps -o output_folder -s number_of_lines_in_target_finder_header -c AllFeaturesPlusIdNB2.model -r file_for_weka_output -w weka.jar`

The file AllFeaturesPlusIdNB2.model is a Weka model file that contains a Naive Bayes classifier for distinguishing true and false positive Target Finder predictions. It is the same classifier created by Cernadas et al (2014) updated to use transcriptional and translational start sites from an annotation file in Wilkins et al (under review).

#References

Cernadas, R.A., Doyle, E.L., Niño-Liu, D.O., Wilkins, K.E., Bancroft, T., Wang, L., Schmidt, C.L., Caldo, R., Yang, B., White, F.F., Nettleton, D., Wise, R.P., and Bogdanove, A.J. (2014). Code-assisted discovery of TAL effector targets in bacterial leaf streak of rice reveals contrast with bacterial blight and a novel susceptibility gene. PLoS Pathogens 10, 1-24. doi: 10.1371/journal.ppat.1003972.

Hall M, Frank E, Holmes G, Pfahringer B, Reutemann P, Witten I. (2009) The WEKA data mining software: an update. SIGKDD Explor Newsl 11, 10–18. doi: 10.1145/1656274.1656278

Kawahara, Y., De La Bastide, M., Hamilton, J., Kanamori, H., Mccombie, W.R., Ouyang, S., Schwartz, D., Tanaka, T., Wu, J., Zhou, S., Childs, K., Davidson, R., Lin, H., Quesada-Ocampo, L., Vaillancourt, B., Sakai, H., Lee, S.S., Kim, J., Numa, H., Itoh, T., Buell, C.R., and Matsumoto, T. (2013). Improvement of the Oryza sativa Nipponbare reference genome using next generation sequence and optical map data. Rice 6, 4. doi: 10.1186/1939-8433-6-4.

Wilkins K, Booher N, Wang L, and Bogdanove A. (under review) TAL effector content and host transcriptional response across diverse strains of the rice bacterial leaf streak pathogen *Xanthomonas oryzae* pv. oryzicola.

Yamamoto YY, Ichida H, Matsui M, Obokata J, Sakurai T, Satou M, Seki M, Shinozaki K, and Abe T. (2007) Identification of plant promoter constituents by analysis of local distribution of short sequences. BMC Genomics 8, 67. doi: 10.1186/1471-2164-8-67

