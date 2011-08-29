#Testing various filters
--bed "../testing/data/refFlat.txt.gz.1" -v "../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1" -t "../testing/data/458_traits.csv" -g "T2D" --gene SAMD11 -p results/basicgraph --nograph --png
--bed "../testing/data/refFlat.txt.gz.1" -v "../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1" -t "../testing/data/458_traits.csv" -g "T2D" --gene SAMD11 -p results/firstgraph --title "Only markers that passed all filters (default)" --nograph --png
--bed "../testing/data/refFlat.txt.gz.1" -v "../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1" -t "../testing/data/458_traits.csv" -g "T2D" --gene SAMD11 -p results/secondgraph --filter SBFilter --title "SBFilter enabled" --nograph --png
--bed "../testing/data/refFlat.txt.gz.1" -v "../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1" -t "../testing/data/458_traits.csv" -g "T2D" --gene SAMD11 -p results/thirdgraph --filter SNPQDFilter --title "SNPQDFitler enabled" --nograph --png
--bed "../testing/data/refFlat.txt.gz.1" -v "../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1" -t "../testing/data/458_traits.csv" -g "T2D" --gene SAMD11 -p results/fourthgraph --filter "SBFilter, SNPQDFilter" --title "SBFilter and SNPQDFilter enabled" --nograph --png

#Testing different gene
--bed "../testing/data/refFlat.txt.gz.1" -v "../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1" -t "../testing/data/458_traits.csv" -g "T2D" --gene SCRIB -p results/newgene --title "New gene, default filter" --nograph --png
#Testing user-provided colors, dimensions, introns
--bed "../testing/data/refFlat.txt.gz.1" -v "../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1" -t "../testing/data/458_traits.csv" -g "T2D" --gene SAMD11 -p results/diffcolors --title "Different colors, dimensions, introns" --nograph --png --color "#00ff00,#ff0000,#0000ff" --dimensions "12,8" --nointrons
#Testing codons, exoncolor, nolegend
--bed "../testing/data/refFlat.txt.gz.1" -v "../testing/data/458_samples_from_bcm_bi_and_washu.annot.vcf.gz.1" -t "../testing/data/458_traits.csv" -g "T2D" --gene SAMD11 -p results/testing --title "codons, exoncolor, nolegend" --nograph --png --exoncolor "#00ff00,#ff0000" --nolegend --codons
