#Generate mutated RNA based on the fasta, bed and vcf file, select only the validated
RNA.mutate <- RNA.validate(fasta = fasta, vcf = vcf, bed = bed)
RNA.mutate <- subset(RNA.mutate, MATCH == TRUE)

#Perform secondary structure prediction
ori.structure <- predict.Structure(name = RNA.mutate$NAME, seq = RNA.mutate$SEQ)
