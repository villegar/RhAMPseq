
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RhAMPseq <img src="https://raw.githubusercontent.com/villegar/RhAMPseq/master/figures/logo.png" alt="logo" align="right" height=200px/>

Analysis of rhAmpSeq data.

## Data description

**Author:**
    J.L.

### Files

  - <span style="color: #03045e; background-color: #E5F1FF; padding: 2px 5px;">hap\_genotype</span>:
    haplotype genotyping data for all samples. This is in a number code.
    The number code corresponds to a specific sequenced haplotype at
    each
    locus.
  - <span style="color: #03045e; background-color: #E5F1FF; padding: 2px 5px;">HaplotypeAllele.fasta</span>:
    each marker is a conserved region of the genome, so the haplotype
    allele may differ only by a single SNP, or it could be many SNPS or
    indels.  
  - <span style="color: #03045e; background-color: #E5F1FF; padding: 2px 5px;">readCountMatrixFile</span>:
    number of read for each amplicon and each sample. This data is also
    provided in the
    <span style="color: #03045e; background-color: #E5F1FF; padding: 2px 5px;">hap\_genotype</span>
    file.

### Example

| Locus              | Haplotypes                                                                        | blank\_\_vDNAfen551A01\_A01 | X160\_271\_303\_\_vDNAfen551A01\_A02 | X160\_271\_311\_\_vDNAfen551A01\_A03 |
| :----------------- | :-------------------------------------------------------------------------------- | :-------------------------- | :----------------------------------- | :----------------------------------- |
| rh\_chr1\_10143990 | 2(0.37);4(0.16);9(0.10);5(0.09);1(0.07);3(0.07);8(0.05);13(0.04);6(0.03);7(0.03); | ./.:0                       | 4/2:34,24                            | 2/2:71                               |

where:

  - `Locus` is the marker.

  - `Haplotypes` lists the frequency of haplotypes at this `Locus` in
    the *ENTIRE* dataset, not just our samples. You can mostly ignore
    this.

  - Next column starts with `blank`, this is a blank control well in the
    plate. Here at this `Locus` it sequenced as `./.`, which means
    `NULL`. and `:0` which is the read depth.

  - Next column is sample `303` from our mapping family, its on the
    `fen551` plate. Here you see we detect `4/2` and `34, 24` for
    haplotype 4 with 34 reads, and 2 with 24 reads. This locus is
    *heterozygous*.

  - Next column is sample `311` from our mapping family. Here you see we
    only detect haplotype `2` at `71` reads. This is *homozygous*.

So, if we wanted to know the DNA sequence of this marker, we would look
up haplotypes 4 and 2 from the
<span style="color: #03045e; background-color: #E5F1FF; padding: 2px 5px;">HaplotypeAllele.fasta</span>
file.

One other thing to note that is very important. the marker name,
`rh_ch1_10143990` stands for rhampseq, chromosome 1, position 10143990.
This position information does not indicate the position of the
polymorphism detected at this marker, it is the position of the primer
sequence. So when we map these markers, they will be very near, but not
exactly at, these positions.

## Initial analysis
