# Pluto:Ultra-Fast Haplotype Phasing and Genotype Refinement Tool for Ultra-Low Coverage DNA sequence reads

We propose an efficient statistical method Pluto that enables haplotype phasing and genotype refinement from individual genomes sequenced in ultra-low coverage with orders of magnitude smaller computational cost compared to alternative methods.

Our method integrates the Positional Burros-Wheeler Transform (PBWT) and Variable Length Markov Chains (VLMC) to build a haplotype graph from reference haplotypes to account for genotype uncertainty. We leverage statistical methods, such as Kolmogorov-Smirnov (KS) test, to accurately and efficiently build haplotype graphs for VLMC from conditional distributions rapidly obtained from PBWT.

Our experiments show that our method can achieve comparable or better genotype accuracy to existing tools under various sequencing depth between 0.1x or 4x when using 1000 Genomes as reference haplotypes. At the same time, the computational speed of Pluto is more than 10x faster than Beagle when evaluated in 4x whole genome sequence reads.

## Installation

  - mkdir build
  - cd build
  - cmake ..
  - make
  - make test

## Usage
```
Pluto index

Available Options
   Shotgun Sequences : --refVCF [Empty]
      Optional Files : --includeUnphasedIDs [], --includePhasedIDs [],
                       --excludeUnphasedIDs [], --excludePhasedIDs []
       Graph Builder : --graphComplexity [1400], --PvalueMatrix [],
                       --calPvalueMatrix [], --geneticDistance [],
                       --seed [123456], --onlyHeterSite
        Output Files : --outPrefix [mach1.out]

```
```
phase --refVcf reference.panel.vcf.gz --unphasedVcf target.vcf.gz --outPrefix target.phased
```

## Example Usage
```
Pluto index --refVcf reference.panel.vcf.gz --PvalueMatrix /Users/fanzhang/Downloads/PlutoTest/PvalueMatrix
```
```
phase --refVcf reference.panel.vcf.gz --unphasedVcf target.vcf.gz --outPrefix target.phased
```

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## Bug Report

fanzhang@umich.edu

## License

This project is licensed under the terms of the MIT license.