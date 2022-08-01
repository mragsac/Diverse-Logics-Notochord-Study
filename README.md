# Diverse logics and grammar encode notochord enhancers

This GitHub repository contains the code used in the analyses for the following manuscript: 

> **Song, Benjamin P**+, **Michelle F. Ragsac**+, Krissie Tellez, Granton A Jindal, Jessica L Grudzien, Sophia H Le, Emma K Farley. "Diverse logics and grammar encode notochord enhancers." *In submission*, 2022. 
> + *These authors contributed equally*

---

## Background Information

For this study, we searched through the *Ciona intestinalis* genome to find regions containing both *Zic* and *Ets* binding sites that could confer notochord activity in the developing embryo. To do this, we first searched through the genome for candidate regions to include within a massively-parallel reporter assay (MPRA) or "enhancer screen." 

This repository contains analyses for the genome search (`01_genome-search/`) and the MPRA performed in whole *Ciona* embryos (`03_bulk-screen`). We call the linking between barcode tags and associated enhancer the "dictionary" analysis (`02_enhancer-dictionary`). If you would like to perform the analyses for the whole embryo MPRA on your own computer, you can find the raw sequencing data under the following SRA ascession identifier: `PRJNA861319`. 

## Additional Data Resources

Within this study, we also used additional data resources provided by the scientific community in our analyses: 

| Resource Name                                                | Resource Usage | Resource Link | 
| ------------------------------------------------------------ | -------------- | ------------- | 
| *Ciona intestinalis* Genome Sequence (`JoinedScaffold.fa`)   | Genome Search  | [Ghost Database](http://ghost.zool.kyoto-u.ac.jp/datas/) | 
| Ets1 Binding Data for *Mus musculus* (`Ets1/Ets1_8mers.txt`) | Genome Search  | [UniPROBE Database](http://thebrain.bwh.harvard.edu/uniprobe/detailsDef.php?id=414) | 

## Script Parameters

For this study, we used several scripts to analyze our data. We wanted to share the parameters that we used for each of the software tools: 

| Software Tool / Script Name        | User-Defined Parameters        |
| ---------------------------------- | ------------------------------ | 
| `genome-search.py` (Custom Script) | Window Size = 30 bp            | 
| `FLASH` Software Tool              | Maximum Mismatch Density = 0.0 | 

---

## Analysis Procedure

1. Perform the genome search within *Ciona intestinalis* for regions with at least one *Zic* site and two *Ets* sites using `01_genome-search/genome-search.py` with a flanking window size of 30 bp. 
2. Select regions from the search to include in this study.
3. Assemble the barcode tag-enhancer sequence dictionary and analyze the sequencing data: 
	1. Determine the sequencing quality with `fastqc` (`02_enhancer-dictionary/01_fastqc`).
	2. Because we have overlapping sequence between our reads, assemble the entire read sequence using `FLASH` (`02_enhancer-dictionary/02_flash`).
	3. Consolidate unique reads by "collapsing" them and tallying up the number of times a unique sequence appears in our data (`02_enhancer-dictionary/03_collapse`).
	4. Analyze and determine the barcode tag-enhancer sequence dictionary associations (`02_enhancer-dictionary/04_assemble`).