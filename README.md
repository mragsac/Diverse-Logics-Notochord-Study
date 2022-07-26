# Diverse logics and grammar encode notochord enhancers

This GitHub repository contains the code used in the analyses for the following manuscript: 

> Song, Benjamin P+, Michelle F. Ragsac+, Krissie Tellez, Granton A Jindal, Jessica L Grudzien, Sophia H Le, Emma K Farley. "Diverse logics and grammar encode notochord enhancers." *In submission*, 2022. 
> + *These authors contributed equally*

---

## Background Information

For this study, we searched through the *Ciona intestinalis* genome to find regions containing both *Zic* and *Ets* binding sites that could confer notochord activity in the developing embryo. To do this, we first searched through the genome for candidate regions to include within a massively-parallel reporter assay (MPRA) or "enhancer screen." 

This repository contains analyses for the genome search (`01_genome-search/`) and the MPRA performed in whole *Ciona* embryos (`02_bulk-screen`). If you would like to perform the analyses for the whole embryo MPRA on your own computer, you can find the raw sequencing data under the following SRA ascession identifier: `PRJNA861319`. 

## Additional Resources

Within this study, we also used additional resources provided by the scientific community in our analyses: 

| Resource Name                                            | Resource Usage | Resource Link | 
| -------------------------------------------------------- | -------------- | ------------- | 
| Ciona intestinalis Genome Sequence (`JoinedScaffold.fa`) | Genome Search  | [Ghost Database](http://ghost.zool.kyoto-u.ac.jp/datas/) | 
| Ets Binidng Data (`Ets1/Ets1_8mers.txt`)                 | Genome Search  | [UniPROBE Database](http://thebrain.bwh.harvard.edu/uniprobe/detailsDef.php?id=414) | 

