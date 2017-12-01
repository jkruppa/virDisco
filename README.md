# viralDetectTools
Modules for a viral detection pipeline. Uses external tools and software. No standalone.



## Installation

The development version from github. Do not do it, if you are not sure what you are doing!

```R
devtools::install_github("jkruppa/viralDetectTools")
```

# Tutorial

```R
  A DNAStringSet instance of length 98
       width seq                                             names
  [1]    375 ACCAACACCGGCAAAGCTAGCA...CGAGGAACTGTTCGTGTACGTG AY949829 AY949829...
  [2]    731 GACTTTGCTAGTCTGTACCCCA...GGGTCGACCTGGTGAGCAAGAC AY949831 AY949831...
  [3]    241 GTAACTCGGTGTACGGGTTCAC...TCTACGGCGACACGGACTCTGT AY949832 AY949832...
  [4]    731 GATTTTGCTAGTTTGTACCCCA...GCGTCGACCTGGTGAGCAAGAC AY952776 AY952776...
  [5]    676 CATTCAAGCCCACAATCTGTGC...ACTGTCCACGGATAAAATCTTA AY952777 AY952777...
  ...    ... ...
 [94]   1726 CTAGGGGTAGAATAACAGATAA...TCATCAATATTACAAAAAACTT X75961 X75961 Dol...
 [95]   2212 AGGAGTAAAGTAATTAAGACTC...ATAAATAAATAATTAAAGAAAA Z30086 Z30086 Dol...
 [96]   1453 AGGATTCAAGACTAAGTTGACT...TGATGCAAATTATTAAATAAAA Z30087 Z30087 Dol...
 [97]   1946 AGGGTGCAAGTTGACCAACCAT...CTCTAGACCCCATTAAGAAAAA Z36978 Z36978 Dol...
 [98]   1655 AGGACCAAAGTCCAAGGAATTG...CAAAATGATCAATTATAAAAAA Z47758 Z47758 Dol...
```


```R
   A AAStringSet instance of length 304
       width seq                                             names
   [1]   125 TNTGKASTSFLFNLKYSSDDLL...LLYRPSTTTRRGAMAEELFVYV AAX47052.1_AY9498...
   [2]   244 DFASLYPSIIQAHNLCYSTLIP...RYIGILSTDKILMKGVDLVSKT AAX47054.1_AY9498...
   [3]    80 NSVYGFTGVAQGLPPCLQIAAT...LLASPPAPPYSIHVIYGDTDSV AAX47055.1_AY9498...
   [4]   244 DFASLYPSIIQAHNLCYSTLIP...RYIGILSTDKILMKGVDLVSKT AAX55676.1_AY9527...
   [5]   225 IQAHNLCYSTLIPDGEMHRHPT...FKCLLLLTKKRYIGILSTDKIL AAX55677.2_AY9527...
   ...   ... ...
 [300]   121 MDSNTVSSFQDILMRMSKMQLG...MQALQLLLEVEQEIRTFSFQLI ALZ47758.1_CY1963...
 [301]   506 MAEEQAYHINKGLECLKSLREN...LNDVKSGKDLGEFYQMVKKIIK CAA87685.1_Z47758...
 [302]   506 MAEEQAYHINKGLECLKSLREN...LNDVKSGKDLGEFYQMVKKIIK CAA87687.1_Z47758...
 [303]   160 MLSKLRKPKLSEARPPAKNQAR...ILPLTGDLLPGLRSRDRLTLRL CAA87686.1_Z47758...
 [304]   843 MRVMGMLRNYQQWWTWGILGFW...RIGRAILNVPTRIRQGFERALL ADZ47758.1_HM0707...
```

```R
# A tibble: 284 x 5
       prot_id genebank_id pos_start pos_end    ind
         <chr>       <chr>     <dbl>   <dbl>  <int>
  1 ABC33905.1    DQ288666         1     725 332091
  2 ABC33906.1    DQ288667         1     725 332092
  3 AOV92933.1    KT964777         1     559 373602
  4 AHB63480.1    KF793824       520   20510 498203
  5 AHB63481.1    KF793824     20446   24927 498204
  6 AHB63482.1    KF793824     24933   25220 498205
  7 AHB63483.1    KF793824     25204   26028 498206
  8 AHB63484.1    KF793824     26031   26447 498207
  9 AHB63485.1    KF793824     26440   26958 498208
 10 AHB63486.1    KF793824     26948   27475 498209
 # ... with 274 more rows
```

```R
# A tibble: 98 x 6
    genebank_id tax_id Accession Length                                     Organism
          <chr>  <chr>     <chr>  <int>                                        <chr>
  1    DQ288666 359948  DQ288666    725             Risso's dolphin gammaherpesvirus
  2    DQ288667 319745  DQ288667    725 Atlantic bottlenose dolphin gammaherpesvirus
  3    EF451565  37131  EF451565    353                        Dolphin morbillivirus
  4    EF469546  37131  EF469546    281                        Dolphin morbillivirus
  5    EU039963  37131  EU039963    310                        Dolphin morbillivirus
  6    EU124652  37131  EU124652    188                        Dolphin morbillivirus
  7    EU886967 577566  EU886967   2613               Bottlenose dolphin enterovirus
  8    FJ890355 645425  FJ890355   3990              Bottlenose dolphin astrovirus 1
  9    GQ258353 107323  GQ258353    356               Bottlenose dolphin herpesvirus
 10    GQ258354 107323  GQ258354    356               Bottlenose dolphin herpesvirus
 # ... with 88 more rows, and 1 more variables: Description <chr>
 ```
