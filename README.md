___
                                            _____                                             
                                         _.'_____`._                                          
                                       .'.-'  12 `-.`.                                        
                                      /,' 11      1 `.\                                       
                                     // 10      /   2 \\                                      
                                    ;;         /       ::                                     
                                    || 9  ----O      3 ||                                     
                                    ::                 ;;                                     
                                     \\ 8           4 //                                      
                                      \`. 7       5 ,'/                                       
                                       '.`-.__6__.-'.'                                        
                                        ((-._____.-))                                         
                                        _))       ((_                                         
                                       '--'       '--'                                        
                                                                                              
    -.       .-.       .-.       .-.       .-.       .-.       .-.       .-.       .-.       .
      \     /   \     /   \     /   \     /   \     /   \     /   \     /   \     /   \     / 
       \   /     \   /     \   /     \   /     \   /     \   /     \   /     \   /     \   /  
        `-~       `-~       `-~       `-~       `-~       `-~       `-~       `-~       `-~   
                                                                                              
                             ██████╗        ██████╗ █████╗ ███╗   ███╗                        
                            ██╔════╝       ██╔════╝██╔══██╗████╗ ████║                        
                            ██║  ██╗ █████╗╚█████╗ ███████║██╔████╔██║                        
                            ██║  ╚██╗╚════╝ ╚═══██╗██╔══██║██║╚██╔╝██║                        
                            ╚██████╔╝      ██████╔╝██║  ██║██║ ╚═╝ ██║                        
                             ╚═════╝       ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝                        
                                                                                              
    -.       .-.       .-.       .-.       .-.       .-.       .-.       .-.       .-.       .
      \     /   \     /   \     /   \     /   \     /   \     /   \     /   \     /   \     / 
       \   /     \   /     \   /     \   /     \   /     \   /     \   /     \   /     \   /  
        `-~       `-~       `-~       `-~       `-~       `-~       `-~       `-~       `-~   
___



# G-SAM RNA pipe-line and analysis code #

Dr. Youri Hoogstrate, 2019 - 2022

---

# Citing recursiveCorPlot

Please cite this paper ([10.1016/j.ccell.2023.02.019](https://doi.org/10.1016/j.ccell.2023.02.019)) when using *recursiveCorPlot* for your
publications:

```
Youri Hoogstrate, Kaspar Draaisma, Santoesha A. Ghisai, Levi van Hijfte, Nastaran Barin, Iris de Heer, Wouter Coppieters,
Thierry P.P. van den Bosch, Anne Bolleboom, Zhenyu Gao, Arnaud J.P.E. Vincent, Latifa Karim, Manon Deckers,
Martin J.B. Taphoorn, Melissa Kerkhof, Astrid Weyerbrock, Marc Sanson, Ann Hoeben, Slávka Lukacova, Giuseppe Lombardi,
 Sieger Leenstra, Monique Hanse, Ruth E.M. Fleischeuer, Colin Watts, Nicos Angelopoulos, Thierry Gorlia, Vassilis Golfinopoulos,
Vincent Bours, Martin J. van den Bent, Pierre A. Robe, Pim J. French,

Transcriptome analysis reveals tumor microenvironment changes in glioblastoma,

Cancer Cell, 2023, ISSN 1535-6108, https://doi.org/10.1016/j.ccell.2023.02.019
```

---

Directory structure:

```
# main R data analyses:

./scripts/load_*.R     [loading of data]
                        - data is typically structured in 3 levels:
                          * per gene
                          * per patient
                          * per resection/sample
./scripts/analysis_*.R [analysis of data, typically exporting cached results to ./cache]
./scripts/vis_*.R      [visualisation of data]


./if/*.ijm             [ImageJ IF analysis scripts]

LICENSE                [copy of MIT license]
```
