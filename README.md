wigreader
=========

Information on using WigReader in your own scripts can be found here: http://ron.gejman.com/blog/2012/07/27/wigreader-a-chip-seq-tool/

Some C programs for handling wiggle files.

The ones I have spent the most time on, in terms of using and testing, are

**WigReader:** Base class for reading wiggle files.

**normalizeWig:** Used to normalize wiggle files to RPKM.

**geneProfile:** Used to output wiggle values along a length of genome (could be a gene or a promoter, or anything you can put into a bed file)

The other tools I have used at various times for my projects, but they have not been tested as rigorously. For example, edge-cases that I haven't tested for may fail. 


These tools are really only going to be useful for bioinformaticians trying to extract data from wiggle files. 

Creating wiggle files
=========
If you have a BAM file, you can create a wiggle file using IGVtools (e.g. `igvtools count -e 100 -w 25 chip.bam "chip.cov.tdf, chip.cov.wig" genome.genome`). However, these wiggle files are not normalized. You have to normalize the output using `normalizeWig` to get wiggle files in RPKM.


Installation
=========
`make`

That's it!

Notes
=========
These scripts are very memory intensive. They are designed to be run on computers with lots of ram (read: > 32 GBs) and to work on very large ChIP-Seq files.

There may be bugs. I am not responsible for any bugs.

The file naming conventions (e.g. output files) are a bit strange and "dictatorial" -- they are setup to spit out files named the way I like as data progresses from one script to the next.

The most important class here is wigreader. You can use wigreader in any of your own C/C++ programs

License
=========
You may use any of these scripts for commericial or non-commercial purposes, as you see fit. If you do use these tools and they result in a figure or data for a publication, please cite me or the tools in some way.

Contact
=========
Please contact me if you find any bugs. If you fix a bug, please push the fix to me and I will accept it!