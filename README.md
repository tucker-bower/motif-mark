# Motif-Marker
#### by Tucker Bower

## Program Summary:
This python script locates nucleotide sequence patterns (motifs) of interest in DNA or RNA sequences from a .fasta file and outputs a vector image (.svg) figure that visualizes the location of each instance of each motif along each sequence, to scale.

Note that 

## Software Dependencies
The details of the bio-conda environment this program was developed and tested within can be found in this [.yml file](https://github.com/tucker-bower/motif-mark/blob/main/environment.yml)

## Executing the Script
To get started, here is a command that runs the program with the example data from this repository:
```
python motif_mark.py -f example_data/Figure_1.fasta -m example_data/Fig_1_motifs.txt
```

### Arguments 

Parameter Name | Explanation
------------ | -------------
--fasta (-f) |  Input .fasta file with sequence introns as lowercase and exons as uppercase
--motifs (-m) | File containing motif sequences of interest (maximum: 5)


## Example Output
![figure1.svg](https://github.com/tucker-bower/motif-mark/blob/main/example_data/Figure_1.svg)
