def translate_that(seq):
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',  #Data table
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    }
    protein = ""
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon]   #return proten which contents are not of a length divisible by 3
    return protein

def read_seq(seq):
    with open(seq, "r") as f: 
        seq = f.read() 
    seq = seq.replace("\n", "")  #read the contents in the DNA textfile change the formatted tect
    seq = seq.replace("\r", "")
    return seq

def mutate():
    seq =read_seq ("DNA.txt")
    seqNormal = seq.replace("a", "A") # chanes "a" to "A"
    seqMutated = seq.replace ("a", "T") #replace "a" with "T"
    
    with open ("mutatedDNA.txt", "w") as f:
        f.write(seqMutated)                                                                                                            
                                               # write new results with changes
    with open ("normalDNA.txt", "w") as f:
        f.write(seqNormal)

print(translate_that("AACT"))

def txtTranslate():
    mutate()
    seqMutated = read_seq("mutatedDNA.txt")
    seqNormal = read_seq("normalDNA.txt")         # read the output sequence of the mutatedDNA and normalDNA
        
    print(translate_that(seqMutated))
    print(translate_that(seqNormal))             #display results
    
txtTranslate()   #calls the translate function so to print out the output
  
    
