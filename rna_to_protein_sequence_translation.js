/* Codewars Kata - http://www.codewars.com/kata/rna-to-protein-sequence-translation/javascript
/  Given a string of RNA, create a funciton which translates the RNA into its protein sequence.
/  Also, added a few extra functions to find molecular weight and isoelectric point of translated protein sequence.
*/

var codons = {
    // Phenylalanine
    'UUC':'F', 'UUU':'F',
    // Leucine
    'UUA':'L', 'UUG':'L', 'CUU':'L', 'CUC':'L','CUA':'L','CUG':'L', 
    // Isoleucine
    'AUU':'I', 'AUC':'I', 'AUA':'I', 
    // Methionine
    'AUG':'M', 
    // Valine
    'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V', 
    // Serine
    'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S', 'AGU':'S', 'AGC':'S', 
    // Proline
    'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 
    // Threonine
    'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    // Alanine
    'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 
    // Tyrosine
    'UAU':'Y', 'UAC':'Y', 
    // Histidine
    'CAU':'H', 'CAC':'H',
    // Glutamine
    'CAA':'Q', 'CAG':'Q', 
    // Asparagine
    'AAU':'N', 'AAC':'N', 
    // Lysine
    'AAA':'K', 'AAG':'K',
    // Aspartic Acid
    'GAU':'D', 'GAC':'D', 
    // Glutamic Acid
    'GAA':'E', 'GAG':'E',
    // Cystine
    'UGU':'C', 'UGC':'C',
    // Tryptophan
    'UGG':'W', 
    // Arginine
    'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R', 
    // Glycine
    'GGU':'G',  'GGC':'G', 'GGA':'G', 'GGG':'G', 
    // Stop codon
    'UAA':'Stop', 'UGA':'Stop', 'UAG':'Stop' 
};

var pI = {
    'G':5.97,  'A':6.02, 'V':5.97,
    'L':5.98,  'I':6.02, 'F':5.48,
    'W':5.88,  'Y':5.65, 'H':7.58,
    'S':5.68,  'T':6.53, 'M':5.75,
    'C':5.14,  'D':2.87, 'E':3.22,
    'N':5.41,  'Q':3.22, 'K':9.74,
    'R':10.76, 'P':6.10,
};

var weight = {
    'A':89.1,  'R':174.2, 'N':132.1,
    'D':133.1, 'C':121.2, 'E':147.1,
    'Q':146.2, 'G':75.1,  'H':155.2,
    'I':131.2, 'L':131.2, 'K':146.2, 
    'M':149.2, 'F':165.2, 'P':115.1,
    'S':105.1, 'T':119.1, 'W':204.2,
    'Y':181.2, 'V':117.1,
}

function translation(rna, callback){
    var protein = "", sequence = "";
    for(var i = 0; i < rna.length; i+=3){
        if (rna[i+3] !== undefined){
            sequence = rna.slice(i, i+3);
            if (codons[sequence] == "Stop") break;
            protein += codons[sequence];
        }
    }
    return typeof callback == "function" ? callback(protein) : protein;
}

//sum of the atomic weights of the atoms in a protein (g/mol)
function molecularWeight(protein){
    var sum = protein.split("").reduce(function(total,amino){
        return total + weight[amino];
    }, 0);
    return sum;
}

// Isoelectric point: 
// pH at which the amino acid in protein do not migrate in an electric field.
// Net charge of protein is zero.
function meanPI(protein, callback){
    var amino_acids = protein.split("")
    var sum = amino_acids.reduce(function(charge,amino){
        return charge + pI[amino];
    }, 0);
    var avg_pI = sum/amino_acids.length;
    return typeof callback == "function" ? callback(avg_pI) : avg_pI;
}

//Determine if isoelectric point makes the protein basic, acidic, or neutral
function acidOrBase(pI){
    var type = "acidic";
    if (pI > 7.00){
        type = "basic";
    }
    else if (pI == 7.00){
        type = "neutral";
    }
    return type;
}
var mrna = "GGUAUCGUGCAAUGUUGCACUUCCAUU"; //insulin mRNA strand from human
var insulin = translation(mrna);
var i_weight = translation(mrna, molecularWeight);
var i_pI = translation(mrna, function(protein){
    return meanPI(protein, acidOrBase);
});
