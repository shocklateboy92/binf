package iHMMuneAlign;


class Codon
{

    char codon[];
    String sequenceString;
    int nuclSeqPos;

    public Codon(char codon[], String sequenceString, int nuclSeqPos)
    {
        this.codon = codon;
        this.sequenceString = sequenceString;
        this.nuclSeqPos = nuclSeqPos;
    }
}
