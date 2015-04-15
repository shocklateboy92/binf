package iHMMuneAlign;


public class CodonHelperFunctions
{

    static final char A_NUCLEOTIDE = 97;
    static final char C_NUCLEOTIDE = 99;
    static final char G_NUCLEOTIDE = 103;
    static final char T_NUCLEOTIDE = 116;
    static final char N_NUCLEOTIDE = 110;
    static final int CODON_LENGTH = 3;

    public CodonHelperFunctions()
    {
    }

    public static boolean isAminoAcidW(char codon[])
    {
        if(!isT(codon[0]))
        {
            return false;
        }
        if(!isG(codon[1]))
        {
            return false;
        }
        return isG(codon[2]);
    }

    public static boolean isAminoAcidG(char codon[])
    {
        if(!isG(codon[0]))
        {
            return false;
        }
        if(!isG(codon[1]))
        {
            return false;
        }
        return isX(codon[2]);
    }

    public static boolean isAminoAcidX(char codon[])
    {
        if(!isX(codon[0]))
        {
            return false;
        }
        if(!isX(codon[1]))
        {
            return false;
        }
        return isX(codon[2]);
    }

    public static boolean isTAG(char codon[])
    {
        if(!isT(codon[0]))
        {
            return false;
        }
        if(!isA(codon[1]))
        {
            return false;
        }
        return isG(codon[2]);
    }

    public static boolean isTAA(char codon[])
    {
        if(!isT(codon[0]))
        {
            return false;
        }
        if(!isA(codon[1]))
        {
            return false;
        }
        return isA(codon[2]);
    }

    public static boolean isTGA(char codon[])
    {
        if(!isT(codon[0]))
        {
            return false;
        }
        if(!isG(codon[1]))
        {
            return false;
        }
        return isA(codon[2]);
    }

    private static boolean isN(char nucleotide)
    {
        return nucleotide == 'n';
    }

    private static boolean isX(char nucleotide)
    {
        return isA(nucleotide) || isC(nucleotide) || isG(nucleotide) || isT(nucleotide);
    }

    private static boolean isA(char nucleotide)
    {
        return Character.toLowerCase(nucleotide) == 'a';
    }

    private static boolean isC(char nucleotide)
    {
        return Character.toLowerCase(nucleotide) == 'c';
    }

    private static boolean isG(char nucleotide)
    {
        return Character.toLowerCase(nucleotide) == 'g';
    }

    private static boolean isT(char nucleotide)
    {
        return Character.toLowerCase(nucleotide) == 't';
    }

    private static boolean isR(char nucleotide)
    {
        return isA(nucleotide) || isG(nucleotide);
    }

    private static boolean isY(char nucleotide)
    {
        return isC(nucleotide) || isT(nucleotide);
    }

    private static boolean isW(char nucleotide)
    {
        return isA(nucleotide) || isT(nucleotide);
    }
}
