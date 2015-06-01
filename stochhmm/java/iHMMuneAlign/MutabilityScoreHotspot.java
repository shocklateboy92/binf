package iHMMuneAlign;

import java.io.PrintStream;

// Referenced classes of package iHMMuneAlign:
//            MutabilityScore

/**
 * this class holds methods for calculating the mutability score of a nuclteotide, based on
 * on whether it is part of a mutational hot spot or not.
 */

public class MutabilityScoreHotspot
    implements MutabilityScore
{

    final char A_NUCLEOTIDE = 'a';
    final char C_NUCLEOTIDE = 'c';
    final char G_NUCLEOTIDE = 'g';
    final char T_NUCLEOTIDE = 't';
    final char N_NUCLEOTIDE = 'n';
    
    // the probability of a mutation at a nt. at pos "p" is already calculated,
    // we only find the mutability score (greater or less than 1)

    // WAN covers 1/8 of all trinucleotides (if assuming a random distribution)
    // but covers 1/3 of all mutations, therefore, the mutatbility score for
    // "a" nts. who are part of a WAN trinucleotide is:
    final double ONE_THIRD = ((double)1/(double)3);
    final double ONE_SIXTH = ((double)1/(double)6);
    final double ONE_EIGTH = ((double)1/(double)8);
    final double ONE_QUARTER = ((double)1/(double)4);
    final double THREE_QUARTER = ((double)3/(double)4);

    // RGYW covers (1/2 * 1/4 * 1/2 * 1/2) 1/32 of all quadnucleotides, but
    // covers 1/6 of all mutations, therefore, the mutability score for a "g" nt.
    // which is part of a RGYW quadnucleotide is:
    final double ONE_OVER_32 = ((double)1/(double)32);
 
    // the mutability score for a non-hotspot area is the fraction of mutations
    // not part of any hotspot: 12/12 - 4/12(WAN) - 2/12(RGYW) - 2/12(WRCY) = 4/12 = 1/3
    // divided by the fraction of nts. not in a hot spot
    // 32/32 (all nucleotides) - 1/8 (WAN) - 1/32 (RGYW) - 1/32 (WRCY) = 26/32
    final double _26_OVER_32 = ((double)26/(double)32);



    //coverage indicates the frequency at which the hotspots
    //occur within the germline, regardless of whether they are mutated or not
    //originals just based on occurance within all possible 3-(WAN) or 4(RGYW,WRCY)-mers
    //final double RGYW_Nucleotide_coverage = ONE_OVER_32;
    //final double WRCY_Nucleotide_coverage = ONE_OVER_32;
    //final double WAN_Nucleotide_coverage = ONE_EIGHT;
    //final double NO_HOTSPOT_Nucleotide_coverage = _26_OVER_32;
    
    //new coverage values are based on the germline frequency of the hotspots
    //for 59982 4mers in all germline IGHV, 2284 RGYW, 2116 WRCY, 4349 WAN and 51233 Non-HS
    final double RGYW_Nucleotide_coverage = 0.0381D;
    final double WRCY_Nucleotide_coverage = 0.0353D;
    final double WAN_Nucleotide_coverage = 0.0725D;
    final double NO_HOTSPOT_Nucleotide_coverage = 0.854D;
    

    //final double RGYW_mutation_prob = 0.13800000000000001D;
    //final double WRCY_mutation_prob = 0.096000000000000002D;
    //final double WAN_mutation_prob = 0.161D;
    //final double NO_HOTSPOT_mutation_prob = 0.60499999999999998D;

    //mutation probabilities calculated from Barington 2007 JI paper
    //4592 muts total, 598 RGYW, 633 WAN, 558 WRCY, 2803 in non-hotspots
    final double RGYW_mutation_prob = 0.130D;
    final double WRCY_mutation_prob = 0.122D;
    final double WAN_mutation_prob = 0.138D;
    final double NO_HOTSPOT_mutation_prob = 0.610D;

    final double NO_HOTSPOT_MUTABILITY_SCORE = (NO_HOTSPOT_mutation_prob/NO_HOTSPOT_Nucleotide_coverage);
    final double WAN_MUTABILITY_SCORE = (WAN_mutation_prob/WAN_Nucleotide_coverage);
    final double NAN_MUTABILITY_SCORE = ((WAN_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
    final double RGYW_MUTABILITY_SCORE = (RGYW_mutation_prob / RGYW_Nucleotide_coverage);
    final double NGYW_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
    final double RGYN_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
    final double RGNN_MUTABILITY_SCORE = ((RGYW_MUTABILITY_SCORE*ONE_QUARTER) + (NO_HOTSPOT_MUTABILITY_SCORE*THREE_QUARTER));
    final double WRCY_MUTABILITY_SCORE = (WRCY_mutation_prob/WRCY_Nucleotide_coverage);
    final double NRCY_MUTABILITY_SCORE = ((WRCY_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
    final double NNCY_MUTABILITY_SCORE = (WRCY_MUTABILITY_SCORE * ONE_QUARTER)+(NO_HOTSPOT_MUTABILITY_SCORE*THREE_QUARTER);
    final double WRCN_MUTABILITY_SCORE = ((WRCY_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2);
    boolean DEBUGGING;

    public MutabilityScoreHotspot()
    {
        DEBUGGING = false;
    }

    public double pentaNucleotideScore(String pentaNucleotideUnknownCase)
    {
        String pentaNucleotide = pentaNucleotideUnknownCase.toLowerCase();
        char centerNucleotide = pentaNucleotide.charAt(2);
        switch(centerNucleotide)
        {
        case 116: // 't'
            if(DEBUGGING)
            {
                System.out.println("Main: No Hot");
            }
            return NO_HOTSPOT_MUTABILITY_SCORE;

        case 97: // 'a'
            return testWAN(pentaNucleotide);

        case 103: // 'g'
            return testRGYW(pentaNucleotide);

        case 99: // 'c'
            return testWRCY(pentaNucleotide);
        default: // if it isn't an a/c/g then returns same as t -> changed because else iHMMune align throws error for seqs with 'n'
            return NO_HOTSPOT_MUTABILITY_SCORE;
        }
        //throw new Error("pentaNucleotideScore(): center nucleotide not a,c,g or t, but : " + centerNucleotide);
    }

    private double testWAN(String pentaNucleotide)
    {
        char secondNucleotide = pentaNucleotide.charAt(1);
        if(isW(secondNucleotide))
        {
            if(DEBUGGING)
            {
                System.out.println("WAN");
            }
            return WAN_MUTABILITY_SCORE;
        }
        if(isN(secondNucleotide))
        {
            if(DEBUGGING)
            {
                System.out.println("NAN");
            }
            return NAN_MUTABILITY_SCORE;
        }
        if(DEBUGGING)
        {
            System.out.println("WAN: No Hot");
        }
        return NO_HOTSPOT_MUTABILITY_SCORE;
    }

    private double testRGYW(String pentaNucleotide)
    {
        char secondNucleotide = pentaNucleotide.charAt(1);
        char fourthNucleotide = pentaNucleotide.charAt(3);
        char fifthNucleotide = pentaNucleotide.charAt(4);
        if(isR(secondNucleotide) && isY(fourthNucleotide) && isW(fifthNucleotide))
        {
            return RGYW_MUTABILITY_SCORE;
        }
        if(isN(secondNucleotide))
        {
            if(isY(fourthNucleotide) && isW(fifthNucleotide))
            {
                if(DEBUGGING)
                {
                    System.out.println("NGYW");
                }
                return NGYW_MUTABILITY_SCORE;
            }
        } else
        if(isN(fourthNucleotide))
        {
            if(isR(secondNucleotide))
            {
                if(DEBUGGING)
                {
                    System.out.println("RGNN");
                }
                return RGNN_MUTABILITY_SCORE;
            }
        } else
        if(isN(fifthNucleotide) && isR(secondNucleotide) && isY(fourthNucleotide))
        {
            if(DEBUGGING)
            {
                System.out.println("RGYN");
            }
            return RGYN_MUTABILITY_SCORE;
        }
        if(DEBUGGING)
        {
            System.out.println("RGYW: No Hot");
        }
        return NO_HOTSPOT_MUTABILITY_SCORE;
    }

    private double testWRCY(String pentaNucleotide)
    {
        char firstNucleotide = pentaNucleotide.charAt(0);
        char secondNucleotide = pentaNucleotide.charAt(1);
        char fourthNucleotide = pentaNucleotide.charAt(3);
        if(isW(firstNucleotide) && isR(secondNucleotide) && isY(fourthNucleotide))
        {
            if(DEBUGGING)
            {
                System.out.println("WRCY");
            }
            return WRCY_MUTABILITY_SCORE;
        }
        if(isN(fourthNucleotide))
        {
            if(isW(firstNucleotide) && isR(secondNucleotide))
            {
                if(DEBUGGING)
                {
                    System.out.println("WRCN");
                }
                return WRCN_MUTABILITY_SCORE;
            }
        } else
        if(isN(secondNucleotide))
        {
            if(isY(fourthNucleotide))
            {
                if(DEBUGGING)
                {
                    System.out.println("NNCY");
                }
                return NNCY_MUTABILITY_SCORE;
            }
        } else
        if(isN(firstNucleotide) && isR(secondNucleotide) && isY(fourthNucleotide))
        {
            if(DEBUGGING)
            {
                System.out.println("NRCY");
            }
            return NRCY_MUTABILITY_SCORE;
        }
        if(DEBUGGING)
        {
            System.out.println("WRCY: No Hot");
        }
        return NO_HOTSPOT_MUTABILITY_SCORE;
    }

    private boolean isN(char nucleotide)
    {
        return nucleotide == 'n';
    }

    private boolean isA(char nucleotide)
    {
        return nucleotide == 'a';
    }

    private boolean isC(char nucleotide)
    {
        return nucleotide == 'c';
    }

    private boolean isG(char nucleotide)
    {
        return nucleotide == 'g';
    }

    private boolean isT(char nucleotide)
    {
        return nucleotide == 't';
    }

    private boolean isR(char nucleotide)
    {
        return isA(nucleotide) || isG(nucleotide);
    }

    private boolean isY(char nucleotide)
    {
        return isC(nucleotide) || isT(nucleotide);
    }

    private boolean isW(char nucleotide)
    {
        return isA(nucleotide) || isT(nucleotide);
    }

    public static void main(String args[])
    {
        MutabilityScoreHotspot object = new MutabilityScoreHotspot();
        System.out.println("dividing two integers: 5 / 2 = 2.5");
        double mutabilityScore = object.pentaNucleotideScore("naacg");
        object.getClass();
        System.out.println("ONE_EIGTH mutability score = " + 0.125D);
        object.getClass();
        System.out.println("WAN mutability score = " + 1.288D);
        System.out.println("mutability score = " + mutabilityScore);
    }
}
