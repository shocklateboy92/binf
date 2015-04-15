package iHMMuneAlign;

import java.io.PrintStream;

// Referenced classes of package iHMMuneAlign:
//            MutabilityScore

public class MutabilityScoreTrinucleotide
    implements MutabilityScore
{

    final char U = 'u';
    final char A = 'a';
    final char C = 'c';
    final char G = 'g';
    final char T = 't';
    final int NUCL_A_NUMBER = 0;
    final int NUCL_C_NUMBER = 1;
    final int NUCL_G_NUMBER = 2;
    final int NUCL_T_NUMBER = 3;
    static final double TRINUCLEOTIDE_MUTABILITY_SCORE_ARRAY[] = {
        0.65000000000000002D, 1.3799999999999999D, 1.1499999999999999D, 2.0099999999999998D, 1.0600000000000001D, 1.1399999999999999D, 1.5D, 1.29D, 0.91000000000000003D, 2.1200000000000001D, 
        0.70999999999999996D, 1.3200000000000001D, 1.46D, 1.3D, 1.1799999999999999D, 1.6499999999999999D, 1.1699999999999999D, 1.1299999999999999D, 1.51D, 1.0900000000000001D, 
        1.0800000000000001D, 0.46000000000000002D, 0.75D, 0.56999999999999995D, 0.64000000000000001D, 0.56000000000000005D, 0.62D, 1.25D, 1.6699999999999999D, 0.59999999999999998D, 
        0.85999999999999999D, 0.64000000000000001D, 1.05D, 0.48999999999999999D, 0.90000000000000002D, 0.87D, 1.73D, 0.48999999999999999D, 0.29999999999999999D, 2.1499999999999999D, 
        0.70999999999999996D, 0.76000000000000001D, 0.70999999999999996D, 0.78000000000000003D, 3.0600000000000001D, 0.65000000000000002D, 0.95999999999999996D, 1.29D, 1.9299999999999999D, 1.76D, 
        2.1299999999999999D, 1.6100000000000001D, 0.91000000000000003D, 0.68999999999999995D, 0.46999999999999997D, 0.60999999999999999D, 0.31D, 0.60999999999999999D, 0.57999999999999996D, 1.01D, 
        1.3200000000000001D, 0.60999999999999999D, 0.28000000000000003D, 0.47999999999999998D
    };

    public static void main(String args[])
    {
        MutabilityScoreTrinucleotide mutability_score_object = new MutabilityScoreTrinucleotide();
        System.out.println("TRINUCLEOTIDE_MUTABILITY_SCORE_ARRAY size = " + TRINUCLEOTIDE_MUTABILITY_SCORE_ARRAY.length);
        double mutability_score = mutability_score_object.pentaNucleotideScore("cgtcc");
        System.out.println("mutability_score = " + mutability_score);
    }

    public MutabilityScoreTrinucleotide()
    {
    }

    public double pentaNucleotideScore(String pentaNucleotide_any_case)
    {
        if(pentaNucleotide_any_case == null)
        {
            throw new Error("penta nucleotide is null");
        }
        String pentaNucleotide = pentaNucleotide_any_case.toLowerCase();
        if(pentaNucleotide.length() != 5)
        {
            throw new Error("pentaNucleotideScore(): pentaNucleotide is not of length 5");
        } else
        {
            int TRINUCLEOTIDES_IN_PENTANUCLEOTIDE = 3;
            String first_trinucleotide = pentaNucleotide.substring(0, 3);
            String middle_trinucleotide = pentaNucleotide.substring(1, 4);
            String last_trinucleotide = pentaNucleotide.substring(2, 5);
            double first_trinucleotide_score = triNucleotideScore(first_trinucleotide);
            double middle_trinucleotide_score = triNucleotideScore(middle_trinucleotide);
            double last_trinucleotide_score = triNucleotideScore(last_trinucleotide);
            double total_trinucleotide_score = first_trinucleotide_score + middle_trinucleotide_score + last_trinucleotide_score;
            double average_trinucleotide_score = total_trinucleotide_score / 3D;
            return average_trinucleotide_score;
        }
    }

    private double triNucleotideScore(String trinucleotide)
    {
        double trinucleotide_score = 4.9406564584124654E-324D;
        if(isNNN(trinucleotide))
        {
            trinucleotide_score = lookupTriNucleotideMutabilityScore(trinucleotide);
        } else
        if(isNNU(trinucleotide))
        {
            trinucleotide_score = NNU_TriNucleotideScore(trinucleotide);
        } else
        if(isNUU(trinucleotide))
        {
            trinucleotide_score = NUU_TriNucleotideScore(trinucleotide);
        } else
        if(isUNN(trinucleotide))
        {
            trinucleotide_score = UNN_TriNucleotideScore(trinucleotide);
        } else
        if(isUUN(trinucleotide))
        {
            trinucleotide_score = UUN_TriNucleotideScore(trinucleotide);
        }
        return trinucleotide_score;
    }

    private double lookupTriNucleotideMutabilityScore(String strict_trinucleotide)
    {
        char charOne = strict_trinucleotide.charAt(0);
        int numberOne = getNucleotideNumber(charOne);
        char charTwo = strict_trinucleotide.charAt(1);
        int numberTwo = getNucleotideNumber(charTwo);
        char charThree = strict_trinucleotide.charAt(2);
        int numberThree = getNucleotideNumber(charThree);
        int trinucleotide_mutatbility_score_array_index = numberOne * 16 + numberTwo * 4 + numberThree;
        double mutability_score = TRINUCLEOTIDE_MUTABILITY_SCORE_ARRAY[trinucleotide_mutatbility_score_array_index];
        return mutability_score;
    }

    private int getNucleotideNumber(char nucleotide)
    {
        switch(nucleotide)
        {
        case 97: // 'a'
            return 0;

        case 99: // 'c'
            return 1;

        case 103: // 'g'
            return 2;

        case 116: // 't'
            return 3;
        }
        throw new Error("NucleotideNumber: nucleotide is not of type a,c,g ot t");
    }

    private char getNucleotideFromNumber(int nucleotide_number)
    {
        switch(nucleotide_number)
        {
        case 0: // '\0'
            return 'a';

        case 1: // '\001'
            return 'c';

        case 2: // '\002'
            return 'g';

        case 3: // '\003'
            return 't';
        }
        throw new Error("getNucleotideFromNumber: nucleotide number is not in range");
    }

    private boolean isNNN(String trinucleotide)
    {
        char charOne = trinucleotide.charAt(0);
        if(!isN(charOne))
        {
            return false;
        }
        char charTwo = trinucleotide.charAt(1);
        if(!isN(charTwo))
        {
            return false;
        }
        char charThree = trinucleotide.charAt(2);
        return isN(charThree);
    }

    private boolean isNNU(String trinucleotide)
    {
        char charOne = trinucleotide.charAt(0);
        if(!isN(charOne))
        {
            return false;
        }
        char charTwo = trinucleotide.charAt(1);
        if(!isN(charTwo))
        {
            return false;
        }
        char charThree = trinucleotide.charAt(2);
        return isU(charThree);
    }

    private boolean isNUU(String trinucleotide)
    {
        char charOne = trinucleotide.charAt(0);
        if(!isN(charOne))
        {
            return false;
        }
        char charTwo = trinucleotide.charAt(1);
        if(!isU(charTwo))
        {
            return false;
        }
        char charThree = trinucleotide.charAt(2);
        return isU(charThree);
    }

    private boolean isUNN(String trinucleotide)
    {
        char charOne = trinucleotide.charAt(0);
        if(!isU(charOne))
        {
            return false;
        }
        char charTwo = trinucleotide.charAt(1);
        if(!isN(charTwo))
        {
            return false;
        }
        char charThree = trinucleotide.charAt(2);
        return isN(charThree);
    }

    private boolean isUUN(String trinucleotide)
    {
        char charOne = trinucleotide.charAt(0);
        if(!isU(charOne))
        {
            return false;
        }
        char charTwo = trinucleotide.charAt(1);
        if(!isU(charTwo))
        {
            return false;
        }
        char charThree = trinucleotide.charAt(2);
        return isN(charThree);
    }

    private boolean isN(char trinucleotide)
    {
        return trinucleotide == 'a' || trinucleotide == 'c' || trinucleotide == 'g' || trinucleotide == 't';
    }

    private boolean isU(char trinucleotide)
    {
        return trinucleotide == 'u';
    }

    private double NNU_TriNucleotideScore(String NNU_trinucleotide)
    {
        String temp_trinucleotide = "";
        double total_trinucleotide_score = 0.0D;
        double TOTAL_TRINUCLEOTIDE_COMBINATIONS = 4D;
        for(int nucl_number = 0; nucl_number <= 3; nucl_number++)
        {
            char replacement_nucl = getNucleotideFromNumber(nucl_number);
            temp_trinucleotide = "" + NNU_trinucleotide.charAt(0) + NNU_trinucleotide.charAt(1) + replacement_nucl;
            double temp_trinucleotide_score = lookupTriNucleotideMutabilityScore(temp_trinucleotide);
            total_trinucleotide_score += temp_trinucleotide_score;
        }

        double average_trinucleotide_score = total_trinucleotide_score / 4D;
        return average_trinucleotide_score;
    }

    private double NUU_TriNucleotideScore(String NUU_trinucleotide)
    {
        int TOTAL_NXU_TRINUCLEOTIDE_COMBINATIONS = 4;
        double total_NXU_trinucleotide_score = 0.0D;
        String temp_NXU_trinucleotide = "";
        double total_trinucleotide_score = 0.0D;
        for(int nucl_number = 0; nucl_number <= 3; nucl_number++)
        {
            char middle_U_replacement_nucl = getNucleotideFromNumber(nucl_number);
            temp_NXU_trinucleotide = "" + NUU_trinucleotide.charAt(0) + middle_U_replacement_nucl + NUU_trinucleotide.charAt(2);
            double temp_NXU_trinucleotide_score = NNU_TriNucleotideScore(temp_NXU_trinucleotide);
            total_NXU_trinucleotide_score += temp_NXU_trinucleotide_score;
        }

        double average_trinucleotide_score = total_NXU_trinucleotide_score / 4D;
        return average_trinucleotide_score;
    }

    private double UNN_TriNucleotideScore(String UNN_trinucleotide)
    {
        String temp_trinucleotide = "";
        double total_trinucleotide_score = 0.0D;
        double TOTAL_TRINUCLEOTIDE_COMBINATIONS = 4D;
        for(int nucl_number = 0; nucl_number <= 3; nucl_number++)
        {
            char replacement_nucl = getNucleotideFromNumber(nucl_number);
            temp_trinucleotide = "" + replacement_nucl + UNN_trinucleotide.charAt(1) + UNN_trinucleotide.charAt(2);
            double temp_trinucleotide_score = lookupTriNucleotideMutabilityScore(temp_trinucleotide);
            total_trinucleotide_score += temp_trinucleotide_score;
        }

        double average_trinucleotide_score = total_trinucleotide_score / 4D;
        return average_trinucleotide_score;
    }

    private double UUN_TriNucleotideScore(String UUN_trinucleotide)
    {
        int TOTAL_UXN_TRINUCLEOTIDE_COMBINATIONS = 4;
        double total_UXN_trinucleotide_score = 0.0D;
        String temp_UXN_trinucleotide = "";
        double total_trinucleotide_score = 0.0D;
        for(int nucl_number = 0; nucl_number <= 3; nucl_number++)
        {
            char middle_U_replacement_nucl = getNucleotideFromNumber(nucl_number);
            temp_UXN_trinucleotide = "" + UUN_trinucleotide.charAt(0) + middle_U_replacement_nucl + UUN_trinucleotide.charAt(2);
            double temp_UXN_trinucleotide_score = UNN_TriNucleotideScore(temp_UXN_trinucleotide);
            total_UXN_trinucleotide_score += temp_UXN_trinucleotide_score;
        }

        double average_trinucleotide_score = total_UXN_trinucleotide_score / 4D;
        return average_trinucleotide_score;
    }

}
