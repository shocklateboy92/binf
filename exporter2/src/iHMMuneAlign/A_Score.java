package iHMMuneAlign;

import java.io.FileWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
/**
 * this class is a container for one or more static methods used to
 * calculate the A value for a given VGene
 * The A-value for an alignment is the probability of an alignment in the VGene
 * at pos 0, this probability decays exponentially as we move away from
 * the start of the VGene
 */

public class A_Score
{
    	/**
     * the A value used for the CDR prediction is found
     * by looking at the number of mutations in the common
     * region of a sequence over the number of nucleotides in
     * this common region multiplied by the A_slope for this
     * common region.
     *
     * The A_slope is previously determined by using linear regression
     * on the A values (Y's) over the mutation ratio for the common region (X's)
     * and thus estimating the relationship between A value and number of 
     * of mutations for common region
     */

    static final char N_NUCLEOTIDE = 110;
    static final char X_NUCLEOTIDE = 120;
    static final int COMMON_AREA_START_POS = 129;
    static final int COMMON_AREA_END_POS = 243;
    static final int COMMON_AREA_NUCL_LENGTH = 115;
    static final double A_SLOPE = 0.0064999999999999997D;
    static final double B_ADDITION = 0.0015D;
    static final boolean DEBUGGING = false;
    
    /**
     * Data for dumping
     */
    
    int number_of_mutations_in_CR;
    
    int number_of_NandX_nts_in_CR;
    
    int VGene_start_offset;
    
    int vGeneLength;
    
    double aProbability;
    	
    
     /**
     * calculate the A probability using the above rules and
     * the number of mutations found in the V-Gene "VGene" when
     * compared to it's mathcing sequence "matching sequence", the
     * "VGene_start_offset" must be taken into account when placing
     * the common region
     * 
     * Pre: 	"VGene" and "matching_sequence" are aligned
     *
     * returns		A probability based on number of mutations in VGene's common area plus above rules
     * 
     * Errors	if VGene does not contain the entire common region, and
     *			if matching sequence is not at least the length of the VGene
     */

    public A_Score()
    {
        this.number_of_mutations_in_CR = 0;
        
        this.number_of_NandX_nts_in_CR = 0;
        
        this.VGene_start_offset = 0;
        
        this.vGeneLength = 0;
        
        this.aProbability = 0;
    }

    public double A_probabilityLong(String VGeneUnknownCase, String matching_sequenceUnknownCase, int VGene_start_offset, PrintWriter pw)
    {
        String VGene = VGeneUnknownCase.toLowerCase();
        String matching_sequence = matching_sequenceUnknownCase.toLowerCase();
        if (DEBUGGING) {
        	System.out.println("VGene length = " + VGene.length());
        	System.out.println("Matching seq. length = " + matching_sequenceUnknownCase.length());
        }
        if(VGene.length() + VGene_start_offset < (COMMON_AREA_END_POS))
        {
            //throw new Error("A_Score: A_probability(): VGene does not contain common region (VGene too short)");
            return 2;
        }
        if(VGene_start_offset > COMMON_AREA_START_POS)
        {
            //throw new Error("A_Score: A_probability(): VGene does not contain common region (VGene offset too great)");
            return 3;
        }
        if(VGene.length() > matching_sequence.length())
        {
            //throw new Error("A_Score: A_probability(): VGene length is greater than its matching sequence length");
            return 4;
        }
        int mutationsInVGene = mutationsInVGene(VGeneUnknownCase, matching_sequenceUnknownCase, VGene_start_offset, pw);
        int vGeneLength = VGeneUnknownCase.length();
        double CRtoVGeneLengthRatio = COMMON_AREA_NUCL_LENGTH / (double)vGeneLength;
        double averageMutationsInCR = (double)mutationsInVGene * CRtoVGeneLengthRatio;
        int number_of_mutations_in_CR = 0;
        int number_of_NandX_nts_in_CR = 0;
        for(int i = (COMMON_AREA_START_POS-1) - VGene_start_offset; i <= ((COMMON_AREA_END_POS-1)) - VGene_start_offset; i++)
        {        //*128*
        	if (DEBUGGING) {
	        	System.out.println(VGene.length() + " vs " + i);
	        }
            char V_nucl = VGene.charAt(i);
            char MS_nucl = matching_sequence.charAt(i);
            if(V_nucl != MS_nucl)
            {
                if(MS_nucl == 'n' || MS_nucl == 'x')
                {
                    number_of_NandX_nts_in_CR++;
                } else
                {
                    number_of_mutations_in_CR++;
                }
            }
        }
		
		if (DEBUGGING) {
	        System.out.println("number of mutations found in common area = " + number_of_mutations_in_CR);
    	    System.out.println("number of n/x nts. found in common area = " + number_of_NandX_nts_in_CR);
    	}
        if(pw != null)
        {
            pw.println(number_of_mutations_in_CR + "\t" + COMMON_AREA_NUCL_LENGTH);
        }
        double probabilityMutationPlusB = (double)number_of_mutations_in_CR * A_SLOPE + B_ADDITION;
        double NandXaddition = (averageMutationsInCR / COMMON_AREA_NUCL_LENGTH) * A_SLOPE * (double)number_of_NandX_nts_in_CR;
        double aProbability = probabilityMutationPlusB + NandXaddition;

        //store: number_of_mutations_in_CR, 
        //number_of_NandX_nts_in_CR, aProbability
        //VGene_start_offset, vGeneLength
        //for write out to ouput file
        
        this.number_of_mutations_in_CR = number_of_mutations_in_CR;
        
        this.number_of_NandX_nts_in_CR = number_of_NandX_nts_in_CR;
        
        this.VGene_start_offset = VGene_start_offset;
        
        this.vGeneLength = vGeneLength;
        
        this.aProbability = aProbability;
        
        return aProbability;
    }
    
     public double A_probabilityShort(String VGeneUnknownCase, String matching_sequenceUnknownCase, int VGene_start_offset, PrintWriter pw)
    {
        String VGene = VGeneUnknownCase.toLowerCase();
        String matching_sequence = matching_sequenceUnknownCase.toLowerCase();
        if (DEBUGGING) {
	        System.out.println("VGene length = " + VGene.length());
    	    System.out.println("Matching seq. length = " + matching_sequenceUnknownCase.length());
    	}
        if (VGene.length() < 100) { //at least 100 nucleotides are length is required to calculate
        	//there is insufficient VGene length for the calculation, so throw an error
        	//throw new Error("A_Score: A_probabilityShort(): There is insufficient VGene length to calculate the A_score (VGene is too short)");
            return 2;
        }
        if(VGene.length() > matching_sequence.length())
        {
            //throw new Error("A_Score: A_probability(): VGene length is greater than its matching sequence length");
            return 4;
        }

        int number_of_mutations_in_CR = 0;
        int number_of_NandX_nts_in_CR = 0;
        //determine the offset for the start of the common region given the VgeneStartOffset for this sequence
        int commonRegionStartOffset = VGene_start_offset - COMMON_AREA_START_POS;
        
        int amendedCommonAreaStart = 0;
        int amendedCommonAreaEnd = 0;
        int amendedCommonAreaNuclLength = COMMON_AREA_NUCL_LENGTH;
        if (commonRegionStartOffset > 0) {
        	//start from the first nucleotide in the truncated common region for the sequence
        	amendedCommonAreaStart = 0;
        	//adjusted the common region end position to include extra length equal to the nucleotides lost from the 5` truncation
        	amendedCommonAreaEnd = COMMON_AREA_END_POS + commonRegionStartOffset - VGene_start_offset;
        	
        	if ((amendedCommonAreaEnd > VGene.length()) && (VGene.length() >= 100)) {
        		//the new common area end is now out-of-bounds of the sequence length
        		//if there are at least 100 nts in the common area, just use the full length of the gene for the calculation
        		amendedCommonAreaEnd = VGene.length();
        		amendedCommonAreaNuclLength = VGene.length();	
        	} 
        }
        
        
        for(int i = amendedCommonAreaStart; i < amendedCommonAreaEnd; i++)
        {        //*128*
        	if (DEBUGGING) {
	        	System.out.println(VGene.length() + " vs " + i);
	        }
            char V_nucl = VGene.charAt(i);
            char MS_nucl = matching_sequence.charAt(i);
            if(V_nucl != MS_nucl)
            {
                if(MS_nucl == 'n' || MS_nucl == 'x')
                {
                    number_of_NandX_nts_in_CR++;
                } else
                {
                    number_of_mutations_in_CR++;
                }
            }
        }
		int mutationsInVGene = mutationsInVGene(VGeneUnknownCase, matching_sequenceUnknownCase, VGene_start_offset, pw);
        int vGeneLength = VGeneUnknownCase.length();
        double CRtoVGeneLengthRatio = amendedCommonAreaNuclLength / (double)vGeneLength;
        double averageMutationsInCR = (double)mutationsInVGene * CRtoVGeneLengthRatio;
        if (DEBUGGING) {
		    System.out.println("number of mutations found in common area = " + number_of_mutations_in_CR);
		    System.out.println("number of n/x nts. found in common area = " + number_of_NandX_nts_in_CR);
        }
        if(pw != null)
        {
            pw.println(number_of_mutations_in_CR + "\t" + amendedCommonAreaNuclLength);
        }
        double probabilityMutationPlusB = (double)number_of_mutations_in_CR * A_SLOPE + B_ADDITION;
        double NandXaddition = (averageMutationsInCR / amendedCommonAreaNuclLength) * A_SLOPE * (double)number_of_NandX_nts_in_CR;
        double aProbability = probabilityMutationPlusB + NandXaddition;

        //store: number_of_mutations_in_CR, 
        //number_of_NandX_nts_in_CR, aProbability
        //VGene_start_offset, vGeneLength
        //for write out to ouput file
        
        this.number_of_mutations_in_CR = number_of_mutations_in_CR;
        
        this.number_of_NandX_nts_in_CR = number_of_NandX_nts_in_CR;
        
        this.VGene_start_offset = VGene_start_offset;
        
        this.vGeneLength = vGeneLength;
        
        this.aProbability = aProbability;
        
        return aProbability;
    }
    
    
    public void dumpInfo(String htmlFileName) {   	
    	PrintWriter outStream = null;
    	try
        {
            outStream = new PrintWriter(new FileWriter(htmlFileName, true));
        }
        catch(Exception ex)
        {
            System.out.println("couldn't open new file: " + ex.getMessage());
        }
        
        outStream.print(";\""+this.number_of_mutations_in_CR+"\";\""+this.number_of_NandX_nts_in_CR+"\";\"");
        outStream.println(this.VGene_start_offset+"\";\""+this.vGeneLength+"\";\""+this.aProbability+"\"");
        
        outStream.close();    
    }

    private static int mutationsInVGene(String VGeneUnknownCase, String matching_sequenceUnknownCase, int VGene_start_offset, PrintWriter pw)
    {
        String VGene = VGeneUnknownCase.toLowerCase();
        String matching_sequence = matching_sequenceUnknownCase.toLowerCase();
        int number_of_mutations_in_VGene = 0;
        int number_of_NandX_nts_in_VGene = 0;
        for(int i = 0; i < VGene.length(); i++)
        {
            char V_nucl = VGene.charAt(i);
            char MS_nucl = matching_sequence.charAt(i);
            if(V_nucl != MS_nucl)
            {
                if(MS_nucl == 'n' || MS_nucl == 'x')
                {
                    number_of_NandX_nts_in_VGene++;
                } else
                {
                    number_of_mutations_in_VGene++;
                }
            }
        }

        int vGeneLength = VGene.length();
        if (DEBUGGING) {
        	System.out.println("number of mutations found in VGene area = " + number_of_mutations_in_VGene);
        	System.out.println("number of n/x nts. found in VGene  area = " + number_of_NandX_nts_in_VGene);
        	System.out.println("VGene area length = " + vGeneLength);
        }
        if(pw != null)
        {
            pw.print(number_of_mutations_in_VGene + "\t" + VGene_start_offset + "\t" + vGeneLength + "\t");
        }
        return number_of_mutations_in_VGene;
    }

    public static void main(String args[])
    {
    	A_Score myA_Score = new A_Score();
        double A_prob = myA_Score.A_probabilityLong("aaacccgggtttaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" +
"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" +
"aaaaaagaatacaa"
, "aaacccaaatttaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" +
"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" +
"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaggggggggggggggggggggggggggggggg"
, 128, null);
        System.out.println("A_prob = " + A_prob);
    }
}
