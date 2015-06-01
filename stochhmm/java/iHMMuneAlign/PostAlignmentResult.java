package iHMMuneAlign;

import java.io.PrintStream;

class PostAlignmentResult
{

    private String UMS_string;
    private String VGene_string;
    private String UMS_name;
    private String VGene_name;
    private int row_offset;
    private int column_offset;
    private int gaps_in_alignment;
    private int similarity;
    //adding a value for the score
    private float align_score;
    //adding a value for the aligned input sequence
    private char[] aligned_input;

    public PostAlignmentResult(String UMS_string, String VGene_string, int row_offset, int column_offset, int gaps_in_alignment, int similarity, String UMS_name, 
            String VGene_name, float align_score, char[] aligned_input)
    {
        this.UMS_string = UMS_string;
        this.VGene_string = VGene_string;
        this.row_offset = row_offset;
        this.column_offset = column_offset;
        this.gaps_in_alignment = gaps_in_alignment;
        this.similarity = similarity;
        this.UMS_name = UMS_name;
        this.VGene_name = VGene_name;
        this.align_score = align_score;
        this.aligned_input = aligned_input;
    }
	
	//add a function to get the similarity
	public int getSimilarity()
	{
		int similarityFromAlignment = similarity;
		return similarityFromAlignment;
	}

    //add a function to get the alingment score
    public float getAlignmentScore()
    {
        float ScoreFromAlignment = align_score;
        return ScoreFromAlignment;
    }

    //add a function to get the gaps
    public int getAlignmentGaps() {
        int GapsFromAlignment = gaps_in_alignment;
        return GapsFromAlignment;
    }
    
    //add a function to get the aligned input sequence as returned from jaligner
    public char[] getAlignedInput() {
        char[] AlignedInputSeq = aligned_input;
        return AlignedInputSeq;
    }

    public String getUMSfromAlignmentString()
    {
        String UMSfromAlignment = UMS_string.substring(row_offset, UMS_string.length());
        return UMSfromAlignment;
    }

    public String getVGeneAlignmentString()
    {
        String VGeneFromAlignment = VGene_string.substring(column_offset, VGene_string.length());
		//System.out.println("V gene offset: " + column_offset);
		//System.out.println("V gene string: " + VGene_string);
		//System.out.println("V gene length: " + VGene_string.length());
        return VGeneFromAlignment;
    }
    
    public String getJGeneAlignmentString() {
    	//return a the germline gene string with any insertions or deletions accounted for
    	//for the trailingJVectorFinder the best gene is not the full sequence as for the BestVGeneFinder and rather just the hseq
    	//therefore just need to replace any '-' from the hseq and return the string
    	String JGeneGermlineFromAlignment = VGene_string;
    	//JGeneGermlineFromAlignment = JGeneGermlineFromAlignment.replace("-","n");
    	    	
    	return JGeneGermlineFromAlignment;
    }
    
    public String getFullLengthUMSString() {
    	return UMS_string;
    }
    
    public String getFullLengthVGeneString() {
    	return VGene_string;
    }

    public String getUMSname()
    {
        return UMS_name;
    }

    public String getVGeneName()
    {
        return VGene_name;
    }

    public int getColOffset()
    {
        return column_offset;
    }

    public int getRowOffset()
    {
        return row_offset;
    }

    public void displayAlignmentInfo()
    {
        System.out.println();
        System.out.println();
        System.out.println();
        System.out.println("*** VGene ALIGNMENT INFO ***");
        System.out.println();
        System.out.println("complete input string length = " + UMS_string.length());
        System.out.println("complete input string:");
        System.out.println(UMS_string);
        System.out.println();
		System.out.println("position from which input aligns to germline (offset): " + row_offset);
        String from_alignment_UMS = getUMSfromAlignmentString();
        System.out.println("length of input that aligns to germline sequence: " + from_alignment_UMS.length());
        System.out.println("portion of input sequence that occurs downstream of first match to germline gen:");
        System.out.println(from_alignment_UMS);
        System.out.println();
        System.out.println("VGene: " + VGene_name);
        System.out.println("position from which germline gene sequence aligns to input (offset): " + column_offset);
		String from_alignment_VGene = getVGeneAlignmentString();
        System.out.println("length of germline gene that aligns against input: " + from_alignment_VGene.length());
        System.out.println("portion of germline gene that aligns against input:");
        System.out.println(from_alignment_VGene);
        System.out.println();
        System.out.println("similarity between input and germline gene: " + similarity);
        System.out.println("Gaps in alignment: " + gaps_in_alignment);
        //add a line to output the score
        System.out.println("Alignment Score: " + align_score);
        System.out.println();
        System.out.println();
    }
}
