package iHMMuneAlign;

import java.io.PrintStream;

class BLASTResult
{

    private String inputStr;
    private String germlineStr;
    private String inputName;
    private String germlineName;
    private int input_offset; //row_offset
    private int germline_offset;
    private int gaps_in_alignment;
    private int similarity;
    //adding a value for the score
    private float align_score;
    //adding a value for the aligned input sequence
    private char[] aligned_input;
	//adding vars for the full BLAST results
	private float bitScore;
	private float evalue;
	private int inputEnd;
	private int germlineEnd;
	private int inputFrame;
	private int identity;
	private int positive;
	private int alignLen;
	private char[] germlineAligned;
	private char[] midLine;
	
	//object representing BLAST result that mirrors the Alignment results object for iHMMune Align
    /*public BLASTResult(String inputStr, String germlineStr, int input_offset, int germline_offset, int gaps_in_alignment, int similarity, String inputName, 
            String germlineName, float align_score, char[] aligned_input)
    {
        this.inputStr = inputStr;
        this.germlineStr = germlineStr;
        this.input_offset = input_offset;
        this.germline_offset = germline_offset;
        this.gaps_in_alignment = gaps_in_alignment;
        this.similarity = similarity;
        this.inputName = inputName;
        this.germlineName = germlineName;
        this.align_score = align_score;
        this.aligned_input = aligned_input;
    }*/
	
	//object representing all BLAST alignment result fields
	public BLASTResult(String inputStr, float bitScore, float score, float evalue, int inputFrom, int inputTo, int hitFrom, int hitTo, int queryFrame, int identity, int positive, int gaps, int alignLen, char[] qSeq, char[] hSeq, char[] midLine, String hit_def) 
	{
		this.inputStr = inputStr;
		this.bitScore = bitScore;
		this.align_score = score;
		this.evalue = evalue;
		this.input_offset = inputFrom;
		this.inputEnd = inputTo;
		this.germline_offset = hitFrom;
		this.germlineEnd = hitTo;
		this.inputFrame = queryFrame;
		this.identity = identity;
		this.positive = positive;
		this.gaps_in_alignment = gaps;
		this.alignLen = alignLen;
		this.aligned_input = qSeq;
		this.germlineAligned = hSeq;
		this.midLine = midLine;
		this.germlineName = hit_def;
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
	
	//add a function to get the bit score
	public float getBitScore()
	{
		return bitScore;
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

    public String getInputStrFromAlignment()
    {
        String tmpStr = aligned_input.toString();
		return tmpStr;
    }

    public String getGermlineStrFromAlignment()
    {
        String tmpStr = germlineAligned.toString();
		return tmpStr;
    }

    public String getInputName()
    {
        return inputName;
    }

    public String getGermlineName()
    {
        return germlineName;
    }

    public int getInputOffset()
    {
        return input_offset;
    }

    public int getGermlineOffset()
    {
        return germline_offset;
    }

    public void displayAlignmentInfo()
    {
        System.out.println();
        System.out.println();
        System.out.println();
        System.out.println("*** BLAST ALIGNMENT INFO ***");
        System.out.println();
        System.out.println("complete input string length = " + inputStr.length());
        System.out.println("complete input string:");
        System.out.println(inputStr);
        System.out.println();
        String from_alignment_input = getInputStrFromAlignment();
        System.out.println("from alignment UMS length: " + from_alignment_input.length());
        System.out.println("input offset: " + input_offset);
        System.out.println("from alignment input string:");
        System.out.println(from_alignment_input);
        System.out.println();
        System.out.println("Germline gene: " + germlineName);
        String from_alignment_germline = getGermlineStrFromAlignment();
        System.out.println("from alignment germline length: " + from_alignment_germline.length());
        System.out.println("Germline offset: " + germline_offset);
        System.out.println("from alignment germline:");
        System.out.println(from_alignment_germline);
        System.out.println();
        System.out.println("similarity: " + similarity);
        System.out.println("Gaps in alignment: " + gaps_in_alignment);
        //add a line to output the score
        System.out.println("Alignment Score: " + align_score);
        System.out.println();
        System.out.println();
    }
}
