package iHMMuneAlign;

import java.io.*;
import java.util.Vector;
import java.util.ArrayList;
//import javax.swing.JTextArea;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.dp.*;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;
import org.biojavax.bio.seq.RichSequence;

// Referenced classes of package iHMMuneAlign:
//            GlobalDefines, HTMLForEXCELWriter, HTMLforWWW_Writer, BestVGeneFinder, 
//            PostAlignmentResult, A_Score, TrailingJGeneVectorFinder, VpnpDpnpJCnoC, 
//            StateInfo, PostAlignProcessing, GeneInfo, Codon, 
//            JGeneReadingFrame, HTML_Writer, RegionInfo, ProbabilityHolder, 
//            MutabilityScore

class AlignmentThread extends Thread
    implements GlobalDefines
{

    RichSequence sequence;
    String txtFileName;
    String htmlFileName;
    ProbabilityHolder probHolder;
    MutabilityScore mutabilityScore;
    File IGHV_FILE;
    File IGHD_FILE;
    File IGHJ_FILE;
    File IGKV_FILE;
    File IGKJ_FILE;
    //addition files for lambda chain
	File IGLV_FILE;
    File IGLJ_FILE;
    int specifiedChainType;
    byte alignmentType;
    int dGeneAcceptanceType;
    HTML_Writer htmlWriter;
    String scoringMatrix;
    int minJLength;
    int alignmentModelType;
    int numOfVgenes;
    int numSuboptimal;
    int Gap_Behaviour;
    int GAPS_ALLOWED = 1;
    int GAPS_ERROR = 2;
    int GAPS_IGNORE_SEQ = 3;
    String blastLocation;
    boolean DEBUGGING = false;
    String mutationSpecFileName;
    String dProbsFilename;

    public AlignmentThread(RichSequence sequence, File IGHV_FILE, File IGHD_FILE, File IGHJ_FILE, File IGKV_FILE, File IGKJ_FILE, int userSpecifiedChainType, String htmlFileName, ProbabilityHolder probHolder,
            MutabilityScore mutabilityScore, char html_writer_type, byte alignmentType, int dGeneAcceptanceType, String scoringMatrix, int minJLength, int alignmentModelType, int numOfVgenes, int numSuboptimal, String blastLocation,
            File IGLV_FILE, File IGLJ_FILE, String mutationSpecFileName, String dProbsFilename)
    {
        htmlWriter = null;
        this.sequence = sequence;
        //this.out = out;
        this.htmlFileName = htmlFileName;
        this.probHolder = probHolder;
        this.mutabilityScore = mutabilityScore;
        this.IGHV_FILE = IGHV_FILE;
        this.IGHD_FILE = IGHD_FILE;
        this.IGHJ_FILE = IGHJ_FILE;
        this.IGKV_FILE = IGKV_FILE;
        this.IGKJ_FILE = IGKJ_FILE;
        //additional files for the IGL repertoires
        this.IGLV_FILE = IGLV_FILE;
        this.IGLJ_FILE = IGLJ_FILE;
        this.specifiedChainType = userSpecifiedChainType;
        this.alignmentType = alignmentType;
        this.dGeneAcceptanceType = dGeneAcceptanceType;
        this.scoringMatrix = scoringMatrix;
        this.minJLength = minJLength;
        this.alignmentModelType = alignmentModelType;
        this.numOfVgenes = numOfVgenes;
        this.numSuboptimal = numSuboptimal;
	this.blastLocation = blastLocation;
        this.mutationSpecFileName = mutationSpecFileName;
        this.dProbsFilename = dProbsFilename;
        htmlWriter = getHtmlWriter(html_writer_type);
        
        Gap_Behaviour = GAPS_ALLOWED;
    }

    private HTML_Writer getHtmlWriter(char html_writer_type)
    {
        switch(html_writer_type)
        {
        case 1: // '\001'
            return new HTMLForEXCELWriter();

        case 2: // '\002'
            return new HTMLforWWW_Writer();

        case 3: // '\003'
            return new HTMLforCSV_Writer();
        }

        throw new Error("No HTML Writer Type was selected");
    }

    private Object startAlignment()
    {
        //String sequenceString = removeFastaStyleWhiteSpace(sequence.seqString());
        String sequenceName = sequence.getName();
        String errorType = "";
        //System.out.println("about to search for the v genes");
		BestVGeneFinder IGHVfinder = new BestVGeneFinder(IGHV_FILE, blastLocation);
        BestVGeneFinder IGKVfinder = new BestVGeneFinder(IGKV_FILE, blastLocation);
        //add another V gene finder for the lambda
        BestVGeneFinder IGLVfinder = new BestVGeneFinder(IGLV_FILE, blastLocation);
        int chainType = 1; // 1 for heavy, 2 for kappa, 3 for lambda <- this is over-written by user specified chain type below
        int MIN_UMS_ALIGNMENT_END_OFFSET = 12;
        int MIN_IGHJ_ALIGNMENT_START_OFFSET = minJLength;
        //set the MaxScore for the identification of the first Vgene to the largest possible number
        float prevMaxScore = java.lang.Float.MAX_VALUE;
        int reverseComplement = 0; // set to 1 if input is reverse complement of original gene sequence
        PostAlignmentResult V_align_result = null;
        PostAlignmentResult IGHV_align_result = null;
        PostAlignmentResult IGKV_align_result = null;
        //for lamdba
        PostAlignmentResult IGLV_align_result = null;
	
		//create arraylists from the IGHD and IGHJ germline repertoire files
		FastaReader fr = new FastaReader();
		ArrayList<RichSequence> dSequences = fr.readFile(IGHD_FILE);
		ArrayList<RichSequence> jhSequences = fr.readFile(IGHJ_FILE);
		ArrayList<RichSequence> jkSequences = fr.readFile(IGKJ_FILE);
		//IGLJ
		ArrayList<RichSequence> jlSequences = fr.readFile(IGLJ_FILE);
		jlSequences = fr.readFile(IGLJ_FILE);
		
        //System.out.println("input seq & rc seq:");
        //System.out.println(sequenceString+"\n"+rcSeqString+"\n");
        //System.exit(0);
        
        //find multiple V genes
        int vIndex;
        for (vIndex=0; vIndex<numOfVgenes; vIndex++)
        {
            // add an empty line to indicate new seq or new v-gene
        	PrintWriter outStream = null;
        	try
            {
                outStream = new PrintWriter(new FileWriter(htmlFileName, true));
            }
            catch(Exception ex)
            {
                System.out.println("couldn't open new file: " + ex.getMessage());
            }
            outStream.println();
            //outStream.println("Alignments for "+sequenceName+" containing the V-gene "+VGene_name+":");
            outStream.close();
            
            /*find the prev scoring alignment if we are getting any subsequent v genes
            //if (vIndex > 0) {
               // if previous alignment successful then get the best score for the last alignment
            //   prevMaxScore = V_align_result.getAlignmentScore();
            }*/

            // do the V gene alignment with IGHV & IGKV repertoires and determine whether the input seq is heavy or light
            if (specifiedChainType == 1) {// heavy chain
                V_align_result = IGHVfinder.getResult(sequence, vIndex);
                chainType = 1;
            } else if (specifiedChainType == 2) {// kappa chain
            	V_align_result = IGKVfinder.getResult(sequence, vIndex);
            	chainType = 2;
            } else if (specifiedChainType == 3) { //lambda chain
            	V_align_result = IGLVfinder.getResult(sequence, vIndex);
             	chainType = 3;
            } else if (specifiedChainType == 0) {// determine automatically
            	//find best matching V genes for each of the repertoires
            	IGHV_align_result = IGHVfinder.getResult(sequence, vIndex);
            	IGKV_align_result = IGKVfinder.getResult(sequence, vIndex);
            	IGLV_align_result = IGLVfinder.getResult(sequence, vIndex);
            	
            	//determine which repertoire provided the highest alignment score
                if ((IGKV_align_result.getAlignmentScore() > IGHV_align_result.getAlignmentScore()) && (IGKV_align_result.getAlignmentScore() > IGLV_align_result.getAlignmentScore()) ){
            	    //kappa repertoire was the source of the highest alignment result
            	    chainType = 2;
                    V_align_result = IGKV_align_result;
                } else if ((IGHV_align_result.getAlignmentScore() > IGKV_align_result.getAlignmentScore()) && (IGHV_align_result.getAlignmentScore() > IGLV_align_result.getAlignmentScore())) {
                	//heavy chain repertoire was the source of the highest alignment result
            	    chainType = 1;
            	    V_align_result = IGHV_align_result;
                } else {
                	//assume that it must be the lambda that achieved the highest alignment score
                	chainType = 3;
                	V_align_result = IGLV_align_result;
                }
            }
            
            // try the v-alignment with the reverse complement seq
            //temp comment out for the 454 run
            /*
            if (vIndex == 0) {
				// create the reverse complement seq of the input seq
				String rcSeqString = getReverseComp(sequence.seqString());
				RichSequence rcSequence = null;
				
				try {
					rcSequence = RichSequence.Tools.createRichSequence((sequence.getName() + "(reverse_comp)"), rcSeqString, DNATools.getDNA());
				} catch (BioException bioe) {
					throw new Error("(AlignmentThread) Could not create rich sequence from reverse complement sequence: " + bioe.getMessage());
				}
				
            	int rcChainType = 1;
            	PostAlignmentResult V_align_result_rc = null;
            	PostAlignmentResult IGHV_align_result_rc = null;
            	PostAlignmentResult IGKV_align_result_rc = null;
            	if (specifiedChainType == 1) {// heavy chain
	                V_align_result_rc = IGHVfinder.getResult(rcSequence, vIndex);
	                rcChainType = 1;
	            } else if (specifiedChainType == 2) {// kappa chain
	            	V_align_result_rc = IGKVfinder.getResult(rcSequence, vIndex);
	                rcChainType = 2;
	            } else if (specifiedChainType == 0) {// determine automatically
	            	IGHV_align_result_rc = IGHVfinder.getResult(rcSequence, vIndex);
	            	IGKV_align_result_rc = IGKVfinder.getResult(rcSequence, vIndex);
	                if (IGHV_align_result_rc.getAlignmentScore() < IGKV_align_result_rc.getAlignmentScore()) {
	            	    rcChainType = 2;
	                    V_align_result_rc = IGKV_align_result_rc;
	                } else {
	            	    rcChainType = 1;
	            	    V_align_result_rc = IGHV_align_result_rc;
	                }
	            }
            	if (V_align_result_rc.getAlignmentScore() > V_align_result.getAlignmentScore()) {
            		reverseComplement = 1;
            		//sequenceString = rcSeqString;
            		V_align_result = V_align_result_rc;
            		chainType = rcChainType;
            	}
            }*/
            
            //get gaps in the best alignment
            int alignment_gap_no = 0;
            alignment_gap_no = V_align_result.getAlignmentGaps();

            //if there are gaps print msg and skip the hmm model creation as gaps not yet implemented
            if((alignment_gap_no > 0) && (Gap_Behaviour == GAPS_ERROR))
            {
                 //throw new Error("VGene finder result equal NULL");

                 //do something more meaningful then just throwing an error
                 //System.out.println("There were gaps in the IGHV gene. iHMMune-align cannot proceed with HMM creation.");

                 //get the IGHV name
                 String VGene_name = V_align_result.getVGeneName();
                 //get the input seq that matches to the IGHV region
                 char[] alignedInput = V_align_result.getAlignedInput();
                 //convert the alignedInput to a string
                 String vOutput = new String(alignedInput);
                 //convert the output to lower case -- note: not showing any match/mismatch details just where the gaps are
                 vOutput = vOutput.toLowerCase();

                 //set the error to return
                 //errorType = "NA - Gaps in IGHV gene";
                 errorType = "NA - Gaps in V gene " + VGene_name + " with " + alignment_gap_no + " gaps";

                 //call the html writer to create the output file - or to append to the output file if already exists and alignment type =4
                 htmlWriter.createErrDocument(htmlFileName ,sequenceName, alignmentType, errorType, vOutput);

                 // close the file once we are done writing the output
                 if(htmlWriter != null)
                 {
                    htmlWriter.closeHTMLFile();
                 }
        
                 //return
                 //return null;   // this will exit the loop and skip subsequent IGHV genes (if there are any remaining to be sought)

            } else if (alignment_gap_no < 0) { //this would flag that no IGHV result was found from the repertoire search
				//throw new Error("VGene finder result equal NULL");
				
				//do something more meaningful then just throwing an error
				//System.out.println("A V gene could not be identified in the input sequence. iHMMune-align cannot proceed with HMM creation.");
				
				//set the error to return
				//errorType = "NA - Gaps in IGHV gene";
				errorType = "NA - No V found in input sequence";
				String vOutput = "NA - not found";
				
				//call the html writer to create the output file - or to append to the output file if already exists and alignment type =4
				htmlWriter.createErrDocument(htmlFileName ,sequenceName, alignmentType, errorType, vOutput);
				
				// close the file once we are done writing the output
				if(htmlWriter != null)
				{
                    htmlWriter.closeHTMLFile();
				}
				
			} else { // no gaps and V gene found, or gaps allowed and v gene found, so ready to try and build the HMM
                String from_alignment_start_UMS_string = V_align_result.getUMSfromAlignmentString();
                String from_alignment_start_VGene_string = V_align_result.getVGeneAlignmentString();
                String UMS_name = V_align_result.getUMSname();
                String VGene_name = V_align_result.getVGeneName();
                RichSequence from_alignment_start_UMS_seq = null;
                RichSequence from_alignment_start_VGene_seq = null;
                try
                {
                    from_alignment_start_UMS_seq = RichSequence.Tools.createRichSequence(UMS_name, from_alignment_start_UMS_string, DNATools.getDNA());
                    from_alignment_start_VGene_seq = RichSequence.Tools.createRichSequence(VGene_name, from_alignment_start_VGene_string, DNATools.getDNA());
                } 
                catch(IllegalSymbolException illse)
                {
                    throw new Error(illse.getMessage());
                } 
				catch (BioException bioe) {
					throw new Error("(AlignmentThread:startAlignment) Error creating RichSequence: " + bioe.getMessage());
				}

                int VGene_start_offset = V_align_result.getColOffset();
                int completeVGeneLength = VGene_start_offset + from_alignment_start_VGene_string.length();
                
                if (DEBUGGING) {
                	System.out.println("(AlignmentThread:startAlignment) input string after V alignment: " + from_alignment_start_UMS_string);
                	System.out.println("(AlignmentThread:startAlignment) germline string after V alignment: " + from_alignment_start_VGene_string);
                	System.out.println("(AlignmentThread:startAlignment) input name: " + UMS_name);
                	System.out.println("(AlignmentThread:startAlignment) germline gene name from V alignment: " + VGene_name);
                	System.out.println("(AlignmentThread:startAlignment) input RichSequence: " + from_alignment_start_UMS_seq.seqString());
                	System.out.println("(AlignmentThread:startAlignment) germline RichSequence: " + from_alignment_start_VGene_seq.seqString());
                	System.out.println("(AlignmentThread:startAlignment) germline start offset: " + VGene_start_offset);
                	System.out.println("(AlignmentThread:startAlignment) germline complete length: " + completeVGeneLength);
                }
                
                A_Score myA_Score = new A_Score();
                double A_probability = 0;
                if (VGene_start_offset < 128) { //128 is the nucleotide before the start of the defined common region
                     //the Voffset indicates that the alignment includes the start of the common region so use the "long" Vgene prob for A_Score		
                     A_probability = myA_Score.A_probabilityLong(from_alignment_start_VGene_string, from_alignment_start_UMS_string, VGene_start_offset, null);
                } else  {
                	// the offset is indicating that the V-REGION starts downstream of the common region defined for the A_score calcs
                	// this is likely the result of the use of FR2 primer sets, use the "short" V-REGION probability for A_score calcs
                	 A_probability = myA_Score.A_probabilityShort(from_alignment_start_VGene_string, from_alignment_start_UMS_string, VGene_start_offset, null);
                }
                
                //System.out.println("A_probability = " + A_probability);
                //if the V was too short for A_prob to be calculated then just output
                if (A_probability > 1) {
                   //call the functions to output results
                   if (A_probability == 2) {
                      //System.out.println("Invalid A_probability due to short IGHV gene");
                   }
                   if (A_probability == 3) {
                      //System.out.println("Invalid A_probability due to too great truncation of the 5` end of the V");
                   }
                   if (A_probability == 4) {
                      //System.out.println("Invalid A_probability due to differing lengths of the V and corresponding input sequence");
                   }
        
                   errorType = "NA - V too short";
                   String vOutput = "NA";
        
                   //call the html writer to create the output file - or to append to the output file if already exists and alignment type =4
                   htmlWriter.createErrDocument(htmlFileName ,sequenceName, alignmentType, errorType, vOutput);
        
                   // close the file once we are done writing the output
                   if(htmlWriter != null)
                   {
                      htmlWriter.closeHTMLFile();
                   }
        
                   //return
                   return null;

                } else { //we have an A_prob so continue with the model construction
					//hacky way to prevent building of full sub-optimal tree when we are not interested in retriving anything other than a single match
						if (numSuboptimal == 0) {
							//determine if there is any sequence downstream of the IGHJ gene
							TrailingJGeneVectorFinder vectorFinder = null;
							if (chainType == 1) {  //heavy chain
							    vectorFinder = new TrailingJGeneVectorFinder(IGHJ_FILE, blastLocation);
							} else if (chainType == 2) {  //kappa chain
							    vectorFinder = new TrailingJGeneVectorFinder(IGKJ_FILE, blastLocation);
							} else if (chainType == 3) { //lambda chain
								vectorFinder = new TrailingJGeneVectorFinder(IGLJ_FILE, blastLocation);
							}
							PostAlignmentResult JAlignResult = vectorFinder.getAlignmentResult(from_alignment_start_UMS_string);
							String UMS_no_C = vectorFinder.getResult(JAlignResult, MIN_IGHJ_ALIGNMENT_START_OFFSET);
							
							if (DEBUGGING) {
								System.out.println("(AlignmentThread) input string without trailing J vector:" + UMS_no_C);
							}
							
							//need to check that the TrailingJGeneVector finder actually found a suitable IGHJ match
							if (UMS_no_C == "NA")  {
								//the min length for the J gene was not satisfied
								
								//want to stop further processing and just give some meaningful output
								if (DEBUGGING) {
									System.out.println("J gene could not be found of a suitable length. Min length specified as " + MIN_IGHJ_ALIGNMENT_START_OFFSET);
								}
								
								//get the IGHV name
								VGene_name = V_align_result.getVGeneName();
								//set the error type
								errorType = "NA - J too short, possible V alignment to " + VGene_name;
								if (DEBUGGING) {
									System.out.println("(AlignmentThread) Throwing error: " + errorType);
								}
								//String vOutput = "NA";
								//amended to still output an IGHV
								//get the input seq that matches to the IGHV region
								char[] alignedInput = V_align_result.getAlignedInput();
								//convert the alignedInput to a string
								String vOutput = new String(alignedInput);
								vOutput = vOutput.toLowerCase();
								
								
								//call the html writer to create the output file - or to append to the output file if already exists and alignment type =4
								htmlWriter.createErrDocument(htmlFileName ,sequenceName, alignmentType, errorType, vOutput);
								
								// close the file once we are done writing the output
								if(htmlWriter != null)
								{
									htmlWriter.closeHTMLFile();
								}
								
								//return
								return null;
							} else if (UMS_no_C == "GAP") {
								//this occurs when the GAP_BEHAVIOUR for the IGHJ ID step of the trailingJVector is set to GAP_ERROR
								//and the best BLAST result contains a gap
								//get the IGHV name
								VGene_name = V_align_result.getVGeneName();
								//set the error type
								errorType = "NA - IGHJ contains gaps, possible V alignment to " + VGene_name;
								if (DEBUGGING) {
									System.out.println("(AlignmentThread) Throwing error: " + errorType);
								}
								//String vOutput = "NA";
								//amended to still output an IGHV
								//get the input seq that matches to the IGHV region
								char[] alignedInput = V_align_result.getAlignedInput();
								//convert the alignedInput to a string
								String vOutput = new String(alignedInput);
								vOutput = vOutput.toLowerCase();
								
								
								//call the html writer to create the output file - or to append to the output file if already exists and alignment type =4
								htmlWriter.createErrDocument(htmlFileName ,sequenceName, alignmentType, errorType, vOutput);
								
								// close the file once we are done writing the output
								if(htmlWriter != null)
								{
									htmlWriter.closeHTMLFile();
								}
								
								//return
								return null;
								
							} //end throwing error for gap in J when gaps not allowed or J found being shorter than J min
							
							
							//if a suitable IGHJ match was found
							//set J sequences based on chain type
							ArrayList<RichSequence> jSequences = null;
							if (chainType == 1) {  //heavy chain
							    jSequences = jhSequences;
							} else if (chainType == 2) {  //kappa chain
							    jSequences = jkSequences;
							    dSequences = null;
							} else if (chainType == 3) {
								jSequences = jlSequences;
								dSequences = null;
							}
							
							//check if there are gaps in the J, if this has not caused an error to be thrown above need to add a new J sequence to the repertoire
							//that incorporates the gaps so the the HMM can build using the sequence
							if (JAlignResult.getAlignmentGaps() > 0) {
								//gaps in the best J found and no error so add a sequence
								String germlineJPlusIndels = JAlignResult.getJGeneAlignmentString();
								String indelJName = JAlignResult.getVGeneName();
								
								//using the string and the name create a rich sequence
								RichSequence indelJGene = null;
								try {
									indelJGene = RichSequence.Tools.createRichSequence(indelJName, germlineJPlusIndels, DNATools.getDNA());
								} catch (BioException bioe) {
									throw new Error("(AlignmentThread:startAlignment) There was a problem created the rich sequence " + indelJName + "\n" + bioe.getMessage());
								}
								
								if (DEBUGGING) {
									System.out.println("(AlignmentThread:startAlignment) Created a new IGHJ entry: " + indelJName + " " + germlineJPlusIndels);
								}
								
								//add the rich sequence to the Jsequence arraylist
								jSequences.add(indelJGene);
							}
							
							String from_alignment_start_UMS_no_C_string = null;
							RichSequence from_alignment_start_UMS_no_C_seq = null;
							//markov model //VpnpDpnpJCnoC vdj = new VpnpDpnpJCnoC(mutabilityScore);
							VnDnJCnoC vdj = new VnDnJCnoC(mutabilityScore, mutationSpecFileName);
							org.biojava.bio.dp.MarkovModel markov_model;
															
							if(UMS_no_C != null) {
								
								if (UMS_no_C == "NA") {
									//we have no match to the IGHJ seq, but we are allowing partial seqs
									markov_model = vdj.createModel(from_alignment_start_VGene_seq, VGene_start_offset, completeVGeneLength, from_alignment_start_UMS_no_C_string.length(), dSequences, jSequences, probHolder, false, A_probability, alignmentModelType, chainType);
									
								} else {
									from_alignment_start_UMS_no_C_string = UMS_no_C;
									try {
										//attempt to fix state path issue with gap char in input sequence
										String inputSeqNoGapChars = UMS_no_C.replace("-","n");
										from_alignment_start_UMS_no_C_seq = RichSequence.Tools.createRichSequence("UMS_name", inputSeqNoGapChars, DNATools.getDNA());
										//from_alignment_start_UMS_no_C_seq = RichSequence.Tools.createRichSequence("UMS_name", UMS_no_C, DNATools.getDNA());
										
										if (DEBUGGING) {
											System.out.println("(AlignmentThread:startAlignment) RichSequence created from: " + from_alignment_start_UMS_no_C_seq.seqString());
										}
									} catch(IllegalSymbolException illse) {
										throw new Error(illse.getMessage());
									} catch (BioException bioe) {
										throw new Error("(AlignmentThread:startAlignment) There was a problem creating RichSequence for " + UMS_no_C + "\n" + bioe.getMessage());
									} 
									
									if (DEBUGGING) {
										System.out.println("(AlignmentThread:startAlignment) About to create model.");
									}
																	
									markov_model = vdj.createModel(from_alignment_start_VGene_seq, VGene_start_offset, completeVGeneLength, from_alignment_start_UMS_no_C_string.length(), dSequences, jSequences, probHolder, true, A_probability, alignmentModelType, chainType);

									//vdj = null;
									
									if (DEBUGGING) {
										System.out.println("(AlignmentThread:startAlignment) Model creation complete.");
									}
								}
							} else {
								//means that the IGHJ length found is equal to or less than full length germline gene
								//possible that there is truncation of the IGHJ
								//System.out.println("NO C_Region Removal");
								from_alignment_start_UMS_no_C_string = from_alignment_start_UMS_string;
								from_alignment_start_UMS_no_C_seq = from_alignment_start_UMS_seq;
								markov_model = vdj.createModel(from_alignment_start_VGene_seq, VGene_start_offset, completeVGeneLength, from_alignment_start_UMS_no_C_string.length(), dSequences, jSequences, probHolder, false, A_probability, alignmentModelType, chainType);
								//vdj = null;
							}

//                            Utils.printModel(markov_model);

// ================================================================================
                            try {
                                FileOutputStream file = new FileOutputStream(new File("model.hmm"));
                                ModelWriter mr = new ModelWriter(markov_model, file);
                                mr.write();
                            } catch (FileNotFoundException e) {
                                e.printStackTrace();
                            }

// ================================================================================

							DP dp = null;
							try {
								
								//create the markov model
								if (DEBUGGING) {
									System.out.println("about to create DP");
								}
							
							
								dp = DPFactory.DEFAULT.createDP(markov_model);
								
								if (DEBUGGING) {
									System.out.println("Dp created");
								}
								
							} catch(BioException bioe) {
								throw new Error(bioe.getMessage());
							}

							
							SymbolList res_array[] = {from_alignment_start_UMS_no_C_seq};
							if(res_array == null) {
								if (DEBUGGING) {
									System.out.println("res_array == null : could not create nt array from input sequence");
								}
								
								System.exit(0);
							} else {
								if (DEBUGGING) {
									System.out.println("res_array != null : created nt array from input sequence");
								}
							}
							
							StatePath v = null;
							try {
								v = dp.viterbi(res_array, ScoreType.PROBABILITY);
								//forwardScore = dp.forward(res_array, ScoreType.PROBABILITY);
							} catch(IllegalAlphabetException illalphe) {
							
								throw new Error(illalphe.getMessage());
								
							} catch(IllegalSymbolException illsymbole) {
							
								throw new Error(illsymbole.getMessage());
								
							} catch(IllegalTransitionException illsymbole) {
							
								throw new Error(illsymbole.getMessage());
								
							} catch(BioException bioex) { 
							
								throw new Error(bioex.getMessage());
								
							} catch(NullPointerException nulle) {
								if (DEBUGGING) {
								
									System.out.println("(AlignmentThread): Null pointer expection in creation of state path.");
									nulle.printStackTrace();
									
								}
								//System.exit(0);
								
								//create an error and return null instead of exiting as there may still be other seqs to partition as part of this batch
								errorType = "State path could not be created";
						        String vOutput = from_alignment_start_UMS_no_C_seq.seqString();
				
						        //call the html writer to create the output file - or to append to the output file if already exists and alignment type =4
						        htmlWriter.createErrDocument(htmlFileName ,sequenceName, alignmentType, errorType, vOutput);
				
						        // close the file once we are done writing the output
						        if(htmlWriter != null)
						        {
						           htmlWriter.closeHTMLFile();
						        }
				
						        //return
						        return null;
								
								
							}
							
							//get state path
							SymbolList viterbi_state_seq = v.symbolListForLabel(StatePath.STATES);
							//get path score
							double pathScore = v.getScore();
							
							displayStateInfo(viterbi_state_seq, from_alignment_start_UMS_no_C_string, sequenceName);
							displayEmissionInfo(viterbi_state_seq, from_alignment_start_UMS_no_C_string, sequenceName, pathScore, htmlFileName, chainType, reverseComplement);
							// output A_score and related data
							myA_Score.dumpInfo(htmlFileName);
							
							//done with model so reset state path, dp & markov model objects
							v = null;
							dp = null;
							//markov_model = vdj.destroyModel();
						
						} else { //build sub-optimal alignments
								
							// the list of n best alternative alignments for this VGene identified
							PathNode[] bestPaths = new PathNode[numSuboptimal];
							// fill the list with empty nodes first
							for (int i = 0; i < numSuboptimal; i++) {
								bestPaths[i] = new PathNode();
							}
							
							// nodes on the curr level & children level of the binary tree
							ArrayList currLevelNodes = new ArrayList(512); //for a maximum of 10 sub-optimal, 2^9=512
							ArrayList nextLevelNodes = new ArrayList(512);
							
							// the root of the tree is the optimal alignment (no removal of D/J repertoire sequences)
							if (chainType == 1) { // heavy chain
								currLevelNodes.add(new PathNode("","", IGHD_FILE, IGHJ_FILE));
							} else if (chainType == 2) { // kappa light chain
								currLevelNodes.add(new PathNode("","", null, IGKJ_FILE));
							} else if (chainType == 3) { // lambda light chain
								currLevelNodes.add(new PathNode("","", null, IGLJ_FILE));
							}
							
							//System.out.println("Size of currLevelNodes: (before while loop) " + currLevelNodes.size());
							
							// loop through each level of the tree until there's no more potential candidates
							while (currLevelNodes.size() != 0) {
								//System.out.println("Size of currLevelNodes: (after while loop) " + currLevelNodes.size());
								
								nextLevelNodes.clear();
								for (int i = 0; i < currLevelNodes.size(); i++) {
									PathNode path = (PathNode)currLevelNodes.get(i);
									File D_FILE = null;
									if (chainType == 1) { // heavy chain
										D_FILE = path.getDFile();
										//output the D filename
										if (DEBUGGING) {
											System.out.println("IGHD from PathNode: " + i + " " + D_FILE.getName());
										}
									}
									File J_FILE = path.getJFile();
									//ouptut the J filename
									if (DEBUGGING) {
										System.out.println("IGHJ from PathNode: " + i + " " + J_FILE.getName());
									}
									
									//determine if there is any sequence downstream of the IGHJ gene
									//TrailingJGeneVectorFinder vectorFinder = new TrailingJGeneVectorFinder(); //changed to actually pass the input j repertoire
									TrailingJGeneVectorFinder vectorFinder = new TrailingJGeneVectorFinder(J_FILE, blastLocation);
									PostAlignmentResult JAlignResult = vectorFinder.getAlignmentResult(from_alignment_start_UMS_string);
									String UMS_no_C = vectorFinder.getResult(JAlignResult, MIN_IGHJ_ALIGNMENT_START_OFFSET);
									//String UMS_no_C = vectorFinder.getResult(from_alignment_start_UMS_string, MIN_IGHJ_ALIGNMENT_START_OFFSET);
					  
									//need to check that the TrailingJGeneVector finder actually found a suitable IGHJ match
									if (UMS_no_C == "NA")  {
										//the min length for the J gene was not satisfied
							
										//want to stop further processing and just give some meaningful output
										//System.out.println("J gene could not be found of a suitable length. Min length specified as " + MIN_IGHJ_ALIGNMENT_START_OFFSET);
							
										//get the IGHV name
										VGene_name = V_align_result.getVGeneName();
										//set the error type
										errorType = "NA - J too short, possible IGHV alignment to " + VGene_name;
										//String vOutput = "NA";
										//amended to still output an IGHV
										//get the input seq that matches to the IGHV region
										char[] alignedInput = V_align_result.getAlignedInput();
										//convert the alignedInput to a string
										String vOutput = new String(alignedInput);
										vOutput = vOutput.toLowerCase();
										

										//call the html writer to create the output file - or to append to the output file if already exists and alignment type =4
										htmlWriter.createErrDocument(htmlFileName ,sequenceName, alignmentType, errorType, vOutput);
							
										// close the file once we are done writing the output
										if(htmlWriter != null)
										{
										   htmlWriter.closeHTMLFile();
										}
							
										//return
										return null;
									} else if (UMS_no_C == "GAP") {
										//this occurs when the GAP_BEHAVIOUR for the IGHJ ID step of the trailingJVector is set to GAP_ERROR
										//and the best BLAST result contains a gap
										//get the IGHV name
										VGene_name = V_align_result.getVGeneName();
										//set the error type
										errorType = "NA - IGHJ contains gaps, possible V alignment to " + VGene_name;
										//String vOutput = "NA";
										//amended to still output an IGHV
										//get the input seq that matches to the IGHV region
										char[] alignedInput = V_align_result.getAlignedInput();
										//convert the alignedInput to a string
										String vOutput = new String(alignedInput);
										vOutput = vOutput.toLowerCase();
								
								
										//call the html writer to create the output file - or to append to the output file if already exists and alignment type =4
										htmlWriter.createErrDocument(htmlFileName ,sequenceName, alignmentType, errorType, vOutput);
								
										// close the file once we are done writing the output
										if(htmlWriter != null)
										{
											htmlWriter.closeHTMLFile();
										}
								
										//return
										return null;
								
									} //end throwing error for gap in J when gaps not allowed or J found being shorter than J min                     
								  
									String from_alignment_start_UMS_no_C_string = null;
									RichSequence from_alignment_start_UMS_no_C_seq = null;
									//VpnpDpnpJCnoC vdj = new VpnpDpnpJCnoC(mutabilityScore);
									VnDnJCnoC vdj = new VnDnJCnoC(mutabilityScore, mutationSpecFileName);
									
									org.biojava.bio.dp.MarkovModel markov_model;
									//set J sequences based on chain type
								    ArrayList<RichSequence> jSequences = null;
								    if (chainType == 1) {  //heavy chain
								         jSequences = jhSequences;
								    } else if (chainType == 2) {  //kappa chain
								        jSequences = jkSequences;
								        dSequences = null;
								    } else if (chainType == 3) { //lambda chain
								    	jSequences = jlSequences;
								    	dSequences = null;
								    }
								    
								    //check if there are gaps in the J, if this has not caused an error to be thrown above need to add a new J sequence to the repertoire
									//that incorporates the gaps so the the HMM can build using the sequence
									if (JAlignResult.getAlignmentGaps() > 0) {
										//gaps in the best J found and no error so add a sequence
										String germlineJPlusIndels = JAlignResult.getJGeneAlignmentString();
										String indelJName = JAlignResult.getVGeneName();
								
										//using the string and the name create a rich sequence
										RichSequence indelJGene = null;
										try {
											indelJGene = RichSequence.Tools.createRichSequence(indelJName, germlineJPlusIndels, DNATools.getDNA());
										} catch (BioException bioe) {
											throw new Error("(AlignmentThread:startAlignment) There was a problem created the rich sequence " + indelJName + "\n" + bioe.getMessage());
										}
								
										//add the rich sequence to the Jsequence arraylist
										jSequences.add(indelJGene);
									}
								       								  
									  if(UMS_no_C != null)
									  {

								          
										  if (UMS_no_C == "NA") {
											//we have no match to the J seq, but we are allowing partial seqs
											  markov_model = vdj.createModel(from_alignment_start_VGene_seq, VGene_start_offset, completeVGeneLength, from_alignment_start_UMS_no_C_string.length(), dSequences, jSequences, probHolder, false, A_probability, alignmentModelType, chainType);
							
										  } else {
											from_alignment_start_UMS_no_C_string = UMS_no_C;
											try
											{
												//attempt to fix state path issues
												String inputSeqNoGapChars = UMS_no_C.replace("-","n");
												from_alignment_start_UMS_no_C_seq = RichSequence.Tools.createRichSequence("UMS_name", inputSeqNoGapChars, DNATools.getDNA());
												//from_alignment_start_UMS_no_C_seq = RichSequence.Tools.createRichSequence("UMS_name", UMS_no_C, DNATools.getDNA());
												
												if (DEBUGGING) {
													System.out.println("(AlignmentThread) RichSequence created for input sequence: " + from_alignment_start_UMS_no_C_seq.seqString());
												}
											}
											catch(IllegalSymbolException illse)
											{
												throw new Error(illse.getMessage());
											}
											catch (BioException bioe) {
												throw new Error("There was a problem creating RichSequence for " + UMS_no_C + "\n" + bioe.getMessage());
											}
											//VpnpDpnpJCnoC vdj = new VpnpDpnpJCnoC(mutabilityScore);
											markov_model = vdj.createModel(from_alignment_start_VGene_seq, VGene_start_offset, completeVGeneLength, from_alignment_start_UMS_no_C_string.length(), dSequences, jSequences, probHolder, true, A_probability, alignmentModelType, chainType);
										  }
									  } else
									  {
										  //means that the IGHJ length found is equal to or less than full length germline gene
										  //possible that there is truncation of the IGHJ
										  //System.out.println("NO C_Region Removal");
										  from_alignment_start_UMS_no_C_string = from_alignment_start_UMS_string;
										  from_alignment_start_UMS_no_C_seq = from_alignment_start_UMS_seq;
										  //VpnpDpnpJCnoC vdj = new VpnpDpnpJCnoC(mutabilityScore);
										  markov_model = vdj.createModel(from_alignment_start_VGene_seq, VGene_start_offset, completeVGeneLength, from_alignment_start_UMS_no_C_string.length(), dSequences, jSequences, probHolder, false, A_probability, alignmentModelType, chainType);
										  //vdj = null;
									  }
					
									  //System.out.println("about to create DP");
									  DP dp = null;
									  try
									  {
										  dp = DPFactory.DEFAULT.createDP(markov_model);
									  }
									  catch(BioException bioe)
									  {
										  throw new Error(bioe.getMessage());
									  }
									 //System.out.println("Dp created");
						
									  SymbolList res_array[] = {from_alignment_start_UMS_no_C_seq};
									  if(res_array == null)
									  {
										  System.out.println("res_array == null: could not create nt array from input sequence.");
										  System.exit(0);
									  }
									  
									  StatePath v = null;
									  try
									  {
										  v = dp.viterbi(res_array, ScoreType.PROBABILITY);
										  //forwardScore = dp.forward(res_array, ScoreType.PROBABILITY);
									  }
									  catch(IllegalAlphabetException illalphe)
									  {
										  throw new Error(illalphe.getMessage());
									  }
									  catch(IllegalSymbolException illsymbole)
									  {
										  throw new Error(illsymbole.getMessage());
									  }
									  catch(IllegalTransitionException illsymbole)
									  {
										  throw new Error(illsymbole.getMessage());
									  }
									  catch(BioException bioex) {
										  throw new Error(bioex.getMessage());
									  }
									  catch(NullPointerException nulle) {
										 	if (DEBUGGING) {
												System.out.println("(AlignmentThread): Null pointer expection in creation of state path.");
												nulle.printStackTrace();
											}
										   //System.exit(0);
								
										   //create an error and return null instead of exiting as there may still be other seqs to partition as part of this batch
											errorType = "NA - State path could not be created";
											String vOutput = from_alignment_start_UMS_no_C_seq.seqString();
				
											//call the html writer to create the output file - or to append to the output file if already exists and alignment type =4
											htmlWriter.createErrDocument(htmlFileName ,sequenceName, alignmentType, errorType, vOutput);
				
											// close the file once we are done writing the output
											if(htmlWriter != null)
											{
											   htmlWriter.closeHTMLFile();
											}
				
											//return
											return null;
									  }
									  //catch(Exception ex)
									  //{
									  //  throw new Error(ex.getMessage());
									  //}
									
									double pathScore = v.getScore(); 
									
									//if curr path has a higher score than the min of the n-best, replace the path which has lowest score in n-best
									if (pathScore > bestPaths[numSuboptimal-1].getScore()) {
										SymbolList viterbi_state_seq = v.symbolListForLabel(StatePath.STATES);
										// print to nowhere, we are just using the codes in displayGene() to get the gene names
										PrintWriter stream = null;
										GeneInfo DGeneInfo = null;
										if (chainType == 1) {
											DGeneInfo = displayGene(from_alignment_start_UMS_no_C_string, viterbi_state_seq, 'D', stream);
										}
										GeneInfo JGeneInfo = displayGene(from_alignment_start_UMS_no_C_string, viterbi_state_seq, 'J', stream);
										String DGeneName = "";
										if (chainType == 1) {
											DGeneName = DGeneInfo.gene_name;
										}
										String JGeneName = JGeneInfo.gene_name;
										
										if (DEBUGGING) {
											System.out.println("curr gene info: "+DGeneName+" "+JGeneName+" "+pathScore); //System.exit(0);
										}
											
										path.setDName(DGeneName);
										path.setJName(JGeneName);
										path.setScore(pathScore);
										path.setStateSeq(viterbi_state_seq);
										path.setUMSNoC(from_alignment_start_UMS_no_C_string);
										bestPaths[numSuboptimal-1] = path;
										for (int k=0; k < bestPaths.length; k++) {
											if (DEBUGGING) {
												System.out.println("best[k]:"+bestPaths[k].getDName()+bestPaths[k].getJName());
											}
										}
										
										// bubble the path that has the lowest score to the end of the bestPaths list
										for (int j = 0; j < numSuboptimal-1; j++) {
											if (bestPaths[j].getScore() < bestPaths[j+1].getScore()) {
												PathNode tmpPath = bestPaths[j];
												bestPaths[j] = bestPaths[j+1];
												bestPaths[j+1] = tmpPath;
											}
										}
										  
										// add potential children node candidates
										if (chainType == 1) { // heavy chain
											nextLevelNodes.add(new PathNode(DGeneName, "", D_FILE, J_FILE));//remove D
											nextLevelNodes.add(new PathNode("", JGeneName, D_FILE, J_FILE));//remove J
										} else { // light chain
											nextLevelNodes.add(new PathNode("", JGeneName, null, J_FILE));//remove top J
										}
									}
									
									//done with model so reset state path, dp & markov model objects
									v = null;
									dp = null;
									//markov_model = vdj.destroyModel();
									vdj = null;
								}
								currLevelNodes.clear();
								//System.out.println("next level:"+nextLevelNodes.size());
								
								currLevelNodes.addAll(nextLevelNodes);
								
							}
							
							//output
							//create a file to store the output in --> this is separate to the selected HTML/CSV/WWW type of output
							//PrintWriter fostream = null;
							//String outputStreamFile = htmlFileName + "fostream";
							//fostream = new PrintWriter (new FileWriter (outputStreamFile, true));
							
							// sort bestPaths, output in order from highest to lowest score
							double maxScore;
							int maxIndex;
							PathNode tmpPath;
							for (int i = numSuboptimal; i > 0; i--) {
								// find the highest scoring alignment
								maxScore = bestPaths[0].getScore();
								maxIndex = 0;
								for (int j = 1; j < i; j++) {
									if (bestPaths[j].getScore() > maxScore) {
										maxScore = bestPaths[j].getScore();
										maxIndex = j;
									}
								}
								// swap with the last active alignment
								if (i > 1) {
									tmpPath = bestPaths[maxIndex];
									bestPaths[maxIndex] = bestPaths[i-1];
									bestPaths[i-1] = tmpPath;
								}
								// output information on this alignment
								displayStateInfo(bestPaths[i-1].getStateSeq(), bestPaths[i-1].getUMSNoC(), sequenceName);
								displayEmissionInfo(bestPaths[i-1].getStateSeq(), bestPaths[i-1].getUMSNoC(), sequenceName, bestPaths[i-1].getScore(), htmlFileName, chainType, reverseComplement);
								// output A_score and related data
								myA_Score.dumpInfo(htmlFileName);
								
							} 
						}//end else for numSubOptimal == 0
				
                    //attempt to invoke garbage collection to deal with mem usages
                    System.gc();
                    //System.out.println("HINTING TO Garbage Collection");

                } //end of else statement for V with A_prob

            } //end else for if IGHV alignment was NULL
          
        }//end of processing for extra IGHV gene options

        //done with the alignments....
        return null;

    }

    private void displayStateInfo(SymbolList viterbi_state_seq, String UMS, String selected_ums_name)
    {
        String stateResult = "";
        String umsResult = "";
        String preEmitResult = "";
        String matchResult = "";
        //System.out.println("Pre Emission");
        //fostream.println("Pre Emission");
        int ums_index = 0;
        for(int i = 1; i <= viterbi_state_seq.length(); i++)
        {
            Annotation annotation = viterbi_state_seq.symbolAt(i).getAnnotation();
            StateInfo si = (StateInfo)annotation.getProperty(null);
            char geneToken = si.geneToken;
            char statePreEmit = si.preEmissionSymbol;
            if(geneToken != 'X')
            {
                if(geneToken == '2' || geneToken == '1')
                {
                    stateResult = stateResult + 'N';
                    preEmitResult = preEmitResult + statePreEmit;
                } else
                if(geneToken == '3' || geneToken == '4' || geneToken == '5' || geneToken == '6')
                {
                    stateResult = stateResult + 'P';
                    preEmitResult = preEmitResult + statePreEmit;
                } else
                {
                    stateResult = stateResult + geneToken;
                    preEmitResult = preEmitResult + statePreEmit;
                }
            }
            //System.out.print(statePreEmit);
            if(geneToken != 'X')
            {
                char UMSEmit = UMS.charAt(ums_index++);
                umsResult = umsResult + UMSEmit;
                if(statePreEmit == '?')
                {
                    matchResult = matchResult + '?';
                } else
                if(statePreEmit == UMSEmit)
                {
                    matchResult = matchResult + '.';
                } else
                {
                    matchResult = matchResult + '|';
                }
            }
        }

        //System.out.println("******** Result of aligning " + selected_ums_name + "  *******" + "\n");
        //System.out.println("\n");
        //System.out.println("emitting state : " + stateResult + "\n");
        //System.out.println("ums region     : " + umsResult + "\n");
        //System.out.println("nucl. in gene  : " + preEmitResult + "\n");
        //System.out.println("match result   : " + matchResult + "\n");
        
        //fostream.println("******** Result of aligning " + selected_ums_name + "  *******" + "\n");
        //fostream.out.println("\n");
        //fostream.out.println("emitting state : " + stateResult + "\n");
        //fostream.out.println("ums region     : " + umsResult + "\n");
        //fostream.out.println("nucl. in gene  : " + preEmitResult + "\n");
        //fostream.out.println("match result   : " + matchResult + "\n");
		
		//clean-up
		stateResult = null;
        umsResult = null;
        preEmitResult = null;
        matchResult = null;
    }

    private void displayEmissionInfo(SymbolList viterbi_state_seq, String UMS, String selected_ums_name, double statePathProbability, String htmlFileName, int chainType, int reverseComplement)
    {
        int state_index = 1;
        int UMS_index = 0;
        PrintWriter stream = null;
        try {
	        if (stream == null) {
	        	String outputFileName = htmlFileName + ".fostream";
	        	boolean append = true;
	        	stream = new PrintWriter (new FileWriter (outputFileName, append));
	        }
        }
        catch(Exception ex)
        {
            System.out.println("couldn't open new file: " + ex.getMessage());
        }
        
        
        //System.out.println("\n*** Details ***\n");
        //System.out.println("\n");
        //System.out.println("State Path Probability = " + statePathProbability);
        //System.out.println("\n");
        //System.out.println("\n");
        stream.println("*** Details ***");
        //stream.println("\n");
        stream.println("State Path Probability = " + statePathProbability);
        //stream.println("\n");
        //stream.println("\n");
        
        GeneInfo VGeneInfo = displayGene(UMS, viterbi_state_seq, 'V', stream);
        RegionInfo n1Region = displayNRegion(UMS, viterbi_state_seq, '1', stream);
        RegionInfo n2Region = null;
        //RegionInfo p1Region = null; 
        //RegionInfo p2Region = null;
        //RegionInfo p3Region = null; 
        //RegionInfo p4Region = null;
        GeneInfo DGeneInfo = null;
        if (chainType == 1) {
        	//p1Region = displayNRegion(UMS, viterbi_state_seq, '3', stream);
            //p2Region = displayNRegion(UMS, viterbi_state_seq, '4', stream);
            DGeneInfo = displayGene(UMS, viterbi_state_seq, 'D', stream);
            //p3Region = displayNRegion(UMS, viterbi_state_seq, '5', stream);
            n2Region = displayNRegion(UMS, viterbi_state_seq, '2', stream);
            //p4Region = displayNRegion(UMS, viterbi_state_seq, '6', stream);
        }
        GeneInfo JGeneInfo = displayGene(UMS, viterbi_state_seq, 'J', stream);
        PostAlignProcessing postAlign = new PostAlignProcessing(dProbsFilename);
        if (chainType == 1) {// heavy chain
        	boolean acceptedDGene = postAlign.acceptDGene(DGeneInfo, dGeneAcceptanceType, n1Region, n2Region);
            if(!acceptedDGene) {
                DGeneInfo.acceptedAlignment = false;
            }
        }
        Vector stopCodons = postAlign.noStopCodons(VGeneInfo, UMS);
        if(stopCodons.size() == 0)
        {
            boolean sequenceHasStopCodons = false;
            //System.out.println("NO STOP CODONS");
            stream.println("NO STOP CODONS");
        } else
        {
            boolean sequenceHasStopCodons = true;
            //System.out.println("/n****/n*****/n**** stop codons\n*****\n");
            stream.println("**** stop codons *****");
            Codon stopCodon = (Codon)stopCodons.get(0);
            //System.out.println("at pos: " + stopCodon.nuclSeqPos);
            //System.out.println("codon: " + new String(stopCodon.codon));
            stream.println("at pos: " + stopCodon.nuclSeqPos);
            stream.println("codon: " + new String(stopCodon.codon));
        }
        JGeneReadingFrame jGeneReadingFrame = new JGeneReadingFrame(UMS, JGeneInfo, VGeneInfo);
        boolean isJGeneInFrame = jGeneReadingFrame.isJGeneInFrame();
        int relativeMotifPosition = jGeneReadingFrame.getWGXGmotifInAlignedJGenePos();
        
        //htmlWriter.createDocument(htmlFileName, selected_ums_name, chainType, VGeneInfo, p1Region, n1Region, p2Region, DGeneInfo, p3Region, n2Region, p4Region, JGeneInfo, selected_ums_name, UMS, statePathProbability, alignmentType, stopCodons, relativeMotifPosition, isJGeneInFrame, reverseComplement);
	htmlWriter.createDocument(htmlFileName, selected_ums_name, chainType, VGeneInfo, n1Region, DGeneInfo, n2Region, JGeneInfo, selected_ums_name, UMS, statePathProbability, alignmentType, stopCodons, relativeMotifPosition, isJGeneInFrame, reverseComplement, dProbsFilename);
        if(htmlWriter != null)
        {
            htmlWriter.closeHTMLFile();
        }
        
        if(stream != null)
        {        
        	//System.out.println("Mutational File Output Closed for editing");
        	stream.close();
        	stream = null;
        	
        }
		
		//clean-up
		VGeneInfo = null;
		n1Region = null;
		n2Region = null;
		DGeneInfo = null;
		JGeneInfo = null;
		postAlign = null;
		stopCodons = null;
		jGeneReadingFrame = null;
    }

    private RegionInfo displayNRegion(String ums, SymbolList viterbi_state_seq, char state_token, PrintWriter stream)
    {
        //System.out.println("in display N region");
        String GeneName = "";
        String UMS = "emitted sequence:  ";
        String regionString = "";
        int ums_start_location = -1;
        int ums_end_location = -1;
        int ums_index = 0;
        int state_index = 1;
        boolean foundStateRegion = false;
        Annotation anno;
        StateInfo si;
        char geneToken;
        do
        {
            anno = viterbi_state_seq.symbolAt(state_index).getAnnotation();
            si = (StateInfo)anno.getProperty(null);
            geneToken = si.geneToken;
            if(geneToken == state_token)
            {
                foundStateRegion = true;
                ums_start_location = ums_index;
            } else
            if(geneToken == 'X')
            {
                state_index++;
            } else
            {
                state_index++;
                ums_index++;
            }
        } while(state_index <= viterbi_state_seq.length() && !foundStateRegion);
        if(!foundStateRegion)
        {
            return null;
        }
        for(; geneToken == state_token; geneToken = si.geneToken)
        {
            UMS = UMS + ums.charAt(ums_index);
            regionString = regionString + ums.charAt(ums_index);
            state_index++;
            ums_index++;
            if(state_index > viterbi_state_seq.length())
            {
                break;
            }
            anno = viterbi_state_seq.symbolAt(state_index).getAnnotation();
            si = (StateInfo)anno.getProperty(null);
        }

        ums_end_location = ums_index - 1;
        anno = viterbi_state_seq.symbolAt(state_index - 1).getAnnotation();
        si = (StateInfo)anno.getProperty(null);
        GeneName = si.geneName + "  input start index:  " + ums_start_location + "  input end index:  " + ums_end_location;
        //System.out.println(GeneName);
        //System.out.println(UMS + "\n");
        //System.out.println("\n");
        stream.println(GeneName);
        stream.println(UMS);
        //stream.println("\n");
        
        RegionInfo tempRegionInfo = new RegionInfo(regionString, ums_start_location, ums_end_location);
        return tempRegionInfo;
    }

    private GeneInfo displayGene(String ums, SymbolList viterbi_state_seq, char state_token, PrintWriter stream)
    {
        String GeneName = "";
        String Nucleotide_Numbers = "";
        String Gene = "germline:   ";
        String geneString = "";
        String UMS = "input:      ";
        String umsString = "";
        String Match = "match:      ";
        int first_nucleotide_number = -1;
        int last_nucleotide_number = -1;
        int ums_index = 0;
        int state_index = 1;
        boolean foundStateRegion = false;
        Annotation anno;
        StateInfo si;
        char geneToken;
        do
        {
            anno = viterbi_state_seq.symbolAt(state_index).getAnnotation();
            si = (StateInfo)anno.getProperty(null);
            geneToken = si.geneToken;
            if(geneToken == state_token)
            {
                foundStateRegion = true;
                first_nucleotide_number = si.getNucleotidePosition();
            } else
            if(geneToken == 'X')
            {
                state_index++;
            } else
            {
                state_index++;
                ums_index++;
            }
        } while(state_index <= viterbi_state_seq.length() && !foundStateRegion);
        if(!foundStateRegion)
        {
            return null;
        }
        for(; geneToken == state_token; geneToken = si.geneToken)
        {
            UMS = UMS + ums.charAt(ums_index);
            umsString = umsString + ums.charAt(ums_index);
            Gene = Gene + si.preEmissionSymbol;
            geneString = geneString + si.preEmissionSymbol;
            if(si.preEmissionSymbol == ums.charAt(ums_index))
            {
                Match = Match + '.';
            } else
            {
                Match = Match + '|';
            }
            state_index++;
            ums_index++;
            if(state_index > viterbi_state_seq.length())
            {
                break;
            }
            anno = viterbi_state_seq.symbolAt(state_index).getAnnotation();
            si = (StateInfo)anno.getProperty(null);
        }

        anno = viterbi_state_seq.symbolAt(state_index - 1).getAnnotation();
        si = (StateInfo)anno.getProperty(null);
        last_nucleotide_number = si.getNucleotidePosition();
        GeneName = si.geneName + "  length of complete gene:  " + si.completeGeneLength + "  start: " + first_nucleotide_number + "  end: " + last_nucleotide_number;
        //System.out.println(GeneName);
        //System.out.println(UMS + "\n");
        //System.out.println(Gene + "\n");
        //System.out.println(Match + "\n");
        //System.out.println("\n");
        
        if (stream != null) {
	        stream.println(GeneName);
	        stream.println(UMS);
	        stream.println(Gene);
	        stream.println(Match);
	        //stream.println("\n");
        }
        int completeGeneLength = si.completeGeneLength;
        GeneInfo tempGeneInfo = new GeneInfo(si.geneName, geneString, umsString, completeGeneLength, first_nucleotide_number, last_nucleotide_number);
        return tempGeneInfo;
    }

    private String removeFastaStyleWhiteSpace(String string)
    {
        String noWhiteSpaceString = "";
        for(int i = 0; i < string.length(); i++)
        {
            char curr = string.charAt(i);
            if(curr != '-')
            {
                noWhiteSpaceString = noWhiteSpaceString + curr;
            }
        }

        return noWhiteSpaceString;
    }
	
	private String getReverseComp(String seqStr) {
		//return the reverse complement of a sequence
		String rcSeqStr = "";
		
		//make the input sequence lowercase
        seqStr.toLowerCase();
        
		//work through the input sequence string and create the reverse comp
		for (int i = (seqStr.length()-1); i >= 0; i--) {
        	char currNt = seqStr.charAt(i);
        	if (currNt == 'a') {
        		rcSeqStr += 't';
        	} else if (currNt == 't') {
        		rcSeqStr += 'a';
        	} else if (currNt == 'c') {
        		rcSeqStr += 'g';
        	} else if (currNt == 'g') {
        		rcSeqStr += 'c';
        	} else {
        		rcSeqStr += currNt;
        	}
        }
		
		//return the rc of the input sequence
		return rcSeqStr;
	}

    public void run()
    {
         //System.out.println("about to search for the v genes");
         
        startAlignment();
    }
}
