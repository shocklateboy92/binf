package iHMMuneAlign;

import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Vector;

// Referenced classes of package iHMMuneAlign:
//            HTML_Writer, GeneInfo, RegionInfo

public class HTMLforCSV_Writer extends HTML_Writer
{

    static boolean DEBUG = false;

    public HTMLforCSV_Writer()
    {
    }

    public void createErrDocument(String fileName, String sequenceName, byte alignmentType, String errorType, String vOutput)
    {
      // creates the output when the model creation hasn't been sucessful
      //for example an A_probability hasn't been able to be calculated due to a short IGHV gene
      //open up the output file
      if(fostream == null)
        {
            createFile(fileName, alignmentType);

            if (DEBUG) {
                System.out.println("(HMTLforCSV_Writer.java:createErrDocument) Creating a file: " + fileName + "for alignment type: " + alignmentType);
            }
        }

        if(fostream == null)
        {
            throw new Error("Fo stream is null ?????????");
        }

      //create place holders for the different outputs
      String geneNamesString = "\"" + errorType + "\";\"NA\";\"NA\"";
      String regionInfoString = "\"" + vOutput + "\";\"NA\";\"NA\";\"NA\";\"NA\"";
      String mutationsInfoString = "\"NA\";\"NA\";\"NA\"";
      String NandXnucleotidesString = "\"NA\"";
      String vGeneStartString = "\"NA\"";
      String noOfStopCodonsString = "\"NA\"";
      String isJGeneInFrameString = "\"NA\"";
      String sequenceNameString = "\"" + sequenceName + "\"";
      String dProb = "\"NA\"";
      String score = "\"NA\"";
      String rc = "\"NA\"";

      // call the function to write the deatils to a file
      //writeToFile(sequenceNameString, regionInfoString, geneNamesString, mutationsInfoString, NandXnucleotidesString, isJGeneInFrameString, vGeneStartString, noOfStopCodonsString, dProb, score, rc);

      // create an output string matching the style of the output format
      String finalOutputString = "";
      finalOutputString = finalOutputString + sequenceNameString + ";" + geneNamesString + ";";
      finalOutputString = finalOutputString + regionInfoString + ";" + mutationsInfoString + ";";
      finalOutputString = finalOutputString + NandXnucleotidesString + ";" + isJGeneInFrameString + ";";
      finalOutputString = finalOutputString + vGeneStartString + ";" + noOfStopCodonsString + ";" + dProb + ";" + score + ";" + rc;

      //output the csv entry for the given sequence
      fostream.println(finalOutputString);
      //System.out.println(finalOutputString);

      //output the csv entry for the given sequence
      //fostream.println("blah" + finalOutputString);
    }

    public void createDocument(String fileName, String sequenceName, int chainType, GeneInfo VGene, RegionInfo P1, RegionInfo N1, RegionInfo P2, GeneInfo DGene,
            RegionInfo P3, RegionInfo N2, RegionInfo P4, GeneInfo JGene, String selectedSequenceName, String selectedSequenceString, double score, byte alignmentType, Vector stopCodons, int relativeMotifPosition, boolean isJGeneInFrame, int reverseComplement, String dProbsFilename)
    {
        if(fostream == null)
        {
            createFile(fileName, alignmentType);

            if (DEBUG) {
                System.out.println("(HTMLforCSV_Writer.java:createDocument) Creating a document with name: " + fileName + " for alignment type: " + alignmentType);
            }
        }

        if(fostream == null)
        {
            throw new Error("Fo stream is null ?????????");
        }

        int completeVGeneLength = VGene.gene_length;
        //String VGeneInfoString = getVGeneString(VGene, completeVGeneLength);  changed to complete V length rather than just end
        String VGeneInfoString = getDorJGeneString(VGene, false, -1);
        boolean acceptedDGene = false;
        String DGeneInfoString = "NA"; //will output NA for kappa and lambda chain alignments for the D
        if (chainType == 1) {// heavy chain
	        acceptedDGene = DGene.acceptedAlignment;
	        DGeneInfoString = getDorJGeneString(DGene, false, -1);

                if (DEBUG) {
                    System.out.println("(HTMLforCSV_Writer.java:createDocument) acceptedDGene set to " + acceptedDGene + " and D gene info string is: " + DGeneInfoString);
                }
        }
	String JGeneInfoString = getDorJGeneString(JGene, true, relativeMotifPosition);
        String unidentifiedRegion = "";
        
        if ((DEBUG) && (!acceptedDGene)) {
            System.out.println("(HTMLforCSV_Writer.java:createDocument) No acceptable D gene, creating an 'unidentifed region': " + unidentifiedRegion);
        }

        String n1String = null;
        String n2String = "NA"; //will output NA for kappa and lambda chain alignments for the N2 region
        String p1String = null;
        String p2String = null;
        String p3String = null;
        String p4String = null;
        if (chainType == 1) {// heavy chain
	        //if the IGHD has been accepted set the P1 region
	        p1String = null;
	        if(acceptedDGene)
	        {
	            if(P1 == null)
	            {
	                p1String = null;
	            } else
	            {
	                p1String = P1.region_string;
	            }
	        } else   //if IGHD not accepted add the P1 to the unidentified region
	        if(P1 == null)
	        {
	            p1String = null;
	        } else
	        {
	            unidentifiedRegion = unidentifiedRegion + P1.region_string;
	            
	            if (DEBUG) {
	            	System.out.println("(HTMLforCSV_Writer.java:createDocument) Added the following to the unidentifed region (P1): " + P1.region_string);
	            }
	        }
         }
	     //if the IGH region is accepted then set the N1 string if there were N1 nts found
         n1String = null;
         if ((chainType != 1) || (acceptedDGene)) { //for IGL & IGK chains, or where D gene found for the IGH alignment
            if(N1 == null) { //if the alignment included no N1 region then the N1 string for output is also set to null
                n1String = null;
            } else { //otherwise alignment includes an N1 region
                n1String = N1.region_string;
            } 
         } else if (N1 == null) {
            n1String = null; 
         } else { // if no acceptable IGHD set any potential N1 nts to be unided juntion nts
            unidentifiedRegion = unidentifiedRegion + N1.region_string;
            
            if (DEBUG) {
	           	System.out.println("(HTMLforCSV_Writer.java:createDocument) Added the following to the unidentifed region (N1): " + N1.region_string);
	        }
         }
        
	     if (chainType == 1) {   
	        //check if there are p2 nts, if D found set the p2 region
	        p2String = null;
	        if(acceptedDGene)
	        {
	            if(P2 == null)
	            {
	                p2String = null;
	            } else
	            {
	                p2String = P2.region_string;
	            }
	        } else //if d not found set any potential p2 nts to be unidentified junction nts
	        if(P2 == null)
	        {
	            p2String = null;
	        } else
	        {
	            unidentifiedRegion = unidentifiedRegion + P2.region_string;
	            
	            if (DEBUG) {
	            	System.out.println("(HTMLforCSV_Writer.java:createDocument) Added the following to the unidentifed region (P2): " + P2.region_string);
	            }
	        }
	
	        //add the ighd gene to the unidentified region (junction)
	        unidentifiedRegion = unidentifiedRegion + DGene.aligned_ums_string;
	        if (DEBUG) {
	            	System.out.println("(HTMLforCSV_Writer.java:createDocument) Added the following to the unidentifed region (D): " + DGene.aligned_ums_string);
	        }
	        
	        //add the p3 nts if there are any and we have an IGHD accepted
	        p3String = null;
	        if(acceptedDGene)
	        {
	            if(P3 == null)
	            {
	                p3String = null;
	            } else
	            {
	                p3String = P3.region_string;
	            }
	        } else  //or add the P3 nts to the junction region if no confident IGHD identification
	        if(P3 == null)
	        {
	            p3String = null;
	        } else
	        {
	            unidentifiedRegion = unidentifiedRegion + P3.region_string;
	            if (DEBUG) {
	            	System.out.println("(HTMLforCSV_Writer.java:createDocument) Added the following to the unidentifed region (P3): " + P3.region_string);
	            }
	        }
	        
	        //if there are N2 nts found add them to the N2 string if an acceptable ighd was found
	        n2String = null;
	        if(acceptedDGene)
	        {
	            if(N2 == null)
	            {
	                n2String = null;
	            } else
	            {
	                n2String = N2.region_string;
	            }
	        } else  // if there is no D found then add the N2 nts to the junction
	        if(N2 == null)
	        {
	            n2String = null;
	        } else
	        {
	            unidentifiedRegion = unidentifiedRegion + N2.region_string;
	            if (DEBUG) {
	            	System.out.println("(HTMLforCSV_Writer.java:createDocument) Added the following to the unidentifed region (N2): " + N2.region_string);
	            }
	        }
	        
	        //deal with the p4 nts
	        p4String = null;
	        if(acceptedDGene)  // if there is a D then added to the p4 region
	        {
	            if(P4 == null)
	            {
	                p4String = null;
	            } else
	            {
	                p4String = P4.region_string;
	            }
	        } else   // if there is no confident IGHD then add the p4 to the unidentified junction
	        if(P4 == null)
	        {
	            p4String = null;
	        } else
	        {
	            unidentifiedRegion = unidentifiedRegion + P4.region_string;
	            if (DEBUG) {
	            	System.out.println("(HTMLforCSV_Writer.java:createDocument) Added the following to the unidentifed region (P4): " + P4.region_string);
	            }
	        }

                if ((DEBUG) && (!acceptedDGene)) {
                    System.out.println("(HTMLforCSV_Writer.java:createDocument) After processing the gene regions across the junction the unidentifed region is: " + unidentifiedRegion);
                }
        }

        //determine the number of mutations in the IGHV
        int VGeneMutations = mutationsInGeneAlignment(VGene.aligned_gene_string, VGene.aligned_ums_string);
        
        //set the number of IGHD mutations to zero
        int DGeneMutations = 0;
        String dProbStr = "\"NA\"";
        if (chainType == 1) {
        	double dProb = 0D;
	        //if an acceptable match to a germline IGHD was found, calculate the number of mutations in the IGHD
	        if(acceptedDGene)
	        {
	            DGeneMutations = mutationsInGeneAlignment(DGene.aligned_gene_string, DGene.aligned_ums_string);
	        }
	        
	        //get the probability associated with the d if dynamic has been used
	        String dString = DGene.aligned_gene_string;
	        int dLength = dString.length();
	        String junction = n1String + n2String + DGene.aligned_gene_string;
	        int junctionLength = junction.length();
	        PostAlignProcessing postAlign = new PostAlignProcessing(dProbsFilename);
	        dProb = postAlign.getAcceptableDynamic(dLength, DGeneMutations, junctionLength, dProbsFilename);
	        
	        dProbStr = "" + dProb;
	
	        if (acceptedDGene) {
	            dProbStr = "\"" + dProbStr + "\"";
	        } else {
	            dProbStr = "\"ua - " + dProbStr + "\"";
	        }
        }

        //determine the number of mutations in the IGHJ alignment
        int JGeneMutations = mutationsInGeneAlignment(JGene.aligned_gene_string, JGene.aligned_ums_string);
        
        //build the output string for the mutation counts for the IGHV, IGHD and IGHJ
        String mutationsInfoString = createMutationsInfoString(VGeneMutations, DGeneMutations, JGeneMutations, chainType);
        
        //counts the number of unknown (x and n) nucleotides present in the input sequence
        int NandXnucleotideCount = getNandXnucleotideCount(selectedSequenceString);
        //formats the output line for the unknown nucleotide count
        String NandXnucleotidesString = createNandXNuclString(NandXnucleotideCount);
        //creates the string that prints whether the IGHJ gene is in frame or not
        String isJGeneInFrameString = createIsJGeneInFrameString(isJGeneInFrame);
        //creates the string that states the number of stop codons present
        String noOfStopCodonsString = "\"" + stopCodons.size() + "\"";

        //builds the output string that contains the gene names
        String geneNamesString;
        //if an acceptable match to the IGHD is found then the string is created with the IGHV,D and J
        if (chainType == 1) {
	        if(acceptedDGene)
	        {
	            geneNamesString = createGeneNamesString(VGene.gene_name, DGene.gene_name, JGene.gene_name);
	        } else  //if an acceptable match to the IGHD was NOT found then IGHV and IGHJ added, and D replaced by 'no dgene alignment'
	        {
	            geneNamesString = createGeneNamesString(VGene.gene_name, "NO_DGENE_ALIGNMENT", JGene.gene_name);
	        }
        } else {
        	geneNamesString = createGeneNamesString(VGene.gene_name, JGene.gene_name);
        }
        //builds the string for the sequence name
        String sequenceNameString = createSequenceNameString(selectedSequenceName);
        
        //builds the string that contains the sequence region information
        String regionInfoString = "";
        //if acceptable IGHD then bring together all the sequence regions
        if ((chainType != 1) || (acceptedDGene))
        {
            regionInfoString = createRegionInfoString(VGeneInfoString, DGeneInfoString, JGeneInfoString, n1String, n2String, p1String, p2String, p3String, p4String);
        } else   //if no confident IGHD the unidentified region that is made up of N and D is passed instead
        {
            regionInfoString = createRegionInfoString(VGeneInfoString, unidentifiedRegion, JGeneInfoString, n1String, n2String, p1String, p2String, p3String, p4String);
        }
        
        //get the start position of matches between the input IGHV and the germline
        //useful later for introducing the imgt gaps
        int vGeneStart = VGene.start_gene_pos;
        String vGeneStartString = "\"" + vGeneStart + "\"";
        
        // string for the score (state path probability)
        String scoreStr = "\"" + score + "\"";

        // string to indicate whether the input seq is a reverse complement of the original gene seq
        String rcStr = "\"" + reverseComplement + "\"";

        //write all the strings that have been built to a file using the writeToFile funstion
        writeToFile(sequenceNameString, regionInfoString, geneNamesString, mutationsInfoString, NandXnucleotidesString, isJGeneInFrameString, vGeneStartString, noOfStopCodonsString, dProbStr, scoreStr, rcStr);
    }

	public void createDocument(String fileName, String sequenceName, int chainType, GeneInfo VGene, RegionInfo N1, GeneInfo DGene, RegionInfo N2, GeneInfo JGene, String selectedSequenceName, String selectedSequenceString, double score, byte alignmentType, Vector stopCodons, int relativeMotifPosition, boolean isJGeneInFrame, int reverseComplement, String dProbsFilename)
	{
		if(fostream == null)
        {
            createFile(fileName, alignmentType);
            
            if (DEBUG) {
                System.out.println("(HTMLforCSV_Writer.java:createDocument) Creating a document with name: " + fileName + " for alignment type: " + alignmentType);
            }
        }
		
        if(fostream == null)
        {
            throw new Error("Fo stream is null ?????????");
        }
		
        int completeVGeneLength = VGene.gene_length;
        //String VGeneInfoString = getVGeneString(VGene, completeVGeneLength);  changed to complete V length rather than just end
        String VGeneInfoString = getDorJGeneString(VGene, false, -1);
        boolean acceptedDGene = false;
        String DGeneInfoString = "NA";
        if (chainType == 1) {// heavy chain
	        acceptedDGene = DGene.acceptedAlignment;
	        DGeneInfoString = getDorJGeneString(DGene, false, -1);
        }
	    String JGeneInfoString = getDorJGeneString(JGene, true, relativeMotifPosition);
        String unidentifiedRegion = "";
        
        if ((DEBUG) && (!acceptedDGene)) {
            System.out.println("(HTMLforCSV_Writer.java:createDocument) No acceptable D gene, creating an 'unidentifed region': " + unidentifiedRegion);
        }
		
        String n1String = null;
        String n2String = "NA";
        //String p1String = null;
        //String p2String = null;
        //String p3String = null;
        //String p4String = null;
        /*if (chainType == 1) {// heavy chain
	        //if the IGHD has been accepted set the P1 region
	        p1String = null;
	        if(acceptedDGene)
	        {
	            if(P1 == null)
	            {
	                p1String = null;
	            } else
	            {
	                p1String = P1.region_string;
	            }
	        } else   //if IGHD not accepted add the P1 to the unidentified region
				if(P1 == null)
				{
					p1String = null;
				} else
				{
					unidentifiedRegion = unidentifiedRegion + P1.region_string;
				}
		}*/
		//if the IGH region is accepted then set the N1 string if there were N1 nts found
		n1String = null;
		if ((chainType != 1) || (acceptedDGene)) {  //not heavy chain OR acceptedD true
            if(N1 == null) {
                n1String = null;
            } else {
                n1String = N1.region_string;
            }
			
		} else  {// if no acceptable IGHD set any potential N1 nts to be unided juntion nts
			if(N1 == null) {
				n1String = null;
			} else  {
				unidentifiedRegion = unidentifiedRegion + N1.region_string;
				//add the ighd gene to the unidentified region (junction)
				if (DEBUG) {
	            	System.out.println("(HTMLforCSV_Writer.java:createDocument) Added the following to the unidentifed region (N1): " + N1.region_string);
	            }
			}
		}
		
        
		if (chainType == 1) {   
	        /*
			//check if there are p2 nts, if D found set the p2 region
	        p2String = null;
	        if(acceptedDGene)
	        {
	            if(P2 == null)
	            {
	                p2String = null;
	            } else
	            {
	                p2String = P2.region_string;
	            }
	        } else //if d not found set any potential p2 nts to be unidentified junction nts
				if(P2 == null)
				{
					p2String = null;
				} else
				{
					unidentifiedRegion = unidentifiedRegion + P2.region_string;
				}
			*/
		

	        
	        //add the p3 nts if there are any and we have an IGHD accepted
	        /*p3String = null;
	        if(acceptedDGene)
	        {
	            if(P3 == null)
	            {
	                p3String = null;
	            } else
	            {
	                p3String = P3.region_string;
	            }
	        } else  //or add the P3 nts to the junction region if no confident IGHD identification
				if(P3 == null)
				{
					p3String = null;
				} else
				{
					unidentifiedRegion = unidentifiedRegion + P3.region_string;
				}
	        */
	        // D
	        if (!acceptedDGene) {
	        	unidentifiedRegion = unidentifiedRegion + DGene.aligned_ums_string;
	        	if (DEBUG) {
	            		System.out.println("(HTMLforCSV_Writer.java:createDocument) Added the following to the unidentifed region (D)): " + DGene.aligned_ums_string);
	            }
	        }
	        
	        //if there are N2 nts found add them to the N2 string if an acceptable ighd was found
	        n2String = null;
	        if(acceptedDGene)
	        {
	            if(N2 == null)
	            {
	                n2String = null;
	            } else
	            {
	                n2String = N2.region_string;
	            }
	        } else  // if there is no D found then add the N2 nts to the junction
				if(N2 == null)
				{
					n2String = null;
				} else
				{
					
					unidentifiedRegion = unidentifiedRegion + N2.region_string;
					if (DEBUG) {
	            		System.out.println("(HTMLforCSV_Writer.java:createDocument) Added the following to the unidentifed region (N2)): " + N2.region_string);
	            	}
				}
			
			if ((DEBUG) && (!acceptedDGene)) {
                   System.out.println("(HTMLforCSV_Writer.java:createDocument) After processing the gene regions across the junction the unidentifed region is: " + unidentifiedRegion);
             }	
	        
	        //deal with the p4 nts
	        /*p4String = null;
	        if(acceptedDGene)  // if there is a D then added to the p4 region
	        {
	            if(P4 == null)
	            {
	                p4String = null;
	            } else
	            {
	                p4String = P4.region_string;
	            }
	        } else   // if there is no confident IGHD then add the p4 to the unidentified junction
				if(P4 == null)
				{
					p4String = null;
				} else
				{
					unidentifiedRegion = unidentifiedRegion + P4.region_string;
				}
			 */
        }
		
        //determine the number of mutations in the IGHV
        int VGeneMutations = mutationsInGeneAlignment(VGene.aligned_gene_string, VGene.aligned_ums_string);
        
        //set the number of IGHD mutations to zero
        int DGeneMutations = 0;
        String dProbStr = "\"NA\"";
        if (chainType == 1) {
        	double dProb = 0D;
	        //if an acceptable match to a germline IGHD was found, calculate the number of mutations in the IGHD
	        if(acceptedDGene)
	        {
	            DGeneMutations = mutationsInGeneAlignment(DGene.aligned_gene_string, DGene.aligned_ums_string);
	        }
	        
	        //get the probability associated with the d if dynamic has been used
	        String dString = DGene.aligned_gene_string;
	        int dLength = dString.length();
	        String junction = n1String + n2String + DGene.aligned_gene_string;
	        int junctionLength = junction.length();
	        PostAlignProcessing postAlign = new PostAlignProcessing(dProbsFilename);
	        dProb = postAlign.getAcceptableDynamic(dLength, DGeneMutations, junctionLength, dProbsFilename);
	        
	        dProbStr = "" + dProb;
			
	        if (acceptedDGene) {
	            dProbStr = "\"" + dProbStr + "\"";
	        } else {
	            dProbStr = "\"ua - " + dProbStr + "\"";
	        }
        }
		
        //determine the number of mutations in the IGHJ alignment
        int JGeneMutations = mutationsInGeneAlignment(JGene.aligned_gene_string, JGene.aligned_ums_string);
        
        //build the output string for the mutation counts for the IGHV, IGHD and IGHJ
        String mutationsInfoString = createMutationsInfoString(VGeneMutations, DGeneMutations, JGeneMutations, chainType);
        
        //counts the number of unknown (x and n) nucleotides present in the input sequence
        int NandXnucleotideCount = getNandXnucleotideCount(selectedSequenceString);
        //formats the output line for the unknown nucleotide count
        String NandXnucleotidesString = createNandXNuclString(NandXnucleotideCount);
        //creates the string that prints whether the IGHJ gene is in frame or not
        String isJGeneInFrameString = createIsJGeneInFrameString(isJGeneInFrame);
        //creates the string that states the number of stop codons present
        String noOfStopCodonsString = "\"" + stopCodons.size() + "\"";
		
        //builds the output string that contains the gene names
        String geneNamesString;
        //if an acceptable match to the IGHD is found then the string is created with the IGHV,D and J
        if (chainType == 1) {
	        if(acceptedDGene)
	        {
	            geneNamesString = createGeneNamesString(VGene.gene_name, DGene.gene_name, JGene.gene_name);
	        } else  //if an acceptable match to the IGHD was NOT found then IGHV and IGHJ added, and D replaced by 'no dgene alignment'
	        {
	            geneNamesString = createGeneNamesString(VGene.gene_name, "NO_DGENE_ALIGNMENT", JGene.gene_name);
	        }
        } else {
        	geneNamesString = createGeneNamesString(VGene.gene_name, JGene.gene_name);
        }
        //builds the string for the sequence name
        String sequenceNameString = createSequenceNameString(selectedSequenceName);
        
        //builds the string that contains the sequence region information
        String regionInfoString = "";
        //if acceptable IGHD then bring together all the sequence regions
        if ((chainType != 1) || (acceptedDGene))
        {
            regionInfoString = createRegionInfoString(VGeneInfoString, DGeneInfoString, JGeneInfoString, n1String, n2String);
        } else   //if no confident IGHD the unidentified region that is made up of N and D is passed instead
        {
            regionInfoString = createRegionInfoString(VGeneInfoString, unidentifiedRegion, JGeneInfoString, n1String, n2String);
        }
        
        //get the start position of matches between the input IGHV and the germline
        //useful later for introducing the imgt gaps
        int vGeneStart = VGene.start_gene_pos;
        String vGeneStartString = "\"" + vGeneStart + "\"";
        
        // string for the score (state path probability)
        String scoreStr = "\"" + score + "\"";
		
        // string to indicate whether the input seq is a reverse complement of the original gene seq
        String rcStr = "\"" + reverseComplement + "\"";
		
        //write all the strings that have been built to a file using the writeToFile funstion
        writeToFile(sequenceNameString, regionInfoString, geneNamesString, mutationsInfoString, NandXnucleotidesString, isJGeneInFrameString, vGeneStartString, noOfStopCodonsString, dProbStr, scoreStr, rcStr);
		
	}
    //new function to bring together all the strings
    private void writeToFile(String sequenceNameString, String regionInfoString, String geneNamesString, String mutationsInfoString, String NandXnucleotidesString, String isJGeneInFrameString, String vGeneStartString, String noOfStopCodonsString, String dProb, String score, String rc)
    { // bring all the output strings together into a single line format
      String finalOutputString = "";
      finalOutputString = finalOutputString + sequenceNameString + ";" + geneNamesString + ";";
      finalOutputString = finalOutputString + regionInfoString + ";" + mutationsInfoString + ";";
      finalOutputString = finalOutputString + NandXnucleotidesString + ";" + isJGeneInFrameString + ";";
      finalOutputString = finalOutputString + vGeneStartString + ";" + noOfStopCodonsString + ";" + dProb + ";" + score + ";" + rc;

      //output the csv entry for the given sequence
      fostream.print(finalOutputString);
    }

    //counts the number of unknown nucleotides in a given sequence
    private int getNandXnucleotideCount(String selectedSequenceString)
    {
        int NandXnucleotideCount = 0;
        for(int i = 0; i < selectedSequenceString.length(); i++)
        {
            char currNucl = selectedSequenceString.charAt(i);
            if(currNucl == 'n' || currNucl == 'x')
            {
                NandXnucleotideCount++;
            }
        }

        return NandXnucleotideCount;
    }

    //counts the number of mutations between the input sequence and the germline gene sequence
    private int mutationsInGeneAlignment(String geneAlignment, String sequenceAlignment)
    {
        int mutationsFound = 0;
        for(int i = 0; i < sequenceAlignment.length(); i++)
        {
            char sequenceNucl = sequenceAlignment.charAt(i);
            char geneNucl = geneAlignment.charAt(i);
            if(geneNucl != sequenceNucl && sequenceNucl != 'x' && sequenceNucl != 'n')
            {
                mutationsFound++;
            }
        }

        return mutationsFound;
    }

    //formats the output line about the IGHJ being in frame or not
    //get rid of the html formatting and add the csv formatting
    private String createIsJGeneInFrameString(boolean isJGeneInFrame)
    {
        String result = "";
        if(isJGeneInFrame)
        {
            result = result + "\"true\"";
        } else
        {
            result = result + "\"false\"";
        }
        //result = result + "</pre>";
        return result;
    }

    //creates the outputstring the describes the number of unknown nucleotides in the sequence
    //get rid of the html formatting and add the csv formatting
    private String createNandXNuclString(int NandXnucleotideCount)
    {
        String NandXnucleotidesString = "";
        NandXnucleotidesString = NandXnucleotidesString + "\"" + NandXnucleotideCount;
        NandXnucleotidesString = NandXnucleotidesString + "\"";
        return NandXnucleotidesString;
    }

    //creates the output string that has the mutation counts for each of the IGHV/D/J
    //remove the html formatting and the header line and add in the csv formatting
    private String createMutationsInfoString(int VGeneMutations, int DGeneMutations, int JGeneMutations, int chainType)
    {
        String mutationsInfoString = "";
        //mutationsInfoString = mutationsInfoString + "Vmut\tDmut\tJmut\n";
        if (chainType == 1) {
	        mutationsInfoString = mutationsInfoString + "\"" + VGeneMutations + "\"" + ";\"" + DGeneMutations + "\"" + ";\"" + JGeneMutations;
	    } else {
	        mutationsInfoString = mutationsInfoString + "\"" + VGeneMutations + "\"" + ";\"" + "NA" + "\"" + ";\"" + JGeneMutations;
		}
        mutationsInfoString = mutationsInfoString + "\"";
        return mutationsInfoString;
    }

    //creates the output string that contains the names of the germline genes
    //also formats the header line (which is not required)
    //removed html and header formatting and add in the csv formatting
    private String createGeneNamesString(String VGeneName, String DGeneName, String JGeneName)
    {
        String GeneNameString = "";
        GeneNameString = GeneNameString + "\"" + VGeneName + "\";\"" + DGeneName + "\";\"" + JGeneName;
        GeneNameString = GeneNameString + "\"";
        return GeneNameString;
    }

    private String createGeneNamesString(String VGeneName, String JGeneName)
    {
        String GeneNameString = "";
        GeneNameString = GeneNameString + "\"" + VGeneName + "\";\"" + "\";\"" + JGeneName;
        GeneNameString = GeneNameString + "\"";
        return GeneNameString;
    }

    //creates the header line that contains the accession number
    //remove the html formatting and add in the csv style
    private String createSequenceNameString(String sequenceName)
    {
        String sequenceNameString = "";
        sequenceNameString = sequenceNameString + "\"" + sequenceName + "\"";
        //sequenceNameString = sequenceNameString + "</pre></h3>";
        return sequenceNameString;
    }

    //formats the different regions sequences so that they are aligned under the relevant headings
    //change so that format is just the regions encapsulated by quotation marks and separated by ;
    private String createRegionInfoString(String VRegionString, String DRegionString, String JRegionString, String N1, String N2, String P1, String P2,
            String P3, String P4)
    {
        String regionInfoString = "";
        regionInfoString = regionInfoString + "\"" + VRegionString + "\";";
        //join the P1 to the start of the N1 region
        if(P1 != null)
        {
            N1 = P1 + N1;
        }
        //join the P2 to the end of the N1 region
        if(P2 != null)
        {
           N1 = N1 + P2;
        }
        if(N1 != null)
        {
            regionInfoString = regionInfoString + "\"" + N1 + "\";";
        } else { //we just need to add an empty value if there is no N1 region
            regionInfoString = regionInfoString + "\"\";";
        }

        //add the D region to the outputstring
        regionInfoString = regionInfoString + "\"" + DRegionString + "\";";

        //add the P3 and P4 nucleotides to the N2 region if they are present
        if(P3 != null)
        {
           N2 = P3 + N2;
        }
        if(P4 != null)
        {
           N2 = N2 + P4;
        }
        //output the N2 region if there is one
        if(N2 != null)
        {
           regionInfoString = regionInfoString + "\"" + N2 + "\";";
        } else { //if there is no N2 region we output empty value for the csv
           regionInfoString = regionInfoString + "\"\";";
        }

        //add the IGHJ
        regionInfoString = regionInfoString + "\"" + JRegionString + "\"";
        //regionInfoString = regionInfoString + "</pre>";
        return regionInfoString;
    }

	//create region info for just the V, D, J, N1 & N2
	private String createRegionInfoString(String VRegionString, String DRegionString, String JRegionString, String N1, String N2)
	{
		String regionInfoString = "";
        regionInfoString = regionInfoString + "\"" + VRegionString + "\";";
        
        if(N1 != null)
        {
            regionInfoString = regionInfoString + "\"" + N1 + "\";";
        } else { //we just need to add an empty value if there is no N1 region
            regionInfoString = regionInfoString + "\"\";";
        }
		
        //add the D region to the outputstring
        regionInfoString = regionInfoString + "\"" + DRegionString + "\";";
		
        //output the N2 region if there is one
        if(N2 != null)
        {
			regionInfoString = regionInfoString + "\"" + N2 + "\";";
        } else { //if there is no N2 region we output empty value for the csv
			regionInfoString = regionInfoString + "\"\";";
        }
		
        //add the IGHJ
        regionInfoString = regionInfoString + "\"" + JRegionString + "\"";
        //regionInfoString = regionInfoString + "</pre>";
        return regionInfoString;
	}
	
    //format the IGHV gene string -> was previously a separate function as only wanted to deal with the V gene end rather than full sequence
    //deal with the full sequence by passing the VGENE_END_SIZE as equal to the full sequence length
    private String getVGeneString(GeneInfo VGene, int VGENE_END_SIZE)
    {
        String name = VGene.gene_name;
        String alignedGeneString = VGene.aligned_gene_string;
        String alignedUMSString = VGene.aligned_ums_string;
        int completeGeneLength = VGene.gene_length;
        int startGenePos = VGene.start_gene_pos; //start of the matches between the IGHV input sequence and germline sequence
        int endGenePos = VGene.end_gene_pos; //end of matches between input and germline
        if(alignedGeneString.length() != alignedUMSString.length())
        {
            throw new Error("aligned VGene String and aligned UMS String not of same length");
        }
        String resultGeneInfo = "";
        int alignedGeneLength = alignedGeneString.length(); //length of the aligned portion of the input sequence
        int endGeneOffset = completeGeneLength - endGenePos;  //the portion of the complete input sequence that is the IGHV
        int alignedEndSize = VGENE_END_SIZE - endGeneOffset; // previously used to get rid of most of the V, now just get any leading mismatches
        String alignedGeneEndString = alignedGeneString.substring(alignedGeneLength - alignedEndSize, alignedGeneLength); //get the substrings that are the IGHV
        String alignedUMSEndString = alignedUMSString.substring(alignedGeneLength - alignedEndSize, alignedGeneLength);   //for the input and germline
        int mutationsInAlignmentEnd = 0;
        //count the mutations in the IGHV gene, now for the full length not just the IGHV end
        for(int i = 0; i < VGENE_END_SIZE; i++)
        {
            if(i < alignedEndSize)
            {
                char currGeneNucl = alignedGeneEndString.charAt(i);
                char currUMSNucl = alignedUMSEndString.charAt(i);
                if(currGeneNucl == currUMSNucl)
                {
                    resultGeneInfo = resultGeneInfo + currGeneNucl;
                } else
                {
                    if(currUMSNucl == 'n' || currUMSNucl == 'x')
                    {
                        resultGeneInfo = resultGeneInfo + currUMSNucl;
                    } else
                    {
                        resultGeneInfo = resultGeneInfo + getMutatedNucl(currUMSNucl);
                    }
                    mutationsInAlignmentEnd++;
                }
            } else
            {
                resultGeneInfo = resultGeneInfo + '.';
            }
        }

        return resultGeneInfo;
    }

    //deals with the IGHD and IGHJ sequences, makes the mutations capitalised
    private String getDorJGeneString(GeneInfo DorJGene, boolean isJGene, int jGeneWGXGmotifPos)
    {
        String name = DorJGene.gene_name;
        String alignedDorJGeneString = DorJGene.aligned_gene_string;
        String alignedUMSString = DorJGene.aligned_ums_string;
        int completeGeneLength = DorJGene.gene_length;
        int startGenePos = DorJGene.start_gene_pos;
        int endDorJGenePos = DorJGene.end_gene_pos;
        if(alignedDorJGeneString.length() != alignedUMSString.length())
        {
            throw new Error("aligned DorJGene String and aligned UMS String not of same length");
        }
        String resultGeneInfo = "";
        int mutationsInAlignmentEnd = 0;
        for(int i = 1; i <= completeGeneLength; i++)
        {
            if(i >= startGenePos && i <= endDorJGenePos) {
                int alignedJGeneIndex = i - startGenePos;
                char currGeneNucl = alignedDorJGeneString.charAt(i - startGenePos);
                char currUMSNucl = alignedUMSString.charAt(i - startGenePos);
                if(currGeneNucl == currUMSNucl)  {
                        resultGeneInfo = resultGeneInfo + currGeneNucl;
                } else  {
                    if(currUMSNucl == 'n' || currUMSNucl == 'x') {
                        resultGeneInfo = resultGeneInfo + currUMSNucl;
                    } else {
                        resultGeneInfo = resultGeneInfo + getMutatedNucl(currUMSNucl);
                    }
                    mutationsInAlignmentEnd++;
                }
            } else {
                resultGeneInfo = resultGeneInfo + '.';
            }
        }

        return resultGeneInfo;
    }

    //write mutation nucleotides as UPPERCASE
    private String getMutatedNucl(char mutNucl)
    {
        String result = "" + mutNucl;
        result = result.toUpperCase();
        return result;
    }


}
