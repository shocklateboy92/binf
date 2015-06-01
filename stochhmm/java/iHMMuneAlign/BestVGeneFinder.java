package iHMMuneAlign;

import java.io.*;
import java.lang.Process; 
import java.util.ArrayList;
import java.util.Random;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.OutputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;


import org.biojava.bio.seq.*;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.seq.DNATools;
import org.biojava.utils.ChangeVetoException;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.io.SymbolTokenization;

import org.biojavax.*;
import org.biojavax.bio.seq.io.*;
import org.biojavax.bio.seq.io.RichStreamReader;
import org.biojavax.bio.SimpleBioEntry;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequence.Tools.*;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.biojavax.bio.seq.RichSequence.IOTools.*;
import org.biojavax.bio.seq.io.RichStreamWriter.*;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

// Referenced classes of package iHMMuneAlign:
//            PostAlignmentResult, FastaReader

public class BestVGeneFinder
{

    //define the gap and gap extension penalties for use with jaligner
    final float OPEN_GAP_COST = 10F;
    final float EXTEND_GAP_COST = 5F;
	final int GAPS_ALLOWED = 1;
    final int GAPS_ERROR = 2;
    final int GAPS_IGNORE_SEQUENCE = 3;
    int Gap_Behaviour;
	String geneRepertoireName;
	String blastLocation;
	boolean DEBUGGING = false;
    //org.biojava.bio.seq.Sequence VGenes[];

    public static void main(String args[])
    {
        try {
			BestVGeneFinder VGene_finder = new BestVGeneFinder(new File("v-regions.fa"), "/usr/bin/");
		} catch (Throwable t) {
			t.printStackTrace();
		}
    }

    //takes the germline repertoire FASTA file and uses it to create a BLAST database for alignment against the input sequence
    public BestVGeneFinder(File fastaFile, String blastPath)
    {
        Gap_Behaviour = GAPS_ALLOWED; //GAPS_ALLOWED; //GAPS_IGNORE_SEQUENCE;  //GAPS_ERROR;
        
		try {
			//get the name from the fastaFile
			String fastaFileName = fastaFile.getAbsolutePath();
			//based on the repertoire name, create name of db file
			File dbFastaFileName = new File (fastaFileName + ".nhr");
			
			//get the location of blast
			this.blastLocation = blastPath;
			
			//check if format db has previously been performed
			if (dbFastaFileName.exists()) {
				//file already present
				//System.out.println("(BestVGeneFinder) BLAST database already exists for repertoire:" + dbFastaFileName);
			} else {
				//no blastdb files present so need to format the germline repertoire db
				
				//create the command string to be executed
				//String cmdStrCreateDB = blastLocation + "makeblastdb -in " + fastaFileName + " -dbtype nucl -parse_seqids";
				String cmdStrCreateDB = "formatdb -i " + fastaFileName + " -n " + fastaFileName + " -t " + fastaFileName + " -p F -o T";  //changed to blastall
				if (DEBUGGING) {
					System.out.println("(BestVGeneFinder) Executing command to create BLAST db from repertoire file: " + cmdStrCreateDB);
				}
				
				//execute the command to create the db files
				Runtime rt = Runtime.getRuntime();
				Process makeblastdb = rt.exec(cmdStrCreateDB);
				makeblastdb.waitFor();
				
				//destroy the process once it is complete
				if (makeblastdb != null) {
		    		close(makeblastdb.getOutputStream());
				    close(makeblastdb.getInputStream());
		    		close(makeblastdb.getErrorStream());
		    		makeblastdb.destroy();
		    		makeblastdb = null;
		    		rt = null;
      			}
				
				if (DEBUGGING) {
					System.out.println("(BestVGeneFinder) The repertoire file was made into a BLAST database file for searching: " + fastaFileName);
				}
			}

			this.geneRepertoireName = fastaFileName;
		}
		catch (Throwable t) {
			t.printStackTrace();
		}

		
    }

    //returns the result of the alignment in the format of a post alignment class
    public PostAlignmentResult getResult(RichSequence inputSeq, int resultNo)
    {
		//align the input sequence against the germline repertoire and get the result
		int resultCount = resultNo + 5;
		PostAlignmentResult result = doAlignment(inputSeq, geneRepertoireName, resultNo, resultCount);
		
		//displays the alignment info as formatted by the postalignment result funtion
		//result.displayAlignmentInfo();
            
		//returns the postalignment result object containing the alignment information
		return result;
    }


    private PostAlignmentResult doAlignment(RichSequence inputSeq, String dbName, int resultNo, int resultCount)
    {
		try {
			// write the input sequence to a temp file with FASTA format	
			String tmpFileName = makeFASTAFileFrmSingleSeq(inputSeq, "tmpSeq.input");
			String tmpFileNameOutput = tmpFileName.replace("input","output");
			
			File tmpFileNameFile = new File(tmpFileName);
			File tmpFileNameOutputFile = new File(tmpFileNameOutput);
			//tmpFileNameFile.deleteOnExit();
			//tmpFileNameOutputFile.deleteOnExit();
			
			//set up gap behaviour for blast alignment
			String blastGapStr = "";
			switch (Gap_Behaviour) {
				case 1: //gaps allowed
					//blastGapStr = " -gapopen 16 -gapextend 4 -penalty -5 -reward 2 -window_size 0"; //change later if want to use custom
					blastGapStr = "-g T -G 16 -E 4 -q -5 -r 2"; //changed to blastall args
					break;
				case 2:
					//use default but need to throw error later during processing
					blastGapStr = "-g T -G 16 -E 4 -q -5 -r 2"; //allow the BLAST search to find gaps but throw an error later
					break;
				case 3: 
					//blastGapStr = " -ungapped  -penalty -5 -reward 2 -window_size 0";
					blastGapStr = "-g F -q -5 -r 2";  //changed to blastall args
					break;
				default:
					break;
			}
			
			//set up scoring for blast alignment
			//String blastScoringStr = "";
			
			//create the command string to execute blast
			//String blastCmdStr = blastLocation + "blastn -query "  + tmpFileName + " -out " + tmpFileNameOutput + " -db " + dbName + " -num_alignments 5 -num_descriptions 5 -outfmt 5" + blastGapStr;
			String blastCmdStr = "blastall -p blastn -i " + tmpFileName + " -o " + tmpFileNameOutput + " -d " + dbName + " -b " + resultCount + " -v " + resultCount + " -m 7 " + blastGapStr;
			if (DEBUGGING) {
				System.out.println("(BestVGeneFinder) BLAST command to be executed: " + blastCmdStr);
			}
			
			//create the blast process, execute blast and wait for the process to finish
			Runtime rt = Runtime.getRuntime();
			Process blastn = rt.exec(blastCmdStr);
			blastn.waitFor();
			
			//destroy the process now that it is complete
			if (blastn != null) {
        		close(blastn.getOutputStream());
		        close(blastn.getInputStream());
        		close(blastn.getErrorStream());
        		blastn.destroy();
        		blastn = null;
        		rt = null;
      		}
			
			//send to result parser to extract result number X --> in case of multi-V alignments can request non-best alignment
			PostAlignmentResult[] results = getBLASTResults(tmpFileNameOutput, resultNo, inputSeq, resultCount);
			
			PostAlignmentResult result = null;
			if (resultNo == 0) {  //only the top scoring V gene have been requested
				//create a new BLAST result that incorps all top scoring (where the alignment score is equal top scoring)
				result = createPseduoBLASTResult(results);
			} else { //looking to return a match other then the top scoring
				//just return the requested result, this means that if equal top scoring will just be outputting the x top scoring
				result = results[resultNo];
			}	
			
			
			
			//delete the temp files
			tmpFileNameFile.delete(); //file that contained single FASTA sequence for blast
			tmpFileNameOutputFile.delete();  //file that contained the BLAST results that has now been saved
			
			return result;
			
		} catch (Throwable t) {
			t.printStackTrace();
			return null;
		}
        
    }
	
	//create a fasta format file from given string
	private String makeFASTAFileFrmSingleSeq (RichSequence seq, String FASTAFileName) {
		try {
			//location for a temp file from the input sequence
			Random rnd = new Random();
			FASTAFileName = FASTAFileName + rnd.nextInt();
			FileOutputStream tempFile = new FileOutputStream(FASTAFileName);
			
			//get default ns for writing to FASTA
			Namespace ns = RichObjectFactory.getDefaultNamespace();
			
			//set a tmp id & create a RichSequence <- used if sequence was passed as String
			//String tmpID = "input_query_seq";
			//RichSequence rs = org.biojavax.bio.seq.RichSequence.Tools.createRichSequence(ns, tmpID, seq, DNATools.getDNA());
			
			//write the sequence object to the temp file
			RichSequenceIterator it = new SingleRichSeqIterator(seq);
			RichSequenceFormat fasta = new FastaFormat();
			RichStreamWriter seqOut = new RichStreamWriter(tempFile, fasta);
			seqOut.writeStream(it,ns);
			
			//close the output stream
			if (tempFile != null) {
				tempFile.close();
			}
			
			if (DEBUGGING) {
				System.out.println("(BestVGeneFinder:makeFASTAFrmSingleSeq) input sequence written to file ready for BLASTN: " + FASTAFileName);
			}
			
			return FASTAFileName;
		} catch (Throwable t) {
			t.printStackTrace();
			return null;
		}
	}
	
	//parse through the XML BLAST output file created and extract the details of alignment of input sequence against the germline repertoire
	private PostAlignmentResult[] getBLASTResults (String blastOutputFile, int resultNo, RichSequence inputSeq, int resultCount) {
		try {
            //array to store results
			//BLASTResult[] results = new BLASTResult[5];
			
			//parse the XML file
			String[][] BLASTResultArray = null;
			BLASTResultArray = parseBLASTXMLFile(blastOutputFile, resultCount);
					
			//default values
			String UMS_string_lowercase = "NA";
			String best_alignment_VGene = "NA";
			int best_alignment_row_offset = 0;
			int best_alignment_col_offset = 0;
			int best_alignment_gaps = -1;
			int best_alignment_similarity = 0;
			String UMS_name = "NA";
			String best_VGene_name = "NA";
			float best_align_score = 0F;
			char[] aligned_input = null;

			PostAlignmentResult[] results = new PostAlignmentResult[resultCount];
			
			for (int resultIndex = 0; resultIndex < resultCount; resultIndex++) {
				if (BLASTResultArray[resultIndex][0] != null) {
					//determine if this is a reverse complement alignment based on the value of HIT-FRAME
					int hitFrame = Integer.parseInt(BLASTResultArray[resultIndex][10]);
				
					if (hitFrame == -1) {
						//the alignment was to the reverse complement/minus/anti-sense strand of the germline gene
					
					} else {
						//alignment was to the plus/sense strand of the germline gene
					
						 //get the 'best' IGHV alignment
						 //adjust result index as to retrieve other result as required
						 //int resultIndex = resultNo;
						 //IGHV name
						 best_VGene_name = BLASTResultArray[resultIndex][0];
			
						 //input string
						 UMS_string_lowercase = (inputSeq.seqString()).toLowerCase();
						 UMS_name = inputSeq.getName();
				
						 //offsets
						 best_alignment_row_offset = Integer.parseInt(BLASTResultArray[resultIndex][6]) - 1; //'from' needs to be adjusted to give offset
						 best_alignment_col_offset = Integer.parseInt(BLASTResultArray[resultIndex][8]) - 1; //adjust to give an offset value
						 //gaps
						 best_alignment_gaps = Integer.parseInt(BLASTResultArray[resultIndex][1]);
						 //similarity
						 best_alignment_similarity = Integer.parseInt(BLASTResultArray[resultIndex][12]);
						 //alignment score
						 best_align_score = Float.parseFloat(BLASTResultArray[resultIndex][4]);
						 //aligned input
						 String alignedStr = BLASTResultArray[resultIndex][14];
						 aligned_input = alignedStr.toCharArray();
						 
						 String hit_id = BLASTResultArray[resultIndex][17];
						 
						 //get the full germline sequence
						 //best_alignment_VGene = getSeqFrmBlastDB(hit_id, geneRepertoireName);
						 best_alignment_VGene = getSeqFrmBlastDB(best_VGene_name, geneRepertoireName);
						 
						 if ((best_alignment_gaps > 0) && (Gap_Behaviour == GAPS_ALLOWED)) {
						 	//allowing gaps so need to make changes to the full length germline gene that is to be used for the PostAlignmentResult
						 	//need to if there has been a deletion relative to germline need to remove to pad the input seq or if there has been an 
						 	//insertion relative to the germline need to pad the germline
		
							//get the aligned V that will have '-' for ins in input, use this with the relevant up and downstream regions to create new germline string
						 	String vGeneHit = BLASTResultArray[resultIndex][15];
			
							//gap in the germline string
							if (vGeneHit.indexOf("-") > 0) {		
								//get the full length sequence of the v gene
							 	String fullLengthVGene = best_alignment_VGene;
							 	
							 	//get the index for the start and end of the hit found to the vgene
							 	int VgeneFromIndex = Integer.parseInt(BLASTResultArray[resultIndex][8]);
							 	int VgeneToIndex = Integer.parseInt(BLASTResultArray[resultIndex][9]);
							 	
							 	//determine what needs to be added up and downstream of the aligned portion of the V to make the full length seq that includes indels
							 	String prependForFullLength = fullLengthVGene.substring(0, (VgeneFromIndex-1));
							 	String postpendForFullLength = "";
							 	if (VgeneToIndex < fullLengthVGene.length()) {
							 		postpendForFullLength = fullLengthVGene.substring(VgeneToIndex, fullLengthVGene.length());
							 	}
							 	

							 	//add the up and downstream seq to create the full length V
							 	String vGeneWithIndels = prependForFullLength + vGeneHit + postpendForFullLength;
							 	if (DEBUGGING) {
							 		System.out.println("(BestVGeneFinder) full length V: " + vGeneWithIndels);
							 	}
							 	
							 	//set the newly created sequence to the V sequence
							 	best_alignment_VGene = (vGeneWithIndels.toLowerCase());
							 	
							 	//adjust the name
							 	best_VGene_name = best_VGene_name + " [indels]";
							}
						 	
						 	//gap in the input string
						 	if(alignedStr.indexOf("-") > 0 ) {
							 	//get the aligned input sequence that will have '-' if there has been a deletion from the germline
							 	String fullLengthInput = UMS_string_lowercase;
							 	int inputFromIndex = Integer.parseInt(BLASTResultArray[resultIndex][6]);
							 	int inputToIndex = Integer.parseInt(BLASTResultArray[resultIndex][7]);
							 	String prependForInput = fullLengthInput.substring(0, (inputFromIndex-1));
							 	
							 	String postpendForInput = "";
							 	if (inputToIndex < fullLengthInput.length()) {
							 		postpendForInput = fullLengthInput.substring(inputToIndex, fullLengthInput.length());
							 	}
							 	String inputWithIndels = prependForInput + alignedStr + postpendForInput;
							 	
							 	//make lower case
							 	UMS_string_lowercase = (inputWithIndels.toLowerCase());
							 	
							 	//adjust the Vgene name to flag to indel
							 	if (best_VGene_name.indexOf("indel") < 0) {
							 		best_VGene_name = best_VGene_name + " [indels]";
							 	}
							 	
							 }
						 }
				   }	
					 
				}
			
				//create the PostAlignmentResult from the extracted data
				PostAlignmentResult result = new PostAlignmentResult(UMS_string_lowercase, best_alignment_VGene, best_alignment_row_offset, best_alignment_col_offset, best_alignment_gaps, best_alignment_similarity, UMS_name, best_VGene_name, best_align_score, aligned_input);
				//save to the array
				results[resultIndex] = result;
			}
			
			return results;
	
		}  catch (Throwable t) {
			t.printStackTrace();
			return null;
		}
		
		
	}
	
	//parse the BLAST XML file and store the hit results in an array
	private String[][] parseBLASTXMLFile(String blastOutputFile, int resultCount){
		//arrays to store the data from the XML
		String[][] BLASTResultsArray = new String[resultCount][20];
		// changed to get hit-frame rather than query frame String[] featNames = {"num","bit-score","score","evalue","query-from","query-to","hit-from","hit-to","query-frame","identity","positive","align-len","qseq","hseq","midline"};
		String[] featNames = {"num","bit-score","score","evalue","query-from","query-to","hit-from","hit-to","hit-frame","identity","positive","align-len","qseq","hseq","midline"};
		try {
			//read the XML output file
			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			dbf.setValidating(false);
			dbf.setFeature("http://xml.org/sax/features/namespaces",false);
			dbf.setFeature("http://xml.org/sax/features/validation", false);
			dbf.setFeature("http://apache.org/xml/features/nonvalidating/load-dtd-grammar", false);
			dbf.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false);
			DocumentBuilder db = dbf.newDocumentBuilder();
			
			//output file
			Document doc = db.parse(blastOutputFile);
			doc.getDocumentElement().normalize();
			NodeList hitLst = doc.getElementsByTagName("Hit");
			
			for (int s = 0; s < hitLst.getLength(); s++) {
				
				Node hitNode = hitLst.item(s);
				
				if (hitNode.getNodeType() == Node.ELEMENT_NODE) {
					
					//name of germline gene
					Element hitElement = (Element) hitNode;
					//NodeList Hit_def = hitElement.getElementsByTagName("Hit_def");
					NodeList Hit_def = hitElement.getElementsByTagName("Hit_accession"); //changed for blastall as hit_def = 'No definition line found'
					Element e = (Element) Hit_def.item(0);
					NodeList nodeLs = e.getChildNodes();

					
					//get the germline gene name from the hit defintion --> separate so can reformat gene name
					String geneName = ((Node) nodeLs.item(0)).getNodeValue();
					BLASTResultsArray[s][0] = geneName;

					//blastall does not output to XML the Hsp_gaps unless the alignment actually included gaps
					//need to determine gaps separately and set to 0 if no hsp_gaps data can be found
					NodeList hspList = hitElement.getElementsByTagName("Hsp_gaps");
					Element hspElement = (Element) hspList.item(0);
					if (hspElement != null) {
						NodeList hspLst = hspElement.getChildNodes();
						
						BLASTResultsArray[s][1] = ((Node) hspLst.item(0)).getNodeValue();
					} else {
						BLASTResultsArray[s][1] = "0";
					}
					
					//work through the features array and get the results
					for (int i = 0; i < featNames.length; i++) {
						hspList = hitElement.getElementsByTagName("Hsp_" + featNames[i]);
						hspElement = (Element) hspList.item(0);
						NodeList hspLst = hspElement.getChildNodes();
						
						BLASTResultsArray[s][(i+2)] = ((Node) hspLst.item(0)).getNodeValue();  //start at 2 as 0 -> gene 1 -> gaps
					}
					
					//get the ID for the hit so the full length germline gene can be retrieved
					hitElement = (Element) hitNode;
					NodeList Hit_id = hitElement.getElementsByTagName("Hit_id");
					e = (Element) Hit_id.item(0);
					nodeLs = e.getChildNodes();
					String hit_id = ((Node) nodeLs.item(0)).getNodeValue();
					BLASTResultsArray[s][17] = hit_id;
				}
			}
		} catch (Throwable t) {
			t.printStackTrace();
			return null;
		}
		
		return BLASTResultsArray;
	}
	
	
	
	private String getSeqFrmBlastDB(String entry, String blastDB) {
		
		try {
			//var to store retrieved sequence
			String geneSeqStr = "";
			Random rnd = new Random();
			String tmpVSeqFilename = "tmpVseq" + rnd.nextInt() + ".fasta";
			
			//command to search the blast db for the sequence based on hit_id
			//String blastdbcmdCmd = blastLocation + "blastdbcmd" + " -db " + blastDB + " -entry lcl|" + entry + "" + " -outfmt \"%s\"";
			String blastdbcmdCmd = "fastacmd -d " + blastDB + " -s " + entry + " -o " + tmpVSeqFilename; //changed for blastall
		
			
			//print details of command string and db entry search
			//System.out.println("(BestVGeneFinder:getSeqFrmBlastDB) Attempting to get sequence for accession " + entry + " from db " + blastDB + " using command:");
			//System.out.println(blastdbcmdCmd);
			
			//create the blast process, execute blast and wait for the process to finish
			Runtime runtime = Runtime.getRuntime();
			Process blastdbcmd = runtime.exec(blastdbcmdCmd);
			
			//read the output from the process to get the sequence (sequence written to the processes stdout)
			/*InputStream is = blastdbcmd.getInputStream();
			InputStream errStream = blastdbcmd.getErrorStream();
			InputStreamReader isr = new InputStreamReader(is);
			InputStreamReader errorISR = new InputStreamReader(errStream);
			BufferedReader br = new BufferedReader(isr);
			BufferedReader errorBR = new BufferedReader(errorISR);
			String tmpline;
			while ((tmpline = br.readLine()) != null) {
				//System.out.println("br:" + tmpline);
				geneSeqStr = geneSeqStr + tmpline; // Prints all standard output to display
			}
			while ((tmpline = errorBR.readLine()) != null) {
				System.out.println("(BestVGeneFinder:getSeqFromBlastDB) <ERROR>" + tmpline + "</ERROR>");
			}
			
			//close the buffers for reading stdout (inputStream) and stderr (errorStream)
			errorBR.close();
			br.close();*/
			//wait for the command to finish
			int exitVal = blastdbcmd.waitFor();
			//destroy the process now that it is completesh
			if (blastdbcmd != null) {
        		close(blastdbcmd.getOutputStream());
		        close(blastdbcmd.getInputStream());
        		close(blastdbcmd.getErrorStream());
        		blastdbcmd.destroy();
        		blastdbcmd = null;
        		runtime = null;
      		}

			 //outputting to file for blastall so don't need to read from stdout
			
			//reformat the germline gene string
			//geneSeqStr = geneSeqStr.replace("\"", "");
			//geneSeqStr = geneSeqStr.toLowerCase();
			
			//need to read from FASTA file and save the gene seq
			File nf = new File(tmpVSeqFilename);
			//nf.deleteOnExit();
			FastaReader fr = new FastaReader();
			ArrayList<RichSequence> seqList = fr.readFile(nf);
			nf.delete();
			RichSequence tempV = seqList.get(0);
			
			geneSeqStr = tempV.seqString();
			//print out the sequence retrieved from stdout using process getInputStream
			//System.out.println("(BestVGeneFinder:getSeqFrmBlastDB) Exit val: " + exitVal + " \nThe sequence retrieved frm blastdb is: " + geneSeqStr);
			
			//return the retrieved sequence
			return geneSeqStr;
			
		} catch (Throwable t) {
			t.printStackTrace();
			return null;
		}
		
	}
	
		//reverse complement a sequence
    public String reverseComplement(String input) {
        //string to hold the reverse and complemented sequence
        String reverseComplementStr = "";
        
        //make the input all lowercase
        String lcInput = input.toLowerCase();
        
        //get the length of the input
        int inputLength = input.length();
        
        //start from the end of the sequence
        for (int i = inputLength; i > 0; i--) {
            //get the char at position i
            int index = i - 1;
            char currChar = lcInput.charAt(index);
            
            //convert char to complementary nt
            char rcChar;
            if (currChar == 'a') {
                rcChar = 't';
            } else if (currChar == 't') {
                rcChar = 'a';
            } else if (currChar == 'g') {
                rcChar = 'c';
            } else if (currChar == 'c') {
                rcChar = 'g';
            } else if (currChar == 'r') {
                rcChar = 'y';
            } else if (currChar == 'y') {
                rcChar = 'r';
            } else if (currChar == 'k') {
                rcChar = 'm';
            } else if (currChar == 'm') {
                rcChar = 'k';
            } else if (currChar == 's') {
                rcChar = 's';
            } else if (currChar == 'w') {
                rcChar = 'w';
            } else if (currChar == 'b') {
                rcChar = 'v';
            } else if (currChar == 'd') {
                rcChar = 'h';
            } else if (currChar == 'h') {
                rcChar = 'd';
            } else if (currChar == 'v') {
                rcChar = 'b';
            } else if (currChar == 'n') {
                rcChar = 'n';
            } else {
                rcChar = currChar;
            }
            
            //add to the string
            reverseComplementStr += rcChar;
            
        }
        
        //return the reverse complement string
        return reverseComplementStr;
    }

 	private static void close(Closeable c) {
		if (c != null) {
			try {
		    	c.close();
		  	} catch (IOException e) {
		    	// ignored
		  	}
    	}
    }
    
    //creates a 'pseudo' BLAST result that includes the gene names of all the genes with equal top score
	private PostAlignmentResult createPseduoBLASTResult(PostAlignmentResult[] results) {
		//var to track the score
		float topScore = Float.parseFloat("0");
		//string to store appended germline names
		String germlineStr = "";
		
		//necessary to create the new post alignment result
		String UMS_string_lowercase = "NA";
		String best_alignment_VGene = "NA";
		int best_alignment_row_offset = 0;

		int best_alignment_col_offset = 0;
		int best_alignment_gaps = -1;
		int best_alignment_similarity = 0;
		String UMS_name = "NA";
		String best_VGene_name = "NA";
		float best_align_score = 0F;
		char[] aligned_input = null;
		
		for (int i = 0; i < results.length; i++) {
			//get a BLASTResult
			PostAlignmentResult currBLASTResult = results[i];
			
			//get the bit score
			float currBitScore = currBLASTResult.getAlignmentScore();
			//System.out.println("BLASTAlignment:createPseudoBLASTResult) The bit-score is: " + currBitScore);
			
			if (currBitScore == topScore) {
				//current score is equal to top score so just want to append the germline gene name
				if (germlineStr.indexOf("or") > 0) {
					germlineStr = germlineStr + " or " + currBLASTResult.getVGeneName();
				} else {
					germlineStr = germlineStr + " [or " + currBLASTResult.getVGeneName();
				}
				//System.out.println("BestVGeneFinder:createPseudoBLASTResult) Germlinestr (curr == top score) " + germlineStr);
			}
			
			
			if (currBitScore > topScore) {
				//current score is greater than prev topScore
				topScore = currBitScore; //make this the new topScore
				
				//start a fresh string for the appended germline gene names
				germlineStr = currBLASTResult.getVGeneName();
				//System.out.println("BestVGeneFinder:createPseudoBLASTResult) Germlinestr (curr > top score) " + germlineStr);
				
				//this is the first of the genes in the list, therefore use all its features for the HMM and displaying the alignment
				UMS_string_lowercase = currBLASTResult.getFullLengthUMSString();
				best_alignment_VGene = currBLASTResult.getFullLengthVGeneString();
				best_alignment_row_offset = currBLASTResult.getRowOffset();
				best_alignment_col_offset = currBLASTResult.getColOffset();
				best_alignment_gaps = currBLASTResult.getAlignmentGaps();
				best_alignment_similarity = currBLASTResult.getSimilarity();
				UMS_name = currBLASTResult.getUMSname();
				best_align_score = currBitScore;
				aligned_input = currBLASTResult.getAlignedInput();
			}
			
		}
		
		//check if a closing bracket needs to be added to the germline gene line
		if (germlineStr.indexOf("or") > 0) {
			germlineStr = germlineStr + "]";
		}
		//set the newly created gene name string to the Vgene name
		best_VGene_name = germlineStr;
	
					
		//create the PostAlignmentResult from the extracted data for the 'first' Vgene plus the newly created pseudo Vgene name
		PostAlignmentResult result = new PostAlignmentResult(UMS_string_lowercase, best_alignment_VGene, best_alignment_row_offset, best_alignment_col_offset, best_alignment_gaps, best_alignment_similarity, UMS_name, best_VGene_name, best_align_score, aligned_input);
			
		//return the pseudo BLAST result
		return result;
	}
    
}
