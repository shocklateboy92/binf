package iHMMuneAlign;

//import gnu.bioinformatics.jaligner.JAligner;
//import gnu.bioinformatics.jaligner.util.Alignment;
import java.io.*;
import java.util.ArrayList;
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
import java.util.Random;


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
//            FastaReader

public class TrailingJGeneVectorFinder
{

    final String DEFAULT_JGENE_FILE_NAME = "j-regions.fa";
    final float OPEN_GAP_COST = 100F;
    final float EXTEND_GAP_COST = 5F;
	final int GAPS_ALLOWED = 1;
    final int GAPS_ERROR = 2;
    final int GAPS_IGNORE_SEQUENCE = 3;
    int Gap_Behaviour;
    String SCORE_MATRIX_NAME;
	String geneRepertoireName = "";
	String blastPath = "";
	boolean DEBUGGING = false;

    public static void main(String args[])
    {
       
    }

    public TrailingJGeneVectorFinder(File fastaFile, String blastPath)
    {
        Gap_Behaviour = GAPS_ALLOWED; //GAPS_ERROR;
		this.blastPath = blastPath;
		
		try {
			//get the name from the fastaFile
			String fastaFileName = fastaFile.getAbsolutePath();
			//build name of the blast db file
			File dbFastaFileName = new File (fastaFileName + ".nhr");
			
			//determine if a pre-existing db file is present
			if (dbFastaFileName.exists()) {
				//found pre-existing
				if (DEBUGGING) {
					System.out.println("(TrailingJGeneVectorFinder) Pre-existing gene repertoire blastdb file found: " + dbFastaFileName);
				}
				dbFastaFileName = null;
			} else {
				//no pre-existing found so need to create
				
				//create the command string to be executed
				//String cmdStrCreateDB = blastPath + "makeblastdb -in " + fastaFileName + " -dbtype nucl -parse_seqids";
				String cmdStrCreateDB = "formatdb -i " + fastaFileName + " -n " + fastaFileName + " -p F -o T";  //changed to blastall
				if (DEBUGGING) {
					System.out.println("(TrailingJGeneVectorFinder) Executing command to create BLAST db from repertoire file: " + cmdStrCreateDB);
				}
				//execute the command to create the db files
				Runtime rt = Runtime.getRuntime();
				Process makeblastdb = rt.exec(cmdStrCreateDB);
				makeblastdb.waitFor();
				
				//process completed so destroy
				if (makeblastdb != null) {
		    		close(makeblastdb.getOutputStream());
				    close(makeblastdb.getInputStream());
		    		close(makeblastdb.getErrorStream());
		    		makeblastdb.destroy();
		    		makeblastdb = null;
		    		rt = null;
		    		dbFastaFileName = null;
      			}
      			
				if (DEBUGGING) {
					System.out.println("(TrailingJGeneVectorFinder) The repertoire file was made into a BLAST database file for searching: " + fastaFileName);
				}
			}
			
			this.geneRepertoireName = fastaFileName;
			
		}
		catch (Throwable t) {
			t.printStackTrace();
		}
        
    } //end TrailingJGeneVectorFinder(File)

    public TrailingJGeneVectorFinder()
    {
        
    } //end TrailingJGeneVectorFinder()
    
    public PostAlignmentResult getAlignmentResult(String inputSeq) {
    	//JAligner jaligner = new JAligner();
        String input_string_uppercase = inputSeq.toUpperCase();
        int input_LENGTH = input_string_uppercase.length();
        if (DEBUGGING) {
        	System.out.println("LENGTH of input string: " + input_LENGTH);
        }
        String best_alignment_JGene = null;
        String best_JGene_name = null;
        PostAlignmentResult alignmentResult = null;
        float best_score = 0F;

		try
		{
			alignmentResult = getAlignment(inputSeq, geneRepertoireName);
		}
		catch(Exception ex)
		{
			throw new Error("(TrailingJGeneVectorFinder) do alignment failed " + ex.getMessage());
		} //end try/catch
		
		return alignmentResult;
    }

    public String getResult(PostAlignmentResult curr_alignment, int MIN_UMS_ALIGNMENT_START_OFFSET)
    {
        String inputSeq = curr_alignment.getFullLengthUMSString();
        String input_string_uppercase = inputSeq.toUpperCase();
        int input_LENGTH = input_string_uppercase.length();
        
         		
		//check for gaps in the returned alignment to the IGHJ repertoire
		int best_alignment_gaps = curr_alignment.getAlignmentGaps();
        if ((best_alignment_gaps > 0) && (Gap_Behaviour == GAPS_ERROR))
        {
            //throw new Error("best JGene alignment has gaps");
            String tempReturnStr = "GAP";
            return tempReturnStr;
        } else if ((best_alignment_gaps > 0) && (Gap_Behaviour == GAPS_ALLOWED)) {
        	//there was a gap found in the best matching IGHJ, but the current gap bahaviour is allowing this so need to create a 'fake' entry in the J repertoire file
        	//the JRepertoire has already been loaded in the AlignmentThread class into an ArrayList <RichSequence>
        	//need to add an additional richsequence to the ArrayList
        }
		
		//extract the alignment results out of the alignment result object
		//similarty between input and best matching IGHJ
        int best_alignment_similarity = curr_alignment.getSimilarity();
		//offsets
        int best_alignment_row_offset = curr_alignment.getRowOffset(); //input sequence offset
        int best_alignment_col_offset = curr_alignment.getColOffset(); //germline repertoire offset
		//length of portion of input sequence that matches against the IGHJ gene
		String best_alignment_match = new String(curr_alignment.getAlignedInput());
		int jLength = best_alignment_match.length();
		int matchingJGenePartLength = jLength;
		//total length of input sequence that is upstream of last nt of IGHJ
        int endJGeneMatchPosition = best_alignment_row_offset + matchingJGenePartLength;
		//alignment score
        float score = curr_alignment.getAlignmentScore();
		
		//print out the details of input sequence alignment against IGHJ repertoire
        if (DEBUGGING) {
        	System.out.println("score = " + score);
        	System.out.println("row offset (input sequence) = " + best_alignment_row_offset);
        	System.out.println("col offset (JGene repertoire) = " + best_alignment_col_offset);
        	System.out.println("minimum J length required = " + MIN_UMS_ALIGNMENT_START_OFFSET);
        	System.out.println("JGene Length = " + jLength);
        	System.out.println("position at which input stops matching to the JGene = " + endJGeneMatchPosition);
        	System.out.println("length of the input sequence = " + input_LENGTH);
        }
         
        //check that the J match found is long enough
        if (jLength < MIN_UMS_ALIGNMENT_START_OFFSET) {
          String tempReturnStr = "NA";
          return tempReturnStr;

        } else {
          //extract the string the corresponds to the portion of the input seq that is upstream of end of IGHJ match
		    String UMSminusTrailingJGenePart;
		    		   
		    if(endJGeneMatchPosition < input_LENGTH) {
		        endJGeneMatchPosition = endJGeneMatchPosition - 1; //need reindex from 0 rather than 1 for use in finding the substr below, else get an extra nt in subs
		        //the input sequence extends beyond the final position of the J alignment, so need to remove the downstream by taking substring up to end of alignment to J
		        UMSminusTrailingJGenePart = inputSeq.substring(0, endJGeneMatchPosition);
		        if (DEBUGGING) {
		        	System.out.println("Trailing JGene Vector or C-region removal");
		        }
		        
		    } else {
		    	//the input sequence is either the same length as the final match to the J or the J finishes after the end of the input sequence
		        UMSminusTrailingJGenePart = inputSeq;
		        if (DEBUGGING) {
		        	System.out.println("NO Trailing JGene Vector or C-region removal");
		        }
		    }
		    
		    if (DEBUGGING) {
		    	System.out.println("UMS minus trailing JGene part = \n" + UMSminusTrailingJGenePart);
				System.out.println("Length of input sequence less trailing J gene part = " + UMSminusTrailingJGenePart.length());
			}
		
		    UMSminusTrailingJGenePart = UMSminusTrailingJGenePart.toLowerCase();
		    return UMSminusTrailingJGenePart;
        
        }
    }//end getResult

    private PostAlignmentResult getAlignment(String inputSeq, String dbName)
    {
		try {
			PostAlignmentResult result = doAlignment(inputSeq, dbName);
			return result;
		}
		catch (Exception ex) {
			throw new Error ("getAlignment failed: " + ex.getMessage());
		}
    } //end getAlignment

    private PostAlignmentResult doAlignment(String inputSeq, String dbName) throws Exception
    {
		try {
			//location for the temp file that holds the input sequence
            Random rnd = new Random();
			String tmpFileName = "temp_input"+rnd.nextInt();
			FileOutputStream tempFile = new FileOutputStream(tmpFileName + ".input");
			File tempFileFile = new File(tmpFileName + ".input");
			//tempFileFile.deleteOnExit();
			File tempFileOutput = new File(tmpFileName + ".output");
			//tempFileOutput.deleteOnExit();

			//grab the default namespace
			Namespace ns = RichObjectFactory.getDefaultNamespace();
			//create a temp id & a rich sequence for the input <- only required if passing strings
			String tmpID = "input_query_seq";
			RichSequence rs = org.biojavax.bio.seq.RichSequence.Tools.createRichSequence(ns, tmpID, inputSeq, DNATools.getDNA());

			//write the sequence object to the temp file
			RichSequenceIterator it = new SingleRichSeqIterator(rs);
			RichSequenceFormat fasta = new FastaFormat();
			RichStreamWriter seqOut = new RichStreamWriter(tempFile, fasta);
			seqOut.writeStream(it,ns);

			//sequence has been written so can close the output stream
			if (tempFile != null) {
				tempFile.close();
				tempFile = null;
			}
			
			if (DEBUGGING) {
				System.out.println("(TrailingJGeneVectorFinder) input sequence written to file ready for BLASTN: " + tmpFileName);
			}
			
			String blastGapStr = "";
			switch (Gap_Behaviour) {
				case 1:
					//blastGapStr = "-gapopen -gapextend"
					blastGapStr = " -G 7 -E 4";
					//specify gap penalties above
					break;
				case 2:
					//should throw error if there is a gap
					//keep default gap behaviour in blast and then throw error if gap in best alignment when processing results
					blastGapStr = " -G 7 -E 4";
					break;
				case 3:
					//ignore alignments that contains gaps
					//blastGapStr = " -ungapped";
					blastGapStr = " -g F";  //changed to blastall args
					break;
			}
			//blast scoring for IGHJ alignments
			//String rewardStr = " -reward 2";
			String rewardStr = " -r 2"; //changed to blastall args
			//String penaltyStr = " -penalty -5";
			String penaltyStr = " -q -5"; //changed to blastall args
			//String wordSizeStr = " -wordsize 7";
			String wordSizeStr = " -W 7"; //changed to blastall args
			
			//create the command string to execute blast
			//String blastCmdStr = blastPath + "blastn -query "  + tmpFileName + ".input -out " + tmpFileName + ".output" + " -db " + dbName + " -num_alignments 5 -num_descriptions 5 -outfmt 5" + blastGapStr + rewardStr + penaltyStr + wordSizeStr;
			String blastCmdStr = "blastall -p blastn -i " + tmpFileName + ".input -o " + tmpFileName + ".output -d " + dbName + " -b 5 -v 5 -m 7 " + blastGapStr + rewardStr + penaltyStr + wordSizeStr;

			if (DEBUGGING) {
				System.out.println("(TrailingJGeneVectorFinder) BLAST command to be executed: " + blastCmdStr);
			}

			//create the blast process, execute blast and wait for the process to finish
			Runtime rt = Runtime.getRuntime();
			Process blastn = rt.exec(blastCmdStr);
			blastn.waitFor();

			int resultNo = 0;
			PostAlignmentResult result = getBLASTResults(tmpFileName + ".output", resultNo, inputSeq);

			//destroy the blastn process
			if (blastn != null) {
        		close(blastn.getOutputStream());
		        close(blastn.getInputStream());
        		close(blastn.getErrorStream());
        		blastn.destroy();
        		blastn = null;
        		rt = null;
      		}
      		
      		//destroy the temp files used during the blast alignment
      		if (tempFileFile != null) {
      			tempFileFile.delete();
      			tempFileFile = null;
      		}
      		if (tempFileOutput != null) {
	      		tempFileOutput.delete();
	      		tempFileOutput = null;
	      	}
      		
			return result;

		} catch (Throwable t) {
			t.printStackTrace();
			return null;
		}
	}//end doAlignment

	//parse through the XML BLAST output file created and extract the details of alignment of input sequence against the germline repertoire
	private PostAlignmentResult getBLASTResults (String blastOutputFile, int resultNo, String inputSeq) {
		try {
			//arrays to store the data from the XML
			String[][] jGeneBLASTResults = new String[5][20];
			String[] featNames = {"num","bit-score","score","evalue","query-from","query-to","hit-from","hit-to","query-frame","identity","positive","align-len","qseq","hseq","midline"};

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
					Element hitElement = (Element) hitNode;
					NodeList Hit_def = hitElement.getElementsByTagName("Hit_accession");
					Element e = (Element) Hit_def.item(0);
					NodeList nodeLs = e.getChildNodes();

					//get the IGHJ name from the hit defintion
					String geneName = ((Node) nodeLs.item(0)).getNodeValue();
					jGeneBLASTResults[s][0] = geneName;
					if (DEBUGGING) {
						System.out.println("output : " + geneName);
					}

					//blastall does not output to XML the Hsp_gaps unless the alignment actually included gaps
					//need to determine gaps separately and set to 0 if no hsp_gaps data can be found
					NodeList hspList = hitElement.getElementsByTagName("Hsp_gaps");
					Element hspElement = (Element) hspList.item(0);
					if (hspElement != null) {
						NodeList hspLst = hspElement.getChildNodes();
						jGeneBLASTResults[s][1] = ((Node) hspLst.item(0)).getNodeValue();
					} else {
						jGeneBLASTResults[s][1] = "0";
					}

					for (int i = 0; i < featNames.length; i++) {
						hspList = hitElement.getElementsByTagName("Hsp_" + featNames[i]);
						hspElement = (Element) hspList.item(0);
						NodeList hspLst = hspElement.getChildNodes();

						jGeneBLASTResults[s][(i+2)] = ((Node) hspLst.item(0)).getNodeValue();
					}//end for loop
				}//end if
			}//end for
			
			//get the 'best' IGHV alignment
			//default values
			String best_Gene_name = "NA - no BLAST hit found";
			String best_alignment_Gene = "NA";
			//input string
			String UMS_string_lowercase = inputSeq.toLowerCase();
			String UMS_name = "input_seq";
			//offsets
			int best_alignment_row_offset = 0;
			int best_alignment_col_offset = 0;
			//gaps
			int best_alignment_gaps = 0;
			//similarity
			int best_alignment_similarity = 0;
			//alignment score
			float best_align_score = 0F;
			//aligned input
			String alignedStr = "";
			char[] aligned_input = alignedStr.toCharArray();
			
			if (jGeneBLASTResults[0][0] != null) {
				//adjust result number required
				int resultIndex = resultNo;
				
				//IGHV name & seq
				best_Gene_name = jGeneBLASTResults[resultIndex][0];
				best_alignment_Gene = jGeneBLASTResults[resultIndex][15];
				//input string
				//UMS_string_lowercase = (inputSeq.seqString()).toLowerCase();
				//UMS_name = seqString.getName();
				//offsets
				best_alignment_row_offset = Integer.parseInt(jGeneBLASTResults[resultIndex][6]);
				best_alignment_col_offset = Integer.parseInt(jGeneBLASTResults[resultIndex][8]);
				//gaps
				best_alignment_gaps = Integer.parseInt(jGeneBLASTResults[resultIndex][1]);
				//similarity
				best_alignment_similarity = Integer.parseInt(jGeneBLASTResults[resultIndex][12]);
				//alignment score
				best_align_score = Float.parseFloat(jGeneBLASTResults[resultIndex][4]);
				//aligned input
				alignedStr = jGeneBLASTResults[resultIndex][14];
				aligned_input = alignedStr.toCharArray();
				
				best_alignment_Gene = getSeqFrmBlastDB(best_Gene_name, geneRepertoireName);

				 if ((best_alignment_gaps > 0) && (Gap_Behaviour == GAPS_ALLOWED)) {
						 	//allowing gaps so need to make changes to the full length germline gene that is to be used for the PostAlignmentResult
						 	//need to if there has been a deletion relative to germline need to remove to pad the input seq or if there has been an 
						 	//insertion relative to the germline need to pad the germline
							
							//get the aligned V that will be missing dels and have '-' for ins in input
						 	String jGeneHit = jGeneBLASTResults[0][15];
							
							//gap in the germline gene string
							if (jGeneHit.indexOf("-") > 0) {
								//get the full length sequence of the v gene
							 	String fullLengthJGene = best_alignment_Gene;
							 	
							 	//get the index for the start and end of the hit found to the vgene
							 	int JgeneFromIndex = Integer.parseInt(jGeneBLASTResults[0][8]);
							 	int JgeneToIndex = Integer.parseInt(jGeneBLASTResults[0][9]);
							 	
							 	//determine what needs to be added up and downstream of the aligned portion of the V to make the full length seq that includes indels
							 	String prependForFullLength = fullLengthJGene.substring(0, (JgeneFromIndex-1));
								String postpendForFullLength = "";
							 	if (JgeneToIndex < fullLengthJGene.length()) {
							 		postpendForFullLength = fullLengthJGene.substring(JgeneToIndex, fullLengthJGene.length());
							 	}
						 	
							 	//add the up and downstream seq to create the full length V
							 	String jGeneWithIndels = prependForFullLength + jGeneHit + postpendForFullLength;
							 	
							 	//set the newly created sequence to the J sequence
							 	best_alignment_Gene = jGeneWithIndels;
							 	
							 	//adjust the name
							 	best_Gene_name = best_Gene_name + " [indels]";
							}
						 	
						 	//gap in the input string
						 	if (alignedStr.indexOf("-") > 0) {
							 	//get the aligned input sequence that will have '-' if there has been a deletion from the germline
							 	String fullLengthInput = UMS_string_lowercase;
							 	int inputFromIndex = Integer.parseInt(jGeneBLASTResults[0][6]);
							 	int inputToIndex = Integer.parseInt(jGeneBLASTResults[0][7]);
							 	String prependForInput = fullLengthInput.substring(0, (inputFromIndex-1));
							 	
							 	String postpendForInput = "";
							 	if (inputToIndex < fullLengthInput.length()) {
							 		postpendForInput = fullLengthInput.substring(inputToIndex, fullLengthInput.length());
							 	}
							 	String inputWithIndels = prependForInput + alignedStr + postpendForInput;
							 	
							 	UMS_string_lowercase = (inputWithIndels.toLowerCase());
							}
				}
			} 


			//create the PostAlignmentResult from the extracted data
			PostAlignmentResult result = new PostAlignmentResult(UMS_string_lowercase, best_alignment_Gene, best_alignment_row_offset, best_alignment_col_offset, best_alignment_gaps, best_alignment_similarity, UMS_name, best_Gene_name, best_align_score, aligned_input);

			return result;

		}  catch (Throwable t) {
			t.printStackTrace();
			return null;
		}

	}//end getBLASTResults

 	private static void close(Closeable c) {
		if (c != null) {
			try {
		    	c.close();
		  	} catch (IOException e) {
		    	// ignored
		  	}
    	}
    }
    
    private String getSeqFrmBlastDB(String entry, String blastDB) {
		
		try {
			//var to store retrieved sequence
			String geneSeqStr = "";
			Random rnd = new Random();
			String tmpJSeqFilename = "tmpJseq" + rnd.nextInt() + ".fasta";
			
			//command to search the blast db for the sequence based on hit_id
			//String blastdbcmdCmd = blastLocation + "blastdbcmd" + " -db " + blastDB + " -entry lcl|" + entry + "" + " -outfmt \"%s\"";
			String blastdbcmdCmd = "fastacmd -d " + blastDB + " -s " + entry + " -o " + tmpJSeqFilename; //changed for blastall
		
			
			//print details of command string and db entry search
			if (DEBUGGING) {
				System.out.println("(BestVGeneFinder:getSeqFrmBlastDB) Attempting to get sequence for accession " + entry + " from db " + blastDB + " using command:");
				System.out.println(blastdbcmdCmd);
			}
			
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
			File nf = new File(tmpJSeqFilename);
			//nf.deleteOnExit();
			FastaReader fr = new FastaReader();
			ArrayList<RichSequence> seqList = fr.readFile(nf);
			nf.delete();
			RichSequence tempV = seqList.get(0);
			
			geneSeqStr = tempV.seqString();
			//print out the sequence retrieved from stdout using process getInputStream
			if (DEBUGGING) {
				System.out.println("(BestVGeneFinder:getSeqFrmBlastDB) Exit val: " + exitVal + " \nThe sequence retrieved frm blastdb is: " + geneSeqStr);
			}
			
			//return the retrieved sequence
			return geneSeqStr;
			
		} catch (Throwable t) {
			t.printStackTrace();
			return null;
		}
		
	}        

    
}//end TrailingJGeneVectorFinder
