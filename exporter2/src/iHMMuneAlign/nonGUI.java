package iHMMuneAlign;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Vector;
import org.biojavax.bio.seq.RichSequence;
import java.util.Random;

// Referenced classes of package iHMMuneAlign:
//            GlobalDefines, MutabilityScoreTrinucleotide, MutabilityScoreHotspot, ProbabilityHolder,
//            FastaReader, OutputWindow, AlignmentThread, MutabilityScore

public class nonGUI
    implements GlobalDefines
{
    ArrayList main_sequence_list;
    File IGHV_FILE;
    File IGHD_FILE;
    File IGHJ_FILE;
    File IGKV_FILE;
    File IGKJ_FILE;
    File IGLV_FILE; //added to allow alignment of IGL sequences
    File IGLJ_FILE; //added to allow alignment of IGL sequences
    int specifiedChainType;
    char html_writer_type;
    int d_gene_acceptance_type;
    MutabilityScoreTrinucleotide mutabilityScoreTrinucleotide;
    MutabilityScoreHotspot mutabilityScoreHotspot;
    MutabilityScore mutabilityScore;
    PrintWriter pw;
    PrintWriter out;
    boolean DEBUGGING = false;

    public nonGUI(String sequenceSetFileName, String IGHV_Filename, String IGHD_Filename, String IGHJ_Filename, String IGKV_Filename, String IGKJ_Filename, int userSpecifiedChainType, String probabilityFilename, String htmlFileName, int scoringSystem, int htmlWriterType, int dAcceptanceLevel, byte alignment_type, int startIndex, int endIndex, String scoringMatrix, int minJLength, int alignmentModelType, int numOfVgenes, int numSuboptimal, String blastLocation, String IGLV_Filename, String IGLJ_Filename, String mutationSpecFileName, String dProbsFilename)
    {
        if (DEBUGGING) {
	        System.out.println("alignment model type from the commmand line " + alignmentModelType);
	    }

        main_sequence_list = null;

        //names of the repertoire files -> get these from the args list
        IGHV_FILE = new File(IGHV_Filename);
        IGHD_FILE = new File(IGHD_Filename);
        IGHJ_FILE = new File(IGHJ_Filename);
        IGKV_FILE = new File(IGKV_Filename);
        IGKJ_FILE = new File(IGKJ_Filename);
        //additional files for lambda
        IGLV_FILE = new File(IGLV_Filename);
        IGLJ_FILE = new File(IGLJ_Filename);
        specifiedChainType = userSpecifiedChainType;

        mutabilityScoreTrinucleotide = null;
        mutabilityScoreHotspot = null;
        mutabilityScore = null;
        pw = null;

        //determine the output type
        if (htmlWriterType == 1) {
          html_writer_type = '\001'; //excel based output
          //use html extension
          htmlFileName = htmlFileName + ".html";
        } else if (htmlWriterType == 2) {
          html_writer_type = '\002';  //web based output
          //use html extension for output file
          htmlFileName = htmlFileName + ".html";
        } else {
          html_writer_type = '\003';
          //change the output file to the text extension rahter than html
          htmlFileName = htmlFileName + ".txt";
        }

        //read in the sequences from the specified file
        ArrayList sequenceList = loadSequenceSequenceSet(sequenceSetFileName);
		setSequenceList(sequenceList);
		
        //load the probabilities from the file -> get the file from the args list
        loadProbabilities(probabilityFilename);
        
		//load the D & J repertoires
		FastaReader fr = new FastaReader();
		//ArrayList<RichSequence> dSequences = fr.readFile(IGHD_FILE);
		//ArrayList<RichSequence> jSequences = fr.readFile(IGHJ_FILE);
		if (DEBUGGING) {
			System.out.println("(nonGUI) D and J repertoires loaded");
		}
        
		//create and set the scoring systems
        mutabilityScoreTrinucleotide = new MutabilityScoreTrinucleotide();
        mutabilityScoreHotspot = new MutabilityScoreHotspot();
		if (scoringSystem == 0) {
		   mutabilityScore = mutabilityScoreTrinucleotide;    //scoring based on trinucleotides
		}
		else {
		   mutabilityScore = mutabilityScoreHotspot;          //scoring based on hotspots
		}
		if (DEBUGGING) {
			System.out.println("(nonGUI) mutability scoring loaded");
		}
		
		//set the IGHD gene acceptance behaviour
		if (dAcceptanceLevel == 5) {
		   d_gene_acceptance_type = 2;     //5mer required
		}
		else if (dAcceptanceLevel == 8) {
		   d_gene_acceptance_type = 1;    //8mer required
		} 
		else if(dAcceptanceLevel == 0) {   //dynamic setting of the d acceptance
		   d_gene_acceptance_type = 3;
		}
		if (DEBUGGING) {
			System.out.println("(nonGUI) d acceptance set");
		}

		//create file for storing details of mutations in alignments
		try
		{
			String mutOutputFilename = null;
			Random rnd = new Random();
			mutOutputFilename = "mutations_in_alignments" + rnd.nextInt() + ".txt";
			File mutOutputFile = new File(mutOutputFilename);
			mutOutputFile.deleteOnExit();
			
			pw = new PrintWriter(new FileWriter(mutOutputFilename));
		}
		catch(Exception ex)
		{
			throw new Error(ex.getMessage());
		}
		if (DEBUGGING) {
			System.out.println("(nonGUI) mutations in alignments file created");
		}
            
		//selectedIndices is a array of integers that lists the indices for all input sequences for processing
		//int indexStart = startIndex;
		//int indexEnd = endIndex;
		int numberOfAlignments = ((endIndex - startIndex) + 1);
		if (DEBUGGING) {
			System.out.println("(nonGUI) startIndex: " + startIndex + " end index: " + endIndex);
			System.out.println("(nonGUI) number of alignment set: " + numberOfAlignments);
		}
		if (numberOfAlignments > 99) {
		  //don't want to be performing more than 100 alignments in one go unless array size is changed for selectedIndicies
		  //so set the indexEnd back to something smaller
		  endIndex = startIndex + 99;
		  numberOfAlignments = 100;
		}
		if (DEBUGGING) {
			System.out.println("(nonGUI) number of alignment set: " + numberOfAlignments);
		}
		
		/*
			  
		int currentIndex = indexStart;
		System.out.println("current index set to " + indexStart);
		int[] selectedIndices = new int[100];  //change to get indices from the args list
		for (int k=0; k < numberOfAlignments; k++) {
		  selectedIndices[k] = currentIndex;
		  currentIndex++;
		  System.out.println("current index in making array = " + currentIndex);
		}*/

		//change to output this stuff to file rather then to the GUI window
		//OutputWindow tempOutputWindow = new OutputWindow("Multiple Sequence Alignment Result");
		//JTextArea outputArea = tempOutputWindow.getTextArea();

		//System.out.println("selected number of sequences = " + selectedIndices.length);
		for(int i = startIndex; i <= endIndex; i++)
		{
			//int currSelectionIndex = selectedIndices[i];
			//RichSequence currSequence = (RichSequence)main_sequence_list.get(currSelectionIndex);
			//get sequence from list
			if (DEBUGGING) {
				System.out.println("(nonGUI) About to get sequence: " + (i+1));
			}
			RichSequence currSequence = (RichSequence)sequenceList.get(i);
			if (DEBUGGING) {
				System.out.println("aligning number i = " + (i+1));
			}
			//send current sequence for partitioning
			alignSequence(currSequence, htmlFileName, specifiedChainType, alignment_type, d_gene_acceptance_type, probabilityFilename, scoringMatrix, minJLength, alignmentModelType, numOfVgenes, numSuboptimal, blastLocation, mutationSpecFileName, dProbsFilename);
			if (DEBUGGING) {
				System.out.println("finished aligning number i = " + (i+1));
			}
		}

		if(pw != null)
		{
			pw.close();
		}
		
		//out.close();
		// end of the mulitple aligment code
    }
	
	//create vectors from input sequence arraylist
    private void setSequenceList(ArrayList sequenceList)
    {
        Vector ums_vector = new Vector();
        Vector ums_name_vector = new Vector();
        for(int i = 0; i < sequenceList.size(); i++)
        {
            RichSequence tempSeq = (RichSequence)sequenceList.get(i);
            String indexed_sequence_name = (i + 1) + "\t: " + tempSeq.getName();
            ums_name_vector.add(indexed_sequence_name);
            ums_vector.add(tempSeq.seqString());
        }

    }


    private void loadProbabilities(String openFileName)
    {
        File openFile = new File(openFileName);
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(openFile));
            String input = "";
            for(String line = br.readLine(); line != null; line = br.readLine())
            {
                input = input + line + " ";
            }
			
			if (DEBUGGING) {
	            System.out.println(input);
	        }
            ProbabilityHolder newProbabilityHolder = ProbabilityHolder.createFromString(input);
        }
        catch(IOException ioe)
        {
            throw new Error("loadProbsButton failed:  " + ioe.getMessage());
        }
    }

    private ArrayList loadSequenceSequenceSet(String openFileName)
    {
        if (DEBUGGING) {
	        System.out.println("loadSequenceSequenceSet function called with fileName: " + openFileName);
	    }
        File openFile = new File(openFileName);
        ArrayList sequenceList = null;
        FastaReader fr = new FastaReader();
        sequenceList = fr.readFile(openFile);
        if(sequenceList == null)
        {
            System.out.println("Error: Could not read Sequence List FASTA file: " + openFile.getName());
            return null;
        } else
        {
            main_sequence_list = sequenceList;
            return sequenceList;
        }
    }

    private void alignSequence(RichSequence sequence, String htmlFileName, int specifiedChainType, byte alignmentType, int dGeneAcceptanceType, String probabilityFilename, String scoringMatrix, int minJLength, int alignmentModelType, int numOfVgenes, int numSuboptimal, String blastLocation, String mutationSpecFileName, String dProbsFilename)
    {
        if(alignmentType != 4)
        {
            String sequenceName = sequence.getName();
            htmlFileName = htmlFileName + " - " + sequenceName;
        }
        //ProbabilityHolder probHolder = ;
        File openProbFile = new File(probabilityFilename);
        try
        {
            BufferedReader br = new BufferedReader(new FileReader(openProbFile));
            String input = "";
            for(String line = br.readLine(); line != null; line = br.readLine())
            {
                input = input + line + " ";
            }
			
			//close the buffered reader
			br.close();

            ProbabilityHolder probHolder = ProbabilityHolder.createFromString(input);
            //create the alignment thread
			AlignmentThread alignmentThread = new AlignmentThread(sequence, IGHV_FILE, IGHD_FILE, IGHJ_FILE, IGKV_FILE, IGKJ_FILE, specifiedChainType, htmlFileName, probHolder, mutabilityScore, html_writer_type, alignmentType, dGeneAcceptanceType, scoringMatrix, minJLength, alignmentModelType, numOfVgenes, numSuboptimal, blastLocation, IGLV_FILE, IGLJ_FILE, mutationSpecFileName, dProbsFilename);
			
			//start the alignment thread
			alignmentThread.start();
			
			//if outputting to single file waits for current thread to end before starting new thread
			if(alignmentType == 4) {
				//if(!alignmentThread.isAlive()); //check if thread alive
				if (alignmentThread.isAlive()) {
					try
					{
						if (DEBUGGING) {
							System.out.println("(nonGUI:alignSequence) Waiting for alignment thread to complete: " + alignmentThread.getName());
							System.out.println("(nonGUI:alignSequence) The number of currently acive threads is: " + alignmentThread.activeCount());
						}
						alignmentThread.join();  //wait for thread to die
						//once thread dead make null
						if (DEBUGGING) {
							System.out.println("(nonGUI:alignSequence) The thread is complete.");
						}
						//attempt to invoke garbage collection to deal with mem usages
						System.gc();
						if (DEBUGGING) {
							System.out.println("HINTING TO Garbage Collection");
						}
					}
					catch(InterruptedException ie)
					{
						throw new Error("alignment thread interrupted:  " + ie.getMessage());
					}
				}
			}
			
			//clean-up
			alignmentThread = null;
			probHolder = null;

        }
        catch(IOException ioe)
        {
            throw new Error("loadProbsButton failed:  " + ioe.getMessage());
        }

     }

    public static void main(String args[])
    {
        String startupSequencesFileName;
        String IGHV_Filename;
        String IGHD_Filename;
        String IGHJ_Filename;
        String IGKV_Filename;
        String IGKJ_Filename;
        String IGLV_Filename;
        String IGLJ_Filename;
        int userSpecifiedChainType;
        String probabilityFilename;
        String htmlFilename;
        String txtFilename;
        int scoringSystem;
        int htmlWriterType;
        int dAcceptanceLevel;
        byte alignment_type;
        int startIndex;
        int endIndex;
        String scoringMatrix;
        int minJLength;
        int alignmentModelType;
        int numOfVgenes;
        int numSuboptimal;
	String blastLocation;
        String mutationSpecFileName;
        String dProbsFilename;

        try {
           startupSequencesFileName = args[0];
           IGHV_Filename = args[1];
           IGHD_Filename = args[2];
           IGHJ_Filename = args[3];
           IGKV_Filename = args[4];
           IGKJ_Filename = args[5];
           userSpecifiedChainType = Integer.parseInt(args[6]);
           probabilityFilename = args[7];
           htmlFilename = args[8];
           scoringSystem = Integer.parseInt(args[9]);
           htmlWriterType = Integer.parseInt(args[10]);
           dAcceptanceLevel = Integer.parseInt(args[11]);
           alignment_type = Byte.parseByte(args[12]);
           startIndex = Integer.parseInt(args[13]);
           endIndex = Integer.parseInt(args[14]);
           scoringMatrix = args[15];
           minJLength = Integer.parseInt(args[16]);
           alignmentModelType = Integer.parseInt(args[17]);
           numOfVgenes = Integer.parseInt(args[18]);
           numSuboptimal = Integer.parseInt(args[19]); // the numSuboptimal is actually the total number of ALTERNATIVE alignments (ie. including the optimal)
	   blastLocation = args[20];
	   //additional args to allow the IGL files
	   IGLV_Filename = args[21];
	   IGLJ_Filename = args[22];
			
           //adjust these as the indexing of the list of sequences start at 0 rather than 1
           startIndex = startIndex -1;
           endIndex = endIndex -1;

           //added var for passing the location of the mutation spectrum file
           mutationSpecFileName = args[23];
           //added var for passing the location of the file containing the probabilties for N-additions being mis-identified as D-REGION nucleotides
           dProbsFilename = args[24];

	} catch (Exception ex) {
	   	   throw new Error("The correct number of parameters were not provided");
	   	   /*startupSequencesFileName = "Andrews Sequences All - altered.txt";
           //filenames for the repertoire files
           IGHV_Filename = "IGHV_Repertoire.fa";
           IGHD_Filename = "IGHD_Repertoire.fa";
           IGHJ_Filename = "IGHJ_Repertoire.fa";
           IGKV_Filename = "IGKV_Repertoire.fa";
           IGKJ_Filename = "IGKJ_Repertoire.fa";
           userSpecifiedChainType = 0;
           probabilityFilename = "Current iHMMune-align Probabilities.PH";
           htmlFilename = "test outputFile";

           //type of scoring required
           scoringSystem = 1; //0 for trinucleotide, 1 for HS
           //the type of output required -> get from the args list
           htmlWriterType = 1;

           //the level of IGHD gene acceptance
           dAcceptanceLevel = 5; //for fivermer can set to 8 for 8mers, defaults to 8mer if nothing supplied

           //single or multiplefile output
           //byte alignment_type = 2; //for individual file output
           alignment_type = 4;      // for single file output
           
           //some default values for the start and end indicies
           startIndex = 0;
           endIndex = 49;
           
           //the file containing the scoring matrix
           scoringMatrix = "NUC.4.4.txt";
           //minimum IGHJ length required
           minJLength = 30;
           
           //type of alignment - full, partial, ighv only
           alignmentModelType = 1;
           
           //number of vgenes to return alignment for
           numOfVgenes = 5;
           numSuboptimal = 5;
			
			//blast location
		   blastLocation = "/usr/local/ncbi/blast/bin/";*/

	}

        nonGUI myNonGUI = new nonGUI(startupSequencesFileName, IGHV_Filename, IGHD_Filename, IGHJ_Filename, IGKV_Filename, IGKJ_Filename, userSpecifiedChainType, probabilityFilename, htmlFilename, scoringSystem, htmlWriterType, dAcceptanceLevel, alignment_type, startIndex, endIndex, scoringMatrix, minJLength, alignmentModelType, numOfVgenes, numSuboptimal, blastLocation, IGLV_Filename, IGLJ_Filename, mutationSpecFileName, dProbsFilename);
        //GUI myGUI = new GUI(startupSequencesFileName);
        //myGUI.show();
    }


}
