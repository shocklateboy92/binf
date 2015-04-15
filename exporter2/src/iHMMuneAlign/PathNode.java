package iHMMuneAlign;

import java.io.*;
import java.util.ArrayList;
import java.util.Random;

import org.biojavax.bio.seq.RichSequence;
import org.biojavax.RichObjectFactory;
import org.biojavax.Namespace;

//import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.SymbolList;

public class PathNode {
   private String DGeneName; // name of chosen D gene in path
   private String JGeneName; // name of chosen J gene in path
   private File DFile; // D file used to build HMM
   private File JFile; // J file used to build HMM
   private double score;
   private SymbolList viterbiStateSeq;
   private String UMS_No_C;
   
   public PathNode(String DNameToRemove, String JNameToRemove, File baseDFile, File baseJFile) { 
	   File newDFile = baseDFile;
	   //newDFile.deleteOnExit(); 
	   File newJFile = baseJFile;
	   //newJFile.deleteOnExit(); 
	   
	   if (!DNameToRemove.equals("")) {
		   try {
			   newDFile = removeSeqFromFile(DNameToRemove, baseDFile);
		   }
		   catch (IOException ioe) {
			   System.out.println("Error creating modified repertoire files");
		   }
	   } 
	   if (!JNameToRemove.equals("")) {
		   try {
			   newJFile = removeSeqFromFile(JNameToRemove, baseJFile);
		   }
		   catch (IOException ioe) {
			   System.out.println("Error creating modified repertoire files");
		   }
	   } 
	   this.DFile = newDFile;
	   this.JFile = newJFile; 
	   //this.isHeavy = isHeavy;
	   //this.score = score;
       //this.viterbiStateSeq = stateSeq;
   }
   
   public PathNode () {
	   this.DGeneName = "";
	   this.JGeneName = "";
	   // initialise to lowest score
	   this.score = 0-java.lang.Float.MAX_VALUE;
	   this.viterbiStateSeq = null;
   }
   
   public String getDName () {
	   return this.DGeneName;
   }
   
   public String getJName () {
	   return this.JGeneName;
   }
   
   public void setDName (String DName) {
	   this.DGeneName = DName;
   }
   
   public void setJName (String JName) {
	   this.JGeneName = JName;
   }
   
   public File getDFile () {
	   return this.DFile;
   }
   
   public File getJFile () {
	   return this.JFile;
   }
   
   public double getScore () {
	   return this.score;
   }
   
   public void setScore (double score) {
	   this.score = score;
   }
   
   public SymbolList getStateSeq() {
	   return this.viterbiStateSeq;
   }
   
   public void setStateSeq(SymbolList stateSeq) {
	   this.viterbiStateSeq = stateSeq; 
   }
   
   public String getUMSNoC() {
	   return this.UMS_No_C;
   }
   
   public void setUMSNoC(String UMSNoC) {
	   this.UMS_No_C = UMSNoC;
   }
   
   private File removeSeqFromFile(String removeName, File baseFile) throws IOException {
	   // generate a random file name
	   Random rnd = new Random();
	   String newFileName = "modifiedRep"+rnd.nextInt()+".fa";
	   FileOutputStream fostream = new FileOutputStream(newFileName);
	   
	   //get default namespace to use for FASTA writing
	   Namespace ns = RichObjectFactory.getDefaultNamespace();
       
       // read in seqs from the baseFile
	   ArrayList<RichSequence> genes = new ArrayList<RichSequence>();  
       try
       {
           FastaReader fr = new FastaReader();
           genes = fr.readFile(baseFile);
       }
       catch(Exception e)
       {
           e.printStackTrace();
           throw new Error(e.getMessage());
       }
       
       // write all seqs to newFile except the one to be removed
       for(int i = 0; i < genes.size(); i++)
       {
           RichSequence seq = (RichSequence)genes.get(i);
           if (!seq.getName().equals(removeName)) {
        	   RichSequence.IOTools.writeFasta(fostream, seq, ns);
           }
       }
       
       File newFile = new File(newFileName);
       // Delete temp file when program exits
       newFile.deleteOnExit();
       return newFile;
   }
}
