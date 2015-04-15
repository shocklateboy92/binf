package iHMMuneAlign;

import java.io.*;
import java.util.ArrayList;
import org.biojava.bio.BioException;
//import org.biojava.bio.seq.Sequence;
//import org.biojava.bio.seq.SequenceIterator;
//import org.biojava.bio.seq.io.SeqIOTools;
import org.biojavax.Namespace;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

public class FastaReader
{

    public FastaReader()
    {
    }

    public static void main(String args[])
    {
        FastaReader fastaReader = new FastaReader();
        try
        {
            ArrayList dGenes = fastaReader.readFile(new File("d-regions.fa"));
            int totalLength = 0;
            int numberOfDGenes = dGenes.size();
            for(int i = 0; i < numberOfDGenes; i++)
            {
                RichSequence currGene = (RichSequence)dGenes.get(i);
                totalLength += currGene.length();
            }

            double averageLength = (double)totalLength / (double)numberOfDGenes;
            System.out.println("average D-Gene length = " + averageLength);
        }
        catch(Exception ex)
        {
            throw new Error(ex.getMessage());
        }
    }

    public ArrayList readFile(File file)
    {
        ArrayList sequences = new ArrayList();
        BufferedReader br = null;
        try
        {
            br = new BufferedReader(new FileReader(file));
        }
        catch(FileNotFoundException fnfe)
        {
            System.out.println("The specified file was not found");
            return null;
        }
		//namespace
		Namespace ns = RichObjectFactory.getDefaultNamespace();
		RichSequenceIterator stream = RichSequence.IOTools.readFastaDNA(br,ns);
        //int counter = 0;
        RichSequence seq = null;
        while(stream.hasNext()) 
        {
            try
            {
                seq = stream.nextRichSequence();
            }
            catch(BioException bioe)
            {
                System.out.println("Bio Exception when attempting to read from fasta file");
                return null;
            }
            sequences.add(seq);
            //counter++;
        }
        return sequences;
    }
}
