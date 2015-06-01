package iHMMuneAlign;

import java.io.PrintStream;
import java.util.Vector;

// Referenced classes of package iHMMuneAlign:
//            Codon, CodonHelperFunctions

public class StopCodonSearch
{

    final int DEFAULT_FRAME_READ_POS = 0;
    String sequenceString;
    Vector stopCodons;
    boolean DEBUGGING = false;

    public StopCodonSearch(String sequenceString, int readingFramePos)
    {
        this.sequenceString = null;
        stopCodons = null;
        this.sequenceString = sequenceString;
        if(readingFramePos == -1)
        {
            stopCodons = findStopCodons(sequenceString, 0);
        } else
        {
            stopCodons = findStopCodons(sequenceString, readingFramePos);
        }
    }

    private Vector findStopCodons(String sequenceString, int codonStartPos)
    {
        if (DEBUGGING) {
        	System.out.println("sequence string: \n" + sequenceString);
        	System.out.println("codon start nucl pos: " + codonStartPos);
        }
        Vector stopCodonPositions = new Vector();
        
        for(int i = codonStartPos; i <= sequenceString.length() - 3; i += 3)
        {
            char currCodonC[] = getCodon(sequenceString, i);
            if (DEBUGGING) {
	            System.out.println("Current Codon: " + new String(currCodonC));
	        }
            if(isStopCodon(currCodonC))
            {
                Codon currCodon = new Codon(currCodonC, sequenceString, i);
                stopCodonPositions.add(currCodon);
            }
        }

        return stopCodonPositions;
    }

    private boolean isStopCodon(char codon[])
    {
        if(CodonHelperFunctions.isTAA(codon))
        {
            return true;
        }
        if(CodonHelperFunctions.isTAG(codon))
        {
            return true;
        }
        return CodonHelperFunctions.isTGA(codon);
    }

    private char[] getCodon(String sequenceString, int codonNuclStartPos)
    {
        int nucleotideOnePos = codonNuclStartPos;
        int nucleotideTwoPos = codonNuclStartPos + 1;
        int nucleotideThreePos = codonNuclStartPos + 2;
        char nucleotideOne = sequenceString.charAt(nucleotideOnePos);
        char nucleotideTwo = sequenceString.charAt(nucleotideTwoPos);
        char nucleotideThree = sequenceString.charAt(nucleotideThreePos);
        char codon[] = {
            nucleotideOne, nucleotideTwo, nucleotideThree
        };
        return codon;
    }

    public boolean hasStopCodons()
    {
        return stopCodons.size() > 0;
    }
}
