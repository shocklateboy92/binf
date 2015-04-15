package iHMMuneAlign;

import java.io.PrintStream;

// Referenced classes of package iHMMuneAlign:
//            GeneInfo, CodonHelperFunctions

public class JGeneReadingFrame
{

    private boolean jGeneInFrame;
    private int WGXGmotifInAlignedJGenePos;
    boolean DEBUGGING;

    public JGeneReadingFrame(String sequence, GeneInfo jGeneInfo, GeneInfo vGeneInfo)
    {
        jGeneInFrame = isJGeneInFrame(sequence, jGeneInfo, vGeneInfo);
        DEBUGGING = false;
    }

    private boolean isJGeneInFrame(String sequence, GeneInfo jGeneInfo, GeneInfo vGeneInfo)
    {
        String alignedJGeneString = jGeneInfo.aligned_gene_string;
        int alignedJGeneWGXGmotifPos = findAlignedJGeneWGXGpos(alignedJGeneString);
        if (DEBUGGING) {
        	System.out.println("Aligned JGene motif position: " + alignedJGeneWGXGmotifPos);
        }
        int jGeneStartPos = jGeneInfo.start_gene_pos;
        int alignedJGeneSequencePos = sequence.indexOf(jGeneInfo.aligned_ums_string);
        if (DEBUGGING) {
        	System.out.println("Aligned alignedJGeneSequencePos: " + alignedJGeneSequencePos);
        }
        int WGXGmotifIndexInSequence = alignedJGeneSequencePos + alignedJGeneWGXGmotifPos;
        int vGeneOffset = vGeneInfo.start_gene_pos - 1;
        int nuclPrecedingWGXGmotif = vGeneOffset + WGXGmotifIndexInSequence;
        if (DEBUGGING) {
        	System.out.println("Found the WGXG-motif at position: " + alignedJGeneWGXGmotifPos);
        }
        WGXGmotifInAlignedJGenePos = alignedJGeneWGXGmotifPos;
        if(nuclPrecedingWGXGmotif % 3 == 0)
        {
            if (DEBUGGING) {
            	System.out.println("JGene in Frame");
            }
            return true;
        } else
        {
            if (DEBUGGING) {
            	System.out.println("JGene NOT in Frame");
            }
            return false;
        }
    }

    private int findAlignedJGeneWGXGpos(String alignedJGeneString)
    {
        char alignedJGeneStringC[] = alignedJGeneString.toCharArray();
        for(int i = 0; i < alignedJGeneStringC.length - 2; i++)
        {
            char currCodon[] = {
                alignedJGeneStringC[i], alignedJGeneStringC[i + 1], alignedJGeneStringC[i + 2]
            };
            if(CodonHelperFunctions.isAminoAcidW(currCodon))
            {
                if (DEBUGGING) {
                	System.out.println("found the W-codon in JGene at pos: " + i);
                }
                if(hasGXGaminoAcids(i + 3, alignedJGeneStringC))
                {
                    return i;
                }
            }
        }

        return -1;
    }

    private boolean hasGXGaminoAcids(int startPos, char alignedJGeneStringC[])
    {
        int alignedJGeneLength = alignedJGeneStringC.length;
        int nucleotidesRequired = 9;
        int remainingNucleotidesInJGene = alignedJGeneLength - startPos;
        if(remainingNucleotidesInJGene < nucleotidesRequired)
        {
            return false;
        }
        char codonOneG[] = {
            alignedJGeneStringC[startPos], alignedJGeneStringC[startPos + 1], alignedJGeneStringC[startPos + 2]
        };
        char codonTwoX[] = {
            alignedJGeneStringC[startPos + 3], alignedJGeneStringC[startPos + 4], alignedJGeneStringC[startPos + 5]
        };
        char codonThreeG[] = {
            alignedJGeneStringC[startPos + 6], alignedJGeneStringC[startPos + 7], alignedJGeneStringC[startPos + 8]
        };
        return CodonHelperFunctions.isAminoAcidG(codonOneG) && CodonHelperFunctions.isAminoAcidX(codonTwoX) && CodonHelperFunctions.isAminoAcidG(codonThreeG);
    }

    public boolean isJGeneInFrame()
    {
        return jGeneInFrame;
    }

    public int getWGXGmotifInAlignedJGenePos()
    {
        return WGXGmotifInAlignedJGenePos;
    }
}
