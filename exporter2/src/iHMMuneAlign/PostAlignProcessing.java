package iHMMuneAlign;

import java.util.Vector;
import java.io.*;
import java.util.*;

// Referenced classes of package iHMMuneAlign:
//            GeneInfo, StopCodonSearch

public class PostAlignProcessing
{
    String dProbsFilename;

    public PostAlignProcessing(String dProbsFilename)
    {
        this.dProbsFilename = dProbsFilename;
    }

    public Vector noStopCodons(GeneInfo vGeneInfo, String sequence)
    {
        int umsVGeneStartPos = vGeneInfo.start_gene_pos;
        int nucleotideCodonStartPos = nucleotideCodonStartPos(umsVGeneStartPos);
        StopCodonSearch stopCodonSearch = new StopCodonSearch(sequence, nucleotideCodonStartPos);
        Vector stopCodons = stopCodonSearch.stopCodons;
        return stopCodons;
    }

    private int nucleotideCodonStartPos(int umsVGeneStartPos)
    {
        //int vGeneNucleotideOffset = umsVGeneStartPos - 1;
        //int startCodonNucleotideOffset = vGeneNucleotideOffset % 3;
        int startCodonNucleotideOffset = umsVGeneStartPos % 3;
        if(startCodonNucleotideOffset == 1)
        {
            return 0;
        }
        if(startCodonNucleotideOffset == 2)
        {
            return 2;
        }
        if(startCodonNucleotideOffset == 0)
        {
            return 1;
        } else
        {
            throw new Error("Invalid result");
        }
    }

    public boolean acceptDGene(GeneInfo dGene, int dGeneAlignmentAcceptanceType, RegionInfo n1Region, RegionInfo n2Region)
    {
        String geneString = dGene.aligned_gene_string;
        String alignmentString = dGene.aligned_ums_string;

        int alignmentLength = geneString.length();
        int mismatches = getMismatches(geneString, alignmentString);
        int consecutiveMatches = getConsecutiveMatches(geneString, alignmentString);
        
        //for minimum of 8 nts
        if(dGeneAlignmentAcceptanceType == 1)
        {
            if(alignmentLength > 12)
            {
                return analyzeTwelvePlusAlignment(geneString, alignmentString);
            } else
            {
                return hasEightMerDGeneAcceptance(alignmentLength, mismatches, consecutiveMatches);
            }
        }

        //for a minimum of 5
        if(dGeneAlignmentAcceptanceType == 2)
        {
            return hasFiveMerDGeneAcceptance(alignmentLength, mismatches, consecutiveMatches);
        }
        
        if (dGeneAlignmentAcceptanceType == 3) //dynamic probs for the minimum D based on the DLength/DMut/JunctionLength
        {
            //get the N1 and N2 region
            String N1 = "";
            String N2 = "";

            //get the junction length
            if (n1Region == null) {
               N1 = "";
            } else {
               N1 = n1Region.region_string;
            }
            
            if (n2Region == null) {
              N2 = "";
            } else {
              N2 = n2Region.region_string;
            }


            String junction = N1 + geneString + N2;
            int junctionLength = junction.length();

            double dProb = getAcceptableDynamic(alignmentLength, mismatches, junctionLength, dProbsFilename);

            //determine what to return based on the value of the d prob
            if (dProb < 0) { //represents an error as the d is too short
               //can't accept the d as too short
               return false;
            } else { //dealing with the any real prob returned
               //use 0.05 as the cut off
               if (dProb < 0.05) { //below the cut-off so accept
                  return true;
               } else { //above the cut-off so don't accept
                  return false;
               }
            }

        }

        //if no valid acceptance type found just return false
        return false;

    }

    public double getAcceptableDynamic(int dLength, int mismatches, int junctionLength, String dProbsFilename) {
      //have the dLength, the d mutations and the junction length
      //need to use these too calculate the prob of misidentification given the junction length

      //need to get the d probs from the file
      String dProbFilename = dProbsFilename;//"Dprobs.txt";
      File dProbFile = new File(dProbFilename);

      String[][] dProbsArray = getDProbs(dProbFile);


      //minimum of 4 nucleotides required for an acceptable D
      if (dLength < 4) {
         //returning the probs so they should be greater than one, use negative number to represent an error
         return -1D;
      } else {
          //need to adjust the index for the row length to correct starting at length 4
          int dRow = dLength - 3;

          //need to adjust the index for the cols as we have a title col at 0
          int dCol = dLength - mismatches + 1;

          //get the value and convert to integer
          String dProbStr = dProbsArray[dRow][dCol];
          double dProb = Double.parseDouble(dProbStr);
          
          //var to hold the final prob
          double finalProb = 0D;

         //need to account for the length of the junction
         //do this by multiplying the prob by (junction length - d length - 1)
          if (dLength == junctionLength) {
             finalProb = dProb;
          } else {
             double factor = junctionLength - dLength + 1;
             finalProb = factor * dProb;
          }
          
          return finalProb;
          //check if the prob is significant
          //if (finalProb < 0.05) {
          //  return true;
          //} else {
          //  return false;
          //}
      }

      //return false;
    } //- hasAcceptableDynamic

    private String[][] getDProbs(File dProbFile) {
       try {
       //create an array list to hold the lines
       ArrayList lines = new ArrayList();
       
       //open up the file
       BufferedReader br = new BufferedReader(new FileReader(dProbFile));

       //var to hol
       String str = null;

       //get all the lines from the file
       while ((str = br.readLine()) != null) {
         lines.add(str);
       }
       
       //2d array for the d-probs
       String[][] dProbs = new String[lines.size()][];

       for (int i = 0; i < dProbs.length; i++) {
          //get a line stored in the lines array
          str = (String) lines.get(i);
          
          //grab the probs from the line
          StringTokenizer st = new StringTokenizer(str, "\t");
          //create a temp array to hold the probs from this line now that they are split
          String[] tempArray = new String[st.countTokens()];

          //for each piece of data on the current line
          for (int j=0; j < tempArray.length; j++) {
             tempArray[j] = st.nextToken();
          }

          dProbs[i] = tempArray;
       }

       return dProbs;
       }  catch(IOException ioe)
       {
            throw new Error("load misidentifcation D probs file failed:  " + ioe.getMessage());
       }
    }

    private boolean hasEightMerDGeneAcceptance(int alignmentLength, int mismatches, int consecutiveMathces)
    {
        if(alignmentLength < 8)
        {
            return false;
        }
        if(alignmentLength < 10)
        {
            return mismatches <= 0;
        }
        if(alignmentLength < 12)
        {
            return mismatches <= 1;
        }
        if(alignmentLength == 12)
        {
            return mismatches <= 2;
        } else
        {
            return false;
        }
    }

    private boolean hasFiveMerDGeneAcceptance(int alignmentLength, int mismatches, int consecutiveMathces)
    {
        if(alignmentLength < 5)
        {
            return false;
        }
        return consecutiveMathces >= 5;
    }

    private int getConsecutiveMatches(String geneString, String alignmentString)
    {
        int consecutiveMatches = 0;
        int longestConsecutiveMatches = 0;
        for(int i = 0; i < geneString.length(); i++)
        {
            char currGeneNucl = geneString.charAt(i);
            char currAlignedNucl = alignmentString.charAt(i);
            if(currGeneNucl == currAlignedNucl)
            {
                consecutiveMatches++;
            } else
            {
                if(consecutiveMatches > longestConsecutiveMatches)
                {
                    longestConsecutiveMatches = consecutiveMatches;
                }
                consecutiveMatches = 0;
            }
        }

        if(consecutiveMatches > longestConsecutiveMatches)
        {
            longestConsecutiveMatches = consecutiveMatches;
        }
        return longestConsecutiveMatches;
    }

    private boolean analyzeTwelvePlusAlignment(String geneString, String alignmentString)
    {
        for(int i = 0; i < geneString.length(); i++)
        {
            int matches = 0;
            int mismatches = 0;
            for(int j = i; j < geneString.length(); j++)
            {
                char geneNucl = geneString.charAt(i);
                char alignmentNucl = alignmentString.charAt(i);
                if(geneNucl != alignmentNucl)
                {
                    if(++mismatches > 2)
                    {
                        break;
                    }
                } else
                if(++matches >= 8)
                {
                    if(mismatches == 0)
                    {
                        return true;
                    }
                    if(matches >= 10)
                    {
                        if(mismatches <= 1)
                        {
                            return true;
                        }
                        if(matches >= 12 && mismatches <= 2)
                        {
                            return true;
                        }
                    }
                }
            }

        }

        return false;
    }

    private int getMismatches(String geneString, String alignmentString)
    {
        int mismatches = 0;
        for(int i = 0; i < geneString.length(); i++)
        {
            char geneNucl = geneString.charAt(i);
            char alignmentNucl = alignmentString.charAt(i);
            if(geneNucl != alignmentNucl)
            {
                mismatches++;
            }
        }

        return mismatches;
    }
}
