package iHMMuneAlign;

import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Vector;

// Referenced classes of package iHMMuneAlign:
//            HTML_Writer, GeneInfo, RegionInfo

public class HTMLforWWW_Writer extends HTML_Writer
{

    public HTMLforWWW_Writer()
    {
    }

    public void createDocument(String fileName, String sequenceName, int chainType, GeneInfo VGene, RegionInfo P1, RegionInfo N1, RegionInfo P2, GeneInfo DGene, 
            RegionInfo P3, RegionInfo N2, RegionInfo P4, GeneInfo JGene, String selectedSequenceName, String selectedSequenceString, double score, byte alignmentType, Vector stopCodons, int relativeMotifPosition, boolean isJGeneInFrame, int reverseComplement, String dProbsFilename)
    {
        if(fostream == null)
        {
            createFile(fileName, alignmentType);
        }
        writeTitle("title");
        if(fostream == null)
        {
            throw new Error("Fo stream is null ?????????");
        }
        boolean acceptedDGene = DGene.acceptedAlignment;
        int VGENE_END_SIZE = 20;
        int completeVGeneLength = VGene.gene_length;
        String VGeneInfoString = getVGeneString(VGene, 20);
        String DGeneInfoString = getDorJGeneString(DGene, false, -1);
        String JGeneInfoString = getDorJGeneString(JGene, true, relativeMotifPosition);
        String unidentifiedRegion = "";
        String p1String = null;
        if(acceptedDGene)
        {
            if(P1 == null)
            {
                p1String = null;
            } else
            {
                p1String = markupNandXNucleotides(minLengthTwo(P1.region_string));
            }
        } else
        if(P1 == null)
        {
            p1String = null;
        } else
        {
            unidentifiedRegion = unidentifiedRegion + P1.region_string;
        }
        String n1String = null;
        if(acceptedDGene)
        {
            if(N1 == null)
            {
                n1String = null;
            } else
            {
                n1String = markupNandXNucleotides(minLengthTwo(N1.region_string));
            }
        } else
        if(N1 == null)
        {
            n1String = null;
        } else
        {
            unidentifiedRegion = unidentifiedRegion + N1.region_string;
        }
        String p2String = null;
        if(acceptedDGene)
        {
            if(P2 == null)
            {
                p2String = null;
            } else
            {
                p2String = markupNandXNucleotides(minLengthTwo(P2.region_string));
            }
        } else
        if(P2 == null)
        {
            p2String = null;
        } else
        {
            unidentifiedRegion = unidentifiedRegion + P2.region_string;
        }
        unidentifiedRegion = unidentifiedRegion + DGene.aligned_ums_string;
        String p3String = null;
        if(acceptedDGene)
        {
            if(P3 == null)
            {
                p3String = null;
            } else
            {
                p3String = markupNandXNucleotides(minLengthTwo(P3.region_string));
            }
        } else
        if(P3 == null)
        {
            p3String = null;
        } else
        {
            unidentifiedRegion = unidentifiedRegion + P3.region_string;
        }
        String n2String = null;
        if(acceptedDGene)
        {
            if(N2 == null)
            {
                n2String = null;
            } else
            {
                n2String = markupNandXNucleotides(minLengthTwo(N2.region_string));
            }
        } else
        if(N2 == null)
        {
            n2String = null;
        } else
        {
            unidentifiedRegion = unidentifiedRegion + N2.region_string;
        }
        String p4String = null;
        if(acceptedDGene)
        {
            if(P4 == null)
            {
                p4String = null;
            } else
            {
                p4String = markupNandXNucleotides(minLengthTwo(P4.region_string));
            }
        } else
        if(P4 == null)
        {
            p4String = null;
        } else
        {
            unidentifiedRegion = unidentifiedRegion + P4.region_string;
        }
        String markedUpUnidentifiedRegion = markupRegion(unidentifiedRegion);
        int VGeneMutations = mutationsInGeneAlignment(VGene.aligned_gene_string, VGene.aligned_ums_string);
        int DGeneMutations = 0;
        if(acceptedDGene)
        {
            DGeneMutations = mutationsInGeneAlignment(DGene.aligned_gene_string, DGene.aligned_ums_string);
        }
        int JGeneMutations = mutationsInGeneAlignment(JGene.aligned_gene_string, JGene.aligned_ums_string);
        String mutationsInfoString = createMutationsInfoString(VGeneMutations, DGeneMutations, JGeneMutations);
        int NandXnucleotideCount = getNandXnucleotideCount(selectedSequenceString);
        String NandXnucleotidesString = createNandXNuclString(NandXnucleotideCount);
        String isJGeneInFrameString = createIsJGeneInFrameString(isJGeneInFrame);
        String stopCodonsInSequenceString = createHasSequenceStopCodonsString(stopCodons);
        String geneNamesString;
        if(acceptedDGene)
        {
            geneNamesString = createGeneNamesString(VGene.gene_name, DGene.gene_name, JGene.gene_name);
        } else
        {
            geneNamesString = createGeneNamesString(VGene.gene_name, "NO DGENE ALIGNMENT", JGene.gene_name);
        }
        String sequenceNameString = createSequenceNameString(selectedSequenceName);
        String regionInfoString;
        if(acceptedDGene)
        {
            regionInfoString = createRegionInfoString(VGeneInfoString, DGeneInfoString, JGeneInfoString, n1String, n2String, p1String, p2String, p3String, p4String);
        } else
        {
            regionInfoString = createRegionInfoString(VGeneInfoString, markedUpUnidentifiedRegion, JGeneInfoString, n1String, n2String, p1String, p2String, p3String, p4String);
        }
        int VRegionLength = 20;
        int DRegionLength;
        if(acceptedDGene)
        {
            DRegionLength = DGene.gene_length;
        } else
        {
            DRegionLength = unidentifiedRegion.length();
        }
        int JRegionLength = JGene.gene_length;
        String regionHeaderString = createRegionHeaderString(VRegionLength, DRegionLength, JRegionLength, n1String, n2String, p1String, p2String, p3String, p4String, acceptedDGene);
        System.out.println("*****Region Header String: " + regionHeaderString);
        writeToFile(sequenceNameString, regionHeaderString, regionInfoString, geneNamesString, mutationsInfoString, NandXnucleotidesString, isJGeneInFrameString, stopCodonsInSequenceString);
    }
    
	public void createDocument(String fileName, String sequenceName, int chainType, GeneInfo VGene, RegionInfo N1, GeneInfo DGene, RegionInfo N2, GeneInfo JGene, String selectedSequenceName, String selectedSequenceString, double score, byte alignmentType, Vector stopCodons, int relativeMotifPosition, boolean isJGeneInFrame, int reverseComplement, String dProbsFilename)
    {
        if(fostream == null)
        {
            createFile(fileName, alignmentType);
        }
        writeTitle("title");
        if(fostream == null)
        {
            throw new Error("Fo stream is null ?????????");
        }
        boolean acceptedDGene = DGene.acceptedAlignment;
        int VGENE_END_SIZE = 20;
        int completeVGeneLength = VGene.gene_length;
        String VGeneInfoString = getVGeneString(VGene, 20);
        String DGeneInfoString = getDorJGeneString(DGene, false, -1);
        String JGeneInfoString = getDorJGeneString(JGene, true, relativeMotifPosition);
        String unidentifiedRegion = "";
		
		String n1String = null;
		if(acceptedDGene) {
            if(N1 != null) {
                n1String = markupNandXNucleotides(minLengthTwo(N1.region_string));
            }
        } else {
			if(N1 != null) {
				unidentifiedRegion = unidentifiedRegion + N1.region_string;
			}
		}
       
        String n2String = null;
        if(acceptedDGene) {
            if(N2 != null) {
                n2String = markupNandXNucleotides(minLengthTwo(N2.region_string));
            }
        } else {
			if(N2 != null) {
				unidentifiedRegion = unidentifiedRegion + DGene.aligned_ums_string;
				unidentifiedRegion = unidentifiedRegion + N2.region_string;
			}
		}
		
		String markedUpUnidentifiedRegion = markupRegion(unidentifiedRegion);
        int VGeneMutations = mutationsInGeneAlignment(VGene.aligned_gene_string, VGene.aligned_ums_string);
        int DGeneMutations = 0;
        if(acceptedDGene)
        {
            DGeneMutations = mutationsInGeneAlignment(DGene.aligned_gene_string, DGene.aligned_ums_string);
        }
        int JGeneMutations = mutationsInGeneAlignment(JGene.aligned_gene_string, JGene.aligned_ums_string);
        String mutationsInfoString = createMutationsInfoString(VGeneMutations, DGeneMutations, JGeneMutations);
        int NandXnucleotideCount = getNandXnucleotideCount(selectedSequenceString);
        String NandXnucleotidesString = createNandXNuclString(NandXnucleotideCount);
        String isJGeneInFrameString = createIsJGeneInFrameString(isJGeneInFrame);
        String stopCodonsInSequenceString = createHasSequenceStopCodonsString(stopCodons);
        String geneNamesString;
        if(acceptedDGene)
        {
            geneNamesString = createGeneNamesString(VGene.gene_name, DGene.gene_name, JGene.gene_name);
        } else
        {
            geneNamesString = createGeneNamesString(VGene.gene_name, "NO DGENE ALIGNMENT", JGene.gene_name);
        }
        String sequenceNameString = createSequenceNameString(selectedSequenceName);
        String regionInfoString;
        if(acceptedDGene)
        {
            regionInfoString = createRegionInfoString(VGeneInfoString, DGeneInfoString, JGeneInfoString, n1String, n2String);
        } else
        {
            regionInfoString = createRegionInfoString(VGeneInfoString, markedUpUnidentifiedRegion, JGeneInfoString, n1String, n2String);
        }
        int VRegionLength = 20;
        int DRegionLength;
        if(acceptedDGene)
        {
            DRegionLength = DGene.gene_length;
        } else
        {
            DRegionLength = unidentifiedRegion.length();
        }
        int JRegionLength = JGene.gene_length;
        String regionHeaderString = createRegionHeaderString(VRegionLength, DRegionLength, JRegionLength, n1String, n2String, acceptedDGene);
        System.out.println("*****Region Header String: " + regionHeaderString);
        writeToFile(sequenceNameString, regionHeaderString, regionInfoString, geneNamesString, mutationsInfoString, NandXnucleotidesString, isJGeneInFrameString, stopCodonsInSequenceString);
    }
    
	
    public void createErrDocument(String fileName, String sequenceName, byte alignmentType, String errorType, String vOutput)
    {
      // creates the output when the model creation hasn't been sucessful
      //for example an A_probability hasn't been able to be calculated due to a short IGHV gene
      //open up the output file
      if(fostream == null)
        {
            createFile(fileName, alignmentType);
        }

        if(fostream == null)
        {
            throw new Error("Fo stream is null ?????????");
        }
      //create place holders for the different outputs
      String regionHeaderString = "" + errorType;
      String regionInfoString = "NA";
      String geneNamesString = "NA";
      String mutationsInfoString = "NA";
      String NandXnucleotidesString = "NA";
      String stopCodonsInSequenceString = "NA";
      String isJGeneInFrameString = "NA";
      String sequenceNameString = sequenceName;

      // write the output for the sequence
      writeToFile(sequenceNameString, regionHeaderString, regionInfoString, geneNamesString, mutationsInfoString, NandXnucleotidesString, isJGeneInFrameString, stopCodonsInSequenceString);

    }

    private void writeToFile(String sequenceNameString, String regionHeaderString, String regionInfoString, String geneNameInfoString, String mutationsInfoString, String NandXnucleotidesString, String isJGeneInFrameString, 
            String stopCodonsInSequenceString)
    {
        fostream.println(sequenceNameString);
        fostream.println(regionHeaderString);
        fostream.println(regionInfoString);
        fostream.println(geneNameInfoString);
        fostream.println(mutationsInfoString);
        fostream.println(NandXnucleotidesString);
        fostream.println(isJGeneInFrameString);
        fostream.println(stopCodonsInSequenceString);
    }

    private String minLengthTwo(String region)
    {
        if(region.length() < 1)
        {
            throw new Error("minLengthTwo(): length found = zero");
        }
        if(region.length() == 1)
        {
            return region + " ";
        } else
        {
            return region;
        }
    }

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

    private String createHasSequenceStopCodonsString(Vector stopCodons)
    {
        String result = "<pre>";
        if(stopCodons.size() == 0)
        {
            result = result + "NO Stop-Codons in Sequence";
        } else
        {
            result = result + "NB! " + stopCodons.size() + " Stop Codons in Sequence";
        }
        result = result + "</pre>";
        return result;
    }

    private String createIsJGeneInFrameString(boolean isJGeneInFrame)
    {
        String result = "<pre>";
        if(isJGeneInFrame)
        {
            result = result + "J-Gene is read IN-FRAME";
        } else
        {
            result = result + "NB! J-Gene is NOT read IN-FRAME";
        }
        result = result + "</pre>";
        return result;
    }

    private String createNandXNuclString(int NandXnucleotideCount)
    {
        String NandXnucleotidesString = "<pre>";
        NandXnucleotidesString = NandXnucleotidesString + "N / X nucleotides in sequence:  " + NandXnucleotideCount;
        NandXnucleotidesString = NandXnucleotidesString + "</pre>";
        return NandXnucleotidesString;
    }

    private String createMutationsInfoString(int VGeneMutations, int DGeneMutations, int JGeneMutations)
    {
        String mutationsInfoString = "<pre>";
        mutationsInfoString = mutationsInfoString + "Vmut\tDmut\tJmut\n";
        mutationsInfoString = mutationsInfoString + " " + VGeneMutations + "\t" + " " + DGeneMutations + "\t" + " " + JGeneMutations;
        mutationsInfoString = mutationsInfoString + "</pre>";
        return mutationsInfoString;
    }

    private String createGeneNamesString(String VGeneName, String DGeneName, String JGeneName)
    {
        String GeneNameString = "<pre>";
        String vNameHeader = "V name";
        String dNameHeader = "D name";
        String jNameHeader = "J name";
        String centeredVNameHeader = getCenteredHeaderString(VGeneName.length(), vNameHeader);
        String centeredDNameHeader = getCenteredHeaderString(DGeneName.length(), dNameHeader);
        String centeredJNameHeader = getCenteredHeaderString(JGeneName.length(), jNameHeader);
        GeneNameString = GeneNameString + centeredVNameHeader + "\t" + centeredDNameHeader + "\t" + centeredJNameHeader + "\n";
        GeneNameString = GeneNameString + VGeneName + "\t" + DGeneName + "\t" + JGeneName;
        GeneNameString = GeneNameString + "</pre>";
        return GeneNameString;
    }

    private String createSequenceNameString(String sequenceName)
    {
        String sequenceNameString = "<pre><h3>";
        sequenceNameString = sequenceNameString + "Input Sequence: " + sequenceName + "\t";
        sequenceNameString = sequenceNameString + "</pre></h3>";
        return sequenceNameString;
    }

    private String createRegionHeaderString(int VRegionLength, int DRegionLength, int JRegionLength, String N1, String N2, String P1, String P2, 
            String P3, String P4, boolean acceptedDGene)
    {
        String headerString = "<pre>";
        headerString = headerString + getCenteredHeaderString(VRegionLength, "VGene (end)") + " ";
        if(P1 != null)
        {
            headerString = headerString + getCenteredHeaderString(P1.length(), "P1") + " ";
        }
        if(N1 != null)
        {
            headerString = headerString + getCenteredHeaderString(N1.length(), "N1") + " ";
        }
        if(P2 != null)
        {
            headerString = headerString + getCenteredHeaderString(P2.length(), "P2") + " ";
        }
        if(acceptedDGene)
        {
            headerString = headerString + getCenteredHeaderString(DRegionLength, "DGene") + " ";
        } else
        {
            headerString = headerString + getCenteredHeaderString(DRegionLength, "UNIDENTIFIED REGION") + " ";
        }
        if(P3 != null)
        {
            headerString = headerString + getCenteredHeaderString(P3.length(), "P3") + " ";
        }
        if(N2 != null)
        {
            headerString = headerString + getCenteredHeaderString(N2.length(), "N2") + " ";
        }
        if(P4 != null)
        {
            headerString = headerString + getCenteredHeaderString(P4.length(), "P4") + " ";
        }
        headerString = headerString + getCenteredHeaderString(JRegionLength, "JGene") + " ";
        headerString = headerString + "</pre>";
        return headerString;
    }
	
	private String createRegionHeaderString(int VRegionLength, int DRegionLength, int JRegionLength, String N1, String N2, boolean acceptedDGene)
    {
        String headerString = "<pre>";
        headerString = headerString + getCenteredHeaderString(VRegionLength, "VGene (end)") + " ";
        if(N1 != null)
        {
            headerString = headerString + getCenteredHeaderString(N1.length(), "N1") + " ";
        }
        if(acceptedDGene)
        {
            headerString = headerString + getCenteredHeaderString(DRegionLength, "DGene") + " ";
        } else
        {
            headerString = headerString + getCenteredHeaderString(DRegionLength, "UNIDENTIFIED REGION") + " ";
        }
        if(N2 != null)
        {
            headerString = headerString + getCenteredHeaderString(N2.length(), "N2") + " ";
        }
        headerString = headerString + getCenteredHeaderString(JRegionLength, "JGene") + " ";
        headerString = headerString + "</pre>";
        return headerString;
    }

    private String createRegionInfoString(String VRegionString, String DRegionString, String JRegionString, String N1, String N2, String P1, String P2, 
            String P3, String P4)
    {
        String regionInfoString = "<pre>";
        regionInfoString = regionInfoString + VRegionString + " ";
        if(P1 != null)
        {
            regionInfoString = regionInfoString + P1 + " ";
        }
        if(N1 != null)
        {
            regionInfoString = regionInfoString + N1 + " ";
        }
        if(P2 != null)
        {
            regionInfoString = regionInfoString + P2 + " ";
        }
        regionInfoString = regionInfoString + DRegionString + " ";
        if(P3 != null)
        {
            regionInfoString = regionInfoString + P3 + " ";
        }
        if(N2 != null)
        {
            regionInfoString = regionInfoString + N2 + " ";
        }
        if(P4 != null)
        {
            regionInfoString = regionInfoString + P4 + " ";
        }
        regionInfoString = regionInfoString + JRegionString + " ";
        regionInfoString = regionInfoString + "</pre>";
        return regionInfoString;
    }

	private String createRegionInfoString(String VRegionString, String DRegionString, String JRegionString, String N1, String N2)
    {
        String regionInfoString = "<pre>";
        regionInfoString = regionInfoString + VRegionString + " ";
        if(N1 != null)
        {
            regionInfoString = regionInfoString + N1 + " ";
        }
        regionInfoString = regionInfoString + DRegionString + " ";
        if(N2 != null)
        {
            regionInfoString = regionInfoString + N2 + " ";
        }
        regionInfoString = regionInfoString + JRegionString + " ";
        regionInfoString = regionInfoString + "</pre>";
        return regionInfoString;
    }
	
    private String getCenteredHeaderString(int regionLength, String header)
    {
        int headerLength = header.length();
        if(headerLength >= regionLength)
        {
            return header;
        }
        int difference = regionLength - headerLength;
        int remainder = difference % 2;
        int startOffset;
        if(remainder != 0)
        {
            startOffset = difference / 2 + 1;
        } else
        {
            startOffset = difference / 2;
        }
        int endOffset = difference / 2;
        String startPadding = spaceString(startOffset);
        String endPadding = spaceString(endOffset);
        String headerString = startPadding + header + endPadding;
        return headerString;
    }

    private String spaceString(int emptyStringLength)
    {
        String result = "";
        for(int i = 0; i < emptyStringLength; i++)
        {
            result = result + " ";
        }

        return result;
    }

    private String getVGeneString(GeneInfo VGene, int VGENE_END_SIZE)
    {
        String name = VGene.gene_name;
        String alignedGeneString = VGene.aligned_gene_string;
        String alignedUMSString = VGene.aligned_ums_string;
        int completeGeneLength = VGene.gene_length;
        int startGenePos = VGene.start_gene_pos;
        int endGenePos = VGene.end_gene_pos;
        if(alignedGeneString.length() != alignedUMSString.length())
        {
            throw new Error("aligned VGene String and aligned UMS String not of same length");
        }
        String resultGeneInfo = "";
        int alignedGeneLength = alignedGeneString.length();
        int endGeneOffset = completeGeneLength - endGenePos;
        int alignedEndSize = VGENE_END_SIZE - endGeneOffset;
        String alignedGeneEndString = alignedGeneString.substring(alignedGeneLength - alignedEndSize, alignedGeneLength);
        String alignedUMSEndString = alignedUMSString.substring(alignedGeneLength - alignedEndSize, alignedGeneLength);
        int mutationsInAlignmentEnd = 0;
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
                        resultGeneInfo = resultGeneInfo + getUnspecifiedNucl(currUMSNucl);
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
            if(i >= startGenePos && i <= endDorJGenePos)
            {
                int alignedJGeneIndex = i - startGenePos;
                boolean isCurrNuclInMotif = false;
                if(isJGene && alignedJGeneIndex >= jGeneWGXGmotifPos && alignedJGeneIndex <= jGeneWGXGmotifPos + 11)
                {
                    isCurrNuclInMotif = true;
                }
                char currGeneNucl = alignedDorJGeneString.charAt(i - startGenePos);
                char currUMSNucl = alignedUMSString.charAt(i - startGenePos);
                if(currGeneNucl == currUMSNucl)
                {
                    if(isCurrNuclInMotif)
                    {
                        resultGeneInfo = resultGeneInfo + getWGXGmotifNucl(currUMSNucl);
                    } else
                    {
                        resultGeneInfo = resultGeneInfo + currGeneNucl;
                    }
                } else
                {
                    if(currUMSNucl == 'n' || currUMSNucl == 'x')
                    {
                        resultGeneInfo = resultGeneInfo + getUnspecifiedNucl(currUMSNucl);
                    } else
                    if(isCurrNuclInMotif)
                    {
                        resultGeneInfo = resultGeneInfo + getMutatedWGXGmotifNucl(currUMSNucl);
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

    private void writeTitle(String title)
    {
        fostream.print("<title>" + title + "</title>");
    }

    private void writeMutatedNucl(char mutNucl)
    {
        fostream.print("<font color=\"FF0000\">" + mutNucl + "</font>");
    }

    private String getMutatedNucl(char mutNucl)
    {
        String result = "<font color=\"FF0000\">" + mutNucl + "</font>";
        return result;
    }

    private String getWGXGmotifNucl(char WGXGnucl)
    {
        String result = "<font color=\"00FF00\">" + WGXGnucl + "</font>";
        return result;
    }

    private String getMutatedWGXGmotifNucl(char mutatedWGXGNucl)
    {
        String result = "<font color=\"FFFF00\">" + mutatedWGXGNucl + "</font>";
        return result;
    }

    private String getUnspecifiedNucl(char unspecifiedNucl)
    {
        String result = "<font color=\"00FF00\">" + unspecifiedNucl + "</font>";
        return result;
    }

    private String getUnidentifiedNucl(char unidentifiedNucl)
    {
        String result = "<font color=\"0000FF\">" + unidentifiedNucl + "</font>";
        return result;
    }

    private String markupNandXNucleotides(String region)
    {
        String markedUpRegion = "";
        for(int i = 0; i < region.length(); i++)
        {
            char currNucl = region.charAt(i);
            if(currNucl == 'n' || currNucl == 'x')
            {
                markedUpRegion = markedUpRegion + getUnspecifiedNucl(currNucl);
            } else
            {
                markedUpRegion = markedUpRegion + currNucl;
            }
        }

        return markedUpRegion;
    }

    private String markupRegion(String nucleotideRegion)
    {
        String markedUpRegion = "";
        for(int i = 0; i < nucleotideRegion.length(); i++)
        {
            char currNucl = nucleotideRegion.charAt(i);
            markedUpRegion = markedUpRegion + getUnidentifiedNucl(currNucl);
        }

        return markedUpRegion;
    }
}
