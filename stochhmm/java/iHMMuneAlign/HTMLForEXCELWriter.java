package iHMMuneAlign;

import java.io.PrintWriter;
import java.util.Vector;

// Referenced classes of package iHMMuneAlign:
//            HTML_Writer, GeneInfo, RegionInfo

public class HTMLForEXCELWriter extends HTML_Writer
{

    final char N_NUCLEOTIDE = 'n';
    final char X_NUCLEOTIDE = 'x';

    public HTMLForEXCELWriter()
    {
    }

    public void createDocument(String fileName, String sequenceName, int chainType, GeneInfo VGene, RegionInfo P1, RegionInfo N1, RegionInfo P2, GeneInfo DGene, 
            RegionInfo P3, RegionInfo N2, RegionInfo P4, GeneInfo JGene, String selectedSequenceName, String selectedSequenceString, double score, byte alignmentType, Vector stopCodons, int relativeMotifPosition, boolean isJGeneInFrame, int reverseComplement, String dProbsFilename)
    {
        if(fostream == null)
        {
            createFile(fileName, alignmentType);
        }
        boolean acceptedDGene;
        if(DGene.acceptedAlignment)
        {
            acceptedDGene = true;
        } else
        {
            acceptedDGene = false;
        }
        String unidentifiedRegion = "";
        String markedUpDGeneString = null;
        String dGeneString = null;
        if(acceptedDGene)
        {
            markedUpDGeneString = markupMutations(DGene.aligned_gene_string, DGene.aligned_ums_string);
        } else
        {
            dGeneString = DGene.aligned_ums_string;
        }
        String VGeneName = VGene.gene_name;
        String DGeneName = DGene.gene_name;
        String JGeneName = JGene.gene_name;
        int p1Length = getRegionLength(P1);
        String p1String = getRegionString(P1);
        int p2Length = getRegionLength(P2);
        String p2String = getRegionString(P2);
        int n1Length = getRegionLength(N1);
        String n1String = getRegionString(N1);
        String pnpONEString = p1String + n1String + p2String;
        int p3Length = getRegionLength(P3);
        String p3String = getRegionString(P3);
        int p4Length = getRegionLength(P4);
        String p4String = getRegionString(P4);
        int n2Length = getRegionLength(N2);
        String n2String = getRegionString(N2);
        String pnpTWOString = p3String + n2String + p4String;
        int DGeneMutations = mutationsInGeneAlignment(DGene.aligned_gene_string, DGene.aligned_ums_string);
        String result = null;
        if(acceptedDGene)
        {
            result = "<table><tr>" + tableData(sequenceName) + tableData(VGeneName) + tableData(score) + tableData(p1Length) + tableData(n1Length) + tableData(p2Length) + tableData(pnpONEString) + tableData(DGeneName) + tableData(markedUpDGeneString) + tableData(DGeneMutations) + tableData(pnpTWOString) + tableData(p3Length) + tableData(n2Length) + tableData(p4Length) + tableData(JGeneName) + "</tr>" + "</table>";
        } else
        {
            unidentifiedRegion = pnpONEString + dGeneString + pnpTWOString;
            String markedUpUnidentifiedRegion = markupRegion(unidentifiedRegion);
            result = "<table><tr>" + tableData(sequenceName) + tableData(VGeneName) + tableData(score) + tableData(0.0D) + tableData(0.0D) + tableData(0.0D) + tableData("") + tableData("NO ACCEPTED D-GENE") + tableData(markedUpUnidentifiedRegion) + tableData("") + tableData("") + tableData(0.0D) + tableData(0.0D) + tableData(0.0D) + tableData(JGeneName) + "</tr>" + "</table>";
        }
        fostream.println(result);
    }
    
	public void createDocument(String fileName, String sequenceName, int chainType, GeneInfo VGene, RegionInfo N1, GeneInfo DGene, RegionInfo N2, GeneInfo JGene, String selectedSequenceName, String selectedSequenceString, double score, byte alignmentType, Vector stopCodons, int relativeMotifPosition, boolean isJGeneInFrame, int reverseComplement, String dProbsFilename)
	{
		if(fostream == null)
        {
            createFile(fileName, alignmentType);
        }
        boolean acceptedDGene;
        if(DGene.acceptedAlignment)
        {
            acceptedDGene = true;
        } else
        {
            acceptedDGene = false;
        }
        String unidentifiedRegion = "";
        String markedUpDGeneString = null;
        String dGeneString = null;
        if(acceptedDGene)
        {
            markedUpDGeneString = markupMutations(DGene.aligned_gene_string, DGene.aligned_ums_string);
        } else
        {
            dGeneString = DGene.aligned_ums_string;
        }
        String VGeneName = VGene.gene_name;
        String DGeneName = DGene.gene_name;
        String JGeneName = JGene.gene_name;
        int n1Length = getRegionLength(N1);
        String n1String = getRegionString(N1);
		int n2Length = getRegionLength(N2);
        String n2String = getRegionString(N2);
        int DGeneMutations = mutationsInGeneAlignment(DGene.aligned_gene_string, DGene.aligned_ums_string);
        String result = null;
        if(acceptedDGene)
        {
            result = "<table><tr>" + tableData(sequenceName) + tableData(VGeneName) + tableData(score) + tableData(n1Length) + tableData(n1String) + tableData(DGeneName) + tableData(markedUpDGeneString) + tableData(DGeneMutations) + tableData(n2String) + tableData(n2Length) + tableData(JGeneName) + "</tr>" + "</table>";
        } else
        {
            unidentifiedRegion = n1String + dGeneString + n2String;
            String markedUpUnidentifiedRegion = markupRegion(unidentifiedRegion);
            result = "<table><tr>" + tableData(sequenceName) + tableData(VGeneName) + tableData(score) + tableData(0.0D) + tableData("") + tableData("NO ACCEPTED D-GENE") + tableData(markedUpUnidentifiedRegion) + tableData("") + tableData("") + tableData(0.0D) + tableData(JGeneName) + "</tr>" + "</table>";
        }
        fostream.println(result);
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
      String VGeneName = "" + errorType;
      String blank = "NA";

      //create the output
      String result = null;
      result = "<table><tr>" + tableData(sequenceName) + tableData(VGeneName) + tableData(blank) + tableData(blank) + tableData(blank) + tableData(blank) + tableData(blank) + tableData(blank) + tableData(blank) + tableData(blank) + tableData(blank) + tableData(blank) + tableData(blank) + tableData(blank) + tableData(blank) + "</tr>" + "</table>";
      
      //write the output
      fostream.println(result);
    }

    private String tableData(String data)
    {
        String result = "<td>" + data + "</td>";
        return result;
    }

    private String tableData(double data)
    {
        String result = "<td>" + data + "</td>";
        return result;
    }

    private int getRegionLength(RegionInfo region)
    {
        if(region == null)
        {
            return 0;
        } else
        {
            return region.region_string.length();
        }
    }

    private String getRegionString(RegionInfo region)
    {
        if(region == null)
        {
            return "";
        } else
        {
            return region.region_string;
        }
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

    private String getMutatedNucl(char mutNucl)
    {
        String result = "<font color=\"FF0000\">" + mutNucl + "</font>";
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

    private String markupMutations(String geneNucleotideString, String sequenceNucleotideString)
    {
        String markedUpRegion = "";
        for(int i = 0; i < sequenceNucleotideString.length(); i++)
        {
            char currGeneNucl = geneNucleotideString.charAt(i);
            char currSequenceNucl = sequenceNucleotideString.charAt(i);
            if(currGeneNucl != currSequenceNucl)
            {
                markedUpRegion = markedUpRegion + getMutatedNucl(currSequenceNucl);
            } else
            {
                markedUpRegion = markedUpRegion + currSequenceNucl;
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
