package iHMMuneAlign;

import java.io.*;
import java.util.Vector;

// Referenced classes of package iHMMuneAlign:
//            GeneInfo, RegionInfo

public abstract class HTML_Writer
{

    final String VHEADER = "VGene (end)";
    final String DHEADER = "DGene";
    final String JHEADER = "JGene";
    final String N1HEADER = "N1";
    final String N2HEADER = "N2";
    final String P1HEADER = "P1";
    final String P2HEADER = "P2";
    final String P3HEADER = "P3";
    final String P4HEADER = "P4";
    final String TABLE_START = "<table>";
    final String TABLE_END = "</table>";
    final String TD = "<td>";
    final String TAB = "\t";
    final String QTAB = "\t\t\t\t";
    final String SPACE = " ";
    final String PARAGRAPH_START = "<pre>";
    final String PARAGRAPH_END = "</pre>";
    final char N_NUCLEOTIDE = 'n';
    final char X_NUCLEOTIDE = 'x';
    PrintWriter fostream;

    public HTML_Writer()
    {
        fostream = null;
    }

    public abstract void createDocument(String s, String s1, int chainType, GeneInfo geneinfo, RegionInfo regioninfo, RegionInfo regioninfo1, RegionInfo regioninfo2, GeneInfo geneinfo1, 
            RegionInfo regioninfo3, RegionInfo regioninfo4, RegionInfo regioninfo5, GeneInfo geneinfo2, String s2, String s3, double d, byte byte0, Vector vector, int i, boolean flag, int reverseComplement, String dProbsFilename);

    public abstract void  createErrDocument(String fileName, String sequenceName, byte alignmentType, String errorType, String vOutput);
	
    public abstract void createDocument(String s, String s1, int chainType, GeneInfo geneInfo, RegionInfo regioninfo, GeneInfo geneinfo1, RegionInfo regioninfo2, GeneInfo geneinfo2, String s2, String s3, double d, byte byte0, Vector vector, int i, boolean flag, int reverseComplement, String dProbsFilename);

    public void createFile(String fileName, byte alignmentType)
    {
        boolean appendToFile = false;
        if(alignmentType == 4)
        {
            appendToFile = true;
        }
        try
        {
            String temp = fileName.replace('|', '_');
            String acceptedFileName = temp;
            //fostream = new PrintWriter(new FileWriter(acceptedFileName + ".html", appendToFile));
            fostream = new PrintWriter(new FileWriter(acceptedFileName, appendToFile));
        }
        catch(Exception ex)
        {
            System.out.println("couldn't open new file: " + ex.getMessage());
        }
    }

    protected void closeHTMLFile()
    {
        System.out.println("HTML file closed in HTMLForExcelWriter");
        fostream.close();
        fostream = null;
    }
}
