package iHMMuneAlign;


public class GeneInfo
{

    int start_gene_pos;
    int end_gene_pos;
    int gene_length;
    String gene_name;
    String aligned_gene_string;
    String aligned_ums_string;
    boolean acceptedAlignment;

    public GeneInfo(String geneName, String geneString, String umsString, int geneLength, int startGenePos, int endGenePos)
    {
        acceptedAlignment = true;
        gene_name = geneName;
        aligned_gene_string = geneString;
        aligned_ums_string = umsString;
        gene_length = geneLength;
        start_gene_pos = startGenePos;
        end_gene_pos = endGenePos;
    }
}
