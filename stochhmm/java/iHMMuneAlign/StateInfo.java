package iHMMuneAlign;


// Referenced classes of package iHMMuneAlign:
//            GlobalDefines

class StateInfo
    implements GlobalDefines
{

    String geneName;
    char geneToken;
    char preEmissionSymbol;
    int geneNumber;
    int nucleotidePosition;
    int completeGeneLength;

    public StateInfo(String geneName, char geneToken, int geneNumber, int nucleotidePosition, char preEmissionSymbol, int completeGeneLength)
    {
        this.geneName = geneName;
        this.geneToken = geneToken;
        this.geneNumber = geneNumber;
        this.nucleotidePosition = nucleotidePosition;
        this.preEmissionSymbol = preEmissionSymbol;
        this.completeGeneLength = completeGeneLength;
    }

    public StateInfo(String geneName, char geneToken)
    {
        this.geneName = geneName;
        this.geneToken = geneToken;
        if(geneToken == 'X')
        {
            preEmissionSymbol = ' ';
        } else
        {
            preEmissionSymbol = '?';
        }
    }

    public String getGeneName()
    {
        return geneName;
    }

    public char getGeneToken()
    {
        return geneToken;
    }

    public int getGeneNumber()
    {
        return geneNumber;
    }

    public int getNucleotidePosition()
    {
        return nucleotidePosition;
    }

    public char getPrepreEmissionSymbol()
    {
        return preEmissionSymbol;
    }
}
