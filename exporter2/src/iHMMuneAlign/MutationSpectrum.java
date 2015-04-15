package iHMMuneAlign;

import java.io.*;
import java.util.StringTokenizer;

public class MutationSpectrum
{
    class KeyPair
    {

        String key;
        double value;

        public KeyPair(String key, double value)
        {
            this.key = key;
            this.value = value;
        }
    }


    final String DIVISION_BY_ZERO = "#DIV/0!";
    final char NUCLEOTIDE_N = 'n';
    final double ONE_THIRD = 0.33333333333333331D;
    final double VALUE_NO_MATCH = -1D;
    final int UNIQUE_NUCLEOTIDE_COUNT = 4;
    final char UNIQUE_NUCLEOTIDE_ARRAY[] = {'a', 'c', 'g', 't'};
    KeyPair key_pairs[];
    int key_pairs_index;
    BufferedReader br;
    boolean DEBUGGING = false;
    //gap char
    final char GAP_NUCLEOTIDE = '-';

    public static void main(String args[])
    {
        MutationSpectrum ms = new MutationSpectrum("Mutation spectrum.txt");
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        String trinucleotide = null;
        String mutation = null;
        for(boolean finishedReading = false; !finishedReading;)
        {
            System.out.println("enter the trinucleotide (\"END\" TO EXIT):  ");
            trinucleotide = getInputLine(br);
            if(trinucleotide.equalsIgnoreCase("end"))
            {
                finishedReading = true;
            } else
            {
                System.out.println("enter the mutation:  ");
                mutation = getInputLine(br);
                System.out.println();
                System.out.println();
                System.out.println("Searching for match");
                double prob = ms.getTNProbability(trinucleotide, mutation.charAt(0));
                System.out.println();
                System.out.println("finished Searching for match");
                System.out.println("prob = " + prob);
                System.out.println();
            }
        }

    }

    private static String getInputLine(BufferedReader br)
    {
        String input = null;
        try
        {
            while(input == null) 
            {
                System.out.println("enter input followed by enter");
                input = br.readLine();
            }
        }
        catch(IOException ioe)
        {
            System.out.println("getInputLine(): could not read from keyboard");
            return null;
        }
        return input;
    }

    public MutationSpectrum(String mutationSpectrumFileName)
    {
        key_pairs = new KeyPair[192];
        key_pairs_index = 0;
        br = null;
        try
        {
            br = new BufferedReader(new FileReader(mutationSpectrumFileName));
        }
        catch(FileNotFoundException fnfe)
        {
            throw new Error("MutationSpectrum():  " + fnfe.getMessage());
        }
        String data = "";
        String temp = "";
        try
        {
            for(temp = br.readLine(); temp != null; temp = br.readLine())
            {
                data = data + temp + " ";
            }

        }
        catch(IOException ioe)
        {
            throw new Error("MutationSpectrum():  " + ioe.getMessage());
        }
        for(StringTokenizer st = new StringTokenizer(data); st.hasMoreTokens();)
        {
            String token = st.nextToken();
            String key_token = token;
            if(!st.hasMoreTokens())
            {
                throw new Error("MutationSpectrum constructor:  Key-Value-pair missing Value");
            }
            token = st.nextToken();
            double value_token;
            if(token.equalsIgnoreCase("#DIV/0!"))
            {
                value_token = 0.33333333333333331D;
            } else
            {
                try
                {
                    value_token = Double.parseDouble(token);
                }
                catch(NumberFormatException nfe)
                {
                    throw new Error("Number Format Exception: " + nfe.getMessage());
                }
            }
            key_pairs[key_pairs_index] = new KeyPair(key_token, value_token);
            key_pairs_index++;
        }

    }

    public double getTNProbability(String trinucleotide, char mutationResult)
    {
        //check for correct trinucleotide length
        if(trinucleotide.length() != 3)
        {
            throw new Error("getTNProbability(): trinucleotide length not equal to THREE");
        }
        //check that a mutation of the middle nt of the trinucleotide has actually occurred
        if(trinucleotide.charAt(1) == mutationResult)
        {
            if (DEBUGGING) {
            	System.out.println("Middle nucleotide has silent mutation, not determining triNt mutation prob");
            }
            return VALUE_NO_MATCH;
        }
        
        double resultProb = VALUE_NO_MATCH;
        
        //determine if the trinucleotide includes any 'n' nucleotides
        if(trinucleotide.charAt(0) == NUCLEOTIDE_N || trinucleotide.charAt(2) == NUCLEOTIDE_N) {
			//mid nt is an 'n', so cannot determine a mutation probability        
            if(trinucleotide.charAt(1) == NUCLEOTIDE_N)
            {
                if (DEBUGGING) {
                	System.out.println("Middle nt of trinucleotide is \"n\", cannot lookup triNt mutation probability");
                }
                //don't want to throw an error as this does occur in some germline genes, just return a 0 prob of mut
                return 0.0D;
            }
            //either the first or the last nts within the trint are 'n', use findEndTriNuclProb to find the average prob for the triNt if n = {a,g,c,t}
            resultProb = findEndTriNuclProb(trinucleotide, mutationResult);
        } else if ( (trinucleotide.charAt(0) == GAP_NUCLEOTIDE) || (trinucleotide.charAt(2) == GAP_NUCLEOTIDE)) {
        	//determine if there are any gap chars within the trinucleotide
        	
        	//if the gap char is at the central nt, no triNt mutation prob can be determined
        	if(trinucleotide.charAt(1) == GAP_NUCLEOTIDE)
            {
                if (DEBUGGING) {
                	System.out.println("Middle nt of trinucleotide is \"-\", cannot lookup triNt mutation probability");
                }
                //return a 0 mutation prob
                return 0.0D;
            }
            
            //if the gap char is at the first or last position of the trinucleotide can still determine a probability of mutation based on averaging all possible triNts
            resultProb = findGappedTriNuclProb(trinucleotide, mutationResult);
        	
        } else {
            //central nt is an 'n' or a '-' so can't calculation a probability of mutation for this nt
            if( (trinucleotide.charAt(1) == NUCLEOTIDE_N) || (trinucleotide.charAt(1) == GAP_NUCLEOTIDE) ){
               // the middle nt is an unknown nucleotide or a gap char, this does occur in some germline genes
               // set the prob for an n target to 0
               if (DEBUGGING) {
                	System.out.println("Middle nt of trinucleotide is \"-\" or \"n\", cannot lookup triNt mutation probability");
                }
               return 0.0D;
            } else {
            	//determine the mutation prob based on the trint
               resultProb = findTriNuclProb(trinucleotide, mutationResult);
            }
        }
        if(resultProb == VALUE_NO_MATCH)
        {
            if (DEBUGGING) {
            	System.out.println("NO match could be found (after lookup)");
            }
            return VALUE_NO_MATCH;
        } else
        {
            return resultProb;
        }
    }

    private double findTriNuclProb(String trinucleotide, char mutationResult)
    {
        String lookup_key = getLookup_Key(trinucleotide, mutationResult);
        double probability = getMatch(lookup_key);
        return probability;
    }

    private double findEndTriNuclProb(String generalTrinucleotide, char mutationResult)
    {
        double totalProbability = 0.0D;
        for(int i = 0; i < 4; i++)
        {
            char unique_nucl = UNIQUE_NUCLEOTIDE_ARRAY[i];
            String uniqueTrinucleotide = generalTrinucleotide.replace(NUCLEOTIDE_N, unique_nucl);
            String lookup_key = getLookup_Key(uniqueTrinucleotide, mutationResult);
            double curr_prob = getMatch(lookup_key);
            totalProbability += curr_prob;
        }

        double averageProbability = totalProbability / 4D;
        return averageProbability;
    }
    
    private double findGappedTriNuclProb(String generalTrinucleotide, char mutationResult)
    {
        double totalProbability = 0.0D;
        for(int i = 0; i < 4; i++)
        {
            char unique_nucl = UNIQUE_NUCLEOTIDE_ARRAY[i];
            String uniqueTrinucleotide = generalTrinucleotide.replace(GAP_NUCLEOTIDE, unique_nucl);
            String lookup_key = getLookup_Key(uniqueTrinucleotide, mutationResult);
            double curr_prob = getMatch(lookup_key);
            totalProbability += curr_prob;
        }

        double averageProbability = totalProbability / 4D;
        return averageProbability;
   }
     

    private String getLookup_Key(String trinucleotide, char mutationResult)
    {
        String lookup_key = trinucleotide.charAt(0) + "(" + trinucleotide.charAt(1) + "->" + mutationResult + ")" + trinucleotide.charAt(2);
        return lookup_key;
    }

    private double getMatch(String lookup_key)
    {
        for(int i = 0; i < key_pairs.length; i++)
        {
            KeyPair keyPair = key_pairs[i];
            if(keyPair.key.equalsIgnoreCase(lookup_key))
            {
                return keyPair.value;
            }
        }

        return VALUE_NO_MATCH;
    }
}
