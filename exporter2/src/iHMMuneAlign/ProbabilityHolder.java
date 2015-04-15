package iHMMuneAlign;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class ProbabilityHolder
{

    private static final String VD_N = "VD_N";
    private static final String DJ_N = "DJ_N";
    private static final String V_END_EXO_MEAN = "V_end_exo_mean";
    private static final String V_END_EXO_STDDEV = "V_end_exo_stdDev";
    private static final String D_START_EXO_MEAN = "D_start_exo_mean";
    private static final String D_START_EXO_STDDEV = "D_start_exo_stdDev";
    private static final String D_END_EXO_MEAN = "D_end_exo_mean";
    private static final String D_END_EXO_STDDEV = "D_end_exo_stdDev";
    private static final String J_START_EXO_MEAN = "J_start_exo_mean";
    private static final String J_START_EXO_STDDEV = "J_start_exo_stdDev";
    private static final String V_END_P = "V_end_P";
    private static final String D_START_P = "D_start_P";
    private static final String D_END_P = "D_end_P";
    private static final String J_START_P = "J_start_P";
    private static final String GENE_MUTATION = "Gene_Mutation";
    double VD_N_Fields[];
    double DJ_N_Fields[];
    double Kappa_N_Fields[];
    //added to allow IGL alignments
    double Lambda_N_Fields[];
    double V_end_exo_mean_Fields[];
    double V_end_exo_stdDev_Fields[];
    double D_start_exo_mean_Fields[];
    double D_start_exo_stdDev_Fields[];
    double D_end_exo_mean_Fields[];
    double D_end_exo_stdDev_Fields[];
    double J_start_exo_mean_Fields[];
    double J_start_exo_stdDev_Fields[];
    double Kappa_V_end_exo_Fields[];
    double Kappa_J_start_exo_Fields[];
    //added for IGL
    double Lambda_V_end_exo_Fields[];
    double Lambda_J_start_exo_Fields[];
    double V_end_P_Fields[];
    double D_start_P_Fields[];
    double D_end_P_Fields[];
    double J_start_P_Fields[];
    double Gene_Mutation_Field;
    //added to allow supression of output unless debugging
    static boolean DEBUGGING = false;

    public ProbabilityHolder()
    {
        VD_N_Fields = null;
        DJ_N_Fields = null;
        Kappa_N_Fields = null;
        Lambda_N_Fields = null;
        V_end_exo_mean_Fields = null;
        V_end_exo_stdDev_Fields = null;
        D_start_exo_mean_Fields = null;
        D_start_exo_stdDev_Fields = null;
        D_end_exo_mean_Fields = null;
        D_end_exo_stdDev_Fields = null;
        J_start_exo_mean_Fields = null;
        J_start_exo_stdDev_Fields = null;
        Kappa_V_end_exo_Fields = null;
        Kappa_J_start_exo_Fields = null;
        Lambda_V_end_exo_Fields = null;
        Lambda_J_start_exo_Fields = null;
        V_end_P_Fields = null;
        D_start_P_Fields = null;
        D_end_P_Fields = null;
        J_start_P_Fields = null;
        Gene_Mutation_Field = 0.0D;
    }

    public double[] getExoProbArray(double mean, double stdDev, double minProbLimit, boolean reverseArray)
    {
        ArrayList temp_result = new ArrayList();
        int i = 0;
        do
        {
            double prob = getNormalDistProb(i, mean, stdDev);
            Double temp_D = new Double(prob);
            if(prob < minProbLimit && (double)i > mean)
            {
                break;
            }
            temp_result.add(temp_D);
            i++;
        } while(true);
        double result_normalized[] = new double[temp_result.size()];
        double sum = 0.0D;
        for(i = 0; i < temp_result.size(); i++)
        {
            Double temp_D = (Double)temp_result.get(i);
            sum += temp_D.doubleValue();
        }

        if(reverseArray)
        {
            int array_index = 0;
            for(i = result_normalized.length - 1; i >= 0; i--)
            {
                Double temp_D = (Double)temp_result.get(i);
                result_normalized[array_index] = temp_D.doubleValue() / sum;
                array_index++;
            }

        } else
        {
            for(i = 0; i < result_normalized.length; i++)
            {
                Double temp_D = (Double)temp_result.get(i);
                result_normalized[i] = temp_D.doubleValue() / sum;
            }

        }
        sum = 0.0D;
        for(i = 0; i < result_normalized.length; i++)
        {
            sum += result_normalized[i];
        }

        return result_normalized;
    }

    private double getNormalDistProb(double x, double mean, double stdDev)
    {
        double temp1 = 1.0D / Math.sqrt(6.2831853071795862D * Math.pow(stdDev, 2D));
        double temp2 = Math.pow(2.7182818284590451D, -0.5D * Math.pow((x - mean) / stdDev, 2D));
        double result = temp1 * temp2;
        return result;
    }

    public static ProbabilityHolder createFromString(String string)
    {
        String END = "end";
        ProbabilityHolder result = new ProbabilityHolder();
        StringTokenizer st = new StringTokenizer(string);
        try
        {
            String token = st.nextToken();
            if(!token.equalsIgnoreCase("Probabilities"))
            {
                throw new Error();
            }
            while(st.hasMoreTokens()) 
            {
                token = st.nextToken();
                if(token.equals("VD_N")) {
                    result.VD_N_Fields = parseDoubleArray(st);
                } else if(token.equals("DJ_N")) {
                    result.DJ_N_Fields = parseDoubleArray(st);
                } else if(token.equals("Kappa_N")) {
                    result.Kappa_N_Fields = parseDoubleArray(st);
                } else if(token.equals("V_end_exo_mean")) {
                    result.V_end_exo_mean_Fields = parseDoubleArray(st);
                } else if(token.equals("V_end_exo_stdDev")) {
                    result.V_end_exo_stdDev_Fields = parseDoubleArray(st);
                } else if(token.equals("D_start_exo_mean")) {
                    result.D_start_exo_mean_Fields = parseDoubleArray(st);
                } else if(token.equals("D_start_exo_stdDev")) {
                    result.D_start_exo_stdDev_Fields = parseDoubleArray(st);
                } else if(token.equals("D_end_exo_mean")) {
                    result.D_end_exo_mean_Fields = parseDoubleArray(st);
                } else if(token.equals("D_end_exo_stdDev")) {
                    result.D_end_exo_stdDev_Fields = parseDoubleArray(st);
                } else if(token.equals("J_start_exo_mean")) {
                    result.J_start_exo_mean_Fields = parseDoubleArray(st);
                } else if(token.equals("J_start_exo_stdDev")) {
                    result.J_start_exo_stdDev_Fields = parseDoubleArray(st);
                } else if(token.equals("Kappa_V_end_exo")) {
                    result.Kappa_V_end_exo_Fields = parseDoubleArray(st);
                } else if(token.equals("Kappa_J_start_exo")) {
                    result.Kappa_J_start_exo_Fields = parseDoubleArray(st);
                } else if(token.equals("V_end_P")) {
                    result.V_end_P_Fields = parseDoubleArray(st);
                } else if(token.equals("D_start_P")) {
                    result.D_start_P_Fields = parseDoubleArray(st);
                } else if(token.equals("D_end_P")) {
                    result.D_end_P_Fields = parseDoubleArray(st);
                } else if(token.equals("J_start_P")) {
                    result.J_start_P_Fields = parseDoubleArray(st);
                } else if(token.equals("Gene_Mutation")) {
                    result.Gene_Mutation_Field = parseDouble(st);
                } else if (token.equals("Lambda_N")) {
                	result.Lambda_N_Fields = parseDoubleArray(st);
                } else if (token.equals("Lambda_V_end_exo")) {
                	result.Lambda_V_end_exo_Fields = parseDoubleArray(st);
                } else if (token.equals("Lambda_J_start_exo")) {
                	result.Lambda_J_start_exo_Fields = parseDoubleArray(st);
                }
                
            }
            if(result.VD_N_Fields == null)
            {
                if (DEBUGGING) {
                	System.out.println(" vd n 53 = null");
                }
            } else
            {
                if (DEBUGGING) {
                	System.out.println(" vd n 53 not equal null");
                }
            }
        }
        catch(Exception e)
        {
            System.out.println("parsing exception");
        }
        return result;
    }

    private static double parseDouble(StringTokenizer st)
    {
        double d = Double.parseDouble(st.nextToken());
        if (DEBUGGING) {
        	System.out.println("parse double = " + d);
        }
        return d;
    }

    private static double[] parseDoubleArray(StringTokenizer st)
    {
        ArrayList temp_result = new ArrayList();
        double result[] = (double[])null;
        String token = st.nextToken();
        if(token.equalsIgnoreCase("Array"))
        {
            do
            {
                token = st.nextToken();
                if(token.equals("end"))
                {
                    break;
                }
                token = st.nextToken();
                temp_result.add(new Double(Double.parseDouble(token)));
            } while(st.hasMoreTokens());
        } else
        {
            return null;
        }
        result = new double[temp_result.size()];
        for(int i = 0; i < result.length; i++)
        {
            Double temp_D = (Double)temp_result.get(i);
            double temp_d = temp_D.doubleValue();
            result[i] = temp_d;
        }

        return result;
    }

    public String toString()
    {
        String N = "\n";
        String END = "end";
        String result = new String("Probabilities" + N + "VD_N" + doubleArrayToString(VD_N_Fields) + END + N + "DJ_N" + doubleArrayToString(DJ_N_Fields) + END + N + "V_end_exo_mean" + doubleArrayToString(V_end_exo_mean_Fields) + END + N + "V_end_exo_stdDev" + doubleArrayToString(V_end_exo_stdDev_Fields) + END + N + "D_start_exo_mean" + doubleArrayToString(D_start_exo_mean_Fields) + END + N + "D_start_exo_stdDev" + doubleArrayToString(D_start_exo_stdDev_Fields) + END + N + "D_end_exo_mean" + doubleArrayToString(D_end_exo_mean_Fields) + END + N + "D_end_exo_stdDev" + doubleArrayToString(D_end_exo_stdDev_Fields) + END + N + "J_start_exo_mean" + doubleArrayToString(J_start_exo_mean_Fields) + END + N + "J_start_exo_stdDev" + doubleArrayToString(J_start_exo_stdDev_Fields) + END + N + "V_end_P" + doubleArrayToString(V_end_P_Fields) + END + N + "D_start_P" + doubleArrayToString(D_start_P_Fields) + END + N + "D_end_P" + doubleArrayToString(D_end_P_Fields) + END + N + "J_start_P" + doubleArrayToString(J_start_P_Fields) + END + N + "Gene_Mutation" + "  " + Gene_Mutation_Field + "  " + END + N);
        return result;
    }

    private String doubleArrayToString(double array[])
    {
        String result = "  Array ";
        for(int i = 0; i < array.length; i++)
        {
            String temp = "[" + i + "]  " + array[i] + "  ";
            result = result + temp;
        }

        return result;
    }
}
