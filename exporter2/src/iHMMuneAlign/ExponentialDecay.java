package iHMMuneAlign;

import java.io.PrintStream;

public class ExponentialDecay
{

    static final double EXPONENTIAL_DECAY_RATE = -0.0023999999999999998D;
    static final int average_DGene_length = 24;

    public ExponentialDecay()
    {
    }

    public static void main(String args[])
    {
        double DGene_nucl_prob = exponentialDecayDGene("acgacgacgacg", 10, 100);
        System.out.println("Dgene nucl prob = " + DGene_nucl_prob);
        double JGene_nucl_prob = exponentialDecayJGene("acgacgacgacg", 10, 100);
        System.out.println("Jgene nucl prob = " + JGene_nucl_prob);
    }

    public static double exponentialDecayVGene(String VGene, int VGene_nucl_position, int VGene_start_offset)
    {
        int sequence_position = VGene_nucl_position + VGene_start_offset;
        double exponential_decay_VGene_prob = exponentialDecay(sequence_position);
        return exponential_decay_VGene_prob;
    }

    public static double exponentialDecayDGene(String DGene, int DGene_nucl_position, int completeVgeneLength)
    {
        int sequence_position = completeVgeneLength + DGene_nucl_position;
        double exponential_decay_DGene_prob = exponentialDecay(sequence_position);
        return exponential_decay_DGene_prob;
    }

    public static double exponentialDecayJGene(String JGene, int JGene_nucl_position, int completeVgeneLength)
    {
        int sequence_position = completeVgeneLength + 24 + JGene_nucl_position;
        double exponential_decay_JGene_prob = exponentialDecay(sequence_position);
        return exponential_decay_JGene_prob;
    }

    private static double exponentialDecay(int position_k)
    {
        double result = Math.exp(-0.0023999999999999998D * (double)position_k);
        return result;
    }
}
