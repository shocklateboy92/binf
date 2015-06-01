package iHMMuneAlign;

import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeVetoException;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by fernie on 10/08/14.
 */
public class DebugMarkovModel extends SimpleMarkovModel {

    List<State> stateList;
    Map<State, List<State>> transitionMap;

    public DebugMarkovModel(int i, FiniteAlphabet dna, String vdjGene) {
        super(i, dna, vdjGene);

        this.stateList = new ArrayList<State>();
    }

    @Override
    public void addState(State toAdd) throws IllegalSymbolException, ChangeVetoException {
        stateList.add(toAdd);
        super.addState(toAdd);
    }

    @Override
    public void createTransition(State from, State to) throws IllegalSymbolException, ChangeVetoException {
        super.createTransition(from, to);
    }
}
