package iHMMuneAlign;

import java.io.PrintStream;
import java.util.Map;
import java.util.Set;
import org.biojava.bio.Annotation;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeType;
//import org.biojava.utils.Changeable;
//import org.biojava.utils.ChangeSupport;

// Referenced classes of package iHMMuneAlign:
//            StateInfo

public class FakeAnnotation implements Annotation
{

    StateInfo stateInfo = null;

    public FakeAnnotation(StateInfo stateInfo)
    {
        //this.stateInfo = null;
        this.stateInfo = stateInfo;
    }

    /*public static void main(String args[])
    {
        StateInfo siTest = new StateInfo("abc", 'D', 10, 8, 'e', -1);
        FakeAnnotation fa = new FakeAnnotation(siTest);
        if(!fa.containsProperty(null))
        {
            fa.setProperty(null, siTest);
        }
        StateInfo si = (StateInfo)fa.getProperty(null);
        System.out.println("si nucloetide number = " + si.getNucleotidePosition());
    }*/

    public Map asMap()
    {
        return null;
    }

    public boolean containsProperty(Object key)
    {
        return stateInfo != null;
    }

    public Object getProperty(Object key)
    {
        return stateInfo;
    }

    public Set keys()
    {
        return null;
    }

    public void removeProperty(Object key)
    {
        stateInfo = null;
    }

    public void setProperty(Object key, Object value)
    {
        stateInfo = (StateInfo)value;
   }

    public void addChangeListener(ChangeListener changelistener){
		//DEPRECATED
		addChangeListener(changelistener, ChangeType.UNKNOWN);
	}

    public void addChangeListener(ChangeListener changelistener, ChangeType changetype) {}

    public void removeChangeListener(ChangeListener changelistener){
		//DEPRECATED
		removeChangeListener(changelistener, ChangeType.UNKNOWN);
	}

    public void removeChangeListener(ChangeListener changelistener, ChangeType changetype){}

    public boolean isUnchanging(ChangeType ct)
    {
        return true;
    }

}
