package COMMoDORE;

/**
 *
 * @author thinkens
 * on 4/4/2018
 * @since 0.1.0
 */
public class Atom {
    
    private double v, p, pe, ke, m;
    private double[] position = {0,0};
    
    
    /**
     * Instantiates a new Atom of type
     * @param atomInfo 
     */
    public Atom(AtomInfo atomInfo) {
        
        this.v = 0;
        this.p = 0;
        this.ke = 0;
        this.pe = 0;
        this.m = atomInfo.m;
    }
    
    public void setPosition(double[] pos) {
        
        this.position = pos;
    }
}
