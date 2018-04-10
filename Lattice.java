package COMMoDORE;

/**
 *
 * @author thinkens
 * on 4/4/2018
 * @since 0.1.0
 */
public class Lattice {
    
    private Atom[] inhabitants; 
    
    
    /**
     * 
     * @param latticeType desired Bravais lattice (SimpleCubic, FCC, BCC, HCP)
     * @param N desired number of atoms to simulate
     */
    public Lattice(CrystalLattice latticeType, int N) {
        
        inhabitants = new Atom[N];
        inhabitants = initializePositions(latticeType);
    }
    
    /**
     * 
     * @return an array of atoms with all of their position set to the correct
     * spots on the desired crystal lattice
     */
    private Atom[] initializePositions(CrystalLattice latticeType) {
        
        double latticeWidth;
        double latticeHeight;
        
        for(Atom atom : inhabitants) {
            switch(latticeType) {
                case SimpleCubic:
                    
                    break;
                case FCC:
                    
                    break;
                case BCC:
                    
                    break;
                case HCP:
                    
                    break;
            }
        }
        return inhabitants;
    }
    
    public Atom[] getInhabitants() {
        
        return inhabitants;
    }
}
