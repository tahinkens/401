package COMMoDORE;

/**
 * Enumeration AtomInfo contains basic properties of atoms such as atomic number,
 * particle mass, equilibrium separation, etc. 
 * 
 * @author thinkens
 * on 4/4/2018
 * @since 0.1.0
 */
public enum AtomInfo {
    
    HELIUM(2,6.646632348e-27,-1,-1,-1,-1,-1),
    ALUMINUM(13,4.48039e-27,1.23e-10,368,4.04e-10,4.04e-10,4.04e-10), //epsilon/m (eq. sep from Totten and Mackenzie)
    COPPER(29,1.0552308203254732e-25,-1,-1,-1,-1,-1),
    XENON(54,2.1802e-25,-1,-1,-1,-1,-1);
    
    public final int Z; //atomic number
    public final double m; //mass per atom/kg
    public final double epsilon; //equilibrium interatomic separation/m
    public final double homopotentialWellDepth; //depth of potential well for two Atoms of same species/J-mol^-1
    public final double a, b, c; //lattice parameters/m
    
    private AtomInfo(int Z, double m, double epsilon, double homopotentialWellDepth, double a, double b, double c) {
        
        this.Z = Z;
        this.m = m;
        this.epsilon = epsilon;
        this.homopotentialWellDepth = homopotentialWellDepth;
        this.a = a;
        this.b = b;
        this.c = c;
    }
    
    
}
