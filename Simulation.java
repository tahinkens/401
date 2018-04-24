package COMMoDORE;

import static COMMoDORE.Main.*;

/**
 *
 * @author thinkens
 */
public class Simulation {
    
    int N;
    
    /**
     * Calculates the scalar kinetic energy of the system. Sums the
     * kinetic energy of each atom individually.
     * 
     * @param atoms all the atoms in the simulation
     * @return a double representing the total kinetic energy of the system
     */
    public static double kineticEnergy(Atom[] atoms) {
        
        double ke = 0;
        double v;
        AtomInfo atomType = atoms[0].atomType;
        
        for(Atom atom : atoms) {
            v = atom.getSpeed();
            ke += v*v * atomType.m;
        }
        return ke;
    }
    
    /**
     * Calculates the potential energy of an atomic system as a function of the
     * positions of the particles. More precisely, the interatomic potential of
     * each atom on every other atom is calculated using a standard 
     * Lennard-Jones potential. This energy can then be differentiated with 
     * respect to a particle's position to find the force on that particle.
     * 
     * @param atoms an array containing all of the atoms in the simulation
     * @return the potential energy 
     */
    protected static double potentialEnergy(Atom[] atoms) {
        
        double potentialEnergy = 0;
        AtomInfo species = atoms[0].atomType; //assume that every Atom in atoms is same species
        
        for(Atom atomi : atoms) {
            for(Atom atomj : atoms) {
                if(!atomj.equals(atomi)) //only sum contributions from different atoms
                    potentialEnergy += atomj.lennardJones(species.homopotentialWellDepth,species.equilibriumSeparation,atomj.getDistanceFromNeighbor(atomi));
            }
        }
        return 0.5 * potentialEnergy;
    }
    
    protected static double temperature(double ke, int N) {
        
        return (2d/3d) * (ke / (N * K_B));
    }
    
    
}
