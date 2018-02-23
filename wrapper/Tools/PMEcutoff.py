#Script to reweight the reaction field energy with PME, to correct for the presence
#of electrostatic finite size artefacts
import os,sys, random
import math
from Sire.Tools.OpenMMMD import *
from Sire.Tools import Parameter, resolveParameters


# Python dependencies
#
try:
    mdtraj = Sire.try_import("mdtraj")
except:
    pass

try:
    numpy = Sire.try_import("numpy")
except:
    pass


short_cutoff_dist = Parameter("short cutoff distance", 10 * angstrom,
                        """The short cutoff distance to use for the non-bonded interactions.""")

long_cutoff_dist = Parameter("long cutoff distance", 10 * angstrom,
                        """The long cutoff distance to use for the non-bonded interactions.""")


use_solute_inter_nrg = Parameter("use solute intermolecular energy",False,"""Whether to use the solute intermolecular energy instead of the total electrostatic potential energy.""")

trajfile = Parameter("trajfile", "traj000000001.dcd",
                    """File name of the trajectory to process.""")
                    
stepframe = Parameter("step_frame",1,"""The number of frames to step to between two succcessive evaluations.""")

skip_frames = Parameter("skip_frames",[ [1,2], [3,4] ],"""List of intervals of frames to skip.""")



def setupRFFF(system, space, cutoff=10* angstrom, cutoff_type="cutoffperiodic"):

    print ("Creating force fields... ")

    solutes = system[MGName("solutes")]
    solute = system[MGName("solute_ref")]
    solute_hard = system[MGName("solute_ref_hard")]
    solute_todummy = system[MGName("solute_ref_todummy")]
    solute_fromdummy = system[MGName("solute_ref_fromdummy")]

    solvent = system[MGName("solvent")]

    #Solvent intramolecular CLJ energy
    solvent_intraclj = IntraCLJFF("solvent_intralj")
    solvent_intraclj.add(solvent)
    if (cutoff_type != "nocutoff"):
        solvent_intraclj.setUseReactionField(True)
        solvent_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solvent_intraclj.setDisableReactionFieldShift(disable_crf.val)

    #import pdb; pdb.set_trace()

    # Solvent-solvent LJ energy
    solventff = InterCLJFF("solvent:solvent")
    solventff.add(solvent)
    if (cutoff_type != "nocutoff"):
        solventff.setUseReactionField(True)
        solventff.setReactionFieldDielectric(rf_dielectric.val)
        solventff.setDisableReactionFieldShift(disable_crf.val)

    # Solute intramolecular LJ energy
    solute_hard_intraclj = IntraCLJFF("solute_hard_intralj")
    solute_hard_intraclj.add(solute_hard)
    if (cutoff_type != "nocutoff"):
        solute_hard_intraclj.setUseReactionField(True)
        solute_hard_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_hard_intraclj.setDisableReactionFieldShift(disable_crf.val)

    solute_todummy_intraclj = IntraSoftCLJFF("solute_todummy_intralj")
    #solute_todummy_intraclj.setShiftDelta(shift_delta.val)
    #solute_todummy_intraclj.setCoulombPower(coulomb_power.val)
    solute_todummy_intraclj.add(solute_todummy)
    if (cutoff_type != "nocutoff"):
        solute_todummy_intraclj.setUseReactionField(True)
        solute_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_todummy_intraclj.setDisableReactionFieldShift(disable_crf.val)

    solute_fromdummy_intraclj = IntraSoftCLJFF("solute_fromdummy_intralj")
    #solute_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    #solute_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    solute_fromdummy_intraclj.add(solute_fromdummy)
    if (cutoff_type != "nocutoff"):
        solute_fromdummy_intraclj.setUseReactionField(True)
        solute_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_fromdummy_intraclj.setDisableReactionFieldShift(disable_crf.val)


    solute_hard_todummy_intraclj = IntraGroupSoftCLJFF("solute_hard:todummy_intralj")
    #solute_hard_todummy_intraclj.setShiftDelta(shift_delta.val)
    #solute_hard_todummy_intraclj.setCoulombPower(coulomb_power.val)
    solute_hard_todummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_todummy_intraclj.add(solute_todummy, MGIdx(1))
    if (cutoff_type != "nocutoff"):
        solute_hard_todummy_intraclj.setUseReactionField(True)
        solute_hard_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_hard_todummy_intraclj.setDisableReactionFieldShift(disable_crf.val)

    solute_hard_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_hard:fromdummy_intralj")
    #solute_hard_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    #solute_hard_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    solute_hard_fromdummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))
    if (cutoff_type != "nocutoff"):
        solute_hard_fromdummy_intraclj.setUseReactionField(True)
        solute_hard_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_hard_fromdummy_intraclj.setDisableReactionFieldShift(disable_crf.val)

    solute_todummy_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_todummy:fromdummy_intralj")
    #solute_todummy_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    #solute_todummy_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    solute_todummy_fromdummy_intraclj.add(solute_todummy, MGIdx(0))
    solute_todummy_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))
    if (cutoff_type != "nocutoff"):
        solute_todummy_fromdummy_intraclj.setUseReactionField(True)
        solute_todummy_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_todummy_fromdummy_intraclj.setDisableReactionFieldShift(disable_crf.val)

    # Solute-solvent LJ energy
    solute_hard_solventff = InterGroupCLJFF("solute_hard:solvent")
    solute_hard_solventff.add(solute_hard, MGIdx(0))
    solute_hard_solventff.add(solvent, MGIdx(1))
    if (cutoff_type != "nocutoff"):
        solute_hard_solventff.setUseReactionField(True)
        solute_hard_solventff.setReactionFieldDielectric(rf_dielectric.val)
        solute_hard_solventff.setDisableReactionFieldShift(disable_crf.val)

    solute_todummy_solventff = InterGroupSoftCLJFF("solute_todummy:solvent")
    solute_todummy_solventff.add(solute_todummy, MGIdx(0))
    solute_todummy_solventff.add(solvent, MGIdx(1))
    if (cutoff_type != "nocutoff"):
        solute_todummy_solventff.setUseReactionField(True)
        solute_todummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
        solute_todummy_solventff.setDisableReactionFieldShift(disable_crf.val)

    solute_fromdummy_solventff = InterGroupSoftCLJFF("solute_fromdummy:solvent")
    solute_fromdummy_solventff.add(solute_fromdummy, MGIdx(0))
    solute_fromdummy_solventff.add(solvent, MGIdx(1))
    if (cutoff_type != "nocutoff"):
        solute_fromdummy_solventff.setUseReactionField(True)
        solute_fromdummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
        solute_fromdummy_solventff.setDisableReactionFieldShift(disable_crf.val)

    # TOTAL
    forcefields = [solute_hard_intraclj, solute_todummy_intraclj, solute_fromdummy_intraclj,
                   solute_hard_todummy_intraclj, solute_hard_fromdummy_intraclj,
                   solute_todummy_fromdummy_intraclj,
                   solventff, solvent_intraclj,
                   solute_hard_solventff, solute_todummy_solventff, solute_fromdummy_solventff]    

    if  (use_solute_inter_nrg.val):
        system.add(solute_hard_solventff)
        system.add(solute_todummy_solventff)
        system.add(solute_fromdummy_solventff)
    else:
        for forcefield in forcefields:
            system.add(forcefield)

    #import pdb; pdb.set_trace()

    system.setProperty("space", space)
    system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff))
    system.setProperty("combiningRules", VariantProperty(combining_rules.val))
    system.setProperty("coulombPower", VariantProperty(coulomb_power.val))
    system.setProperty("shiftDelta", VariantProperty(shift_delta.val))


    total_nrg = solute_hard_intraclj.components().coulomb() + \
                solute_todummy_intraclj.components().coulomb(0) + solute_fromdummy_intraclj.components().coulomb(0) + \
                solute_hard_todummy_intraclj.components().coulomb(0) + solute_hard_fromdummy_intraclj.components().coulomb(0) + \
                solute_todummy_fromdummy_intraclj.components().coulomb(0) + \
                solventff.components().coulomb() + \
                solvent_intraclj.components().coulomb() + \
                solute_hard_solventff.components().coulomb() + \
                solute_todummy_solventff.components().coulomb(0) + \
                solute_fromdummy_solventff.components().coulomb(0)



    e_total = system.totalComponent()

    lam = Symbol("lambda")
    system.setComponent(e_total, total_nrg)
    system.setConstant(lam, 0.0)
    system.add(PerturbationConstraint(solutes))

    # NON BONDED Alpha constraints for the soft force fields
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy:solvent"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy:solvent"), 1 - lam))
    if (not use_solute_inter_nrg.val):
        system.add(PropertyConstraint("alpha0", FFName("solute_todummy_intralj"), lam))
        system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy_intralj"), 1 - lam))
        system.add(PropertyConstraint("alpha0", FFName("solute_hard:todummy_intralj"), lam))
        system.add(PropertyConstraint("alpha0", FFName("solute_hard:fromdummy_intralj"), 1 - lam))
        system.add(PropertyConstraint("alpha0", FFName("solute_todummy:fromdummy_intralj"), Max(lam, 1 - lam)))




    system.setComponent(lam, lambda_val.val)

    return system


def setupPME(system, space, cutoff=10* angstrom ,cutoff_type="pme"):

    print("Creating force fields... ")

    all = system[MGName("all")]
    molecules = system[MGName("molecules")]
    ions = system[MGName("ions")]

    # - first solvent-solvent coulomb/LJ (CLJ) energy
    internonbondedff = InterCLJFF("molecules:molecules")
    internonbondedff.setToleranceEwaldPME(0.00001)
    internonbondedff.add(molecules)

    inter_ions_nonbondedff = InterCLJFF("ions:ions")
    inter_ions_nonbondedff.setToleranceEwaldPME(0.00001)
    inter_ions_nonbondedff.add(ions)

    inter_ions_molecules_nonbondedff = InterGroupCLJFF("ions:molecules")
    inter_ions_molecules_nonbondedff.setToleranceEwaldPME(0.00001)
    inter_ions_molecules_nonbondedff.add(ions, MGIdx(0))
    inter_ions_molecules_nonbondedff.add(molecules, MGIdx(1))

    # Now solute bond, angle, dihedral energy
    intrabondedff = InternalFF("molecules-intrabonded")
    intrabondedff.setToleranceEwaldPME(0.00001)
    intrabondedff.add(molecules)

    # Now solute intramolecular CLJ energy
    intranonbondedff = IntraCLJFF("molecules-intranonbonded")
    intranonbondedff.setToleranceEwaldPME(0.00001)
    intranonbondedff.add(molecules)


    # Here is the list of all forcefields
    forcefields = [internonbondedff, intrabondedff, intranonbondedff, inter_ions_nonbondedff,
                   inter_ions_molecules_nonbondedff, restraintff]

    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff))
    system.setProperty("combiningRules", VariantProperty(combining_rules.val))

    total_nrg = internonbondedff.components().total() + \
                intranonbondedff.components().total() + intrabondedff.components().total() + \
                inter_ions_nonbondedff.components().total() + inter_ions_molecules_nonbondedff.components().total() + \
                restraintff.components().total()

    e_total = system.totalComponent()

    system.setComponent(e_total, total_nrg)

    # Add a monitor that calculates the average total energy and average energy
    # deltas - we will collect both a mean average and an zwanzig average
    system.add("total_energy", MonitorComponent(e_total, Average()))

    return system


def updateSystemfromTraj(system, frame_xyz, cell_lengths, cell_angles):
    #print("Here we are processing traj")
    traj_coordinates = frame_xyz[0]

    traj_box_x = cell_lengths[0][0].tolist()
    traj_box_y = cell_lengths[0][1].tolist()
    traj_box_z = cell_lengths[0][2].tolist()

    traj_natoms = len(traj_coordinates)

    # Sire does not support non rectangular boxes
    newmols_coords = {}

    traj_index = 0
    mol_index = 0

    molnums = system.molNums()
    molnums.sort()
    for molnum in molnums:
        #print(molnum)
        mol = system.molecule(molnum)[0].molecule() #-1 to take all the solute atoms
        #print(mol)
        molatoms = mol.atoms()
        molnatoms = mol.nAtoms()
        # Create an empty coord group using molecule so we get the correct layout
        newmol_coords = AtomCoords( mol.property("coordinates") )
        for x in range(0,molnatoms):
            tmparray = traj_coordinates[traj_index]
            atom_coord = Vector( tmparray[0].tolist() , tmparray[1].tolist() , tmparray[2].tolist() )
            atom = molatoms[x]
            cgatomidx = atom.cgAtomIdx()
            newmol_coords.set( cgatomidx, atom_coord)
            traj_index += 1

        newmols_coords[molnum] = newmol_coords
        mol_index += 1

    if traj_natoms != traj_index:
        print ("The number of atoms in the system is not equal to the number of atoms in the trajectory file ! Aborting.")
        sys.exit(-1)
        
    changedmols = MoleculeGroup("changedmols")
    mol_index = 0
    for molnum in molnums:
        mol = system.molecule(molnum)[0].molecule()
        newmol_coords = newmols_coords[molnum]
        mol = mol.edit().setProperty("coordinates", newmol_coords).commit()
        changedmols.add(mol)
        
    system.update(changedmols)
    space = PeriodicBox(Vector( traj_box_x, traj_box_y, traj_box_z ) )
    system.setProperty("space",space)

    return system

def getFreeEnergy(delta_nrgs):
    free_nrg = FreeEnergyAverage(temperature.val)
    for nrg in delta_nrgs:
        free_nrg.accumulate(nrg.value())
    deltaG = free_nrg.average() * kcal_per_mol
    return deltaG

def getMean(delta_nrgs):
    mean_nrg = Average()
    for nrg in delta_nrgs:
        mean_nrg.accumulate(nrg.value())
    deltaU = mean_nrg.average() * kcal_per_mol
    return deltaU

def resample(values):
    nvals = len(values)
    new_values = []
    for x in range(0,nvals):
        i = random.randint(0,nvals-1)
        new_values.append(values[i])
    return new_values

def switch(R,Rc,Rs):
    if (R>Rc):
        return 0
    elif (R<Rs):
        return 1
    else:
        x=(R-Rs)/(Rc-Rs)
        S = 1 - 6*math.pow(x,5)+15*math.pow(x,4)-10*math.pow(x,3)
        return S

def electrostatic(system):
    'Experiments with messing up electrostatics Returns the Barker-Watts switched energy of the solute with the rest of the system'

    solute = system[MGName("solute_ref")]
    solventNums = system[MGName("solvent")].molNums()
    sol_atoms = solute.first().molecule().atoms()

    space = system.property("space")

    Rcoul = system.property("switchingFunction").electrostaticCutoffDistance().value()
    Rcoul2 = Rcoul**2
    Rs = Rcoul - 0.5

    Rcut = 10.0
    Rcut2 = Rcut**2

    # slow !!!
    krf = (1/Rcoul**3) * ( ( rf_dielectric.val -1 ) / ( 2 * rf_dielectric.val +1 ) )
    crf = (1/Rcoul) * ( (3 * rf_dielectric.val)/(2*rf_dielectric.val+1))
    #crf = 0.0
    cnrg = 0.0
    for sol_atom in sol_atoms:
        print (sol_atom)
        ri = sol_atom.property("coordinates")
        qi = sol_atom.property("charge").value()
        inner_charge = qi
        for nummol in solventNums:
            slv = system.molecule(nummol)[0]
            for mol_atoms in slv.atoms():
                for mol_atom in mol_atoms:
                    rj = mol_atom.property("coordinates")
                    dij2 = space.calcDist2(ri,rj)
                    if dij2 > Rcut2:
                        continue
                    #print ("WITHIN CUTOFF")
                    qj = mol_atom.property("charge").value()
                    inner_charge += qj
                    dij = math.sqrt(dij2)
                    nrg = qi*qj*(1/dij+krf*dij**2-crf)
                    cnrg += nrg
                    #cnrg += nrg*switch(dij,Rc,Rs)
        print ("Inner charge: %8.5f" % (inner_charge))
        #ccor = (-inner_charge*qi)*(1/Rc+krf*Rc**2)
        #cnrg += ccor
    cnrg *= one_over_four_pi_eps0 * kcal_per_mol
    print ("cnrg: %s " % (cnrg))
    #import pdb; pdb.set_trace()
    return cnrg

def compute_histogram(delta_nrgs):
    values = []
    for delta in delta_nrgs:
        val = delta.value()
        values.append(val)
    hist,bins = numpy.histogram(values,bins=50)
    width=0.7*(bins[1]-bins[0])
    center = (bins[:-1] + bins[1:])/2
    #print("~~~~~~~~~~~~~~~~~~~~~~~len~~~~~~~~~~~~~~~~~~~~~~~")
    #print(len(hist),len(bins))
    #print("~~~~~~~~~~~~~~~~~~~~~~~BINS~~~~~~~~~~~~~~~~~~~~~~~")
    #print(bins)
    fig,ax = plt.subplots()
    ax.bar(center,hist,align="center",width=width)
    #x_axis_name = str(key) + " angle"
    #plt.xlabel(x_axis_name)
    #plt.ylabel("Occurrence")
    name ="histogram.png"
    plt.savefig(name)
    plt.clf()
    plt.close("all")
    #save bins and hist
    ofile = open("histogram.dat","w")
    counter= 0 
    for i,bin in enumerate(bins,0):
        #write frequency-hist and bin val
        if counter==len(bins)-1:
            break
        else:
            ofile.write("%.4f,%4f\n" % (hist[i],bin))
            counter+=1
    ofile.close()


@resolveParameters
def runLambda():
    try:
        host = os.environ['HOSTNAME']
    except KeyError:
        host = "unknown"
    print("### Running RF tail correction calculation on %s ###" % host)
    if verbose.val:
        print("###================= Simulation Parameters=====================###")
        Parameter.printAll()
        print ("###===========================================================###\n")
    print("lambda is %s" % lambda_val.val)

    if os.path.exists(s3file.val):
        (molecules, space) = Sire.Stream.load(s3file.val)
    else:
        amber = Amber()
        (molecules, space) = amber.readCrdTop(crdfile.val, topfile.val)
        Sire.Stream.save((molecules, space), s3file.val)

    #create the main system
    system = createSystemFreeEnergy(molecules)
    #now create a system with "short" cutoff, which is the Reaction Field
    #system
    system_shortc = System()
    system_shortc.copy(system)
    # Construct a forcefield to compute electrostatic energies
    system_shortc = setupRFFF(system_shortc, space, \
                              cutoff=short_cutoff_dist.val)
   

    #now create an MD system, which supports PME and use the "long" cutoff
    system_longc = createSystem(molecules)
    # Update this to setupRFFF
    system_longc = setupPME(system_longc, space, \
                             cutoff=long_cutoff_dist.val)

    sys.exit(-1)
    start_frame = 1
    end_frame = 1000000000
    step_frame = stepframe.val

    mdtraj_trajfile = mdtraj.open(trajfile.val,'r')
    nframes = len(mdtraj_trajfile)
    if end_frame > (nframes - 1):
        end_frame = nframes - 1
    mdtraj_trajfile.seek(start_frame)
    current_frame = start_frame


    # Now loop over snapshots in dcd and accumulate energies
    delta_nrgs = []
    ofile = open("electrostatics.dat","w")
    en_file = open("energies.dat","w")
    while (current_frame <= end_frame):
        #check if the current_frame is to be skip
        for interval in skip_frames.val:
            if (current_frame > interval[0] and current_frame < interval[1]):
                #create a variable that is true once  the condition is hold
                convergence = True
                break
            else:
                convergence = False
            
        #once we have cycled through all the intervals check the convergence
        if convergence:
            print("Discarding frame %s" % current_frame)
            current_frame += step_frame
            mdtraj_trajfile.seek(current_frame)
            print("Updating frame to %s" % current_frame)
            
        else:
            print ("Processing frame %s " % current_frame)
            #print ("CURRENT POSITION %s " % mdtraj_trajfile.tell() )
            frames_xyz, cell_lengths, cell_angles = mdtraj_trajfile.read(n_frames=1)
            system_shortc = updateSystemfromTraj(system_shortc, frames_xyz, cell_lengths, cell_angles)
            system_longc = updateSystemfromTraj(system_longc, frames_xyz, cell_lengths, cell_angles)

            delta_nrg = (system_longc.energy() - system_shortc.energy())

            delta_nrgs.append(delta_nrg)
            en_file.write("%d,%.4f\n" % (current_frame, delta_nrg.value()))
            current_frame += step_frame
            mdtraj_trajfile.seek(current_frame)
        #import pdb; pdb.set_trace()
    #print (delta_nrgs)
    #compute the histogram
    #compute_histogram(delta_nrgs)
    # Now compute free energy change
    deltaG = getFreeEnergy(delta_nrgs)
    #print (deltaG)
    nbootstrap = 100
    deltaG_bootstrap = np.zeros(nbootstrap)
    for x in range(0,nbootstrap):
        resampled_nrgs = resample(delta_nrgs)
        dG = getFreeEnergy(resampled_nrgs)
        deltaG_bootstrap[x] = dG.value()
    devG = deltaG_bootstrap.std()
    print ("Rcshort %.2f --> Rclong %.2f Angstrom , DG_ELEC = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (
        short_cutoff_dist.val.value(),long_cutoff_dist.val.value(),deltaG.value(), devG))
    deltaU = getMean(delta_nrgs)
    deltaU_bootstrap = np.zeros(nbootstrap)
    for x in range(0,nbootstrap):
        resampled_nrgs = resample(delta_nrgs)
        dU = getMean(resampled_nrgs)
        deltaU_bootstrap[x] = dU.value()
    dev = deltaU_bootstrap.std()
    print ("deltaU %8.5f +/- %8.5f kcal/mol (1 sigma)" % (deltaU.value(),dev))
    print((short_cutoff_dist.val.value(),long_cutoff_dist.val.value(),deltaG.value(), devG,deltaU.value(),dev))
    ofile.write("%.2f,%.2f,%8.5f,%8.5f,%8.5f,%8.5f\n" %(short_cutoff_dist.val.value(),long_cutoff_dist.val.value(),deltaG.value(), devG,deltaU.value(),dev))
    ofile.close()
