import openawsem
import open3SPN2

from functools import partial

from openawsem.functionTerms import *
from openawsem.helperFunctions.myFunctions import *

#Location of the AWSEM information folder, including fragment memories
AWSEM_folder = "/path/to/awsem/directory/for/protein/part"

#File of fragment memory to be used
fragment = f"single_frags.mem"

#Native (or other reference) file
reference = f"{AWSEM_folder}/crystal_structure-openmmawsem.pdb"

def set_up_forces(s,protein, dna, computeQ, AWSEM = AWSEM_folder, fragment = fragment):
    # apply forces
    forces = {}

    for i in range(s.getNumForces()):
        force = s.getForce(i)
        force_name="CMMotionRemover"

    #Add 3SPN2 forces
    for force_name in open3SPN2.forces:
        print(force_name)
        force = open3SPN2.forces[force_name](dna)
        if force_name in ['BasePair','CrossStacking']:
            force.addForce(s)
        else:
            s.addForce(force)
        forces.update({force_name:force})

    #Add AWSEM forces. Fragment memories are in the protein residue only AWSEM-created folder
    frags_dir = AWSEM
    openAWSEMforces = dict(Connectivity=openawsem.functionTerms.basicTerms.con_term,
                        Chain=openawsem.functionTerms.basicTerms.chain_term,
                        Chi=openawsem.functionTerms.basicTerms.chi_term,
                        Excl=openawsem.functionTerms.basicTerms.excl_term,
                        rama=openawsem.functionTerms.basicTerms.rama_term,
                        rama_pro=openawsem.functionTerms.basicTerms.rama_proline_term,
                        rama_ssweight = partial(openawsem.functionTerms.basicTerms.rama_ssweight_term,
                                        k_rama_ssweight=2*8.368,
                                        ssweight_file=f"{frags_dir}/ssweight"),
                        contact=openawsem.functionTerms.contactTerms.contact_term,
                        frag  = partial(openawsem.functionTerms.templateTerms.fragment_memory_term, 
                                        frag_file_list_file=f"{frags_dir}/{fragment}", 
                                        UseSavedFragTable=False, 
                                        k_fm=0.04184),
                        beta1 = openawsem.functionTerms.hydrogenBondTerms.beta_term_1,
                        beta2 = openawsem.functionTerms.hydrogenBondTerms.beta_term_2,
                        beta3 = openawsem.functionTerms.hydrogenBondTerms.beta_term_3,
                        pap1 = partial(openawsem.functionTerms.hydrogenBondTerms.pap_term_1,
                                        ssweightFileName=f"{frags_dir}/ssweight"),
                        pap2 = partial(openawsem.functionTerms.hydrogenBondTerms.pap_term_2,
                                        ssweightFileName=f"{frags_dir}/ssweight"),
                        DH = partial(openawsem.functionTerms.debyeHuckelTerms.debye_huckel_term, 
                                        chargeFile=f"{frags_dir}/charge.txt")
                        )

    protein.setup_virtual_sites(s)

    #Add DNA-protein interaction forces
    for force_name in open3SPN2.protein_dna_forces:
        print(force_name)
        force = open3SPN2.protein_dna_forces[force_name](dna,protein)
        s.addForce(force)
        forces.update({force_name: force})  
        
    #OpenAWSEM forces with exclusions
    for force_name in openAWSEMforces:
        print(force_name)
        if force_name in ['contact']:
            force = openAWSEMforces[force_name](protein, withExclusion=False,periodic=False)
            print(force_name, "pre-add #Exclusions", force.getNumExclusions())
            open3SPN2.addNonBondedExclusions(dna,force)
            print(force_name, "post-add #Exclusions", force.getNumExclusions())
        elif force_name in ['Excl']:
            force = openAWSEMforces[force_name](protein)
            print(force_name, "pre-add #Exclusions", force.getNumExclusions())
            open3SPN2.addNonBondedExclusions(dna,force)
            print(force_name, "post-add #Exclusions", force.getNumExclusions())
        #continue
        else:
            force = openAWSEMforces[force_name](protein)
        s.addForce(force)
        forces.update({force_name: force})


    if computeQ:
        analysis = dict(Rg = openawsem.functionTerms.biasTerms.rg_term,
                        Q = partial(openawsem.functionTerms.biasTerms.q_value, reference_pdb_file = reference, forceGroup=1))
        for force_name in analysis:
            print(force_name)
            force = analysis[force_name](protein)
            s.addForce(force)
            forces.update({force_name: force})
    return forces
