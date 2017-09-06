#!/usr/bin/env python
import re, libsbml, numpy
from libsbml import SBMLNamespaces
import sys, os
import operator
import numpy
import re
from numpy import savetxt
from math import isnan

#===== Constructors Definitions =====
class SBML_Object(object):
    """
    This class define some variables used by all the other classes
    describing a SBML-object
    """
    def __init__(self,obj):
        self.ID = obj.getId()
        self.meta_ID = obj.getMetaId()
        self.name = obj.getName()

class species(SBML_Object):
    """
    This class describes the chemical species defined in the
    model. Inherits form SBML_object
    """
    def __init__(self,spe):
        super(species,self).__init__(spe)

        self.boundaryCondition = spe.getBoundaryCondition()
        self.compartment = spe.getCompartment()
        self.spec_con = spe.getConstant()
        #Da verificare se lasciare le units qui o tenerla nella classe con
        #il metodo per lanciare il tutto
        self.units = spe.getSubstanceUnits()
        self.conc = spe.getInitialConcentration()
        self.amount= spe.getInitialAmount()

class parameter(SBML_Object):
    """
    This constructor describes reaction constants (k)
    Inherits from SBML_Object
    """
    def __init__(self,par):
        super(parameter,self).__init__(par)
        self.par_const = par.getConstant()
        self.value = par.getValue()
        self.sbo = par.getSBOTermID()

class compartment(SBML_Object):
    """
    This Constructor describes a compartment
    Inherits from SBML_Object
    """
    def __init__(self,comp):
        super(compartment,self).__init__(comp)
        self.comp_const = comp.getConstant()
        self.size = comp.getSize()
        self.dims = comp.getSpatialDimensions()

class kin_law(SBML_Object):
    """
    This Constructor describes a compartment
    Inherits from SBML_Object
    """
    def __init__(self,kl):
        super(kin_law,self).__init__(kl)
        self.comp_const = kl.getConstant()
        self.size = kl.getSize()
        self.dims = kl.getSpatialDimensions()
        
class reaction(SBML_Object):
    """
    This constructor describes reactions.  To describe a reaction you
    need istances of the previus classes
    Inherits from SBML_Object
    """
    def __init__(self,rt):
        super(reaction,self).__init__(rt)
        self.kin_law = rt.getKineticLaw()
        self.react_list = rt.getListOfReactants()
        self.prod_list = rt.getListOfProducts()
        self.rev = rt.getReversible()

#===== END =====    

#==== Useful Functions ====
def separator():
    print ("")
    print ("="*100)
    print ("")
#==== END ====

#==== PRINCIPAL CLASS ====

class SBML2BSW():
    def __init__(self,model,OUTPATH):
        self.model=model
        self.out=OUTPATH

    wd=os.path.dirname(os.path.abspath(__file__))

    def create_folder(self):
        """
        Creates the specified output folder
        """
        try: os.mkdir(self.out)
        except: print ("WARNING: directory", self.out, "already exists")
            
        
    def react(self,model,use_fba=False):

        #list needed
        self.ALPHABET = []
        self.Species_list = []
        self.Species_ID = []
        self.IN_AMOUNT = []
        self.M_FEED = []
        self.LEFT = []
        self.RIGHT = []
        self.PARAMS = []
        self.REACT_NAME = []
        self.BOUNDIARIES = []
        self.FLUX_BOUND = []
        
        #-- Am I using fba?
        #NB: comes directly from the old version
        if use_fba:
            # FBA stuff
            fbc = self.model.getPlugin("fbc")
            print (" * Detected objectives:", fbc.getNumObjectives())
            for o in xrange(fbc.getNumObjectives()):
                OF = ObjectiveFunction()
                SBOf = fbc.getObjective(o)
                OF.name = SBOf.getName()
                OF.ID = SBOf.getId()
                OF.type = SBOf.getType()
                print (" * Objective function detected:", OF.type, OF.ID)
                # for ofl_ in range(SBOf.getNumFluxObjectives()):
                #    print ofl_
                for fo in xrange(SBOf.getNumFluxObjectives()):
                    fluxobj = SBOf.getFluxObjective(fo)
                    print (" * Coefficient", fo, "for reaction", fluxobj.getReaction(), "is equal to", fluxobj.getCoefficient())
        
        
        for chem in model.getListOfSpecies():
            a=species(chem)
            self.Species_list.append(a)
            self.Species_ID.append(a.ID)
            #--prepare alphabet
            if a.name != "":
                Alph_element = a.name+"_in_"+a.compartment
            else:
                Alph_element = a.ID+"_in_"+a.compartment

            self.ALPHABET.append(Alph_element)

            #-- prepare initial amount
            
            
            amount=0
            if isnan(a.amount):
                amount=float(a.conc)
            else:
                amount=float(a.amount)

            self.IN_AMOUNT.append(amount)
                
            #--prepare M_feed
            
            if a.spec_con:
                self.M_FEED.append(a.spec_con)
            else:
                self.M_FEED.append(int(0))

        #----Maniputlating the reactants of each reaction
        
        for react in model.getListOfReactions():
            
            if use_fba:
                RFBCplg = react.getPlugin('fbc')
                lfbid = RFBCplg.getLowerFluxBound()
                ufbid = RFBCplg.getUpperFluxBound()
                self.FLUX_BOUND.append( [PARAM_D[lfbid]['value'], PARAM_D[ufbid]['value']] )

            
            tmp_reactants =[0]*(len(self.Species_list))
            tmp_products = [0]*(len(self.Species_list))
            create_reverse = False
            rc = reaction(react)
            self.REACT_NAME.append(rc.ID)
            for single_react in rc.react_list:
                react_species = single_react.getSpecies()
                index = self.Species_ID.index(react_species)
                if react_species in self.Species_ID:
                    spe  = self.Species_list[index]
                    sto = single_react.getStoichiometry()
                    tmp_reactants[index] = int(sto)
                # else:
                #     tmp_reactants[index] = 0
                        
            self.LEFT.append(tmp_reactants)

            #----Maniputlating the reactants of each reactio
            for prods in rc.prod_list:
                prod_species = prods.getSpecies()
                index = self.Species_ID.index(prod_species)
                if prod_species in self.Species_ID:
                    spe = self.Species_list[index]
                    sto = prods.getStoichiometry()
                    tmp_products[index] = int(sto)
            self.RIGHT.append(tmp_products)

            #---- Create the constant vector
            if rc.kin_law != None:
                if len(rc.kin_law.getListOfParameters())==0:
                    nomi_parametri = [par.getId() for par in self.model.getListOfParameters() ]
                    parameters_in_kineticaw = re.findall(r"[\w']+", rc.kin_law.getFormula())
                    parameters_in_kineticaw = [x.strip() for x in parameters_in_kineticaw ]
                    #filter object in python3 don't have len attribute, converting to string
                    parameters_in_kineticaw = list(filter( lambda x: x in nomi_parametri, parameters_in_kineticaw))
                    
                    if len(parameters_in_kineticaw)==0:
                        print ("parameters_in_kineticaw==0")
                        print ("ERROR: can't find any kinetic parameters for reaction", rc.ID)
                        exit(-3)
                    elif rc.rev:
                        print ("WARNING: detected reversible reaction by getReversible", rc.ID)
                        create_reverse = True
                        #self.rev_par=1/(float(self.model.getParameter(parameters_in_kineticaw[0]).getValue()))
                    elif len(list(parameters_in_kineticaw))==2:
                        print ("WARNING: detected two parameters in kinetic law of reaction", rc.ID, ", assuming reversible reaction")
                        create_reverse = True
                        #self.rev_par=self.model.getParameter(parameters_in_kineticaw[1]).getValue()
                    elif len(list(parameters_in_kineticaw))==1:
                        pass
                    else:
                        print ("ERROR: too many parameters in kinetic law, aborting")
                        print (list(parameters_in_kineticaw))
                        exit(-3)
                    
                    for el in parameters_in_kineticaw:                   
                        p=parameter(self.model.getParameter(el))
                        temp=0
                        if p.value==0:
                            if not p.par_const:
                                print ("reaction",rc.name)
                                print ("WARNING: non constant parameter, assignment rule?")
                                #-- To Be Verified--------
                                if model.getListOfRules().get(1).isParameter():
                                    print ("vvvv=======")
                                    print ("*Rule for parameter",p.name,"detected")
                                    print ("^^^^=======")
                                    tokenized_rule = model.getListOfRules().get(1).getFormula()
                                    if tokenized_rule[0:8] == 'stepfunc':
                                        tokenized_rule = tokenized_rule.replace("stepfunc(", "")
                                        tokenized_rule = tokenized_rule.replace(")", "")
                                    print ("* tokenized_rule =" ,tokenized_rule)
                                    tokenized_rule = tokenized_rule.replace(",", "")
                                    tokenized_rule =  tokenized_rule.split()
                                    temp = 0
                                    for token in tokenized_rule:
                                        try:
                                            temp = float(token)
                                            if temp>0:
                                                break
                                        except:
                                            pass
                                    print ("token =",temp)
                                    self.PARAMS.append(temp)
                                #-- To Be Verified-------
                            else:
                                print ("WARNING: constant value set to 0, parameter:", p.name)
                                self.PARAMS.append(temp)
                        else:
                            print (float(p.value))
                            self.PARAMS.append(p.value)
                else:
                    for p in rc.kin_law.getListOfParameters():
                        self.PARAMS.append(p.value)

                if create_reverse:
                    self.REACT_NAME.append(rc.ID+" (reverse)")
                    self.LEFT.append(tmp_products)
                    self.RIGHT.append(tmp_reactants)
                    # verify if reverse parameters are needed
                    # print ("Reverse constant = ",self.rev_par)
                    # self.PARAMS.append(self.rev_par)


    def save(self):
        os.chdir(self.out)
        numpy.savetxt("alphabet", self.ALPHABET,fmt="%s", delimiter="\t",newline="\t")
        numpy.savetxt("M_0",self.IN_AMOUNT,fmt="%e", delimiter="\t",newline="\t")
        numpy.savetxt("M_feed",self.M_FEED,fmt="%d", delimiter="\t",newline="\t")
        numpy.savetxt("left_side", self.LEFT, fmt="%d", delimiter="\t")
        numpy.savetxt("right_side", numpy.array(self.RIGHT), fmt="%d", delimiter="\t")
        numpy.savetxt("c_vector", numpy.array(self.PARAMS),fmt="%e", delimiter="\t")
        numpy.savetxt("boundaries",self.FLUX_BOUND,fmt="%e",delimiter="\t")
        os.chdir(SBML2BSW.wd)

#==== END ====



#==== END ====

if __name__ == '__main__':

    REACT = None

    OUTPUT_FOLDER = "./output"

    if len(sys.argv)>1: INPUT_FILE = sys.argv[1]

    sbml = libsbml.SBMLReader().readSBML(INPUT_FILE)
    REACT = sbml.getModel()

    SB=SBML2BSW(REACT,OUTPUT_FOLDER)
    SB.create_folder()
    SB.react(REACT)
    SB.save()

    #-- screen output --
    separator()
    print("Reaction Names",SB.REACT_NAME)
    separator()
    print("PSA chemical species",SB.ALPHABET)
    separator()
    print("Chemicals Initial Amount",SB.IN_AMOUNT)
    separator()
    print("Feed Species",SB.M_FEED)
    separator()
    print("Reactants:",SB.LEFT)
    separator()
    print("Products:",SB.RIGHT)
    separator()
    print("Parameters Vector",SB.PARAMS)
    separator()
    #-- END  --
