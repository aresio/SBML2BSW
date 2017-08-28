#!/usr/bin/env python
import re, libsbml, numpy
from libsbml import SBMLNamespaces
import sys, os
import operator
import numpy
import re
from numpy import savetxt
from math import isnan

"""
This is the object Oriented version of SBML2BSW.
It's written in python 3.

THIS VERSION IS UNDER DEVELOPEMENT!
"""

#FIX the REACTACT MATRICE

#===This part is to conduct the test===
#create a model object

model = None
#sbml = libsbml.SBMLReader().readSBML("Dismutase.xml")
#model = sbml.getModel()
#=== END ===


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

#==== Function definitions ====
def create_folder(OUTPATH):
    """
    Creates the specified output folder
    """
    wd = os.path.dirname(__file__)
    try: os.mkdir(OUTPATH)
    except: print ("WARNING: directory", OUTPATH, "already exists")

#For now functions, for testing purpose
def react(model,OUTPATH):

    #----Lists needed for the output---
    LEFT = []
    RIGHT = []
    ALPHABET = []
    C_VECT = []
    IN_AMOUNT = []
    M_FEED = []
    PARAMS = []
    #---- end ----
    
    #list of all species for now, a species object is created for each
    #object in the model. With models with a lot of reactants it's not
    #good.
    
    Species_list=[]
    Species_ID=[]
    
    for chem in model.getListOfSpecies():
        a=species(chem)
        Species_list.append(a)
        Species_ID.append(a.ID)
        
        #--prepare alphabet
        if a.name != "":
            Alph_element = a.name+"_in_"+a.compartment
        else:
            Alph_element = a.ID+"_in_"+a.compartment

        ALPHABET.append(Alph_element)

        #--prepare initial amount

        amount=0
        if isnan(a.amount):
            print("isnan")
            amount=float(a.conc)
            print(amount)
        else:
            print("is a num")
            amount=float(a.amount)
            print(amount)
            
        IN_AMOUNT.append(amount)

        #--prepare M_feed
        if a.spec_con:            
            M_FEED.append(a.spec_con)
        else:
            M_FEED.append(a.spec_con)

    #----Maniputlating the reactants of each reaction
    for j in model.getListOfReactions():
        tmp_reactants = numpy.zeros(len(Species_list))
        rc = reaction(j)
        for i in rc.react_list:
            alias = i.getSpecies()
            index = Species_ID.index(alias)
            if alias in Species_ID:
                spe  = Species_list[index]
                sto = i.getStoichiometry()
                tmp_reactants[index] = int(sto)
            else:
                tmp_reactants[index] = 0

        LEFT.append(tmp_reactants)
    
    #----Maniputlating the products of each reaction
    
        tmp_products =  numpy.zeros(len(Species_list))
        if rc.prod_list != "":
            for i in rc.prod_list:
                alias = i.getSpecies()
                index = Species_ID.index(alias)
                if alias in Species_ID:
                    spe  = Species_list[index]
                    if spe.compartment != "":
                        sto = i.getStoichiometry()
                        tmp_products[index] = int(sto)
                    else:
                        tmp_products[index] = 0
                else:
                    for index in len(Species_ID)-1:
                        tmp_products[index] = 0

            RIGHT.append(tmp_products)

        if rc.kin_law != None:
            if len(rc.kin_law.getListOfParameters())==0:
                nomi_parametri = [par.getId() for par in model.getListOfParameters() ]
                parameters_in_kineticaw = re.findall(r"[\w']+", rc.kin_law.getFormula())
                parameters_in_kineticaw = [x.strip() for x in parameters_in_kineticaw ]
                #filter object in python3 don't have len attribute, converting to string
                parameters_in_kineticaw = list(filter( lambda x: x in nomi_parametri, parameters_in_kineticaw))
                print(parameters_in_kineticaw)
                if len(parameters_in_kineticaw)==0:
                    print ("parameters_in_kineticaw==0")
                    print ("ERROR: can't find any kinetic parameters for reaction", rc.name)
                    exit(-3)
                elif rc.rev:
                    #print "WARNING: detected reversible reaction by getReversible", reaction_name
                    create_reverse = True
                elif len(list(parameters_in_kineticaw))==2:
                    print ("WARNING: detected two parameters in kinetic law of reaction", rc.name, ", assuming reversible reaction")
                    create_reverse = True
                elif len(list(parameters_in_kineticaw))==1:
                    pass
                else:
                    print ("ERROR: too many parameters in kinetic law, aborting")
                    exit(-3)
                    
                for el in parameters_in_kineticaw:                   
                    p = model.getParameter(el)
                    
                    print (len(parameters_in_kineticaw), p)
                    
                    if p.getValue()==0:
                        if not p.constant:
                            print ("WARNING: non constant parameter, assignment rule?")
                            if model.getListOfRules().get(p.getName()).isParameter():
                                # print " * Rule for parameter", p.getName(), "detected"
                                # print " * Rule implemented as", self.model.getListOfRules().get("k1").getFormula()
                                tokenized_rule = model.getListOfRules().get("k1").getFormula()
                                if tokenized_rule[0:8] == 'stepfunc':
                                    tokenized_rule = tokenized_rule.replace("stepfunc(", "")
                                    tokenized_rule = tokenized_rule.replace(")", "")
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
                                    #print token, "is NAN"
                                    
                                PARAMS.append(temp)                
                                #print p.getName(), float(token)
                        else:
                            #print "WARNING: constant value set to 0, parameter:", p.getName()
                            PARAMS.append(temp)           
                else:
                    PARAMS.append(p.getValue())
            else:
                for p in rc.kin_law.getListOfParameters():
                    PARAMS.append(p.getValue())
    
    os.chdir(OUTPUT_FOLDER)
    numpy.savetxt("c_vector", numpy.array(PARAMS),fmt="%e", delimiter="\t")
    numpy.savetxt("alphabet", ALPHABET,fmt="%s", delimiter="\t",newline="\t")
    numpy.savetxt("M_0",IN_AMOUNT,fmt="%e", delimiter="\t",newline="\t")
    numpy.savetxt("M_feed",M_FEED,fmt="%d", delimiter="\t",newline="\t")
    numpy.savetxt("left_side", LEFT, fmt="%d", delimiter="\t") 
    numpy.savetxt("right_side", RIGHT, fmt="%d", delimiter="\t") 
#==== END ====

if __name__ == '__main__':
    OUTPUT_FOLDER = "./output"

    if len(sys.argv)>1: INPUT_FILE = sys.argv[1]

    sbml = libsbml.SBMLReader().readSBML(INPUT_FILE)
    model = sbml.getModel()

    create_folder(OUTPUT_FOLDER)
    react(model,OUTPUT_FOLDER)
