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
# def create_folder(OUTPATH):
#     """
#     Creates the specified output folder
#     """
#     wd = os.path.dirname(__file__)
#     try: os.mkdir(OUTPATH)
#     except: print ("WARNING: directory", OUTPATH, "already exists")

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
            
        
    def react(self,model):

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

        #----Maniputlating the reactants of each reaction
        
        for react in model.getListOfReactions():
            
            tmp_reactants = numpy.zeros(len(self.Species_list))
            tmp_products = numpy.zeros(len(self.Species_list))
            
            rc = reaction(react)
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
                        p = self.model.getParameter(el)
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

                                   self.PARAMS.append(temp)                
                                   #print p.getName(), float(token)
                        else:
                            #print "WARNING: constant value set to 0, parameter:", p.getName()
                            self.PARAMS.append(temp)           
                    else:
                        self.PARAMS.append(p.getValue())
                else:
                    for p in rc.kin_law.getListOfParameters():
                        self.PARAMS.append(p.getValue())


    def save(self):
        os.chdir(self.out)
        numpy.savetxt("alphabet", self.ALPHABET,fmt="%s", delimiter="\t",newline="\t")
        numpy.savetxt("M_0",self.IN_AMOUNT,fmt="%e", delimiter="\t",newline="\t")
        numpy.savetxt("M_feed",self.M_FEED,fmt="%d", delimiter="\t",newline="\t")
        numpy.savetxt("left_side", self.LEFT, fmt="%d", delimiter="\t")
        numpy.savetxt("right_side", self.RIGHT, fmt="%d", delimiter="\t")
        numpy.savetxt("c_vector", numpy.array(self.PARAMS),fmt="%e", delimiter="\t")
        os.chdir(SBML2BSW.wd)

#==== END ====



#==== END ====

if __name__ == '__main__':

    AA = None

    OUTPUT_FOLDER = "./output"

    if len(sys.argv)>1: INPUT_FILE = sys.argv[1]

    sbml = libsbml.SBMLReader().readSBML(INPUT_FILE)
    AA = sbml.getModel()

    SB=SBML2BSW(AA,OUTPUT_FOLDER)
    SB.create_folder()
    SB.react(AA)
    SB.save()
    
    separator()
    
    print("PSA chemical species",SB.ALPHABET)

    separator()

    print("Chemicals Initial Amount",SB.IN_AMOUNT)

    separator()

    print(SB.M_FEED)

    separator()

    print("Reactants:",SB.LEFT)

    separator()

    print("Products:",SB.RIGHT)

    separator()

    print(SB.PARAMS)

    separator()
