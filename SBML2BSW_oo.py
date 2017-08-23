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
sbml = libsbml.SBMLReader().readSBML("Dismutase.xml")
model = sbml.getModel()
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
        
#===== END =====    

#==== Function definitions ====
def create_folder(OUTPATH):
    """
    Creates the specified output folder
    """
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

        #prepare alphabet
        Alph_element = a.name+"_in_"+a.compartment
        ALPHABET.append(Alph_element)

        #prepare initial amount
        IN_AMOUNT.append(a.conc)

        #prepare M_feed
        if a.spec_con:
            M_FEED.append(1)
            print("1")
        else:
            M_FEED.append(0)
            print("0")
    print (M_FEED)

    for par in model.getListOfParameters():
        Parameter=parameter(par)
#        print (Parameter.value)
        C_VECT.append(Parameter.value)

    
    reactants_vetctor=[]
    #Maniputlating the reactants of each reaction
    for j in model.getListOfReactions():
        tmp_reactants = [0] * len(Species_list)
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
    
    products_vetctor=[]
    #Maniputlating the products of each reaction
    for i in model.getListOfReactions():
        tmp_products = [0] * len(Species_list)
        rc = reaction(i)
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
            
    os.chdir(OUTPUT_FOLDER)
    numpy.savetxt("alphabet", ALPHABET,fmt="%s", delimiter="\t",newline="\t")
    numpy.savetxt("M_0",IN_AMOUNT,fmt="%e", delimiter="\t",newline="\t")
    numpy.savetxt("M_feed",M_FEED,fmt="%d", delimiter="\t",newline="\t")
    numpy.savetxt("c_vector",C_VECT,fmt="%e",delimiter="\t")
    numpy.savetxt("left_side", LEFT, fmt="%d", delimiter="\t") 
    numpy.savetxt("right_side", RIGHT, fmt="%d", delimiter="\t") 
#==== END ====

OUTPUT_FOLDER = "./output"

create_folder(OUTPUT_FOLDER)
react(model,OUTPUT_FOLDER)
