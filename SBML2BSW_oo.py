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
#For now functions, for testing purpose
def react(model):

    LEFT = []
    RIGHT = []
    
    #list of all species for now, a species object is created for each
    #object in the model. With models with a lot of reactants it's not
    #good.
    
    Species_list=[]
    Species_ID=[]
    for chem in model.getListOfSpecies():
        a=species(chem)
        Species_list.append(a)
        Species_ID.append(a.ID)

    reactants_vetctor=[]
    #Maniputlating the reactants of each reaction
    for i in model.getListOfReactions():
        tmp_reactants = [0] * len(Species_list)
        rc = reaction(i)
        for i in rc.react_list:
            alias = i.getSpecies()
            index = Species_ID.index(alias)
            if alias in Species_ID:
                spe  = Species_list[index]
                fullname = spe.name
                if fullname == "":
                    fullname = spe.D
                if spe.compartment != "":
                    fullname = fullname + "_in_"+spe.compartment
                sto = i.getStoichiometry()
                tmp_reactants[index] = int(sto)
            else:
                tmp_reactants[index] = 0
                    
            print (tmp_reactants)

            LEFT.append(tmp_reactants)
    print (LEFT)
    
    numpy.savetxt("left_side", LEFT, fmt="%d", delimiter="\t") 
    
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
                    fullname = spe.name
                    if fullname == "":
                        fullname = spe.D
                    if spe.compartment != "":
                        fullname = fullname + "_in_"+spe.compartment
                        sto = i.getStoichiometry()
                        tmp_products[index] = int(sto)
                    else:
                        tmp_products[index] = 0
                else:
                    for index in len(Species_ID)-1:
                        tmp_products[index] = 0
                    
            print (tmp_products)

            RIGHT.append(tmp_products)
    print (RIGHT)
    numpy.savetxt("right_side", RIGHT, fmt="%d", delimiter="\t") 
#==== END ====

react(model)
