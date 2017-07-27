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
        
#===== END =====    

#==== Function definitions ====
#For now functions, for testing purpose
def react(model):
    #dictionary id-name
    id2name = {}
    
    #list of all species
    Species_list=[]
    for chem in model.getListOfSpecies():
        Species_list.append(chem)
   
    #Maniputlating the reactants of each reaction
    for i in model.getListOfReactions():
        rc = reaction(i)
        for i in rc.react_list:
            alias = i.getSpecies()
            #vvv==Messed up, to be cleaned
            #according to the previous versions
            for j in Species_list:
                sp=species(j)
                if str(sp.ID) == alias:
                    fullname = sp.name
                    if fullname == "":
                        fullname = sp.ID
                    if sp.compartment != "":
                        fullname = fullname + "_in_"+sp.compartment
                    id2name[sp.ID] = fullname
    print (id2name)
                    
#==== END ====

react(model)
