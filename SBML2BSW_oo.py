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

        self.boundiaryCondition = getBoundaryCondition()
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
#===== END =====    

#==== Function definitions ====
#For now functions, for testing purpose
def react(model):

    #Sistemare name (riga 167 versione python3)
    #sfruttare classe specie

    for i in model.getListOfReactions():

        rc = reaction(i)
        if rc.kin_law != None:

            print ("vvv===========")
            print (rc.kin_law.getListOfParameters())
            print ("^^^===========")

            if len(rc.kin_law.getListOfParameters())==0:
                nomi_parametri = [name.getId() for rc.name in model.getListOfParameters()]
                parameters_in_kineticaw = re.findall(r"[\w']+", react.getKineticLaw().getFormula())
                parameters_in_kineticaw = [x.strip() for x in parameters_in_kineticaw ]
                parameters_in_kineticaw = list(filter( lambda x: x in nomi_parametri, parameters_in_kineticaw ))

                print (react.getName(), "PIKL:", parameters_in_kineticaw)
#==== END ====

react(model)
