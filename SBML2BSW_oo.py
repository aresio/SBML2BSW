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
        self.react_comp = rt.getCompartment()
#===== END =====    

#==== Useful Functions ====
def separator():
    print ("")
    print ("="*100)
    print ("")
#==== END ====

#==== PRINCIPAL CLASS ====

class SBML2BSW():
    """
    This Class uses all the previous ones to calculate and write the
    BSW files
    """
    def __init__(self,model,OUTPATH,lvl):
        self.model=model
        self.out=OUTPATH
        self.lvl=lvl
    wd=os.path.dirname(os.path.abspath(__file__))

    def create_folder(self):
        """
        Creates the specified output folder
        """
        try: os.mkdir(self.out)
        except: print ("WARNING: directory", self.out, "already exists")
            
        
    def react(self,model,verbose=False):

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

        #usefull dictionaries
        self.id2name = {}
        self.Dictionary = {}

        
        for chem in model.getListOfSpecies():
            a=species(chem)
            self.Species_list.append(a)
            self.Species_ID.append(a.ID)
            #--prepare alphabet

            #Pay attention: name is just the label i use to recognise
            #different species.
            name=""
            if a.name != "":
                name=a.name
            else:
                name=a.ID

            if a.compartment != "":
                Alph_element = name+"_in_"+a.compartment

            if verbose:
                print ("* parsing species", Alph_element)
            

            #To avoid double entries
            if Alph_element in self.ALPHABET:
                pass
            else:
                self.ALPHABET.append(Alph_element)

            #dictionary of univocal id and elements in compartments
            #(wich can be double)
            self.id2name[a.ID]=Alph_element

            #-- prepare initial amount
            amount=0

            if self.lvl==3:
                if isnan(a.amount):
                    amount = float(a.conc)
                else:
                    amount=a.amount
            elif self.lvl==2:
                if not a.units:
                    amount=a.amount
                else:
                    amount=a.conc
                
            alias = self.id2name[a.ID]
            index = self.ALPHABET.index(alias)
            
            try:
                self.IN_AMOUNT[index] = amount
            except:
                self.IN_AMOUNT.append(amount)
            
            
            #--prepare M_feed
            if a.spec_con:
                try:
                    self.M_FEED[index]=a.spec_con
                except:
                    self.M_FEED.append(a.spec_con)
            else:
                try:
                    self.M_FEED[index]=0
                except:
                    self.M_FEED.append(0)

            
        #----Invert id2name dictionary
        self.Dictionary = {v: k for k, v in self.id2name.items()}

        #----Maniputlating the reactants of each reaction
        
        for react in model.getListOfReactions():           
            tmp_reactants =[0]*(len(self.ALPHABET))
            tmp_products = [0]*(len(self.ALPHABET))
            create_reverse = False
            
            rc = reaction(react)
            self.REACT_NAME.append(rc.ID)
            
            for reactant in rc.react_list:
                sto=reactant.getStoichiometry()
                alias = self.id2name[reactant.getSpecies()]
                index = self.ALPHABET.index(alias)
                tmp_reactants[index]=int(sto)
            self.LEFT.append(tmp_reactants)

            #----Maniputlating the reactants of each reaction
            for prods in rc.prod_list:
                sto=prods.getStoichiometry()
                alias = self.id2name[prods.getSpecies()]
                index = self.ALPHABET.index(alias)
                tmp_products[index]=int(sto)
            self.RIGHT.append(tmp_products)

            #---- Create the constant vector
            if rc.kin_law != None:
#                print("* direct", rc.ID)
                if len(rc.kin_law.getListOfParameters())==0:
                    nomi_parametri = [par.getId() for par in self.model.getListOfParameters()]
                    parameters_in_kineticlaw = re.findall(r"[\w']+", rc.kin_law.getFormula())
                    parameters_in_kineticlaw = [x.strip() for x in parameters_in_kineticlaw ]
                    #filter object in python3 don't have len attribute, converting to string
                    parameters_in_kineticlaw = list(filter( lambda x: x in nomi_parametri, parameters_in_kineticlaw))
                    print (parameters_in_kineticlaw)
                    if len(parameters_in_kineticlaw)==0:
                        print ("parameters_in_kineticlaw==0")
                        print ("ERROR: can't find any kinetic parameters for reaction", rc.ID)
                        exit(-3)
                        
                    elif len(parameters_in_kineticlaw)==1:
                        p=parameter(self.model.getParameter(parameters_in_kineticlaw[0]))
                        if p.value==0:
                            temp=0
                            if not p_par_const:
                                print ("WARNING: non constant parameter, assignment rule?")
                                if self.model.getListOfRules().get(p.name).isParameter:
                                    tokenized_rule = model.getListOfRules().get(p.name).getFormula()
                                    if tokenized_rule[0:8] == 'stepfunc':
                                        tokenized_rule = tokenized_rule.replace("stepfunc(", "")
                                        tokenized_rule = tokenized_rule.replace(")", "")
                                    tokenized_rule = tokenized_rule.replace(",", "")
                                    tokenized_rule =  tokenized_rule.split()
                                    for token in tokenized_rule:
                                        try:
                                            temp = float(token)
                                            if temp>0:
                                                break
                                        except:
                                            pass
                                        self.PARAMS.append(temp)

                            else:
                                print ("WARNING: constant value set to 0, parameter:", p.name," ",temp)
                                self.PARAMS.append(temp)
                        else:
                            self.PARAMS.append(p.value)
                        if rc.rev:
                            print ("WARNING: detected reversible reaction by getReversible", rc.ID,"but only one parameter defined")
                            print ("Assuming Reverse reaction can't take place (constant=0)")
                            create_reverse = True
                            self.PARAMS.append(0)
                        
                    elif len(list(parameters_in_kineticlaw))==2:
                        if rc.rev:
                            create_reverse = True
                        else:
                            print ("WARNING: detected two parameters in kinetic law of reaction", rc.ID, ", assuming reversible reaction")
                            create_reverse = True
                        for el in parameters_in_kineticlaw:
                            p=parameter(self.model.getParameter(el))
                            if p.value==0:
                                temp = 0
                                if not p.par_const:
                                    print ("WARNING: non constant parameter, assignment rule?")
                                    if self.model.getListOfRules().get(p.name).isParameter:
                                        tokenized_rule = model.getListOfRules().get(p.name).getFormula()
                                        if tokenized_rule[0:8] == 'stepfunc':
                                            tokenized_rule = tokenized_rule.replace("stepfunc(", "")
                                            tokenized_rule = tokenized_rule.replace(")", "")
                                        tokenized_rule = tokenized_rule.replace(",", "")
                                        tokenized_rule =  tokenized_rule.split()
                                        for token in tokenized_rule:
                                            try:
                                                temp = float(token)
                                                if temp>0:
                                                    break
                                            except:
                                                pass

                                            self.PARAMS.append(temp)

                                else:
                                    print ("WARNING: constant value set to 0, parameter:", p.name," ",temp)
                                    self.PARAMS.append(temp)

                            else:
                                print(p.value)
                                self.PARAMS.append(p.value)

                    else:
                        print ("ERROR: too many parameters in kinetic law, aborting")
                        print (list(parameters_in_kineticlaw))
                        print (len(list(parameters_in_kineticlaw)))
                        exit(-3)
#if the length of rc.kin_law.getListOfParameters>0 we
#still can have all the cases af before!

                elif len(rc.kin_law.getListOfParameters())==1:
                    p=parameter(rc.kin_law.getListOfParameters()[0])
                    self.PARAMS.append(p.value)
                    if rc.rev:
                        print ("WARNING: detected reversible reaction by getReversible", rc.ID,"but only one parameter defined")
                        print ("Assuming Reverse reaction can't take place (constant=0)")
                        self.PARAMS.append(0)
                        create_reverse=True
                        
                elif len(rc.kin_law.getListOfParameters()) == 2:
                    if rc.rev:
                        create_reverse=True
                        
                    if len(rc.kin_law.getListOfParameters()) == 2:
                        create_reverse=True
                        
                        
                    for el in rc.kin_law.getListOfParameters():
                        p=parameter(el)
                        if p.value==0:
                            temp = 0
                            if not p.par_const:
                                print ("WARNING: non constant parameter, assignment rule?")
                                if self.model.getListOfRules().get(p.name).isParameter:
                                    tokenized_rule = model.getListOfRules().get(p.name).getFormula()
                                    if tokenized_rule[0:8] == 'stepfunc':
                                        tokenized_rule = tokenized_rule.replace("stepfunc(", "")
                                        tokenized_rule = tokenized_rule.replace(")", "")
                                    tokenized_rule = tokenized_rule.replace(",", "")
                                    tokenized_rule =  tokenized_rule.split()
                                    for token in tokenized_rule:
                                        try:
                                            temp = float(token)
                                            if temp>0:
                                                break
                                        except:
                                            pass
                                        
                                
                                    print ("token =",temp)
                                    self.PARAMS.append(temp)

                            else:   
                                print ("WARNING: constant value set to 0, parameter:", p.name," ",temp)
                                self.PARAMS.append(temp)

                        else:
                            print(p.name," in reaction ",rc.name)
                            self.PARAMS.append(p.value)
                else:
                    print("ERROR: too many parameters")
                    print("Aborting")
                    exit(-3)
                if create_reverse:
#                    print("* reverse",rc.ID)
                    self.REACT_NAME.append(rc.ID+" (reverse)")
                    self.LEFT.append(tmp_products)
                    self.RIGHT.append(tmp_reactants)

    def save(self,verbose=False):
        os.chdir(self.out)
        
        if verbose:
            print(" * Creating 'aphabet' with names of species")
        with open("alphabet", "w") as fo:
            for i in self.ALPHABET:
                fo.write(i+"\t")
        if verbose: print("DONE")
        
        if verbose:
            print(" * Creating 'M_0' with initial state")
        numpy.savetxt("M_0",[self.IN_AMOUNT],fmt="%e",delimiter="\t")
        if verbose:print("DONE")

        if verbose:
            print(" * Creating 'M_feed' file")       
        numpy.savetxt("M_feed",[self.M_FEED],fmt="%d",delimiter="\t")
        if verbose: print("DONE")
        
        if verbose:
            print (" * Creating 'left_side' matrix with reactants")
        numpy.savetxt("left_side", numpy.array(self.LEFT), fmt="%d",newline="\n",delimiter="\t")
        if verbose: print("DONE")
        
        if verbose:
            print (" * Creating 'right_side' matrix with reactants")
        numpy.savetxt("right_side", self.RIGHT, fmt="%d",newline="\n", delimiter="\t")
        if verbose: print("DONE")

        if verbose:
            print (" * Creating 'c_vector' vector with constants")
        numpy.savetxt("c_vector", numpy.array(self.PARAMS),fmt="%e", delimiter="\t")
        if verbose: print("DONE")
        
        if verbose:
            print (" * Creating 'boundaries' matrix for FBA fluxes limits")
        numpy.savetxt("boundaries",self.FLUX_BOUND,fmt="%e",delimiter="\t")
        if verbose: print("DONE")

#==== END ====

if __name__ == '__main__':

    REACT = None

    if len(sys.argv)>1:
        INPUT_FILE = sys.argv[1]
        try:
            OUTPUT_FOLDER = sys.argv[2]
        except:
            OUTPUT_FOLDER = "./output"
    else:
        error_string = """
        No input file specified 
        Positional Argument 1: path to input file
        Positional Argument 2: output folder (optional, default './output')
        """
        print(error_string)
        exit(1)
    
    sbml = libsbml.SBMLReader().readSBML(INPUT_FILE)
    level=sbml.getLevel()
    REACT = sbml.getModel()

    SB=SBML2BSW(REACT,OUTPUT_FOLDER,level)
    SB.create_folder()
    SB.react(REACT)
    SB.save()

    # #-- screen output --
    # separator()
    print("Reaction Names",SB.REACT_NAME)
    print("Reaction Names",len(SB.REACT_NAME))
    # separator()
    # print("PSA chemical species",SB.Dictionary)
    # separator()
    # print("Chemicals Initial Amount",SB.IN_AMOUNT)
    # separator()
    # print("Feed Species",SB.M_FEED)
    # separator()
    print("Parameters Vector",SB.PARAMS)
    print("Parameters Vector",len(SB.PARAMS))
    # separator()
    # #-- END  --
