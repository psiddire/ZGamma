import ROOT
import array

TMVA_tools = ROOT.TMVA.Tools.Instance()

class FunctorFromMVA(object):
    def __init__(self, name, xml_filename, *variables, **kwargs):
        self.reader    = ROOT.TMVA.Reader( "!Color:Silent=%s:Verbose=%s" % (kwargs.get('silent','T'), kwargs.get('verbose','F')))
        self.var_map   = {}
        self.name      = name
        self.variables = variables
        self.xml_filename = xml_filename
        for var in variables:
            self.var_map[var] = array.array('f',[0]) 
            self.reader.AddVariable(var, self.var_map[var])
        self.reader.BookMVA(name, xml_filename)

    def evaluate_(self):
        return self.reader.EvaluateMVA(self.name)

    #@memo_last
    def __call__(self, **kvars):
        if not ( 
                all(name in self.variables for name in kvars.keys()) and \
                all(name in kvars.keys() for name in self.variables)
        ):
            raise Exception("Wrong variable names. Available variables: %s" % self.variables.__repr__())
        for name, val in kvars.iteritems():
            self.var_map[name][0] = val
        retval = self.evaluate_()
        return retval

