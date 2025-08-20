class Properties:
        def __init__(self,propfilename):
                self.dict={}
                infile=open(propfilename,'r')
                for line in infile:
                        if (line) and (line[0]!='#') and ('=' in line):
                                (key,val)=line.split('=')
                                self.dict[key]=val.strip()
                return

        def getport(self):
                if 'port' in self.dict:
                        return self.dict['port']

        def getusername(self): 
                if 'username' in self.dict:
                        return self.dict['username']

        def getpassword(self): 
                if 'password' in self.dict:
                        return self.dict['password']

        def getdatabase(self): 
                if 'database' in self.dict:
                        return self.dict['database']

        def getipaddress(self): 
                if 'ip' in self.dict:
                        return self.dict['ip']

