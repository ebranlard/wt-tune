import os

class WrongFormatError(Exception):
    pass

class EmptyFileError(Exception):
    pass

class BrokenFormatError(Exception):
    pass

try: #Python3
    FileNotFoundError=FileNotFoundError
except NameError: # Python2
    FileNotFoundError = IOError

class File(object):
    def __init__(self,filename=None):
        if filename:
            ### If there is a new filename, replace the object variable
            self.filename = filename
            ### If the filename is provided, read the file
            self.read()
        else:
            self.filename = None


    def read(self, filename=None):
        if filename:
            self.filename = filename

        if self.filename:
            if not os.path.isfile(self.filename):
                raise OSError(2,'File not found:',self.filename)
            if os.stat(self.filename).st_size == 0:
                raise EmptyFileError('File is empty:',self.filename)
            # Calling children function
            self._read()
        else:  
            raise Exception('No filename provided')

    def write(self, filename=None):
        if filename:
            self.filename = filename

        if self.filename:
            # Calling children function
            self._write()
        else:
            raise Exception('No filename provided')

    def toDataFrame(self):
        return self._toDataFrame()

    
    def _read(self):
        raise NotImplementedError("Method must be implemented in the subclass")

    def _write(self):
        raise NotImplementedError("Method must be implemented in the subclass")

    def _toDataFrame(self):
        raise NotImplementedError("Method must be implemented in the subclass")

    @staticmethod
    def defaultExtension():
        raise NotImplementedError("Method must be implemented in the subclass")

    @staticmethod
    def formatName():
        raise NotImplementedError("Method must be implemented in the subclass")

    @classmethod
    def isRightFormat(cls,filename):
        """ Tries to open a file, return true and the file if it succeeds """
        #raise NotImplementedError("Method must be implemented in the subclass")
        try:
            F=cls(filename=filename)
            return True,F
        except MemoryError:
            raise
        except WrongFormatError:
            return False,None
        except:
            raise




    def test_ascii(self):
        # compare ourselves (assuming read has occured) with what we write
        # --- First re-read original as ascii
        # comparison is done ignoring multiple white spaces for now
        with open(self.filename, 'r') as f1:
            lines1 = f1.read().splitlines();
            lines1 = '|'.join([l.replace('\t',' ').strip() for l in lines1])
            lines1 = ' '.join(lines1.split())

        # --- Then test write function (assuming read)
        try:
            filename_out = self.filename+'_TMP'
            self.write(filename_out)
        except Exception as e:
            raise Exception('Error writing what we read\n'+e.args[0])

        # --- Third read what we wrote as ascii
        with open(filename_out, 'r') as f2:
            lines2 = f2.read().splitlines();
            lines2 = '|'.join([l.replace('\t',' ').strip() for l in lines2])
            lines2 = ' '.join(lines2.split())

        # --- Fourth re-read what we wrote
        try:
            self.read(filename_out)
        except Exception as e:
            raise Exception('Error reading what we wrote\n'+e.args[0])

        # Last, we perform the ascii comparison
        if lines1 == lines2:
            print('[ OK ] '+self.filename)
            os.remove(filename_out)
        else:
            print('[FAIL] '+self.filename)


# --------------------------------------------------------------------------------}
# --- Helper functions
# --------------------------------------------------------------------------------{
def isBinary(filename):
    with open(filename, 'r') as f:
        try:
            f.readline()
            return False
        except UnicodeDecodeError:
            return True
